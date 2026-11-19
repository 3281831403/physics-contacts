"""Contact detection and impulse helpers for stacked rigid bodies.

The public API is deliberately dict-based so every intermediate value is
visible to the tests.  A body is an inline dict with the fields ``mass``,
``inv_mass``, ``pos``, ``vel``, ``radius``, ``friction`` and
``restitution``.  A contact is an inline dict with ``a``, ``b``, ``depth``
and ``normal``.

Only the Python standard library and pytest are used by this repository.
"""

import math
from dataclasses import dataclass

POSITION_TOLERANCE = 1e-6
VELOCITY_TOLERANCE = 1e-6
DEFAULT_DT = 1.0 / 60.0


@dataclass(frozen=True)
class Tolerances:
    """Numeric tolerances shared by the contact solver."""

    position: float = POSITION_TOLERANCE
    velocity: float = VELOCITY_TOLERANCE


TOLERANCES = Tolerances()


class PhysicsError(Exception):
    """Raised when a physics operation cannot produce finite, valid state."""


def _clone_body(body) -> dict:
    """Return a shallow-enough copy of a body with fresh vector lists."""
    clone = dict(body)
    if "pos" in clone:
        clone["pos"] = list(clone["pos"])
    if "vel" in clone:
        clone["vel"] = list(clone["vel"])
    return clone


def _clone_contact(contact) -> dict:
    """Return a fresh contact dict with a fresh normal list."""
    clone = dict(contact)
    clone["normal"] = list(clone["normal"])
    return clone


def _is_finite(value) -> bool:
    """True for every finite number, recursing into lists/dicts."""
    if isinstance(value, bool):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(value)
    if isinstance(value, list):
        return all(_is_finite(item) for item in value)
    if isinstance(value, dict):
        return all(_is_finite(item) for item in value.values())
    return True


def detect(contacts, bodies, tolerance) -> list[dict]:
    """Filter and order candidate contacts.

    ``contacts`` contains mappings with integer ``a``/``b`` body indices,
    float ``depth`` and a 2D ``normal`` pointing from ``a`` to ``b``.
    Contacts with ``depth < tolerance`` are discarded.  The returned list
    contains fresh contact dicts sorted by depth descending; equal depths
    keep their input order.  Neither input list is mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    kept = [
        _clone_contact(contact)
        for contact in contacts
        if contact["depth"] >= tolerance
    ]
    # A stable sort by descending depth preserves input order on ties, so
    # the result no longer depends on the order detection emitted contacts.
    kept.sort(key=lambda contact: contact["depth"], reverse=True)
    return kept


def resolve_penetration(bodies, contact, slop) -> dict:
    """Remove penetration from one contact.

    Returns ``{"bodies": ..., "applied": bool, "correction": float}``.
    For ``depth <= slop`` nothing moves.  Otherwise ``correction = depth -
    slop`` and it is split along ``normal`` in inverse-mass proportion, so
    the mass-weighted centre of mass does not move.  A body with
    ``inv_mass == 0`` is immovable.  Velocities are not changed.  Inputs
    are not mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    new_bodies = [_clone_body(body) for body in bodies]
    depth = contact["depth"]

    if depth <= slop:
        return {"bodies": new_bodies, "applied": False, "correction": 0.0}

    correction = depth - slop
    index_a = contact["a"]
    index_b = contact["b"]
    inv_mass_a = new_bodies[index_a]["inv_mass"]
    inv_mass_b = new_bodies[index_b]["inv_mass"]
    total_inv_mass = inv_mass_a + inv_mass_b

    if total_inv_mass <= 0.0:
        # Both bodies are immovable: the penetration cannot be removed.
        return {"bodies": new_bodies, "applied": False, "correction": 0.0}

    normal = contact["normal"]
    # Inverse-mass weighting keeps the mass-weighted centre of mass fixed:
    # a moves against the normal, b moves along it, and m_a*d_a + m_b*d_b = 0.
    share_a = inv_mass_a / total_inv_mass
    share_b = inv_mass_b / total_inv_mass
    for axis in range(len(normal)):
        component = normal[axis] * correction
        new_bodies[index_a]["pos"][axis] -= component * share_a
        new_bodies[index_b]["pos"][axis] += component * share_b

    return {"bodies": new_bodies, "applied": True, "correction": correction}


def apply_friction(body, normal, tangent_impulse, mu) -> dict:
    """Apply a clamped tangent impulse to a body.

    ``normal`` is ``{"direction": [nx, ny], "jn": float}``.  The tangent is
    the unit normal rotated by +90 degrees, ``[-ny, nx]``.  The requested
    scalar ``tangent_impulse`` is clamped to ``[-mu * jn, mu * jn]`` before
    the velocity change ``tangent * impulse * inv_mass`` is added.  Returns
    ``{"body": new_body, "tangent_impulse": clamped, "clamped": bool}``.
    Raises ``PhysicsError`` for a zero normal direction or negative
    ``mu``/``jn``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    direction = normal["direction"]
    jn = normal["jn"]
    nx, ny = direction
    length = math.hypot(nx, ny)
    if length == 0.0:
        raise PhysicsError("normal direction must not be zero")
    if mu < 0.0:
        raise PhysicsError("friction coefficient mu must not be negative")
    if jn < 0.0:
        raise PhysicsError("normal impulse jn must not be negative")

    limit = mu * jn
    clamped_impulse = max(-limit, min(limit, tangent_impulse))
    was_clamped = clamped_impulse != tangent_impulse

    tangent = [-ny, nx]
    new_body = _clone_body(body)
    inv_mass = new_body["inv_mass"]
    for axis in range(2):
        new_body["vel"][axis] += tangent[axis] * clamped_impulse * inv_mass

    return {
        "body": new_body,
        "tangent_impulse": clamped_impulse,
        "clamped": was_clamped,
    }


def combine_materials(a, b) -> dict:
    """Combine two materials.

    ``friction = sqrt(max(0, a.friction) * max(0, b.friction))`` and
    ``restitution = max(a.restitution, b.restitution)``.  Returns
    ``{"friction": ..., "restitution": ...}``.  Inputs are not mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    friction = math.sqrt(max(0.0, a["friction"]) * max(0.0, b["friction"]))
    restitution = max(a["restitution"], b["restitution"])
    return {"friction": friction, "restitution": restitution}


def solve_iterations(bodies, contacts, iterations) -> dict:
    """Run sequential impulse iterations for the supplied contacts.

    Each pass resolves normal velocities, then returns a fresh body list and
    ``{"bodies": ..., "residual": ..., "iterations": ..., "converged": ...}``.
    ``residual`` is the largest remaining closing speed after a pass.
    Iteration stops when ``residual <= 1e-6`` or ``iterations`` passes have
    been made.  The input bodies are not mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    new_bodies = [_clone_body(body) for body in bodies]
    # A depth-sorted, stable contact order makes the solver independent of
    # the order contacts happened to be detected in.
    ordered = sorted(contacts, key=lambda contact: contact["depth"], reverse=True)

    def closing_speed(contact) -> float:
        """Approaching speed along the normal; positive means closing.

        The normal points from ``a`` to ``b``, so bodies approach when
        ``a`` moves along the normal faster than ``b`` does.
        """
        index_a = contact["a"]
        index_b = contact["b"]
        normal = contact["normal"]
        velocity_a = new_bodies[index_a]["vel"]
        velocity_b = new_bodies[index_b]["vel"]
        relative = (
            velocity_a[0] - velocity_b[0],
            velocity_a[1] - velocity_b[1],
        )
        return relative[0] * normal[0] + relative[1] * normal[1]

    residual = max(
        (closing_speed(contact) for contact in ordered), default=0.0
    )
    residual = max(0.0, residual)
    passes = 0
    converged = residual <= VELOCITY_TOLERANCE

    while passes < iterations and not converged:
        for contact in ordered:
            speed = closing_speed(contact)
            if speed <= 0.0:
                # Already separating: an impulse would pull the bodies back.
                continue
            index_a = contact["a"]
            index_b = contact["b"]
            inv_mass_a = new_bodies[index_a]["inv_mass"]
            inv_mass_b = new_bodies[index_b]["inv_mass"]
            total_inv_mass = inv_mass_a + inv_mass_b
            if total_inv_mass <= 0.0:
                continue
            # Normal impulse that exactly cancels the closing velocity.
            impulse = speed / total_inv_mass
            normal = contact["normal"]
            for axis in range(2):
                new_bodies[index_a]["vel"][axis] -= (
                    normal[axis] * impulse * inv_mass_a
                )
                new_bodies[index_b]["vel"][axis] += (
                    normal[axis] * impulse * inv_mass_b
                )

        passes += 1
        residual = max(
            (max(0.0, closing_speed(contact)) for contact in ordered),
            default=0.0,
        )
        if residual <= VELOCITY_TOLERANCE:
            converged = True

    return {
        "bodies": new_bodies,
        "residual": residual,
        "iterations": passes,
        "converged": converged,
    }


def sleep_check(body, lin_threshold, ang_threshold, frames) -> bool:
    """Track consecutive low-motion frames and report sleep.

    A body is low-motion when ``abs(vel) <= lin_threshold`` and
    ``abs(ang_vel) <= ang_threshold`` (both inclusive).  Each call updates
    ``body["sleep_frames"]``: increment when low-motion, reset to zero
    otherwise.  Return True once the counter reaches ``frames``.
    ``frames`` must be at least 1; otherwise raise ``PhysicsError``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    if frames < 1:
        raise PhysicsError("frames must be at least 1")

    velocity = body["vel"]
    angular_velocity = body.get("ang_vel", 0.0)
    linear_low = all(abs(component) <= lin_threshold for component in velocity)
    angular_low = abs(angular_velocity) <= ang_threshold

    if linear_low and angular_low:
        body["sleep_frames"] = body.get("sleep_frames", 0) + 1
    else:
        body["sleep_frames"] = 0

    return body["sleep_frames"] >= frames


def wake(body, impulse) -> dict:
    """Wake a sleeping body with an impulse.

    Returns a new body dict with ``sleep_frames == 0`` and
    ``sleeping == False``.  The existing velocity is kept and the impulse
    changes it by ``impulse * inv_mass``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    new_body = _clone_body(body)
    inv_mass = new_body["inv_mass"]
    for axis in range(len(impulse)):
        new_body["vel"][axis] += impulse[axis] * inv_mass
    new_body["sleep_frames"] = 0
    new_body["sleeping"] = False
    return new_body


def step_batch(bodies, contacts, substeps) -> dict:
    """Advance a batch by one fixed timestep.

    Total step duration is 1/60 s, split into ``substeps`` equal slices of
    ``1/60 / substeps``.  Bodies are integrated in place.  Returns
    ``{"bodies": bodies, "dt": 1/60, "substeps": substeps}``.
    ``substeps <= 0`` raises ``PhysicsError``.  Every intermediate value
    must be finite; if any body or contact produces NaN/inf, restore all
    bodies to their exact pre-call state and raise ``PhysicsError``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    if substeps <= 0:
        raise PhysicsError("substeps must be a positive integer")

    # Snapshot every body before touching anything so a non-finite value at
    # any point rolls the whole batch back to its exact pre-call state.
    snapshot = [_clone_body(body) for body in bodies]

    def rollback_and_raise() -> None:
        for original, current in zip(snapshot, bodies):
            current.clear()
            current.update(original)
        raise PhysicsError("non-finite value encountered; batch rolled back")

    if not all(_is_finite(body) for body in bodies):
        rollback_and_raise()
    if not all(_is_finite(contact) for contact in contacts):
        rollback_and_raise()

    substep_dt = DEFAULT_DT / substeps
    try:
        for _ in range(substeps):
            for body in bodies:
                for axis in range(len(body["vel"])):
                    body["pos"][axis] += body["vel"][axis] * substep_dt
                if not _is_finite(body):
                    rollback_and_raise()
    except PhysicsError:
        raise
    except Exception:
        rollback_and_raise()

    return {"bodies": bodies, "dt": DEFAULT_DT, "substeps": substeps}
