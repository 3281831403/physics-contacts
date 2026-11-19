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


def detect(contacts, bodies, tolerance) -> list[dict]:
    """Filter and order candidate contacts.

    ``contacts`` contains mappings with integer ``a``/``b`` body indices,
    float ``depth`` and a 2D ``normal`` pointing from ``a`` to ``b``.
    Contacts with ``depth < tolerance`` are discarded.  The returned list
    contains fresh contact dicts sorted by depth descending; equal depths
    keep their input order.  Neither input list is mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    raise NotImplementedError("detect is not implemented")


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
    raise NotImplementedError("resolve_penetration is not implemented")


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
    raise NotImplementedError("apply_friction is not implemented")


def combine_materials(a, b) -> dict:
    """Combine two materials.

    ``friction = sqrt(max(0, a.friction) * max(0, b.friction))`` and
    ``restitution = max(a.restitution, b.restitution)``.  Returns
    ``{"friction": ..., "restitution": ...}``.  Inputs are not mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    raise NotImplementedError("combine_materials is not implemented")


def solve_iterations(bodies, contacts, iterations) -> dict:
    """Run sequential impulse iterations for the supplied contacts.

    Each pass resolves normal velocities, then returns a fresh body list and
    ``{"bodies": ..., "residual": ..., "iterations": ..., "converged": ...}``.
    ``residual`` is the largest remaining closing speed after a pass.
    Iteration stops when ``residual <= 1e-6`` or ``iterations`` passes have
    been made.  The input bodies are not mutated.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    raise NotImplementedError("solve_iterations is not implemented")


def sleep_check(body, lin_threshold, ang_threshold, frames) -> bool:
    """Track consecutive low-motion frames and report sleep.

    A body is low-motion when ``abs(vel) <= lin_threshold`` and
    ``abs(ang_vel) <= ang_threshold`` (both inclusive).  Each call updates
    ``body["sleep_frames"]``: increment when low-motion, reset to zero
    otherwise.  Return True once the counter reaches ``frames``.
    ``frames`` must be at least 1; otherwise raise ``PhysicsError``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    raise NotImplementedError("sleep_check is not implemented")


def wake(body, impulse) -> dict:
    """Wake a sleeping body with an impulse.

    Returns a new body dict with ``sleep_frames == 0`` and
    ``sleeping == False``.  The existing velocity is kept and the impulse
    changes it by ``impulse * inv_mass``.

    Tolerance: position 1e-6, velocity 1e-6.
    """
    raise NotImplementedError("wake is not implemented")


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
    raise NotImplementedError("step_batch is not implemented")

