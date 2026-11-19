# physics-contacts

## 题面

堆叠箱子一多，物理就开始抖：穿透修正只推一个物体，双质量场景动量不守恒；忽略 slop，静止物体持续被微推；摩擦冲量不钳制，切向冲量超过 mu*jn 后物体加速；材质合成用相加，摩擦系数超过 1；接触排序按检测顺序，求解结果依赖输入顺序；迭代不判收敛，固定跑满次数且残差不下降；休眠只看速度不看角速度，旋转物体被错误休眠；唤醒时把速度清零，物体被击中后不动；子步 dt 不除 substeps，模拟速度随子步数变化；批量出现 NaN 时留下部分更新，物理状态被污染。physics/contacts.py 是核心，tests/test_contacts.py 15个红。要修的是：按逆质量分配穿透修正、尊重 slop、钳制摩擦锥、摩擦取几何平均且弹性取最大、按深度排序接触、迭代判收敛、休眠同时检查角速度、唤醒保留速度、子步缩放 dt、NaN 时整批回滚。

## 运行方式

在仓库根目录执行：

```bash
python -m pytest tests -q
```

初始状态必须恰好为 `15 failed`，且没有 collection、import 或 syntax error。所有失败都来自 `physics/contacts.py` 里 8 个尚未实现的桩函数。

## 字段表

| 字段 | 取值 |
| --- | --- |
| 所属领域 | 游戏、图形与多媒体 |
| 子方向 | 图形渲染 |
| 任务类型 | Bugfix |
| 语言/框架 | Python 3, pytest |
| 环境可复现等级 | 无外部依赖 |

## 约束

- `tests/test_contacts.py` 不可修改。
- 实现只使用 Python 标准库（`math`、`dataclasses`）和 pytest。
- 离线运行，无随机，不调用 `sleep`；固定时间步 `dt = 1/60`。
- 刚体使用内联 dict：`mass`、`inv_mass`、`pos`、`vel`、`radius`、`friction`、`restitution`。
- 位置比较容差 `1e-6`，速度比较容差 `1e-6`；异常类型为 `PhysicsError`。
