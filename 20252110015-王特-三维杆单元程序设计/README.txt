三维桁架杆单元有限元程序
运行环境：Python3 + numpy
安装依赖：pip install numpy
运行方式：PyCharm打开truss3d.py，右键Run执行
功能说明：
1. truss3d_element_stiffness：输入两点坐标、弹性模量E、截面积A，输出杆长、方向余弦、6×6全局刚度矩阵；两点重合会抛出报错，防止退化单元。
2. truss3d_element_stress：输入坐标、材料、节点位移向量，输出杆件轴向应变、应力、轴力。
内置两组验证算例，自动输出结果，包含刚度矩阵对称性、刚体位移、特征值、刚度物理意义校验。