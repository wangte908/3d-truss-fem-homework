import numpy as np

def truss3d_element_stiffness(x1, x2, E, A):
    x1 = np.array(x1, dtype=float)
    x2 = np.array(x2, dtype=float)
    dx = x2 - x1
    L = np.linalg.norm(dx)
    # 退化单元判断
    if L < 1e-12:
        raise ValueError("错误：两个节点重合，单元退化，无法计算。")
    cx, cy, cz = dx / L
    cc = np.outer(np.array([cx, cy, cz]), np.array([cx, cy, cz]))
    Ke = (E * A / L) * np.block([[cc, -cc], [-cc, cc]])
    return L, (cx, cy, cz), Ke

def truss3d_element_stress(x1, x2, E, A, de):
    L, c, _ = truss3d_element_stiffness(x1, x2, E, A)
    cx, cy, cz = c
    de = np.array(de, dtype=float)
    B = (1.0 / L) * np.array([-cx, -cy, -cz, cx, cy, cz])
    epsilon = np.dot(B, de)
    sigma = E * epsilon
    N = sigma * A
    return epsilon, sigma, N

if __name__ == "__main__":
    # 算例1：沿X轴一维杆单元
    print("==========算例1：沿X轴三维杆单元==========")
    x1_1 = [0, 0, 0]
    x2_1 = [2, 0, 0]
    E1 = 200e9
    A1 = 1e-4
    de1 = [0, 0, 0, 1e-3, 0, 0]
    L1, c1, Ke1 = truss3d_element_stiffness(x1_1, x2_1, E1, A1)
    eps1, sig1, N1 = truss3d_element_stress(x1_1, x2_1, E1, A1, de1)
    print(f"杆长L = {L1:.6f} m")
    print(f"方向余弦 cx,cy,cz = {c1}")
    print("单元刚度矩阵Ke：")
    print(Ke1)
    print(f"应变ε={eps1:.6e}, 应力σ={sig1/1e6:.2f} MPa, 轴力N={N1:.2f} N\n")

    # 算例2：空间任意斜杆单元
    print("==========算例2：空间任意方向杆单元==========")
    x1_2 = [0, 0, 0]
    x2_2 = [1, 2, 2]
    E2 = 210e9
    A2 = 2e-4
    de2 = [0, 0, 0, 1e-3, 2e-3, 2e-3]
    L2, c2, Ke2 = truss3d_element_stiffness(x1_2, x2_2, E2, A2)
    eps2, sig2, N2 = truss3d_element_stress(x1_2, x2_2, E2, A2, de2)
    print(f"杆长L = {L2:.6f} m")
    print(f"方向余弦 cx,cy,cz = {c2}")
    print(f"应变ε={eps2:.6e}, 应力σ={sig2/1e6:.2f} MPa, 轴力N={N2:.2f} N")

    # 刚度矩阵对称性验证
    sym_diff = np.max(np.abs(Ke2 - Ke2.T))
    print(f"\n刚度矩阵对称最大误差：{sym_diff:.2e}")
    # 特征值验证半正定、奇异性
    eig_vals = np.linalg.eigvals(Ke2)
    print(f"刚度矩阵最小特征值：{np.min(eig_vals):.2e}")
    # 刚体位移检验
    de_rigid = np.array([1,1,1,1,1,1])
    eps_rigid, _, _ = truss3d_element_stress(x1_2, x2_2, E2, A2, de_rigid)
    print(f"刚体平移下轴向应变：{eps_rigid:.2e}\n")

    # 刚度矩阵物理意义验证
    print("==========刚度矩阵物理意义验证==========")
    de_unit = np.array([1,0,0,0,0,0])
    Fe = Ke2 @ de_unit
    print("仅自由度1单位位移对应的节点力列阵：")
    print(Fe)
    print("该结果与Ke第一列完全相等。")