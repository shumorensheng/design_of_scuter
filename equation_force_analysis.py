from equation_motion_analysis import *
import numpy as np
import os
# 压强列表
p_list = np.array(
    [
        [0, 140],
        [-1, 140],
        [-1, 72.5],
        [-1, 50.5],
        [-1, 32.5],
        [-1, 15],
        [-1, 7.5],
        [-1, 2.5],
        [-1, 2.5],
        [-1, 1],
        [1.5, 1],
        [5, 1],
        [8, 1],
        [15, 1],
        [43.5, 1],
    ]
)


def solve_cordiration(tangential_angel_of_OA, tangential_angel_of_AB, e, a, b):
    """
    计算A,B点的坐标系下的坐标
    """
    A = np.float32(
        np.array(
            list(
                zip(
                    a * np.cos(tangential_angel_of_OA),
                    a * np.sin(tangential_angel_of_OA),
                )
            )
        )
    )
    B = np.float32(
        np.array(
            list(
                zip(
                    e * np.ones_like(tangential_angel_of_AB),
                    
                    a* np.sin(tangential_angel_of_OA) + b * np.sin(tangential_angel_of_AB)
                    ,
                )
            )
        )
    )
    return A, B


def solve_force(
    ac_x,
    ac_y,
    a_B,
    alpha_list,
    tangential_angel_of_OA,
    tangential_angel_of_AB,
    e,
    a,
    b,
    d,
    i,
    m1,
    m2,
    m3,
    pc,
    g=9.81,
):
    A, B = solve_cordiration(tangential_angel_of_OA, tangential_angel_of_AB, e, a, b)
    # print(tangential_angel_of_AB)
    
    # print(B)

    c_x = np.float32(A[:, 0].T + b * i * np.cos(tangential_angel_of_AB))
    c_y = np.float32(A[:, 1].T + b * i * np.sin(tangential_angel_of_AB))
    C = np.float32(np.vstack((c_x, c_y)).T)  # C点坐标
    # print(C)
    jc = m2 * pc * b**2/g#转动惯量

    #p假设向上为正
    area = np.float32( np.pi * (d / 2) ** 2)
    p1 = np.float32( p_list * 2 * area)  # 气体压力,压强取向下为正
    p2 = np.float32(m3 * a_B /  g).reshape(
        -1, 1
    )  # 惯性力 加速度是向下为正，为了校准方向，所以对重力加速度再取负
    
    p_force = -p1 + p2 - m3  # 作用在B滑块上的压力
    # 连杆的惯性力与惯性力矩
    #  这里的加速度再计算的时候是按照x正方向,y负方向计算的，只需要按照惯性力取负即可
    
    c_force_x = np.float32( -1*ac_x * m2 / g)
    c_force_y = np.float32( ac_y * m2 / g)
    # print(np.sqrt(ac_x**2 + ac_y**2))
    c_force = np.float32(np.sqrt(c_force_x**2 + c_force_y**2))
    c_torque = np.float32(-1 * jc * alpha_list)#c_torque是顺时针为正方向
    
    # 根据对B点的力矩计算R12t,力矩以逆时针方向为正
    r_bc = np.float32(C - B)
    
    M_gravity = -1*m2 * r_bc[:, 0]
    M_acceleration = c_force_y * r_bc[:, 0] - c_force_x * r_bc[:, 1]
    
    R12t = np.float32((c_torque - M_gravity - M_acceleration) / b)
    # print(M_gravity)
    # print(M_acceleration)
    R12t_x = np.float32(R12t * np.sin(tangential_angel_of_AB))#r12t_x是x轴正方向
    R12t_y = np.float32(R12t * np.cos(tangential_angel_of_AB))#r12t_y方向向下

    # 求解R12n
    R12n = np.float32(
        (p_force - R12t_y.reshape(-1, 1) + c_force_y.reshape(-1, 1) - m2)
        / np.sin(tangential_angel_of_AB).reshape(-1, 1)
    )

    # 求解R03
    R03 = np.float32(
        R12n * np.cos(tangential_angel_of_AB).reshape(-1, 1)
        - R12t_x.reshape(-1, 1)
        - c_force_x.reshape(-1, 1)
    )#R12n * np.cos(tangential_angel_of_AB)方向是x轴负方向，

    R12t_x = np.hstack((R12t_x, R12t_x))
    R12t_y = np.hstack((R12t_y, R12t_y))
    tangential_angel_of_AB = np.hstack((tangential_angel_of_AB, tangential_angel_of_AB))
    p_force = np.hstack((p_force[:, 0], p_force[:, 1]))
    R12n = np.hstack((R12n[:, 0], R12n[:, 1]))
    # 求解R12
    R12t=np.hstack((R12t,R12t))
    R12 = np.sqrt(R12t**2 + R12n**2)

    # 求解Mb
    R21_x = np.float32(R12n * np.cos(tangential_angel_of_AB)-R12t_x)#假设正方向为x正方向，为什么直接使用R12_x呢，因为R12的正负是为他们的正方向服务的
    R21_y = np.float32(R12t_y + R12n * np.sin(tangential_angel_of_AB))
    A_x = np.float32(np.hstack((A[:, 0], A[:, 0])))
    A_y = np.float32(np.hstack((A[:, 1], A[:, 1])))
    

    # print(A_x.shape, A_y.shape, R12_x.shape, R12_y.shape,tangential_angel_of_AB.shape)
    Mb = -1*np.float32(A_x * R21_y - A_y * R21_x)#逆时针为正，加负号转换为顺时针，同曲柄转向一致

    R03 = np.hstack((R03[:, 0], R03[:, 1]))
    R23 = np.float32(np.sqrt(R03**2 + p_force**2))
    c_force = np.hstack((c_force, c_force))
    c_torque = np.hstack((c_torque, c_torque))
    # R12t=np.hstack((R12t,R12t))

    return p_force, R12n, R12t, c_force, c_torque, R12, R03, Mb


def force(e, a, b, w, d, i, m1, m2, m3, pc):
    (
        sigma,
        shift,
        veltime,
        acctime,
        x,
        x_vals,
        shift_vals,
        veltime_vals,
        acctime_vals,
    ) = nihao(e, a, b, w,enable_singel_graph=False)
    (   angel_list,
        alpha_list,
        tangential_angel_of_OA,
        tangential_angel_of_AB,
        ac_x,
        ac_y,
        a_B,
    ) = output_result_v_and_a(sigma, e, a, b, w, i, print=False, output=True)
    
    os.makedirs('result', exist_ok=True)
    p_force, R12n, R12t, c_force, c_torque, R12, R03, Mb = solve_force(
        ac_x,
        ac_y,
        a_B,
        alpha_list,
        tangential_angel_of_OA,
        tangential_angel_of_AB,
        e,
        a,
        b,
        d,
        i,
        m1,
        m2,
        m3,
        pc,
        g=9.81,
    )
    angel_list=np.float32(np.hstack((angel_list,angel_list)))
    force_list =np.array(list(zip(angel_list, p_force, R12n, R12t, c_force, c_torque, R12, R03, -1*R12,Mb)))

    # print(len(angel_list),len(p_force),len(R12n),len(R12t),len(c_force),len(c_torque),len(R12),len(R03),len(Mb))
    with open('result/力分析.txt','w',encoding='utf-8') as f:
        f.write("第4行右边的垂直的角度、第8行是最低点的角度、第12行左边的垂直角度\n\
            和压力是向上为正方向，r12t是逆时针为正方向，ro3向右为正方向，、\n\
            ro1只有大小，Mi以及Mb是顺时针为正方向，\n")
        f.write(
                "{:<10}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n".format(
                    "angel", "和压力p(N)", "R12n(N)", "R12t(N)", 
                    "PI2(N)", "MI2(N·m)", "R12(N)", "R03(N)", "R01(N)", "Mb(N·m)"
                )
            )
        for i in range(len(force_list)):
            f.write(
                "{:<10.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}\n".format(
                    float(force_list[i][0]),
                    float(force_list[i][1]),
                    float(force_list[i][2]),
                    float(force_list[i][3]),
                    float(force_list[i][4]), 
                    float(force_list[i][5]), 
                    float(force_list[i][6]),
                    float(force_list[i][7]), 
                    float(force_list[i][8]),
                    float(force_list[i][9]),
                )
            )


if __name__ == "__main__":
    e = 0.060 # 曲柄偏心距单位m
    K = 1.05  # 行程速比系数
    H = 0.270  # 滑块行程单位m
    I = 0.36  # ab杆质心到a端距离系数
    D = 22  # 滑块直径 ，单位cm
    M1 = 170  # 曲柄质量单位N
    M2 = 135  # 连杆质量
    M3 = 210  # 滑块质量
    pc = 0.165  # 转动半径系数
    a, b = solve_longth_of_stick(K, H, e)  # 曲柄长度, 连杆长度
    print("曲柄长度: ", a)
    print("连杆长度: ", b)
    n = 610  # 转速rpm
    w = n * 2 * sp.pi / 60  # 转速转化为弧度/s

    force(e, a, b, w,D, I, M1, M2, M3,pc)
