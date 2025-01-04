import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号无法显示的问题
def plot_graphs(x_vals, shift_vals, veltime_vals, acctime_vals):
    # 创建一个 1 行 3 列的子图
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # 绘制位移
    axs[0].plot(x_vals, shift_vals, label="shift")
    axs[0].set_title("shift")
    axs[0].set_xlabel("angel")
    axs[0].set_ylabel("s(m)")
    axs[0].legend()

    # 绘制速度
    axs[1].plot(x_vals, veltime_vals, label="velocity", color="orange")
    axs[1].set_title("velocity")
    axs[1].set_xlabel("angel")
    axs[1].set_ylabel("v(m/s)")
    axs[1].legend()

    # 绘制加速度
    axs[2].plot(x_vals, acctime_vals, label="acceleration", color="green")
    axs[2].set_title("acceleration")
    axs[2].set_xlabel("angel")
    axs[2].set_ylabel("a(m/s²)")
    axs[2].legend()

    # 自动调整子图参数和保存总图
    plt.tight_layout()
    plt.savefig("result/运动方程解析.png")  # 保存整个图
    plt.close()
    # 显示图像
    # #plt.show(block=False)  # 不阻塞后续代码执行

    # # 暂停2秒
    # plt.pause(2)

    # # 关闭图像
    # plt.close()


 
def save_singel_graph(s, v, a, x):
    s=sp.lambdify(x,s,"numpy")
    v=sp.lambdify(x,v,"numpy")
    a=sp.lambdify(x,a,"numpy")
    x_vals = np.linspace(0, 360, 1000)
    s_range = s(x_vals)
    v_range = v(x_vals)
    a_range = a(x_vals)
    plt.figure(1)
    plt.plot(x_vals, s_range,color="red")
    plt.title("位移")
    plt.xlabel("角度（°）")
    plt.ylabel("s(m)")
    plt.savefig("result/位移.png")
    plt.close()
    
    plt.figure(2)
    plt.plot(x_vals, v_range, color="orange")
    plt.title("速度")
    plt.xlabel("角度（°）")
    plt.ylabel("v(m/s)")
    plt.savefig("result/速度.png")
    plt.close()
    
    plt.figure(3)
    plt.plot(x_vals, a_range, color="green")
    plt.title("加速度")
    plt.xlabel("角度（°）")
    plt.ylabel("a(m/s²)")
    plt.savefig("result/加速度.png")    
    plt.close()

def find_extra_angle(sigma, a1, b, e):
    extra_angle = []
    angle_for_s_max = (sigma + sp.pi - sp.acos(e / (b - a1))) * 180 / sp.pi
    vertical_angel1 = (
        (sigma - sp.acos(e / sp.sqrt((a1**2 + b**2))) + sp.atan(b / a1)) * 180 / sp.pi
    )
    vertical_angel2 = (
        (-sigma + 2 * sp.pi + sp.acos(e /sp.sqrt (a1**2 + b**2)) - sp.atan(b / a1))
        * 180
        / sp.pi
    )
    extra_angle.append(angle_for_s_max), extra_angle.append(
        vertical_angel1
    ), extra_angle.append(vertical_angel2)
    extra_angle = [i.evalf() for i in extra_angle]
    return extra_angle


def output_result(sigma, a1, b, e, s, v, a, x, filename="result/辅助绘图.txt"):

    extra_angle = find_extra_angle(sigma, a1, b, e)

    angle_list = [i for i in range(0, 360, 30)]
    angle_list.extend(extra_angle)
    angle_list=np.sort(np.float32(angle_list))
    

    s_list = [s.subs(x, angles).evalf() for angles in angle_list]
    # s_list.extend([s.subs(x, angels).evalf() for angels in extra_angle])

    v_list = [v.subs(x, angles).evalf() for angles in angle_list]
    # v_list.extend([v.subs(x, angels).evalf() for angels in extra_angle])

    a_list = [a.subs(x, angles).evalf() for angles in angle_list]
    # a_list.extend([a.subs(x, angles).evalf() for angles in extra_angle])
    result = list(zip(angle_list, s_list, v_list, a_list))

    with open(filename, "w", encoding="utf-8") as f:
        f.write("第4行右边的垂直的角度、第8行是最低点的角度、第12行左边的垂直角度\n")
        f.write("{:<10}{:<15}{:<15}{:<15}\n".format("angel", "s(m)", "v(m/s)", "a(m/s²)"))  # 表头对齐

        for i in range(len(result)):
            f.write(
                "{:<10.2f}{:<15.4f}{:<15.4f}{:<15.4f}\n".format(
                    float(result[i][0]),  # 确保是float
                    float(result[i][1]),  # 确保是float
                    float(result[i][2]),  # 确保是float
                    float(result[i][3]),  # 确保是float
                )
            )



def tangential_angel_of_OA_AB(e, a, b, sigma, angel_list: np.ndarray):
    """
    这些切向量夹角都是与y轴夹角
    """
    # 计算垂直于 AB 方向的角度

    angel_list = np.radians(angel_list)  # 将角度转换为弧度
    tangential_angel_of_OA = sigma - angel_list
    tangential_angel_of_OA = np.float32(tangential_angel_of_OA)
    tangential_angel_of_AB = np.arccos((e - a * np.cos(tangential_angel_of_OA)) / b)
    tangential_angel_of_AB = np.float32(tangential_angel_of_AB)
    return tangential_angel_of_OA, tangential_angel_of_AB


def solve_velocity(tangential_angel_of_OA, tangential_angel_of_AB, v, b,i):
    orientation_of_B = tangential_angel_of_AB - tangential_angel_of_OA
    VBA_list = np.sin(tangential_angel_of_OA) * v / np.sin(tangential_angel_of_AB)
    VB_list = np.sin(orientation_of_B) * v / np.sin(tangential_angel_of_AB)
    W2_list = VBA_list / b
    W2_list = np.float32(W2_list)
    v = np.float32(v)
    VC2_list = np.sqrt(
        v**2
        + (W2_list *b*i) ** 2
        - 2*v *b*i* W2_list * np.cos(tangential_angel_of_AB - tangential_angel_of_OA)
    )

    return VB_list, VBA_list, W2_list, VC2_list  # 返回计算结果


def solve_acceleration(tangential_angel_of_OA, tangential_angel_of_AB, w_list, a, b, w,i):
    acceleration_of_A = np.float32(a * w**2)
    an_BA = np.float32(b * w_list**2)
    at_BA = np.float32(
        (
            acceleration_of_A * np.cos(tangential_angel_of_OA)
            + an_BA * np.cos(tangential_angel_of_AB)
        )
        / np.sin(tangential_angel_of_AB)
    )
    alpha_list = np.float32(at_BA / b)
    acceleration_of_B = np.float32(
        at_BA * np.cos(tangential_angel_of_AB)
        + acceleration_of_A * np.sin(tangential_angel_of_OA)
        + an_BA * np.sin(tangential_angel_of_AB)
    )
    a_BA = np.float32(np.sqrt(at_BA**2 + an_BA**2))
    
    
    an_CA=np.float32(i *b* w_list**2)
    at_CA=np.float32(i*b*alpha_list)
    ac_x=np.float32(at_CA*np.sin(tangential_angel_of_AB)-an_CA*np.cos(tangential_angel_of_AB)-acceleration_of_A*np.cos(tangential_angel_of_OA))
    ac_y=np.float32(at_CA*np.cos(tangential_angel_of_AB)+an_CA*np.sin(tangential_angel_of_AB)+acceleration_of_A*np.sin(tangential_angel_of_OA))
    ac_B = np.float32(np.sqrt(ac_x**2 + ac_y**2))
    
    return an_BA, at_BA, a_BA, alpha_list,ac_x,ac_y, ac_B, acceleration_of_B  # 返回计算结果


def output_result_v_and_a(
    sigma, e, a, b, w ,i, filename1="result/表一_速度.txt", filename2="result/表二_加速度.txt",print:bool=True,output:bool=False
):
    v = w * a
    extra_angle = find_extra_angle(sigma, a, b, e)

    angle_list = [i for i in range(0, 360, 30)]
    angle_list.extend(extra_angle)
    angle_list =np.sort(np.float32(angle_list)) 

    tangential_angel_of_OA, tangential_angel_of_AB = tangential_angel_of_OA_AB(
        e, a, b, sigma, angle_list
    )

    VB_list, VBA_list, W2_list, VC2_list = solve_velocity(
        tangential_angel_of_OA, tangential_angel_of_AB, v, b,i
    )
    an_BA, at_BA, a_BA, alpha_list,ac_x,ac_y, ac_B, acceleration_of_B = solve_acceleration(
        tangential_angel_of_OA, tangential_angel_of_AB, W2_list, a, b, w,i
    )

    velocity_list = list(zip(angle_list, VB_list, VBA_list, VC2_list, W2_list))
    acceleration_list = list(
        zip(angle_list, an_BA, at_BA, a_BA, alpha_list, ac_B, acceleration_of_B)
    )
    if print:
        with open(filename1, "w", encoding="utf-8") as f:
            # 写入文件头
            f.write("第4行右边的垂直的角度、第8行是最低点的角度、第12行左边的垂直角度\n")

            f.write(
                "{:<10}{:<15}{:<15}{:<15}{:<15}\n".format(
                    "angel", "VB(m/s)", "VBA(m/s)", "VC2(m/s)", "W2(rad/s)"
                )
            )

            # 写入数据，并将 SymPy 对象转换为 float
            for i in range(len(velocity_list)):
                f.write(
                    "{:<10.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}\n".format(
                        float(velocity_list[i][0]),
                        float(velocity_list[i][1]),
                        float(velocity_list[i][2]),
                        float(velocity_list[i][3]),
                        float(velocity_list[i][4]),
                    )
                )


    # 文件写入部分
        with open(filename2, "w",encoding="utf-8") as f:
            # 写入表头
            f.write("第4行右边的垂直的角度、第8行是最低点的角度、第12行左边的垂直角度\n")

            
            f.write(
                "{:<10}{:<15}{:<15}{:<15}{:<15}{:<15}{:<25}\n".format(
                    "angel", "an_BA(m/s²)", "at_BA(m/s²)", "a_BA(m/s²)", 
                    "alpha(rad/s²)", "ac_B(m/s²)", "aB(m/s²)"
                )
            )

            # 写入数据
            for i in range(len(acceleration_list)):
                f.write(
                    "{:<10.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<25.2f}\n".format(
                        acceleration_list[i][0],
                        acceleration_list[i][1],
                        acceleration_list[i][2],
                        acceleration_list[i][3],
                        acceleration_list[i][4],
                        acceleration_list[i][5],
                        acceleration_list[i][6],
                    )
                )

    if output:
        return angle_list,alpha_list,tangential_angel_of_OA, tangential_angel_of_AB,ac_x,ac_y,acceleration_of_B

def solve_longth_of_stick(k, h, e):
    """
    计算连杆长度
    参数：
    k:行程速比系数
    h:滑块行程
    返回：曲柄以及连杆长度
    """
    angle = 180 * (k - 1) / (k + 1)  # 极位偏角

    angle = np.radians(angle)
    R = h * 0.5 / np.sin(angle)
    CD = np.sqrt(R**2 - (h / 2) ** 2)
    OE = np.sqrt(R**2 - (CD - e) ** 2)
    OF = OE + h / 2

    zeta = np.arcsin(OE / R) - angle

    CONST1 = np.sqrt(OF**2 + e**2)
    CONST2 = 2 * R * np.sin(zeta / 2)

    a = (CONST1 - CONST2) / 2
    b = (CONST1 + CONST2) / 2

    return a, b


def nihao(e, a, b, w,enable_singel_graph=True):
    """
    运 动方程解析函数:
    参数：
    e: 曲柄偏心距
    a: 曲柄长度
    b: 连杆长度
    w: 转速
    """
    # 定义符号变量和表达式
    h = sp.sqrt(((a + b) ** 2 - e**2))  # 最高终止点到转动中心距离
    sigma = sp.acos(e / (a + b))  # 连杆与水平面的夹角

    x = sp.symbols("x")  # x 为转过的角度
    shift = (
        h
        - sp.sqrt(b**2 - (a * sp.cos(sigma - x * sp.pi / 180) - e) ** 2)
        - a * sp.sin(sigma - x * sp.pi / 180)
    )  # 运动方程表达式
    veltime = sp.diff(w * 180 / sp.pi * shift, x)  # 运动速度表达式
    acctime = sp.diff(w * 180 / sp.pi * veltime, x)  # 运动加速度表达式

    # veltime = w*veltime
    # acctime = w**2*acctime
    if enable_singel_graph:
        save_singel_graph(shift, veltime, acctime, x)

    # 生成 x 的数值范围（0 到 360度）
    x_vals = np.linspace(0, 360, 360)

    # 替换 x 的值计算运动方程、速度和加速度
    shift_vals = [shift.subs(x, angle).evalf() for angle in x_vals]
    veltime_vals = [veltime.subs(x, angle).evalf() for angle in x_vals]
    acctime_vals = [acctime.subs(x, angle).evalf() for angle in x_vals]
    return sigma,shift,veltime,acctime,x,x_vals, shift_vals, veltime_vals, acctime_vals
    
def motion(e, a, b, w,i):
    os.makedirs('result', exist_ok=True)

    sigma,shift,veltime,acctime,x,x_vals, shift_vals, veltime_vals, acctime_vals\
         = nihao(e, a, b, w)
    # 调用绘图函数
    #     # 调用绘图函数
    output_result(sigma, a, b, e, shift, veltime, acctime, x)
    output_result_v_and_a(sigma, e, a, b, w ,i)
    plot_graphs(x_vals, shift_vals, veltime_vals, acctime_vals)
    

if __name__ == "__main__":
    e = 0.07  # 曲柄偏心距
    K = 1.08  # 行程速比系数
    H = 0.32  # 滑块行程
    I= 0.4    # ab杆质心到a端距离系数
    a, b = solve_longth_of_stick(K, H, e)  # 曲柄长度, 连杆长度
    print("曲柄长度: ", a)
    print("连杆长度: ", b)
    n = 590  # 转速rpm
    v = n * 2 * sp.pi / 60  # 转速转化为弧度/s

    motion(e, a, b, v,I)
    
