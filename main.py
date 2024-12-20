from cam_analysis import cam
from  equation_force_analysis import force
from equation_motion_analysis import motion,solve_longth_of_stick
import sympy as sp
import matplotlib.pyplot as plt

# 设置字体，确保你的环境中有该字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号无法显示的问题

def main(e0, a, b,  I,  w, D, M1, M2, M3, pc, h, theta1, theta2, theta3, theta4, e1, r, r0, h2, e2, r2):
    
    motion(e0, a, b, w,I)
    force(e0, a, b, w,D, I, M1, M2, M3,pc)
    cam(h,theta1,theta2,theta3,theta4,e1,r,r0,name='凸轮1',filename1='result/凸轮1-s-δ曲线.png',filename2='result/凸轮1-轮廓线.png')
    cam(h2,theta1,theta2,theta3,theta4,e2,r2,r0,name='凸轮2',filename1='result/凸轮2-s-δ曲线.png',filename2='result/凸轮2-轮廓线.png')
    
    
    
if __name__ == '__main__':
    e0 = 0.07  # 曲柄偏心距单位m
    K = 1.08  # 行程速比系数
    H = 0.32  # 滑块行程单位m
    I = 0.38  # ab杆质心到a端距离系数
    D = 23  # 滑块直径 ，单位cm
    M1 = 190  # 曲柄质量单位N
    M2 = 140  # 连杆质量
    M3 = 230  # 滑块质量
    pc = 0.17  # 转动半径系数
    a, b = solve_longth_of_stick(K, H, e0)  # 曲柄长度, 连杆长度
    print("曲柄长度: ", a)
    print("连杆长度: ", b)
    n = 590  # 转速rpm
    w = n * 2 * sp.pi / 60  # 转速转化为弧度/s
    h=9 #行程 单位mm
    theta1=60 #升程角
    theta3=10 #远休角
    theta2=60 #回程角
    theta4=360-theta1-theta2-theta3 #近休角
    e1=3 #偏心距 单位mm
    r=55 #基圆半径 单位mm
    
    h2=10
    e2=0
    r2=60
    r0=10 #凸轮滚子半径 单位mm

    main(e0, a, b,  I,  w, D, M1, M2, M3, pc, h,\
        theta1, theta2, theta3, theta4, e1, r, r0, h2, e2, r2)
