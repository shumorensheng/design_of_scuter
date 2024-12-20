import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import os
# 设置字体，确保你的环境中有该字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号无法显示的问题

def cam(h,theta1,theta2,theta3,theta4,e,r,r0,name,filename1,filename2):
    os.makedirs('result',exist_ok=True)
    x=sp.symbols('x')
    s_pression=sp.Piecewise((h*(x/theta1-sp.sin(2*sp.pi*x/theta1)/(2*sp.pi)),(x>=0)&(x<theta1)),\
                            (h,(x>=theta1)&(x<theta1+theta3)),\
                            (h*(1-(x-theta1-theta3)/theta2+sp.sin(2*sp.pi*(x-theta1-theta3)/theta2)/(2*sp.pi)),(x>=theta1+theta3)&(x<theta1+theta3+theta2)),\
                            (0,(x>=theta1+theta3+theta2)&(x<theta1+theta3+theta2+theta4)))
    s_pression=sp.lambdify(x,s_pression,'numpy')
    x_range=np.linspace(0,theta1+theta3+theta2+theta4,1000)
    y_range=s_pression(x_range)
    plt.figure(1)
    plt.plot(x_range,y_range)
    plt.xlabel('δ(°)')
    plt.ylabel('s(mm)')
    plt.title('%ss-δ曲线'%name)
    plt.savefig(filename1)
    # # 显示图像
    # plt.show(block=False)  # 不阻塞后续代码执行

    # # 暂停2秒
    # plt.pause(2)

    # # 关闭图像
    plt.close()

    
    
    s0=np.sqrt(r**2-e**2)
    zeta=np.linspace(0,360,1000)
    x=(s0+s_pression(zeta))*np.sin(np.radians(zeta))+e*np.cos(np.radians(zeta))
    y=(s0+s_pression(zeta))*np.cos(np.radians(zeta))-e*np.sin(np.radians(zeta))
    sin_theta=y/np.sqrt(x**2+y**2)
    cos_theta=x/np.sqrt(x**2+y**2)
    x_2=x-r0*cos_theta
    y_2=y-r0*sin_theta
    
    
    plt.figure(2)#解析法设计凸轮轮廓
    plt.plot(x,y,'r')
    plt.plot(x_2,y_2,'b')
    plt.legend(['理论廓线','工作廓线'])
    plt.title('%s轮廓'%name)
    plt.axis('off')
    # 显示图像
    plt.savefig(filename2)
    # plt.show(block=False)  # 不阻塞后续代码执行

    # # 暂停2秒
    # plt.pause(2)

    # # 关闭图像
    plt.close()
    
    

if __name__ == '__main__':
    h=9 #行程 单位mm
    theta1=60 #升程角
    theta3=10 #远休角
    theta2=60 #回程角
    theta4=360-theta1-theta2-theta3 #近休角
    e=3 #偏心距 单位mm
    r=55 #基圆半径 单位mm
    
    h2=10
    e2=0
    r2=60
    r0=10 #凸轮滚子半径 单位mm
    cam(h,theta1,theta2,theta3,theta4,e,r,r0,name='凸轮1',filename1='result/凸轮1-s-δ曲线.png',filename2='result/凸轮1-轮廓线.png')
    cam(h2,theta1,theta2,theta3,theta4,e2,r2,r0,name='凸轮2',filename1='result/凸轮2-s-δ曲线.png',filename2='result/凸轮2-轮廓线.png')