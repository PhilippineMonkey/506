import numpy as np
import matplotlib.pyplot as plt

# 材料参数
k = 0.5    # 热导率 (W/m·K)
rho = 980  # 密度 (kg/m³)
c = 3500   # 比热容 (J/kg·K)
alpha = k / (rho * c)  # 热扩散率


R_list = [0.05, 0.075, 0.1]  # 鸡块半径 (m)
Nr = 100     # 空间网格数
dt = 0.1    # 时间步长 (秒)
total_time = 1200  # 总模拟时间 (秒)

# 环境参数
T_env = 180  # 环境温度 (°C)
h = 25       # 对流换热系数 (W/m²·K)



def solve_heat_equation(R):
    dr = R / (Nr - 1)
    r = np.linspace(0, R, Nr)
    T = np.ones(Nr) * 20 
    snapshots = {}
    
    for step in range(int(total_time/dt) + 1):
        t = step * dt  
        T_new = T.copy()
        for i in range(1, Nr-1):
            term1 = (T[i+1] - 2*T[i] + T[i-1]) / dr**2
            term2 = (T[i+1] - T[i-1]) / (2*dr * r[i]) if r[i] != 0 else 0
            T_new[i] = T[i] + alpha * dt * (term1 + 2*term2)
        
        T_new[0] = T_new[1]  
        T_new[-1] = (h * T_env * dr + k * T_new[-2]) / (h * dr + k)
        T = T_new
        if t in [300, 600, 1200]:
            snapshots[t] = T.copy()

    return r, snapshots

for idx, R in enumerate(R_list):
    r, snapshots = solve_heat_equation(R)
    
    plt.figure(figsize=(8, 6))
    for time, data in snapshots.items():
        plt.plot(r, data, label=f'{time//60} min')
    
    plt.title(f'')
    plt.xlabel('Radius (m)')
    plt.grid(False)
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.savefig(f'./heat_distribution{idx+1}.png')  
    plt.close()
