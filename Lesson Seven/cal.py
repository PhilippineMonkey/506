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
total_time = 2000  # 总模拟时间 (秒)

# 环境参数
T_env = 180  # 环境温度 (°C)
h = 25       # 对流换热系数 (W/m²·K)


#configuration
#半径范围
config_R_list =[]
config_R_list.append(0.06)
config_R_list.append(0.08)
config_R_list.append(0.12)
config_R_list.append(0.14)
config_R_list.append(0.16)
config_R_list.append(0.18)
config_R_list.extend(np.arange(0.20,1.0,0.05))
#时间变化的R值
selected_R = 0.05




def solve_heat_equation(R):
    dr = R / (Nr - 1)
    r = np.linspace(0, R, Nr)
    T = np.ones(Nr) * 20 
    snapshots = {}
    surface_temps = []  
    time_points = []   
    
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
        surface_temps.append(T[-1])  
        time_points.append(t)
        
        if any(abs(t - target) < dt/2 for target in [300, 600, 1200]):
            snapshots[t] = T.copy()

    return r, snapshots, np.array(surface_temps), np.array(time_points) 


def solve_heat_time(R, target_time):
    dr = R / (Nr - 1)
    r = np.linspace(0, R, Nr)
    T = np.ones(Nr) * 20
    total_steps = int(target_time / dt) + 1
    
    for step in range(total_steps):
        t = step * dt
        T_new = T.copy()
        for i in range(1, Nr-1):
            term1 = (T[i+1] - 2*T[i] + T[i-1]) / dr**2
            term2 = (T[i+1] - T[i-1]) / (2*dr * r[i]) if r[i] != 0 else 0
            T_new[i] = T[i] + alpha * dt * (term1 + 2*term2)
        
        T_new[0] = T_new[1]
        surface_temp = (h * T_env * dr + k * T_new[-2]) / (h * dr + k)
        T_new[-1] = surface_temp
        
        T = T_new
        if abs(t - target_time) < dt/2:
            return surface_temp  
    
    return T[-1]  


r, snapshots, t, T_surface = solve_heat_equation(selected_R)

plt.figure(figsize=(10, 6))
plt.plot(t/60, T_surface, '-', lw=1, alpha=0.7) 
plt.scatter(t[::100]/60, T_surface[::100], s=20, c='r',  
            label=f'Δt={dt*100}s') 

plt.xlabel('Time (minutes)')
plt.ylabel('Surface Temperature (°C)')
plt.title(f'')
plt.axhline(180, c='gray', ls='--', label='Critical Temp (180°C)')
plt.legend()
plt.grid()
plt.savefig(f'./surface_temp_R{int(selected_R*100)}cm.png')
plt.close()


def plot_radius_impact():
    plt.figure(figsize=(8, 6))
    T_list = []
    for R in config_R_list:
        temp = solve_heat_time(R, 600)
        T_list.append(temp)
    plt.plot(config_R_list, T_list,label='Surface Temperature', marker='o', color='blue')
    plt.title('Temperature vs Radius at t=600s')
    plt.xlabel('Radius (m)')
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.grid()
    plt.savefig('./radius_impact.png')
    plt.close()

if __name__ == "__main__":
    plot_radius_impact()

