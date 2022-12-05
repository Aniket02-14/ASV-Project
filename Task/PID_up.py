import  numpy as np
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

def PID(b,c,a,e,i,del_t,xf,psi,r,x,y,t):
        j = len(x) - 1
        while x[-1] <= xf:
            j += 1
            k1 = b*e+c*i-a*r[j-1]
            k2 = b*(e-r[j-1]*del_t/2)+c*(i-e*del_t/2)-a*(r[j-1]+del_t*k1/2)
            k3 = b*(e-r[j-1]*del_t/2)+c*(i-e*del_t/2)-a*(r[j-1]+del_t*k2/2)
            k4 = b*(e-r[j-1]*del_t)+c*(i-e*del_t)-a*(r[j-1]+del_t*k3)
            k_o = (k1+2*k2+2*k3+k4)/6
            R = r[j-1] + k_o
            r = np.append(r,R)
            e = e - r[j]*del_t
            PSI = psi[j-1] + r[j]*del_t
            psi = np.append(psi,PSI)
            #print(PSI)
            i += e*del_t
            X = x[j-1] + v*np.sin(psi[j-1])*del_t
            x = np.append(x,X)
            Y = y[j-1] + v*np.cos(psi[j-1])*del_t
            y = np.append(y,Y)
            T = t[j-1] + del_t
            t = np.append(t,T)
        return x,y,psi,r,t

xsi = 500
ysi = 100
ini_x = np.array([xsi,ysi])
xf = np.array(xsi)
yf = np.array(ysi)

xai = 0.0
yai = 0.0
ini_a = np.array([xai,yai])

xaf = 1250
yaf = 751
fin_a = np.array([xaf,yaf])

cte = np.abs(np.linalg.norm(np.cross(fin_a-ini_a, ini_a-ini_x)))/np.linalg.norm(fin_a-ini_a)
m = yaf/xaf
m_p = -1/m
theta = np.pi/2-np.arctan(m_p)
x_inter = np.array([xsi + cte*np.cos(theta), ysi + cte*np.sin(theta)])
delta = np.linalg.norm(x_inter - fin_a)

alpha = np.pi/2 - np.arctan(m)

T = 107
k = 0.2
kp = 1
kd = 10
ki = 0.01
v = 10
del_t = 0.1
psi_d = alpha - np.arctan(cte/delta)
t = np.zeros(1)
ka = (1+k*kd)/T
kb = k*kp/T
kc = k*ki/T
ew = psi_d
iw = 0
d = psi_d*k*ki


psi = np.zeros(1)
r = np.zeros(1)
x = np.zeros(1)
x[0] = xsi
y = np.zeros(1)
y[0] = ysi
[x,y,psi,r,t] = PID(kb,kc,ka,ew,iw,del_t,1250,psi,r,x,y,t)
#print(x)
# x_main1 = x
# y_main1 = m*x_main1

xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

xbf = 3445
ybf = 946


ew = np.pi/2 - np.arctan((ybf-y[-1])/(xbf-x[-1])) - psi[-1]
iw = 0
xa = np.zeros(1)
ya = np.zeros(1)
xa[0] = x[-1]
ya[0] = y[-1]
ra = np.zeros(1)
psia = np.zeros(1)                                                                                                                   
psia[0] = psi[-1]
ta = np.zeros(1)
ta[0] = t[-1]
# print(psi)
#print(psi[-1])
#print(psia)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,3445,psia,ra,xa,ya,ta)
# x_main2 = xa
# y_main2 = ((ybf-y[-1])/(xbf-x[-1]))*x_main2
x = np.concatenate((x,xa))
y = np.concatenate((y,ya))
psi = np.concatenate((psi,psia))
r = np.concatenate((r,ra))
t = np.concatenate((t,ta))

xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xcf = 5847
ycf = -48


ew = 3*np.pi/2 - (np.pi + np.arctan((ycf-y[-1])/(xcf-x[-1]))) - psi[-1]
iw = 0
xa = np.zeros(1)
ya = np.zeros(1)
xa[0] = x[-1]
ya[0] = y[-1]
ra = np.zeros(1)
psia = np.zeros(1)
psia[0] = psi[-1]
ta = np.zeros(1)
ta[0] = t[-1]
# print(psi)
#print(psi[-1])
#print(psia)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xcf,psia,ra,xa,ya,ta)
# x_main3 = xa
# y_main3 = ((ycf-y[-1])/(xcf-x[-1]))*x_main3
# print(x)
x = np.concatenate((x,xa))
y = np.concatenate((y,ya))
psi = np.concatenate((psi,psia))
r = np.concatenate((r,ra))
t = np.concatenate((t,ta))


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])
# x_main = np.concatenate((x_main1,x_main2,x_main3))
# y_main = np.concatenate((y_main1,y_main2,y_main3))
# fig = plt.figure() 
# axis = plt.axes(xlim =(0, 3200),
#                 ylim =(0, 1500)) 
  
# line, = axis.plot([], [], lw = 2)


# def init(): 
#     line.set_data([], []) 
#     return line, 
   
# # initializing empty values
# # for x and y co-ordinates
# xdata, ydata = [], [] 
   
# # animation function 
# def animate(i):
       
#     xdat = x[i]
#     ydat = y[i]
       
#     xdata.append(xdat) 
#     ydata.append(ydat) 
#     line.set_data(xdata, ydata) 
      
#     return line,
       
# anim = FuncAnimation(fig, animate, init_func = init, 
#                                frames = None, interval = 1, blit = True)

xi =[0,1250,3445,5847]
yi = [0,751,946,-48]

plt.figure(1)
plt.plot(x,y)
plt.plot(xi,yi)
plt.scatter(xf,yf)


plt.figure(2)
plt.plot(t,psi)

plt.figure(3)
plt.plot(t,r)
plt.show()