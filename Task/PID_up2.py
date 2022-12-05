import  numpy as np
import matplotlib.pyplot as plt

# defining pid function
def PID(b,c,a,e,i,del_t,xf,yf,psi,r,x,y,t):
        j = len(x) - 1
        if yf > 0 or (yf == 0 and xf > 0):
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
        else:
            while x[-1] >= xf:
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


# defining parameters
def point(x_cord, y_cord, t, psi):
    x_i = np.zeros(1)
    y_i = np.zeros(1)
    t_i = np.zeros(1)
    psi_i = np.zeros(1)
    r_i = np.zeros(1)

    x_i[0] = x_cord[-1]
    y_i[0] = y_cord[-1]
    t_i[0] = t[-1]
    psi_i[0] = psi[-1]

    return x_i,y_i,t_i,psi_i,r_i


# updating parameters
def update(x_cord, y_cord, t, psi, r, x_i, y_i, t_i, psi_i, r_i):
    x_cord = np.concatenate((x_cord,x_i))
    y_cord = np.concatenate((y_cord,y_i))
    psi = np.concatenate((psi,psi_i))
    r = np.concatenate((r,r_i))
    t = np.concatenate((t,t_i))

    return x_cord,y_cord,psi,r,t


# slope finding
def slope(xi,yi,xf,yf,psi):
    m = (yf-yi)/(xf-xi)
    theta = np.arctan(m)
    if m > 0:
        psid = np.pi/2 - theta
    else:
        psid = -theta + np.pi/2

    if (yi >= 0 and yf >= 0):
        ew = psid-psi
    else:
        ew = np.pi + psid-psi
    iw = 0
    return ew, iw

xsi = -2000
ysi = 0
ini_x = np.array([xsi,ysi])
xf = np.array(xsi)
yf = np.array(ysi)

xai = -2000.0
yai = 0.0
ini_a = np.array([xai,yai])

xaf = -1000*np.sqrt(3)
yaf = 1000
fin_a = np.array([xaf,yaf])

# cte = np.abs(np.linalg.norm(np.cross(fin_a-ini_a, ini_a-ini_x)))/np.linalg.norm(fin_a-ini_a)
# m = (yaf-yai)/(xaf-xai)
# m_p = -1/m
# theta = np.pi/2-np.arctan(m_p)
# x_inter = np.array([xsi + cte*np.cos(theta), ysi + cte*np.sin(theta)])
# delta = np.linalg.norm(x_inter - fin_a)

# alpha = np.pi/2 - np.arctan(m)

T = 107
k = 0.2
kp = 1
kd = 10
ki = 0.02
v = 10
del_t = 0.1
# # psi_d = alpha - np.arctan(cte/delta)
# psi_d = slope(xai,yai,xaf,yaf)
t = np.zeros(1)
ka = (1+k*kd)/T
kb = k*kp/T
kc = k*ki/T


psi = np.zeros(1)
r = np.zeros(1)
x = np.zeros(1)
x[0] = xsi
y = np.zeros(1)
y[0] = ysi
[ew,iw]= slope(xai,yai,xaf,yaf,psi[0])
[x,y,psi,r,t] = PID(kb,kc,ka,ew,iw,del_t,xaf,yaf,psi,r,x,y,t)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

xbf = -1000
ybf = 1000*np.sqrt(3)

[ew,iw]= slope(x[-1],y[-1],xbf,ybf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xbf,ybf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)
xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xcf = 0
ycf = 2000


[ew,iw]= slope(x[-1],y[-1],xcf,ycf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xcf,ycf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

# xi =[0,1250,3445,5847]
# yi = [0,751,946,-48]

xdf = 1000
ydf = 1000*np.sqrt(3)

[ew,iw]= slope(x[-1],y[-1],xdf,ydf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xdf,ydf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])



xef = 1000*np.sqrt(3)
yef = 1000

[ew,iw]= slope(x[-1],y[-1],xef,yef,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xef,yef,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xff = 2000
yff = 0

[ew,iw]= slope(x[-1],y[-1],xff,yff,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xff,yff,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xgf = 1000*np.sqrt(3)
ygf = -1000

[ew,iw]= slope(x[-1],y[-1],xgf,ygf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xgf,ygf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

xhf = 1000
yhf = -1000*np.sqrt(3)

[ew,iw]= slope(x[-1],y[-1],xhf,yhf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xhf,yhf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

xif = 0
yif = -2000

[ew,iw]= slope(x[-1],y[-1],xif,yif,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xif,yif,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)


xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xjf = -1000
yjf = -1000*np.sqrt(3)

[ew,iw]= slope(x[-1],y[-1],xjf,yjf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xjf,yjf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)
xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])

xkf = -1000*np.sqrt(3)
ykf = -1000
[ew,iw]= slope(x[-1],y[-1],xkf,ykf,psi[-1])
[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xkf,ykf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)

xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])


xlf = -2000
ylf = 0

[ew,iw]= slope(x[-1],y[-1],xlf,ylf,psi[-1])

[xa,ya,ta,psia,ra] = point(x,y,t,psi)
[xa,ya,psia,ra,ta] = PID(kb,kc,ka,ew,iw,del_t,xlf,ylf,psia,ra,xa,ya,ta)
[x,y,psi,r,t] = update(x,y,t,psi,r,xa,ya,ta,psia,ra)
print(xa)

xf = np.append(xf,x[-1])
yf = np.append(yf,y[-1])



plt.figure(1)
plt.plot(x,y)
#plt.plot(xi,yi)
plt.scatter(xf,yf)


plt.figure(2)
plt.plot(t,psi)

plt.figure(3)
plt.plot(t,r)
plt.show()