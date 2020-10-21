import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import numpy as np
import math
from scipy.optimize import Bounds
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FuncFormatter, MultipleLocator
from datetime import date
#Data frames
df = pd.read_csv("time_series_covid19_confirmed_global.csv")
df1 = pd.read_csv("time_series_covid19_recovered_global.csv")
df2 = pd.read_csv("time_series_covid19_deaths_global.csv")

#Datos Guatemala --------------------------------------------------------
guatei = df.loc[df['Country/Region'] == 'Guatemala']
guatei = guatei.iloc[0,56:117]
guater = df1.loc[df1['Country/Region'] == 'Guatemala']
guater = guater.iloc[0,56:117]
guated = df2.loc[df2['Country/Region'] == 'Guatemala']
guated = guated.iloc[0,56:117]
guates = 17000000 - guater - guatei - guated
#Datos Italia ------------------------------------------------------------
italiai = df.loc[df['Country/Region'] == 'Italy']
italiai = italiai.iloc[0,15:117]
italiar = df1.loc[df1['Country/Region'] == 'Italy']
italiar = italiar.iloc[0,15:117]
italiad = df2.loc[df2['Country/Region'] == 'Italy']
italiad =italiad.iloc[0,15:117]
italias = 60.36e6- italiar - italiai - italiad
#Datos China --------------------------------------------------------------
chinai = df.loc[df['Province/State']=='Hubei']
chinai = chinai.iloc[0,4:117]
chinar = df1.loc[df1['Province/State']=='Hubei']
chinar = chinar.iloc[0,4:117]
chinad = df2.loc[df2['Province/State']=='Hubei']
chinad = chinad.iloc[0,4:117]
chinas = 11e6 - chinar - chinai-chinad

caso = "China" 

if(caso == "Guatemala"):
    N =  17000000
    ndatos = 61
elif(caso == "Italia"):
    N = 60.36e6
    ndatos = 102
elif(caso == "China"):
    N = 11e6
    ndatos = 113
#Functions -----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
def sird(t,y,beta,gamma,sigma):
        S,I,R,D = y
        return  [-beta * S * I/N, beta * S * I/N  - gamma * I - sigma* I, gamma * I, sigma * I]
#control--------------------------------------------------------------------------------------------------------------------------------
def beta_t(t,t0,b0,alpha):
    if(t <= t0):
        return b0
    elif(t > t0):
        return  b0*(1 - alpha*(t - t0)/t)

def sird_control(t,y,beta,gamma,sigma):
        S,I,R,D = y
        bt = beta_t(t,50,beta,0.5)
        return  [-bt * S * I/N, bt * S * I/N  - gamma * I - sigma* I, gamma * I, sigma * I]
    
def f(w):
    beta,gamma,sigma = w 
    I0, R0, D0  = 1,0,0
    S0 = N - I0 - R0 -D0 
    y0 = [S0, I0, R0,D0]
    t = np.linspace(0,ndatos,ndatos)
    solu= solve_ivp(sird, [0, ndatos], y0, method ='RK45' , t_eval=np.arange(0, ndatos, 1),dense_output = True, vectorized=False, args = (beta,gamma,sigma))
    S,I,R,D = solu.sol(t)
    error = 0
    if(caso == "Guatemala"):
        error1 = (np.mean((guates - S)**2)) 
        error2 = (np.mean((guatei - I)**2))
        error3 = (np.mean((guater - R)**2)) 
        error4 = (np.mean((guated - D)**2)) 
    elif(caso == "Italia"):
        error1 = (np.mean((italias - S)**2)) 
        error2 = (np.mean((italiai - I)**2))
        error3 = (np.mean((italiar - R)**2)) 
        error4 = (np.mean((italiad - D)**2)) 
 
    elif(caso == "China"):
        error1 = (np.mean((chinas - S)**2)) 
        error2 = (np.mean((chinai - I)**2))
        error3 = (np.mean((chinar - R)**2)) 
        error4 = (np.mean((chinad - D)**2)) 
    return error1  + error2  + error3  + error4

#Optimal parameters 
parametros = minimize(f,[0.0001, 0.0001, 0.0001],method='L-BFGS-B',bounds=[(0.00000001, 0.3), (0.00000001, 0.2),(0.00000001, 0.2)])
beta,gamma,sigma = parametros.x
print(parametros.x)

# initial conditions---------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
if(caso == "Guatemala"):
    I0, R0,D0  = guatei[0],guater[0], guated[0] 
    S0 = N - I0 - R0 - D0 
    y0 = [S0, I0, R0,D0]
elif(caso == "Italia"):
    I0, R0,D0  = italiai[0],italiar[0], italiad[0]
    S0 = N - I0 - R0 -D0 
    y0 = [S0, I0, R0, D0]
    print(y0)
elif(caso == "China"):
    I0, R0,D0 =chinai[0] , chinar[0],chinad[0]
    S0 = N - I0 - R0  -D0
    y0 = [S0, I0, R0,D0]
    print(y0)

#solution with optimized parameters------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
ndatos =300 
t = np.linspace(0,ndatos,ndatos)

#sol = odeint(sird,y0,t,args = (beta,gamma))
solu = solve_ivp(sird, [0, ndatos], y0, method ='RK45' , t_eval=np.arange(0, ndatos, 1),dense_output = True, vectorized=False, args = (beta,gamma,sigma))
S,I,R,D = solu.sol(t)

soluc = solve_ivp(sird_control, [0, ndatos], y0, method ='RK45' , t_eval=np.arange(0, ndatos, 1),dense_output = True, vectorized=False, args = (beta,gamma,sigma))
S_c,I_c,R_c,D_c= soluc.sol(t)



# Plot the data on four separate curves for S(t), I(t) and R(t)--------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
fig,a = plt.subplots(2,2)

#guates.plot(kind ='line',ax =a[0,0],label='data',logy = False)
a[0][0].plot(t,S ,label="susceptibles")
a[0][0].plot(t,S_c,label="susceptibles control")
#a[0][0].set_yscale("log")
a[0][0].set_title('Susceptibles')
a[0][0].legend(loc = 'upper right')
a[0][0].set_xlabel('tiempo (días)')
a[0][0].set_ylabel('Numero de susceptibles')

#guatei.plot(kind ='line',ax =a[1,0],label='data',logy = False)
a[1][0].plot(t,I ,label="infected")
a[1][0].plot(t,I_c,label="infected control")
#a[1][0].set_yscale("log")
a[1][0].set_title('Infectados')
a[1][0].legend(loc = 'upper right')
a[1][0].set_xlabel('tiempo (días)')
a[1][0].set_ylabel('Numero de infectados')

#guater.plot(kind ='line',ax =a[0,1],label='data',logy = False)
a[0][1].plot(t,R ,label="Recuperados")
a[0][1].plot(t,R_c,label="Recuperados control")
#a[0][1].set_yscale("log")
a[0][1].set_title('Recuperados')
a[0][1].legend(loc = 'upper right')
a[0][1].set_xlabel('tiempo (días)')
a[0][1].set_ylabel('Numero de recuperados')

#guated.plot(kind ='line',ax =a[1,1],label='data',logy = False)
a[1][1].plot(t,D ,label="Muertos")
a[1][1].plot(t,D_c,label="Muertos control")
#a[1][1].set_yscale("log")
a[1][1].set_title('Muertos')
a[1][1].legend(loc = 'upper right')
a[1][1].set_xlabel('tiempo (días)')
a[1][1].set_ylabel('Numero de muertos ')
plt.show()

def data_gen(t=0):
    cnt = 0
    while cnt < 10000:
        cnt += 1
        t += 1
        yield t, S_c[t],I_c[t],R_c[t],D_c[t] #np.sin(2*np.pi*t) * np.exp(-t/10.)


def init():
    ax.set_ylim(-1, N)
    ax.set_xlim(0, 150)
    ax.set_title(caso)
    ax.set_xlabel('tiempo (días)')
    ax.set_ylabel('Cantidad de personas')
    
    del xdata[:]
    del ydata[:]
    del y2data[:]
    del y3data[:]
    del y4data[:]
    line.set_data(xdata, ydata)
    line1.set_data(xdata, y2data)
    line2.set_data(xdata, y3data)
    line3.set_data(xdata, y4data)

    return line,

fig, ax = plt.subplots()
fig.set_size_inches(20,20)
line, = ax.plot([], [], lw=2)
line1, = ax.plot([],[], lw = 2)
line2, = ax.plot([],[], lw = 2)
line3, = ax.plot([],[], lw = 2)
annotation_s = ax.annotate('Susceptibles (0)', xy=(date(2020, 2, 22), 0),xytext=(date(2020, 2, 22),0))
annotation_i = ax.annotate('Infectados (0)', xy=(date(2020, 2, 22), 0),xytext=(date(2020, 2, 22),0))
annotation_r = ax.annotate('Recuperados (0)', xy=(date(2020, 2, 22), 0),xytext=(date(2020, 2, 22),0))
annotation_d = ax.annotate('Muertos (0)', xy=(date(2020, 2, 22), 0),xytext=(date(2020, 2, 22),0))
ax.grid()
xdata, ydata ,y2data,y3data,y4data = [], [], [],[],[]


def run(data):
    # update the data
    t, y,y2,y3,y4= data
    xdata.append(t)
    ydata.append(y)
    y2data.append(y2)
    y3data.append(y3)
    y4data.append(y4)
    xmin, xmax = ax.get_xlim()

    if t >= xmax:
        ax.set_xlim(xmin, 2*xmax)
        ax.figure.canvas.draw()
    line.set_data(xdata, ydata)
    line1.set_data(xdata,y2data)
    line2.set_data(xdata,y3data)
    line3.set_data(xdata,y4data)
    
    annotation_s.set_position((t,y))
    annotation_s.xy = (t,y)
    annotation_s.set_text('Susceptibles (%d)' % y)

    annotation_i.set_position((t,y2))
    annotation_i.xy = (t,y2)
    annotation_i.set_text('Infectados (%d)' % y2)
    
    annotation_r.set_position((t,y3))
    annotation_r.xy = (t,y3)
    annotation_r.set_text('Recuperados (%d)' % y3)

    annotation_d.set_position((t,y4))
    annotation_d.xy = (t,y4)
    annotation_d.set_text('Muertos (%d)' % y3)

    return line,

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,
                              repeat=False, init_func=init)
plt.show()
