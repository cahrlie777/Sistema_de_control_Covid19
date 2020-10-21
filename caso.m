%% Clear workspace
clc;
clear all; 
close all; 
%% Import Data
%%%% https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases
glb_conf = readtable('time_series_covid19_confirmed_global.csv'); 
glb_dead = readtable('time_series_covid19_deaths_global.csv');
glb_recv = readtable('time_series_covid19_recovered_global.csv');
%%% Hubei 55
%%% Guate 118
%%% Italia 132
%%% Japón 141
%%% USA 227
%%% España 203
%% Get coutry data
pais = 'guate'; 
corte=20;
switch pais
    case 'china'
        inf = glb_conf(64, corte:end);
        dead = glb_dead(64, corte:end);
        recov = glb_recv(55,corte:end); 
    case 'guate'
        inf = glb_conf(125, corte:end);
        dead = glb_dead(125, corte:end);
        recov = glb_recv(118,corte:end); 
    case 'italia'
        inf = glb_conf(139, corte:end);
        dead = glb_dead(139, corte:end);
        recov = glb_recv(133,corte:end);
    case 'japon'
        inf = glb_conf(141, corte:end);
        dead = glb_dead(141, corte:end);
        recov = glb_recv(135,corte:end);
    case 'usa'
        inf = glb_conf(227, corte:end);
        dead = glb_dead(227, corte:end);
        recov = glb_recv(227,corte:end);
    case 'spain'
        inf = glb_conf(203, corte:end);
        dead = glb_dead(203, corte:end);
        recov = glb_recv(201,corte:end);
        size(recov)
end
inf_double = str2double((table2array(inf)));
INF = [inf_double]'; 
recov_double = str2double((table2array(recov)));
RECOV = [recov_double]'; 
dead_double = str2double((table2array(dead)));
DEAD = [dead_double]'; 
switch pais 
    case 'china'
        POP = 58500000;
    case 'guate'
        POP = 17000000; 
    case 'italia'
        POP = 60360000; 
    case 'japon'
        POP = 126600000; 
    case 'usa'
        POP = 328200000; 
    case 'spain'
        POP = 46940000; 
end
%%
%% PARAMETROS
% lb = [0.1,0.001]; 
% ub = [0.5,0.15];
lb = [0,0,0];
ub = [10,10,10];
GammaR = 0.0101;
GammaD = 0.0101; 
Beta = 0.2027; 
GammaR = 0;
GammaD = 0; 
Beta = 0.01;
y0 = [POP, 1, 3, 1];
%%
%% PARAMETROS
%% Prepare Data
SUSC = POP-INF-DEAD-RECOV; 
DB_SIRD = [SUSC,INF,RECOV,DEAD];
time = size(DB_SIRD,1);
t_db = 1:1:time; 
%% Optimización 
x0=[Beta,GammaR,GammaD];
A = [];
b = [];
Aeq = []; 
beq = []; 
fun = @(x)func_sird(x(1),x(2),x(3),DB_SIRD, POP,time,y0);
RES = fmincon(fun,x0,A,b,Aeq,beq,lb,ub); 
%% Resolver ecuación diferencial utilizando las nuevas constantes
[t,y] = ode45(@(t,y) SIRD(POP, RES(1), RES(2), RES(3), y), t_db, y0);
graficar(1,1,3,1,t_db,DB_SIRD,2,"(I)(DB)","(I)(modelo)",t,y,time,RES);
graficar(1,1,3,2,t_db,DB_SIRD,3,"(R)(DB)","(R)(modelo)",t,y,time,RES);
graficar(1,1,3,3,t_db,DB_SIRD,4,"(D)(DB)","(D)(modelo)",t,y,time,RES);
function graficar(fignum,fila,columna,subplotx,vector_t,vector_db,columna_graf,strLdb,strLy1,t,y,time,RES)
figure(fignum);
subplot(fila,columna,subplotx);
db1 = plot(vector_t,vector_db(:,columna_graf)); 
hold on;
y1 = plot(t,y(:,columna_graf));
xlim([0 time]);
ylim([0 max([max(vector_db(:,columna_graf)) max(y(:,columna_graf))])]);
xlabel('Tiempo');
ylabel('Personas');
beta_f=num2str(RES(1,1));
gammaR_f=num2str(RES(1,2));
gammaD_f=num2str(RES(1,3));
title("SIRD:" + strLdb + " vs " + strLy1 + newline + "Beta=" + beta_f + newline + "gammaR=" + gammaR_f + newline + "gammaD=" + gammaD_f);
grid on;
grid minor;
legend([db1; y1], strLdb, strLy1);
end
function escalar = func_sird(beta,gamma_r, gamma_d,DB, pop,span,y0)
tspan = 1:1:span;
[t,y] = ode45(@(t,y) SIRD(pop, beta, gamma_r, gamma_d, y), tspan, y0);
escalar = 0; 
L = size(tspan,2);  
for m = 1:L
    escalar = escalar + norm(y(m,:) - DB(m,:))^2; 
end
end
function x_prima = SIRD(N, beta, gamma_r, gamma_d,x)
S = x(1);
I = x(2); 
dS = -(beta/N)*I.*S;
dI = (beta/N).*S.*I - (gamma_r+gamma_d).*I; 
dR = gamma_r*I;
dD = gamma_d.*I; 
x_prima = [dS,dI,dR, dD]'; 
end