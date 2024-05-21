clc
clear
close all

%Initialization

beta = 0.1; %how fast it transmits
D = 14; %average time spent sick
gamma = 1/D; %recovery rate
N = 5; %total number of population in millions
R_zero = beta/gamma; % how many people are infected by one person
R_zero

dt = 0.1; %time variable in days
time = 100; % roughly two years
t = 0:dt:time; %increment

Nt = length(t); % return the length of t/all time locations incremented by dt

%Store the results, in the following array
S = zeros(1,Nt);
I = zeros(1,Nt);
R = zeros(1,Nt);

I_initial = 0.0001; % initial population infected
I(1) = I_initial;
S(1) = N-I_initial; % initial suseptible population, given that R0=0

for i = 1:Nt-1
  %use the given formulae
  dS = -beta*I(i)*S(i); % change in S over change in time
  S(i+1) = S(i) + dS*dt; % S over time
  
  S_analytical(i+1) = S(i) * exp((-1.4*R(i+1))/N);
  
  dI = (beta*I(i)*S(i)) - gamma*I(i); % change in I over in time 
  I(i+1) = I(i)+dI*dt; % I over time
 
  dR = gamma*I(i); %change in R over change in time
  R(i+1) = R(i) + dR*dt; % R over time
 
end
%calculate the last value for the suseptible population
%this is because the for loop goes up until Nt-1



% Newton Raphsons Method:
% x_n+1 = x_n - (f(x_n)/f'(x_n))
%Starting guess and associated function values
x_old = 0.0001 ;
f_x_old = my_func(x_old) ;
f_deriv_x_old = my_func_deriv(x_old) ;


%How accurate do you want answer to be?
threshold = 1e-5 ;


%Keep going until f(x_old) is very small
while (abs(f_x_old) > threshold )
 
%Newton Raphson update
x_new = x_old - (f_x_old/f_deriv_x_old) ;

x_old = x_new ;
f_x_old = my_func(x_old) ;
% error = x_new - answer ;
f_deriv_x_old = my_func_deriv(x_old) ;
x_old
end

%plot the graph
hold on,
plot(t, x_old, t, S, t, I, t, R, t, S_analytical)
xlabel('Time in weeks')
ylabel('Population in millions')
legend('Rinf','Susceptible','Infected','Recovered', 'S_a')
title('Spread of COVID-19 (Ireland)')
grid on;
grid minor;
set(gca, 'FontSize', 26)