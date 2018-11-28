clear; close; clc;
%% constants
W = 0.05; %width
H = 0.005; %height
L = 0.1; %length
k = 120; %conductivity
h = 10; %heat transfer coefficient
Rc = 0.15; %contact resistance between base & 1st element
Qb = 10; %base heat flux [W]
T_inf = 20; %ambient temperature
n = 1000; %number of fin elements (when n=70, Tb=121)
c = 1e-5; %covergence criteria (abs((Tnew-Told)/Tnew))
%% variables
dx = L/n;
As = 2*(H+W)*dx;
Ac = H*W;
R1 = dx/(k*Ac);
R2 = 1/(h*As);
R_tip = 1/(h*Ac);
%% calculation
A(n,n) = 0;
B(n,1) = 0;
A(1,1:2)=[(1/R1+1/R2) -1/R1];
A(n,(n-1):n)=[-1/R1 (1/R1+1/R2+1/R_tip)];
B([1 n],:)=[(Qb+T_inf/R2) T_inf*(1/R_tip+1/R2)];
for i = 2:(n-1)
    A(i,(i-1):(i+1))=[-1/R1 (2/R1+1/R2) -1/R1];
    B(i,:)=T_inf/R2;
end
T = A\B; %temperature distribution column
for i = 1:(n-1)
    if abs((T(i+1)-T(i))/T(i+1)) >= c %convergence criteria
        T(i) = T(i+1);
    end
end
Tb = Qb*(Rc+0.5*R1)+T(1); %base temperature
Q = (sum(T)-n*T_inf)/R2;
Q_max = h*(n*As+2*Ac)*(Tb-T_inf);
e = Q/Q_max; %fin efficiency
x(n) = 0; Q_conv(n) = 0;
for i = 1:n
    x(i) = dx*i;
    Q_conv(i) = (T(i)-T_inf)/R2;
end
%% plot
figure(1);
plot(x,T/Tb);
title('Temperature Plot n = 1000');
xlabel('Dimensionless position x*=x/L');
ylabel('Dimensionless temperature T(x*)/Tb');
figure(2);
plot(x,Q_conv/Qb);
title('Heat Transfer Plot');
xlabel('Dimensionless position x*=x/L');
ylabel('Dimensionless heat transfer Q_conv(x*)/Qb');