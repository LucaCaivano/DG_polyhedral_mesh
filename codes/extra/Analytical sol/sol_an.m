%soluzione analitica
clear all
syms x
syms y
% syms Young
% syms Poisson
% syms mu
% syms lambda

Young=10^7;
Poisson=0.48;
%Cohesion=450;
%Friction_angle=pi/9; 

mu          = Young/(2*(1+Poisson));                                               % 
lambda       = Young*Poisson/((1-2*Poisson)*(1+Poisson)); 
U=[10.0 + sin(pi*x)*exp(pi*y); 10.0 + sin(pi*x)*sin(pi*y)];
A=diff(U(1), x);
B=diff(U(2), y);
C=(diff(U(1), y)+diff(U(2), x))/2;
Eps=[A C; C B]; 
Sigma=2*mu*Eps+lambda*(A+B)*[1 0; 0 1];
Div_sigma=mydivergence(Sigma, [x, y]);
f=-Div_sigma;


%Neumann 
%           g3          
%       ---------
%    g4 -       -  g2
%       -       -
%       ---------
%           g1


n1=[0; -1];
n2=[1;  0];
n3=[0;  1];
n4=[-1; 0];

g1=Sigma*n1;
g2=Sigma*n2;
g3=Sigma*n3;
g4=Sigma*n4;


