clear all;
clc;
kp=100;
tp=20;
tt=0.5;
tg=0.4;
R=3;
delp=0.01;
u=0;
A=[-1/tp  kp/tp  0
    0 -1/tt 1/tt
    -1/(R*tg) 0 -1/tg];
B=[0
    0
    1/tg];
Tau=[-kp/tp
    0
    0];
C=[1 0 0];
D=0;

delt=1*10^-3;

Ad=expm(A*delt);

I=eye(length(A));

Bd=inv(A)*(Ad-I)*B;
Cd=C;
Dd=D;
Taud=inv(A)*(Ad-I)*Tau;

X0=[0
    0
    0];

t=0;
it=1;
y(1)=0;
time(it)=t;

while t<20
    it=it+1;
    X=Ad*X0+Bd*u+Taud*delp;
    y(it)=Cd*X;
    X0=X;
    t=t+delt;
    time(it)=t;
end

plot(time,y)