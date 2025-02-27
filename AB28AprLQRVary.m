% clear all;
% clc;
delt=10^-3;
Tg1=0.08;
Tt1=0.30;
Tr1=10.0;
Kr1=0.3;
TRH1=41.600;
TR1=5.000;
TGH1=0.513;
TW1=1.000;
X1=0.60;
Y1=1.00;
b1=0.05;
c1=1.0;
TF1=0.23;
TCR1=0.01;
TCD1=0.20;
TWD1=5.0;
T12=0.0433;
R1=2.4;
bt1=0.425;
al1=-1.0;
KP1=120;
TP1=20;
Kt1=0.60;
Kh1=0.20;
Kg1=0.10;
KW1=0.10;

a=zeros(25,25);

a(1,1)=-1/TP1;
a(1,2)=-KP1/TP1;
a(1,4)=KP1/TP1;
a(1,7)=KP1/TP1;
a(1,10)=KP1/TP1;

a(1,24)=KP1/TP1;
a(2,1)=2*pi*T12;
a(2,3)=-2*pi*T12;
a(3,2)=-(al1*KP1)/TP1;
a(3,3)=-1/TP1;
a(3,14)=KP1/TP1;
a(3,17)=KP1/TP1;
a(3,20)=KP1/TP1;
a(3,25)=KP1/TP1;

a(4,4)=-1/Tr1;
a(4,5)=(Kt1/Tr1)-(Kr1*Kt1/Tt1);
a(4,6)=Kr1*Kt1/Tt1;
a(5,5)=-1/Tt1;
a(5,6)=1/Tt1;
a(6,1)=-1/(R1*Tg1);
a(6,6)=-1/Tg1;

a(7,1)=(2*Kh1*TR1)/(TGH1*R1*TRH1);
a(7,7)=-2/TW1;
a(7,8)=(2*Kh1/TW1)+(2*Kh1/TGH1);
a(7,9)=((2*Kh1*TR1)/(TGH1*TRH1))-(2*Kh1/TGH1);

a(8,1)=-TR1/(TGH1*R1*TRH1);
a(8,8)=-1/TGH1;
a(8,9)=(1/TGH1)-(TR1/(TGH1*TRH1));
a(9,1)=-1/(R1*TRH1);
a(9,9)=-1/TRH1;

a(10,10)=-1/TCD1;
a(10,11)=Kg1/TCD1;
a(10,12)=-(Kg1*TCR1)/(TF1*TCD1);
a(11,11)=-1/TF1;
a(11,12)=1/TF1+(TCR1/(TF1*TF1));

a(12,1)=-X1/(b1*R1*Y1);
a(12,12)=-c1/b1;
a(12,13)=1/b1;

a(13,1)=(X1/(R1*Y1*Y1))-(1/(R1*Y1));
a(13,13)=-1/Y1;

a(14,14)=-1/Tr1;
a(14,15)=(Kt1/Tr1)-(Kr1*Kt1/Tt1);
a(14,16)=Kr1*Kt1/Tt1;
a(15,15)=-1/Tt1;
a(15,16)=1/Tt1;
a(16,3)=-1/(R1*Tg1);
a(16,16)=-1/Tg1;

a(17,3)=(2*Kh1*TR1)/(TGH1*R1*TRH1);
a(17,17)=-2/TW1;
a(17,18)=(2*Kh1/TW1)+(2*Kh1/TGH1);
a(17,19)=((2*Kh1*TR1)/(TGH1*TRH1))-(2*Kh1/TGH1);

a(18,3)=-TR1/(TGH1*R1*TRH1);
a(18,18)=-1/TGH1;
a(18,19)=(1/TGH1)-(TR1/(TGH1*TRH1));
a(19,3)=-1/(R1*TRH1);
a(19,19)=-1/TRH1;

a(20,20)=-1/TCD1;
a(20,21)=Kg1/TCD1;
a(20,22)=-(Kg1*TCR1)/(TF1*TCD1);
a(21,21)=-1/TF1;
a(21,22)=1/TF1+(TCR1/(TF1*TF1));

a(22,3)=-X1/(b1*R1*Y1);
a(22,22)=-c1/b1;
a(22,23)=1/b1;

a(23,3)=(X1/(R1*Y1*Y1))-(1/(R1*Y1));
a(23,23)=-1/Y1;

a(24,24)=-1/TWD1;
a(25,25)=-1/TWD1;

b=zeros(25,2);

b(6,1)=1/Tg1;
b(7,1)=(-2*Kh1*TR1)/(TGH1*TRH1);
b(8,1)=TR1/(TGH1*TRH1);
b(9,1)=1/TRH1;

b(12,1)=X1/(b1*Y1);
b(13,1)=(1/Y1)-(X1/(Y1*Y1));

b(16,2)=1/Tg1;
b(17,2)=(-2*Kh1*TR1)/(TGH1*TRH1);
b(18,2)=TR1/(TGH1*TRH1);
b(19,2)=1/TRH1;

b(22,2)=X1/(b1*Y1);
b(23,2)=(1/Y1)-(X1/(Y1*Y1));

tau=zeros(25,4);

tau(1,1)=-KP1/TP1;
tau(3,2)=-KP1/TP1;
tau(24,3)=1/TWD1;
tau(25,4)=1/TWD1;

det(a);
eig(a);



[ad,bd]=c2d(a,b,delt);
[ad1,taud]=c2d(a,tau,delt);

xa=zeros(25,1);
u=[0
    0];
delp=[0.01
    0
    0
    0];


c=zeros(1,25);
c(1,1)=1;
c(3,3)=1;

t=0;
it=1;

y(1)=[0];
time(it)=t;
cd=c;

% Q MATRIX START
v1=100;
Q=zeros(25,25);
for ii=1:1:25
    for jj=1:1:25
        if ii==jj
            Q(ii,jj)=v1;
        end
    end
end
% Q MATRIX END

% R MATRIX START
v2=50;
R=zeros(2,2);
R(1,1)=v2;
R(2,2)=v2;
% R MATRIX END
N=zeros(25,2);
[gain,sy,er]=lqr(a,b,Q,R,N);


while t<40
    it=it+1;
    x=ad*xa+bd*u+taud*delp;
    y(it)=cd*x;
    xa=x;
    t=t+delt;
    time(it)=t;
    
    for ii=1:1:25
        for jj=1:1:25
            if ii==jj
                Q(ii,jj)=Q(ii,jj)+10;
            end
        end
    end
    
    R(1,1)=R(1,1)+5;
    R(2,2)=R(2,2)+5;
    u=-gain*xa;
end
plot(time,y);

lqr_resp=[transpose(time) transpose(y)];
plot(transpose(time),transpose(y));
filename='AB28AprLQRVaryxl';
xlswrite(filename,lqr_resp);