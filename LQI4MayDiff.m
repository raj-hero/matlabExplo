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

xa=zeros(25,1);
u=[0
    0];
delp=[0.01
    0
    0
    0];


c=zeros(1,25);
% c(1,1)=1;
c(1,3)=1;

t=0;
it=1;

y(1)=[0];
time(it)=t;
cd=c;

Q=100*eye(25);
Q(26,:)=0;
Q(:,26)=0;
Q(26,26)=100;
v2=[50 50];
R=diag(v2);
N=zeros(26,2);

% fun=@(y) -y;
% q=integral(fun,0,1);
diff(q,y)=-y;
xak=[xa;q];
ck=[c 0];
bk=[b;0 0];
tauk=[tau;0 0 0 0];
ve1=zeros(25,1);
a=[a ve1];
c=[c 0];
ak=[a;-c];

[gain,sy,er]=lqr(ak,bk,Q,R,N);
[ad,bd]=c2d(ak,bk,delt);
[ad1,taud]=c2d(ak,tauk,delt);
while t<40
    it=it+1;
    x=ad*xak+bd*u+taud*delp;
    y(it)=ck*x;
    xak=x;
    t=t+delt;
    time(it)=t;
    u=-gain*xak;
end
plot(time,y);
