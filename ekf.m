load('data.mat');

%inicjalizacja
e0 = 1;
e1 = 0;
e2 = 0;
e3 = 0;
pb = 0; %bias p �yroskopu
qb = 0; %bias q �yroskopu
rb = 0; %bias r zyroskopu

%macierze zale�ne od poziomu zaszumienia czujnik�w
%mo�na nimi regulowa�
 P = zeros(7,7);
 Q = diag([[1 1 1 1] * 0.00005, [1 1 1] * 0.000001] .^ 2);
 R = diag([[1 1 1] * 0.045, [1 1 1] * 0.015]);
        
 x = [e0 e1 e2 e3 pb qb rb]'; %inicjalizacja wektora stanu

for i=2:length(time)
 
%czas pr�bkowania
dt = time(i) - time(i-1);
    
%%%%%%%%%%%
%PREDYKCJA%
%%%%%%%%%%%

%dane z �yroskopu
p = gx(i)*pi/180; q = gy(i)*pi/180; r = gz(i)*pi/180;

%wektor wej��
u = [p q r pb qb rb]';

%macierz przej�cia
F = 0.5*[-e1 -e2 -e3 e1 e2 e3;
         e0 -e3  e2 -e0 e3 -e2;
         e3  e0 -e1 -e3 -e0 e1;
        -e2  e1  e0 e2 -e1 -e0;
         0   0   0   0  0   0;
         0   0   0   0  0   0;
         0   0   0   0  0   0];
 
%estymacja wektora stanu
x = x + dt*F*u;
%aktualizacja kwaternionion�w:
e0 = x(1);
e1 = x(2);
e2 = x(3);
e3 = x(4);

%normalizacja kwaternion�w 
norm = sqrt(e0^2 + e1^2 + e2^2 + e3^2);
e0 = e0 / norm;
e1 = e1 / norm;
e2 = e2 / norm;
e3 = e3 / norm;
x(1) = e0;
x(2) = e1;
x(3) = e2;
x(4) = e3;

%Jakobian A - pochodne cz�stkowe dF/du 
A = 0.5*[0 -(p-pb) -(q-qb) -(r-rb) e1 e2 e3;
        (p-pb) 0 (r-rb) -(q-qb)  -e0 e3 -e2;
        (q-qb) -(r-rb) 0 (p-pb)  -e3 -e0 e1;
        (r-rb) (q-qb) -(p-pb) 0   e2 -e1 -e0;
        0        0       0    0   0  0   0;
        0        0       0    0   0  0   0;
        0        0       0    0   0  0   0];

%estymacja macierzy kowariancji
P = P + dt*(A*P+P*A'+Q);

%modelowanie magnetometru
dm = 0; %k�t deklinacji magnetycznej - przyjmuj� 0 bo nie wiem sk�d pochodz� dane
msin=sind(dm); 
mcos=cosd(dm);

m =[msin*(2*e0*e3+2*e1*e2)-mcos*(2*e2*e2+2*e3*e3-1)...
    -mcos*(2*e0*e3-2*e1*e2)-msin*(2*e1*e1+2*e3*e3-1)...
    mcos*(2*e0*e2+2*e1*e3)-msin*(2*e0*e1-2*e2*e3)];

%modelowanie akcelerometru - jedynie wektora grawitacji
a = [2*(e1*e3-e0*e2) 2*(e0*e1+e2*e3) 1-2*(e1^2+e2^2)];

%macierz modeli
z = [a m]';

%normalizacja modeli
norm = sqrt(z(1)^2 + z(2)^2 + z(3)^2);
z(1) = z(1) / norm;
z(2) = z(2) / norm;
z(3) = z(3) / norm;
norm = sqrt(z(4)^2 + z(5)^2 + z(6)^2);
z(4) = z(4) / norm;
z(5) = z(5) / norm;
z(6) = z(6) / norm;

%pomiar z akcelerometru i magnetometru
y = [ax(i) ay(i) az(i) mx(i) my(i) mz(i)]';

%normalizacja pomiaru z akcelerometru i magnetometru
norm = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
y(1) = y(1) / norm;
y(2) = y(2) / norm;
y(3) = y(3) / norm;
norm = sqrt(y(4)^2 + y(5)^2 + y(6)^2);
y(4) = y(4) / norm;
y(5) = y(5) / norm;
y(6) = y(6) / norm;

%%%%%%%%%%       
%KOREKCJA%
%%%%%%%%%%


%Jakobian H - pochodne cz�stkowe dz/dx
H = 2  * [e2 -e3 e0 -e1 0 0 0;
         -e1 -e0 -e3 -e2 0 0 0;
          0  2*e1 2*e2 0 0 0 0;
          e3*msin, e2*msin, e1*msin-2*e2*mcos, e0*msin-2*e3*mcos, 0,0,0;
          -e3*mcos, e2*mcos-2*e1*msin,  e1*mcos, -e0*mcos-2*e3*msin, 0,0,0;
           e2*mcos-e1*msin, e3*mcos-e0*msin, e0*mcos+e3*msin, e1*mcos+e2*msin, 0,0,0];


%wzmocnienie K
K = P*H'/(H*P*H' + R);
%macierz kowariancji
P = (eye(7,7) - K*H)*P;
%wyj�cie filtru
x = x + K*(y - z);

%aktualizacja wektora stanu 
e0 = x(1);
e1 = x(2);
e2 = x(3);
e3 = x(4);
pb = x(5);
qb = x(6);
rb = x(7);

%normalizacja kwaternion�w 
norm = sqrt(e0^2 + e1^2 + e2^2 + e3^2);
e0 = e0 / norm;
e1 = e1 / norm;
e2 = e2 / norm;
e3 = e3 / norm;
x(1) = e0;
x(2) = e1;
x(3) = e2;
x(4) = e3;

%przej�cie na k�ty Eulera
phi(i) = atan2((2*(e0*e1+e3*e2)),1-2*(e1^2+e2^2))*180/pi;
theta(i) = asin(2*(e0*e2-e3*e1))*180/pi;
psi(i) = atan2((2*(e0*e3+e1*e2)),1-2*(e2^2+e3^2))*180/pi;

%k�ty korekcji - do ewentualnych wykres�w
phi_acc(i) = atan2(ay(i), az(i));
theta_acc(i) = atan2(-ax(i), ay(i) * sin(phi_acc(i)) + az(i) * cos(phi_acc(i)));

phi_acc(i) = phi_acc(i) * 180/pi;
theta_acc(i) = theta_acc(i) * 180/pi;
psi_mag(i) = atan2(-my(i),mx(i))* 180/pi;

end

%wykresy
f1 = figure;
figure(f1);
subplot(3,1,1);
plot(time,theta_komp);
grid;
hold on;
plot(time,theta);
title('theta')
legend('filtr komplementarny','filtr Kalmana')

subplot(3,1,2)
plot(time,phi_komp);
grid;
hold on;
plot(time,phi);
title('phi')
legend('filtr komplementarny','filtr Kalmana')

subplot(3,1,3)
plot(time,psi_komp);
grid;
hold on;
plot(time,psi);
title('psi')
legend('filtr komplementarny','filtr Kalmana')

f2 = figure;
figure(f2);
subplot(3,1,1);
plot(time,gx);
grid;
hold on;
plot(time,gy);
plot(time,gz);
title('dane z �yroskop�w');
legend('p','q','r');

subplot(3,1,2)
plot(time,ax);
grid;
hold on;
plot(time,ay);
plot(time,az);
title('dane z akcelerometru')
legend('ax','ay','az')

subplot(3,1,3)
plot(time,ax);
grid;
hold on;
plot(time,ay);
plot(time,az);
title('dane z magnetometru')
legend('mx','my','mz')



