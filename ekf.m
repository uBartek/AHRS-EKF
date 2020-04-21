load('data.mat');

%%%%%%%%%%%%%%%%
%INITIALIZATION%
%%%%%%%%%%%%%%%%

%quaternions
e0 = 1;
e1 = 0;
e2 = 0;
e3 = 0;

pb = 0; %bias p 
qb = 0; %bias q 
rb = 0; %bias r 

%covariance matrix
P = zeros(7,7);
%process noise matrix
Q = diag([[1 1 1 1] * 0.00005, [1 1 1] * 0.000001] .^ 2);
%measurement noise matrix
R = diag([[1 1 1] * 0.045, [1 1 1] * 0.015]);
        
%state space init
x = [e0 e1 e2 e3 pb qb rb]';

for i=2:length(time)
 
%sample time
dt = time(i) - time(i-1);
    
%%%%%%%%%%%%%%%%%
%PREDICTION STEP%
%%%%%%%%%%%%%%%%%

%read data from gyro
p = gx(i)*pi/180; q = gy(i)*pi/180; r = gz(i)*pi/180;

%input vector
u = [p q r pb qb rb]';

%transition matrix
F = 0.5*[-e1 -e2 -e3 e1 e2 e3;
         e0 -e3  e2 -e0 e3 -e2;
         e3  e0 -e1 -e3 -e0 e1;
        -e2  e1  e0 e2 -e1 -e0;
         0   0   0   0  0   0;
         0   0   0   0  0   0;
         0   0   0   0  0   0];
 
%state space estimate
x = x + dt*F*u;

%update quaternions value
e0 = x(1);
e1 = x(2);
e2 = x(3);
e3 = x(4);

%normalise quaternions
norm = sqrt(e0^2 + e1^2 + e2^2 + e3^2);
e0 = e0 / norm;
e1 = e1 / norm;
e2 = e2 / norm;
e3 = e3 / norm;
x(1) = e0;
x(2) = e1;
x(3) = e2;
x(4) = e3;

%Jacobian matrix A - partial derivatives dF/du 
A = 0.5*[0 -(p-pb) -(q-qb) -(r-rb) e1 e2 e3;
        (p-pb) 0 (r-rb) -(q-qb)  -e0 e3 -e2;
        (q-qb) -(r-rb) 0 (p-pb)  -e3 -e0 e1;
        (r-rb) (q-qb) -(p-pb) 0   e2 -e1 -e0;
        0        0       0    0   0  0   0;
        0        0       0    0   0  0   0;
        0        0       0    0   0  0   0];

%covariance matrix estimate
P = P + dt*(A*P+P*A' + Q);

%magnetometer model estimation
dm = 0; %magnetic declination angle
msin=sind(dm); 
mcos=cosd(dm);

m =[msin*(2*e0*e3+2*e1*e2)-mcos*(2*e2*e2+2*e3*e3-1)...
    -mcos*(2*e0*e3-2*e1*e2)-msin*(2*e1*e1+2*e3*e3-1)...
    mcos*(2*e0*e2+2*e1*e3)-msin*(2*e0*e1-2*e2*e3)];

%accelerometer model estimation
a = -[2*(e1*e3-e0*e2) 2*(e0*e1+e2*e3) 1-2*(e1^2+e2^2)];

%models matrix
z = [a m]';

%measure from acc and mag
y = [-ax(i) -ay(i) -az(i) mx(i) my(i) mz(i)]';

%normalise measurements
norm = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
y(1) = y(1) / norm;
y(2) = y(2) / norm;
y(3) = y(3) / norm;
norm = sqrt(y(4)^2 + y(5)^2 + y(6)^2);
y(4) = y(4) / norm;
y(5) = y(5) / norm;
y(6) = y(6) / norm;

%%%%%%%%%%%%%     
%UPDATE STEP%
%%%%%%%%%%%%%

%Jacobian matrix H - partial derivatives dy/dx 
H = 2  * [e2 -e3 e0 -e1 0 0 0;
         -e1 -e0 -e3 -e2 0 0 0;
          0  2*e1 2*e2 0 0 0 0;
          e3*msin, e2*msin, e1*msin-2*e2*mcos, e0*msin-2*e3*mcos, 0,0,0;
          -e3*mcos, e2*mcos-2*e1*msin,  e1*mcos, -e0*mcos-2*e3*msin, 0,0,0;
           e2*mcos-e1*msin, e3*mcos-e0*msin, e0*mcos+e3*msin, e1*mcos+e2*msin, 0,0,0];

%gain 
K = P*H'/(H*P*H' + R);
%covariance matrix
P = (eye(7,7) - K*H)*P;
%state space
x = x + K*(y - z);

%update quaternions and biases
e0 = x(1);
e1 = x(2);
e2 = x(3);
e3 = x(4);
pb = x(5);
qb = x(6);
rb = x(7);

%normalise quaternions
norm = sqrt(e0^2 + e1^2 + e2^2 + e3^2);
e0 = e0 / norm;
e1 = e1 / norm;
e2 = e2 / norm;
e3 = e3 / norm;
x(1) = e0;
x(2) = e1;
x(3) = e2;
x(4) = e3;

%Euler angles
phi(i) = atan2((2*(e0*e1+e3*e2)),1-2*(e1^2+e2^2))*180/pi;
theta(i) = asin(2*(e0*e2-e3*e1))*180/pi;
psi(i) = atan2((2*(e0*e3+e1*e2)),1-2*(e2^2+e3^2))*180/pi;

end

%plots
f1 = figure;
figure(f1);
subplot(3,1,1);
plot(time,theta_komp);
grid;
hold on;
plot(time,theta);
title('theta')
legend('complementary filter','Kalman filter')
ylabel('angle [deg]')

subplot(3,1,2)
plot(time,phi_komp);
grid;
hold on;
plot(time,phi);
title('phi')
legend('complementary filter','Kalman filter')
ylabel('angle [deg]')

subplot(3,1,3)
plot(time,psi_komp);
grid;
hold on;
plot(time,psi);
title('psi')
legend('complementary filter','Kalman filter')
ylabel('angle [deg]')

f2 = figure;
figure(f2);
subplot(3,1,1);
plot(time,gx);
grid;
hold on;
plot(time,gy);
plot(time,gz);
title('gyro');
legend('p','q','r');
ylabel('angular velocity [deg/s]')

subplot(3,1,2)
plot(time,ax);
grid;
hold on;
plot(time,ay);
plot(time,az);
title('accelerometer')
legend('ax','ay','az')
ylabel('acceleration [g]')

subplot(3,1,3)
plot(time,ax);
grid;
hold on;
plot(time,ay);
plot(time,az);
title('magnetometer')
legend('mx','my','mz')
xlabel('time [s]')
ylabel('flux [G]')



