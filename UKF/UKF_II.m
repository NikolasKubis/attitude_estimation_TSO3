% The quaternion Unscented Kalman Filter.


clear all;
% inertia tensor
inertia=diag([6*10^(0);7*10^(0);9*10^(0)]); 

dt=0.01;

simulation_time=100;

samples=simulation_time/dt;


% tracking directions are chosen to be constant. However, the filter works
% with time varying reference directions as well.
tracking_direction1=[0 0 1]';

tracking_direction2=[0 1 0]';

% establishing a deterministic model error
ll=1;
 for i=0:dt:simulation_time

model_error(:,ll)=0.1*[cos(((2*pi)/5)*i);sin(((2*pi)/5)*i);cos(((2*pi)/5)*i)*sin(((2*pi)/5)*i)];


ll=ll+1;
 end
% motion initialization
w=zeros(3,samples);

q=zeros(4,samples);

theta=2.3;

l=[1 1 1];

ll=l/norm(l);

l1=ll(1);

l2=ll(2);

l3=ll(3);

q(:,1)=[cos(theta/2) l1*sin(theta/2) l2*sin(theta/2) l3*sin(theta/2)];

y1=zeros(3,samples);

y2=zeros(3,samples);

% initial angular rate in rad/sec
w(:,1)=[0.3 0.2 0.1]';

y=zeros(6,samples);

d=0.001; % measurement noise variance (the two sensors are considered uncorelated)

k=1;

for i=0:dt:simulation_time
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
    
    % Lie-verlet integration scheme
    w_imp=w(:,k)+((dt/2)*inv(inertia)*((skew(inertia*w(:,k))*w(:,k))+nominal_input+model_error(:,k)));
    
    omega=0.5*[0 -w_imp'; w_imp -skew(w_imp)];
    
    q(:,k+1)=expm(dt*omega)*q(:,k);
    
    % newton solver
    norma=100;
    
    w_n=w_imp;
    
    while norma>(10^(-4))
        
        fw_n=-w_n+w_imp+(dt/2)*inv(inertia)*(skew(inertia*w_n)*w_n+nominal_input+model_error(:,k));
        
        jw_n= -eye(3)+((dt/2)*inv(inertia)*(skew(inertia*w_n)-skew(w_n)*inertia));
        
        w_n_1=w_n-inv(jw_n)*fw_n;
        
        w_n=w_n_1;
        
        norma=norm(fw_n,1);
    end
    
    w(:,k+1)=w_n_1;
    
    y1(:,k)=output_matrix(q(:,k))'*tracking_direction1+d*randn(3,1);
    
    y2(:,k)=output_matrix(q(:,k))'*tracking_direction2+d*randn(3,1);
    
    y(:,k)=[y1(:,k);y2(:,k)];
    
    k=k+1;
end

% filter's initialization
x(:,1)=[1 0 0 0 0 0 0]'; % state vector initialization

q_est(:,1)=[1 0 0 0]'; % quaternion initialization

w_est(:,1)=[0 0 0]'; % rate initialization

P(:,:,1)=0.1*eye(7); % conditional state covariance initialization

process_noise_covariance=0.00000001*eye(7);

output_noise_covariance=0.0000000001*eye(6);

errorw=zeros(3,samples); % rate error initialization

errorR=zeros(1,samples); % otientation error initialization

y1_est=zeros(3,samples);

y1_est(:,1)=output_matrix([1 0 0 0])'*tracking_direction1; % initial estimated output_1

y2_est=zeros(3,samples);

y2_est(:,1)=output_matrix([1 0 0 0])'*tracking_direction2; % initial estimated output_1

y_est(:,1)=[y1_est(:,1)' y2_est(:,1)']';

L=7;

k=1;

for i=0:dt:simulation_time
    
    y_est(:,k)=[y1_est(:,k)' y2_est(:,k)']';
     
    % sigma points calulation
    [point,weight_c,weight_m]=sigma_points(x(:,k),P(:,:,k),L,process_noise_covariance);
    
    nominal_input=[sin(((2*pi)/3)*i);cos(((2*pi)/1)*i);sin(((2*pi)/5)*i)];
        
    predicted_sigma_points = zeros(7,(2*L));
    
    
    % propagation of predicted sigma points 
    for l=1:(2*L)
        
        predicted_sigma_points(:,l)=prediction(point(:,l),nominal_input,inertia,dt); 
        
    end
    % determine the predicted state with conventional average
    predicted_state=prediction_state(predicted_sigma_points,weight_m,L);
    
    % determine the conditional state covariance
    state_covariance=covx(weight_c,predicted_sigma_points,predicted_state,process_noise_covariance,L);
    
    output_sigma_points = zeros(6,(2*L));
    
    % output sigma points 
    for l=1:(2*L)
        
        output_sigma_points(:,l)=out(predicted_sigma_points(:,l),tracking_direction1,tracking_direction2); % exei ena reprojection mesa
        
    end
    
    % determine the predicted ooutput with conventional average
    predicted_output=prediction_output(output_sigma_points,weight_m,L);
    
    y1_est(:,k+1)=predicted_output(1:3);
    
    y2_est(:,k+1)=predicted_output(4:6);
    
    % determine the output covariance matrix
    output_covariance=covy(weight_c,output_sigma_points,predicted_output,output_noise_covariance,L);
    
    % determine the cross covariance matrix
    cross_covariance_matrix=cross(weight_c,predicted_sigma_points,output_sigma_points,...
        predicted_state,predicted_output,L);
    
    % determine the UKF's gain matrix
    gain=cross_covariance_matrix/output_covariance;
    
    % state correction
    x(:,k+1)=predicted_state+gain*(y(:,k)-y_est(:,k));
    
    % quaterion re-projection
    x(1:4,k+1)=x(1:4,k+1)/norm(x(1:4,k+1));
    
    % conditional state covariance correction after new observations become
    % available
    P(:,:,k+1)=state_covariance-gain*output_covariance*gain';
    
    % corrected quaterion
    q_est(:,k+1)=x(1:4,k+1);
    
    %corrected rate
    w_est(:,k+1)=x(5:7,k+1);
    
    % rate error
    errorw(:,k)=(w_est(:,k)-w(:,k));
   
    % //calculation or orientation error
    DD=output_matrix(q(:,k));
    
    DD_est=output_matrix(q_est(:,k));
    
    tr=trace(eye(3)-(DD'*DD_est));
    
    errorR(k)=acosd(1-(tr/2));
    %//
    
 k=k+1;
    
end




figure(1)

t=linspace(0,simulation_time,samples);
plot(t,movmean(errorR(1:samples),1),'-','Linewidth',0.5,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_R (deg)$','interpreter','latex')
set(gca,'FontSize',16)
grid on


figure(2)
subplot(3,1,1)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(1,1:samples),1).*57.2957795,'-','Linewidth',0.5,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_{\omega}^x \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on


subplot(3,1,2)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(2,1:samples),1).*57.2957795,'-','Linewidth',0.5,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
ylabel('$e_{\omega}^y \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on

subplot(3,1,3)
t=linspace(0,simulation_time,samples);
plot(t,movmean(errorw(3,1:samples),1).*57.2957795,'-','Linewidth',0.5,'color','blue')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',16)
xlabel('$ time (sec)$','interpreter','latex')
ylabel('$e_{\omega}^z \ (deg/sec)$','interpreter','latex')
set(gca,'FontSize',16)
grid on






