% if 1
clear all;
%Note: the channel model is given by
%h=H_1*Phi*h_2*beta
%where:
%h:cascaded User-RIS-BS channel
%H_1: BS-RIS channel, rank-1, H_1=a*b' with a,b denoting steering vectors
%   at BS and RIS, respectively
%Phi: diagnoal matrix, RIS phase shift
%h_2: User-RIS channel for the single-antenna user (only 1 user in the system)
%beta: User-RIS-BS cascaded path loss
%   loss
%The direct User-BS channel is assumed to be blocked.




M=10; %# of antennas;
Nlist=[40,200,400]; %# of LIS components;
sigma=10^(-11); %noise power,-110dBW
p=5; % transmit power,5W
x=[0,0,25]; %BS location (x,y,z)
y=[0,80,40]; %RIS location (x,y,z)
du=100;
z=[10,du,1.5]; %user location (x,y,z)

% beta1=(sqrt(norm(y-z)^2))^(-3.76)*10^(-3.53); %User-RIS path loss , no use
% beta2=(sqrt(norm(x-y)^2))^(-3.76)*10^(-3.53); %BS-RIS path loss, no use
beta=(sqrt(norm(x-y)^2+norm(y-z)^2))^(-3.76)*10^(-3.53); %Cascaded path loss

theta_los_1=atan(80/15);
theta_los_2=pi-theta_los_1;
phi_los1=pi/2;
phi_los2=3*pi/2;
a=exp(1i*pi*(0:M-1)'*sin(theta_los_1)*sin(phi_los1)); % RIS-BS channel BS steering vector
theta=atan(sqrt(10^2+(du-80)^2)/(40-1.5));
phi=3*pi/2;
r_angle=sin(theta)*sin(phi);

taulist=-40:5:0; %normalzied MSE in dB, the x-axis in the figure
bitlist=[-1,0]; %0 means perfect phase shift, -1 means no RIS
ntest=1000; %number of Monte Carlo trials
Received_SNR=nan(ntest,length(taulist),length(bitlist),length(Nlist)); %store the result for plot

for nn=1:length(Nlist)
    N=Nlist(nn); %RIS size
    b=exp(1i*pi*(0:N-1)'*sin(theta_los_2)*sin(phi_los2)); %RIS steering vector
    H_1=a*b'; %rank-one channel for BS-RIS
    R=eye(N,N); 
    R_est=R; % assume an i.i.d. channel the channel covariance is identity matrix
    for i=1:ntest
        h_2=sqrt(0.5)*(1i*randn(N,1)+randn(N,1)); % small scale fading coefficients
        for j=1:length(taulist)
            MSE=10^(taulist(j)/10); % channel error
            % error channel model h_est=sqrt(1-\tau)h+\sqrt{\tau}e for
            % Guassian nosie e, and \tau is used to control the MSE level
            tau=2*MSE-MSE^2; 
            h_2_est=sqrt(1-tau)*h_2+sqrt(tau)*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            % same channel estimation error model for b
            b_est=sqrt(1-tau)*b+sqrt(tau)*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            b_est=exp(1i*angle(b_est));
            %optimal LIS design
            %note that both active/passive beamforming are computed via
            %estimated CSI (in the existance of error)
            % active beamforming: maximum ratio transmission
            g_bar=((diag(h_2_est')*sqrtm(R_est)).')*b_est;
            %tentative phase shift vector
            v_star0=exp(1i*angle(g_bar));
            for bn=1:length(bitlist)
                % account for discrete phase shifts, if any
                % no effect in this code
                % v_star=v_star0
                bit=bitlist(bn);
                v_star=proj_lis(v_star0,bit);
                %phase shift matrix
                Phi=diag(v_star);
                %the true cascaded channel
                h=H_1*Phi*sqrtm(R)*h_2*sqrt(beta);
                %normalize active beamforming to have unit norm
                g=h/norm(h);
                %compute receive SNR
                SNR=abs(h'*g)^2*p/sigma;
                %store the result
                Received_SNR(i,j,bn,nn)=(SNR);
                %         fprintf('Trial %d Rate %f\n',i,Received_SNR(i,j));
            end
        end
    end
end
% end
% compute rate
Result=(mean(log2(1+Received_SNR)));

%Store the result
Tri_N1=Result(1,:,1,1);
Infi_N1=Result(1,:,2,1);
Tri_N2=Result(1,:,1,2);
Infi_N2=Result(1,:,2,2);
Tri_N3=Result(1,:,1,3);
Infi_N3=Result(1,:,2,3);

% fit curves
Tri_cur_N1 = fit(taulist',Tri_N1','linearinterp');
Infi_cur_N1 = fit(taulist',Infi_N1','smoothingspline');
Tri_cur_N2 = fit(taulist',Tri_N2','smoothingspline');
Infi_cur_N2 = fit(taulist',Infi_N2','smoothingspline');
Tri_cur_N3 = fit(taulist',Tri_N3','smoothingspline');
Infi_cur_N3 = fit(taulist',Infi_N3','smoothingspline');

%plot
close all
figure( 'Position', [20 20 700 500])
hold on
plot(taulist,Infi_cur_N1(taulist(:)),'+-','Color','blue','LineWidth',1.5,'MarkerSize',6);
plot(taulist,Tri_cur_N1(taulist(:)),':','Color','blue','LineWidth',1.5,'MarkerSize',6);
plot(taulist,Infi_cur_N2(taulist(:)),'>-','Color','red','LineWidth',1.5,'MarkerSize',6);
plot(taulist,Tri_cur_N2(taulist(:)),':','Color','red','LineWidth',1.5,'MarkerSize',6);
plot(taulist,Infi_cur_N3(taulist(:)),'o-','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',6);
plot(taulist,Tri_cur_N3(taulist(:)),':','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',6);
axis([-40 0 0 25])
grid on;
xlabel({'Normalized MSE (dB)'},'FontSize',15)
ylabel('Acheiveable Rate (bit/s/Hz)','FontSize',15);
set(gca,'FontSize',15);
legend({'Optimal Phase Shift, N=40','Random Phase Shift, N=40','Optimal Phase Shift, N=200','Random Phase Shift, N=200','Optimal Phase Shift, N=400','Random Phase Shift, N=400'},'FontSize',10,'Location','southeast')