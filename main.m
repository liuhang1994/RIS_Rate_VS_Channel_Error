% if 1
clear all;
M=10; %# of antennas;
Nlist=[40,200,400]; %# of LIS components;

sigma=10^(-11); %noise power;
p=5;
x=[0,0,25]; %BS location
y=[0,80,40]; %LIS location
du=100;
z=[10,du,1.5]; %user location
beta1=(sqrt(norm(y-z)^2))^(-3.76)*10^(-3.53);
beta2=(sqrt(norm(x-y)^2))^(-3.76)*10^(-3.53);
% beta=beta1*beta2;
beta=(sqrt(norm(x-y)^2+norm(y-z)^2))^(-3.76)*10^(-3.53);
theta_los_1=atan(80/15);
theta_los_2=pi-theta_los_1;
phi_los1=pi/2;
phi_los2=3*pi/2;
a=exp(1i*pi*(0:M-1)'*sin(theta_los_1)*sin(phi_los1));
theta=atan(sqrt(10^2+(du-80)^2)/(40-1.5));
phi=3*pi/2;
r_angle=sin(theta)*sin(phi);
taulist=0:0.1:1;
taulist=-40:5:0;
bitlist=[-1,0];
ntest=1000;
Received_SNR=nan(ntest,length(taulist),length(bitlist),length(Nlist));

for nn=1:length(Nlist)
    N=Nlist(nn);
    b=exp(1i*pi*(0:N-1)'*sin(theta_los_2)*sin(phi_los2));
    H_1=a*b';
    R=eye(N,N); %sqrtm
    R_est=R;
    % for i=1:N
    %     for j=1:N
    %         R(i,j)=exp((i-j)*r_angle*1i*pi);
    %     end
    % end
    for i=1:ntest
        h_2=sqrt(0.5)*(1i*randn(N,1)+randn(N,1));
        for j=1:length(taulist)
            MSE=10^(taulist(j)/10);
            tau=2*MSE-MSE^2;
            h_2_est=sqrt(1-tau)*h_2+sqrt(tau)*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            b_est=sqrt(1-tau)*b+sqrt(tau)*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            %         h_2_est=h_2+sqrt(tau^2/(tau^2+1))*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            %         b_est=b+sqrt(tau^2/(tau^2+1))*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            b_est=exp(1i*angle(b_est));
            %                     H_1_est=sqrt(1-tau^2)*H_1+tau*sqrt(0.5)*(randn(M,N)+1i*randn(M,N));
            %                     [U,S,V]=svd(H_1_est);
            %                     b_est=V(:,1);
            %                     b_est=b_est*sqrt(N);
            %                     b_est=exp(-1i*angle(b_est(1)))*b_est;
            %optimal LIS design
            g_bar=((diag(h_2_est')*sqrtm(R_est)).')*b_est;
            v_star0=exp(1i*angle(g_bar));
            for bn=1:length(bitlist)
                bit=bitlist(bn);
                v_star=proj_lis(v_star0,bit);
                Phi=diag(v_star);
                h=H_1*Phi*sqrtm(R)*h_2*sqrt(beta);
                %h_est=H_1_est*Phi*sqrtm(R_est)*h_2_est*sqrt(beta);
                %             h_est=h;
                g=h/norm(h);
                SNR=abs(h'*g)^2*p/sigma;
                Received_SNR(i,j,bn,nn)=(SNR);
                %         fprintf('Trial %d Rate %f\n',i,Received_SNR(i,j));
            end
        end
    end
end
% end
Result=(mean(log2(1+Received_SNR)));

% Result=mean(log2(1+Received_SNR)>13)
Tri_N1=Result(1,:,1,1);
Infi_N1=Result(1,:,2,1);
Tri_N2=Result(1,:,1,2);
Infi_N2=Result(1,:,2,2);
Tri_N3=Result(1,:,1,3);
Infi_N3=Result(1,:,2,3);
% OneBit=Result(1,:,3)/Result(1,1,3);
% TwoBit=Result(1,:,4)/Result(1,1,4);
Tri_cur_N1 = fit(taulist',Tri_N1','linearinterp');
Infi_cur_N1 = fit(taulist',Infi_N1','smoothingspline');
Tri_cur_N2 = fit(taulist',Tri_N2','smoothingspline');
Infi_cur_N2 = fit(taulist',Infi_N2','smoothingspline');
Tri_cur_N3 = fit(taulist',Tri_N3','smoothingspline');
Infi_cur_N3 = fit(taulist',Infi_N3','smoothingspline');
% OneBit_cur = fit(taulist',OneBit','smoothingspline');
% TwoBit_cur = fit(taulist',TwoBit','smoothingspline');
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