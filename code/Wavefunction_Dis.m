%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavefunction distribution_lattice_Model_OPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
sigma_x=[0,1;1,0];sigma_y=1j*[0,-1;1,0];sigma_z=[1,0;0,-1];

sigma_zeeman=kron(eye(4),sigma_z);
% spin_operator=kron(eye(2),sigma_z);
% orbit_operator=kron(sigma_z,eye(2));

kx_list=[2.2,2.4,3.9,4.1];n_kx=numel(kx_list);

t=1;a=1;lambda=0.18;tso=0.2;
Lx=a;
% sqrt(3)*a;
n_y=101;n_DOF=8;Ham=zeros(n_y*n_DOF,n_y*n_DOF);
Ene=zeros(n_kx,n_y*n_DOF);
position=zeros(n_kx,n_y*n_DOF);

y_operator_temp=kron(linspace(-1,1,n_y),ones(1,8));
y_operator=diag(y_operator_temp,0);
mat_temp=[1,1,0,0,1,1,0,0;0,0,1,1,0,0,1,1];
wav_pos_operator=kron(eye(n_y),mat_temp);
wav_pos=zeros(n_kx,n_y*2);
pos_list=linspace(1,n_y*2,n_y*2);

for i_kx=1:n_kx
    
    kx=kx_list(i_kx);
%     same block
    temp=[1+exp(-1j*kx*Lx),1,1+exp(1j*kx*Lx)];
    h_hopping=diag(temp,1);
    h_hopping=h_hopping+h_hopping';
    H_hopping=kron(h_hopping,eye(2));
    temp=[-1j,0,0];xi_y_temp=diag(temp,1);xi_y_temp=xi_y_temp+xi_y_temp';
    temp=[1,0,0];xi_x_temp=diag(temp,1);xi_x_temp=xi_x_temp+xi_x_temp';
%     xi_{ij}=sigma_i xi_j
    xi_xy=kron(xi_y_temp,sigma_x);xi_yy=kron(xi_y_temp,sigma_y);
    xi_xx=kron(xi_x_temp,sigma_x);xi_yx=kron(xi_x_temp,sigma_y);
%     eta_{ij}=sigma_i eta_j
    temp=[0,0,-1j];eta_y_temp=diag(temp,1);eta_y_temp=eta_y_temp+eta_y_temp';
    temp=[0,0,1];eta_x_temp=diag(temp,1);eta_x_temp=eta_x_temp+eta_x_temp';
    eta_xy=kron(eta_y_temp,sigma_x);eta_yy=kron(eta_y_temp,sigma_y);
    eta_xx=kron(eta_x_temp,sigma_x);eta_yx=kron(eta_x_temp,sigma_y);
%     23 rashba
    temp=[0,-1j,0];matrix_23_temp=diag(temp,1);matrix_23_temp=matrix_23_temp+matrix_23_temp';
    matrix_23=kron(matrix_23_temp,sigma_x);
    for i_y=1:n_y
        Ham(n_DOF*(i_y-1)+1:n_DOF*(i_y),n_DOF*(i_y-1)+1:n_DOF*(i_y))=...
        -t*H_hopping+lambda*sigma_zeeman...
        +tso/2*(1+cos(kx*Lx))*xi_xy...
        -tso*sqrt(3)/2*(cos(kx*Lx)-1)*xi_yy...
        -tso/2*sin(kx*Lx)*xi_xx...
        +tso*sqrt(3)/2*sin(kx*Lx)*xi_yx...
        +tso/2*(1+cos(kx*Lx))*eta_xy...
        -tso*sqrt(3)/2*(cos(kx*Lx)+1)*eta_yy...
        +tso/2*sin(kx*Lx)*eta_xx...
        +tso*sqrt(3)/2*sin(kx*Lx)*eta_yx...
        +tso*matrix_23;
    end
    
%     different block
    h_hopping=diag([1],3);
    H_hopping=kron(h_hopping,eye(2));
    matrix_14_temp=diag([1j],3);
    matrix_14=kron(matrix_14_temp,sigma_x);
    i_y=1;
    Ham(n_DOF*(i_y-1)+1:n_DOF*(i_y),n_DOF*(i_y)+1:n_DOF*(i_y+1))=...
            -t*H_hopping+tso*matrix_14;
    for i_y=2:n_y-1
        Ham(n_DOF*(i_y-1)+1:n_DOF*(i_y),n_DOF*(i_y)+1:n_DOF*(i_y+1))=...
            -t*H_hopping+tso*matrix_14;
        Ham(n_DOF*(i_y-1)+1:n_DOF*(i_y),n_DOF*(i_y-2)+1:n_DOF*(i_y-1))=...
            Ham(n_DOF*(i_y-2)+1:n_DOF*(i_y-1),n_DOF*(i_y-1)+1:n_DOF*(i_y))';
    end
    i_y=n_y;
    Ham(n_DOF*(i_y-1)+1:n_DOF*(i_y),n_DOF*(i_y-2)+1:n_DOF*(i_y-1))=...
        Ham(n_DOF*(i_y-2)+1:n_DOF*(i_y-1),n_DOF*(i_y-1)+1:n_DOF*(i_y))';   
    [vec,ene_temp]=eig(Ham);
    wav=abs(vec(:,floor(n_y*n_DOF/2)+1));
    
    wav_pos(i_kx,:)=wav_pos_operator*wav;
    Ene(i_kx,:)=diag(ene_temp);
    
end
% figure
% plot(kx_list,Ene,'b');ylim([-1,1])

%r:K;b:K';<: negative velocity; >:positive velocity;
figure
plot(pos_list,wav_pos(1,:),'r<-','Linewidth',1,'MarkerSize',8);hold on;
plot(pos_list,wav_pos(2,:),'r>-.','Linewidth',1,'MarkerSize',8);
plot(pos_list,wav_pos(3,:),'b<-','Linewidth',1,'MarkerSize',8);
plot(pos_list,wav_pos(4,:),'b>-.','Linewidth',1,'MarkerSize',8);
xlim([1,n_y*2]);
ylabel('$\vert \psi\vert^2$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
xlabel('$y$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
legend('kx=2.2','kx=2.4','kx=3.9','kx=4.1');
hold off;
