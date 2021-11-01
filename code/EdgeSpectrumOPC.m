%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EdgeSpectrum_lattice_Model_OPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
sigma_x=[0,1;1,0];sigma_y=1j*[0,-1;1,0];sigma_z=[1,0;0,-1];

sigma_zeeman=kron(eye(4),sigma_z);
stagger=diag([1,-1,1,-1],0);
sigma_zeeman_stagger=kron(stagger,sigma_z);
% spin_operator=kron(eye(2),sigma_z);
% orbit_operator=kron(sigma_z,eye(2));

kx_min=0;kx_max=2*pi;n_kx=201;
kx_list=linspace(kx_min,kx_max,n_kx);

t=1;a=1;lambda=0.18;tso=0.2;V=0.;lambda_stagger=0.;
Lx=a;sqrt(3)*a;
n_y=51;n_DOF=8;Ham=zeros(n_y*n_DOF,n_y*n_DOF);
Ene=zeros(n_kx,n_y*n_DOF);
position=zeros(n_kx,n_y*n_DOF);
sublattice=zeros(n_kx,n_y*n_DOF);

y_operator_temp=kron(linspace(-1,1,n_y),ones(1,8));
y_operator=diag(y_operator_temp,0);
sublattice_temp1=kron([1,-1,1,-1],ones(1,2));
sublattice_temp=kron(ones(1,n_y),sublattice_temp1);
sublattice_operator=diag(sublattice_temp,0);
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
    matrix_potential=diag([1,1,-1,-1,1,1,-1,-1],0);
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
        +tso*matrix_23+V*matrix_potential+lambda_stagger*sigma_zeeman_stagger;
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
    position(i_kx,:)=real(diag(vec'*y_operator*vec));
    sublattice(i_kx,:)=diag(vec'*sublattice_operator*vec);
    Ene(i_kx,:)=diag(ene_temp);
end
% figure
% plot(kx_list,Ene,'b');ylim([-1,1])

Kx_list_temp=kron(ones(1,n_y*n_DOF),kx_list');
Kx_list=reshape(Kx_list_temp,1,numel(Kx_list_temp));
Ene_scatter=reshape(Ene,1,numel(Ene));
position_scatter=reshape(position,1,numel(position));

sublattice_scatter=reshape(sublattice,1,numel(sublattice));

% figure
subplot(1,2,1)
scatter(Kx_list,Ene_scatter,[],position_scatter,'filled');colorbar;
c=colorbar;colormap('cool');c.Label.String='\langle y\rangle';
caxis([-1,1]);
xlim([kx_min,kx_max]);ylim([-1,1]);
set(c,'Fontsize',16,'Fontname', 'Times New Roman');
set(gca,'fontsize',18,'Fontname', 'Times New Roman');
xlabel('$k_x$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
ylabel('$E$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
grid on;hold off

subplot(1,2,2)
scatter(Kx_list,Ene_scatter,[],sublattice_scatter,'filled');colorbar;
c=colorbar;colormap('cool');c.Label.String='\langle \tau_z\rangle';
caxis([-1,1]);
xlim([kx_min,kx_max]);ylim([-1,1]);
set(c,'Fontsize',16,'Fontname', 'Times New Roman');
set(gca,'fontsize',18,'Fontname', 'Times New Roman');
xlabel('$k_x$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
ylabel('$E$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
grid on;hold off