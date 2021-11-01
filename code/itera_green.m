%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative Green Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
sigma_x=[0,1;1,0];sigma_y=1j*[0,-1;1,0];sigma_z=[1,0;0,-1];

sigma_zeeman=kron(eye(4),sigma_z);
% spin_operator=kron(eye(2),sigma_z);
% orbit_operator=kron(sigma_z,eye(2));

kx_min=0;kx_max=2*pi;n_kx=151;kx_list=linspace(kx_min,kx_max,n_kx);

t=1;a=1;lambda=0.18;tso=0.2;n_DOF=8;
E_min=-1.0;E_max=1.0;n_E=301;E_list=linspace(E_min,E_max,n_E);
E_imag_temp=0.3e-1;E_list=E_list+1j*E_imag_temp;
n_iter=10;

Lx=a;

DOS=zeros(n_kx,n_E);


for i_kx=1:n_kx
    kx=kx_list(i_kx);

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
    
    Ham_block=...
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
%  H_R0_R1  c_{j+1}^\dagger c_j
    h_hopping=diag([1],3);
    H_hopping=kron(h_hopping,eye(2));
    matrix_14_temp=diag([1j],3);
    matrix_14=kron(matrix_14_temp,sigma_x);
    
    Ham_linking=-t*H_hopping+tso*matrix_14;
    
    for i_E=1:n_E
        E_imag=E_list(i_E);
        t0=inv(E_imag*eye(8)-Ham_block)*Ham_linking';
        t0_tilted=inv(E_imag*eye(8)-Ham_block)*Ham_linking;
        ti=t0;ti_tilted=t0_tilted;
        T=ti;T_i=ti_tilted;
        for i_iter=1:n_iter
            t_temp=pinv(eye(8)-ti*ti_tilted-ti_tilted*ti)*ti*ti;
            t_temp_iter=pinv(eye(8)-ti*ti_tilted-ti_tilted*ti)*ti_tilted*ti_tilted;
            ti=t_temp;
            ti_iter=t_temp_iter;
            T=T+T_i*ti;
            T_i=T_i*ti_tilted;
        end
        g_R0=inv(E_imag*eye(8)-Ham_block-Ham_linking*T);
        sum_R=Ham_linking*g_R0*Ham_linking';
        G_S=inv(E_imag*eye(8)-Ham_block-sum_R);
        DOS(i_kx,i_E)=-1/pi*imag(sum(diag(G_S)));
    end
end
% figure
subplot(1,2,2)
surf(kx_list,real(E_list),exp(DOS'));shading('interp');
view(0,90);xlim([kx_min,kx_max]);ylim([E_min,E_max]);
c=colorbar;
caxis([0,15]);c.Label.String='e^{DOS}';
set(c,'Fontsize',16,'Fontname', 'Times New Roman');
set(gca,'fontsize',18,'Fontname', 'Times New Roman');
xlabel('$k_x$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
ylabel('$E$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
colormap('jet');