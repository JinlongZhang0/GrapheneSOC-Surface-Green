%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Band Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
sigma_x=[0,1;1,0];sigma_y=1j*[0,-1;1,0];sigma_z=[1,0;0,-1];
% sigma_{ij}=tau_i sigma_j
sigma_x0=kron(sigma_x,eye(2));sigma_y0=kron(sigma_y,eye(2));
sigma_z0=kron(sigma_z,eye(2));
sigma_0x=kron(eye(2),sigma_x);sigma_0y=kron(eye(2),sigma_y);
sigma_0z=kron(eye(2),sigma_z);
sigma_zx=kron(sigma_z,sigma_x);sigma_zy=kron(sigma_z,sigma_y);
sigma_xx=kron(sigma_x,sigma_x);sigma_xy=kron(sigma_x,sigma_y);
sigma_yx=kron(sigma_y,sigma_x);sigma_yy=kron(sigma_y,sigma_y);
sigma_z0=kron(sigma_z,eye(2));sigma_zz=kron(sigma_z,sigma_z);
spin_operator=kron(eye(2),sigma_z);
orbit_operator=kron(sigma_z,eye(2));

kx_min=-pi;kx_max=pi;n_kx=211;
ky_min=-pi;ky_max=pi;n_ky=171;
kx_list=linspace(kx_min,kx_max,n_kx);
ky_list=linspace(ky_min,ky_max,n_ky);
Ham=zeros(4,4);
t=1;a=1;lambda=0.;tso=0.;V=0.;lambda_stagger=0.;
Ene=zeros(n_kx,n_ky,4);
Spin=zeros(n_kx,n_ky,4);
orbit=zeros(n_kx,n_ky,4);
for i_kx=1:n_kx
    kx=kx_list(i_kx);
    for i_ky=1:n_ky
        ky=ky_list(i_ky);
        Ham=-t*(2*cos(sqrt(3)/2*kx*a)*cos(ky/2*a)+cos(ky*a))*sigma_x0...
            -t*(-2*cos(sqrt(3)/2*kx*a)*sin(ky/2*a)+sin(ky*a))*sigma_y0...
            +lambda*sigma_0z...
            -tso*(cos(ky*a/2)*cos(sqrt(3)/2*kx*a)-cos(ky*a))*sigma_yx...
            -tso*sqrt(3)*sin(ky*a/2)*sin(sqrt(3)/2*kx*a)*sigma_yy...
            -tso*(sin(ky*a/2)*cos(sqrt(3)/2*kx*a)+sin(ky*a))*sigma_xx...
            -tso*sqrt(3)*(-cos(ky/2*a)*sin(sqrt(3)/2*kx*a))*sigma_xy+V*sigma_z0...
            +lambda_stagger*sigma_zz;
        [vec,ene_temp]=eig(Ham);
        Ene(i_kx,i_ky,:)=diag(ene_temp);
        Spin(i_kx,i_ky,:)=diag(vec'*spin_operator*vec);
        orbit(i_kx,i_ky,:)=diag(vec'*orbit_operator*vec);
    end
end
figure
for i_band=1:4
    mesh(kx_list,ky_list,Ene(:,:,i_band)');hold on;
end
colorbar;xlabel('k_x');ylabel('k_y');

figure
subplot(1,2,1)
Kx_list_temp=kron(ones(1,4),kx_list');
Kx_list=reshape(Kx_list_temp,1,numel(Kx_list_temp));
temp=Ene(:,floor(n_ky/2)+1,:);Ene_re=reshape(temp,numel(kx_list),4);
Ene_Re=reshape(Ene_re,1,numel(Ene_re));
temp=Spin(:,floor(n_ky/2)+1,:);Spin_re=reshape(temp,numel(kx_list),4);
Spin_Re=reshape(Spin_re,1,numel(Spin_re));
scatter(Kx_list,Ene_Re,[],Spin_Re,'filled');
c=colorbar;colormap('cool');
caxis([-1,1]);
c.Label.String='\langle\sigma_z\rangle';
set(c,'Fontsize',16,'Fontname', 'Times New Roman');
xlabel('$k_x$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
ylabel('$E$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
set(gca,'fontsize',18,'Fontname', 'Times New Roman');
grid on;
subplot(1,2,2)
temp=orbit(:,floor(n_ky/2)+1,:);orbit_re=reshape(temp,numel(kx_list),4);
orbit_Re=reshape(orbit_re,1,numel(orbit_re));
scatter(Kx_list,Ene_Re,[],orbit_Re,'filled');
c=colorbar;colormap('cool');c.Label.String='\langle \tau_z\rangle';
caxis([-1,1]);
set(c,'Fontsize',16,'Fontname', 'Times New Roman');
set(gca,'fontsize',18,'Fontname', 'Times New Roman');
xlabel('$k_x$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
ylabel('$E$','Interpreter','latex','fontsize',21,'Fontname', 'Times New Roman');
grid on;
% 
% % figure
% % temp=Ene(:,floor(n_ky/2)+1,:);
% % Ene_re=reshape(temp,numel(kx_list),4);
% % plot(kx_list,Ene_re,'b') 
% 
% 
% % 
% % figure
% % temp=Ene(floor(n_kx/2)+1,:,:);
% % Ene_re=reshape(temp,numel(ky_list),4);
% % plot(ky_list,Ene_re)