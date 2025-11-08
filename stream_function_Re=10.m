function plot_Re10(psi)
Re=10;
%%%%% ci-theta grid %%%%%
n=size(psi, 1); m=size(psi,2);
N=n-1; M=m-1;
h=pi/M; µ grid spacing
xi=(0:N)*h; theta=(0:M)*h;
[XI, THETA] = meshgrid(xi,theta);
%%%%%% x-y grid %%%%%%%
nx=640; ny=480/2; % number of pixels in x and half of y
xmin=-1.5; xmax=2.5; ymax=(xmax-xmin)*ny/nx; ymin=-ymax;
x=linspace(xmin,xmax,nx+1); y=linspace(0,ymax,ny+1);
[X,Y]=meshgrid(x,y);
%%%%% construct interpolation points %%%%%%
xi_i=0.5*log(X.^2+Y.^2);
theta_i=wrapTo2Pi(atan2(Y,X));
%%%%%% interpolate %%%%%%
psi_xy=interp2(XI,THETA,psi',xi_i,theta_i);
%%%%% set psi zero inside cylinder %%%%%%
psi_xy(xi_i<0)=0;
%%%%%% scale contour levels %%%%%%
%%%%% negative values have same range as positive values %%%%%%
psi_min=min(psi_xy(:));
psi_max=max(psi_xi(:));
psi_xy(psi_xy<0)=psi_xy(psi_xy<0)/abs(psi_min);
psi_xy(psi_xy>0)=psi_xy(psi_xy>0)/abs(psi_max);
%%%%%% set colormap for contours %%%%%
levels=linspace(-1,1,1000);
cmap=flipud(jet(length(levels)));
colormap(cmap);
%%%%% plot color contours %%%%%
imagesc(x,y,psi_xy); hold on;
imagesc(x,-y,-psi_xy); % negative values of y
%%%%%% plot contour lines %%%%%
v=[-0.9:0.2:-0.1,0,0.0005,0.001,0.002:0.004:0.01,0.02:0.04:0.9];
contour(X,Y,psi_xy,v,'LineColor','k');
contour(X,-Y,-psi_xy,-v,'LineColor','k');
%%%%% draw black circle for cylinder %%%%%
t=linspace(0,2*pi,1000);
a=cos(t); b=sin(t);
fill(a,b,[0 0 0];
%%%%%% neaten plot %%%%%%
set(gca,'YDir','normal');
axis([xmin xmax ymin ymax]); set(gcf,'color','w'); axis odd; axis equal;
text(xmin+0.75*(xmax-xmin),ymin+0.08*(ymax-ymin);...
  ['Re = ', num2str(Re, '%3.0f')],'FontSize', 22,'Color','k');


function [psi, omega] = flow_around_cylinder_steady
Re=10; 
%%%%% define the grid %%%%%
n=101; m=101; % number of grid points
N=n-1; M=m-1; % number of grid intervals
h=pi/M; % grid spacing based on theta variable
xi=(0:N)*h; theta=(0:M)*h; % xi and theta variables on the grid
%%%%% Initialize the flow fields %%%%%
psi=zeros(n,m);
omega=zeros(n,m);
psi(n,:)=exp(xi(n))*sin(theta(:)); % Write the free stream bc here
%%%%% Set relax params, tol, extra variables %%%%%
r_psi=1.8; % Set the relaxation parameter here, psi equation
r_omega=0.9; % Set the relaxation parameter here, omega equation
delta=1.e-08; % error tolerance
error=2*delta; % initialize error variable
%%%%% Add any additional variable definitions here %%%%%
...
...
%%%%% Main SOR Loop %%%%%
while (error > delta)
    psi_old = psi; omega_old = omega;
    for i=2:n-1
        for j=2:m-1
            psi(i,j)=(1-r_psi)*psi(i,j)+(r_psi/4)*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+h^2*exp(2*xi(i))*omega(i,j)); % Write psi equation here
        end
    end
    error_psi=max(abs(psi(:)-psi_old(:)));
    omega(1,:)=(1/(2*h^2))*(psi(3,:)-8*psi(2,:)); % Write the boundary condition here
    for i=2:n-1
        for j=2:m-1
            f(i,j)=((psi(i+1,j)-psi(i-1,j))*(omega(i,j+1)-omega(i,j-1)))-((psi(i,j+1)-psi(i,j-1))*(omega(i+1,j)-omega(i-1,j)));
            omega(i,j)= (1-r_omega)*omega(i,j)+(r_omega/4)*(omega(i+1,j)+omega(i-1,j)+omega(i,j+1)+omega(i,j-1)+((Re/8)*f(i,j)));% Write omega equation here
        end
    end
    error_omega=max(abs(omega(:)-omega_old(:)));
    error=max(error_psi, error_omega);
end
plot_Re10(psi);
