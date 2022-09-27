% Plot the grid
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.06
% Author:       Wei Zhang
% Email:        zhangwei@sustech.edu.nc
% Date:         2021.08.23

clear all;
close all;
clc;
addmypath

% 2D plot slice all points

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';
%media_type = 'ac_iso';
media_type = 'el_iso';
% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';

% varable to plot 
%  'Vp', 'Vs', 'rho', 'lambda', 'mu' etc...
varnm = 'Vp';
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,z]=gather_coord(parfnm,output_dir);

switch varnm
    case 'Vp'
        rho=gather_media(parfnm,output_dir,'rho');
        if strcmp(media_type,'ac_iso') == 1
          kappa=gather_media(parfnm,output_dir,'kappa');
          v=( (kappa)./rho ).^0.5;
        elseif strcmp(media_type,'el_iso') == 1
          mu=gather_media(parfnm,output_dir,'mu');
          lambda=gather_media(parfnm,output_dir,'lambda');
          v=( (lambda+2*mu)./rho ).^0.5;
        end
        v=v/1e3;
    case 'Vs'
        rho=gather_media(parfnm,output_dir,'rho');
        mu=gather_media(parfnm,output_dir,'mu');
        v=( mu./rho ).^0.5;
        v=v/1e3;
    case 'rho'
        v=gather_media(parfnm,output_dir,varnm);
        v=v/1e3;
    otherwise
        v=gather_media(parfnm,output_dir,varnm);
end

%- set coord unit
flag_km     = 1;
if flag_km
   x=x/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------

% figure plot
hid=figure;
set(hid,'BackingStore','on');

pcolor(x,z,v);

xlabel(['X axis (' str_unit ')']);
ylabel(['Z axis (' str_unit ')']);

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% shading
% shading interp;
shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
end
% title
if flag_title
    title(varnm,'interpreter','none');
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm '.png'],'-dpng');
end
