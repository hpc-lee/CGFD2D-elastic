% Plot the grid
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.06
% Author:       Wei Zhang
% Email:        zhangwei@sustech.edu.nc
% Date:         2021.08.23

clear all;

grid_nc  = '/home/zhangw/work/cgfd2dwave_2lay/01stg/output/coord.nc'

% in nc, dim order: time, z, x
phy_start  = [ 0    , 0 ]
%phy_count  = [ 101, 300]
phy_count  = [ inf, inf]
phy_stride = [ 5, 5]

fig_dir       = './fig'

%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

% to local index with ghost for grid and media
grid_start  = nc_attget(grid_nc,-1,'local_index_of_first_physical_points');
grid_start  = fliplr(double(grid_start)) + phy_start;
grid_count  = phy_count;
grid_stride = phy_stride;

%-- load coord 
xvar_info = ncinfo(grid_nc,'x');
num_dim_x = length(xvar_info.Dimensions);
if num_dim_x == 1 % cart grid
  x = nc_varget(grid_nc, 'x', grid_start(2), grid_count(2), grid_stride(2));
  z = nc_varget(grid_nc, 'z', grid_start(1), grid_count(1), grid_stride(1));
  [x2d,z2d] = meshgrid(x,z);
else
  x2d = nc_varget(grid_nc, 'x', grid_start, grid_count, grid_stride);
  z2d = nc_varget(grid_nc, 'z', grid_start, grid_count, grid_stride);
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

% set others
my_daspect    = [1 1 1];
my_colormap   = 'parula';
my_pause_time = 0.5;
%my_caxis      = [-1,1] * 1e-2;

hid = figure;
set(hid,'BackingStore','on');

plot(x2d,z2d,'k-');
hold on
plot(x2d',z2d','k-');
  
xlabel(['X axis (' str_unit ')']);
ylabel(['Z axis (' str_unit ')']);

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');
  
  %- colorbar range/scale
  if exist('my_caxis','var')
      caxis(my_caxis);
  end
  %- axis daspect
  if exist('my_daspect','var')
      daspect(my_daspect);
  end
  %- colormap and colorbar
  if exist('my_colormap','var')
      colormap(my_colormap);
  end

  axis tight;

  gridtitle='XOZ-Grid';
  title(gridtitle);

% save and print figure
flag_print  = 1;
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    ou_fname = [fig_dir, '/', 'grid.png']
    print(gcf, ou_fname, '-dpng');
end


