% Draw metric
% Author:       Wei Zhang
% Email:        zhangwei@sustech.edu.nc
% Date:         2021.08.23

clear all;

metric_nc = '/home/zhangw/work/cgfd2dwave_2lay/00col/output/metric.nc'
grid_nc   = '/home/zhangw/work/cgfd2dwave_2lay/00col/output/coord.nc'

% in nc, dim order: time, z, x
phy_start  = [ 0    , 0 ]
%phy_count  = [ 101, 300]
phy_count  = [ inf, inf]
phy_stride = [ 1, 1]

%-- metric parameters list
vars_list = { 'jac', 'xi_x', 'xi_z', 'zeta_x', 'zeta_z' };
fig_dir   = './fig'

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
else
  x = nc_varget(grid_nc, 'x', grid_start, grid_count, grid_stride);
  z = nc_varget(grid_nc, 'z', grid_start, grid_count, grid_stride);
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

num_var = length(vars_list);

for ivar = 1 : num_var

  vars_list{ivar}

  var = nc_varget(metric_nc, vars_list{ivar}, grid_start, grid_count, grid_stride);

  %-- figure plot
  hid(ivar) = figure;
  set(hid(ivar),'BackingStore','on');
  
  %-- show
  pcolor(x, z, var);
  
  xlabel(['X axis (' str_unit ')']);
  ylabel(['Z axis (' str_unit ')']);
     
  set(gca,'layer','top');
  set(gcf,'color','white','renderer','painters');
  
  % shading
  % shading interp;
  shading flat;
  
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
  
  axis tight
  
  cid = colorbar('vert');
  %cid.Label.String='(km/s)';
  %cid.Label.String='g/cm^3';
  
  title(vars_list{ivar});

  drawnow
  
  % save and print figure
  flag_print  = 1;
  if flag_print
      width= 500;
      height=500;
      set(gcf,'paperpositionmode','manual');
      set(gcf,'paperunits','points');
      set(gcf,'papersize',[width,height]);
      set(gcf,'paperposition',[0,0,width,height]);
      ou_fname = [fig_dir, '/', vars_list{ivar}, '.png']
      print(gcf, ou_fname, '-dpng');
  end
  
end % ivar
