% Draw seismic wavefield snapshot 
% Author:       Wei Zhang
% Email:        zhangwei@sustech.edu.nc
% Date:         2021.08.23

clear all;

snap_nc  = './00/output/volume_vel.nc'
grid_nc  = './00/output/coord.nc'

% variable and time to plot
var_name   ='Vz';
tdim_start = 50;
tdim_skip  = 50;
tdim_end   = 1500;

% in nc, dim order: time, z, x
snap_start  = [ 0    , 0 ]
snap_count  = [ 101, 300]
snap_stride = [ 1, 1]

%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

% to local index with ghost for grid and media
grid_start  = nc_attget(snap_nc,-1,'first_index_in_this_thread_with_ghosts');
grid_start  = fliplr(double(grid_start))
grid_count  = snap_count;
grid_stride = snap_stride;

%-- load coord and media
xvar_info = ncinfo(grid_nc,'x');
num_dim_x = length(xvar_info.Dimensions);
if num_dim_x == 1 % cart grid
  x = nc_varget(grid_nc, 'x', grid_start(2), grid_count(2), grid_stride(2));
  z = nc_varget(grid_nc, 'z', grid_start(1), grid_count(1), grid_stride(1));
else
  x = nc_varget(grid_nc, 'x', grid_start, grid_count, grid_stride);
  z = nc_varget(grid_nc, 'z', grid_start, grid_count, grid_stride);
end

%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------

% set coord unit
flag_km     = 1;
if flag_km
   x=x/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

% set others
my_daspect    = [1 1 1];
my_colormap   = 'parula';
my_pause_time = 0.5;
%my_caxis      = [-1,1] * 1e-2;
fig_dir       = './fig'

%-- loop to load and show snapshot

hid=figure;
set(hid,'BackingStore','on');

for it = tdim_start : tdim_skip : tdim_end
    
    dim_start  = [ it, snap_start ];
    dim_count  = [ 1 , snap_count ];
    dim_stride = [ 1 , snap_stride ];
    var  = nc_varget(snap_nc, var_name, dim_start, dim_count, dim_stride);
    t    = nc_varget(snap_nc, 'time', [it], [1], [1]);
    
    disp([ 'time step (t=' num2str(t) ')']);

    pcolor(x, z, var);
    
    set(gca,'layer','top');
    set(gcf,'color','white','renderer','painters');

    % axis image
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
    colorbar('vert');
    
    %title
    titlestr=['Snapshot of ' var_name ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(my_pause_time);
    
    % save and print figure
    flag_print  = 0;
    if flag_print==1
        width= 500;
        height=500;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        fnm_out=[fig_dir, '/', varnm, '_tdim_',num2str(it,'%5.5i')]
        print(gcf,[fnm_out '.png'],'-dpng');
    end
    
end


