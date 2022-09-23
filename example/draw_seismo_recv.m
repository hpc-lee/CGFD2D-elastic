% Plot the seismograms on a line
% Author:   Yuanhang Huo, Wei Zhang
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.08.24

clear all;

%-- sac file
%sac_file = '/home/zhangw/work/cgfd2dwave_recv/02col/output/evt_by_par.recv_i200_k099.Vz.sac'
%sac_file = '/home/zhangw/work/cgfd2dwave_recv/03stg/output/evt_by_par.recv_i200_k099.Vz.sac'

%sac_file = '/home/zhangw/work/cgfd2dwave_2lay/00col/output/evt_by_par.recv_i200_k099.Vz.sac'
sac_file = './00/output/evt_by_par.recv_i200_k099.Vz.sac'

%-- get fname only
[filepath,name,ext] = fileparts(sac_file);
fname = [name, ext];

%-- read file
sacdata=rsac(sac_file);

var_val =sacdata(:,2);
var_time=sacdata(:,1);
    
%-- plot 
figure;

plot(var_time, var_val, 'k-', 'linewidth', 1.0);

xlabel('Time (s)');
ylabel('Amplitude');
title(fname, 'interpreter', 'none');
set(gcf,'color','white','renderer','painters');

%-- save and print figure
flag_print = 1;

fig_dir = './fig'

if flag_print

    mkdir(fig_dir);

    width= 800;
    height=400;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[fig_dir, '/', name, '.png'],'-dpng');
end

