% Draw seismic wavefield snapshot by using pcolor style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.05.31

clear all;
close all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json'
output_dir='../project/output'

% which snapshot to plot
id=1;

% variable and time to plot
varnm='Vz';
ns=10;
ne=600;
nt=10;

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
savegif = 0;

% scl_caxis=[-1.0 1.0];
filename1 = ['Vz2.gif'];
scl_daspect =[1 1 1];
clrmp       = 'jetwr';
taut=0.5;
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,z]=gather_coord(parfnm,output_dir);

% coordinate unit
str_unit='m';
if flag_km
   x=x/1e3;
   z=z/1e3;
   str_unit='km';
end

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% snapshot show
for nlayer=ns:nt:ne
    
    [v,t]=gather_snap(parfnm,output_dir,nlayer,varnm);
    
    disp([ '  draw ' num2str(nlayer) 'th time step (t=' num2str(t) ')']);
    
    pcolor(x,z,v)
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
        
    
    set(gca,'layer','top');
    set(gcf,'color','white','renderer','painters');

    % axis image
    % shading
    % shading interp;
    shading flat;
    % colorbar range/scale
    if exist('scl_caxis')
        caxis(scl_caxis);
    end
    % axis daspect
    if exist('scl_daspect')
        daspect(scl_daspect);
    end
    % colormap and colorbar
    if exist('clrmp')
        colormap(clrmp);
    end
    colorbar('vert');
    
    %title
    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(taut);
    %save gif
    if savegif
      im=frame2im(getframe(gcf));
      [imind,map]=rgb2ind(im,256);
      if nlayer==ns
        imwrite(imind,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.5);
      else
        imwrite(imind,map,filename1,'gif','WriteMode','append','DelayTime',0.5);
      end
    end
    
    % save and print figure
    if flag_print==1
        width= 500;
        height=500;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        fnm_out=[varnm '_ndim_',num2str(nlayer,'%5.5i')];
        print(gcf,[fnm_out '.png'],'-dpng');
    end
    
end



