clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
subs=[1,1];     % start from index '1'
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[8,8];

% figure control parameters
flag_km     = 0;
flag_emlast = 1;
flag_print  = 0;
flag_title  = 0;
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

[x,z]=gather_coord(output_dir,subs,subc,subt);

%- set coord unit
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
hid = figure;
set(hid,'BackingStore','on');

plot(x,z,'r-');
hold on
plot(x',z','r-');
  
xlabel(['X axis (' str_unit ')'],FontSize=35);
ylabel(['Y axis (' str_unit ')'],FontSize=35);
  
set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

axis tight;

% title
if flag_title
    gridtitle='XOZ-Grid';
    title(gridtitle);
end

set(gca,'FontSize',20);
% set(gcf,'position',[0,0,1920,1200]);


% print(gcf,['grid_model4.png'],'-r300','-dpng');
% save and print figure

