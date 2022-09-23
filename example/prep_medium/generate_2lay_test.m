% Generate a modified Can4 model (from E2VP)
% reference: Chaljub et al, 2015, 3D numerical simulation of earthquake ground motion
%   in sedimentary basins: testing accuracy through stringent models
clear all;
close all;

% 3lay header: 
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij
media_type = 'elastic_isotropic'
if_write = 1;

NI = 4;

H(1) =     0;
H(2) = - 173*3;
H(3) = - 725*3;
H(4) = -1156*3;

% media info od the first layer
rho(1) = 1000;
vp(1)  = 1500;
vs(1)  = 1000;
par_grad = 0.0;
par_pow  = 1; 


for ni = 2:NI
    rho(ni) = rho(ni-1) * 1.5;
    vp(ni)  =  vp(ni-1) * 1.5;
    vs(ni)  =  vs(ni-1) * 1.5;
end

% x,z
layer1 = [-3000, 0; 8000, 0];
layer2 = [-3000, -519;  1000, -519;  2500, 0; 8000, 0];
layer3 = [-3000, -2175; 1000, -2175; 2500, 0; 8000, 0];
layer4 = [-3000, -3468; 1000, -3468; 2500, 0; 8000, 0];

% plot
figure;
for ni = 1:NI
    plot(eval(['layer',num2str(ni),'(:,1)']),eval(['layer',num2str(ni),'(:,2)']));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The build model');
end

if if_write
	fid = fopen('test.md2lay','w');
	fprintf(fid, '%s\n',media_type);
	fprintf(fid, '%d\n', NI);
    for ni = 1:NI
        layer = eval(['layer',num2str(ni)]);
        npoint = length(layer);
        fprintf(fid, '%d\n', npoint);
        for i = 1:npoint
            fprintf(fid, '%13.6f %13.6f ',layer(i,1), layer(i, 2));
            fprintf(fid, '%13.6f %6.2f %6.2f ', rho(ni), par_grad, par_pow);
            fprintf(fid, '%13.6f %6.2f %6.2f ', vp(ni), par_grad, par_pow);
            fprintf(fid, '%13.6f %6.2f %6.2f\n', vs(ni), par_grad, par_pow);
        end
    end
    fclose(fid);
end

