function [v] = gather_metric(output_dir,varnm,subs,subc,subt)

% load
fnm_metric=[output_dir,'/','metric_px0_pz0.nc'];
    
if ~ exist(fnm_metric,'file')
   error([mfilename ': file ' fnm_metric 'does not exist']);
end

xs = subs(1); 
zs = subs(2); 

xzc = ncreadatt(fnm_metric,'/','count_of_physical_points');
xzc = double(xzc);

if(subc(1) == -1)
  xc = ceil((xzc(1)-subs(1)+1)/subt(1));
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = ceil((xzc(2)-subs(2)+1)/subt(2));
else
  zc = subc(2);
end
%stride
xt = subt(1);
zt = subt(2);

i1 = 1;
i2 = i1 + xc - 1;
k1 = 1;
k2 = k1 + zc - 1;

v(i1:i2,k1:k2)=ncread(fnm_metric,varnm,[xs,zs],[xc,zc],[xt,zt]);

end
