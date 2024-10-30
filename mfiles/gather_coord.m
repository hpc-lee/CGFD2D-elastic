function [x,z] = gather_coord(output_dir,subs,subc,subt)

% load
fnm_coord=[output_dir,'/','coord_px0_pz0.nc'];
    
if ~ exist(fnm_coord,'file')
   error([mfilename ': file ' fnm_coord 'does not exist']);
end

xs = subs(1); 
zs = subs(2); 

xzc = ncreadatt(fnm_coord,'/','count_of_physical_points');
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

x(i1:i2,k1:k2)=ncread(fnm_coord,'x',[xs,zs],[xc,zc],[xt,zt]);

z(i1:i2,k1:k2)=ncread(fnm_coord,'z',[xs,zs],[xc,zc],[xt,zt]);
    
end
