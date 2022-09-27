function [v] = gather_media(parfnm,output_dir,varnm)

% load
fnm_media=[output_dir,'/','media.nc'];
    
if ~ exist(fnm_media,'file')
   error([mfilename ': file ' fnm_media 'does not exist']);
end

xzs = nc_attget(fnm_media,nc_global,'local_index_of_first_physical_points');
xs = double(xzs(1)); % = 3
zs = double(xzs(2)); % = 3

xzc = nc_attget(fnm_media,nc_global,'count_of_physical_points');
xc = double(xzc(1));
zc = double(xzc(2));

i1 = 1;
i2 = i1 + xc - 1;
k1 = 1;
k2 = k1 + zc - 1;

v(k1:k2,i1:i2)=nc_varget(fnm_media,varnm,[zs,xs],[zc,xc],[1,1]);

v=v';

end
