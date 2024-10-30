function [v,t]=gather_snap(parfnm,output_dir,nlayer,varnm,subs,subc,subt)

% load
fnm_snap=[output_dir,'/','volume_vel.nc'];

% check dir exists
if ~ exist(fnm_snap,'file')
    error([mfilename ': file ' fnm_snap ' does not exist']);
end

info_snap=ncinfo(fnm_snap);
tdim=info_snap.Dimensions(1);
if tdim.Length==0 || (nlayer-1)-1>=tdim.Length
   error([num2str(nlayer) 'th layer is beyond current time dim (' ...
        num2str(tdim.Length) ') in ' fnm_snap]);
end

% read parameters file
jsontext=fileread(parfnm);
par=jsondecode(jsontext);
snap_subc=par.snapshot(1).grid_index_count;
snap_subc = double(snap_subc);

xs = subs(1); 
zs = subs(2); 
if(subc(1) == -1)
  xc = ceil((snap_subc(1)-subs(1)+1)/subt(1));
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = ceil((snap_subc(2)-subs(2)+1)/subt(2));
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
% get data
v(i1:i2,k1:k2)=ncread(fnm_snap,varnm, ...
      [xs,zs,nlayer],[xc,zc,1],[xt,zt,1]);
t=ncread(fnm_snap,'time',[nlayer],[1]);

%v=v';

end
