% DEMO for loading PIV data
load('input.mat',pfn,pxtom,pairsep);
h = pivGetHeader([pfn '.def.piv']); % load PIV file

tcut = floor(h.nt/2); % split fluid data at tcut and allocate separately
status = cat(3,uint8(zeros(h.ny, h.nx, tcut)),uint8(zeros(h.ny, h.nx, h.nt-tcut)));
u = cat(3,single(zeros(h.ny, h.nx,tcut)),single(zeros(h.ny, h.nx, h.nt-tcut)));
v = cat(3,single(zeros(h.ny, h.nx, tcut)),single(zeros(h.ny, h.nx, h.nt-tcut)));

for i = 0:(h.nt-1)
  data = pivGetFrame(h,i);
  
  status(:,:,i+1) = uint8(data{1});
  u(:,:,i+1) = single(data{2});
  v(:,:,i+1) = single(data{3});
end
clear data

% flip status, u, v 
status = flipud(status);
u = flipud(u);
v = -flipud(v);

% Convert u, v to physical units
u = u*pxtom/pairsep;
v = v*pxtom/pairsep;

% Remove rejected vectors (vectors for which status =/= 1)
u(status ~= 1) = nan;
v(status ~= 1) = nan;

% convert arrays into correct form for clustering script
U = permute(u,[3 2 1]); clear u
V = permute(v,[3 2 1]); clear v
status = permute(status,[3 2 1]);

% save fluid data
save(['fluid_vel_' datasetid '.mat'],'U','V','status');
