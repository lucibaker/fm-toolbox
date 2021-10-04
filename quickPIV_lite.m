% Quickly load & parse a .piv file
clear
% close all
addpath('Z:\luci\MATLAB\pivio-matlab')

% Load PIV header
load input.mat
% pfn = 'run_a.corr.mask.0024';
% pxtom = 1;
% pairsep = 1;
h = pivGetHeader([pfn '.def.piv']);
fprintf(['\n' pfn '\n'])

tlim =  4000; %h.nt; % limit number of frames used to calculate statistics if memory is an issue


%% Load velocity and status data from piv file

% status
status = uint8(zeros(h.ny, h.nx, h.nt));
u = single(zeros(h.ny, h.nx, h.nt));
v = single(zeros(h.ny, h.nx, h.nt));
parfor i = 0:(h.nt-1)
    data = pivGetFrame(h,i);  
    status(:,:,i+1) = uint8(data{1});
    u(:,:,i+1) = single(data{2});
    v(:,:,i+1) = single(data{3});
end
% flip status, u, v 
status = flipud(status);
u = flipud(u);
v = -flipud(v);
u = u*pxtom/pairsep;
v = v*pxtom/pairsep;

save('fluidvel.mat','u','v','status','-v7.3')


%%
% x,y coordinate
y = pxtom * (h.y0 : h.dy : (h.y0 + h.dy*(h.ny-1)))'-mean([y0_left,y0_right])*pxtom;
x = pxtom * (h.x0 : h.dx : (h.x0 + h.dx*(h.nx-1)))';


%% Data quality
% Rejected PIV vectors
u(status ~= 1) = nan;
v(status ~= 1) = nan;

%% Rejected vectors
rej_vec = status~=1 & status~=5;  % true if not rejected or masked
rej = sum(rej_vec(:))/numel(status); % percent rejected vectors [%]
fprintf([num2str(rej*100) '%% nonvalid vectors\n'])
rej_map = sum(rej_vec,3)/h.nt;
figure; pcolor(x*1000,y(y>0)*1000,rej_map(y>0,:)); colormap hot; shading flat;
xlabel('x [mm]'); ylabel('y [mm]'); c = colorbar; c.Limits = [0 1]; ylabel(c,'Rejected vector fraction'); axis equal tight

clear rej_vec

%% Displacement PDFs 
nbins = 200;
figure; subplot(211); histogram(u(:,:,1:tlim)/(pxtom/pairsep),nbins,'normalization','pdf');%,'BinLimits',[11 23]); 
xlabel('U [px]'); ylabel('PDF'); title('Displacement PDF');
subplot(212); hu = histogram(u(:,:,1:tlim)/(pxtom/pairsep) - round(u(:,:,1:tlim)/(pxtom/pairsep)),nbins,'normalization','pdf'); 
xlabel('U [px]'); ylabel('PDF'); title('Fractional U PDF'); xlim([-.5 .5])
%goodplot([5 7])

Cu = 1 - min(hu.Values)/max(hu.Values);
fprintf(['C_u = ' num2str(Cu) '\n'])

%% Mean velocity profile
ut = 0.0183; del_nu = nu/ut;
u_mean = nanmean(nanmean(u,2),3);
v_mean = nanmean(nanmean(v,2),3);
U0 = max(u_mean);

figure; %subplot(121);
semilogy(u_mean/ut,y/del_nu,'.'); xlabel('u^+'); ylabel('y^+')
k = 0.41; B = 5.5; rhs = 1/k*log(y/del_nu)+B; % assume ideal log law behavior
hold on; semilogy(rhs,y/del_nu,'k--'); axis([10 25 1e1 1e3])

% subplot(122); plot(vmean,y); xlabel('v [m/s]'); 

% save('input.mat','ut','-append')

%% Reynolds stresses
u_fluct = u(:,:,1:tlim) - repmat(u_mean,1,length(x),tlim);
u_var = rms(squeeze(rms(permute(u_fluct,[3 2 1]),'omitnan')),'omitnan').^2; 
v_fluct = v(:,:,1:tlim) - repmat(v_mean,1,length(x),tlim);
v_var = rms(squeeze(rms(permute(v_fluct,[3 2 1]),'omitnan')),'omitnan').^2; 
uv = nanmean(nanmean(u_fluct.*v_fluct,3),2);
%%
figure; %subplot(121);
semilogy(u_var/ut^2,y/del_nu,'.',v_var/ut^2,y/del_nu,'.',-uv/ut^2,y/del_nu,'.'); xlabel('u''^+^2'); ylabel('y^+')
% subplot(122); plot(vmean,y); xlabel('v [m/s]'); 
legend('u''^2','v''^2','u''v''^2')
