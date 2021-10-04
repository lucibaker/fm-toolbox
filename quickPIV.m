% Analyze water channel PIV data
% Outputs: Data quality, mean & RMS velocity profiles
% function quickPIV
clear
close all
addpath('Z:\luci\MATLAB\pivio-matlab')

% Experiment conditions
load('input.mat')
% dt = 1/30; %1200e-6; % s
% pairsep = dt;
% lm = 1; lp = 584.30; % calibration distance, m and px
% pxtom = lm/lp; %0.0583e-3; % m/pixel
% % dt = 1; pxtom = 1; pairsep = 1;
% H = 1; %0.15; % channel height, m
% % nu = 1e-6; % viscosity, m^2/s

% Load PIV header
% pfn = 'run-b.bsub.0064';
h = pivGetHeader([pfn '.def.piv']);
fprintf(['\n' pfn '\n'])

% h.nt = 3000;
tlim = 4000; % h.nt; %  limit number of frames used to calculate statistics if memory is an issue
tcut = floor(h.nt/4); % split fluid data at tcut and allocate separately


%% Load velocity and status data from piv file

% status
status = cat(3,uint8(zeros(h.ny, h.nx, tcut)),uint8(zeros(h.ny, h.nx, tcut)), ...
    uint8(zeros(h.ny, h.nx, tcut)),uint8(zeros(h.ny, h.nx, h.nt-3*tcut)));
parfor i = 0:(h.nt-1)
    try
  data = pivGetFrame(h,i);  
  status(:,:,i+1) = uint8(data{1});
    catch
    end
end
% flip status, u, v 
status = flipud(status);
save('fluidvel.mat','status','-v7.3');
clear status

%% u
u = cat(3,single(zeros(h.ny, h.nx,tcut)),single(zeros(h.ny, h.nx,tcut)), ...
    single(zeros(h.ny, h.nx,tcut)),single(zeros(h.ny, h.nx, h.nt-3*tcut)));
parfor i = 0:(h.nt-1)
    try
  data = pivGetFrame(h,i);
  u(:,:,i+1) = single(data{2});
    catch
    end
end
% flip status, u, v 
u = flipud(u);
u = u*pxtom/pairsep;
save('fluidvel.mat','u','-v7.3','-append');
clear u

% v
v = cat(3,single(zeros(h.ny, h.nx,tcut)),single(zeros(h.ny, h.nx,tcut)), ...
    single(zeros(h.ny, h.nx,tcut)),single(zeros(h.ny, h.nx, h.nt-3*tcut)));
parfor i = 0:(h.nt-1)
    try
  data = pivGetFrame(h,i);
  v(:,:,i+1) = single(data{3});
    catch
    end
end
% flip status, u, v 
v = -flipud(v);
v = v*pxtom/pairsep;
save('fluidvel.mat','v','-v7.3','-append');
clear v

% %% write to text file
% f = fopen([pfn '.txt'],'wt');
% fprintf(f,'PIV results - zoomed (invalid vectors are replaced by NaN)\n');
% fprintf(f,'frame# x[px] y[px] u[px/frame] v[px/frame]\n\n');
% 
% x_px = h.x0:h.dx:(h.nx*h.dx);
% y_px = h.y0:h.dy:(h.ny*h.dy);
% for k = 1:tlim
%     for i = 1:h.nx
%         for j = 1:h.ny
%             fprintf(f,'%4u %4u %4u %1.6f %1.6f\n',k,x_px(i),y_px(j),u(j,i,k),v(j,i,k));
%         end
%     end
% end
% 
% fclose(f);

% keyboard

%% reload data
load('fluidvel.mat','u');
u = u(:,:,1:tlim);
% load('fluidvel.mat','v');
load('fluidvel.mat','status');
status = status(:,:,1:tlim);

%%
% x,y coordinate
img = imgRead([ifn '.corr.mask.img']);
[~,y0] = max(mean(flipud(imgGetFrame(img,1000)),2));  % floor height
y0_piv = max([round((y0-h.y0)/h.dy)+1, 1]); 

x = pxtom*(h.x0 : h.dx : (h.x0 + h.dx*(h.nx-1)));
y = pxtom*((h.y0 : h.dy : (h.y0 + h.dy*(h.ny-1))) - y0);

% % mask
% if exist('mask.bmp','file')
%     mask = flipud(imread('mask.bmp'));
%     mask_piv = mask(h.y0:h.dy:(h.iy-1),h.x0:h.dx:(h.ix-1));
%     status(repmat(mask_piv,1,1,size(status,3)) == 0) = 5;
%     save('fluidvel.mat','u','v','status','-v7.3');
% end

%% Data quality
% Rejected PIV vectors
u(status ~= 1) = nan;
% v(status ~= 1) = nan;

%% Rejected vectors
rej_vec = status(:,:,1:tlim)~=1 & status(:,:,1:tlim)~=5;  % true if not rejected or masked
rej = sum(rej_vec(:))/numel(status(:,:,1:tlim)); % percent rejected vectors [%]
fprintf([num2str(rej*100) '%% nonvalid vectors\n'])
rej_map = sum(rej_vec,3)/tlim;
figure; pcolor(x*1000,y*1000,rej_map); colormap hot; shading flat;
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

%% Velocity profiles
% Mean velocity profiles
u_mean = nanmean(nanmean(u(:,:,1:tlim),3),2);
% v_mean = nanmean(nanmean(v(:,:,1:tlim),3),2);

% RMS velocity profiles
u_rms = rms(squeeze(rms(permute(u(:,:,1:tlim)-u_mean,[3 2 1]),'omitnan')),'omitnan'); 
% v_rms = rms(squeeze(rms(permute(v(:,:,1:tlim)-v_mean,[3 2 1]),'omitnan')),'omitnan');

% Freestream velocity 
U0 = 1; %max(umean);

% plot mean velocity
figure; subplot(121); 
plot(u_mean/U0,y/H,'b^'); 
xlabel('\langleu\rangle'); ylabel('y/H'); title('Mean streamwise')
% subplot(122); plot(v_mean/U0,y/H,'b^'); %goodplot([7 6])
% xlabel('\langlev\rangle'); ylabel('y/H'); title('Mean wall-normal')

% plot RMS velocity
% figure; 
subplot(122); plot(u_rms/U0,y/H,'b^');
xlabel('u_{rms}'); ylabel('y/H'); title('RMS streamwise')

% subplot(122); plot(v_rms/U0,y/H,'b^'); %goodplot([7 6])
% xlabel('v_{rms}'); ylabel('y/H'); title('RMS wall-normal')

%% Log layer profile
ut = 0.0126; %0.0148; %0.0181; % friction velocity estimate
del_nu = nu/ut;    % viscous length 
k = .41; B = 5.5; %5.2;  % assume canonical log law coefficients
rhs = 1/k*log(y/del_nu)+B;
% rhs2 = y/del_nu;

figure; semilogx(y(y0_piv:end)/del_nu,u_mean(y0_piv:end)/ut,'k.',y/del_nu,rhs,'r--'); hold on
% plot(y(3:5)/del_nu,rhs2(3:5),'b--')
xlabel('y^+'); ylabel('u^+'); title('Mean velocity profile'); grid on; %goodplot([7 6]);
fprintf(['u_tau = ' num2str(ut) ' m/s\n'])

%%
keyboard
%%
save('input.mat','ut','-append')

%% plot
u_mean_field = nanmean(u(:,:,1:tlim),3);
% v_mean_field = nanmean(v(:,:,1:tlim),3);
figure; pcolor_img(u_mean_field); axis equal; colormap jet; colorbar
hold on; quiver(u_mean_field,v_mean_field,'color','k')

%% Reynolds number
U_b = trapz(y(1:end-1),[0;u_mean(2:end-1)])/(y(end-1)+pxtom*h.dy/2);
Re_b = U_b*H/nu;
fprintf('Re_b = %2.f\n',Re_b);

%% degradation over time
% N = 4;
% N_edges = round(linspace(1,h.nt,N+1));
% u_mean_time = zeros(N,length(y));
% l = cell(N,1);
% c = {'r','b','k','g'};
% figure; 
% for i = 1:N
%     u_mean_time(i,:) = nanmean(nanmean(u(:,:,N_edges(i):N_edges(i+1)),3),2);
%     plot(u_mean_time(i,:),y/H,'linewidth',1.5); hold on
%     l{i} = [num2str(N_edges(i)) ' - ' num2str(N_edges(i+1))];
% end
% ylabel('y/H');xlabel('\langleu\rangle [m/s]')
% legend(l,'location','nw')
% axis auto