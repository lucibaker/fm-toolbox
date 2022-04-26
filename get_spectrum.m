function [E, f] = get_spectrum(q_fluct, fs)
% get power spectral density E over frequency range f from fluctuating
% profile quantity q and sampling frequency fs using Welch (1967) method
% (with 50% overlap)

% interpolate NaNs
q_fluct = naninterp(q_fluct);

% signal divisions
% length
if length(q_fluct) > 5*2048 && length(q_fluct) < 10*2048
    Ldiv = 2048;
elseif length(q_fluct) >= 10*2048
    Ldiv = 8192;
else
    Ldiv = 2^nextpow2(length(q_fluct)/5);
end

% number of divisions 
Ndiv = floor(2*size(q_fluct,1)/Ldiv) - 1;

if ~mod(Ldiv,2)
    P1 = zeros(Ldiv/2+1,size(q_fluct,2),Ndiv);
else
    P1 = zeros((Ldiv+1)/2,size(q_fluct,2),Ndiv);
end

% window and compute fft over each division
hannwin = hann(Ldiv);
for j = 1:Ndiv
    idx_start = (j-1)*Ldiv/2 + 1;
    idx_end = idx_start + Ldiv - 1;
    
    Y = fft(hannwin.*q_fluct( idx_start:idx_end, :));
    
    P2 = (Y/Ldiv).*conj(Y/Ldiv); 
    if ~mod(Ldiv,2)
        P1(:,:,j) = P2(1:Ldiv/2+1,:);
    else
        P1(:,:,j) = P2(1:(Ldiv+1)/2,:);
    end
end

% double all values except that at f=0 to convert 2-sided to 1-sided spectrum
P1(2:end-1,:,:) = 2*P1(2:end-1,:,:);

% frequency range
f = fs*(0:(Ldiv/2))/Ldiv;

% average over divisions and profile
E = mean(mean(P1,3),2);

end