function hybrid = combineHRIR(head, hrir, delay_t, fs)

delay_t = delay_t + 20 / fs; % add a 20 sample delay to all results

l = length(hrir);

Hc = fft(hrir, 2^12);
Hh = fft(head, 2^12);
Ac = (abs(Hc));
Ah = (abs(Hh));
L = length(Ah);
nfreq = floor(L/2);      % DC is not counted
omega = ((0:nfreq)*(fs/(2*nfreq)))';
omega_low = 250;
omega_high = 1000;

Ac_cut = Ac(1:nfreq+1);
Ah_cut = Ah(1:nfreq+1);

[~, idx_low] = min(abs(omega-omega_low));
[~, idx_high] =  min(abs(omega-omega_high));

% construct desired log magnitude
As = [Ah_cut(1:idx_low); zeros(idx_high-idx_low-1, 1); Ac_cut(idx_high:end)];
for i = idx_low+1:idx_high-1
    As(i) = Ah_cut(i) + (Ac_cut(i)-Ah_cut(i))/(omega(idx_high)-omega(idx_low));
end
         
As = [As; flipud(conj(As(2:end-1,:)))];

% minimum phase reconstruction   [0:L-1].*fs/L
As_mps = mps(As);
As_mps_phase = angle(As_mps);
f = ((0:length(As_mps)-1)*(fs/(length(As_mps))))';
tau_phi = -2 * pi * delay_t * f;
As_mps_phase = As_mps_phase + tau_phi;
Hs = abs(As_mps).*exp(1i*As_mps_phase);

hybrid = ifft(Hs, 'symmetric');

hybrid = hybrid(1:l);

end

