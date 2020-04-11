function [hrir_left, hrir_right] = getGenericHRTF(source, microphone, fs)

yes = 0; % plot

% load cipic kemar small pinna (subject 165)
cipic_small = load('C:\Users\jtissi16\Documents\MATLAB\TH\CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_021\hrir_final.mat');
% angle sampling
azimuths = [-80 -65 -55 -45:5:45 55 65 80];
elevations = -45 + 5.626*(0:49);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get angle position of source and listener and get hrtf from database

% get azimuth and elevation relative to source and mic positions and heading
azimuth = microphone.getAzimuthInteraural(source);
elevation = microphone.getElevationInteraural(source);

% % adopt elevation angle for cipic library [-90;270]
% if ( (rad2deg(elevation) < -90) && (elevation > -180) )
%     elevation = deg2rad(360+rad2deg(elevation));
% end

[az_err, az_idx] = min(abs(azimuths-rad2deg(azimuth)));
[el_err, el_idx] = min(abs(elevations-rad2deg(elevation)));

% get hrtf from cipic kemar small pinna database (vertical-polar)
hrir_left = squeeze(cipic_small.hrir_l(az_idx, el_idx, :));
hrir_right = squeeze(cipic_small.hrir_r(az_idx, el_idx, :));


% change sampling rate
if fs ~= 44100
    hrir_left = changeSamplingRate(hrir_left, 44100, fs);
    hrir_right = changeSamplingRate(hrir_right, 44100, fs);
end

% plots
if yes
    % compute FFTs of HRIR (CIPIC function)
    [hrtf_left, sampfreq_l] = freq_resp(hrir_left, 0, Inf); % [hrir_left; zeros(2^15,1)]
    [hrtf_right, sampfreq_r] = freq_resp(hrir_right, 0, Inf); % [hrir_right; zeros(2^15,1)]
    figure
    subplot(2,1,1); plot(hrir_left); hold on; plot(hrir_right); grid on
    subplot(2,1,2); plot(sampfreq_l, hrtf_left); hold on; plot(sampfreq_r, hrtf_right); grid on
end

end