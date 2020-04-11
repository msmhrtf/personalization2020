function out = torsoShadow(in, thetaT, thetaHS, thetaHT, a, b, alfa_min, theta_min, frac, fs)
% torso shadow
% 3 filters in cascade

out = head_shadow(in, thetaT, b, alfa_min, theta_min, fs); % torso shadow
out = head_shadow(out, thetaHS, a, alfa_min, theta_min, fs); % head shadow
out = head_itd(out, thetaHT, a, fs, frac); % head ITD

end

