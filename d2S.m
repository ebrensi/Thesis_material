function d = d2S(mu,s)
smin =  min(abs(s));
smax =  max(abs(s));
d = zeros(length(mu),1);

for i = 1:length(mu)
	imu = abs(imag(mu));  rmu = real(mu);
	if isinf(mu(i))
		d(i) = inf;
	elseif imu > smax  %
		d(i) = norm([rmu  imu-smax]);
	elseif imu < smin
		d(i) = norm([rmu  smin-imu]);
	else
		d(i) = abs(real(mu(i)));
	end
end

end
