function s0str =  s0string(s0)

s0factor = s0/(2*pi*1e9);
if isreal(s0)
	s0str = sprintf('2%s10^9(%g)','\pi',s0factor);
else
	s0str = sprintf('2%s10^9(%g + %gi)','\pi',real(s0factor), imag(s0factor));
end
