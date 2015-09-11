function d = chord(z1,z2)
% d = chord(z1,z2)
%  chordal metric for complex numbers.  z1, z2 can be vectors or scalars
lz1 = length(z1); lz2 = length(z2);
d = zeros(max(lz1,lz2),1);
if lz1 == 1
    z1 = repmat(z1,lz2,1);
elseif lz2 == 1
    z2 = repmat(z2,lz1,1);
end
inf1 = isinf(z1);
d(inf1) = 2 ./ s(z2(inf1));

inf2 = isinf(z2) & ~inf1;
d(inf2) = 2 ./ s(z1(inf2));

fin = ~(inf1 | inf2);
d(fin) = 2*abs(z1(fin)-z2(fin)) ./ (s(z1(fin)) .* s(z2(fin)));

    function y = s(x)
        y = sqrt(1 + abs(x).^2);
    end
end
