function c = csqrt(z)
% Complex square root
%   with the result having both real and imaginary values positive

a = abs(z);
x = real(z);
c = sqrt(0.5*(x+a)) + 1i*sqrt(0.5*(-x+a));

end
