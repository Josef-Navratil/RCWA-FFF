function epsM = epsMatrix(epsW, epsB, w, nMax)
%Toeplitz matrix of permittivity coefficients
%   of a symmetric grating with w=fillfactor

nDim = 2*nMax+1;
kV = 1:nMax;

colV = zeros(nDim,1);
colV(1,1) = w*epsW + (1-w)*epsB;
colV(2:(nMax+1),1) = (1i./(2*pi*kV.')).*(exp(-1i*pi*w*kV.')-exp(1i*pi*w*kV.'))*(epsW-epsB);

rowV = zeros(1,nDim);
rowV(1,1) = w*epsW + (1-w)*epsB;
rowV(1,2:(nMax+1)) = (-1i./(2*pi*kV)).*(exp(1i*pi*w*kV)-exp(-1i*pi*w*kV))*(epsW-epsB);

epsM = toeplitz(colV,rowV);

end
