function [cosM,sinM] = generateSinCosMat (fcos, fsin, Lam, nDim)

%% Inputs are functions of tangent and normal to profile, outputs are Topelitz matrices generated from its fourier coefficients 
cosM = F_series_gen(fcos,12,Lam,nDim); %generates vector of fourier coefficients of cos(angle to normal)
cosM = toeplitz(cosM);
sinM = F_series_gen(fsin,12,Lam,nDim); %generates vector of fourier coefficients of cos(angle to normal)
sinM = toeplitz(sinM);
end
