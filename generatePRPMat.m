function [Z_p0,Z_n0,PP_p,PP_n] = generatePRPMat(CP,nDim,k0,d1)
%% Inputs are CP matrix, truncation number, wavenumber and layer thickness
%% Outputs are admittance and propagation matrices

    [GP,sV] = eig(CP);
    sV      = diag(sV);
    sVPsort = (real(sV)+imag(sV)>0);
    sVNsort = ~sVPsort;
    sVPos   = sV(sVPsort);
    sVNeg   = sV(sVNsort);
    GPE     = GP(1:nDim,sVPsort);
    GNE     = GP(1:nDim,sVNsort);
    GPH     = GP(nDim+1:2*nDim,sVPsort);
    GNH     = GP(nDim+1:2*nDim,sVNsort);
    Z_p0    = GPE/GPH;
    Z_n0    = GNE/GNH;
    PP_p    = GPH*diag(exp(1i*k0*sVPos*d1))/GPH;
    PP_n    = GNH*diag(exp(-1i*k0*sVNeg*d1))/GNH;
    
end