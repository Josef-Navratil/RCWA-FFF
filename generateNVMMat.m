%% Inputs are toeplitz matrices of sin and cos functions to profile and 
% toeplitz matrices of permittivities, 
%% Outputs are auxiliary matrices A,B,C,D of NVM
function [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn)
        A=cosM*(etaMn\cosM) + sinM*epsMn*sinM; %Auxiliary factorization matrices
        B=sinM*epsMn*cosM - cosM*(etaMn\sinM);
        C=cosM*epsMn*sinM - sinM*(etaMn\cosM);
        D=sinM*(etaMn\sinM) + cosM*epsMn*cosM;
end

