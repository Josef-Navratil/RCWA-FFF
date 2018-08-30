function [RP,RS,s0V] = computeScatMatNVM (lam,thI,epsB,Lam,d,epsW,epsS,fsin,fcos,nMax,N)
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate


nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sin(thI);
I=eye(nDim);


%% Relief profile

d1 = d/N;
zV = linspace(d1/2,d-d1/2,N);
wV = (1/pi)*acos(1-2*zV/d); % vector of fill factors ... sinusoidal profile
%wV = 0.5*ones(1,N);


nV=-nMax:nMax;

%% Prepare Factorization Matrices
[cosM,sinM] = generateSinCosMat(fcos,fsin, Lam, nDim);

    qV  = q0+nV*q;
    s0V = sqrt(1-qV.^2);
    sSubV = sqrt(epsS-qV.^2);
    Y0  = diag(s0V);
 
%% FFF Superstrate medium
    epsMn = epsMatrix(epsW,epsB,0,nMax);  % "n"="next"
    etaMn = epsMatrix(1/epsW,1/epsB,0,nMax);
    kX0   = diag(qV);
    [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);

%% Assemble the main matrix
    CP0 = [-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
%% Generate propagation matrices in superstrate, p-polarization
    [Z_p0,Z_n0,~,~] = generatePRPMat(CP0,nDim,k0,0);

%% 1st layer
    epsMn = epsMatrix(epsW,epsB,wV(1),nMax);  % "n"="next"
    etaMn = epsMatrix(1/epsW,1/epsB,wV(1),nMax);
    CSn   = epsMn - diag(qV.^2);
    [GSn,muSn] = eig(CSn);
    sVSn  = csqrt(diag(muSn));
    YSn   = GSn*diag(sVSn)/GSn;
   
    [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    % Assemble the main matrix
    CP = [-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
    % Generate propagation, admittance and reflection matrices
    % s-polarization, interface 0/1
    [Z_pii,Z_nii,PP_p,PP_n]   = generatePRPMat(CP,nDim,k0,d1);
    [RP_01,TP_01,RP_10,TP_10] = generateREFMat(Z_pii,Z_p0,Z_nii,Z_n0,nDim);
    % Mark as old reflection matrices
    RP   = RP_01;
    TP   = TP_01;
    RP_i = RP_10;
    TP_i = TP_10;
    % Generate reflection matrices s-polarization, interface 0/1
    [RSn,TSn,RSin,TSin] = scatInterface(Y0,YSn);



    
    for iL=1:N
        RS  = RSn; 
        TS  = TSn;
        RSi = RSin;
        TSi = TSin; % interface 0/iL
        RP_old  = RP;
        RPi_old = RP_i;
        TP_old  = TP;
        TPi_old = TP_i;
        YS    = YSn;
        Z_pi  = Z_pii;
        Z_ni  = Z_nii;
        GS    = GSn;
        sVS   = sVSn;
        PS    = GS*diag(exp(1i*k0*sVS*d1))/GS;
        PP_pi = PP_p;
        PP_ni = PP_n;        

        if iL<N
            epsMn = epsMatrix(epsW,epsB,wV(iL+1),nMax);
            etaMn = epsMatrix(1/epsW,1/epsB,wV(iL+1),nMax);
            
            % Generate admittance matrices for p-polarization
            [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
            CP = [-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
            [Z_pii,Z_nii,PP_p,PP_n] = generatePRPMat(CP,nDim,k0,d1);
            
            % Generate impendance matrices for s-polarization
            CSn  = epsMn - diag(qV.^2);
            [GSn,muSn] = eig(CSn);
            sVSn = csqrt(diag(muSn));
            YSn  = GSn*diag(sVSn)/GSn;
 
        else % Last layer
            YSn   = diag(sSubV);
            Z_pii = diag(sSubV)/epsS;
            Z_nii = -Z_pii;
        end
   
%% Scattering Layer s-polarization
        [RSn1,TSn1,RSin1,TSin1] = scatInterface(YS,YSn);
        [RSn,TSn,RSin,TSin]     = scatLayer(RS,TS,RSi,TSi, RSn1,TSn1,RSin1,TSin1, PS);

%% Scattering Layer p-polarization
        [RP_ijj,TP_ijj,RP_jji,TP_jji] = generateREFMat(Z_pii,Z_pi,Z_nii,Z_ni,nDim);
        [RP,TP,RP_i,TP_i]             = generateSCTMat(RPi_old,TP_old,RP_old,TPi_old,PP_ni,PP_pi,RP_ijj,RP_jji,TP_ijj,TP_jji,nDim,iL);

    end
end

