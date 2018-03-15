% planar diffraction by a sinusoidal relif grating
clear all
lam = 300e-9;
%permW=1.0;
%permW = 1.5^2;
%permW = 3^2;
permW=(1+5i)^2;
thI = 15*(pi/180);
epsB = 1;
Lam = 600e-9;
d=600e-9;

nMax = 100;
nDim = 2*nMax+1;

% relief profile
N = 400; % number of slices
d1 = d/N;
zV = linspace(d1/2,d-d1/2,N);
wV = (1/pi)*acos(1-2*zV/d); % vector of fill factors ... sinusoidal profile
%wV = 0.5*ones(1,N);

q0=sin(thI);
nV=-nMax:nMax;
I=eye(nDim);

%% Priprava FFF
K1=2*pi/Lam;
%fcos=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %function cos(phi(x)) of tangent component
%fsin=@(x) (((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %function sin(phi(x)) of tangent component
%fcos=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce cos(phi(x))
%fsin=@(x) -(((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce sin(phi(x))
fsin=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce cos(phi(x))
fcos=@(x) (((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce sin(phi(x))
%fcos=@(x) 1+0.*x;
%fsin=@(x) 0.*x;
cosM=F_series_gen(fcos,12,Lam,nDim);
cosM=toeplitz(cosM);
sinM=F_series_gen(fsin,12,Lam,nDim);
sinM=toeplitz(sinM);


clear fcos fsin K1

    
    epsW=permW;
    epsS=epsW;
    k0=2*pi/lam;
    q=lam/Lam;
    qV=q0+nV*q;
    s0V = sqrt(1-qV.^2);
    sSubV = sqrt(epsS-qV.^2);
    Y0 = diag(s0V);
 
  %% FFF Superstrate medium
    epsMn = epsMatrix(epsW,epsB,0,nMax);  % "n"="next"
    etaMn = epsMatrix(1/epsW,1/epsB,0,nMax);
        A=cosM*(etaMn\cosM) + sinM*epsMn*sinM; %Auxiliary factorization matrices
        B=sinM*epsMn*cosM - cosM*(etaMn\sinM);
        C=cosM*epsMn*sinM - sinM*(etaMn\cosM);
        D=sinM*(etaMn\sinM) + cosM*epsMn*cosM;
     kX0=diag(qV);
    %Assemble the main matrix
    CP0=[-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
%% Generate propagation matrices in superstrate    
    [GP0,sV]=eig(CP0);
    sV=diag(sV);
    sVPsort=(real(sV)+imag(sV)>0);
    sVNsort=~sVPsort;
    GPE0=GP0(1:nDim,sVPsort);
    GNE0=GP0(1:nDim,sVNsort);
    GPH0=GP0(nDim+1:2*nDim,sVPsort);
    GNH0=GP0(nDim+1:2*nDim,sVNsort);
    Z_p0 = GPE0/GPH0;
    Z_n0 = GNE0/GNH0;

%% 1st layer
    epsMn = epsMatrix(epsW,epsB,wV(1),nMax);  % "n"="next"
    etaMn = epsMatrix(1/epsW,1/epsB,wV(1),nMax);
    CSn = epsMn - diag(qV.^2);
    CPn = etaMn\(I - diag(qV)*(epsMn\diag(qV)));
    [GSn,muSn] = eig(CSn);
    [GPn,muPn] = eig(CPn);
    sVPn=sqrt(diag(muPn));
    sVSn=csqrt(diag(muSn));
    YSn = GSn*diag(sVSn)/GSn;
    YPn = etaMn*(GPn*diag(sVPn)/GPn);
   
        A=cosM*(etaMn\cosM) + sinM*epsMn*sinM; %Auxiliary factorization matrices
        B=sinM*epsMn*cosM - cosM*(etaMn\sinM);
        C=cosM*epsMn*sinM - sinM*(etaMn\cosM);
        D=sinM*(etaMn\sinM) + cosM*epsMn*cosM;
    %Assemble the main matrix
    CP=[-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
    %Generate propagation matrices
    [GVec,sV]=eig(CP);
    sV=diag(sV);
    sVPsort=(real(sV)+imag(sV)>0);
    sVNsort=~sVPsort;
    sVPos=sV(sVPsort);
    sVNeg=sV(sVNsort);
    GPE=GVec(1:nDim,sVPsort);
    GNE=GVec(1:nDim,sVNsort);
    GPH=GVec(nDim+1:2*nDim,sVPsort);
    GNH=GVec(nDim+1:2*nDim,sVNsort);
    Z_pii = GPE/GPH;
    Z_nii = GNE/GNH;
    PP_p=(GPH*diag(exp(1i*k0*sVPos*d1))/GPH);
    PP_n=(GNH*diag(exp(-1i*k0*sVNeg*d1))/GNH);


    [RSn,TSn,RSin,TSin] = scatInterface(Y0,YSn); % interface 0/1
    [RPn,TPn,RPin,TPin] = scatInterface(Y0,YPn);

    RP_01 = (Z_pii-Z_n0)\(Z_p0-Z_pii);
    TP_01 = I+RP_01;
    RP_10 = (Z_n0-Z_pii)\(Z_nii-Z_n0);
    TP_10 = I+RP_10;

    %Auxiliary block preparing matrices for the cycle
%    RP_ijj=RP_01; %RP_i,i+1
%    RP_jji=RP_10; %RP_i+1,i
%    TP_ijj=TP_01; %TP_i,i+1
%    TP_jji=TP_10; %TP_i+1,i
    RP=RP_01;
    TP=TP_01;
    RP_i=RP_10;
    TP_i=TP_10;

    
    for iL=1:N
        RS=RSn; TS=TSn; RSi=RSin; TSi=TSin; % interface 0/iL
        RP_old=RP; RPi_old=RP_i;
        TP_old=TP; TPi_old=TP_i;
        YS = YSn; Z_pi = Z_pii; Z_ni = Z_nii;
        
        RPP=RPn; TPP=TPn; RPi=RPin; TPi=TPin; GP=GPn; sVP=sVPn; YP = YPn;
        PP_pi=(GPH*diag(exp(1i*k0*sVPos*d1))/GPH);
        PP_ni=(GNH*diag(exp(-1i*k0*sVNeg*d1))/GNH);        
        GS=GSn;
        sVS=sVSn;
        PS = GS*diag(exp(1i*k0*sVS*d1))/GS;
        PP = GP*diag(exp(1i*k0*sVP*d1))/GP;
        if iL<N
            epsMn = epsMatrix(epsW,epsB,wV(iL+1),nMax);
            etaMn = epsMatrix(1/epsW,1/epsB,wV(iL+1),nMax);
            CSn = epsMn - diag(qV.^2);
            CPn = etaMn\(I - diag(qV)*(epsMn\diag(qV)));
            [GSn,muSn] = eig(CSn);
            [GPn,muPn] = eig(CPn);
            sVSn=csqrt(diag(muSn));
            sVPn=sqrt(diag(muPn));
            YSn = GSn*diag(sVSn)/GSn;
            YPn = etaMn*(GPn*diag(sVPn)/GPn);
            %Generate F Matrices
                A=cosM*(etaMn\cosM) + sinM*epsMn*sinM; %Auxiliary factorization matrices
                B=sinM*epsMn*cosM - cosM*(etaMn\sinM);
                C=cosM*epsMn*sinM - sinM*(etaMn\cosM);
                D=sinM*(etaMn\sinM) + cosM*epsMn*cosM;
            %Assemble the main matrix
            CP=[-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];
            %Generate propagation matrices
            [GVec,sV]=eig(CP);
                sV=diag(sV);
                sVPsort=(real(sV)+imag(sV)>0);
                sVNsort=~sVPsort;
                sVPos=sV(sVPsort);
                sVNeg=sV(sVNsort);
                GPE=GVec(1:nDim,sVPsort);
                GNE=GVec(1:nDim,sVNsort);
                GPH=GVec(nDim+1:2*nDim,sVPsort);
                GNH=GVec(nDim+1:2*nDim,sVNsort);
            Z_pii = GPE/GPH;
            Z_nii = GNE/GNH;
        else
            YSn=diag(sSubV);
            Z_pii=diag(sSubV)/epsS;
            Z_nii=-Z_pii;
            YPn=diag(sSubV)/epsS;
        end
        [RSn1,TSn1,RSin1,TSin1] = scatInterface(YS,YSn);
        [RPn1,TPn1,RPin1,TPin1] = scatInterface(YP,YPn);
        RP_ijj = (Z_pii-Z_ni)\(Z_pi-Z_pii);
        TP_ijj = I+RP_ijj;
        RP_jji = (Z_ni-Z_pii)\(Z_nii-Z_ni);
        TP_jji = I+RP_jji;

%        [RPn1,TPn1,RPin1,TPin1] = scatInterfaceP(YPPos,YPNeg,YPPosN,YPNegN);

        [RSn,TSn,RSin,TSin] = scatLayer(RS,TS,RSi,TSi, RSn1,TSn1,RSin1,TSin1, PS);
        [RPn,TPn,RPin,TPin] = scatLayer(RPP,TPP,RPi,TPi, RPn1,TPn1,RPin1,TPin1, PP);
%% Scattering Layer p-polarization
        QP = RPi_old*PP_ni*RP_ijj*PP_pi;
%        QP = RPi_old*PP_pi*RP_ijj*PP_ni;
        tmp = (I-QP)\TP_old;
        disp([iL cond(I-QP)])
        RP = RP_old + TPi_old*PP_ni*RP_ijj*PP_pi*tmp;
%        RP = RP_old + TPi_old*PP_pi*RP_ijj*PP_ni*tmp;
        TP = TP_ijj*PP_pi*tmp;

        QN = RP_ijj*PP_pi*RPi_old*PP_ni;
%        QN = RP_ijj*PP_ni*RPi_old*PP_pi;
        tmp = (I-QN)\TP_jji;
        RP_i = RP_jji + TP_ijj*PP_pi*RPi_old*PP_ni*tmp;
%        RP_i = RP_jji + TP_ijj*PP_ni*RPi_old*PP_pi*tmp;
        TP_i = TPi_old*PP_ni*tmp;
                

    end
    
    % RESULTS
    RPvec=abs(diag(RP).^2);
    RS=abs(diag(RSn).^2);
    RS02 = abs(RSn(nMax+1,nMax+1))^2;
    RP02 = abs(RP(nMax+1,nMax+1))^2;
    RPP02 = abs(RPn(nMax+1,nMax+1))^2;
    TS02 = (abs(TSn(nMax+1,nMax+1))^2)*sSubV(nMax+1)/s0V(nMax+1);
    TP02 = abs(TP(nMax+1,nMax+1))^2./sqrt(epsS);
    TPP02 = abs(TPn(nMax+1,nMax+1))^2./sqrt(epsS);