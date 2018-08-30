% planar diffraction by a sinusoidal relif grating
clear all
lam    = 300e-9;
%permW = 1.0;
%permW = 1.5^2;
%permW = 3^2;
permW  = (1+5i)^2;
thI    = 15*(pi/180);
epsB   = 1;
Lam    = 600e-9;
d      = 600e-9;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

%% Functions of sin and cos tangent to the profile

K1=2*pi/Lam;

fsin=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce cos(phi(x))
fcos=@(x) (((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce sin(phi(x))

% No factorization

%fcos=@(x) 1+0.*x;
%fsin=@(x) 0.*x;

%% Preallocate fields for reflection amplitudes and errors

RPm2=zeros(1,1);
RPm1=zeros(1,1);
RP0=zeros(1,1);
RP1=zeros(1,1);
RSm2=zeros(1,1);
RSm1=zeros(1,1);
RS0=zeros(1,1);
RS1=zeros(1,1);
errorS=zeros(1,1);
errorP=zeros(1,1);


%% Reference values for sinusiodal grating, computed by C-Method
RP_ref=[0.0674606805767746 0.00319796649824055 0.303105896615171 0.124877473586179];
RS_ref=[0.327219701682201 0.207465501464766 0.117286155151613 0.0792741802547489];


%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    = 60;         % lowest number of modes
nMax_u    = 60;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    = 600;           % lower number of layers
N_u    = 600;           % upper number of layers
N_step = 1;             % length of the step - layers

tic                   % start time count

for nMax =  nMax_l:nMax_step:nMax_u
    
    for N = N_l:N_step:N_u        
        [RP,RS,s0V] = computeScatMatNVM (lam,thI,epsB,Lam,d,epsW,epsS,fsin,fcos,nMax,N);  
    
        % RESULTS
        RPvec = (abs((RP(nMax-1:nMax+2,nMax+1)').^2).*s0V(1,nMax-1:nMax+2))./s0V(1,nMax+1);
        RSvec = abs((RS(nMax-1:nMax+2,nMax+1)').^2).*s0V(1,nMax-1:nMax+2)./s0V(1,nMax+1);
        RPm2((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = RPvec(1);
        RPm1((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = RPvec(2);
        RP0((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)  = RPvec(3);
        RP1((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)  = RPvec(4);
        RSm2((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = RSvec(1);
        RSm1((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = RSvec(2);
        RS0((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)  = RSvec(3);
        RS1((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)  = RSvec(4);
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
        errorS((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = (1/4)*(((RSvec(1)-RS_ref(1))/RS_ref(1)).^2+((RSvec(2)-RS_ref(2))/RS_ref(2)).^2+((RSvec(3)-RS_ref(3))/RS_ref(3)).^2+((RSvec(4)-RS_ref(4))/RS_ref(4)).^2);
        errorP((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = (1/4)*(((RPvec(1)-RP_ref(1))/RP_ref(1)).^2+((RPvec(2)-RP_ref(2))/RP_ref(2)).^2+((RPvec(3)-RP_ref(3))/RP_ref(3)).^2+((RPvec(4)-RP_ref(4))/RP_ref(4)).^2);
    end
end

%% Save results to a file
%filename = 'test.mat';
%save(filename)
%exit
  