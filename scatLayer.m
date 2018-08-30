function [R02,T02,R20,T20] = scatLayer(R01,T01,R10,T10, R12,T12,R21,T21, P)
%Scattering by a layer between two interfaces
%

I = eye(size(R01));

Q = R10*P*R12*P;
tmp = (I-Q)\T01;
R02 = R01 + T10*P*R12*P*tmp;
T02 = T12*P*tmp;

Q = R12*P*R10*P;
tmp = (I-Q)\T21;
R20 = R21 + T12*P*R10*P*tmp;
T20 = T10*P*tmp;

end
