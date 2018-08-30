function [RP,TP,RP_i,TP_i] = generateSCTMat(RPi_old,TP_old,RP_old,TPi_old,PP_ni,PP_pi,RP_ijj,RP_jji,TP_ijj,TP_jji,nDim, iL)

%% Inputs are reflection and transmission matrices between two layers, propagation matrices, size of dimension and number of interation
%% Outputs are reflection and transmission matrices for one layer below
    I   = eye(nDim);
    QP  = RPi_old*PP_ni*RP_ijj*PP_pi;
    tmp = (I-QP)\TP_old;
    RP  = RP_old + TPi_old*PP_ni*RP_ijj*PP_pi*tmp;
    TP  = TP_ijj*PP_pi*tmp;
    disp([iL cond(I-QP)]) % Show number of interation and condition number of QP matrices - for testing purposes

    QN   = RP_ijj*PP_pi*RPi_old*PP_ni;
    tmp  = (I-QN)\TP_jji;
    RP_i = RP_jji + TP_ijj*PP_pi*RPi_old*PP_ni*tmp;
    TP_i = TPi_old*PP_ni*tmp;

end

