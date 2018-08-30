function [RP_ijj,TP_ijj,RP_jji,TP_jji] = generateREFMat(Z_pii,Z_pi,Z_nii,Z_ni,nDim)
%% Inputs are admittance matrices and number of modes, p-polarization
%% Outputs are reflection and transmission matrices for one interface, p-polarization
        I = eye(nDim);
        RP_ijj = (Z_pii-Z_ni)\(Z_pi-Z_pii); % R_i,i+1
        TP_ijj = I+RP_ijj;                  % T_i,i+1 
        RP_jji = (Z_ni-Z_pii)\(Z_nii-Z_ni); % R_j+1,j
        TP_jji = I+RP_jji;                  % T_j+1,j

end

