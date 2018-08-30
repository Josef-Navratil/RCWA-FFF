function [R,T,Ri,Ti] = scatInterface(Y0,Y1)
%
%   Reflection and transmission matrices between two interfaces
    
    I = eye(size(Y0));
    R = (Y0+Y1)\(Y0-Y1);
    T = I + R;
    
    Ri = -R;
    Ti= I + Ri;

end

