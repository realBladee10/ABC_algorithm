function [modeshape,k_term] = TangentialMode(n1,n2,l1,l2,rA,rB,r0A,r0B,k)
%-----------------------------------------------------
% Compute the Tangential Mode in a room.
% Syntax: [modeshape,k_term] = 
%          TangentialMode(n1,n2,l1,l2,rA,rB,r0A,r0B,k)
%-----------------------------------------------------
global S V c  bta


tau_m = (3*V)/(5*c*S*bta); % time constant of mth mode
An = sqrt(4); %normalized constant

% Initialisation of variables

% Calculate tangential modes shapes and wavenumber of modes

    modeshape = An.*cos(n1.*pi*rA/l1).*cos(n2'*pi*rB/l2).*...
                     An.*cos(n1.*pi*r0A/l1).*cos(n2'*pi*r0B/l2);
    k_m = sqrt((n1.*pi/l1).^2 + (n2'*pi/l2).^2);

% Denominator term for Green function
k_term = k^2 - k_m.^2 - 1i*k/(tau_m*c);
