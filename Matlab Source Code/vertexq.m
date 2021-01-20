%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% Function: vertexq
% 
% Generates the momentum dependence of the electron-phonon coupling
%
% vertex |g(q)|^2 = g_0^2*f(q).
%   f(q) = vq. The normalization is
%   [1/N_q]\sum_{q_x,q_y} f(q) = 1 in discrete form
%   or [1/(2\pi)^2]\int_{q\in BZ} dq_x dq_y f(q) = 1 in continuous form.
%   funtype: 'exp', 'gauss' for forward scattering coupling
%   with a peak at q=0, and an optional width q0.
%
%   For cuprates, see Devereaux/Virosztek/Zawadowski PRB1995;
%   Bulut/Scalapino PRB1996; Sandvik/Scalapino/Bickers PRB2004.
%             O(3)
%             |
%   Cu--O(2)--Cu--O(2)
%   |         |
%   O(3)      O(3)      y
%   |         |         |
%   Cu--O(2)--Cu        |-->x
%   (0) apical oxygen mode: O(4)-Cu where O(4) is right above Cu.
%   (1) breathing mode: O(2)-Cu, O(3)-Cu (half breathing: only O(2) or O(3)).
%   (2) buckling mode: Cu-O(2)-Cu,  Cu-O(3)-Cu. A1g: O(2)/O(3) in phase;
%   B1g: O(2)/O(3) out of phase.
%   From Bulut/Scalapino PRB1996:
%   (0) |g(q)|^2 = g_0^2;
%   (1) |g(q)|^2 = g_0^2*(sin(qx/2)^2 + sin(qy/2)^2);
%   (2) |g(q)|^2 = g_0^2*(cos(qx/2)^2 + cos(qy/2)^2);
%   where g_0^2 = (2*)|g|^2/(2M\omega_ph) [factor 2 for case (1) and
%   (2) because of two oxygens].
%   funtype:
%   'br': vq = sin(qx/2)^2 + sin(qy/2)^2
%   'bk': vq = cos(qx/2)^2 + cos(qy/2)^2
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function vq = vertexq(KX,KY,varargin)

    funtype = varargin{1};
    if strcmp(funtype,'exp')
        % Exponential forward scattering
        q0 = varargin{2};
        nk = size(KX,1)/2;
        MQ = 2*sqrt(sin(KX/2).^2 + sin(KY/2).^2); %=|q|
        vq = 4*4*(pi^2)*exp(-2*MQ/q0)/(2*pi*q0*q0);
    elseif strcmp(funtype,'frohlich')
        % Frohlich coupling (1/q)
        MQ = 2*sqrt(sin(KX/2).^2 + sin(KY/2).^2); %=|q|
        idx_mask = (MQ~=0);
        nk = size(KX,1)/2;
        vq(1:2*nk,1:2*nk) = 0;
        vq(idx_mask) = 1.0./MQ(idx_mask);
    elseif strcmp(funtype,'gauss')
        % Gaussian forward scattering
        q0 = varargin{2};
        nk = size(KX,1)/2;
        MQ = (KX-pi).^2 + (KY-pi).^2;
        MQ = MQ([nk+1:end 1:nk], [nk+1:end 1:nk]);
        vq = 4*(pi^2)*exp(-MQ/(2*q0*q0))/(2*pi*q0*q0);
    elseif strcmp(funtype,'uniform')
        % Uniform (standard BCS)
        vq = ones(size(KX));    
    elseif strcmp(funtype,'br')
        % Sine
        vq = sin(KX/2).^2 + sin(KY/2).^2;
    elseif strcmp(funtype,'bk')
        % Cosine
        vq = cos(KX/2).^2 + cos(KY/2).^2;
    else
        error([funtype, ' is not implemented!'])
    end
