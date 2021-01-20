%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% >>Function: getNCwt
% Get the Newton-Cotes weights for numerical integration.
%
% >>Usage:
% Compute the weights before you call fourier_bwd, where the weights are
% used to Fourier transform from (tau) to (i omega_n).
% Note that the getNCwt function relies on the function 'fI' defined below.
%
% >> Input:
% beta:     beta = 1/T, inverse temperature.
% Nc:       2*Nc = number of Matsubara frequencies.
% NCorder:  the order of Newton-Cotes formulae to calculate the Fourier
%           integral from \tau-space to \omega_n-space. See Press et al., 
%           Numerical Recipes. Specifically, the method used here is due 
%           to Filon (1928); nowadays alternative is Clenshaw-Curtis 
%           quadrature (QUADPACK FORTRAN90 library).
% isfermion: frequency grid shift is different for fermion and boson.
%
% >> Output:
% Wt:       weights for the corresponding Newton-Cotes order.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function Wt = getNCwt(beta,Nc,NCorder,isfermion)

s = (isfermion==1);
nwtot = 2*Nc;
del = beta/nwtot;
tau = (0:nwtot-1)*del;
vwn = (2*(-Nc:Nc-1) + s)*pi/beta;
idx = (vwn~=0);

% Select which NCorder you want (0/1/2)
if NCorder == 0
    Wt{1}(1:nwtot) = 0;
    Wt{1}(idx) = (exp(1i*vwn(idx)*del)-1)./(1i*vwn(idx));
    %the wn-->0 limit of the above
    Wt{1}(~idx) = del;
    Wt = Wt{1}(:);
elseif NCorder == 1
    for ii = 1:2
        Wt{ii}(1:nwtot) = 0;
    end
    Wt{1}(idx) = (2/del)*(1-cos(del*vwn(idx)))./(vwn(idx).^2);
    Wt{2}(idx) = (-1/del)*(1-1i*del*vwn(idx)-exp(-1i*del*vwn(idx)))./(vwn(idx).^2);
    %the wn-->0 limit of the above
    Wt{1}(~idx) = del;
    Wt{2}(~idx) = -del/2;
    Wt = [Wt{1}(:) Wt{2}(:)];
elseif NCorder == 2
    for ii = 1:7
        Wt{ii}(1:nwtot) = 0;
    end
    c1 = (del/4)*(fI(vwn,2,2,del) - 3*fI(vwn,2,1,del) + 2*fI(vwn,2,0,del));
    c2 = (-del/2)*exp(-1i*vwn*del).*(fI(vwn,2,2,del) - 2*fI(vwn,2,1,del));
    c3 = (del/4)*exp(-2i*vwn*del).*(fI(vwn,2,2,del) - fI(vwn,2,1,del));

    s1 = (del/4)*(fI(vwn,1,2,del) - 3*fI(vwn,1,1,del) + 2*fI(vwn,1,0,del));
    s2 = (-del/2)*exp(-1i*vwn*del).*(fI(vwn,1,2,del) - 2*fI(vwn,1,1,del));
    s3 = (del/4)*exp(-2i*vwn*del).*(fI(vwn,1,2,del) - fI(vwn,1,1,del));

    p1 = (del/4)*exp(1i*vwn*del).*(fI(vwn,1,2,del) - fI(vwn,1,1,del));
    p2 = (-del/2)*(fI(vwn,1,2,del) - fI(vwn,1,0,del));
    p3 = (del/4)*exp(-1i*vwn*del).*(fI(vwn,1,2,del) + fI(vwn,1,1,del));

    Wt{1} = c1 + c2 + c3;
    Wt{2} = -(c3+p3);
    Wt{3} = (s1-c2+p3).*exp(1i*vwn*tau(1));
    Wt{4} = (s2-c3).*exp(1i*vwn*tau(2));
    Wt{5} = (s3).*exp(1i*vwn*tau(3));
    Wt{6} = (p1).*exp(1i*vwn*tau(end-1));
    Wt{7} = (p2-c1).*exp(1i*vwn*tau(end));
    tmp = [];
    for ii = 1:7
        tmp = [tmp Wt{ii}(:)];
    end
    Wt = tmp;
else
    error('Wrong Newton-Cotes order!')
end

function out = fI(x,u,ord,dta)
indx = (x~=0);
out = zeros(size(x));
x = 1i*dta*x;
if ord == 0
    out(indx) = (exp(u*x(indx))-1)./x(indx);
    out(~indx) = u;
elseif ord == 1
    out(indx) = u*exp(u*x(indx))./x(indx) - (exp(u*x(indx))-1)./(x(indx).^2);
    out(~indx) = (u^2)/2;
elseif ord == 2
    out(indx) = u^2*exp(u*x(indx))./x(indx) - 2*u*exp(u*x(indx))./(x(indx).^2) ...
        + 2*(exp(u*x(indx))-1)./(x(indx).^3);
    out(~indx) = (u^3)/3;
else
    error('Wrong order of integral!')
end
