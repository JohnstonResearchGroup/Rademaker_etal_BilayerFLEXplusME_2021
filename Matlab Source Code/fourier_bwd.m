%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 26 March 2020
%
% >>Function: fourier_bwd
% Fourier transform from (r,\tau) to (k,i\omega_n).
%
% >>Usage:
% gin = fourier_bwd(gin,dg,NCorder,Wt,beta,Nk,Nc,isfermion,useSymmetry)
%
% >>Input:
% gin:      function in (r,\tau) space.
% dg:       discontinuity of Green's function at \tau=0. dg = f(0+)-f(0-).
% NCorder:  the order of Newton-Cotes formulae to calculate the Fourier
%           integral from \tau-space to \omega_n-space. See Press et al., 
%           Numerical Recipes. Specifically, the method used here is due 
%           to Filon (1928); nowadays alternative is Clenshaw-Curtis 
%           quadrature (QUADPACK FORTRAN90 library).
% Wt:       weights for the corresponding Newton-Cotes order.
% beta:     beta = 1/T, inverse temperature.
% Nk:       Nk = [Nk(1) Nk(2)], k-grid and r-grid.
% Nc:       2*Nc = number of Matsubara frequencies.
% isfermion: frequency grid shift is different for fermion and boson.
% useSymmetry: use C_2 symmetry for k-grid. Nk(1) and Nk(2) must both be
%           even number if useSymmetry=1. Set useSymmetry=0 and 
%           Nk = [1/2 1/2] to calculate one k-point (\Gamma point).
% varargin: whether to run in serial (1) or parallel (else).
%
% >>Output:
% gin: same variable name as input gin to allow in-place operations without
% allocating more memory for the output. Without explicit pointers, it
% requires careful coding to make Matlab recognize the in-place operations:
% (1) output/input share the same variable name; (2) the in-place function
% must itself be called from another function; (3) changing the size of gin
% might not be allowed within the in-place function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function gin = fourier_bwd(gin,dg,NCorder,Wt,beta,Nk,Nc,isfermion,useSymmetry,varargin)

% Run serial or parallel?
if isempty(varargin)
    runInSerial = 1;
else
    runInSerial = varargin{1};
end


% Determine and check size of the matrix
siz_g = size(gin);
if ~isequal(siz_g,[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array gin has a wrong size.')
end
nktot = 4*Nk(1)*Nk(2);  %# of sites
nwtot = 2*Nc;           %# of tau slices

% ifft on r grid ----------------------------------------------------------
% 1/nktot normalization is included in direct fft to match the definition 
% of Green's function in real space; thus the 1/nktot in Matlab's ifft
% definition must be removed by "*nktot".

imeps = 1e-100i;
% Nota Bene: imeps is used to keep Matlab from deallocating pointers for
% the imaginary part when they equal zero. (The strange behavior may depend
% on the version of Matlab. Most recent versions do not need this awkward
% trick.)

% inverse Fourier transform from (r,tau) to (k,tau)
if runInSerial == 1
    for nn = 1:nwtot
        gin(:,:,nn) = ifft2(gin(:,:,nn))*nktot + imeps;
    end
else
    % Parallel
    parfor nn = 1:nwtot
        gin(:,:,nn) = ifft2(gin(:,:,nn))*nktot + imeps;
    end
end
gin = gin - imeps;


% ifft on tau grid --------------------------------------------------------

% If useSymmetry, only to tau-wn Fourier transform on 1/4 of the momentum
% points
if useSymmetry
    gin = gin(1:Nk(1)+1,1:Nk(2)+1,1:nwtot);
    nks = (Nk(1)+1)*(Nk(2)+1);
else
    nks = nktot;
end

% Fermions or bosons?
s = (isfermion==1);         %fermion: s=1; boson: s=0;

% Make list of Fourier modes exp(i wn tau)
pre = exp(1i*pi*s/nwtot*(0:nwtot-1));

% Make gin (nwtot,nks)-array to avoid potential memory spikes
gin = transpose(reshape(gin,nks,nwtot));


% For parfor loops, when you index into a sliced variables,
% restrictions are placed on the first-level variable indices.
% This allows parfor to easily distribute the right part of the
% variable to the right workers. One of these first-level
% indices must be the loop counter variable or the counter
% variable plus or minus a constant. Every other first-level
% index must be a constant, a non-loop counter variable, a
% colon, or an end. So gin(:,ii) = gin([Nc+1:end, 1:Nc],ii) is not
% allowed in parfor loops.
%

% Note thatthe inverse Fourier transform from tau to omega_n is an integral
% and therefore we can either use a simple summation (NCorder==0) or
% Newton-Cotes methods (1,2). The weights are given externally.
if NCorder == 0
    %parfor
    if runInSerial == 1
        for ii = 1:nks
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gin(:,ii) = gin([Nc+1:end, 1:Nc],ii);
            gin(:,ii) = gin(:,ii).*Wt(:,1);
        end
    else
        parfor ii = 1:nks
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gtmp = gin(:,ii);
            gin(:,ii) = gtmp([Nc+1:end, 1:Nc],1).*Wt(:,1);
        end
    end
elseif NCorder == 1
    %parfor
    if runInSerial == 1
        for ii = 1:nks
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gin(:,ii) = gin([Nc+1:end, 1:Nc],ii);
            gin(:,ii) = gin(:,ii).*Wt(:,1) + dg*Wt(:,2);
        end
    else
        parfor ii = 1:nks
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gtmp = gin(:,ii);
            gin(:,ii) = gtmp([Nc+1:end, 1:Nc],1).*Wt(:,1) + dg*Wt(:,2);
        end
    end
elseif NCorder == 2
    %parfor
    if runInSerial == 1
        for ii = 1:nks
            gtmp = gin([1 2 3 end-1 end],ii);
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gin(:,ii) = gin([Nc+1:end, 1:Nc],ii);
            gin(:,ii) = gin(:,ii).*Wt(:,1) + dg*Wt(:,2) + ...
                Wt(:,3:end)*gtmp;
        end
    else
        parfor ii = 1:nks           
            gtmpold = gin(:,ii);
            gin(:,ii) = ifft(gin(:,ii).*pre(:))*nwtot;       %remove normalization from ifft
            gtmp = gin(:,ii);
            gin(:,ii) = gtmp([Nc+1:end, 1:Nc],1).*Wt(:,1) + dg*Wt(:,2) + ...
                Wt(:,3:end)*gtmpold([1 2 3 end-1 end],1);
        end
    end
else
    error('Wrong Newton-Cotes order!')
end

% If we used symmetry, change gin back to (2*Nk(1), 2*Nk(2), nwtot)-array
if useSymmetry
    tmp = reshape(transpose(gin),Nk(1)+1,Nk(2)+1,nwtot);
    gin = tmp([1:end, end-1:-1:2],[1:end, end-1:-1:2],:);
else
    gin = reshape(transpose(gin),2*Nk(1),2*Nk(2),nwtot);
end
