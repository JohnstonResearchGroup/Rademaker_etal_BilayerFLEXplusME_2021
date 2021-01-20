%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 26 March 2020
%
% >>Function: fourier_fws
% Fourier transform from (k,i\omega_n) to (r,\tau).
%
% >>Usage:
% gin = fourier_fwd(gin,ek,beta,Nk,Nc,addtail,isfermion,useSymmetry)
%
% >> Input:
% gin:      function in (k,i\omega_n) space.
% ek:       dispersion on Nk(1)*Nk(2) k-grid, used only for adding tail to
%           Green's function.
% beta:     beta = 1/T, inverse temperature.
% Nk:       Nk = [Nk(1) Nk(2)], k-grid and r-grid.
% Nc:       2*Nc = number of Matsubara frequencies.
% addtail:  1 (add tail) / 0(no tail).
% isfermion: frequency grid shift is different for fermion and boson.
% useSymmetry: use C_2 symmetry for k-grid. Nk(1) and Nk(2) must both be
%           even number if useSymmetry=1. Set useSymmetry=0 and 
%           Nk = [1/2 1/2] to calculate one k-point (\Gamma point).
% varargin: whether to run in serial (1) or parallel (else).
%
% >> Output:
% gin: same variable name as input gin to allow in-place operations without
% allocating more memory for the output. Without explicit pointers, it
% requires careful coding to make Matlab recognize the in-place operations:
% (1) output/input share the same variable name; (2) the in-place function
% must itself be called from another function; (3) changing the size of gin
% might not be allowed within the in-place function.
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gin = fourier_fwd(gin,ek,beta,Nk,Nc,addtail,isfermion,useSymmetry,varargin)

% Run in parallel or serial?
if isempty(varargin)
    runInSerial = 1;
else
    runInSerial = varargin{1};
end

% Determine and check size of the input matrix
siz_g = size(gin);
siz_e = size(ek);
if ~isequal(siz_g,[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array gin has a wrong size.')
end
if addtail && ~isequal(siz_e,[2*Nk(1) 2*Nk(2)])
    error(' Input array ek has a wrong size.')
end
nktot = 4*Nk(1)*Nk(2);
nwtot = 2*Nc;

% fft on wn_n/nu_m grid ---------------------------------------------------
% If we useSymmetry, only take 1/4 of the momentum points
if useSymmetry
    gin = gin(1:Nk(1)+1,1:Nk(2)+1,1:nwtot);
    nks = (Nk(1)+1)*(Nk(2)+1);
else
    nks = nktot;
end

s = (isfermion==1);                 %fermion: s=1; boson: s=0;
del = beta/nwtot;
tau = (0:nwtot-1)*del;
vwn = (2*(-Nc:Nc-1)+1)*pi/beta;     %fermion grid
pre = exp(1i*pi*(nwtot-s)/nwtot*(0:nwtot-1))/beta;

% Make gin (nwtot,nks)-array to avoid potential memory spikes
gin = transpose(reshape(gin,nks,nwtot));

% Compute the high-frequency tail here
% This tail gives us the G(tau) for the noninteracting system
if addtail
    if useSymmetry
        ek = ek(1:Nk(1)+1,1:Nk(2)+1);
    end
    ek = reshape(ek,1,nks);
    % Subtract noninteracting Greens function 1/(i wn - ek)
    gin = gin - 1./(repmat(1i*vwn(:),1,nks) - repmat(ek,nwtot,1));

    % Compute the tail
    tail(1:nwtot,1:nks) = 0;
    % Serial or parallel forloop over momenta
    if runInSerial == 1
        % Serial 
        for ii = 1:nks
            %gin(:,ii) = gin(:,ii) - 1./(1i*vwn(:)-ek(ii));
            %Nota Bene: The above is very slow in the special case when
            %gin(:,ii)-1./(1i*vwn(:)-ek(ii))=0. Matlab will deallocate the
            %memory for the imaginary part and convert the whole array into
            %real. The deallocations take a very long time if they are looped
            %over a large nks point. The current solution is to move it out
            %side the loop. TODO: look for a better solution. (The strange
            %behavior may depend on the version of Matlab.)

            %set up tail; must prevent overflow (important at low T/large beta)
            if ek(ii)>0
                tail(:,ii) = exp(-tau(:)*ek(ii))/(exp(-beta*ek(ii)) + 1);
            else
                tail(:,ii) = exp((beta-tau(:))*ek(ii))/(exp(beta*ek(ii)) + 1);
            end
        end
    else
        % Parallel
        parfor ii = 1:nks
            %set up tail; must prevent overflow (important at low T/large beta)
            if ek(ii)>0
                tail(:,ii) = exp(-tau(:)*ek(ii))/(exp(-beta*ek(ii)) + 1);
            else
                tail(:,ii) = exp((beta-tau(:))*ek(ii))/(exp(beta*ek(ii)) + 1);
            end
        end        
    end
end

% Fast Fourier transform from (k,i wn) to (k,tau)
if runInSerial == 1
    for ii = 1:nks
        gin(:,ii) = fft(gin(:,ii)).*pre(:);
    end
else
    parfor ii = 1:nks
        gin(:,ii) = fft(gin(:,ii)).*pre(:);
    end    
end

% Subtract the tail is necessary
if addtail
    gin = gin - tail;
end

% Change gin back to (2*Nk(1), 2*Nk(2), nwtot)-array
if useSymmetry
    tmp = reshape(transpose(gin),Nk(1)+1,Nk(2)+1,nwtot);
    gin = tmp([1:end, end-1:-1:2],[1:end, end-1:-1:2],:);
else
    gin = reshape(transpose(gin),2*Nk(1),2*Nk(2),nwtot);
end

% FFT on k/q grid ---------------------------------------------------------
% From G(k, tau) to G(r, tau)
% 1/nktot normalization is included in direct fft to match the definition
% of Green's function in real space
% Serial or Parallel
if runInSerial == 1
    for nn = 1:nwtot
        gin(:,:,nn) = fft2(gin(:,:,nn))/nktot;
    end
else
    parfor nn = 1:nwtot
        gin(:,:,nn) = fft2(gin(:,:,nn))/nktot;
    end
end
