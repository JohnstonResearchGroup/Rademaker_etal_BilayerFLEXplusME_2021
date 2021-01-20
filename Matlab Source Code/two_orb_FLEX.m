%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bilayer FLEX + Phonons
%       Y. Wang, L. Rademaker, K. Nakatsukasa, G. Alvarez-Suchini, and S. Johnston           
%       Last Update: 10 August 2020
%
%Script: two_orb_FLEX
%
% This is the main script that calculates the FLEX calculations for a
% two-orbital model, looped over a list of temperature (Tlist) and U-values
% (Ulist).
%
% If you want to change from ladders to 2d lattices, make sure you include
% both hopping in the y direction (ty) as well as change the K-grid
% accordingly.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment the following line if you want to convert this script into a
% function:
% function two_orb_FLEX(filling0,U,ty,tz,Tlist,eps,maxchi_cutoff)

% Clear Matlab
clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contains both model parameters and computational parameters

maxchi_cutoff = 0.995;  % Artificial cut-off for the polarization 
                        % to avoid problems with convergence
filling0 = 2.1;         % Initial filling
                        % NOTE: filling0 = 2 n 
U = 7;                  % initial Hubbard U  (units of t)
Ulist = U;
%Ulist = [0.5:0.5:10];              % Make a list of U values
ty = 1;                 % Hopping in the y-direction 
                        % (set to zero for ladders)
tz = 2.36;              % Hopping in the z-direction (interlayer)
T = 0.01;               % initial temperature T (units of t)
Tlist = T;              % Make a list of temperatures (or just one)
%Tlist = [0.0075:0.0025:0.0125];
eps = 1e-6;             % Precision

% "profile" can be used to track the calculation, uncomment if desired
% profile on

%=1: add Hartree-Fock self-energies; =0: do not add.
AddHFSelfE = 1;

% What do to with the self-energies?
runInSerial = 1;    % switch between for-loop and parfor-loop
saveSelfE = 1;      % save self-energy
loadSelfE = 1;      % load self-energy at different T or U as the input
flySelfE = 0;       % load self-energy calculated on the fly
if (loadSelfE == 1) && (saveSelfE ~= 1)
    flySelfE = 1;
end


checkFreqSymm = 0;  % check frequency even/odd symmetry
%Nota Bene: 
% Eliashberg equation for normal self-energy S = 1i*wn*(1-Z) + X is used.
% Although there is no explicit frequency even/odd symmetry, the even/odd
% symmetry for real(S)/imag(S) is imposed every iteration. Most error
% (<1e-14) related to the symmetry might come from the fft transformation,
% which depends on the version of fftw (C library) called by Matlab, hence
% the Matlab version as well. Another error related to the fft is the
% approximation of Fourier integral with fft in the fourier_bwd when
% transforming from \tau-space to \omega_n-space. This can be reduced by
% using a higher Newton-Cotes formula.


useSymmetry = 1;    % C_2(inversion) symmetry is used
NCorder = 2;        % Newton-Cotes integration order
Power2NumFreq = 1;  % Whether the number of frequencies is a power of 2,
                    % setting it to '1' means better use of FFT


%--------------------------------------------------------------------------
% FORWARD ELECTRON-PHONON COUPLING

% Type of electron-phonon coupling
varargin = {'exp'};
% Options for varargin are:
% 'bk': buckling mode  ~cos(q/2)^2 [has CDW at q*=(0,0)]
% 'br': breathing mode ~sin(q/2)^2 [has CDW at q*=(pi,pi)]
% 'exp': forward scattering coupling with a peak at q=0, and width q0.
% 'gauss': gaussian forward scattering with peak at q=0, and width q0
% 'frohlich': 1/q coupling

% Properties of the coupling
wph = 4/3;          % Phonon frequency, in units of t
                    % Note that t=75 meV in FeSe, so (4/3)t=100 meV
q0 = 0.1;           % Forward scattering decay length, units: 1/a


% Strength of electron-phonon coupling
% Set this to zero to switch off e-ph interactions
lambda_mass = 0.191323;  
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Define the momentum grid and all of the constants that we need.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the number of k-points ([1 96] is ladder, [64 64] is 2d grid)
Nk = [32 32];
% Total number of k-points
nktot = 4*Nk(1)*Nk(2);
% Define the momentum grid
for ii = 1:2
    K{ii} = (0:1:2*Nk(ii)-1)*pi/Nk(ii);
end
[KX,KY] = meshgrid(K{2},K{1});      %row: y, column: x

% Define the Boltzmann constant (either 1 or actual value in eV/K)
kb = 1.0;
%kb = 8.6173303e-5;     %Source: 2014 CODATA, eV/K

% Not sure what this is, is also not used later
getlim = @(x) [min(x(:)) max(x(:))];

% Define the orbitals and the dispersion
Norb = 2;   % Nr of orbitals, has to be two
Nelm = 2;   % Nelm: number of elements stored for the matrix in orbital space
% Indexmap for the orbitals
Idxmap_sub = [1 1; 1 2];
Idxmap_ind = [1 2];        % Idxmap_ind = sub2ind([Norb Norb],Idxmap_sub(:,1),Idxmap_sub(:,2))'
% Energy unit is hopping in x-direction
t = 1.0;        %energy unit 1*t = 0.05eV
if Norb == 2
    % Relation between energy unit and eV
    t2eV = 1;    
    
    % Choose one of the possible Bloch Hamiltonians;
    % hk is only used to compute ek
    hk{1,1} = @(kx,ky) -2*t*cos(kx) - 2*ty*cos(ky);
    hk{2,1} = @(kx,ky) -tz*ones(size(kx));
    
    % Option with next-nearest neighbor hopping
    %hk{1,1} = @(kx,ky) -2*t*cos(kx) - 2*ty*cos(ky) + 0.6*cos(kx).*cos(ky);    
    %hk{2,1} = @(kx,ky) -0.25*tz*(cos(kx) - cos(ky)).^2;
    
    % Option with flipped sign of 
    %hk{1,1} = @(kx,ky) 2*t*cos(kx) + 2*ty*cos(ky);
    %hk{2,1} = @(kx,ky) -tz*ones(size(kx)); 
    
    % If necessary, symmetrize
    %hk{1,2} = hk{2,1};
    %hk{2,2} = hk{1,1};
    
    % Dispersion in orbital space
    ek(1:2*Nk(1),1:2*Nk(2),1:2) = 0;
    ek(:,:,1) = hk{1,1}(KX,KY);
    ek(:,:,2) = hk{2,1}(KX,KY);
    
    % Dispersion in band space
    % Is only used to estimate the density of states and fillings later,
    % not really needed for FLEX calculation
    xik(:,:,1) = ek(:,:,1) + ek(:,:,2);     % band: bonding, A
    xik(:,:,2) = ek(:,:,1) - ek(:,:,2);     % band: anti-bonding, B
                                            %for tz>0, xik(:,:,1) < xik(:,:,2)
    % Frequency cut-off in units of the bandwidth
    % Ncut = w_c/band_width
	Ncut = 5;       

else
    error('Number of orbitals must be 2!')
end

% Range of chemical potentials to check
muLim = [min(xik(:)) max(xik(:))];

% Frequency range
Wband = muLim(2) - muLim(1);
% Matsubara freq. cut-off w_c = Ncut*(band width)
nuscale = Ncut*Wband/2;

% Computational parameters
sig = 2*t*pi/Nk(1);     % precision, ~v_F*\Delta k ~ 4*t*pi/nk
mixing = 1;             % Should we include mixing?
wgt = 0.1;              % Mixing weight (new = wgt * new + (1-wgt) * old)
%mix_method: 0, mix self-energy; 1, mix Green's function; 
%            2, mix Green's function and use Anderson acceleration
mix_method = 0;

% When to stop the self-consistent iterations
maxiter = 20000;

%--------------------------------------------------------------------------
% FORWARD ELECTRON-PHONON COUPLING

% Set up the actual electron-phonon coupling strength
lam0 = lambda_mass*(2*wph);

% Read the variable argument
if isempty(varargin)
    error('Wrong forward coupling input!')
    vertexinput = {'uniform'} % default
else
    vertex_type = varargin{1};
    if strcmp(vertex_type,'br') || strcmp(vertex_type,'bk') || ...
            strcmp(vertex_type,'uniform') || strcmp(vertex_type,'frohlich')
        vertexinput = {vertex_type};       %'bk' or 'br' coupling
    elseif strcmp(vertex_type,'exp') || strcmp(vertex_type,'gauss')
        vertexinput = {vertex_type,q0};     %forward scattering coupling
    else
        error('Wrong type of coupling!')
        vertexinput = {'uniform'} % default
    end
end

% Get the momentum -dependence of the electron-phonon coupling
fqph = vertexq(KX,KY,vertexinput{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Setup the file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileDir = './';
filenamestr = ['_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) '_U=' num2str(U) ...
    '_lam=' num2str(lambda_mass) '_q0=' num2str(q0) ...
    '.dat'];
fileout = ['out' filenamestr];
filegaps = ['gaps' filenamestr];
filelamb = ['calc_lambda_mass' filenamestr];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Start outputting to the screen and file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'diary' saves all command line input and output to an external file
diary([fileDir, fileout]);
% start outputting
fprintf('\n');
fprintf('Gap calculation\n');
fprintf('  grid [nky nkx] = [%4d,%4d], convergence criterion = %g [t], smearing = %g [t]\n',...
    Nk(1),Nk(2),eps,sig)
fprintf('  U = %g [t]\n',U)
fprintf('  lambda_mass = %g, wph = %g [t], q0 = %g [1/a]\n',lambda_mass, wph, q0)
fprintf('  ty = %g [t], tz = %g [t]\n',ty,tz)
fprintf('  [xik_min, xik_max, Wband] = [%g, %g, %g] [t]\n',muLim(1),muLim(2),Wband)
fprintf('  Ncut = %g [t]\n',Ncut)
fprintf('  U*chi_max_cutoff = %g \n',maxchi_cutoff)
fprintf('  mixing = %d, weight = %g\n',mixing,wgt)
fprintf('  maxiter = %d\n',maxiter)
fprintf('  useSymmetry = %1d\n',useSymmetry)
fprintf('  NCorder = %1d\n',NCorder)
fprintf('  [loadSelfE, saveSelfE, flySelfE] = [%1d, %1d, %1d]\n',loadSelfE,saveSelfE,flySelfE)
fprintf('  runInSerial = %1d\n',runInSerial)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Start looping over T and U 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open output file for gaps vs T, U
fidgaps = fopen([fileDir filegaps],'a');
fidlamb = fopen([fileDir filelamb],'a');

% Set initial SC gap on the two bands to nonzero
gap_input = [0.01 0.01]*20;

% Clear temperature
T = [];

% Start temperature loop
for nt = 1:numel(Tlist)
    % Start U loop
    for iiu = 1:numel(Ulist)
        % Timing
        tic;
        
        % Set temperature
        Told = T;
        Uold = U;
        T = Tlist(nt);
        U = Ulist(iiu);
        
        % Define frequency scale
        beta = 1/T;
        wn = pi/beta;
        numwi = round(nuscale/wn);
        % Rescale, if necessary, the number of frequencies to be 
        %  a power of 2
        if Power2NumFreq == 1
            %numwi = 2^(ceil(log2(numwi)));
            numwi = max(numwi,256);
        end
        
        % Define parameter that helps calculating density in get_filling
        xifill0 = T*log(2*Norb/filling0 - 1);
        
        % Output of model parameters
        fprintf('\n')
        fprintf('\n')
        fprintf('**************************************************************\n')
        fprintf('**************************************************************\n')
        fprintf('  T = %g [t], U = %g [t]\n',T,U)
        fprintf('  xifill0 = %g [t]\n',xifill0)
        fprintf('  Ncut = %g [t] --> numwi =%6d, Power2NumFreq = %1d\n',Ncut,round(nuscale/wn),Power2NumFreq)
        fprintf('  w(n=0) = %12.8f [t], number of wn>0 numwi =%6d\n',wn,numwi)
        
        % If not the first temperature, store the old frequency grid
        if nt>1
            WNold = WN;
        end
        % Define the freqency grid for the current temperature
        WN = pi*(2*(-numwi:numwi-1)+1)/beta;    % Fermions
        WNU = pi*2*(-numwi:numwi-1)/beta;       % Bosons
        
        % Loading/saving self-energies
        if saveSelfE == 1              
            fileSelfE{nt,iiu} = ['selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2)) ...
                '_U=' num2str(U) '_lam=' num2str(lambda_mass) '_q0=' num2str(q0) ...
                '_T=' num2str(T) '.mat'];
            fprintf('  self-energy %s%s will be saved.\n',fileDir,fileSelfE{nt,iiu})
        end
        if (saveSelfE == 1) && exist([fileDir fileSelfE{nt,iiu}],'file')
            fprintf('  self-energy %s%s exists. Go to next T and U.\n',fileDir,fileSelfE{nt,iiu})
            continue
        end
        if (loadSelfE == 1) && (nt>1 || iiu>1)
            if numel(Tlist)>1 && numel(Ulist)>1 && saveSelfE~=1
                error('  set saveSelfE = 1!')
            end
            if numel(Tlist) == 1 && numel(Ulist) > 1
                if flySelfE ~= 1
                    load([fileDir fileSelfE{nt,iiu-1}],'S','P','WN','mu','-mat')
                end
            elseif  numel(Tlist) > 1 && numel(Ulist) == 1
                if flySelfE == 1
                    [S P mu] = interpSelfE(WN,1,S,P,WNold,mu);
                else
                    [S P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,1}]);
                end
            else
                if (Told~=T)
                    [S P mu] = interpSelfE(WN,0,[fileDir fileSelfE{nt-1,iiu}]);
                elseif (Uold~=U)
                    load([fileDir fileSelfE{nt,iiu-1}],'S','P','WN','mu','-mat')
                end
            end
            wn = pi/beta;
            for jj = 1:Nelm
                gap(jj) = max(max( abs(real(P(:,:,numwi+1,jj))) ));
                if abs(gap(jj)) < 0.001
                    for nn = 1:length(WN)
                        P(:,:,nn,jj) = gap_input(jj);
                    end
                end
            end
        else
            % setup the arrays to store the self-energies on the imaginary
            % axis and initialize the self-energies
            clear('S','P')
            S(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
            P(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
            for jj = 1:Nelm
                for nn = 1:2*numwi
                    % Set initial self-energy zero or small Z
                    %S(:,:,nn,jj) = -1i*sign(WN(nn))*0.05;
                    S(:,:,nn,jj) = 0;
                    
                    % Note: set s-wave or d-wave initial suggestion for gap
                    P(:,:,nn,jj) = gap_input(jj);
                    %P(:,:,nn,jj) = gap_input(jj)*(cos(KX)-cos(KY));
                end
            end
            mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
        end
        
        % Get iniital chemical potential
        mu0 = get_mu(zeros(size(S)),zeros(size(P)),WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
        % Initial filling for current 
        filling = get_filling(S,P,WN,ek,mu,xifill0,Norb,beta,useSymmetry);
        
        % Check the different ways to compute the fillings to the initial
        % provided filling
        fk = fermi(xik(:)-mu0,1/T);
        fill_fm = 2*sum(fk(:))/nktot;  %2 for two spins        
        fprintf('  fill_fm = %12.8f from fermi function sum\n',fill_fm)
        fprintf('  filling0= %12.8f, mu0= %12.8f [t]\n',filling0,mu0)        
        fprintf('  filling = %12.8f, mu = %12.8f [t]\n',filling,mu)

        % Compute estimates of density of states
        dk = gauss(xik-mu,sig);
        Nf0 = sum(dk(:))/nktot;        
        dfk = -beta*exp(beta*(xik-mu))./((exp(beta*(xik-mu))+1).^2);
        Nf = -sum(dfk(:))/nktot;
        fprintf('  N(E_f,T =%12.8f [t]) = %12.8f (1/[t]/per spin/Volume)\n',sig,Nf0)
        fprintf('  N(E_f,T =%12.8f [t]) = %12.8f (1/[t]/per spin/Volume)\n',T,Nf)
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The actual calculation is done here
        calculate_im_axis;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % The rest of this loop is output to file
        
        % choose the lowest Matsubara frequency
        wn = pi/beta;
        
        % If no interactions set maxchi to zero
        if (U == 0), maxchi = 0; end
        % Check the gap
        if (maxchi < 1)
            % Calculate the gap
            for jj = 1:Nelm
                gapwn = real(P(:,:,numwi+1,jj));
                gapmax(jj) = max(gapwn(:));
                gapmin(jj) = min(gapwn(:));
            end
            % Output the gaps to file
            if Nelm == 2
                %ekmax_h = max(max(ek(:,:,idxh) + X(:,:,numwi+1,idxh) - mu));
                %ekmin_e = min(min(ek(:,:,idxe) + X(:,:,numwi+1,idxe) - mu));
                fprintf(fidgaps,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %g\n', ...
                    T,mu,gapmax(1),gapmin(1),gapmax(2),gapmin(2),eps);
            elseif Nelm == 3
                fprintf(fidgaps,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %g\n', ...
                    T,mu,gapmax(1),gapmin(1),gapmax(2),gapmin(2),gapmax(3),gapmin(3),eps);
            else
                error('Number of orbitals must be 2 or 3!')
            end
        end
        if (U==0), clear('maxchi'); end
        
        % Save self-energy to file
        if saveSelfE == 1
            save([fileDir fileSelfE{nt,iiu}],'S','P','Vnm','Van','VnmMomentum','VanMomentum','WN','mu','-mat')
        end
        
        % Timing
        vtm = toc;
        
        % Save the polarizations and self-energies to file
        % First prepare the arrays
        % save(['T=' num2str(T) '_U=' num2str(U) '.mat'],'nk','numwi','T','mu','WN','Z','X','P',...
        %     'WNU','chis0','chic0','Vsf','Vcf','-mat')
        % q dependence
        chis00 = permute(chis0(:,:,numwi+1,:),[1 2 4 3]);
        chic00 = permute(chic0(:,:,numwi+1,:),[1 2 4 3]);
        chisR0 = permute(chisR(:,:,numwi+1,:),[1 2 4 3]);
        chicR0 = permute(chicR(:,:,numwi+1,:),[1 2 4 3]);
        Vnm0 = permute(Vnm(:,:,numwi+1,:),[1 2 4 3]);
        Vnm0Momentum = permute(VnmMomentum(:,:,numwi+1,:),[1 2 4 3]);
        Van0Momentum = permute(VanMomentum(:,:,numwi+1,:),[1 2 4 3]);
        S0 = permute(S(:,:,numwi+1,:),[1 2 4 3]);
        P0 = permute(P(:,:,numwi+1,:),[1 2 4 3]);
        % w dependence q = (0 0)
        chis0w1 = reshape(chis0(1,1,:,:),1,[]);
        chic0w1 = reshape(chic0(1,1,:,:),1,[]);
        chisRw1 = reshape(chisR(1,1,:,:),1,[]);
        Vnmw1 = reshape(Vnm(1,1,:,:),1,[]);
        Sw1 = reshape(S(1,1,:,:),1,[]);
        Pw1 = reshape(P(1,1,:,:),1,[]);
        % w dependence q = (pi pi)
        chis0w2 = reshape(chis0(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        chic0w2 = reshape(chic0(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        chisRw2 = reshape(chisR(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Vnmw2 = reshape(Vnm(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Sw2 = reshape(S(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        Pw2 = reshape(P(Nk(1)+1,Nk(2)+1,:,:),1,[]);
        % w dependence q = (0 pi)
        chis0w3 = reshape(chis0(Nk(1)+1,1,:,:),1,[]);
        chic0w3 = reshape(chic0(Nk(1)+1,1,:,:),1,[]);
        chisRw3 = reshape(chisR(Nk(1)+1,1,:,:),1,[]);
        Vnmw3 = reshape(Vnm(Nk(1)+1,1,:,:),1,[]);
        Sw3 = reshape(S(Nk(1)+1,1,:,:),1,[]);
        Pw3 = reshape(P(Nk(1)+1,1,:,:),1,[]);
        % Output to file
        save([fileDir 'selfE_Nk=' num2str(Nk(1)) '_' num2str(Nk(2))  '_T=' num2str(T) '_U=' num2str(U) '_freqFFT.mat'],...
            'Nk','numwi','T','mu','WN','WNU',...
            'S0','P0','chi*00','chi*R0','Vnm0','Vnm0Momentum','Van0Momentum',...
            'Sw*','Pw*','chi*0w*','chisRw*','Vnmw*','ek','-mat')
        
        % Uncomment the following if you want to use 'profile'
        % profile viewer
        % p = profile('info');
        % profsave(p,'profile_results')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate lambda_m based on fermi surface avergae: (done with each
        % band)
        S0pp = S0(:,:,1) + S0(:,:,2);
        S0ppC = conj(S0(:,:,1)) + conj(S0(:,:,2));
        Zpp = ones(2*Nk(1),2*Nk(2)) - ((S0pp - S0ppC)/(2*1i*pi*T));
        S0mm = S0(:,:,1) - S0(:,:,2);
        S0mmC = conj(S0(:,:,1)) - conj(S0(:,:,2));
        Zmm = ones(2*Nk(1),2*Nk(2)) - ((S0mm - S0mmC)/(2*1i*pi*T));
        sigma = 0.005;
        eplus = real(((ek(:,:,1) + ek(:,:,2)) - (mu*ones(size(S0(:,:,1))))+(S0(:,:,1)+S0(:,:,2))));
        eminus = real(((ek(:,:,1) - ek(:,:,2)) - (mu*ones(size(S0(:,:,1))))+(S0(:,:,1)-S0(:,:,2))));
        delta_sum1 = 0;
        for nkx = 1:Nk(1)+1
            kx = KX(nkx);
            for nky = 1:Nk(2) +1
                delta_sum1 = delta_sum1 + gauss2(eplus(nkx,nky),sigma);
            end
        end
        delta_sum1 = delta_sum1/(4*Nk(1)*Nk(2));
        
        
        lam1 = 0;
        for nkx = 1:Nk(1)+1
            kx = KX(nkx);
            for nky = 1:Nk(2) +1
                lam1 = lam1 + (Zpp(nkx,nky)-1)*gauss2(eplus(nkx,nky),sigma);
            end 
        end
        lam1 = lam1/(4*Nk(1)*Nk(2));
        
        lambda_m1 = lam1/delta_sum1;
        
        delta_sum2 = 0;
        for nkx = 1:Nk(1)+1
            kx = KX(nkx);
            for nky = 1:Nk(2) +1
                delta_sum2 = delta_sum2 + gauss2(eminus(nkx,nky),sigma);
            end
        end
        delta_sum2 = delta_sum2/(4*Nk(1)*Nk(2));
        
        
        lam2 = 0;
        for nkx = 1:Nk(1)+1
            kx = KX(nkx);
            for nky = 1:Nk(2) +1
                lam2 = lam2 + (Zmm(nkx,nky)-1)*gauss2(eminus(nkx,nky),sigma);
            end 
        end
        lam2 = lam2/(4*Nk(1)*Nk(2));
        
        lambda_m2 = lam2/delta_sum2;
        fprintf(fidlamb,'%12.8f %12.8f %12.8f  %g\n', ...
                    T,mu,lambda_m1,lambda_m2);
        %%%%
        fprintf('Done gap calculation. Total Time = %.2f s\n',sum(vtm));
        
    end % End U loop
end % End T loop

% Close output file for gaps and diary
fclose(fidgaps);
fclose(fidlamb);
diary off
