%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% Script: calculate_im_axis
% 
% This contains the self-consistent FLEX calculation. In the current
% implementation, it requires the set-up script of two_orb_FLEX to prepare,
% and it calls all other functions.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User output
fprintf('\n')
fprintf('**************************************************************\n')
fprintf('Begin imaginary axis calculation.\n')

% Check Matsubara frequencies
if (2*numwi) ~= numel(WN)
    error('Wrong Matsubara frequency grid!')
end
% Compute Newton-Cotes weights
fmWt = getNCwt(beta,numwi,NCorder,1);   %fermion freq weight
bsWt = getNCwt(beta,numwi,NCorder,0);   %boson freq weight

% Integration measure
pre = 1/(nktot*beta);

% Compute initial filling
filling = get_filling(S,P,WN,ek,mu,xifill0,Norb,beta,useSymmetry);
fprintf('  filling  = %12.8f, mu = %12.8f [t]\n',filling,mu)

% Clear the arrays changing size (numwi) at different temperatures
clear('Gn','Fn','chis0','chic0','chisR','chicR','Vnm','Van','Id','Gbeta','Fbeta')

% Define arrays for Fourier transforms
Gn(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
Fn(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
chis0(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
chic0(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
chisR(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
chicR(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
Vnm(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
Van(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
VnmMomentum(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
VanMomentum(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:Nelm) = 0;
nkwpt = numel(chis0)/Nelm;

% Define identity array (depends on whether or not using symmetry)
if useSymmetry
    Id(1:(Nk(1)+1),1:(Nk(2)+1),1:2*numwi,1:(Norb^2)) = 0;
else
    Id(1:2*Nk(1),1:2*Nk(2),1:2*numwi,1:(Norb^2)) = 0;
end
idxones = diag(reshape(1:(Norb^2),[Norb Norb]));
Id(:,:,:,idxones) = 1;
Id = reshape(Id,[],Norb^2);

% Array for end point (in tau-space)
Gbeta(1:2*Nk(1),1:2*Nk(2)) = 0;
Fbeta(1:2*Nk(1),1:2*Nk(2)) = 0;

%--------------------------------------------------------------------------
% FORWARD ELECTRON-PHONON SCATTERING

% Boson frequency grid WNU is be defined in two_orb_FLEX code
% Note that invD0 and fqphw do NOT carry an orbital index

% Inverse of phonon Green's function D0
invD0 = reshape(repmat(1 + (WNU/wph).^2, [4*Nk(1)*Nk(2), 1]), [2*Nk(1), 2*Nk(2), 2*numwi]);
% Extend matrix fqph to bosonic freqeuncy dimension
fqphw = repmat(fqph, [1, 1, 2*numwi]);
%--------------------------------------------------------------------------


% Initialise the self-consistent loop
loop = 1;
iter = 0;
fprintf(['  Iter.  _______dS_______  _______dP_______',...
    '  _______mu_______  _______<n>______',...
    '  _______gap______\n'])



% Prepare Anderson acceleration
if (mixing == 1) && (mix_method == 2)
    if useSymmetry
        nptmatG = (Nk(1)+1)*(Nk(2)+1)*2*numwi*Nelm;
    else
        nptmatG = 2*Nk(1)*2*Nk(2)*2*numwi*Nelm;
    end
    xx(1:nptmatG*2,1) = 0;  %store (Gnold, Fnold)<--*2 [or (Sold, Pold)] in one column
    mMax = 5;       %mMax>0, 5
    % maximum number of stored residuals (non-negative integer); should not
    % exceed number of rows of xx vectors (so we have a tall/thin matrix in
    % QR decompostion)
    droptol = 1.e10;
    % tolerance for dropping stored residual vectors to improve
    % conditioning: If droptol > 0, drop residuals if the condition number
    % exceeds droptol; if droptol <= 0, do not drop residuals.
    %dampbt = @(x)0.2*(x<20) + 0.2*(x>=20 && x<25) + 0.5*(x>=25); %0.2
    %dampbt = @(x)1*(x<15) + 0.4*(x>=15);
    dampbt = 1;
    % damping factor: If dampbt > 0 (and dampbt ~= 1), then the step is damped
    % by dampbt; otherwise, the step is not damped. NOTE: dampbt can be a
    % function handle; form dampbt(iter), where iter is the iteration number
    % and 0 < dampbt(iter) <= 1.
    AAstart = 1;    %AAstart>=1, 5
    % acceleration delay factor: If AAstart > 0, start acceleration when
    % iter = AAstart.
    res_hist = [];
    % residual history matrix (iteration numbers and residual norms).
    DG = [];
    % Storage of g-value differences.
    R = []; Q = [];
    mAA = 0;
    % Initialize the number of stored residuals.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the selfconsistent loop
while (loop == 1)
    % Count the iterations
    iter = iter + 1;
       
    % Update old total self-energies (=S+S^{HF}) (with mixing)
    if (mixing == 1) && (iter>1)
        if (mix_method == 0)
            Sold = wgt*S + (1-wgt)*Sold;
            Pold = wgt*P + (1-wgt)*Pold;
            %update Green's function G^{+,+} = Gn, G^{+,-} = Fn;
            [Gn Fn] = solve_dyson(Sold,Pold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
        elseif mix_method == 1
            %mix Green's function Gn and Fn
            %(a) update Gn/Fn with pre-mixed self-energies
            [Gn Fn] = solve_dyson(S,P,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
            %(b) mix the Green's function
            Gnold = wgt*Gn + (1-wgt)*Gnold;
            Fnold = wgt*Fn + (1-wgt)*Fnold;
            %(c) find the self-energies from the mixed Green's function
            [Sold Pold] = solve_dysonbwd(Gnold,Fnold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
            Gn = Gnold;
            Fn = Fnold;
        else
            if mix_method ~= 2
                error('mix_method = %d is not supported!',mix_method)
            end
        end
        %
        %Recalculate chemical potential mu after mixing. Note: this
        %may cause convergence problem if the mixed Green's
        %function or self-energies overshoot too much so it is
        %commented now.
        %mu = get_mu(Zold,Xold,Pold,WN,ek,Nk,beta,filling0);
        %
    else
        Sold = S;
        Pold = P;
        %update Green's function G^{+,+} = Gn, G^{+,-} = Fn;
        [Gn Fn] = solve_dyson(Sold,Pold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
        %save old Green's function Gn(k,i\omega_n) and Fn(k,i\omega_n) for mixing 
        Gnold = Gn;
        Fnold = Fn;
    end

    % Update Hartree-Fock self-energies -----------------------------------
    % It should not be folded into FLEX interaction because U is 
    % instantaneous.
    S_HF(1:Nelm) = 0;
    P_HF(1:Nelm) = 0;
    if AddHFSelfE == 1
        if Nelm == 2
            P_HF(1) = U*sum(reshape(Fn(:,:,:,1),[],1))*pre;
            %S_HF(1)=0 because orb1 equivalent to orb2
            %S_HF(2)=P_HF(2)=0 because U'=J=J'=0
            %1=(1,1)_orb, 2=(2,1)_orb.
        elseif Nelm == 3
            for jj = [1 3]
                S_HF(jj) = U*sum(reshape(Gn(:,:,:,jj),[],1))*pre;
                P_HF(jj) = U*sum(reshape(Fn(:,:,:,jj),[],1))*pre;
            end
            S_HF([1 3]) = S_HF([1 3]) - (S_HF(1)+S_HF(3))/2;
            %S_HF(2)=P_HF(2)=0 because U'=J=J'=0
            %1=(1,1)_orb, 2=(2,1)_orb, 3=(2,2)_orb.
        else
            error('The code currently only deals special two-orbital case!')
        end
    end
    
    % ---------------------------------------------------------------------
    % Compute the self-energy ---------------------------------------------
    % Following steps:
    % 1. Forward FT Greens function
    % 2. Compute susceptibilities chi*0
    % 3. Add RPA 
    % 4. Construct Vnm (FLEX + Forward)
    % 5. Fourier transform Vnm
    % 6. Inverse Fourier transform Vnm*Gn

    % forward Fourier transform Green's function (k,i\omega_n) to (r,\tau)
    % G(r,\tau) and F(r,\tau) are real.
    if Nelm==2   
        idxdiag = 1;
    elseif Nelm==3
        idxdiag = [1 3];
    else
        error('The code currently only deals special two-orbital case!')
    end    
    for jj = 1:Nelm
        if (jj==1) || (jj==3)
            Gn(:,:,:,jj) = real(fourier_fwd(Gn(:,:,:,jj),xifill0(ones(2*Nk)),beta,Nk,numwi,1,1,useSymmetry,runInSerial));
        else
            Gn(:,:,:,jj) = real(fourier_fwd(Gn(:,:,:,jj),[],beta,Nk,numwi,0,1,useSymmetry,runInSerial));
        end
        Fn(:,:,:,jj) = real(fourier_fwd(Fn(:,:,:,jj),[],beta,Nk,numwi,0,1,useSymmetry,runInSerial));
    end
    
    % Compute the spin (chis0) & charge (chic0) susceptibilities
    % diagonal elements: chi_{aa,aa}, chi_{bb,bb}
    for jj = idxdiag
        Gbeta = -Gn(:,:,1,jj);
        Gbeta(1,1) = Gbeta(1,1) - 1;
        chis0(:,:,1,jj) = Gn(:,:,1,jj).*Gbeta - Fn(:,:,1,jj).*Fn(:,:,1,jj);
        chic0(:,:,1,jj) = Gn(:,:,1,jj).*Gbeta + Fn(:,:,1,jj).*Fn(:,:,1,jj);
        chis0(:,:,2:end,jj) = Gn(:,:,2:end,jj).*Gn(:,:,end:-1:2,jj) - Fn(:,:,2:end,jj).*Fn(:,:,2:end,jj);
        chic0(:,:,2:end,jj) = Gn(:,:,2:end,jj).*Gn(:,:,end:-1:2,jj) + Fn(:,:,2:end,jj).*Fn(:,:,2:end,jj);

    end
    % off-diagonal elements: chi_{aa,bb}=chi_{bb,aa}
    for jj = 2
        Gbeta = -Gn(:,:,1,jj);
        Fbeta = -Fn(:,:,1,jj);
        chis0(:,:,1,jj) = Gn(:,:,1,jj).*Gbeta + Fn(:,:,1,jj).*Fbeta;
        chic0(:,:,1,jj) = Gn(:,:,1,jj).*Gbeta - Fn(:,:,1,jj).*Fbeta;
        chis0(:,:,2:end,jj) = Gn(:,:,2:end,jj).*Gn(:,:,end:-1:2,jj) + Fn(:,:,2:end,jj).*Fn(:,:,end:-1:2,jj);
        chic0(:,:,2:end,jj) = Gn(:,:,2:end,jj).*Gn(:,:,end:-1:2,jj) - Fn(:,:,2:end,jj).*Fn(:,:,end:-1:2,jj);
    end

    % Backward Fourier transform susceptibility (r,\tau) to (k,i\omega_n)
    for jj = 1:Nelm
        chis0(:,:,:,jj) = fourier_bwd(chis0(:,:,:,jj),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
        chic0(:,:,:,jj) = fourier_bwd(chic0(:,:,:,jj),0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial);
    end
   
    % Cut-off \chi(q) beyond the Stoner instability
    if (U ~= 0)
        if Nelm==3
            error('orb1 =/= orb2 case not coded yet!')
        end
        rechis0a = real(chis0(:,:,:,1));
        rechis0b = real(chis0(:,:,:,2));
        snchis0b = sign(rechis0b);
        rechis0b = abs(rechis0b);
        maxchi = U*max(max(rechis0a(:,:,numwi+1)+rechis0b(:,:,numwi+1)));
        if (maxchi > maxchi_cutoff)
            fprintf('  Warning: max U*chis0 = %g,  decrease U*chis0 to U*chi_max = %g\n',maxchi,maxchi_cutoff)
            %{
            idx = find( (rechis0a+rechis0b) > maxchi_cutoff/U);
            chis0(idx) = rechis0a(idx)./(rechis0a(idx) + rechis0b(idx))*(maxchi_cutoff/U);
            chis0(idx+nkwpt) = (maxchi_cutoff/U - chis0(idx)).*snchis0b(idx);
            %}
        end
    end

    % Compute the full spin/charge-fluctuation propagator 
    %   (RPA interactions)
    idxmap = [1 2 2 1];       %orb1 is equivalent to orb2
    %idxmap = [1 2 2 3];       %orb1 is not equivalent to orb2
    if useSymmetry
        chis0 = reshape(chis0(1:(Nk(1)+1),1:(Nk(2)+1),:,idxmap),[],Norb^2);
        chic0 = reshape(chic0(1:(Nk(1)+1),1:(Nk(2)+1),:,idxmap),[],Norb^2);
        siz = [(Nk(1)+1), (Nk(2)+1), 2*numwi, Nelm];
    else
        chis0 = reshape(chis0(:,:,:,idxmap),[],Norb^2);
        chic0 = reshape(chic0(:,:,:,idxmap),[],Norb^2);
        siz = [2*Nk(1), 2*Nk(2), 2*numwi, Nelm];
    end
    %Note: the frequency grid of chi is {-N_c,...,N_c-1}, so for frequency
    % symmetry, the indices [1:(numwi+1)] should be used. Due to this
    % difference, we currently do not use the boson frequency sysmmetry.

    chisR = mul2(inv2(Id - U*chis0),chis0);
    chicR = mul2(inv2(Id + U*chic0),chic0);

    chis0 = reshape(chis0(:,[1 2]),siz);
    chic0 = reshape(chic0(:,[1 2]),siz);
    chisR = reshape(chisR(:,[1 2]),siz);
    chicR = reshape(chicR(:,[1 2]),siz);
    if useSymmetry
        chis0 = chis0([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
        chic0 = chic0([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
        chisR = chisR([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
        chicR = chicR([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    end
    
    %----------------------------------------------------------------------
    % The FLEX equations define the effective interaction as follows,
    % depending on the RPA redefined susceptibilities (minus a correction
    % from the bare susceptibilities)
    Vnm = (1.5*U*U)*chisR + (0.5*U*U)*chicR - (0.5*U*U)*(chis0 + chic0);
    Van = (1.5*U*U)*chisR - (0.5*U*U)*chicR - (0.5*U*U)*(chis0 - chic0);
    
    
    %----------------------------------------------------------------------
    % FORWARD ELECTRON-PHONON SCATTERING
    
    % The electron-phonon coupling consists of:
    % - 'fqphw': the momentum structure
    % - 'invD0': the bare phonon propagator
    % It only couples to the '1' component of the orbital space
    Vnm(:,:,:,1) = Vnm(:,:,:,1) + lam0*fqphw./invD0;
    Vnm(:,:,:,2) = Vnm(:,:,:,2) + lam0*fqphw./invD0;
    % It couples with a MINUS SIGN on the anomalous self-energy
    Van(:,:,:,1) = Van(:,:,:,1) - lam0*fqphw./invD0;
    Van(:,:,:,2) = Van(:,:,:,2) - lam0*fqphw./invD0;
    
    %----------------------------------------------------------------------
    % Forward Fourier transform effective interaction 
    %        (k,i\omega_n) to (r,\tau)
    
    
    for jj = 1:Nelm
        VnmMomentum(:,:,:,jj) = Vnm(:,:,:,jj);
        VanMomentum(:,:,:,jj) = Van(:,:,:,jj);
        Vnm(:,:,:,jj) = fourier_fwd(Vnm(:,:,:,jj),[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
        Van(:,:,:,jj) = fourier_fwd(Van(:,:,:,jj),[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
    end
    % Check whether the effective interaction in (r,\tau) is real
    tmp = abs(imag(Vnm)); err1 = max(tmp(:));
    tmp = abs(imag(Van)); err2 = max(tmp(:));
    if max([err1,err2]) > 1e-10
        fprintf('  Warning: imag(Vnm) is too large! err_imag = %g, %g\n',err1,err2)
    end
    % Impose that they are real
    Vnm = real(Vnm);
    Van = real(Van);
    VnmMomentum = real(VnmMomentum);
    VanMomentum = real(VanMomentum);
    
    
    %----------------------------------------------------------------------
    % Here we compute the self-energy: the normal interaction Vnm and
    % the anomalous interaction Van are convoluted with Gn and Fn,
    % respectively
    % For two-orbital case, with U'=J=J'=0
    S = Vnm.*Gn;
    P = Van.*Fn;
    
    % Backward Fourier transform self-energies (r,\tau) to (k,i\omega_n)
    for jj = 1:Nelm
        % Normal self-energy
        if (jj==1) || (jj==3)
            S(:,:,:,jj) = fourier_bwd(S(:,:,:,jj),-Vnm(1,1,jj),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
        else
            S(:,:,:,jj) = fourier_bwd(S(:,:,:,jj),0,NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
        end
        % Anomalous self-energy
        P(:,:,:,jj) = fourier_bwd(P(:,:,:,jj),0,NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
        % Add Hartree-Fock term
        if (S_HF(jj)~=0)
            S(:,:,:,jj) = S(:,:,:,jj) + S_HF(jj);
        end
        if (P_HF(jj)~=0)
            P(:,:,:,jj) = P(:,:,:,jj) + P_HF(jj);
        end
    end

    %----------------------------------------------------------------------
    
    % If we need to check frequency symmetry
    if checkFreqSymm == 1
        tmp = abs(S(:,:,1:numwi,:) - conj(S(:,:,end:-1:(numwi+1),:))); err1 = max(tmp(:));
        tmp = abs(P(:,:,1:numwi,:) - conj(P(:,:,end:-1:(numwi+1),:))); err2 = max(tmp(:));
        if max([err1,err2]) > 1e-14
            fprintf('  Error of frequency symmetry is too large! [errS errP] =[%g %g]\n',err1,err2)
        end
    end
    
    % Compute the new chemical potential mu and check the filling
    try
        mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
    catch
        errmsg = lasterror; 
        fprintf(errmsg.message)
        fprintf(' Noninteracting mu0 used!\n')
        mu = mu0;
    end

    % Compute the filling again
    filling = get_filling(S,P,WN,ek,mu,xifill0,Norb,beta,useSymmetry);    

    % Compute the gap
    wn = pi/beta;
    for jj = 1:Nelm
        tmp = real(P(:,:,numwi+1,jj));
        [mgap midx] = max(abs(tmp(:)));
        gap(jj) = tmp(midx);
    end

    %----------------------------------------------------------------------
    % Decide if we have to exit
    diffS = 0; diffP = 0;
    tmp = abs(S - Sold);
    diffS = max(tmp(:));
    tmp = abs(P - Pold);
    diffP = max(tmp(:));

    if (diffS < eps && diffP < eps), loop = 0; end;
    fprintf('  %5d  %16.12f  %16.12f',iter,diffS,diffP)
    fprintf('  %16.12f  %16.12f',mu,filling)
    fprintf('  %16.12f  %16.12f\n',gap(1),gap(2))

    % Decide if we need to exit the loop.
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check convergence!\n',maxiter)
        %continue
    end
    
    %----------------------------------------------------------------------
    
    %Anderson acceleration
    if (loop ~= 0) && (mixing == 1) && (mix_method == 2)
        % Compute the current residual norm.
        %xx = real([Sold(:); Pold(:)]);  gval = real([S(:); P(:)]);
        [Gn Fn] = solve_dyson(S,P,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
        if useSymmetry
            xx = [reshape(Gnold(1:(Nk(1)+1),1:(Nk(2)+1),:,:),[],1); ...
                reshape(Fnold(1:(Nk(1)+1),1:(Nk(2)+1),:,:),[],1)];
            gval = [reshape(Gn(1:(Nk(1)+1),1:(Nk(2)+1),:,:),[],1); ...
                reshape(Fn(1:(Nk(1)+1),1:(Nk(2)+1),:,:),[],1)];
        else
            xx = [Gnold(:); Fnold(:)];
            gval = [Gn(:); Fn(:)];
        end
        fval = gval - xx;
        res_norm = max(abs(fval));   % = norm(A,inf)
        %fprintf('     -- %d %g \n', iter, res_norm);
        res_hist = [res_hist;[iter,res_norm]];

        if iter < AAstart %|| mMax == 0 %note we set mMax>0
            % Without acceleration, update x <- g(x) to obtain the next
            % approximate solution.
            Gnold = wgt*Gn + (1-wgt)*Gnold;
            Fnold = wgt*Fn + (1-wgt)*Fnold;
            %find the self-energies from the mixed Green's function           
            [Sold Pold] = solve_dysonbwd(Gnold,Fnold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
            Gn = Gnold;
            Fn = Fnold;
            %xx = gval;
            %recalculate mu after mixing, commented now, see note above
            %mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
        else
            % Apply Anderson acceleration.

            % Update the df vector and the DG array.
            if iter > AAstart
                df = fval-f_old;
                if mAA < mMax
                    DG = [DG gval-g_old];
                else
                    DG = [DG(:,2:mAA) gval-g_old];
                end
                mAA = mAA + 1;
            end
            f_old = fval;
            g_old = gval;

            if mAA == 0     %iter == AAstart
                % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
                fprintf('  %d, start Anderson acceleration and store %d steps of self-energies\n', iter, mMax);
                Gnold = Gn;
                Fnold = Fn;
                %find the self-energies from the mixed Green's function
                [Sold Pold] = solve_dysonbwd(Gnold,Fnold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
                Gn = Gnold;
                Fn = Fnold;
                %xx = gval;
                %recalculate mu after mixing, commented now, see note above
                %mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
            else
                % If mAA > 0, solve the least-squares problem and update the solution.
                if mAA == 1
                    % If mAA == 1, form the initial QR decomposition.
                    R(1,1) = norm(df);
                    Q = R(1,1)\df;
                else
                    % If mAA > 1, update the QR decomposition.
                    if mAA > mMax
                        % If the column dimension of Q is mMax, delete the first column and
                        % update the decomposition.
                        [Q,R] = qrdelete(Q,R,1);
                        mAA = mAA - 1;
                        if size(R,1) ~= size(R,2),
                            Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                        end
                    end
                    % Now update the QR decomposition to incorporate the new column.
                    for j = 1:mAA - 1
                        R(j,mAA) = Q(:,j)'*df;
                        df = df - R(j,mAA)*Q(:,j);
                    end
                    R(mAA,mAA) = norm(df);
                    Q = [Q, R(mAA,mAA)\df];
                end
                if droptol > 0
                    % Drop residuals to improve conditioning if necessary.
                    condDF = cond(R);
                    while condDF > droptol && mAA > 1
                        fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                        [Q,R] = qrdelete(Q,R,1);
                        DG = DG(:,2:mAA);
                        mAA = mAA - 1;
                        % The following treats the qrdelete quirk described above.
                        if size(R,1) ~= size(R,2)
                            Q = Q(:,1:mAA); R = R(1:mAA,:);
                        end
                        condDF = cond(R);
                    end
                end
                % Solve the least-squares problem.
                gamma = R\(Q'*fval);
                % Update the approximate solution.
                xx = gval - DG*gamma;
                % Apply damping if dampbt is a function handle or if dampbt > 0
                % (and dampbt ~= 1).
                if isa(dampbt,'function_handle')
                    xx = xx - (1-dampbt(iter))*(fval - Q*R*gamma);
                else
                    if dampbt > 0 && dampbt ~= 1
                        xx = xx - (1-dampbt)*(fval - Q*R*gamma);
                    end
                end
                if useSymmetry
                    sizg = [(Nk(1)+1), (Nk(2)+1), 2*numwi, Nelm];
                    Gnold = reshape(xx(1:nptmatG),sizg);
                    Gnold = Gnold([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
                    Fnold = reshape(xx(1+nptmatG:nptmatG*2),sizg);
                    Fnold = Fnold([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
                else
                    Gnold(:) = xx(1:nptmatG);
                    Fnold(:) = xx(1+nptmatG:nptmatG*2);
                end
                %find the self-energies from the mixed Green's function
                [Sold Pold] = solve_dysonbwd(Gnold,Fnold,WN,ek,mu,Norb,Nk,numwi,useSymmetry);
                Gn = Gnold;
                Fn = Fnold;                
                %recalculate mu after mixing, commented now, see note above
                %mu = get_mu(S,P,WN,ek,muLim,xifill0,Norb,beta,useSymmetry);
            end
        end %if iter < AAstart
    end %if (mixing == 1) && (mix_method == 2)
end %self-consistency loop
