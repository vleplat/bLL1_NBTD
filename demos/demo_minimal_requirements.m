%      Nonnegative block-term decomposition with the beta-divergence:     %
%            joint data fusion and blind spectral unmixing                %
%-------------------------------------------------------------------------%

% Copyright (c) 2022 Clemence Prevost, Valentin Leplat
% https://github.com/cprevost4/bLL1_NBTD
% Contact: clemence.prevost@univ-lille.fr

% This software reproduces the results from the preprint called:
% "Nonnegative block-term decomposition with the beta-divergence:
% joint data fusion and blind spectral unmixing" - C.Prévost, V. Leplat

%-------------------------------------------------------------------------%
%                   DEMO WITH MINIMAL REQUIREMENTS                        %
%-------------------------------------------------------------------------%

clear all
close all
clc

%% Load data

load('end4.mat'); clear cood
M = M(1:173,:);
Z =A'*M'; Z = reshape(Z,[100 100 size(M,1)]);

load('SRF_S2'); P3 = Pm; clear Pm
d1 = 4; d2 = d1; q = 9;
[P1,P2] = spatial_deg(Z, q, d1, d2);

X_1 = tmprod(Z,{P1,P2},[1,2]); X_1 = max(0,X_1);
X_2 = tmprod(Z,P3,3); X_2 = max(0,X_2);

C_ref = M; S_ref = A'; clear A M

R = 4; opts.R = R;

%% Add noise

prompt = "What kind of noise ? 1:Add. Gaussian, 2: Poisson, 3% Gamma multi.";
noise_type = input(prompt);

%-- Gamma noise
noiseGamma_Sigma = 0.05;
% variance for Gamma distribution, mean is set to 1
b = noiseGamma_Sigma^2; 
a = 1/b; 
%-- Gaussian and Poisson noise
options.SNR_Y1 = 30; %dB 
options.SNR_Y2 = 30; %dB


if noise_type == 1 || noise_type == 2
    if noise_type == 1 
        %Gaussian noise - Y_1
        Ngaus = randn(size(X_1));
        Noise = Ngaus/norm(Ngaus(:),'fro');
    elseif noise_type == 2
        % Poisson noise - Y_1
        Npois= poissrnd(X_1);
        Noise = Npois/norm(Npois(:),'fro');
    end
    propnoise = 1/10^(options.SNR_Y1/20);
    nX = norm(X_1(:),'fro');
    epsi1=nX*propnoise*Noise;
    Y_1 = max(0,X_1+epsi1);
    
    if noise_type == 1 
        sigma_h = 10^(-options.SNR_Y1/10); sigma_m = 10^(-options.SNR_Y2/10); 
        opts.lambda = (sigma_h^2)./(sigma_m^2);
        % Gaussian noise - Y_2
        Ngaus = randn(size(X_2));
        Noise = Ngaus/norm(Ngaus(:),'fro');
    elseif noise_type == 2
        % Poisson noise - Y_2
        Npois= poissrnd(X_2);
        Noise = Npois/norm(Npois(:),'fro');
    end
    propnoise = 1/10^(options.SNR_Y2/20);
    nX = norm(X_2(:),'fro');
    epsi2=nX*propnoise*Noise;
    Y_2 = max(0,X_2+epsi2);

elseif noise_type == 3
    opts.lambda = 1;
    epsi1 = gamrnd(a,b,size(X_1));
    epsi2 = gamrnd(a,b,size(X_2));
    Y_1 = max(0,X_1.*epsi1);
    Y_2 = max(0,X_2.*epsi2);
else
    warning('wrong choice for the noise statistics')
    quit
end

P3 = max(0,P3); P1 = max(0,P1); P2 = max(0,P2);

%% Initialization

C0 = vca(tens2mat(Y_1,[],3)',R); C0 = max(C0,eps);
C0_tilde = P3*C0;
S0 = tens2mat(Y_2,[],3)/C0_tilde';

if noise_type == 1 %0 = IS, 1 = KL, 2 = Euclidean 
    options.beta = 3/2;
elseif noise_type == 2
    options.beta = 3/5;
else
    options.beta = 1/4;
end

options.kappa = 1e-10;
options.nIter = 500;
if options.beta <1
    options.gamma = 1/(2-options.beta);
elseif options.beta>2
    options.gamma = 1/(options.beta-1);
else
    options.gamma = 1;
end

L = [24 24 24 24];
for r=1:R
    X0 = reshape(S0(:,r),[size(Z,1) size(Z,2)]);
    X0(X0<0) = 0;
    Ainit = rand(size(Y_2,1),L(r));
    Binit = rand(size(Y_2,2),L(r));
    [A0{r},B0{r},cost] = mu_nmf(X0,Ainit,Binit,options);
%     figure(r)
%     imagesc(A0{r}*B0{r}');
% %         rank(A0{r}*B0{r}')
end

A00 = [A0{1} A0{2} A0{3} A0{4}];
B00 = [B0{1} B0{2} B0{3} B0{4}];

%% Algorithm

options.lambda = 1; %balance parameter between observations
options.kappa = 1e-7;
options.nIter = 1000;
options.verbose = 1;

tic;
[A,B,C,cost] = MU_beta_LL1_1L(Y_1,Y_2,A00,B00,C0,L,P1,P2,P3,options);
toc

%% Results for fusion

S = pw_vecL(A,B,R,L);
Zhat = pw_vecL(A,B,R,L)*C'; Zhat = reshape(Zhat,size(Z));

fprintf('CC: %d \n',cc(Z,Zhat))
fprintf('ERGAS: %d \n',ergas(Z,Zhat,d1,d1))
fprintf('SAM: %d \n',sam(Z,Zhat))
fprintf('RMSE: %d \n',frob(Z-Zhat,'squared')/frob(Z,'squared'))
fprintf('PSNR: %d \n',r_snr(Z,Zhat))

figure(1)
img_ref = Z(:,:,20);
subplot(1,2,1); imagesc(Z(:,:,20)); caxis([min(img_ref(:)) max(img_ref(:))]); colorbar
title('Reference')
subplot(1,2,2); imagesc(Zhat(:,:,20)); caxis([min(img_ref(:)) max(img_ref(:))]); colorbar
title('Estimated')

%% Results for unmixing

for r=1:R
    C_ref(:,r) = C_ref(:,r)/norm(C_ref(:,r));
end

[C, ind] = sort_endmembers_to_ref(C_ref,C);  sam = SAM(C_ref,C);

figure(3)
plot(C_ref,'k--','Linewidth',1); hold on; plot(C,'r.','Markersize',10); xlim([1 173])
set(gca,'XTick',[]); set(gca,'FontSize',18); 
title(sprintf('SAD = %g',sam))


S = S(:,ind);
for r=1:R
    S_ref(:,r) = S_ref(:,r)/norm(S_ref(:,r));
    S(:,r) = S(:,r)/norm(S(:,r));
end

rmse = RMSE(S_ref,S); 

figure(4)
for r=1:R
    subplot(2,2,r); imagesc(reshape(S(:,r),[size(Z,1) size(Z,2)]));
    set(gca,'XTick',[]); set(gca,'YTick',[]);
end
subplot(2,2,1); set(gca,'FontSize',18); title(sprintf('RMSE = %f',rmse))

