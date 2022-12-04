clear all
close all
clc

%% Load data

load('samson_1'); clear nBand nCol nRow
V = reshape(V,[156,95,95]); V = permute(V,[2,3,1]);
Z = V; clear V

load('end3'); 
Cref = M; Sref = A'; clear A M cood

load('SRF_S2'); P3 = Pm(:,1:156); clear Pm
d1 = 5; d2 = d1;
[P1,P2] = spatial_deg(Z, 9, d1, d2);

P1 = max(P1,eps); P2 = max(P2,eps); P3 = max(P3,eps);

X_1 = tmprod(Z,{P1,P2},[1,2]); X_1 = max(0,X_1);
X_2 = tmprod(Z,P3,3); X_2 = max(0,X_2);

R = 3;

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
        sigma_h = 10^(-options.SNR_Y1/10); sigma_m = 10^(-options.SNR_Y2/10); 
        opts.lambda = (sigma_h^2)./(sigma_m^2);
        % Poisson noise - Y_2
        Npois= poissrnd(X_2);
        Noise = Npois/norm(Npois(:),'fro');
    end
    propnoise = 1/10^(options.SNR_Y2/20);
    nX = norm(X_2(:),'fro');
    epsi2=nX*propnoise*Noise;
    Y_2 = max(0,X_2+epsi2);

elseif noise_type == 3
    epsi1 = gamrnd(a,b,size(X_1));
    epsi2 = gamrnd(a,b,size(X_2));
    Y_1 = max(0,X_1.*epsi1);
    Y_2 = max(0,X_2.*epsi2);
    opts.lambda = 1;
else
    warning('wrong choice for the noise statistics')
end

P3 = max(0,P3); P1 = max(0,P1); P2 = max(0,P2);

%%

C0 = vca(tens2mat(Y_1,[],3)',R); C0 = sort_endmembers_to_ref(Cref,C0);
C0 = max(C0,eps);
C0_tilde = P3*C0;
S0 = tens2mat(Y_2,[],3)/C0_tilde';


options.beta = 1; %0 = IS, 1 = KL, 2 = Euclidean 
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

prompt = "fprintf('Enter %i values for Lr',R)";
L = input(prompt);
if sum(L)>min(size(Y_1,1),size(Y_2,2))
    warning('The Lr are too large')
end

for r=1:R
    X0 = reshape(S0(:,r),[size(Z,1) size(Z,2)]);
    X0(X0<0) = 1e-25;
    Ainit = rand(size(Y_2,1),L(r));
    Binit = rand(size(Y_2,2),L(r));
    [A0{r},B0{r},cost] = mu_nmf(X0,Ainit,Binit,options);
    figure(r)
    imagesc(A0{r}*B0{r}');
end

A00 = [A0{1} A0{2} A0{3}];
B00 = [B0{1} B0{2} B0{3}];


P10 = rand(size(Y_1,1),size(Y_2,1));
P20 = rand(size(Y_1,1),size(Y_2,1));

for i=1:2
    
    % Update for P1
    V = A00*pw_kronL(C0,P2*B00,R,ones(1,R),L)';
    num = (((P10*V).^(options.beta-2)).*tens2mat(Y_1,1,[]))*V';
    denum = ((P10*V).^(options.beta-1))*V';
    P10 = P10.*((num./denum).^(options.gamma));
    P10 = max(P10,eps);

    % Update for P2
    V = B00*pw_kronL(C0,P1*A00,R,ones(1,R),L)';
    num = (((P20*V).^(options.beta-2)).*tens2mat(Y_1,2,[]))*V';
    denum = ((P20*V).^(options.beta-1))*V';
    P20 = P20.*((num./denum).^(options.gamma));
    P20 = max(P20,eps);

end


%% Algorithm parameters

options.kappa = 1e-10;
options.nIter = 1000;
options.verbose = 1;
options.lambda = 1;

tic;
[A,B,C,P1_est,P2_est,cost] = MU_beta_LL1_2L(Y_1,Y_2,A00,B00,C0,L,P10,P20,P3,options);
t = toc;

%% Show results 1

figure(1)
semilogy(cost)

figure(2)
plot(C)

S = pw_vecL(A,B,R,L);
figure(3)
for r=1:R
    S(:,r) = S(:,r)/norm(S(:,r));
    subplot(1,R,r); imagesc(reshape(S(:,r),[size(Z,1) size(Z,2)]));% caxis([0 1])
end

Zhat = pw_vecL(A,B,R,L)*C'; Zhat = reshape(Zhat,size(Z));

fusion_quality = QualityIndices(Zhat,Z,d1);

frob(P1-P1_est,'squared')/frob(P1,'squared')
frob(P2-P2_est,'squared')/frob(P2,'squared')

frob(Y_1-tmprod(Z,{P1_est,P2_est},[1,2]),'squared')/frob(Y_1,'squared')

figure(1)
subplot(1,3,1); imagesc(P1); 
subplot(1,3,2); imagesc(P1_est);
% subplot(1,3,3); imagesc(P1-P1_est); set(gca,'ColorScale','log')

for r=1:R
    Cref(:,r) = Cref(:,r)/norm(Cref(:,r));
    C(:,r) = C(:,r)/norm(C(:,r));
end
C = sort_endmembers_to_ref(Cref,C);  sam = SAM(Cref,C);
figure(2)
plot(Cref,'k--','Linewidth',1); hold on; plot(C,'r.','Markersize',10); xlim([1 173])
set(gca,'XTick',[]); set(gca,'FontSize',18); 
title(sprintf('SAD = %g',sam))

figure(3)
for r=1:R
    subplot(1,3,r); imagesc(reshape(S(:,r),[size(Z,1) size(Z,2)]));
    set(gca,'XTick',[]); set(gca,'YTick',[]);
end


%% Comparaison méthodes blind

% Blind-STEREO
tic;
[A_hat,B_hat,C,A_tilde,B_tilde,C_tilde]=Blind_TenRec(Y_2,tens2mat(Y_1,[],3),25,50);
[A_hat,B_hat,C_hat,cost] = Blind_STEREO(Y_1,Y_2,P3,10,A_hat,B_hat,A_tilde,B_tilde,C_tilde,1);
t_cp = toc;
Z_cpd = cpdgen({A_hat,B_hat,C});
fusion_quality_cpd = QualityIndices(Z_cpd,Z,d1);

% BSCOTT
%opts.Nblocks = [4 4];
tic;
Z_tuck = bscott(Y_2,Y_1,P3,[45 45 3]);
t_tuck = toc;
fusion_quality_tuck = QualityIndices(Z_tuck,Z,d1);

% BSCOTT
%opts.Nblocks = [4 4];
tic;
Z_tuck = bscott(Y_2,Y_1,P3,[15 15 3]);
t_tuck2 = toc;
fusion_quality_tuck2 = QualityIndices(Z_tuck,Z,d1);

%% Tableau

metrics = cell2mat(struct2cell(fusion_quality))';
metrics_cpd  = cell2mat(struct2cell(fusion_quality_cpd))';
metrics_tuck  = cell2mat(struct2cell(fusion_quality_tuck))';
metrics_tuck2  = cell2mat(struct2cell(fusion_quality_tuck))';

table_blind = [metrics t; metrics_cpd t_cp; metrics_tuck t_tuck; metrics_tuck2 t_tuck2]


