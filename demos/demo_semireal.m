clear all
close all
clc

%% Generate data

load('end4.mat'); clear cood
M = M(1:173,:);
Z =A'*M'; Z = reshape(Z,[100 100 size(M,1)]);

load('SRF_S2'); P3 = Pm; clear Pm
d1 = 4; d2 = d1; q = 9;
[P1,P2] = spatial_deg(Z, q, d1, d2);

X_1 = tmprod(Z,{P1,P2},[1,2]); X_1 = max(0,X_1);
X_2 = tmprod(Z,P3,3); X_2 = max(0,X_2);

C_ref = M; S_ref = A'; clear A M

R = 4;

nTrials = 1;

prompt = "What kind of noise ? 1:Add. Gaussian, 2: Poisson, 3% Gamma multi.";
noise_type = input(prompt);


for l=1:nTrials
   
    
    %% Add noise

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
            options.lambda = (sigma_h^2)./(sigma_m^2);
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
        epsi1 = gamrnd(a,b,size(X_1));
        epsi2 = gamrnd(a,b,size(X_2));
        Y_1 = max(0,X_1.*epsi1);
        Y_2 = max(0,X_2.*epsi2);
    else
        warning('wrong choice for the noise statistics')
        quit
    end

    P3 = max(0,P3); P1 = max(0,P1); P2 = max(0,P2);
    %opts.lambda = options.lambda;
    
    %% Proposed

    fprintf('Running proposed approach ... \n')
    C0 = vca(tens2mat(Y_1,[],3)',R); C0 = sort_endmembers_to_ref(C_ref,C0);
    C0 = max(C0,eps);
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

    L = [24 24 24 24 ];
    for r=1:R
        X0 = reshape(S0(:,r),[size(Z,1) size(Z,2)]);
        X0(X0<0) = 1e-25;
        Ainit = rand(size(Y_2,1),L(r));
        Binit = rand(size(Y_2,2),L(r));
        [A0{r},B0{r},cost] = mu_nmf(X0,Ainit,Binit,options);
    end

    A00 = [A0{1} A0{2} A0{3} A0{4}];
    B00 = [B0{1} B0{2} B0{3} B0{4}];

    options.kappa = 1e-7;
    options.nIter = 1000; options.verbose = 1;
    options.lambda = 1;

    prompt = "What kind of processing unit ? 1:cpu, 2: gpu.";
    unit_type = input(prompt);
    if unit_type == 1
        tic;
        [A,B,C,cost] = MU_beta_LL1_1L(Y_1,Y_2,A00,B00,C0,L,P1,P2,P3,options);
        time = toc;
    elseif unit_type == 2
        tic;
        [A,B,C,cost] = MU_beta_LL1_1L_gpu(Y_1,Y_2,A00,B00,C0,L,P1,P2,P3,options);
        time = toc;
    end

    S = pw_vecL(A,B,R,L);
    Zhat = pw_vecL(A,B,R,L)*C'; Zhat = reshape(Zhat,size(Z));

    metrics = QualityIndices(Zhat,Z,d1);
    table_proposed = [transpose(struct2cell(metrics)) time];
    
    all_proposed(:,:,l) = cell2mat(table_proposed);
    
    %% Matrix benchmark

    fprintf('Running matrix benchmark ... \n')
    opts.R = 4;
    opts.P = 30;
    opts.q = q; opts.d1 = d1;

    results = run_matrix_benchmark(Y_1,Y_2,P1,P2,P3,opts);

    table_HSMSFusion = [];
    nMethods = length(cell2mat(results.times));

    for n=1:nMethods
        estimate = results.estimates{n};
        eval(sprintf('metrics%d = QualityIndices(estimate,Z,d1);',n))
        table_HSMSFusion = [table_HSMSFusion; eval(sprintf('[transpose(struct2cell(metrics%d)) results.times{%d}]',n,n))];
    end

    all_matrix(:,:,l) = cell2mat(table_HSMSFusion);
    
    %% Tensor benchmark

    fprintf('Running matrix benchmark ... \n')
    opts.CPrank = 50; opts.maxit = 25; opts.MAXIT = 10; 

    opts.scottrank = [40 40 6]; opts.Nblocks = [4 4]; opts.bscottrank = [8 8 6];

    opts.CT = [6 6 6]; opts.CB = opts.scottrank; opts.Rpsi = [4 4 2];

    opts.nIter = 20; opts.innerIter = 5; opts.rho = 1e-4; opts.mu = 1; opts.mu2 = 1;

    opts.L = 3; 

    opts.th = 1e-4; opts.maxiter = 10; opts.theta = 1e-4; opts.eta = 5*1e-3;


    results = run_tensor_benchmark(Y_1,Y_2,P1,P2,P3,opts);

    table_TensorFusion = [];
    nMethods = length(cell2mat(results.times));

    for n=1:nMethods
        estimate = results.estimates{n};
        eval(sprintf('metrics%d = QualityIndices(estimate,Z,d1);',n))
        table_TensorFusion = [table_TensorFusion; eval(sprintf('[transpose(struct2cell(metrics%d)) results.times{%d}]',n,n))];
    end
    
    all_tensor(:,:,l) = cell2mat(table_TensorFusion);
   
end

%% Fusion results

mean_proposed = mean(all_proposed,3);

mean_matrix = mean(all_matrix,3);

mean_tensor = mean(all_tensor,3);

table_fusion = [mean_proposed; mean_matrix; mean_tensor];
table_fusion = [["Proposed"; "CNMF"; "FUSE"; "HySure"; "SFIM"; "STEREO"; "Blind-STEREO"; "SCOTT"; "BSCOTT"; "SCLL1"; "CT-STAR"; "CB-STAR"; "CNN-BTD-Var"] table_fusion];
table_fusion = [["Method" "CC" "SAM" "RMSE" "ERGAS" "PSNR" "Time"]; table_fusion]


%% Results of unmixing


for r=1:R
    C_ref(:,r) = C_ref(:,r)/norm(C_ref(:,r));
    C(:,r) = C(:,r)/norm(C(:,r));
end

[C, ind] = sort_endmembers_to_ref(C_ref,C);  sam = SAM(C_ref,C);

figure(1)
 plot(C_ref,'k--'); hold on; plot(C,'r.'); xlim([1 173])
set(gca,'XTick',[]); set(gca,'FontSize',18); xlabel('(a)','Interpreter','latex')
title(sprintf('SAD = %g',sam))

S = pw_vecL(A,B,R,L); 
S = S(:,ind); 
for r=1:R
    S_ref(:,r) = S_ref(:,r)/norm(S_ref(:,r));
    S(:,r) = S(:,r)/norm(S(:,r));
end

rmse = RMSE(S_ref,S); 

figure(2)
for r=1:R
    subplot(2,2,r); imagesc(reshape(S(:,r),[size(Z,1) size(Z,2)]));
    set(gca,'XTick',[]); set(gca,'YTick',[]);
end
subplot(2,2,1); set(gca,'FontSize',18); title(sprintf('RMSE = %f',rmse))