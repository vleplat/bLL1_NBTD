function results = run_tensor_benchmark(Y_1,Y_2,P1,P2,P3,opts)


% Initialize the LL1-based algorithms
C0 = vca(tens2mat(Y_1,[],3)',opts.R); C0 = max(C0,eps); C0_tilde = P3*C0;
S0 = tens2mat(Y_2,[],3)/C0_tilde';

options.beta = 1;
options.kappa = 1e-10;
options.nIter = 500;
options.gamma = 1;

L = opts.L;
A00 = []; B00 = []; S00 = [];
for r=1:opts.R
    X0 = reshape(S0(:,r),[size(Y_2,1) size(Y_2,2)]);
    X0(X0<0) = 0;
    Ainit = rand(size(Y_2,1),L);
    Binit = rand(size(Y_2,2),L);
    [A0{r},B0{r}] = mu_nmf(X0,Ainit,Binit,options);
    Abun{r} = A0{r}*B0{r}';
    A00 = [A00 A0{r}]; B00 = [B00 B0{r}]; S00 = [S00 Abun{r}(:)];
end


% STEREO
tic;
[A,B,C,~,~,C_tilde] = TenRec(Y_2,tens2mat(Y_1,[],3),opts.maxit,opts.CPrank,P1,P2);
[A,B,C] = STEREO(Y_1,Y_2,P1,P2,P3,opts.MAXIT,opts.lambda,A,B,C,C_tilde);
I1 = cpdgen({A,B,C}); t1 = toc;

% Blind-STEREO
tic;
[A,B,~,A_tilde,B_tilde,C_tilde]=Blind_TenRec(Y_2,tens2mat(Y_1,[],3),opts.maxit,opts.CPrank);
[A,B,C] = Blind_STEREO(Y_1,Y_2,P3,opts.MAXIT,A,B,A_tilde,B_tilde,C_tilde,opts.lambda);
I2 = cpdgen({A,B,C}); t2 = toc;

% SCOTT
tic; I3 = scott2(Y_1,Y_2, P1, P2, P3, opts.scottrank); t3 = toc;

%BSCOTT
tic; I4 = bscott(Y_2,Y_1,P3,opts.bscottrank,opts); t4 = toc;

%SCLL1
tic;
[S5, C5] = SCLL1(Y_1,Y_2, P1, P2, P3, opts, S00, C0);
I5 = reshape(S5*C5',[size(Y_2,1) size(Y_2,2) size(Y_1,3)]); t5 = toc;

% CT-STAR
tic; I6 = hosvdvar(Y_2,Y_1,opts.CT,opts.Rpsi,P1,P2,P3); t6 = toc;

%CB-STAR
tic; I7 = hosvdvaropt(Y_2,Y_1,P1,P2,P3,opts.CB,opts.Rpsi,1); t7 = toc;

% BTD-Var
tic;
[~,~,S8,C8] = BTD_Var(Y_1,Y_2,P1,P2,P3,opts.R,B00,C0,P3*C0,opts.nIter,opts.lambda);
I8 = reshape(S8*C8',[size(Y_2,1) size(Y_2,2) size(Y_1,3)]); t8 = toc;


results.times = {t1 t2 t3 t4 t5 t6 t7 t8};
results.estimates = {I1 I2 I3 I4 I5 I6 I7 I8};
results.spectra = {C5 C8};
results.abundances = {S5 S8};

end

