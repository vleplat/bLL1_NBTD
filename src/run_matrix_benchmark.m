function results = run_matrix_benchmark(Y_1,Y_2,P1,P2,P3,opts)

%CNMF
tic; [I1,C1,~,~,S1] = CNMF_fusion(Y_1,Y_2,opts.R); t1 = toc;

% FUSE
tic; I2 = FUSE_wrapper(Y_1,Y_2,opts.d1,opts.lambda,P3,gauss_kernel(opts.q)); t2 = toc;

% HySure
tic; I3 = hysure_adaptor(Y_1,Y_2, P1, P2, P3, opts.R, opts); t3 = toc;

%SFIM
tic; I4 = SFIM(Y_1,Y_2,1); t4 = toc;

results.times = {t1 t2 t3 t4};
results.estimates = {I1 I2 I3 I4};
results.spectra = {C1};
results.abundances = {S1};

end

