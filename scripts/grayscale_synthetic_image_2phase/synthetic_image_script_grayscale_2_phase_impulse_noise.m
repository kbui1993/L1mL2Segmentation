%%This scripts performs two-phase image segmentation using the AITV models
%%on a synthetic grayscale image with impulse noise.

%% read image
fui8 = imread('shape.png');
fui8 = rgb2gray(fui8);
f = double(fui8);
[N,M] = size(f);

%create matrix to store DICE indies
result = zeros(7, 11);
result(:,1) = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7];

for i = 1:7
    %set seed
    rng(1234);
    
    %impulse noise, 0 for salt & pepper noise; 1 for random-valued
    fnoise = impulsenoise(f,result(i,1), 1);
    
    %rescale image
    fg = double(fnoise);
    fg = rescale_image(fg);
    
    %% Chan Vese Models
    %set parameters
    pm.outer_iter = 20;
    pm.alpha = 1.0;
    
    pm.lambda = 2.0;
    pm.c = 1e-8;
    pm.inner_iter = 300;
    pm.tau = 1/8;
    pm.sigma = 1/8;
    pm.method = 'PDHGLS';
    pm.beta = 1;

    %preinitialize segmentation
    u = make_circle(M,N,10);

    u = double(u);

    %%CV method
    %L1-1.0L2
    tic;
    L1_L2_u1 = L1L2_two_phase(fg, u, pm);
    toc

    %L1 - 0.75L2
    pm.alpha = 0.75;
    tic;
    L1_0pt75_L2_u1 = L1L2_two_phase(fg, u, pm);
    toc

    %L1 - 0.5L2
    pm.alpha = 0.5;
    tic;
    L1_0pt5_L2_u1 = L1L2_two_phase(fg, u, pm);
    toc

    %L1 - 0.25L2
    pm.alpha = 0.25;
    tic;
    L1_0pt25_L2_u1 = L1L2_two_phase(fg, u, pm);
    toc

    %anisotropic CV
    pm.alpha =0;
    tic;
    L1_u1 = L1L2_two_phase(fg, u, pm);
    toc

    %% run fuzzy region competition
    pm2.outer_iter = 40;
    pm2.alpha = 1.0;
    
    pm2.lambda = 2.0;
    pm2.c = 1e-8;
    pm2.inner_iter = 300;
    pm2.tau = 1/8;
    pm2.sigma = 1/8;
    pm2.method = 'PDHGLS';
    pm2.nu = 10;
    pm2.beta = 1;
    
    %fuzzy initialization
    u_initial{1} = u;
    u_initial{2} = ones(M,N)-u;
    
    
    %L1-1.0L2
    tic;
    L1_L2_f = fuzzy_L1L2(fg, u_initial, pm2, 2);
    toc
    
    %L1-0.75L2
    pm2.alpha = 0.75;
    tic;
    L1_0pt75_L2_f = fuzzy_L1L2(fg, u_initial, pm2, 2);
    toc
    
    %L1-0.5L2
    pm2.alpha = 0.5;
    tic;
    L1_0pt5_L2_f = fuzzy_L1L2(fg, u_initial, pm2, 2);
    toc
    
    %L1-0.25L2
    pm2.alpha = 0.25;
    tic;
    L1_0pt25_L2_f = fuzzy_L1L2(fg, u_initial, pm2, 2);
    toc
    
    %L1
    pm2.alpha = 0;
    tic;
    L1_f = fuzzy_L1L2(fg, u_initial, pm2, 2);
    toc

    fg1 = rescale_image(f);

    %% compute dice
    result(i,2:11) = [dice(double(L1_L2_u1>0.5), fg1),...
    dice(double(L1_0pt75_L2_u1>0.5), fg1),...
    dice(double(L1_0pt5_L2_u1>0.5), fg1),...
    dice(double(L1_0pt25_L2_u1>0.5), fg1),...
    dice(double(L1_u1>0.5), fg1),...
    dice(double(L1_L2_f{1}>L1_L2_f{2}), fg1),...
    dice(double(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{2}), fg1),...
    dice(double(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{2}), fg1),...
    dice(double(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{2}), fg1),...
    dice(double(L1_f{1}>L1_f{2}), fg1)];
end