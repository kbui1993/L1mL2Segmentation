%%This scripts performs two-phase image segmentation using the AITV models
%%on a synthetic grayscale image with impuse noise.

%% read image
color_synthetic_image2;
f = double(M);
[N,M,~] = size(f);

%store rescaled image
fg = rescale_color_image(f);
fg = double(fg);


result = zeros(5, 11);
result(:,1) = [0.1; 0.2; 0.3; 0.4; 0.5];

for i = 1:5
    
    %set seed
    rng(1234);
    
    %impulse noise; 0 for salt and pepper noise and 1 for random-valued
    %noised
    fnoise = impulsenoise(fg, result(i,1), 1);
    
    %% CV methods
    %set parameters
    pm.outer_iter = 40;
    pm.alpha = 1.0;
    pm.lambda = 2.25;
    pm.c = 1e-8;
    pm.inner_iter = 1000;
    pm.tau = 1/8;
    pm.beta = 1.0;


    %set image segmentation initialization
    u1 = make_circle_shift_x(M,N,30, -5);
    u2 = make_circle_shift_x(M,N, 30, 5);
    u1 = double(u1);
    u2 = double(u2);
    
    %L1-L2
    tic;
    [L1L2_U1,L1L2_U2, L1L2_c1, L1L2_c2, L1L2_c3, L1L2_c4] = L1L2_color_four_phase(fnoise, u1, u2, pm);
    toc
    
    L1L2_c1_im = ones(N,M,3);
    L1L2_c1_im(:,:,1) = ones(N,M)*(L1L2_c1(1)>0.5);
    L1L2_c1_im(:,:,2) = ones(N,M)*(L1L2_c1(2)>0.5);
    L1L2_c1_im(:,:,3) = ones(N,M)*(L1L2_c1(3)>0.5);

    L1L2_c2_im = ones(N,M,3);
    L1L2_c2_im(:,:,1) = ones(N,M)*(L1L2_c2(1)>0.5);
    L1L2_c2_im(:,:,2) = ones(N,M)*(L1L2_c2(2)>0.5);
    L1L2_c2_im(:,:,3) = ones(N,M)*(L1L2_c2(3)>0.5);

    L1L2_c3_im = ones(N,M,3);
    L1L2_c3_im(:,:,1) = ones(N,M)*(L1L2_c3(1)>0.5);
    L1L2_c3_im(:,:,2) = ones(N,M)*(L1L2_c3(2)>0.5);
    L1L2_c3_im(:,:,3) = ones(N,M)*(L1L2_c3(3)>0.5);

    L1L2_c4_im = ones(N,M,3);
    L1L2_c4_im(:,:,1) = ones(N,M)*(L1L2_c4(1)>0.5);
    L1L2_c4_im(:,:,2) = ones(N,M)*(L1L2_c4(2)>0.5);
    L1L2_c4_im(:,:,3) = ones(N,M)*(L1L2_c4(3)>0.5);

    L1L2_approx_im = double(L1L2_U1>0.5).*double(L1L2_U2>0.5).*L1L2_c1_im + double(L1L2_U1>0.5).*double(L1L2_U2<=0.5).*L1L2_c2_im +...
        double(L1L2_U1<=0.5).*double(L1L2_U2>0.5).*L1L2_c3_im + double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5).*L1L2_c4_im;




    %L1-0.75L2
    pm.alpha = 0.75;
    tic;
    [L1L2_75_U1,L1L2_75_U2, L1_0pt75_L2_c1, L1_0pt75_L2_c2, L1_0pt75_L2_c3, L1_0pt75_L2_c4] = L1L2_color_four_phase(fnoise, u1, u2, pm);
    toc

    L1_0pt75_L2_c1_im = ones(N,M,3);
    L1_0pt75_L2_c1_im(:,:,1) = ones(N,M)*(L1_0pt75_L2_c1(1)>0.5);
    L1_0pt75_L2_c1_im(:,:,2) = ones(N,M)*(L1_0pt75_L2_c1(2)>0.5);
    L1_0pt75_L2_c1_im(:,:,3) = ones(N,M)*(L1_0pt75_L2_c1(3)>0.5);

    L1_0pt75_L2_c2_im = ones(N,M,3);
    L1_0pt75_L2_c2_im(:,:,1) = ones(N,M)*(L1_0pt75_L2_c2(1)>0.5);
    L1_0pt75_L2_c2_im(:,:,2) = ones(N,M)*(L1_0pt75_L2_c2(2)>0.5);
    L1_0pt75_L2_c2_im(:,:,3) = ones(N,M)*(L1_0pt75_L2_c2(3)>0.5);

    L1_0pt75_L2_c3_im = ones(N,M,3);
    L1_0pt75_L2_c3_im(:,:,1) = ones(N,M)*(L1_0pt75_L2_c3(1)>0.5);
    L1_0pt75_L2_c3_im(:,:,2) = ones(N,M)*(L1_0pt75_L2_c3(2)>0.5);
    L1_0pt75_L2_c3_im(:,:,3) = ones(N,M)*(L1_0pt75_L2_c3(3)>0.5);

    L1_0pt75_L2_c4_im = ones(N,M,3);
    L1_0pt75_L2_c4_im(:,:,1) = ones(N,M)*(L1_0pt75_L2_c4(1)>0.5);
    L1_0pt75_L2_c4_im(:,:,2) = ones(N,M)*(L1_0pt75_L2_c4(2)>0.5);
    L1_0pt75_L2_c4_im(:,:,3) = ones(N,M)*(L1_0pt75_L2_c4(3)>0.5);

    L1_0pt75_L2_approx_im = double(L1L2_75_U1>0.5).*double(L1L2_75_U2>0.5).*L1_0pt75_L2_c1_im + double(L1L2_75_U1>0.5).*double(L1L2_75_U2<=0.5).*L1_0pt75_L2_c2_im +...
        double(L1L2_75_U1<=0.5).*double(L1L2_75_U2>0.5).*L1_0pt75_L2_c3_im + double(L1L2_75_U1<=0.5).*double(L1L2_75_U2<=0.5).*L1_0pt75_L2_c4_im;

    
    
    %L1-0.5L2
    pm.alpha = 0.5;
    tic;
    [L1L2_05_U1,L1L2_05_U2, L1_0pt5_L2_c1, L1_0pt5_L2_c2, L1_0pt5_L2_c3, L1_0pt5_L2_c4] = L1L2_color_four_phase(fnoise, u1, u2, pm);
    toc

    L1_0pt5_L2_c1_im = ones(N,M,3);
    L1_0pt5_L2_c1_im(:,:,1) = ones(N,M)*(L1_0pt5_L2_c1(1)>0.5);
    L1_0pt5_L2_c1_im(:,:,2) = ones(N,M)*(L1_0pt5_L2_c1(2)>0.5);
    L1_0pt5_L2_c1_im(:,:,3) = ones(N,M)*(L1_0pt5_L2_c1(3)>0.5);

    L1_0pt5_L2_c2_im = ones(N,M,3);
    L1_0pt5_L2_c2_im(:,:,1) = ones(N,M)*(L1_0pt5_L2_c2(1)>0.5);
    L1_0pt5_L2_c2_im(:,:,2) = ones(N,M)*(L1_0pt5_L2_c2(2)>0.5);
    L1_0pt5_L2_c2_im(:,:,3) = ones(N,M)*(L1_0pt5_L2_c2(3)>0.5);

    L1_0pt5_L2_c3_im = ones(N,M,3);
    L1_0pt5_L2_c3_im(:,:,1) = ones(N,M)*(L1_0pt5_L2_c3(1)>0.5);
    L1_0pt5_L2_c3_im(:,:,2) = ones(N,M)*(L1_0pt5_L2_c3(2)>0.5);
    L1_0pt5_L2_c3_im(:,:,3) = ones(N,M)*(L1_0pt5_L2_c3(3)>0.5);

    L1_0pt5_L2_c4_im = ones(N,M,3);
    L1_0pt5_L2_c4_im(:,:,1) = ones(N,M)*(L1_0pt5_L2_c4(1)>0.5);
    L1_0pt5_L2_c4_im(:,:,2) = ones(N,M)*(L1_0pt5_L2_c4(2)>0.5);
    L1_0pt5_L2_c4_im(:,:,3) = ones(N,M)*(L1_0pt5_L2_c4(3)>0.5);

    L1_0pt5_L2_approx_im = double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5).*L1_0pt5_L2_c1_im + double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5).*L1_0pt5_L2_c2_im +...
        double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5).*L1_0pt5_L2_c3_im + double(L1L2_05_U1<=0.5).*double(L1L2_05_U2<=0.5).*L1_0pt5_L2_c4_im;


    %L1-0.25L2
    pm.alpha = 0.25;
    tic;
    [L1L2_25_U1,L1L2_25_U2, L1_0pt25_L2_c1, L1_0pt25_L2_c2, L1_0pt25_L2_c3, L1_0pt25_L2_c4] = L1L2_color_four_phase(fnoise, u1, u2, pm);
    toc

    L1_0pt25_L2_c1_im = ones(N,M,3);
    L1_0pt25_L2_c1_im(:,:,1) = ones(N,M)*(L1_0pt25_L2_c1(1)>0.5);
    L1_0pt25_L2_c1_im(:,:,2) = ones(N,M)*(L1_0pt25_L2_c1(2)>0.5);
    L1_0pt25_L2_c1_im(:,:,3) = ones(N,M)*(L1_0pt25_L2_c1(3)>0.5);

    L1_0pt25_L2_c2_im = ones(N,M,3);
    L1_0pt25_L2_c2_im(:,:,1) = ones(N,M)*(L1_0pt25_L2_c2(1)>0.5);
    L1_0pt25_L2_c2_im(:,:,2) = ones(N,M)*(L1_0pt25_L2_c2(2)>0.5);
    L1_0pt25_L2_c2_im(:,:,3) = ones(N,M)*(L1_0pt25_L2_c2(3)>0.5);

    L1_0pt25_L2_c3_im = ones(N,M,3);
    L1_0pt25_L2_c3_im(:,:,1) = ones(N,M)*(L1_0pt25_L2_c3(1)>0.5);
    L1_0pt25_L2_c3_im(:,:,2) = ones(N,M)*(L1_0pt25_L2_c3(2)>0.5);
    L1_0pt25_L2_c3_im(:,:,3) = ones(N,M)*(L1_0pt25_L2_c3(3)>0.5);

    L1_0pt25_L2_c4_im = ones(N,M,3);
    L1_0pt25_L2_c4_im(:,:,1) = ones(N,M)*(L1_0pt25_L2_c4(1)>0.5);
    L1_0pt25_L2_c4_im(:,:,2) = ones(N,M)*(L1_0pt25_L2_c4(2)>0.5);
    L1_0pt25_L2_c4_im(:,:,3) = ones(N,M)*(L1_0pt25_L2_c4(3)>0.5);

    L1_0pt25_L2_approx_im = double(L1L2_25_U1>0.5).*double(L1L2_25_U2>0.5).*L1_0pt25_L2_c1_im + double(L1L2_25_U1>0.5).*double(L1L2_25_U2<=0.5).*L1_0pt25_L2_c2_im +...
        double(L1L2_25_U1<=0.5).*double(L1L2_25_U2>0.5).*L1_0pt25_L2_c3_im + double(L1L2_25_U1<=0.5).*double(L1L2_25_U2<=0.5).*L1_0pt25_L2_c4_im;


    %anisotropic
    pm.alpha = 0;
    tic;
    [ani_U1,ani_U2, ani_c1, ani_c2, ani_c3, ani_c4] = L1L2_color_four_phase(fnoise, u1, u2, pm);
    toc

    L1_c1_im = ones(N,M,3);
    L1_c1_im(:,:,1) = ones(N,M)*(ani_c1(1)>0.5);
    L1_c1_im(:,:,2) = ones(N,M)*(ani_c1(2)>0.5);
    L1_c1_im(:,:,3) = ones(N,M)*(ani_c1(3)>0.5);

    L1_c2_im = ones(N,M,3);
    L1_c2_im(:,:,1) = ones(N,M)*(ani_c2(1)>0.5);
    L1_c2_im(:,:,2) = ones(N,M)*(ani_c2(2)>0.5);
    L1_c2_im(:,:,3) = ones(N,M)*(ani_c2(3)>0.5);

    L1_c3_im = ones(N,M,3);
    L1_c3_im(:,:,1) = ones(N,M)*(ani_c3(1)>0.5);
    L1_c3_im(:,:,2) = ones(N,M)*(ani_c3(2)>0.5);
    L1_c3_im(:,:,3) = ones(N,M)*(ani_c3(3)>0.5);

    L1_c4_im = ones(N,M,3);
    L1_c4_im(:,:,1) = ones(N,M)*(ani_c4(1)>0.5);
    L1_c4_im(:,:,2) = ones(N,M)*(ani_c4(2)>0.5);
    L1_c4_im(:,:,3) = ones(N,M)*(ani_c4(3)>0.5);

    L1_approx_im = double(ani_U1>0.5).*double(ani_U2>0.5).*L1_c1_im + double(ani_U1>0.5).*double(ani_U2<=0.5).*L1_c2_im +...
        double(ani_U1<=0.5).*double(ani_U2>0.5).*L1_c3_im + double(ani_U1<=0.5).*double(ani_U2<=0.5).*L1_c4_im;

    

    %% fuzzy method
    %set parameters
    pm2.outer_iter = 160;
    pm2.alpha = 1.0;
    pm2.lambda = 2.25;
    pm2.c = 1e-8;
    pm2.inner_iter = 1000;
    pm2.tau = 1/8;
    pm2.beta = 1.0;
    pm2.nu = 5.0;

    u_initial{1} = u1.*u2;
    u_initial{2} = u1.*(1-u2);
    u_initial{3} = (1-u1).*u2;
    u_initial{4} = (1-u1).*(1-u2);



    %L1-1.0L2
    tic;
    [L1_L2_f, L1_L2_c] = fuzzy_color_L1L2(fnoise, u_initial, pm2, 4);
    toc
    L1_L2_c{1} = L1_L2_c{1}>0.5;
    L1_L2_c{2} = L1_L2_c{2}>0.5;
    L1_L2_c{3} = L1_L2_c{3}>0.5;
    L1_L2_c{4} = L1_L2_c{4}>0.5;
    
    L1_L2_approx1 =...
    double((L1_L2_f{1}>L1_L2_f{2}).*(L1_L2_f{1}>L1_L2_f{3}).*(L1_L2_f{1}>L1_L2_f{4})*L1_L2_c{1}(1))+...
    double((L1_L2_f{2}>L1_L2_f{1}).*(L1_L2_f{2}>L1_L2_f{3}).*(L1_L2_f{2}>L1_L2_f{4})*L1_L2_c{2}(1))+...
    double((L1_L2_f{3}>L1_L2_f{2}).*(L1_L2_f{3}>L1_L2_f{1}).*(L1_L2_f{3}>L1_L2_f{4})*L1_L2_c{3}(1))+...
    double((L1_L2_f{4}>L1_L2_f{2}).*(L1_L2_f{4}>L1_L2_f{3}).*(L1_L2_f{4}>L1_L2_f{1})*L1_L2_c{4}(1));
    L1_L2_approx2 =...
        double((L1_L2_f{1}>L1_L2_f{2}).*(L1_L2_f{1}>L1_L2_f{3}).*(L1_L2_f{1}>L1_L2_f{4})*L1_L2_c{1}(2))+...
        double((L1_L2_f{2}>L1_L2_f{1}).*(L1_L2_f{2}>L1_L2_f{3}).*(L1_L2_f{2}>L1_L2_f{4})*L1_L2_c{2}(2))+...
        double((L1_L2_f{3}>L1_L2_f{2}).*(L1_L2_f{3}>L1_L2_f{1}).*(L1_L2_f{3}>L1_L2_f{4})*L1_L2_c{3}(2))+...
        double((L1_L2_f{4}>L1_L2_f{2}).*(L1_L2_f{4}>L1_L2_f{3}).*(L1_L2_f{4}>L1_L2_f{1})*L1_L2_c{4}(2));
    L1_L2_approx3 =...
        double((L1_L2_f{1}>L1_L2_f{2}).*(L1_L2_f{1}>L1_L2_f{3}).*(L1_L2_f{1}>L1_L2_f{4})*L1_L2_c{1}(3))+...
        double((L1_L2_f{2}>L1_L2_f{1}).*(L1_L2_f{2}>L1_L2_f{3}).*(L1_L2_f{2}>L1_L2_f{4})*L1_L2_c{2}(3))+...
        double((L1_L2_f{3}>L1_L2_f{2}).*(L1_L2_f{3}>L1_L2_f{1}).*(L1_L2_f{3}>L1_L2_f{4})*L1_L2_c{3}(3))+...
        double((L1_L2_f{4}>L1_L2_f{2}).*(L1_L2_f{4}>L1_L2_f{3}).*(L1_L2_f{4}>L1_L2_f{1})*L1_L2_c{4}(3));
    L1_L2_approx = zeros(N,M,3);
    L1_L2_approx(:,:,1) = L1_L2_approx1;
    L1_L2_approx(:,:,2) = L1_L2_approx2;
    L1_L2_approx(:,:,3) = L1_L2_approx3;
    

    %L1-0.75L2
    pm2.alpha = 0.75;
    tic;
    [L1_0pt75_L2_f, L1_0pt75_L2_c] = fuzzy_color_L1L2(fnoise, u_initial, pm2, 4);
    toc
    
    L1_0pt75_L2_c{1} = L1_0pt75_L2_c{1}>0.5;
    L1_0pt75_L2_c{2} = L1_0pt75_L2_c{2}>0.5;
    L1_0pt75_L2_c{3} = L1_0pt75_L2_c{3}>0.5;
    L1_0pt75_L2_c{4} = L1_0pt75_L2_c{4}>0.5;
    
    L1_0pt75_L2_approx1 =...
    double((L1_0pt75_L2_f{1}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{1}(1))+...
    double((L1_0pt75_L2_f{2}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{2}(1))+...
    double((L1_0pt75_L2_f{3}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{3}(1))+...
    double((L1_0pt75_L2_f{4}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{1})*L1_0pt75_L2_c{4}(1));
    L1_0pt75_L2_approx2 =...
        double((L1_0pt75_L2_f{1}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{1}(2))+...
        double((L1_0pt75_L2_f{2}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{2}(2))+...
        double((L1_0pt75_L2_f{3}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{3}(2))+...
        double((L1_0pt75_L2_f{4}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{1})*L1_0pt75_L2_c{4}(2));
    L1_0pt75_L2_approx3 =...
        double((L1_0pt75_L2_f{1}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{1}(3))+...
        double((L1_0pt75_L2_f{2}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{2}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{2}(3))+...
        double((L1_0pt75_L2_f{3}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{3}>L1_0pt75_L2_f{4})*L1_0pt75_L2_c{3}(3))+...
        double((L1_0pt75_L2_f{4}>L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{4}>L1_0pt75_L2_f{1})*L1_0pt75_L2_c{4}(3));

    L1_0pt75_L2_approx = zeros(N,M,3);
    L1_0pt75_L2_approx(:,:,1) = L1_0pt75_L2_approx1;
    L1_0pt75_L2_approx(:,:,2) = L1_0pt75_L2_approx2;
    L1_0pt75_L2_approx(:,:,3) = L1_0pt75_L2_approx3;

    %L1-0.5L2
    pm2.alpha = 0.5;
    tic;
    [L1_0pt5_L2_f, L1_0pt5_L2_c]  = fuzzy_color_L1L2(fnoise, u_initial, pm2, 4);
    toc
    
        L1_0pt5_L2_c{1} = L1_0pt5_L2_c{1}>0.5;
    L1_0pt5_L2_c{2} = L1_0pt5_L2_c{2}>0.5;
    L1_0pt5_L2_c{3} = L1_0pt5_L2_c{3}>0.5;
    L1_0pt5_L2_c{4} = L1_0pt5_L2_c{4}>0.5;
    
    
    L1_0pt5_L2_approx1 =...
    double((L1_0pt5_L2_f{1}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{1}(1))+...
    double((L1_0pt5_L2_f{2}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{2}(1))+...
    double((L1_0pt5_L2_f{3}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{3}(1))+...
    double((L1_0pt5_L2_f{4}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{1})*L1_0pt5_L2_c{4}(1));
    L1_0pt5_L2_approx2 =...
        double((L1_0pt5_L2_f{1}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{1}(2))+...
        double((L1_0pt5_L2_f{2}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{2}(2))+...
        double((L1_0pt5_L2_f{3}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{3}(2))+...
        double((L1_0pt5_L2_f{4}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{1})*L1_0pt5_L2_c{4}(2));
    L1_0pt5_L2_approx3 =...
        double((L1_0pt5_L2_f{1}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{1}(3))+...
        double((L1_0pt5_L2_f{2}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{2}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{2}(3))+...
        double((L1_0pt5_L2_f{3}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{3}>L1_0pt5_L2_f{4})*L1_0pt5_L2_c{3}(3))+...
        double((L1_0pt5_L2_f{4}>L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{4}>L1_0pt5_L2_f{1})*L1_0pt5_L2_c{4}(3));
    L1_0pt5_L2_approx = zeros(N,M,3);
    L1_0pt5_L2_approx(:,:,1) = L1_0pt5_L2_approx1;
    L1_0pt5_L2_approx(:,:,2) = L1_0pt5_L2_approx2;
    L1_0pt5_L2_approx(:,:,3) = L1_0pt5_L2_approx3;

    
    %L1-0.25L2
    pm2.alpha = 0.25;
    tic;
    [L1_0pt25_L2_f, L1_0pt25_L2_c] = fuzzy_color_L1L2(fnoise, u_initial, pm2, 4);
    toc
    
        L1_0pt25_L2_c{1} = L1_0pt25_L2_c{1}>0.5;
    L1_0pt25_L2_c{2} = L1_0pt25_L2_c{2}>0.5;
    L1_0pt25_L2_c{3} = L1_0pt25_L2_c{3}>0.5;
    L1_0pt25_L2_c{4} = L1_0pt25_L2_c{4}>0.5;
    
    L1_0pt25_L2_approx1 =...
    double((L1_0pt25_L2_f{1}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{1}(1))+...
    double((L1_0pt25_L2_f{2}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{2}(1))+...
    double((L1_0pt25_L2_f{3}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{3}(1))+...
    double((L1_0pt25_L2_f{4}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{1})*L1_0pt25_L2_c{4}(1));
    L1_0pt25_L2_approx2 =...
        double((L1_0pt25_L2_f{1}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{1}(2))+...
        double((L1_0pt25_L2_f{2}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{2}(2))+...
        double((L1_0pt25_L2_f{3}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{3}(2))+...
        double((L1_0pt25_L2_f{4}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{1})*L1_0pt25_L2_c{4}(2));
    L1_0pt25_L2_approx3 =...
        double((L1_0pt25_L2_f{1}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{1}(3))+...
        double((L1_0pt25_L2_f{2}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{2}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{2}(3))+...
        double((L1_0pt25_L2_f{3}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{3}>L1_0pt25_L2_f{4})*L1_0pt25_L2_c{3}(3))+...
        double((L1_0pt25_L2_f{4}>L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{4}>L1_0pt25_L2_f{1})*L1_0pt25_L2_c{4}(3));

    L1_0pt25_L2_approx = zeros(N,M,3);
    L1_0pt25_L2_approx(:,:,1) = L1_0pt25_L2_approx1;
    L1_0pt25_L2_approx(:,:,2) = L1_0pt25_L2_approx2;
    L1_0pt25_L2_approx(:,:,3) = L1_0pt25_L2_approx3;

    %L1
    pm2.alpha = 0;
    tic;
    [L1_f, L1_c] = fuzzy_color_L1L2(fnoise, u_initial, pm2, 4);
    toc
    
           L1_c{1} = L1_c{1}>0.5;
    L1_c{2} = L1_c{2}>0.5;
    L1_c{3} = L1_c{3}>0.5;
    L1_c{4} = L1_c{4}>0.5;

    L1_approx1 =...
    double((L1_f{1}>L1_f{2}).*(L1_f{1}>L1_f{3}).*(L1_f{1}>L1_f{4})*L1_c{1}(1))+...
    double((L1_f{2}>L1_f{1}).*(L1_f{2}>L1_f{3}).*(L1_f{2}>L1_f{4})*L1_c{2}(1))+...
    double((L1_f{3}>L1_f{2}).*(L1_f{3}>L1_f{1}).*(L1_f{3}>L1_f{4})*L1_c{3}(1))+...
    double((L1_f{4}>L1_f{2}).*(L1_f{4}>L1_f{3}).*(L1_f{4}>L1_f{1})*L1_c{4}(1));
    L1_approx2 =...
        double((L1_f{1}>L1_f{2}).*(L1_f{1}>L1_f{3}).*(L1_f{1}>L1_f{4})*L1_c{1}(2))+...
        double((L1_f{2}>L1_f{1}).*(L1_f{2}>L1_f{3}).*(L1_f{2}>L1_f{4})*L1_c{2}(2))+...
        double((L1_f{3}>L1_f{2}).*(L1_f{3}>L1_f{1}).*(L1_f{3}>L1_f{4})*L1_c{3}(2))+...
        double((L1_f{4}>L1_f{2}).*(L1_f{4}>L1_f{3}).*(L1_f{4}>L1_f{1})*L1_c{4}(2));
    L1_approx3 =...
        double((L1_f{1}>L1_f{2}).*(L1_f{1}>L1_f{3}).*(L1_f{1}>L1_f{4})*L1_c{1}(3))+...
        double((L1_f{2}>L1_f{1}).*(L1_f{2}>L1_f{3}).*(L1_f{2}>L1_f{4})*L1_c{2}(3))+...
        double((L1_f{3}>L1_f{2}).*(L1_f{3}>L1_f{1}).*(L1_f{3}>L1_f{4})*L1_c{3}(3))+...
        double((L1_f{4}>L1_f{2}).*(L1_f{4}>L1_f{3}).*(L1_f{4}>L1_f{1})*L1_c{4}(3));
    L1_approx = zeros(N,M,3);
    L1_approx(:,:,1) = L1_approx1;
    L1_approx(:,:,2) = L1_approx2;
    L1_approx(:,:,3) = L1_approx3;

    %% create classification image
    fg1 = ones(n,n)+double(fg(:,:,1)==1)+2*double(fg(:,:,2)==1)+3*double(fg(:,:,3)==1);
    fg1 = uint8(fg1);
    fg1 = double(fg1);

    %L1-L2
    L1_L2_M = ones(n,n)+3*(L1L2_U1.*L1L2_U2)+(L1L2_U1.*(1-L1L2_U2))+2*(1-L1L2_U1).*(L1L2_U2);
    L1_L2_M = uint8(L1_L2_M);
    L1_L2_M = double(L1_L2_M);

    %L1-0.75L2
    L1_0pt75_L2_M = ones(n,n)+3*(L1L2_75_U1.*L1L2_75_U2)+(L1L2_75_U1.*(1-L1L2_75_U2))+2*(1-L1L2_75_U1).*(L1L2_75_U2);
    L1_0pt75_L2_M = uint8(L1_0pt75_L2_M);
    L1_0pt75_L2_M = double(L1_0pt75_L2_M);

    %L1-0.5L2
    L1_0pt5_L2_M = ones(n,n)+3*(L1L2_05_U1.*L1L2_05_U2)+(L1L2_05_U1.*(1-L1L2_05_U2))+2*(1-L1L2_05_U1).*(L1L2_05_U2);
    L1_0pt5_L2_M = uint8(L1_0pt5_L2_M);
    L1_0pt5_L2_M = double(L1_0pt5_L2_M);

    %L1-0.25L2
    L1_0pt25_L2_M = ones(n,n)+3*(L1L2_25_U1.*L1L2_25_U2)+(L1L2_25_U1.*(1-L1L2_25_U2))+2*(1-L1L2_25_U1).*(L1L2_25_U2);
    L1_0pt25_L2_M = uint8(L1_0pt25_L2_M);
    L1_0pt25_L2_M = double(L1_0pt25_L2_M);

    %L1
    L1_M = ones(n,n)+3*(ani_U1.*ani_U2)+(ani_U1.*(1-ani_U2))+2*(1-ani_U1).*(ani_U2);
    L1_M = uint8(L1_M);
    L1_M = double(L1_M);

    %fuzzy L1-L2
    fuzzy_L1_L2_M = ones(n,n)+3*double((L1_L2_f{1} >= L1_L2_f{2}).*(L1_L2_f{1} >= L1_L2_f{3}).*(L1_L2_f{1} >= L1_L2_f{4}))+...
        +double((L1_L2_f{2} >= L1_L2_f{1}).*(L1_L2_f{2} >= L1_L2_f{3}).*(L1_L2_f{2} >= L1_L2_f{4}))+...
        +2*double((L1_L2_f{3} >= L1_L2_f{1}).*(L1_L2_f{3} >= L1_L2_f{2}).*(L1_L2_f{3} >= L1_L2_f{4}));

    %fuzzy L1-0.75L2
    fuzzy_L1_0pt75_L2_M = ones(n,n)+3*double((L1_0pt75_L2_f{1} >= L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{1} >= L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{1} >= L1_0pt75_L2_f{4}))+...
        +double((L1_0pt75_L2_f{2} >= L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{2} >= L1_0pt75_L2_f{3}).*(L1_0pt75_L2_f{2} >= L1_0pt75_L2_f{4}))+...
        +2*double((L1_0pt75_L2_f{3} >= L1_0pt75_L2_f{1}).*(L1_0pt75_L2_f{3} >= L1_0pt75_L2_f{2}).*(L1_0pt75_L2_f{3} >= L1_0pt75_L2_f{4}));

    %fuzzy L1-0.5L2
    fuzzy_L1_0pt5_L2_M = ones(n,n)+3*double((L1_0pt5_L2_f{1} >= L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{1} >= L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{1} >= L1_0pt5_L2_f{4}))+...
        +double((L1_0pt5_L2_f{2} >= L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{2} >= L1_0pt5_L2_f{3}).*(L1_0pt5_L2_f{2} >= L1_0pt5_L2_f{4}))+...
        +2*double((L1_0pt5_L2_f{3} >= L1_0pt5_L2_f{1}).*(L1_0pt5_L2_f{3} >= L1_0pt5_L2_f{2}).*(L1_0pt5_L2_f{3} >= L1_0pt5_L2_f{4}));

    %fuzzy L1-0.25L2
    fuzzy_L1_0pt25_L2_M = ones(n,n)+3*double((L1_0pt25_L2_f{1} >= L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{1} >= L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{1} >= L1_0pt25_L2_f{4}))+...
        +double((L1_0pt25_L2_f{2} >= L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{2} >= L1_0pt25_L2_f{3}).*(L1_0pt25_L2_f{2} >= L1_0pt25_L2_f{4}))+...
        +2*double((L1_0pt25_L2_f{3} >= L1_0pt25_L2_f{1}).*(L1_0pt25_L2_f{3} >= L1_0pt25_L2_f{2}).*(L1_0pt25_L2_f{3} >= L1_0pt25_L2_f{4}));

    %fuzzy L1
    fuzzy_L1_M = ones(n,n)+3*double((L1_f{1} >= L1_f{2}).*(L1_f{1} >= L1_f{3}).*(L1_f{1} >= L1_f{4}))+...
        +double((L1_f{2} >= L1_f{1}).*(L1_f{2} >= L1_f{3}).*(L1_f{2} >= L1_f{4}))+...
        +2*double((L1_f{3} >= L1_f{1}).*(L1_f{3} >= L1_f{2}).*(L1_f{3} >= L1_f{4}));
    
    %% compute DICE
    mean(dice(L1_L2_M, fg1))
    mean(dice(L1_0pt75_L2_M, fg1))
    mean(dice(L1_0pt5_L2_M, fg1))
    mean(dice(L1_0pt25_L2_M, fg1))
    mean(dice(L1_M, fg1))
    mean(dice(fuzzy_L1_L2_M, fg1))
    mean(dice(fuzzy_L1_0pt75_L2_M, fg1))
    mean(dice(fuzzy_L1_0pt5_L2_M, fg1))
    mean(dice(fuzzy_L1_0pt25_L2_M, fg1))
    mean(dice(fuzzy_L1_M, fg1))
    
    result(i,2:11) = [mean(dice(L1_L2_M, fg1)),...
    mean(dice(L1_0pt75_L2_M, fg1)),...
    mean(dice(L1_0pt5_L2_M, fg1)),...
    mean(dice(L1_0pt25_L2_M, fg1)),...
    mean(dice(L1_M, fg1)),...
    mean(dice(fuzzy_L1_L2_M, fg1)),...
    mean(dice(fuzzy_L1_0pt75_L2_M, fg1)),...
    mean(dice(fuzzy_L1_0pt5_L2_M, fg1)),...
    mean(dice(fuzzy_L1_0pt25_L2_M, fg1)),...
    mean(dice(fuzzy_L1_M, fg1))];
    
end
