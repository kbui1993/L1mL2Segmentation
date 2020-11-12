%%This script performs four-phase segmentation on a color image of a
%%butterfly.

%read butterfly image
fui8 = imread('butterfly.jpg');

%convert to double
f = double(fui8);

%get image size
[N,M,~] = size(f);

%rescale image
fg = rescale_color_image(f);

%% CV models
%set parameters
pm.outer_iter = 40;
pm.alpha = 1.0;
pm.lambda = 1000;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.beta = 1.0;

%set image segmentation initialization
u1 = make_circle_shift_x(M,N,10, -5);
u2 = make_circle_shift_x(M,N,10, 5);
u1 = double(u1);
u2 = double(u2);

%L1-L2 segmentation
tic;
[L1L2_U1,L1L2_U2, L1L2_c1, L1L2_c2, L1L2_c3, L1L2_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

L1L2_c1_im = ones(N,M,3);
L1L2_c1_im(:,:,1) = ones(N,M)*L1L2_c1(1);
L1L2_c1_im(:,:,2) = ones(N,M)*L1L2_c1(2);
L1L2_c1_im(:,:,3) = ones(N,M)*L1L2_c1(3);

L1L2_c2_im = ones(N,M,3);
L1L2_c2_im(:,:,1) = ones(N,M)*L1L2_c2(1);
L1L2_c2_im(:,:,2) = ones(N,M)*L1L2_c2(2);
L1L2_c2_im(:,:,3) = ones(N,M)*L1L2_c2(3);

L1L2_c3_im = ones(N,M,3);
L1L2_c3_im(:,:,1) = ones(N,M)*L1L2_c3(1);
L1L2_c3_im(:,:,2) = ones(N,M)*L1L2_c3(2);
L1L2_c3_im(:,:,3) = ones(N,M)*L1L2_c3(3);

L1L2_c4_im = ones(N,M,3);
L1L2_c4_im(:,:,1) = ones(N,M)*L1L2_c4(1);
L1L2_c4_im(:,:,2) = ones(N,M)*L1L2_c4(2);
L1L2_c4_im(:,:,3) = ones(N,M)*L1L2_c4(3);

L1L2_approx_im = double(L1L2_U1>0.5).*double(L1L2_U2>0.5).*L1L2_c1_im + double(L1L2_U1>0.5).*double(L1L2_U2<=0.5).*L1L2_c2_im +...
    double(L1L2_U1<=0.5).*double(L1L2_U2>0.5).*L1L2_c3_im + double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5).*L1L2_c4_im;


%L1-0.75L2 segmentation
pm.alpha = 0.75;
tic;
[L1L2_75_U1,L1L2_75_U2, L1_0pt75_L2_c1, L1_0pt75_L2_c2, L1_0pt75_L2_c3, L1_0pt75_L2_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

L1_0pt75_L2_c1_im = ones(N,M,3);
L1_0pt75_L2_c1_im(:,:,1) = ones(N,M)*L1_0pt75_L2_c1(1);
L1_0pt75_L2_c1_im(:,:,2) = ones(N,M)*L1_0pt75_L2_c1(2);
L1_0pt75_L2_c1_im(:,:,3) = ones(N,M)*L1_0pt75_L2_c1(3);

L1_0pt75_L2_c2_im = ones(N,M,3);
L1_0pt75_L2_c2_im(:,:,1) = ones(N,M)*L1_0pt75_L2_c2(1);
L1_0pt75_L2_c2_im(:,:,2) = ones(N,M)*L1_0pt75_L2_c2(2);
L1_0pt75_L2_c2_im(:,:,3) = ones(N,M)*L1_0pt75_L2_c2(3);

L1_0pt75_L2_c3_im = ones(N,M,3);
L1_0pt75_L2_c3_im(:,:,1) = ones(N,M)*L1_0pt75_L2_c3(1);
L1_0pt75_L2_c3_im(:,:,2) = ones(N,M)*L1_0pt75_L2_c3(2);
L1_0pt75_L2_c3_im(:,:,3) = ones(N,M)*L1_0pt75_L2_c3(3);

L1_0pt75_L2_c4_im = ones(N,M,3);
L1_0pt75_L2_c4_im(:,:,1) = ones(N,M)*L1_0pt75_L2_c4(1);
L1_0pt75_L2_c4_im(:,:,2) = ones(N,M)*L1_0pt75_L2_c4(2);
L1_0pt75_L2_c4_im(:,:,3) = ones(N,M)*L1_0pt75_L2_c4(3);

L1_0pt75_L2_approx_im = double(L1L2_75_U1>0.5).*double(L1L2_75_U2>0.5).*L1_0pt75_L2_c1_im + double(L1L2_75_U1>0.5).*double(L1L2_75_U2<=0.5).*L1_0pt75_L2_c2_im +...
    double(L1L2_75_U1<=0.5).*double(L1L2_75_U2>0.5).*L1_0pt75_L2_c3_im + double(L1L2_75_U1<=0.5).*double(L1L2_75_U2<=0.5).*L1_0pt75_L2_c4_im;


%L1-0.5L2 segmentation
pm.alpha = 0.5;
tic;
[L1L2_05_U1,L1L2_05_U2, L1_0pt5_L2_c1, L1_0pt5_L2_c2, L1_0pt5_L2_c3, L1_0pt5_L2_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

L1_0pt5_L2_c1_im = ones(N,M,3);
L1_0pt5_L2_c1_im(:,:,1) = ones(N,M)*L1_0pt5_L2_c1(1);
L1_0pt5_L2_c1_im(:,:,2) = ones(N,M)*L1_0pt5_L2_c1(2);
L1_0pt5_L2_c1_im(:,:,3) = ones(N,M)*L1_0pt5_L2_c1(3);

L1_0pt5_L2_c2_im = ones(N,M,3);
L1_0pt5_L2_c2_im(:,:,1) = ones(N,M)*L1_0pt5_L2_c2(1);
L1_0pt5_L2_c2_im(:,:,2) = ones(N,M)*L1_0pt5_L2_c2(2);
L1_0pt5_L2_c2_im(:,:,3) = ones(N,M)*L1_0pt5_L2_c2(3);

L1_0pt5_L2_c3_im = ones(N,M,3);
L1_0pt5_L2_c3_im(:,:,1) = ones(N,M)*L1_0pt5_L2_c3(1);
L1_0pt5_L2_c3_im(:,:,2) = ones(N,M)*L1_0pt5_L2_c3(2);
L1_0pt5_L2_c3_im(:,:,3) = ones(N,M)*L1_0pt5_L2_c3(3);

L1_0pt5_L2_c4_im = ones(N,M,3);
L1_0pt5_L2_c4_im(:,:,1) = ones(N,M)*L1_0pt5_L2_c4(1);
L1_0pt5_L2_c4_im(:,:,2) = ones(N,M)*L1_0pt5_L2_c4(2);
L1_0pt5_L2_c4_im(:,:,3) = ones(N,M)*L1_0pt5_L2_c4(3);

L1_0pt5_L2_approx_im = double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5).*L1_0pt5_L2_c1_im + double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5).*L1_0pt5_L2_c2_im +...
    double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5).*L1_0pt5_L2_c3_im + double(L1L2_05_U1<=0.5).*double(L1L2_05_U2<=0.5).*L1_0pt5_L2_c4_im;


%L1-0.25L2 segmentation
pm.alpha = 0.25;
tic;
[L1L2_25_U1,L1L2_25_U2, L1_0pt25_L2_c1, L1_0pt25_L2_c2, L1_0pt25_L2_c3, L1_0pt25_L2_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

L1_0pt25_L2_c1_im = ones(N,M,3);
L1_0pt25_L2_c1_im(:,:,1) = ones(N,M)*L1_0pt25_L2_c1(1);
L1_0pt25_L2_c1_im(:,:,2) = ones(N,M)*L1_0pt25_L2_c1(2);
L1_0pt25_L2_c1_im(:,:,3) = ones(N,M)*L1_0pt25_L2_c1(3);

L1_0pt25_L2_c2_im = ones(N,M,3);
L1_0pt25_L2_c2_im(:,:,1) = ones(N,M)*L1_0pt25_L2_c2(1);
L1_0pt25_L2_c2_im(:,:,2) = ones(N,M)*L1_0pt25_L2_c2(2);
L1_0pt25_L2_c2_im(:,:,3) = ones(N,M)*L1_0pt25_L2_c2(3);

L1_0pt25_L2_c3_im = ones(N,M,3);
L1_0pt25_L2_c3_im(:,:,1) = ones(N,M)*L1_0pt25_L2_c3(1);
L1_0pt25_L2_c3_im(:,:,2) = ones(N,M)*L1_0pt25_L2_c3(2);
L1_0pt25_L2_c3_im(:,:,3) = ones(N,M)*L1_0pt25_L2_c3(3);

L1_0pt25_L2_c4_im = ones(N,M,3);
L1_0pt25_L2_c4_im(:,:,1) = ones(N,M)*L1_0pt25_L2_c4(1);
L1_0pt25_L2_c4_im(:,:,2) = ones(N,M)*L1_0pt25_L2_c4(2);
L1_0pt25_L2_c4_im(:,:,3) = ones(N,M)*L1_0pt25_L2_c4(3);

L1_0pt25_L2_approx_im = double(L1L2_25_U1>0.5).*double(L1L2_25_U2>0.5).*L1_0pt25_L2_c1_im + double(L1L2_25_U1>0.5).*double(L1L2_25_U2<=0.5).*L1_0pt25_L2_c2_im +...
    double(L1L2_25_U1<=0.5).*double(L1L2_25_U2>0.5).*L1_0pt25_L2_c3_im + double(L1L2_25_U1<=0.5).*double(L1L2_25_U2<=0.5).*L1_0pt25_L2_c4_im;


%anisotropic segmentation
pm.alpha = 0;
tic;
[ani_U1,ani_U2, ani_c1, ani_c2, ani_c3, ani_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

L1_c1_im = ones(N,M,3);
L1_c1_im(:,:,1) = ones(N,M)*ani_c1(1);
L1_c1_im(:,:,2) = ones(N,M)*ani_c1(2);
L1_c1_im(:,:,3) = ones(N,M)*ani_c1(3);

L1_c2_im = ones(N,M,3);
L1_c2_im(:,:,1) = ones(N,M)*ani_c2(1);
L1_c2_im(:,:,2) = ones(N,M)*ani_c2(2);
L1_c2_im(:,:,3) = ones(N,M)*ani_c2(3);

L1_c3_im = ones(N,M,3);
L1_c3_im(:,:,1) = ones(N,M)*ani_c3(1);
L1_c3_im(:,:,2) = ones(N,M)*ani_c3(2);
L1_c3_im(:,:,3) = ones(N,M)*ani_c3(3);

L1_c4_im = ones(N,M,3);
L1_c4_im(:,:,1) = ones(N,M)*ani_c4(1);
L1_c4_im(:,:,2) = ones(N,M)*ani_c4(2);
L1_c4_im(:,:,3) = ones(N,M)*ani_c4(3);

L1_approx_im = double(ani_U1>0.5).*double(ani_U2>0.5).*L1_c1_im + double(ani_U1>0.5).*double(ani_U2<=0.5).*L1_c2_im +...
    double(ani_U1<=0.5).*double(ani_U2>0.5).*L1_c3_im + double(ani_U1<=0.5).*double(ani_U2<=0.5).*L1_c4_im;

%% fuzzy method
%set parameters
pm2.outer_iter = 80;
pm2.alpha = 1.0;
pm2.lambda = 1000;
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.beta = 1.0;
pm2.nu = 650;

rng(1234);
u_initial{1} = rand(N,M);
u_initial{2} = rand(N,M);
u_initial{3} = rand(N,M);
u_initial{4} = rand(N,M);

sum_u = sum(cat(3,u_initial{:}),3);

u_initial{1} = u_initial{1}./sum_u;
u_initial{2} = u_initial{2}./sum_u;
u_initial{3} = u_initial{3}./sum_u;
u_initial{4} = u_initial{4}./sum_u;

%L1-L2
tic;
[fuzzy_L1_L2_f, fuzzy_L1_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

fuzzy_L1_L2_approx1 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{1}(1))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{2}(1))+...
    double((fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{3}(1))+...
    double((fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{4}(1));
fuzzy_L1_L2_approx2 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{1}(2))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{2}(2))+...
    double((fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{3}(2))+...
    double((fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{4}(2));
fuzzy_L1_L2_approx3 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{1}(3))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{2}(3))+...
    double((fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{1}).*(fuzzy_L1_L2_f{3}>fuzzy_L1_L2_f{4})*fuzzy_L1_L2_c{3}(3))+...
    double((fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{2}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{3}).*(fuzzy_L1_L2_f{4}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{4}(3));
fuzzy_L1_L2_approx = zeros(N,M,3);
fuzzy_L1_L2_approx(:,:,1) = fuzzy_L1_L2_approx1;
fuzzy_L1_L2_approx(:,:,2) = fuzzy_L1_L2_approx2;
fuzzy_L1_L2_approx(:,:,3) = fuzzy_L1_L2_approx3;

%L1-0.75L2
pm2.alpha = 0.75;
tic;
[fuzzy_L1_0pt75_L2_f, fuzzy_L1_0pt75_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

fuzzy_L1_0pt75_L2_approx1 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{1}(1))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{2}(1))+...
    double((fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{3}(1))+...
    double((fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{4}(1));
fuzzy_L1_0pt75_L2_approx2 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{1}(2))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{2}(2))+...
    double((fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{3}(2))+...
    double((fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{4}(2));
fuzzy_L1_0pt75_L2_approx3 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{1}(3))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{2}(3))+...
    double((fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{1}).*(fuzzy_L1_0pt75_L2_f{3}>fuzzy_L1_0pt75_L2_f{4})*fuzzy_L1_0pt75_L2_c{3}(3))+...
    double((fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{2}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{3}).*(fuzzy_L1_0pt75_L2_f{4}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{4}(3));

fuzzy_L1_0pt75_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt75_L2_approx(:,:,1) = fuzzy_L1_0pt75_L2_approx1;
fuzzy_L1_0pt75_L2_approx(:,:,2) = fuzzy_L1_0pt75_L2_approx2;
fuzzy_L1_0pt75_L2_approx(:,:,3) = fuzzy_L1_0pt75_L2_approx3;


%L1-0.5L2
pm2.alpha = 0.5;
tic;
[fuzzy_L1_0pt5_L2_f, fuzzy_L1_0pt5_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

fuzzy_L1_0pt5_L2_approx1 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{1}(1))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{2}(1))+...
    double((fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{3}(1))+...
    double((fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{4}(1));
fuzzy_L1_0pt5_L2_approx2 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{1}(2))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{2}(2))+...
    double((fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{3}(2))+...
    double((fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{4}(2));
fuzzy_L1_0pt5_L2_approx3 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{1}(3))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{2}(3))+...
    double((fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{1}).*(fuzzy_L1_0pt5_L2_f{3}>fuzzy_L1_0pt5_L2_f{4})*fuzzy_L1_0pt5_L2_c{3}(3))+...
    double((fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{2}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{3}).*(fuzzy_L1_0pt5_L2_f{4}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{4}(3));
fuzzy_L1_0pt5_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt5_L2_approx(:,:,1) = fuzzy_L1_0pt5_L2_approx1;
fuzzy_L1_0pt5_L2_approx(:,:,2) = fuzzy_L1_0pt5_L2_approx2;
fuzzy_L1_0pt5_L2_approx(:,:,3) = fuzzy_L1_0pt5_L2_approx3;


%L1-0.25L2
pm2.alpha = 0.25;
tic;
[fuzzy_L1_0pt25_L2_f, fuzzy_L1_0pt25_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

fuzzy_L1_0pt25_L2_approx1 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{1}(1))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{2}(1))+...
    double((fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{3}(1))+...
    double((fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{4}(1));
fuzzy_L1_0pt25_L2_approx2 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{1}(2))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{2}(2))+...
    double((fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{3}(2))+...
    double((fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{4}(2));
fuzzy_L1_0pt25_L2_approx3 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{1}(3))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{2}(3))+...
    double((fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{1}).*(fuzzy_L1_0pt25_L2_f{3}>fuzzy_L1_0pt25_L2_f{4})*fuzzy_L1_0pt25_L2_c{3}(3))+...
    double((fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{2}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{3}).*(fuzzy_L1_0pt25_L2_f{4}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{4}(3));

fuzzy_L1_0pt25_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt25_L2_approx(:,:,1) = fuzzy_L1_0pt25_L2_approx1;
fuzzy_L1_0pt25_L2_approx(:,:,2) = fuzzy_L1_0pt25_L2_approx2;
fuzzy_L1_0pt25_L2_approx(:,:,3) = fuzzy_L1_0pt25_L2_approx3;

%L1
pm2.alpha = 0.0;
tic;
[fuzzy_L1_f, fuzzy_L1_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

fuzzy_L1_approx1 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2}).*(fuzzy_L1_f{1}>fuzzy_L1_f{3}).*(fuzzy_L1_f{1}>fuzzy_L1_f{4})*fuzzy_L1_c{1}(1))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1}).*(fuzzy_L1_f{2}>fuzzy_L1_f{3}).*(fuzzy_L1_f{2}>fuzzy_L1_f{4})*fuzzy_L1_c{2}(1))+...
    double((fuzzy_L1_f{3}>fuzzy_L1_f{2}).*(fuzzy_L1_f{3}>fuzzy_L1_f{1}).*(fuzzy_L1_f{3}>fuzzy_L1_f{4})*fuzzy_L1_c{3}(1))+...
    double((fuzzy_L1_f{4}>fuzzy_L1_f{2}).*(fuzzy_L1_f{4}>fuzzy_L1_f{3}).*(fuzzy_L1_f{4}>fuzzy_L1_f{1})*fuzzy_L1_c{4}(1));
fuzzy_L1_approx2 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2}).*(fuzzy_L1_f{1}>fuzzy_L1_f{3}).*(fuzzy_L1_f{1}>fuzzy_L1_f{4})*fuzzy_L1_c{1}(2))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1}).*(fuzzy_L1_f{2}>fuzzy_L1_f{3}).*(fuzzy_L1_f{2}>fuzzy_L1_f{4})*fuzzy_L1_c{2}(2))+...
    double((fuzzy_L1_f{3}>fuzzy_L1_f{2}).*(fuzzy_L1_f{3}>fuzzy_L1_f{1}).*(fuzzy_L1_f{3}>fuzzy_L1_f{4})*fuzzy_L1_c{3}(2))+...
    double((fuzzy_L1_f{4}>fuzzy_L1_f{2}).*(fuzzy_L1_f{4}>fuzzy_L1_f{3}).*(fuzzy_L1_f{4}>fuzzy_L1_f{1})*fuzzy_L1_c{4}(2));
fuzzy_L1_approx3 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2}).*(fuzzy_L1_f{1}>fuzzy_L1_f{3}).*(fuzzy_L1_f{1}>fuzzy_L1_f{4})*fuzzy_L1_c{1}(3))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1}).*(fuzzy_L1_f{2}>fuzzy_L1_f{3}).*(fuzzy_L1_f{2}>fuzzy_L1_f{4})*fuzzy_L1_c{2}(3))+...
    double((fuzzy_L1_f{3}>fuzzy_L1_f{2}).*(fuzzy_L1_f{3}>fuzzy_L1_f{1}).*(fuzzy_L1_f{3}>fuzzy_L1_f{4})*fuzzy_L1_c{3}(3))+...
    double((fuzzy_L1_f{4}>fuzzy_L1_f{2}).*(fuzzy_L1_f{4}>fuzzy_L1_f{3}).*(fuzzy_L1_f{4}>fuzzy_L1_f{1})*fuzzy_L1_c{4}(3));
fuzzy_L1_approx = zeros(N,M,3);
fuzzy_L1_approx(:,:,1) = fuzzy_L1_approx1;
fuzzy_L1_approx(:,:,2) = fuzzy_L1_approx2;
fuzzy_L1_approx(:,:,3) = fuzzy_L1_approx3;

%% compute psnr
psnr_result = [psnr(L1L2_approx_im, fg),...
psnr(L1_0pt75_L2_approx_im, fg),...
psnr(L1_0pt5_L2_approx_im, fg),...
psnr(L1_0pt25_L2_approx_im, fg),...
psnr(L1_approx_im, fg),...
psnr(fuzzy_L1_L2_approx, fg),...
psnr(fuzzy_L1_0pt75_L2_approx, fg),...
psnr(fuzzy_L1_0pt5_L2_approx, fg),...
psnr(fuzzy_L1_0pt25_L2_approx, fg),...
psnr(fuzzy_L1_approx, fg)];

%% plot segmentation
figure;
subplot(2,6,1); imagesc(fg); axis off; axis square; title('Original');
subplot(2,6,2); imagesc(L1L2_approx_im); axis off; axis square; title('L1-L2 CV');
subplot(2,6,3); imagesc(L1_0pt75_L2_approx_im); axis off; axis square; title('L1-0.75L2 CV');
subplot(2,6,4); imagesc(L1_0pt5_L2_approx_im); axis off; axis square; title('L1-0.5L2 CV');
subplot(2,6,5); imagesc(L1_0pt25_L2_approx_im); axis off; axis square; title('L1-0.25L2 CV');
subplot(2,6,6); imagesc(L1_approx_im); axis off; axis square; title('L1 CV');
subplot(2,6,8); imagesc(fuzzy_L1_L2_approx); axis off; axis square; title('L1-L2 fuzzy');
subplot(2,6,9); imagesc(fuzzy_L1_0pt75_L2_approx); axis off; axis square; title('L1-0.75L2 fuzzy');
subplot(2,6,10); imagesc(fuzzy_L1_0pt5_L2_approx); axis off; axis square; title('L1-0.5L2 fuzzy');
subplot(2,6,11); imagesc(fuzzy_L1_0pt25_L2_approx); axis off; axis square; title('L1-0.25L2 fuzzy');
subplot(2,6,12); imagesc(fuzzy_L1_approx); axis off; axis square; title('L1 fuzzy');