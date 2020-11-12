%%This script performs two-phase segmentation on a color image of a
%%hawk.

%read hawk image
fui8 = imread('hawk.jpg');

%convert to double
f = double(fui8);

%get size of image
[N,M,~] = size(f);

%rescale image
fg = rescale_color_image(f);

%%Chan-Vese models
%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 10.0;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.beta = 1.0;


%set segmentation initialization
u = make_circle(M,N,10);
u = double(u);

%L1-L2
tic;
[L1_L2_u1, L1_L2_c1, L1_L2_c2] = L1L2_color_two_phase(fg, u, pm);
time = toc

L1_L2_approx1 = (L1_L2_u1>0.5)*L1_L2_c1(1)+(L1_L2_u1<0.5)*L1_L2_c2(1);
L1_L2_approx2 = (L1_L2_u1>0.5)*L1_L2_c1(2)+(L1_L2_u1<0.5)*L1_L2_c2(2);
L1_L2_approx3 = (L1_L2_u1>0.5)*L1_L2_c1(3)+(L1_L2_u1<0.5)*L1_L2_c2(3);

L1_L2_approx = zeros(N,M,3);
L1_L2_approx(:,:,1) = L1_L2_approx1;
L1_L2_approx(:,:,2) = L1_L2_approx2;
L1_L2_approx(:,:,3) = L1_L2_approx3;

%L1-0.75L2
pm.alpha =0.75;
tic;
[L1_0pt75_L2_u1, L1_0pt75_L2_c1, L1_0pt75_L2_c2] = L1L2_color_two_phase(fg, u, pm);
time = toc

L1_0pt75_L2_approx1 = (L1_0pt75_L2_u1>0.5)*L1_0pt75_L2_c1(1)+(L1_0pt75_L2_u1<0.5)*L1_0pt75_L2_c2(1);
L1_0pt75_L2_approx2 = (L1_0pt75_L2_u1>0.5)*L1_0pt75_L2_c1(2)+(L1_0pt75_L2_u1<0.5)*L1_0pt75_L2_c2(2);
L1_0pt75_L2_approx3 = (L1_0pt75_L2_u1>0.5)*L1_0pt75_L2_c1(3)+(L1_0pt75_L2_u1<0.5)*L1_0pt75_L2_c2(3);

L1_0pt75_L2_approx = zeros(N,M,3);
L1_0pt75_L2_approx(:,:,1) = L1_0pt75_L2_approx1;
L1_0pt75_L2_approx(:,:,2) = L1_0pt75_L2_approx2;
L1_0pt75_L2_approx(:,:,3) = L1_0pt75_L2_approx3;


%L1-0.5L2
pm.alpha =0.5;
tic;
[L1_0pt5_L2_u1, L1_0pt5_L2_c1, L1_0pt5_L2_c2] = L1L2_color_two_phase(fg, u, pm);
time = toc

L1_0pt5_L2_approx1 = (L1_0pt5_L2_u1>0.5)*L1_0pt5_L2_c1(1)+(L1_0pt5_L2_u1<0.5)*L1_0pt5_L2_c2(1);
L1_0pt5_L2_approx2 = (L1_0pt5_L2_u1>0.5)*L1_0pt5_L2_c1(2)+(L1_0pt5_L2_u1<0.5)*L1_0pt5_L2_c2(2);
L1_0pt5_L2_approx3 = (L1_0pt5_L2_u1>0.5)*L1_0pt5_L2_c1(3)+(L1_0pt5_L2_u1<0.5)*L1_0pt5_L2_c2(3);

L1_0pt5_L2_approx = zeros(N,M,3);
L1_0pt5_L2_approx(:,:,1) = L1_0pt5_L2_approx1;
L1_0pt5_L2_approx(:,:,2) = L1_0pt5_L2_approx2;
L1_0pt5_L2_approx(:,:,3) = L1_0pt5_L2_approx3;


%L1-0.25L2
pm.alpha =0.25;
tic;
[L1_0pt25_L2_u1, L1_0pt25_L2_c1, L1_0pt25_L2_c2] = L1L2_color_two_phase(fg, u, pm);
time = toc


L1_0pt25_L2_approx1 = (L1_0pt25_L2_u1>0.5)*L1_0pt25_L2_c1(1)+(L1_0pt25_L2_u1<0.5)*L1_0pt25_L2_c2(1);
L1_0pt25_L2_approx2 = (L1_0pt25_L2_u1>0.5)*L1_0pt25_L2_c1(2)+(L1_0pt25_L2_u1<0.5)*L1_0pt25_L2_c2(2);
L1_0pt25_L2_approx3 = (L1_0pt25_L2_u1>0.5)*L1_0pt25_L2_c1(3)+(L1_0pt25_L2_u1<0.5)*L1_0pt25_L2_c2(3);

L1_0pt25_L2_approx = zeros(N,M,3);
L1_0pt25_L2_approx(:,:,1) = L1_0pt25_L2_approx1;
L1_0pt25_L2_approx(:,:,2) = L1_0pt25_L2_approx2;
L1_0pt25_L2_approx(:,:,3) = L1_0pt25_L2_approx3;

%anisotropic CV
pm.alpha = 0;
tic;
[L1_u1, L1_c1, L1_c2] = L1L2_color_two_phase(fg, u, pm);
time = toc

L1_approx1 = (L1_u1>0.5)*L1_c1(1)+(L1_u1<0.5)*L1_c2(1);
L1_approx2 = (L1_u1>0.5)*L1_c1(2)+(L1_u1<0.5)*L1_c2(2);
L1_approx3 = (L1_u1>0.5)*L1_c1(3)+(L1_u1<0.5)*L1_c2(3);

L1_approx = zeros(N,M,3);
L1_approx(:,:,1) = L1_approx1;
L1_approx(:,:,2) = L1_approx2;
L1_approx(:,:,3) = L1_approx3;

%%fuzzy method
%set parameters
pm2.outer_iter = 40;
pm2.alpha = 1.0;
pm2.lambda = 10.0;
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.beta = 1.0;
pm2.nu = 10.0;

u_initial{1} = u;
u_initial{2} = ones(N,M)-u;

%L1-L2
tic;
[fuzzy_L1_L2_f, fuzzy_L1_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 2);
toc

fuzzy_L1_L2_approx1 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2})*fuzzy_L1_L2_c{1}(1))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{2}(1));
fuzzy_L1_L2_approx2 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2})*fuzzy_L1_L2_c{1}(2))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{2}(2));
fuzzy_L1_L2_approx3 =...
    double((fuzzy_L1_L2_f{1}>fuzzy_L1_L2_f{2})*fuzzy_L1_L2_c{1}(3))+...
    double((fuzzy_L1_L2_f{2}>fuzzy_L1_L2_f{1})*fuzzy_L1_L2_c{2}(3));
fuzzy_L1_L2_approx = zeros(N,M,3);
fuzzy_L1_L2_approx(:,:,1) = fuzzy_L1_L2_approx1;
fuzzy_L1_L2_approx(:,:,2) = fuzzy_L1_L2_approx2;
fuzzy_L1_L2_approx(:,:,3) = fuzzy_L1_L2_approx3;

%L1-0.75L2
pm2.alpha = 0.75;
tic;
[fuzzy_L1_0pt75_L2_f, fuzzy_L1_0pt75_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 2);
toc

fuzzy_L1_0pt75_L2_approx1 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2})*fuzzy_L1_0pt75_L2_c{1}(1))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{2}(1));
fuzzy_L1_0pt75_L2_approx2 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2})*fuzzy_L1_0pt75_L2_c{1}(2))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{2}(2));
fuzzy_L1_0pt75_L2_approx3 =...
    double((fuzzy_L1_0pt75_L2_f{1}>fuzzy_L1_0pt75_L2_f{2})*fuzzy_L1_0pt75_L2_c{1}(3))+...
    double((fuzzy_L1_0pt75_L2_f{2}>fuzzy_L1_0pt75_L2_f{1})*fuzzy_L1_0pt75_L2_c{2}(3));
fuzzy_L1_0pt75_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt75_L2_approx(:,:,1) = fuzzy_L1_0pt75_L2_approx1;
fuzzy_L1_0pt75_L2_approx(:,:,2) = fuzzy_L1_0pt75_L2_approx2;
fuzzy_L1_0pt75_L2_approx(:,:,3) = fuzzy_L1_0pt75_L2_approx3;

%L1-0.5L2
pm2.alpha = 0.5;
tic;
[fuzzy_L1_0pt5_L2_f, fuzzy_L1_0pt5_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 2);
toc

fuzzy_L1_0pt5_L2_approx1 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2})*fuzzy_L1_0pt5_L2_c{1}(1))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{2}(1));
fuzzy_L1_0pt5_L2_approx2 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2})*fuzzy_L1_0pt5_L2_c{1}(2))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{2}(2));
fuzzy_L1_0pt5_L2_approx3 =...
    double((fuzzy_L1_0pt5_L2_f{1}>fuzzy_L1_0pt5_L2_f{2})*fuzzy_L1_0pt5_L2_c{1}(3))+...
    double((fuzzy_L1_0pt5_L2_f{2}>fuzzy_L1_0pt5_L2_f{1})*fuzzy_L1_0pt5_L2_c{2}(3));
fuzzy_L1_0pt5_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt5_L2_approx(:,:,1) = fuzzy_L1_0pt5_L2_approx1;
fuzzy_L1_0pt5_L2_approx(:,:,2) = fuzzy_L1_0pt5_L2_approx2;
fuzzy_L1_0pt5_L2_approx(:,:,3) = fuzzy_L1_0pt5_L2_approx3;

%L1-0.25L2
pm2.alpha = 0.25;
tic;
[fuzzy_L1_0pt25_L2_f, fuzzy_L1_0pt25_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 2);
toc

fuzzy_L1_0pt25_L2_approx1 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2})*fuzzy_L1_0pt25_L2_c{1}(1))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{2}(1));
fuzzy_L1_0pt25_L2_approx2 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2})*fuzzy_L1_0pt25_L2_c{1}(2))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{2}(2));
fuzzy_L1_0pt25_L2_approx3 =...
    double((fuzzy_L1_0pt25_L2_f{1}>fuzzy_L1_0pt25_L2_f{2})*fuzzy_L1_0pt25_L2_c{1}(3))+...
    double((fuzzy_L1_0pt25_L2_f{2}>fuzzy_L1_0pt25_L2_f{1})*fuzzy_L1_0pt25_L2_c{2}(3));
fuzzy_L1_0pt25_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt25_L2_approx(:,:,1) = fuzzy_L1_0pt25_L2_approx1;
fuzzy_L1_0pt25_L2_approx(:,:,2) = fuzzy_L1_0pt25_L2_approx2;
fuzzy_L1_0pt25_L2_approx(:,:,3) = fuzzy_L1_0pt25_L2_approx3;

%L1
pm2.alpha = 0.0;
tic;
[fuzzy_L1_f, fuzzy_L1_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 2);
toc

fuzzy_L1_approx1 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2})*fuzzy_L1_c{1}(1))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1})*fuzzy_L1_c{2}(1));
fuzzy_L1_approx2 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2})*fuzzy_L1_c{1}(2))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1})*fuzzy_L1_c{2}(2));
fuzzy_L1_approx3 =...
    double((fuzzy_L1_f{1}>fuzzy_L1_f{2})*fuzzy_L1_c{1}(3))+...
    double((fuzzy_L1_f{2}>fuzzy_L1_f{1})*fuzzy_L1_c{2}(3));
fuzzy_L1_approx = zeros(N,M,3);
fuzzy_L1_approx(:,:,1) = fuzzy_L1_approx1;
fuzzy_L1_approx(:,:,2) = fuzzy_L1_approx2;
fuzzy_L1_approx(:,:,3) = fuzzy_L1_approx3;

psnr_result = [psnr(L1_L2_approx, fg),...
psnr(L1_0pt75_L2_approx, fg),...
psnr(L1_0pt5_L2_approx, fg),...
psnr(L1_0pt25_L2_approx, fg),...
psnr(L1_approx, fg),...
psnr(fuzzy_L1_L2_approx, fg),...
psnr(fuzzy_L1_0pt75_L2_approx, fg),...
psnr(fuzzy_L1_0pt5_L2_approx, fg),...
psnr(fuzzy_L1_0pt25_L2_approx, fg),...
psnr(fuzzy_L1_approx, fg)];

%plot segmentation
figure;
subplot(2,6,1); imagesc(fg); axis off; axis square; title('Original');
subplot(2,6,2); imagesc(L1_L2_approx); axis off; axis square; title('L1-L2 CV');
subplot(2,6,3); imagesc(L1_0pt75_L2_approx); axis off; axis square; title('L1-0.75L2 CV');
subplot(2,6,4); imagesc(L1_0pt5_L2_approx); axis off; axis square; title('L1-0.5L2 CV');
subplot(2,6,5); imagesc(L1_0pt25_L2_approx); axis off; axis square; title('L1-0.25L2 CV');
subplot(2,6,6); imagesc(L1_approx); axis off; axis square; title('L1 CV');
subplot(2,6,8); imagesc(fuzzy_L1_L2_approx); axis off; axis square; title('L1-L2 fuzzy');
subplot(2,6,9); imagesc(fuzzy_L1_0pt75_L2_approx); axis off; axis square; title('L1-0.75L2 fuzzy');
subplot(2,6,10); imagesc(fuzzy_L1_0pt5_L2_approx); axis off; axis square; title('L1-0.5L2 fuzzy');
subplot(2,6,11); imagesc(fuzzy_L1_0pt25_L2_approx); axis off; axis square; title('L1-0.25L2 fuzzy');
subplot(2,6,12); imagesc(fuzzy_L1_approx); axis off; axis square; title('L1 fuzzy');