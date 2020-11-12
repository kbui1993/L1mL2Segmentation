%%This scripts performs two-phase image segmentation using the AITV models
%%on a synthetic grayscale image with no impulse noise.

%% read image
fui8 = imread('shape.png');
fui8 = rgb2gray(fui8);
f = double(fui8);
[N,M] = size(f);
fg = rescale_image(f);


%% CV method
%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 2.0;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.beta = 1;
pm.nu = 10;

%preinitialize segmentation
u = make_circle(M,N,10);

u = double(u);

u_initial{1} = u;
u_initial{2} = ones(M,N)-u;

%L1-1.0L2
tic;
L1_L2_u1 = L1L2_two_phase(fg, u, pm);
toc

L1_L2_u1_f = fuzzy_L1L2(fg, u_initial, pm, 2);

%L1 - 0.75L2
pm.alpha = 0.75;
tic;
L1_0pt75_L2_u1 = L1L2_two_phase(fg, u, pm);
toc

L1_0pt75_L2_u1_f = fuzzy_L1L2(fg, u_initial, pm, 2);

%L1 - 0.5L2
pm.alpha = 0.5;
tic;
L1_0pt5_L2_u1 = L1L2_two_phase(fg, u, pm);
toc

%L1 - 0.25L2
pm.alpha = 0.5;
tic;
L1_0pt25_L2_u1 = L1L2_two_phase(fg, u, pm);
toc

%anisotropic CV
pm.alpha =0;
tic;
L1_u1 = L1L2_two_phase(fg, u, pm);
toc

%% fuzzy competition models
pm2.outer_iter = 40;
pm2.alpha = 1.0;

pm2.lambda = 2.0;
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
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

%% compute dice
dice(double(L1_L2_u1>0.5), fg)
dice(double(L1_0pt75_L2_u1>0.5), fg)
dice(double(L1_0pt5_L2_u1>0.5), fg)
dice(double(L1_0pt25_L2_u1>0.5), fg)
dice(double(L1_u1>0.5), fg)
dice(double(L1_L2_f{1}>L1_L2_f{2}), fg)
dice(double(L1_0pt75_L2_f{1}>L1_0pt75_L2_f{2}), fg)
dice(double(L1_0pt5_L2_f{1}>L1_0pt5_L2_f{2}), fg)
dice(double(L1_0pt25_L2_f{1}>L1_0pt25_L2_f{2}), fg)
dice(double(L1_f{1}>L1_f{2}), fg)

%% plot segmentation
figure;
subplot(2,6,1); imagesc(fg); colormap(gray); axis off; axis square; title('Original');
subplot(2,6,2); imagesc(L1_L2_u1>0.5); colormap(gray); axis off; axis square; title('L1-L2 CV');
subplot(2,6,3); imagesc(L1_0pt75_L2_u1>0.5); colormap(gray); axis off; axis square; title('L1-0.75L2 CV');
subplot(2,6,4); imagesc(L1_0pt5_L2_u1>0.5);  colormap(gray); axis off; axis square; title('L1-0.5L2 CV');
subplot(2,6,5); imagesc(L1_0pt25_L2_u1>0.5);  colormap(gray); axis off; axis square; title('L1-0.25L2 CV');
subplot(2,6,6); imagesc(L1_u1>0.5);colormap(gray); axis off; axis square; title('L1 CV');
subplot(2,6,8); imagesc(L1_L2_f{1} > L1_L2_f{2}); colormap(gray); axis off; axis square; title('L1-L2 fuzzy');
subplot(2,6,9); imagesc(L1_0pt75_L2_f{1} > L1_0pt75_L2_f{2}); colormap(gray); axis off; axis square; title('L1-0.75L2 fuzzy');
subplot(2,6,10); imagesc(L1_0pt5_L2_f{1} > L1_0pt5_L2_f{2}); colormap(gray); axis off; axis square; title('L1-0.5L2 fuzzy');
subplot(2,6,11); imagesc(L1_0pt25_L2_f{1} > L1_0pt25_L2_f{2}); colormap(gray); axis off; axis square; title('L1-0.25L2 fuzzy');
subplot(2,6,12); imagesc(L1_f{1} > L1_f{2}); colormap(gray); axis off; axis square; title('L1 fuzzy');
