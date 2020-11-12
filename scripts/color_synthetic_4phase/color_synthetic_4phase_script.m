%%This scripts performs four-phase image segmentation using the AITV models
%%on a synthetic color image with no impulse noise.

%% generate synthetic image
color_synthetic_image2;
f = double(M);
[N,M,~] = size(f);

%store rescaled image
fg = rescale_color_image(f);
fg = double(fg);

%%CV methods
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

%%CV method
%L1-L2
tic;
[L1L2_U1,L1L2_U2, L1L2_c1, L1L2_c2, L1L2_c3, L1L2_c4] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%L1-0.75L2
pm.alpha = 0.75;
tic;
[L1L2_75_U1,L1L2_75_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%L1-0.5L2
pm.alpha = 0.5;
tic;
[L1L2_50_U1,L1L2_50_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc;

%L1-0.25L2
pm.alpha = 0.25;
tic;
[L1L2_25_U1,L1L2_25_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%anisotropic
pm.alpha = 0;
pm.c = 0;
tic;
[L1_U1,L1_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

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
L1_L2_f = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

%L1-0.75L2
pm2.alpha = 0.75;
tic;
L1_0pt75_L2_f = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

%L1-0.5L2
pm2.alpha = 0.5;
tic;
L1_0pt5_L2_f = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

%L1-0.25L2
pm2.alpha = 0.25;
tic;
L1_0pt25_L2_f = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

%L1
pm2.alpha = 0;
tic;
L1_f = fuzzy_color_L1L2(fg, u_initial, pm2, 4);
toc

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
L1_0pt5_L2_M = ones(n,n)+3*(L1L2_50_U1.*L1L2_50_U2)+(L1L2_50_U1.*(1-L1L2_50_U2))+2*(1-L1L2_50_U1).*(L1L2_50_U2);
L1_0pt5_L2_M = uint8(L1_0pt5_L2_M);
L1_0pt5_L2_M = double(L1_0pt5_L2_M);

%L1-0.25L2
L1_0pt25_L2_M = ones(n,n)+3*(L1L2_25_U1.*L1L2_25_U2)+(L1L2_25_U1.*(1-L1L2_25_U2))+2*(1-L1L2_25_U1).*(L1L2_25_U2);
L1_0pt25_L2_M = uint8(L1_0pt25_L2_M);
L1_0pt25_L2_M = double(L1_0pt25_L2_M);

%L1
L1_M = ones(n,n)+3*(L1_U1.*L1_U2)+(L1_U1.*(1-L1_U2))+2*(1-L1_U1).*(L1_U2);
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


