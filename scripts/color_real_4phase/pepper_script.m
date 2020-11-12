%%This script performs seven-phase segmentation on a color image of a
%%pepper.

%% read pepper image
fui8 = imread('pepper.jpg');
%convert to double
f = double(fui8);
%get image size
[N,M,~] = size(f);

%rescale image
fg = rescale_color_image(f);


%% fuzzy method
%set parameters
pm2.outer_iter = 160;
pm2.alpha = 1.0;
pm2.lambda = 500;
pm2.c = 1e-8;
pm2.inner_iter = 300;
pm2.tau = 1/8;
pm2.beta = 1.0;
pm2.nu = 400;


rng(1234);

u_initial{1} = rand(N,M);
u_initial{2} = rand(N,M);
u_initial{3} = rand(N,M);
u_initial{4} = rand(N,M);
u_initial{5} = rand(N,M);
u_initial{6} = rand(N,M);
u_initial{7} = rand(N,M);


sum_u = sum(cat(3,u_initial{:}),3);

u_initial{1} = u_initial{1}./sum_u;
u_initial{2} = u_initial{2}./sum_u;
u_initial{3} = u_initial{3}./sum_u;
u_initial{4} = u_initial{4}./sum_u;
u_initial{5} = u_initial{5}./sum_u;
u_initial{6} = u_initial{6}./sum_u;
u_initial{7} = u_initial{7}./sum_u;




%L1-L2
tic;
[fuzzy_L1_L2_f, fuzzy_L1_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 7);
toc

for i=1:7
    fuzzy_L1_L2(:,:,i) = fuzzy_L1_L2_f{i};
end

[~,I] = max(fuzzy_L1_L2, [], 3);

fuzzy_L1_L2_approx1 = ...
    double(I==1)*fuzzy_L1_L2_c{1}(1)+...
    double(I==2)*fuzzy_L1_L2_c{2}(1)+...
    double(I==3)*fuzzy_L1_L2_c{3}(1)+...
    double(I==4)*fuzzy_L1_L2_c{4}(1)+...
    double(I==5)*fuzzy_L1_L2_c{5}(1)+...
    double(I==6)*fuzzy_L1_L2_c{6}(1)+...
    double(I==7)*fuzzy_L1_L2_c{7}(1);

fuzzy_L1_L2_approx2 = ...
    double(I==1)*fuzzy_L1_L2_c{1}(2)+...
    double(I==2)*fuzzy_L1_L2_c{2}(2)+...
    double(I==3)*fuzzy_L1_L2_c{3}(2)+...
    double(I==4)*fuzzy_L1_L2_c{4}(2)+...
    double(I==5)*fuzzy_L1_L2_c{5}(2)+...
    double(I==6)*fuzzy_L1_L2_c{6}(2)+...
    double(I==7)*fuzzy_L1_L2_c{7}(2);

fuzzy_L1_L2_approx3 = ...
    double(I==1)*fuzzy_L1_L2_c{1}(3)+...
    double(I==2)*fuzzy_L1_L2_c{2}(3)+...
    double(I==3)*fuzzy_L1_L2_c{3}(3)+...
    double(I==4)*fuzzy_L1_L2_c{4}(3)+...
    double(I==5)*fuzzy_L1_L2_c{5}(3)+...
    double(I==6)*fuzzy_L1_L2_c{6}(3)+...
    double(I==7)*fuzzy_L1_L2_c{7}(3);


fuzzy_L1_L2_approx = zeros(N,M,3);
fuzzy_L1_L2_approx(:,:,1) = fuzzy_L1_L2_approx1;
fuzzy_L1_L2_approx(:,:,2) = fuzzy_L1_L2_approx2;
fuzzy_L1_L2_approx(:,:,3) = fuzzy_L1_L2_approx3;


%L1-0.75L2
pm2.alpha = 0.75;
tic;
[fuzzy_L1_0pt75_L2_f, fuzzy_L1_0pt75_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 7);
toc

for i=1:7
    fuzzy_L1_0pt75_L2(:,:,i) = fuzzy_L1_0pt75_L2_f{i};
end

[~,I] = max(fuzzy_L1_0pt75_L2, [], 3);

fuzzy_L1_0pt75_L2_approx1 = ...
    double(I==1)*fuzzy_L1_0pt75_L2_c{1}(1)+...
    double(I==2)*fuzzy_L1_0pt75_L2_c{2}(1)+...
    double(I==3)*fuzzy_L1_0pt75_L2_c{3}(1)+...
    double(I==4)*fuzzy_L1_0pt75_L2_c{4}(1)+...
    double(I==5)*fuzzy_L1_0pt75_L2_c{5}(1)+...
    double(I==6)*fuzzy_L1_0pt75_L2_c{6}(1)+...
    double(I==7)*fuzzy_L1_0pt75_L2_c{7}(1);

fuzzy_L1_0pt75_L2_approx2 = ...
    double(I==1)*fuzzy_L1_0pt75_L2_c{1}(2)+...
    double(I==2)*fuzzy_L1_0pt75_L2_c{2}(2)+...
    double(I==3)*fuzzy_L1_0pt75_L2_c{3}(2)+...
    double(I==4)*fuzzy_L1_0pt75_L2_c{4}(2)+...
    double(I==5)*fuzzy_L1_0pt75_L2_c{5}(2)+...
    double(I==6)*fuzzy_L1_0pt75_L2_c{6}(2)+...
    double(I==7)*fuzzy_L1_0pt75_L2_c{7}(2);

fuzzy_L1_0pt75_L2_approx3 = ...
    double(I==1)*fuzzy_L1_0pt75_L2_c{1}(3)+...
    double(I==2)*fuzzy_L1_0pt75_L2_c{2}(3)+...
    double(I==3)*fuzzy_L1_0pt75_L2_c{3}(3)+...
    double(I==4)*fuzzy_L1_0pt75_L2_c{4}(3)+...
    double(I==5)*fuzzy_L1_0pt75_L2_c{5}(3)+...
    double(I==6)*fuzzy_L1_0pt75_L2_c{6}(3)+...
    double(I==7)*fuzzy_L1_0pt75_L2_c{7}(3);


fuzzy_L1_0pt75_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt75_L2_approx(:,:,1) = fuzzy_L1_0pt75_L2_approx1;
fuzzy_L1_0pt75_L2_approx(:,:,2) = fuzzy_L1_0pt75_L2_approx2;
fuzzy_L1_0pt75_L2_approx(:,:,3) = fuzzy_L1_0pt75_L2_approx3;


%L1-0.5L2
pm2.alpha = 0.5;
tic;
[fuzzy_L1_0pt5_L2_f, fuzzy_L1_0pt5_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 7);
toc


for i=1:7
    fuzzy_L1_0pt5_L2(:,:,i) = fuzzy_L1_0pt5_L2_f{i};
end

[~,I] = max(fuzzy_L1_0pt5_L2, [], 3);

fuzzy_L1_0pt5_L2_approx1 = ...
    double(I==1)*fuzzy_L1_0pt5_L2_c{1}(1)+...
    double(I==2)*fuzzy_L1_0pt5_L2_c{2}(1)+...
    double(I==3)*fuzzy_L1_0pt5_L2_c{3}(1)+...
    double(I==4)*fuzzy_L1_0pt5_L2_c{4}(1)+...
    double(I==5)*fuzzy_L1_0pt5_L2_c{5}(1)+...
    double(I==6)*fuzzy_L1_0pt5_L2_c{6}(1)+...
    double(I==7)*fuzzy_L1_0pt5_L2_c{7}(1);

fuzzy_L1_0pt5_L2_approx2 = ...
    double(I==1)*fuzzy_L1_0pt5_L2_c{1}(2)+...
    double(I==2)*fuzzy_L1_0pt5_L2_c{2}(2)+...
    double(I==3)*fuzzy_L1_0pt5_L2_c{3}(2)+...
    double(I==4)*fuzzy_L1_0pt5_L2_c{4}(2)+...
    double(I==5)*fuzzy_L1_0pt5_L2_c{5}(2)+...
    double(I==6)*fuzzy_L1_0pt5_L2_c{6}(2)+...
    double(I==7)*fuzzy_L1_0pt5_L2_c{7}(2);

fuzzy_L1_0pt5_L2_approx3 = ...
    double(I==1)*fuzzy_L1_0pt5_L2_c{1}(3)+...
    double(I==2)*fuzzy_L1_0pt5_L2_c{2}(3)+...
    double(I==3)*fuzzy_L1_0pt5_L2_c{3}(3)+...
    double(I==4)*fuzzy_L1_0pt5_L2_c{4}(3)+...
    double(I==5)*fuzzy_L1_0pt5_L2_c{5}(3)+...
    double(I==6)*fuzzy_L1_0pt5_L2_c{6}(3)+...
    double(I==7)*fuzzy_L1_0pt5_L2_c{7}(3);


fuzzy_L1_0pt5_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt5_L2_approx(:,:,1) = fuzzy_L1_0pt5_L2_approx1;
fuzzy_L1_0pt5_L2_approx(:,:,2) = fuzzy_L1_0pt5_L2_approx2;
fuzzy_L1_0pt5_L2_approx(:,:,3) = fuzzy_L1_0pt5_L2_approx3;



%L1-0.25L2
pm2.alpha = 0.25;
tic;
[fuzzy_L1_0pt25_L2_f, fuzzy_L1_0pt25_L2_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 7);
toc

for i=1:7
    fuzzy_L1_0pt25_L2(:,:,i) = fuzzy_L1_0pt25_L2_f{i};
end

[~,I] = max(fuzzy_L1_0pt25_L2, [], 3);

fuzzy_L1_0pt25_L2_approx1 = ...
    double(I==1)*fuzzy_L1_0pt25_L2_c{1}(1)+...
    double(I==2)*fuzzy_L1_0pt25_L2_c{2}(1)+...
    double(I==3)*fuzzy_L1_0pt25_L2_c{3}(1)+...
    double(I==4)*fuzzy_L1_0pt25_L2_c{4}(1)+...
    double(I==5)*fuzzy_L1_0pt25_L2_c{5}(1)+...
    double(I==6)*fuzzy_L1_0pt25_L2_c{6}(1)+...
    double(I==7)*fuzzy_L1_0pt25_L2_c{7}(1);

fuzzy_L1_0pt25_L2_approx2 = ...
    double(I==1)*fuzzy_L1_0pt25_L2_c{1}(2)+...
    double(I==2)*fuzzy_L1_0pt25_L2_c{2}(2)+...
    double(I==3)*fuzzy_L1_0pt25_L2_c{3}(2)+...
    double(I==4)*fuzzy_L1_0pt25_L2_c{4}(2)+...
    double(I==5)*fuzzy_L1_0pt25_L2_c{5}(2)+...
    double(I==6)*fuzzy_L1_0pt25_L2_c{6}(2)+...
    double(I==7)*fuzzy_L1_0pt25_L2_c{7}(2);

fuzzy_L1_0pt25_L2_approx3 = ...
    double(I==1)*fuzzy_L1_0pt25_L2_c{1}(3)+...
    double(I==2)*fuzzy_L1_0pt25_L2_c{2}(3)+...
    double(I==3)*fuzzy_L1_0pt25_L2_c{3}(3)+...
    double(I==4)*fuzzy_L1_0pt25_L2_c{4}(3)+...
    double(I==5)*fuzzy_L1_0pt25_L2_c{5}(3)+...
    double(I==6)*fuzzy_L1_0pt25_L2_c{6}(3)+...
    double(I==7)*fuzzy_L1_0pt25_L2_c{7}(3);


fuzzy_L1_0pt25_L2_approx = zeros(N,M,3);
fuzzy_L1_0pt25_L2_approx(:,:,1) = fuzzy_L1_0pt25_L2_approx1;
fuzzy_L1_0pt25_L2_approx(:,:,2) = fuzzy_L1_0pt25_L2_approx2;
fuzzy_L1_0pt25_L2_approx(:,:,3) = fuzzy_L1_0pt25_L2_approx3;


%L1
pm2.alpha = 0.0;
tic;
[fuzzy_L1_f, fuzzy_L1_c] = fuzzy_color_L1L2(fg, u_initial, pm2, 7);
toc

for i=1:7
    fuzzy_L1(:,:,i) = fuzzy_L1_f{i};
end

[~,I] = max(fuzzy_L1, [], 3);

fuzzy_L1_approx1 = ...
    double(I==1)*fuzzy_L1_c{1}(1)+...
    double(I==2)*fuzzy_L1_c{2}(1)+...
    double(I==3)*fuzzy_L1_c{3}(1)+...
    double(I==4)*fuzzy_L1_c{4}(1)+...
    double(I==5)*fuzzy_L1_c{5}(1)+...
    double(I==6)*fuzzy_L1_c{6}(1)+...
    double(I==7)*fuzzy_L1_c{7}(1);

fuzzy_L1_approx2 = ...
    double(I==1)*fuzzy_L1_c{1}(2)+...
    double(I==2)*fuzzy_L1_c{2}(2)+...
    double(I==3)*fuzzy_L1_c{3}(2)+...
    double(I==4)*fuzzy_L1_c{4}(2)+...
    double(I==5)*fuzzy_L1_c{5}(2)+...
    double(I==6)*fuzzy_L1_c{6}(2)+...
    double(I==7)*fuzzy_L1_c{7}(2);

fuzzy_L1_approx3 = ...
    double(I==1)*fuzzy_L1_c{1}(3)+...
    double(I==2)*fuzzy_L1_c{2}(3)+...
    double(I==3)*fuzzy_L1_c{3}(3)+...
    double(I==4)*fuzzy_L1_c{4}(3)+...
    double(I==5)*fuzzy_L1_c{5}(3)+...
    double(I==6)*fuzzy_L1_c{6}(3)+...
    double(I==7)*fuzzy_L1_c{7}(3);


fuzzy_L1_approx = zeros(N,M,3);
fuzzy_L1_approx(:,:,1) = fuzzy_L1_approx1;
fuzzy_L1_approx(:,:,2) = fuzzy_L1_approx2;
fuzzy_L1_approx(:,:,3) = fuzzy_L1_approx3;

%% compute psnr

psnr_result = [
psnr(fuzzy_L1_L2_approx, fg),...
psnr(fuzzy_L1_0pt75_L2_approx, fg),...
psnr(fuzzy_L1_0pt5_L2_approx, fg),...
psnr(fuzzy_L1_0pt25_L2_approx, fg),...
psnr(fuzzy_L1_approx, fg)];


%% plot figure
figure;
subplot(1,6,1); imagesc(fg); axis off; axis square; title('Original');
subplot(1,6,2); imagesc(fuzzy_L1_L2_approx); axis off; axis square; title('L1-L2 fuzzy');
subplot(1,6,3); imagesc(fuzzy_L1_0pt75_L2_approx); axis off; axis square; title('L1-0.75L2 fuzzy');
subplot(1,6,4); imagesc(fuzzy_L1_0pt5_L2_approx); axis off; axis square; title('L1-0.5L2 fuzzy');
subplot(1,6,5); imagesc(fuzzy_L1_0pt25_L2_approx); axis off; axis square; title('L1-0.25L2 fuzzy');
subplot(1,6,6); imagesc(fuzzy_L1_approx); axis off; axis square; title('L1 fuzzy');