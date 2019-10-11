%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script creates a color synthetic image for 2 phase image
%segmentation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image size 100 x 100
n=100;

%preinitialize image channel
M1=1+zeros(n,n);
M2=zeros(n,n);
M3=zeros(n,n);

%image location information
r1=15; c=[30,30]; r2=7;  c1=[75,50]; r3=10; 

%generate synthetic image
for i=1:n
    for j=1:n
        
        %create circle
        if sqrt((i-c(1))^2+(j-c(2))^2)<r1 && sqrt((i-c(1))^2+(j-c(2))^2)> r2     
            M1(i,j)=0;
            M2(i,j)=90;
            M3(i,j)=150;
        end
        
        %create weird shape
        if i>55 && i < 85 && j>25 && j <75
            if sqrt((i-c1(1))^2+(j-c1(2))^2) > r3 
                M1(i,j)=0;
                M2(i,j)=90;
                M3(i,j)=150;
            end
        end
    end
end

%create triangle
for i=20:51
    M1(i,i+35:86)=0;
    M2(i,i+35:86)=90;
    M3(i,i+35:86)=150;
end

%form the color image altogether from the channel image
M = zeros(n,n,3);
M(:,:,1) = M1;
M(:,:,2) = M2;
M(:,:,3) = M3;