%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generates a synthetic image for grayscale 4 phase image
%segmentation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image size 100 x100
n=100;

%preinitialize image
M=1+zeros(n,n);     

%image location information
r1=15; c=[30,30]; r2=7;  c1=[75,50]; r3=10; 

%create synthetic image
for i=1:n
    for j=1:n
        
        %create circle
        if sqrt((i-c(1))^2+(j-c(2))^2)<r1 && sqrt((i-c(1))^2+(j-c(2))^2)> r2     
            M(i,j)=0.6;
        end
        
        %create weird shape
        if i>55 && i < 85 && j>25 && j <75
            if sqrt((i-c1(1))^2+(j-c1(2))^2) > r3 
                M(i,j)=0.9;     
            end
        end
    end
end

%create triangle
for i=20:51
    M(i,i+35:86)=0.3;                    
end

%rescale image
M=rescale(M);