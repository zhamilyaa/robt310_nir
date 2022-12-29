% ROBT310: Project – Final Report
% Name: Nurdaulet Zhuzbay
% Name: Zhamilya Saparova
% Name: Mukhamedzhan Nurmukhamed
%%

%
Z = (imread('IR3.png'));
Z=im2double(Z);
IR1 = histeq(Z);
% [J, rect] = imcrop(Z);
%%

Colorim1=im2double(imread('RGB3_m.jpg'));

% [height, width, dim] = size(Colorim1);
% IR1 = imresize(J,[height width]);

%%
figure(1)
%imshow(B)
imshow(medianfilter(sharpening(histeq(Colorim1))));
title('[Input]');
% 
figure(2)
imshow(histeq(colorization(medianfilter(segmentation(sharpening(IR1), Colorim1)), IR1)));
title('[Output]');



function M = medianfilter(I)

[height, width, dim] = size(I);
%assiging an output median filtered image
M = double(zeros(height, width, dim));
     
%applying median filter
for k = 1:dim
   for i = 1:height
       for j = 1:width
                     
            if i~= 1 && j~=1 && i~=height && j~=width 
                     
            X = I(i,j,k);
                  
            A = I(i-1,j-1,k);
            B = I(i-1,j,k);
            C = I(i-1,j+1,k);
            D = I(i,j-1,k);
            E = I(i,j+1,k);
            F = I(i+1,j-1,k);
            G = I(i+1,j,k);
            H = I(i+1,j+1,k); 
                     
            array = [ X A B C D E F G H];
            sorted_array = sort(array);
                     
                     
            median = sorted_array(5);
                     
            M(i,j,k) = median;
            else 
                M(i,j,k) = I(i,j,k);
                  
            end
       end
   end
end

end


function S = sharpening(I)

[height, width, dim] = size(I);

%blurring the image
I_blurred = imgaussfilt(I,1);

%calculating g_mask
G_mask = double(zeros(height, width, dim));
S = double(zeros(height, width, dim));

coef = 1;

for k = 1:dim
   for i = 1:height
       for j = 1:width
           
           G_mask(i,j,k) = I(i,j,k) - I_blurred(i,j,k);
           S(i,j,k) = I(i,j,k) + coef*G_mask(i,j,k);
           
           
       end
   end
end
end



function ColorRegion = segmentation(I, Colorim)

[height, width, dim] = size(I);

Red=0;
Green=0;
Blue=0;
simPixNum=0;
counter = 1;

I_segmented = uint32(zeros(height, width, dim));
ColorRegion = double(zeros(height, width, 3));


for i = 1: height
    for j = 1:width
        
        x = i;
        y = j;

        if I_segmented(x,y) == 0
            
            segment = regiongrowing(I,x,y,0.006); 
            
            for a = 1:height
                for b = 1:width
                    
                    if segment(a,b) == 1
                        
                        if I_segmented(a,b) == 0
                            I_segmented(a,b) = counter;
                            simPixNum=simPixNum+1;
                            Red=Red+Colorim(a,b,1);
                            Green=Green+Colorim(a,b,2);
                            Blue=Blue+Colorim(a,b,3);
                            
                        end
                    end
                    
                end
            end

            for a = 1:height
                for b = 1:width
                    
                    
                    if I_segmented(a,b)==counter
                        ColorRegion(a,b,1)=Red/simPixNum;
                        ColorRegion(a,b,2)=Green/simPixNum;
                        ColorRegion(a,b,3)=Blue/simPixNum;

                    end
                end
            end
            
            counter = counter + 1;
            Red=0;
            Green=0;
            Blue=0;
            simPixNum=0;

        end
    end
end
end


function C = colorization(I, Z)

[height, width,dim] = size(I);

C = double(zeros(height,width,dim));

irMax=max(max(Z));

RedMax=max(max(I(:,:,1)));
GreenMax=max(max(I(:,:,2)));
BlueMax=max(max(I(:,:,3)));

for k = 1:3
 for i = 1:height
     for j = 1:width
         
         x = Z(i,j);
                   
         if (k == 1)
             C(i,j,k) = x/irMax * I(i,j,k)/RedMax; 
         end
         
         if (k == 2)
             C(i,j,k) = x/irMax * I(i,j,k)/GreenMax;%/0.9804;
         end
         
         if (k == 3)
             C(i,j,k) = x/irMax * I(i,j,k)/BlueMax;%/0.9804;
         end 
      end
 end
end

end

function J=regiongrowing(I,x,y,maxdist)

[height,width,dim] = size(I);

J = zeros(height,width); 
M = I(x,y);  %The mean value of the segmented region
NP = 1;      %number of pixels
nf = 10000; 
np=0;
nl = zeros(nf,3); 
Pd = 0;  %pixel distance
N=[-1 0;     %neighbour pixels
    1 0; 
   0 -1;
    0 1];
while( (Pd < maxdist) && (NP < numel(I)))
    for j=1:4
        xn = x + N(j,1);
        yn = y + N(j,2); %coordinates of neighbour
        
        
    if( (xn >= 1) && (yn >= 1) && (xn <= height) && (yn <= width) && (J(xn,yn)==0) ) 
                np = np+1;
                nl(np,:) = [xn yn I(xn,yn)]; 
                J(xn,yn)=1;
    end
    end
    
    if( (np+10) > nf)
        nf=nf+10000; 
        nl((np+1):nf,:)=0;
    end
    
    dist = abs(nl(1:np,3)-M);
    [Pd, index] = min(dist);
    J(x,y) = 2; 
    NP = NP + 1;
    M = (M*NP + nl(index,3))/(NP+1); %new mean value
    
    x = nl(index,1); 
    y = nl(index,2);
    
    nl(index,:) = nl(np,:); 
    np = np-1;
end
J=J>1;
end


