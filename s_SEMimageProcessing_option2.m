%%%% This script is used to process SEM images (gray scale) - 2nd option
%%%% Useful when EDX has been acquired

%% 1. Load the image
clear;

folder = 'C:\Users\spet5905\Desktop\Valcun-B370LON\Backscatter-EDX\layer1\centre'; 
fileImg = 'centre-1-300.tif'; 
fileEDX = 'centre-1-300\Map_001_CountMap_Fe-K.csv';

img = imread([folder,'\',fileImg]); 
img = im2double(img); % Image with 3 channels. For gray-scale SEM image, the 3 channels are identical
map = readmatrix([folder,'\',fileEDX]);

if size(img,3)==4
   img(:,:,4) =[];
    
end    
%figure, imagesc(img); set(gca,'Visible','off'); axis image; % Display the original image
figure, imshow(img);

% Input the image pixel size 

prompt = {'Enter the image pixel size in \mum'};
dlgtitle  = 'Image pixel size';
dims = [1 40];
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,{'1'},opts);
pixsize = str2num(answer{1});  % in microns


%% 2.1 Process the image if SEM image (gray scale)

img = img(1:1920,:); % crop the infomation bar at the bottom of the image

inpaintFlag = 0;

inor = normalise(img(:,:,1),0.45,0.03);
%inor = imadjust(img(:,:,1));


if inpaintFlag
figure, imshow(inor); axis image;    
[ip bw] = improcess;

bwComb = bw(:,:,1);
  if size(bw,3)>1
    for n = 2:size(bw,3)
    
        bwComb = ~(~bwComb.*~bw(:,:,n));

    end
  end
else  
  ip = inor;
end

ib = bilateralFilterxray(ip,ip,12,0.1,12,0.1);
figure, imshow(ib); axis image;   


%% 2.2 Image segmentation of SEM image (gray scale), in this case using EDX map as a mask
map = imresize(map,xxxx);  % To be filled
% B = imresize(A,scale) returns image B that is scale times the size of image A. 
% The input image A can be a grayscale, RGB, binary, or categorical image.
% If A has more than two dimensions, then imresize only resizes the first two dimensions. 
% If scale is between 0 and 1, then B is smaller than A. If scale is greater than 1, then B is larger than A. By default, imresize uses bicubic interpolation.
map_nor = normalise(map); %normalise the map between [0,1]
% ?
bw = map_nor>0.2;
bw2 = bwareaopen(bw,400); % round(1/pixsize^2) Remove small objects smaller than 1 um^2
% Remove small areas of foreground (white regions) with bwareaopen. 
% All areas of foreground made up of fewer than n pixels are removed.
% Open is erosion followed by dilation using the same structuring element
% for both operations.
bw3 = imdilate(bw2,strel('square',15)); % Dilate the image 
% Take the maximum in a neighbourhood
L = imsegkmeans(im2single(bw3.*ib),3); % image clustering using K means
bw4 = L==3; % create a mask using the segmented phase; In this case, it is Fe-rich phase


s = regionprops(bw4,'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList');

IMCcent = cat(1,s.Centroid);
IMCarea = cat(1,s.Area);
IMClength = cat(1,s.MajorAxisLength);
IMCwidth = cat(1,s.MinorAxisLength);
IMCasprat = IMClength./IMCwidth;
IMCareafrac = sum(IMCarea)/numel(bw2);

figure, imshow(bw4);

%Nested function
function [ip bw] = improcess

ha=gca;
hnor = findobj(ha,'Type','Image');

inor = hnor.CData;

bw = roipoly;


ip = inpaintExemplar(inor,bw(:,:,end));

imshow(ip); ha=gca; % update the handle

hf  = figure('color','w','position',[1200 600 250 70]);


multiroipoly(ha);

uiwait(hf); % Stops the code until the user draws all the lines

function multiroipoly(ha)

 %hf  = figure('color','w','position',[1200 600 250 70]);
    hf.MenuBar = 'none';
    hf.ToolBar = 'none';
    hf.Name = 'Draw poly';
 addpoly = uicontrol(hf,'Style','pushbutton',...
                    'FontSize',14,...
                    'fontweight','bold',...
                    'String','Add poly',...
                    'Callback',@addpoly_Callback);
                
 hOK = uicontrol(hf,'Style','pushbutton',...
                    'FontSize',14,...
                    'fontweight','bold',...
                    'String','OK',...
                    'Callback',@okbutton_Callback);
                
                addpoly.Units = 'pixels';
                addpoly.Position = [30 10 100 50 ];
                hOK.Units = 'pixels';
                hOK.Position = [150 10 70 50];
                
                % UI Callbacks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                function addpoly_Callback(source,callbackdata)
                    axes(ha);
                    bw(:,:,end+1) = roipoly;
                    
                    ip = inpaintExemplar(ip,bw(:,:,end));
                    imshow(ip); ha = gca; % update the handle
                end % addpoly_Callback
                
                function okbutton_Callback(source,callbackdata)
                    % Display surf plot of the currently selected data.
                    
                    close(hf);
                end % okbutton_Callback

end


end
