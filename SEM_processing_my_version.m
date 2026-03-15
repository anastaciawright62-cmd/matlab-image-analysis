%%%% This script is used to process SEM images (gray scale)
%%%% 

%% 1. Load the image
clear;

folder = 'C:\Users\Anastacia\OneDrive - Nexus365\part II\part ii micrographs\SEM\evo\2nd feb'; 
file = ['Pure 6063 sample_01.tif']; 
  
img = imread([folder,'\',file]); 
img = im2double(img); % Image with 3 channels. For gray-scale SEM image, the 3 channels are identical

if size(img,3)==4 
    % szdim = size(A,dim) returns the length of dimension dim when dim is a positive integer scalar. 
    % You can also specify dim as a vector of positive integers to query multiple dimension lengths at a time.
   img(:,:,4) =[];
    
end    
%figure, imagesc(img); set(gca,'Visible','off'); axis image; % Display the original image
figure, imshow(img);

% Input the image pixel size 

prompt = {'Enter the image pixel size in \mum'};
dlgtitle  = 'Image pixel size';
dims = [1 40];
opts.Interpreter = 'tex';
% specifies that the dialog box is resizeable in the horizontal direction when opts is set to 'on'. 
% When opts is a structure, it specifies whether the dialog box is resizeable in the horizontal
% direction, whether it is modal, and whether the prompt text is interpreted. Here 'tex' as the interpreter?
answer = inputdlg(prompt,dlgtitle,dims,{'1'},opts);
% creates a modal dialog box containing one or more text edit fields and returns the values entered by the user. 
% The return values are elements of a cell array of character vectors. The first element of the cell array corresponds 
% to the response in the edit field at the top of the dialog box. The second element corresponds to the next edit field response, and so on.
pixsize = str2num(answer{1});  % in microns
% converts a character array or string scalar (first inupt of answer)) to a numeric matrix. The input can include spaces, commas, 
% and semicolons to indicate separate elements. If str2num cannot parse the input as numeric values, then it returns an empty matrix.
% Users can enter scalar or vector values into inputdlg text edit fields. 
% MATLAB® stores the input as a cell array of character vectors. 
% Convert a member of the input cell array to a number, using str2num.
%% 2.1 Process the image if SEM image (gray scale)

img = img(1:700,:); % crop the infomation bar at the bottom of the image - was 1920 before I edited

%figure, imshow(img);

inpaintFlag = 0; %just a flag 

% inor = mat2gray(img(:,:,1));%0.45,0.03);
inor = normalise(img(:,:,1),0.45,0.03); %'normalise' is a separate function to the built-in 'normalize' function in Matlab so if you get an error message, request the normalise function from Enzo/the previous owner
% inor = imadjust(img(:,:,1));
% Offsets and rescales image so that the minimum value is 0 and the maximum value is 1.  Result is returned in n.  
% If the image is colour the image is converted to HSV and the value/intensity component is normalised to 0-1 before being converted back to RGB.
% n = normalise(im, reqmean, reqvar)
% Arguments:  im      - A grey-level input image.
            % reqmean - The required mean value of the image.
            % reqvar  - The required variance of the image.
% Offsets and rescales image so that it has mean reqmean and variance reqvar.  Colour images cannot be normalised in this manner.

if inpaintFlag
figure, imshow(inor); axis image;    
[ip bw] = improcess;

bwComb = bw(:,:,1);
  if size(bw,3)>1
    for n = 2:size(bw,3)
    % szdim = size(A,dim) returns the length of dimension dim when dim is a positive integer scalar. 
    % You can also specify dim as a vector of positive integers to query multiple dimension lengths at a time. 
    % For example, size(A,[2 3]) returns the lengths of the second and third dimensions of A in the 1-by-2 row vector szdim.
        bwComb = ~(~bwComb.*~bw(:,:,n));

    end
  end
else  
  ip = inor;
end

ib = bilateralFilterxray(ip,ip,12,0.1,12,0.1);%also external script, this is basically a gaussian filter but they blurr edges so bilaterial separates the edges and rest and doesn't smooth for the edge
% output = bilateralFilterxray( data, edge, sigmaSpatial, sigmaRange, ... samplingSpatial, samplingRange )
% Bilaterally filters the image 'data' using the edges in the image 'edge'.
% If 'data' == 'edge', then it's the normal bilateral filter. Else, then it is the "cross" or "joint" bilateral filter.
% data and edge should be of the same size and greyscale. (i.e. they should be ( height x width x 1 matrices ))
% data is the only required argument
% By default:
% edge = data
% sigmaSpatial = samplingSpatial = min( width, height ) / 16;
% sigmaRange = samplingRange = ( max( edge( : ) ) - min( edge( : ) ) ) / 10
figure, imshow(ib); axis image;   


%% 2.2 Image segmentation of SEM image (gray scale)

L = imsegkmeans(im2single(ib),4); % image clustering using K means
% L = imsegkmeans(I,k) segments image I into k clusters by performing k-means clustering 
% and returns the segmented labeled output in L.
% J = im2single(I) converts the grayscale, RGB, or binary image I to single, rescaling or offsetting the data as necessary.
bw = L==2; % create a mask using the segmented phase; In this case, it is Fe-rich phase
bw2 = bwareaopen(bw,36); % round(1/pixsize^2) Remove small objects smaller than 1 um^2
% BW2 = bwareaopen(BW,P) removes all connected components (objects) that have fewer than P pixels from the binary image BW, producing another binary image, BW2. 
% This operation is known as an area opening.
s = regionprops(bw2,'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList');
% stats = regionprops(BW,properties) measures properties for each object in the binary image BW.
IMCcent = cat(1,s.Centroid);
% C = cat(dim,A,B) concatenates B to the end of A along dimension dim when A and B have compatible sizes 
% (the lengths of the dimensions match except for the operating dimension dim).
IMCarea = cat(1,s.Area);
IMClength = cat(1,s.MajorAxisLength);
IMCwidth = cat(1,s.MinorAxisLength);
IMCasprat = IMClength./IMCwidth;
IMCareafrac = sum(IMCarea)/numel(bw2);

figure, imshow(bw2);

%Nested function
function [ip bw] = improcess
    ha = gca;
    hnor = findobj(ha,'Type','Image');
    inor = hnor.CData;

    % Ensure image gets mouse events
    figImg = ancestor(ha,'figure');
    figure(figImg); drawnow;
    hImg = findobj(ha,'Type','Image');
    set(hImg,'PickableParts','all','HitTest','on');

    % initialize
    bw = false(size(inor,1), size(inor,2), 0);
    ip = inor;
    imshow(ip); ha = gca; drawnow;

    % Let user draw one-or-more polygons. Finish by pressing Enter (or right-click)
    keepDrawing = true;
    while keepDrawing
        figure(figImg); axes(ha); drawnow;               % ensure correct focus
        roi = drawpolygon(ha);                           % interactive draw
        if isempty(roi.Position)                         % user canceled (Esc)
            break
        end
        mask = createMask(roi, hImg);                    % mask same size as image
        bw(:,:,end+1) = mask;
        ip = inpaintExemplar(ip, mask);                  % update image
        imshow(ip); ha = gca; drawnow;

        % Ask whether to add another polygon
        choice = questdlg('Add another polygon?', 'Continue', 'Yes', 'No', 'No');
        keepDrawing = strcmp(choice,'Yes');
    end

    % Create control panel AFTER drawing is finished
    hf = figure('Color','w','Position',[1200 600 250 70], ...
                'MenuBar','none','ToolBar','none','Name','Draw poly');
    uicontrol(hf,'Style','pushbutton','FontSize',12,'FontWeight','bold', ...
             'String','Process Masks','Position',[15 10 110 50], ...
             'Callback',@(~,~) processMasks(bw, ip));
    uicontrol(hf,'Style','pushbutton','FontSize',12,'FontWeight','bold', ...
             'String','Close','Position',[135 10 100 50], ...
             'Callback',@(~,~) close(hf));

    % --- local helper to process masks (edit as needed)
    function processMasks(masks, currentIp)
        % Example: display number of masks and show combined mask
        n = size(masks,3);
        disp(['Number of masks: ' num2str(n)]);
        if n>0
            combined = any(masks,3);
            figure; imshow(combined); title('Combined mask');
        end
        % Insert your processing here (save/export/apply more inpainting, etc.)
    end
end

