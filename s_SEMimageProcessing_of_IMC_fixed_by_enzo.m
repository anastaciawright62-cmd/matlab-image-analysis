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

inpaintFlag = 1;

% inor = mat2gray(img(:,:,1));%0.45,0.03);
inor = normalise(img(:,:,1),0.45,0.03);
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

ib = bilateralFilterxray(ip,ip,12,0.1,12,0.1);
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

method = 'thresholding';% 'kmeans'

if strcmpi(method,'kmeans')

    L = imsegkmeans(im2single(ib),4); % image clustering using K means
    % L = imsegkmeans(I,k) segments image I into k clusters by performing
    % k-means clustering and returns the segmented labeled output in L. J =
    % im2single(I) converts the grayscale, RGB, or binary image I to single,
    % rescaling or offsetting the data as necessary.
    bw = L==2; % create a mask using the segmented phase; In this case, it is Fe-rich phase
    bw2 = bwareaopen(bw,3); % round(1/pixsize^2) Remove small objects smaller than 1 um^2
    % BW2 = bwareaopen(BW,P) removes all connected components (objects) that
    % have fewer than P pixels from the binary image BW, producing another
    % binary image, BW2. This operation is known as an area opening.
elseif strcmpi(method,'thresholding')
    
    L = imbinarize(ib,0.71);
    bw2 = bwareaopen(bw,3); % round(1/pixsize^2) Remove small objects smaller than 1 um^2
    % BW2 = bwareaopen(BW,P) removes all connected components (objects) that
    % have fewer than P pixels from the binary image BW, producing another
    % binary image, BW2. This operation is known as an area opening.
else
    disp('Method not recognised');
end

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

%% Saving results 
%Nested function
function [ip bw] = improcess

ha=gca;
% ax = gca returns the current axes (or standalone visualization) in the current figure. Use ax to get and set properties of the current axes. 
% If there are no axes or charts in the current figure, then gca creates a Cartesian axes object.
hnor = findobj(ha,'Type','Image');
% h = findobj returns the graphics root object and all of its descendants.
% h = findobj(prop,value) returns all objects in the hierarchy that have their property prop set to value.
% h = findobj('-not',prop,value) returns all objects whose specified property is not set to the specified value.
inor = hnor.CData;

bw = false(size(inor,1), size(inor,2), 0);   % empty 3-D logical stack
bw = roipoly;


ip = inpaintExemplar(inor,bw(:,:,end));
% J = inpaintExemplar(I,mask) fills specific regions in the input image using the exemplar-based inpainting method. 
% mask is a logical image that denotes the target regions in the image to be filled using inpainting.
imshow(ip); ha=gca; % update the handle

hf  = figure('color','w','position',[1200 600 250 70]);
% figure creates a new figure window using default property values. The resulting figure is the current figure.
% figure(Name,Value) modifies properties of the figure using one or more name-value pair arguments. For example, figure('Color','white') sets the background color to white.
% f = figure(___) returns the Figure object. Use f to query or modify properties of the figure after it is created.
% the figure window is positioned 1200 pixels to the right and 600 pixels above the bottom left corner of the primary display, and is 250 pixels wide and 70 pixels tall.

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
                % c = uicontrol creates a push button (the default user interface control) in the current figure, and returns the UIControl object. 
                % If a figure does not exist, then MATLAB® calls the figure function to create one.
                % c = uicontrol(Name,Value) creates a user interface control with property values 
                % specified using one or more name-value pair arguments. For example, 'Style','checkbox' creates a check box.
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
