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
pixsize = str2double(answer{1});  % in microns
% converts a character array or string scalar (first inupt of answer)) to a numeric matrix. The input can include spaces, commas, 
% and semicolons to indicate separate elements. If str2num cannot parse the input as numeric values, then it returns an empty matrix.
% Users can enter scalar or vector values into inputdlg text edit fields. 
% MATLAB® stores the input as a cell array of character vectors. 
% Convert a member of the input cell array to a number, using str2num.
%% 2.1 Process the image if SEM image (gray scale)

img = img(1:700,:); % crop the infomation bar at the bottom of the image

choice = questdlg('Do you want to manually remove features from segmentation?', ...
                  'Inpainting','Yes','No','Yes');

inpaintFlag = strcmp(choice,'Yes');

inor = normalise(img(:,:,1),0.45,0.03);
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

    % Defensive: bw may be empty (no polygons), 2-D (one mask), or 3-D stack.
    if isempty(bw)
        % No polygons drawn: inform the user and continue with no inpaint masks
        msg = 'No polygon masks were drawn. Continuing with segmentation only.';
        warning(msg);
        try
            msgbox(msg,'No masks','warn');
        catch
            % ignore if running headless
        end
        bwComb = false(size(inor));    % safe all-false mask
    else
        bw = logical(bw);
        if ndims(bw) == 2
            bwComb = bw;               % single mask
        else
            bwComb = any(bw,3);        % union of all mask layers
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

method = questdlg('Choose segmentation method','Segmentation','Thresholding','Kmeans','Thresholding');
if isempty(method)
    method = 'thresholding'; % default if user closes dialog
end

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
    bw2 = bwareaopen(L,3); % round(1/pixsize^2) Remove small objects smaller than 1 um^2
    % BW2 = bwareaopen(BW,P) removes all connected components (objects) that
    % have fewer than P pixels from the binary image BW, producing another
    % binary image, BW2. This operation is known as an area opening.
else
    disp('Method not recognised');
end

bw2 = imclearborder(bw2);%not sure if I should remove them

s = regionprops(bw2,'Centroid','Area','MajorAxisLength','MinorAxisLength','Perimeter','Orientation','PixelIdxList','Solidity');

if isempty(s)
    warning('No objects found in bw2.');
    IMCcent = [];
    IMCarea = [];
    IMClength = [];
    IMCwidth = [];
    IMCasprat = [];
    IMCareafrac = 0;
    IMCperimeter = [];
    IMCcircularity = [];
    IMCorient_deg = [];
    IMCsolidity = [];
    IMCECD = [];
    IMCcount = 0;
    imageArea_um2 = numel(bw2) * pixsize^2;
    IMCdensity_um2 = 0;
else
    %Extract raw pixel properties
    IMCcent = cat(1,s.Centroid);           % Nx2 (px)
    IMCarea = cat(1,s.Area);               % px
    IMClength = cat(1,s.MajorAxisLength);    % px
    IMCwidth = cat(1,s.MinorAxisLength);    % px
    IMCperimeter = cat(1,s.Perimeter);          % px
    IMCorient_deg = cat(1,s.Orientation);        % degrees
    IMCsolidity   = cat(1,s.Solidity);          % unitless
    %remove objects with 0 width or perimeter
    valid = (IMCwidth > 0) & (IMCperimeter > 0);   % remove degenerate objects
    if ~any(valid)
        warning('All detected objects removed by validity filter.');
        % produce empty outputs consistently
        IMCcent = []; IMCarea = []; IMClength = []; IMCwidth = []; IMCperimeter = [];
        IMCasprat = []; IMCcircularity = []; IMCorient_deg = []; IMCsolidity = [];
        IMCECD = []; IMCcount = 0;
        imageArea_um2 = numel(bw2) * pixsize^2;
        IMCdensity_um2 = 0;
        mean_ip = []; mean_ib = []; mean_orig = [];
    else
        % apply filter to arrays and to s so PixelIdxList stays aligned
        s = s(valid);
        IMCcent = IMCcent(valid,:);
        IMCarea = IMCarea(valid);
        IMClength = IMClength(valid);
        IMCwidth = IMCwidth(valid);
        IMCperimeter = IMCperimeter(valid);
        IMCorient_deg = IMCorient_deg(valid);
        IMCsolidity = IMCsolidity(valid);
        %now computing the rest
        IMCasprat = IMClength./IMCwidth;
        IMCcircularity = 4*pi*IMCarea ./ (IMCperimeter.^2);
        IMCECD = sqrt(4 * IMCarea / pi);%equivalent circular diameter (good size measurement)
        IMCareafrac = sum(IMCarea)/numel(bw2);%numerator excludes border-touching objects
        IMCcount = numel(IMCarea);
        % Sanity check: total particle area cannot exceed image area
        assert(sum(IMCarea) <= numel(bw2) + 1e-6, 'Total IMC area (px) > image area (px) — segmentation bug?');
    end
    %conversions to microns
    IMCarea_um2 = IMCarea * (pixsize^2);
    IMCECD_um = IMCECD * pixsize; % if ECD computed in px -> convert to µm
    imageArea_um2 = numel(bw2) * pixsize^2;
    IMCdensity_um2 = IMCcount / imageArea_um2;
    IMCdensity_per_mm2 = IMCdensity_um2 * 1e6;
    % ---- per-object means: require ib, compute means ----
    if ~exist('ib','var') || ~isequal(size(ib), size(bw2))
        error('Filtered image ''ib'' is missing or wrong size: mean_ib cannot be computed.');
    end
    
    % original image must exist (let it error if missing)
    orig = img;
    if size(orig,3) > 1
        orig = im2gray(orig);
    end
    
    nobj = numel(s);
    mean_ib = nan(nobj,1);
    mean_ip = nan(nobj,1);
    mean_orig = nan(nobj,1);
    
    for k = 1:nobj
        idx = s(k).PixelIdxList;
        if ~isempty(idx)
            mean_ib(k) = mean(ib(idx));                % filtered
            mean_ip(k) = mean(ip(idx));
            mean_orig(k) = mean(orig(idx));            % raw original
        end
    end
    
end

figure, imshow(bw2);

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

% Draw polygon (may return empty if user cancels)
bw = roipoly;

% If no masks drawn, skip inpainting and keep original image
if isempty(bw)
    ip = inor;                     % no change to image
else
    % If bw is 2-D (one mask) it's fine; if 3-D handle accordingly
    bw = logical(bw);
    if ndims(bw) == 2
        ip = inpaintExemplar(inor, bw);       % single mask
    else
        ip = inpaintExemplar(inor, bw(:,:,end)); % most recent mask
    end
end

% redraw into same axes so ha remains valid
imshow(ip, 'Parent', ha);
drawnow;
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
%% --- Save everything shown in screenshot + related outputs ---

% Build output folder
fullImagePath = fullfile(folder, file);
[imgFolder, imgName, ~] = fileparts(fullImagePath);
outFolder = fullfile(imgFolder, [imgName '_results']);
if ~exist(outFolder,'dir'), mkdir(outFolder); end

% Helper: safe save variable if exists
function safeSaveVar(varname, outpath)
    if evalin('caller', ['exist(''' varname ''',''var'')'])
        v = evalin('caller', varname);
        save(outpath, varname, '-v7.3'); % -v7.3 for large arrays
    else
        warning('Variable %s not found; skipping save %s', varname, outpath);
    end
end

%save s.mat
fullImagePath = fullfile(folder, file);
[imgFolder, imgName, ~] = fileparts(fullImagePath);
outFolder = fullfile(imgFolder, [imgName '_results']);
if ~exist(outFolder,'dir'), mkdir(outFolder); end
% save s struct (use -v7.3 if large). This is the exact geometry of every
% detected particle
save(fullfile(outFolder, [imgName '_s.mat']), 's', '-v7.3');

% Save a logical mask MAT too (smaller than csv)
if exist('bw2','var'), save(fullfile(outFolder,[imgName '_mask.mat']), 'bw2'); end

% --- Save CSV versions for numeric vectors / tables where sensible ---
% IMCarea, IMClength, IMCwidth are vectors -> CSV
if exist('IMCarea','var')
    try writematrix(IMCarea(:), fullfile(outFolder,[imgName '_IMCarea.csv'])); catch; warning('Could not write IMCarea CSV'); end
end
if exist('IMClength','var')
    try writematrix(IMClength(:), fullfile(outFolder,[imgName '_IMClength.csv'])); catch; warning('Could not write IMClength CSV'); end
end
if exist('IMCwidth','var')
    try writematrix(IMCwidth(:), fullfile(outFolder,[imgName '_IMCwidth.csv'])); catch; warning('Could not write IMCwidth CSV'); end
end
if exist('IMCasprat','var')
    try writematrix(IMCasprat(:), fullfile(outFolder,[imgName '_IMCasprat.csv'])); catch; warning('Could not write IMCasprat CSV'); end
end
if exist('IMCareafrac','var')
    try writematrix(IMCareafrac, fullfile(outFolder,[imgName '_IMCareafrac.csv'])); catch; warning('Could not write IMCareafrac CSV'); end
end
if exist('IMCcent','var')
    try writematrix(IMCcent, fullfile(outFolder,[imgName '_IMCcent.csv'])); catch; warning('Could not write IMCcent CSV'); end
end

% --- Create and save histograms/figures shown in screenshot ---
% helper to save invisible figure
function fh = makeHiddenFig()
    fh = figure('Visible','off','Units','pixels','Position',[100 100 800 600]);
end

% hist(IMCarea)
if exist('IMCarea','var') && ~isempty(IMCarea)
    fh = makeHiddenFig(); histogram(IMCarea,40); xlabel('IMC area (px)'); title('hist(IMCarea)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCarea).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCarea).fig'])); close(fh);
end

% hist(IMCareafrac)  -- ordinarily single scalar, histogram can be trivial; make a small bar/figure
if exist('IMCareafrac','var')
    fh = makeHiddenFig(); bar(IMCareafrac); ylabel('area fraction'); title('hist(IMCareafrac)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCareafrac).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCareafrac).fig'])); close(fh);
end

% hist(IMCasprat)
if exist('IMCasprat','var') && ~isempty(IMCasprat)
    fh = makeHiddenFig(); histogram(IMCasprat,40); xlabel('aspect ratio'); title('hist(IMCasprat)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCasprat).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCasprat).fig'])); close(fh);
end

% hist(IMCcent) -- screenshot shows a figure named hist(IMCcent); we save a scatter of centroids (x vs y)
if exist('IMCcent','var') && ~isempty(IMCcent)
    fh = makeHiddenFig(); scatter(IMCcent(:,1), IMCcent(:,2), 6, 'filled'); axis ij; xlabel('x (px)'); ylabel('y (px)'); title('IMC centroids (px)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCcent).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCcent).fig'])); close(fh);
end

% hist(IMClength)
if exist('IMClength','var') && ~isempty(IMClength)
    fh = makeHiddenFig(); histogram(IMClength,40); xlabel('Major axis length (px)'); title('hist(IMClength)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMClength).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMClength).fig'])); close(fh);
end

% hist(IMCwidth)
if exist('IMCwidth','var') && ~isempty(IMCwidth)
    fh = makeHiddenFig(); histogram(IMCwidth,40); xlabel('Minor axis length (px)'); title('hist(IMCwidth)'); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCwidth).png'])); saveas(fh, fullfile(outFolder,[imgName '_hist(IMCwidth).fig'])); close(fh);
end

% Save 'ib' data figure (as .fig and PNG normalized)
if exist('ib','var')
    try
        fh = makeHiddenFig(); imshow(mat2gray(ib)); title('ib (normalized for display)'); saveas(fh, fullfile(outFolder,[imgName '_ib.png'])); saveas(fh, fullfile(outFolder,[imgName '_ib.fig'])); close(fh);
    catch, warning('Could not save ib image preview'); end
end

% Save 'ip' data figure (as .fig and PNG normalized)
if exist('ip','var')
    try
        fh = makeHiddenFig(); imshow(mat2gray(ip)); title('ip (normalized for display)'); saveas(fh, fullfile(outFolder,[imgName '_ip.png'])); saveas(fh, fullfile(outFolder,[imgName '_ip.fig'])); close(fh);
    catch, warning('Could not save ip image preview'); end
end

% Save 'normalised' image (inor) as figure
if exist('inor','var')
    try
        fh = makeHiddenFig(); imshow(mat2gray(inor)); title('normalised'); saveas(fh, fullfile(outFolder,[imgName '_normalised.png'])); saveas(fh, fullfile(outFolder,[imgName '_normalised.fig'])); close(fh);
    catch, warning('Could not save normalised image'); end
end

% Save 'segmented' figure: overlay bw2 on original (if available)
if exist('bw2','var')
    try
        fh = makeHiddenFig(); imshow(mat2gray(img)); hold on;
        h = imshow(bw2); set(h,'AlphaData',0.3); title('segmented'); saveas(fh, fullfile(outFolder,[imgName '_segmented.png'])); saveas(fh, fullfile(outFolder,[imgName '_segmented.fig'])); close(fh);
    catch, warning('Could not save segmented figure'); end
end

% --- Save Excel file if you want a sheet (e.g. layer62-junction... style) ---
% create an Excel workbook with the per-object table and a summary sheet
if exist('IMCarea','var')
    try
        T = table((1:numel(IMCarea))', IMCcent(:,1), IMCcent(:,2), IMCarea(:), IMClength(:), IMCwidth(:), IMCasprat(:), IMCECD(:), mean_ib(:), ...
            'VariableNames', {'id','centroid_x_px','centroid_y_px','area_px','majorAxis_px','minorAxis_px','aspectRatio','ECD_px','mean_ib'});
        writetable(T, fullfile(outFolder, [imgName '_IMC_table.xlsx']), 'Sheet', 'PerObject');
        % summary sheet
        summaryT = table({imgName}, IMCcount, IMCareafrac, IMCdensity_um2, pixsize, min_area_um2, 'VariableNames', {'image','n_particles','area_fraction','density_per_um2','pixsize_um','min_area_um2'});
        writetable(summaryT, fullfile(outFolder, [imgName '_IMC_table.xlsx']), 'Sheet', 'Summary');
    catch, warning('Could not write Excel workbook (maybe file locked).'); end
end

% --- Optionally save full ib/ip to CSV (disabled by default) ---
save_full_image_csv = false;
if save_full_image_csv
    if exist('ib','var')
        try writematrix(ib, fullfile(outFolder,[imgName '_ib.csv'])); catch, warning('Could not write ib CSV'); end
    end
    if exist('ip','var')
        try writematrix(ip, fullfile(outFolder,[imgName '_ip.csv'])); catch, warning('Could not write ip CSV'); end
    end
end

% Save the single scalar summary CSV + JSON (same as earlier block) ---
summaryStruct = struct();
summaryStruct.image = imgName;
summaryStruct.date = datestr(now,'yyyy-mm-dd HH:MM:SS');
summaryStruct.n_particles = exist('IMCcount','var') * IMCcount; % will be 0 if missing
if exist('IMCareafrac','var'), summaryStruct.area_fraction = IMCareafrac; else summaryStruct.area_fraction = NaN; end
if exist('IMCdensity_um2','var'), summaryStruct.density_per_um2 = IMCdensity_um2; summaryStruct.density_per_mm2 = IMCdensity_um2*1e6; else summaryStruct.density_per_um2 = NaN; summaryStruct.density_per_mm2 = NaN; end
summaryStruct.pixel_size_um = exist('pixsize','var') * pixsize;
summaryStruct.min_area_um2 = exist('min_area_um2','var') * min_area_um2;
summaryStruct.notes = 'Saved per-object tables, histograms, mask and images.';

% CSV summary
Snames = fieldnames(summaryStruct);
Svals = cellfun(@(f) summaryStruct.(f), Snames, 'uni', false);
summaryTable = cell2table([Snames, Svals],'VariableNames',{'Metric','Value'});
writetable(summaryTable, fullfile(outFolder, [imgName '_summary.csv']),'WriteVariableNames',false);

% JSON summary
try
    jsonText = jsonencode(summaryStruct);
    fid = fopen(fullfile(outFolder, [imgName '_summary.json']),'w');
    if fid ~= -1
        fwrite(fid, jsonText, 'char');
        fclose(fid);
    end
catch
    warning('Could not write JSON summary (jsonencode maybe unavailable).');
end

fprintf('Saved full set of outputs to folder:\n%s\n', outFolder);