clear all
close all
clc

%% ============================================================
% Ecad_GATA6_S100A2_Analysis_public.m
%
% Example code for cell-level E-cadherin, GATA6, and S100A2
% analysis in a single ROI.
%
% Required inputs:
%   1. H&E image
%   2. Three-channel fluorescence image
%      Channel 1 = S100A2
%      Channel 2 = E-cadherin
%      Channel 3 = DAPI
%   3. GATA6 fluorescence image
%   4. HoverNet nuclear segmentation (.mat)
%   5. Trained ResNet50 model for E-cadherin classification
%
% E-cadherin:
%   ResNet50-based classification into four patterns:
%   P1, P2, P3, and P4.
%
% GATA6:
%   Nuclear GATA6 expression is classified as negative, weak,
%   or strong using predefined fluorescence intensity and
%   positive-area thresholds.
%
% S100A2:
%   Nuclear S100A2 expression is classified as negative or
%   positive using predefined fluorescence intensity and
%   positive-area thresholds.
%
% Fluorescence intensity thresholds used in this example are
% predefined representative values. The threshold determination
% procedure is described below.
%
% Classification results are displayed in the MATLAB Command
% Window. Annotated images are saved in the "output" subfolder.
%% ============================================================

%% Input files
dataDir=fileparts(mfilename('fullpath'));

HEfile=fullfile(dataDir,'example_HE.jpg');
FLfile=fullfile(dataDir,'example_FL.png');
GATAfile=fullfile(dataDir,'example_GATA6.tif');
HoverNetFile=fullfile(dataDir,'example_HoverNet.mat');
modelFile=fullfile(dataDir,'example_Ecad_ResNet50.mat');

%% Output directory
outputDir=fullfile(dataDir,'output');

if ~exist(outputDir,'dir')
    mkdir(outputDir)
end

%% Load trained E-cadherin classifier
load(modelFile,'net');

%% Analysis parameters
patchsize=128;
scalesize=[224 224];

dilsize=5;
seNuc=strel('diamond',dilsize);

% Cell-selection criteria:
% Area_th(1) = E-cadherin-positive area / nuclear area (%)
% Area_th(2) = DAPI-positive area / nuclear area (%)
% Area_th(3) = minimum nuclear area (pixels)
% Area_th(4) = minimum bounding-box width (pixels)
% Area_th(5) = minimum bounding-box height (pixels)
Area_th=[33.7 16.5 552 20 20];

% GATA6-positive area threshold
gata6_area_Th=0.2;

% S100A2-positive area threshold
s100_area_Th=0.5;

%% ============================================================
% Representative fluorescence intensity thresholds
%
% In the original analysis, intensity thresholds were determined
% independently for each TMA core. Otsu's method (graythresh)
% was applied to S100A2, E-cadherin, DAPI, and GATA6.
%
% For GATA6, multithresh was additionally applied to determine
% the strong-positive intensity threshold.
%
% The median threshold across all TMA cores from the same slide
% was used as the representative threshold for subsequent
% cell-level analysis.
%
% For this single-ROI example, representative threshold values
% determined by the above procedure are provided directly.
%% ============================================================

S100A2_th=0.2804;
Ecad_th=0.1373;
DAPI_th=0.0706;
GATA6_th=0.0667;
GATA6_strong_th=0.5333;

%% Read images and HoverNet segmentation
HE=imread(HEfile);
FL=imread(FLfile);
gata6_img=imread(GATAfile);
hovermat=load(HoverNetFile);

if ~isfield(hovermat,'inst_uid') || ~isfield(hovermat,'inst_map') || ...
        ~isfield(hovermat,'inst_centroid')
    error('Required HoverNet variables are missing.')
end

labelid=hovermat.inst_uid;
CyCx=hovermat.inst_centroid;

H_src=size(FL,1);
W_src=size(FL,2);

%% Check input data
if size(FL,3)~=3
    error('The fluorescence image must contain three channels.')
end

if size(HE,1)~=H_src || size(HE,2)~=W_src
    error('H&E and fluorescence image sizes do not match.')
end

if size(gata6_img,1)~=H_src || size(gata6_img,2)~=W_src
    error('Fluorescence and GATA6 image sizes do not match.')
end

if size(hovermat.inst_map,1)~=H_src || size(hovermat.inst_map,2)~=W_src
    error('Fluorescence image and HoverNet segmentation sizes do not match.')
end

%% Fluorescence masks
s100_bw=imbinarize(FL(:,:,1),S100A2_th);
ecad_bw=imbinarize(FL(:,:,2),Ecad_th);
nuc_img=imbinarize(FL(:,:,3),DAPI_th);

ecad_bw=imdilate(ecad_bw,seNuc);
ecad_bw=imfill(ecad_bw,'holes');

% Restrict S100A2-positive signal to the DAPI-positive region
s100_bw=s100_bw & nuc_img;

%% GATA6 intensity classification
% imquantize with two thresholds:
%   1 = below positive threshold
%   2 = weak-positive intensity range
%   3 = strong-positive intensity range
gata6_level=[GATA6_th GATA6_strong_th]*255;
gata6_label=imquantize(gata6_img,gata6_level);

%% Result buffers
NumCandidate=length(labelid);

% Column 1 = predicted E-cadherin class
% Columns 2-5 = ResNet scores for P1-P4
Ecadtype=zeros(NumCandidate,5);

% GATA6: 0=negative, 1=weak, 2=strong
gata6type=zeros(NumCandidate,1);

% S100A2: 0=negative, 1=positive
s100type=zeros(NumCandidate,1);

% Bounding boxes of successfully analyzed nuclei
BBox=zeros(NumCandidate,4);

ci=1;

%% ============================================================
% Cell-level analysis
%% ============================================================

for o=1:NumCandidate

    %% HoverNet nuclear mask
    BW=(hovermat.inst_map==labelid(o));
    BB=regionprops(BW,'BoundingBox');

    if isempty(BB)
        continue
    end

    bbox=round(BB(1).BoundingBox);
    bbox=checkbbox_v2(bbox,H_src,W_src);

    %% Nuclear-centered patch
    Sx=CyCx(o,2)-patchsize/2;
    Sy=CyCx(o,1)-patchsize/2;

    if Sx<1, Sx=1; end
    if Sy<1, Sy=1; end
    if Sx+patchsize>W_src, Sx=W_src-patchsize; end
    if Sy+patchsize>H_src, Sy=H_src-patchsize; end

    Sx=round(Sx);
    Sy=round(Sy);

    FL_patch=FL(Sy:Sy+patchsize,Sx:Sx+patchsize,:);
    nucimg_patch=nuc_img(Sy:Sy+patchsize,Sx:Sx+patchsize);
    ecad_bw_patch=ecad_bw(Sy:Sy+patchsize,Sx:Sx+patchsize);
    gata6_label_patch=gata6_label(Sy:Sy+patchsize,Sx:Sx+patchsize);
    s100_bw_patch=s100_bw(Sy:Sy+patchsize,Sx:Sx+patchsize);

    %% Local bounding box
    bbox_local=bbox;
    bbox_local(1)=bbox_local(1)-Sx;
    bbox_local(2)=bbox_local(2)-Sy;
    bbox_local=checkbbox_v2(bbox_local,patchsize,patchsize);

    %% Nuclear area
    HE_nuc_area=sum(BW(:));

    if HE_nuc_area<=0
        continue
    end

    %% DAPI occupancy
    sig=nucimg_patch(...
        bbox_local(2):bbox_local(2)+bbox_local(4)-1,...
        bbox_local(1):bbox_local(1)+bbox_local(3)-1);

    dapi_area=sum(sig(:));
    dapi_occupancy_rate=dapi_area/HE_nuc_area;

    %% E-cadherin occupancy
    sig=ecad_bw_patch(...
        bbox_local(2):bbox_local(2)+bbox_local(4)-1,...
        bbox_local(1):bbox_local(1)+bbox_local(3)-1);

    ecad_area=sum(sig(:));
    ecad_occupancy_rate=ecad_area/HE_nuc_area;

    %% GATA6 signal
    sig=gata6_label_patch(...
        bbox_local(2):bbox_local(2)+bbox_local(4)-1,...
        bbox_local(1):bbox_local(1)+bbox_local(3)-1);

    gata6_weak_area=sum(sig(:)==2);
    gata6_strong_area=sum(sig(:)==3);

    gata6_weak_area_rate=gata6_weak_area/HE_nuc_area;
    gata6_strong_area_rate=gata6_strong_area/HE_nuc_area;

    %% S100A2 signal
    sig=s100_bw_patch(...
        bbox_local(2):bbox_local(2)+bbox_local(4)-1,...
        bbox_local(1):bbox_local(1)+bbox_local(3)-1);

    s100_area=sum(sig(:));
    s100_occupancy_rate=s100_area/HE_nuc_area;

    %% Cell-selection criteria
    if ~(ecad_occupancy_rate*100>Area_th(1) && ...
            dapi_occupancy_rate*100>Area_th(2) && ...
            HE_nuc_area>Area_th(3) && ...
            bbox_local(3)>=Area_th(4) && ...
            bbox_local(4)>=Area_th(5))
        continue
    end

    %% E-cadherin pattern classification
    % ResNet50 output:
    %   1 = P1
    %   2 = P2
    %   3 = P3
    %   4 = P4
    FL_patch_input=imresize(FL_patch,scalesize);

    [Y,scores]=classify(net,FL_patch_input);
    predY=str2double(cellstr(Y));

    Ecadtype(ci,:)=[predY scores];

    %% GATA6 classification
    % Cells with >20% GATA6-positive nuclear area are classified
    % as positive. Weak/strong status is determined by the
    % dominant GATA6 intensity category.
    if (gata6_weak_area_rate+gata6_strong_area_rate)>gata6_area_Th
        if gata6_weak_area_rate<gata6_strong_area_rate
            gata6type(ci)=2;
        else
            gata6type(ci)=1;
        end
    else
        gata6type(ci)=0;
    end

    %% S100A2 classification
    % Cells with >50% S100A2-positive nuclear area are classified
    % as S100A2-positive.
    if s100_occupancy_rate>s100_area_Th
        s100type(ci)=1;
    else
        s100type(ci)=0;
    end

    BBox(ci,:)=bbox;
    ci=ci+1;
end

%% Keep successfully analyzed cells
Numofcell=ci-1;

Ecadtype=Ecadtype(1:Numofcell,:);
gata6type=gata6type(1:Numofcell,:);
s100type=s100type(1:Numofcell,:);
BBox=BBox(1:Numofcell,:);

%% E-cadherin pattern counts
EcadCount=zeros(1,4);

for ei=1:4
    EcadCount(ei)=sum(Ecadtype(:,1)==ei);
end

TotalEcad=sum(EcadCount);

%% E-cadherin x GATA6 x S100A2 counts
GataNega=(gata6type==0);
GataWeak=(gata6type==1);
GataStrong=(gata6type==2);

s100Posi=(s100type==1);
s100Nega=~s100Posi;

% Dimensions:
%   1 = E-cadherin pattern (P1-P4)
%   2 = GATA6 status (negative / weak / strong)
%   3 = S100A2 status (negative / positive)
counts=zeros(4,3,2);

for ei=1:4
    flg_ecad=(Ecadtype(:,1)==ei);

    counts(ei,1,1)=sum(flg_ecad&GataNega&s100Nega);
    counts(ei,1,2)=sum(flg_ecad&GataNega&s100Posi);
    counts(ei,2,1)=sum(flg_ecad&GataWeak&s100Nega);
    counts(ei,2,2)=sum(flg_ecad&GataWeak&s100Posi);
    counts(ei,3,1)=sum(flg_ecad&GataStrong&s100Nega);
    counts(ei,3,2)=sum(flg_ecad&GataStrong&s100Posi);
end

%% Display classification results
EcadName={'P1','P2','P3','P4'};
GataName={'GATA6-neg','GATA6-weak','GATA6-strong'};
S100Name={'S100A2-neg','S100A2-pos'};

fprintf('\n=============================================\n');
fprintf('E-CADHERIN PATTERN CLASSIFICATION\n');
fprintf('=============================================\n');

for ei=1:4
    if TotalEcad>0
        percentage=100*EcadCount(ei)/TotalEcad;
    else
        percentage=0;
    end

    fprintf('%s : %5d (%5.1f%%)\n',EcadName{ei},EcadCount(ei),percentage);
end

fprintf('---------------------------------------------\n');
fprintf('TOTAL : %5d\n',TotalEcad);

fprintf('\n=============================================\n');
fprintf('E-CADHERIN x GATA6 x S100A2\n');
fprintf('=============================================\n');

for ei=1:4
    fprintf('\n%s\n',EcadName{ei});

    for gi=1:3
        for si=1:2
            fprintf('  %-12s %-11s : %5d\n',...
                GataName{gi},S100Name{si},counts(ei,gi,si));
        end
    end
end

%% ============================================================
% Save classification images
%% ============================================================

%% E-cadherin
% P1=green, P2=magenta, P3=cyan, P4=yellow
cmap_Ecad=[
    0   255   0;
    255   0 255;
    0   255 255;
    255 255   0];

colors_Ecad=zeros(Numofcell,3);

for pi=1:4
    flg=(Ecadtype(:,1)==pi);
    colors_Ecad(flg,:)=repmat(cmap_Ecad(pi,:),sum(flg),1);
end

RGB_Ecad_HE=insertShape(HE,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_Ecad);

RGB_Ecad_FL=insertShape(FL,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_Ecad);

imwrite(RGB_Ecad_HE,fullfile(outputDir,'example_Ecad_HE.jpg'));
imwrite(RGB_Ecad_FL,fullfile(outputDir,'example_Ecad_FL.jpg'));

%% GATA6
% negative=cyan, weak=yellow, strong=red
cmap_GATA6=[
    0   255 255;
    255 255   0;
    255   0   0];

colors_GATA6=zeros(Numofcell,3);

for gi=0:2
    flg=(gata6type==gi);
    colors_GATA6(flg,:)=repmat(cmap_GATA6(gi+1,:),sum(flg),1);
end

RGB_GATA6_HE=insertShape(HE,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_GATA6);

RGB_GATA6=insertShape(gata6_img,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_GATA6);

imwrite(RGB_GATA6_HE,fullfile(outputDir,'example_GATA6_HE.jpg'));
imwrite(RGB_GATA6,fullfile(outputDir,'example_GATA6.jpg'));

%% S100A2
% negative=cyan, positive=yellow
cmap_S100A2=[
    0   255 255;
    255 255   0];

colors_S100A2=zeros(Numofcell,3);

for si=0:1
    flg=(s100type==si);
    colors_S100A2(flg,:)=repmat(cmap_S100A2(si+1,:),sum(flg),1);
end

RGB_S100A2_HE=insertShape(HE,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_S100A2);

RGB_S100A2_FL=insertShape(FL,'Rectangle',BBox,...
    'LineWidth',2,'Color',colors_S100A2);

imwrite(RGB_S100A2_HE,fullfile(outputDir,'example_S100A2_HE.jpg'));
imwrite(RGB_S100A2_FL,fullfile(outputDir,'example_S100A2_FL.jpg'));