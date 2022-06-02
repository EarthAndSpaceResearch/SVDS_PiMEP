%% Triple Point Collocation Example Use Case:
%  Monthly, L3 SMAP, SMOS, and Argo data
%
% j. anderson (janderson@esr.org)
% last edit 8 March 2022

%% Prepare workspace
clear all
close all
clc

%% Define variables -- all user defined variables here

% % Add paths of code and packages to be used if not already on path
% % Remote machine
% addpath(genpath('/vol0/janderson/scripts'));
% addpath(genpath('/vol0/janderson/triplepoint_usecasedata'));

% Output directory, Local
outputDir = '/Users/jesse/Documents/ESR/PiMEP/Data/TriplePointCollocation';

% Input data directories, Local
smosDir = '/Volumes/T7/ESR/Datasets/SMOS/SMOS_SSS_L3_v317_CATDS_CPDC/Ocean_products/GRIDDED/MonthlyOnly'; % location of smos data
smapDir = '/Volumes/T7/ESR/Datasets/SMAP/L3/RSS/V4/monthly/SCI/AllData'; % location of smap data
argoDir = '/Volumes/T7/ESR/Datasets/Argo/RG_ArgoClimatology'; % location of argo data

% % Output directory, Remote
% outputDir = '/vol0/janderson/analysis_output';
% 
% % Input data directories, Remote
% smosDir = '/vol0/janderson/triplepoint_usecasedata/smos/SMOS_SSS_L3_v317_CATDS_CPDC/Ocean_products/GRIDDED/MonthlyOnly'; % location of smos data
% smapDir = '/vol0/janderson/triplepoint_usecasedata/smap/L3/RSS/V4/monthly/SCI/AllData'; % location of smap data
% argoDir = '/vol0/janderson/triplepoint_usecasedata/argo/RG_ArgoClimatology'; % location of argo data

% Date to start and end analysis
% PiMEP integration should have error checking for valid dates for each dataset
startDate = datevec('4-1-2015'); % For this use case, SMAP has the latest start of 1 April 2015
endDate = datevec('1-31-2022'); % Use case maximum date 1-31-2022

% Set large difference flag
excludeLargeDiff = 1; %1 (true) yes exclude large difference (>5PSU), 0 for no, include all salinity data differences regardless of size

%% Parse inputs

% Create start year and month variables
startYear=startDate(1);
startMonth=startDate(2);

% Create end year and month variables
endYear=endDate(1);
endMonth=endDate(2);

% Determine how many months to run analysis
year1Months = (12-startMonth)+1;
yearLastMonths = endMonth; 
yearsBetweenMonths = (endYear-startYear-1)*12;
numMonths = year1Months+yearLastMonths+yearsBetweenMonths;
clear year1Months yearLastMonths yearsBetweenMonths

% Create list of month/year combinations
dateList = datevec(datetime(startDate) + calmonths(0:numMonths-1)); % Create list of year, months, days

%% Load First Dataset - SMOS

% Go to location of SMOS data
cd(smosDir);

% Get list of all data files
tempFileList = dir(['SM_RE07','*', '.nc']); % Get information on all delayed mode files in directory 
tempFileList=char(tempFileList.name); % Keep only the filenames
tempFileListOper = dir(['SM_OPER','*', '.nc']); %  Get information on all operational data files in directory
tempFileList=cat(1,tempFileList,tempFileListOper.name); % Keep only the filenames and add to end of file list
clear tempFileListOper

% Read in static variables from first file
fnTemp = tempFileList(1,:); % get first filename
smosLon = ncread(fnTemp,'lon'); % 694x1 degrees_east -179.8703 to 179.6109 0.5187 degree FillValue = 9.969209968386869e+36
smosLat = ncread(fnTemp,'lat'); % 292x1 degrees_north -83.5171 81.9831 FillValue = 9.969209968386869e+36
clear fnTemp

% Initialize and preallocate dynamic variables with matrices of NaNs
smosDateTime = nan(size(tempFileList,1),1);
smosS = nan(length(smosLon), length(smosLat),size(tempFileList,1));

% Load and write data from monthly files to dynamic variables
for fileLoop = 1:size(tempFileList,1)
    fnTemp = tempFileList(fileLoop,:); % Get name of next file
    
    % Get date and time from filename
    smosDateTime(fileLoop,1) = datenum([str2num(fnTemp(20:23)) str2num(fnTemp(24:25)) ...
        str2num(fnTemp(26:27)) str2num(fnTemp(29:30)) ...
        str2num(fnTemp(31:32)) str2num(fnTemp(33:34))]);
    
    % Get salinity data
    smosS(:,:,fileLoop) = ncread(fnTemp,'Mean_Sea_Surface_Salinity_1'); % 694x292 lon, lat 'Mean sea surface salinity' FillValue = -32767 PSU add_offset = 30 scale_factor = 0.001
    %smosS(:,:,fileLoop) = ncread(fnTemp,'Mean_Sea_Surface_Salinity_2')); % 694x292 lon, lat 'Mean sea surface salinity (rain corrected)' FillValue = -32767 PSU add_offset = 30 scale_factor = 0.001
    
    clear fnTemp
end

clear tempFileList fileLoop

% Replace fill values with NaN
smosS(smosS == -32767) = NaN;

% % Quick plot to ensure looks right
% [x,y]=meshgrid(smosLat, smosLon);
% figure(1)
% contourf(y, x, smosS(:,:,1), 'LevelStep', 0.2,'LineStyle','none')

% Standardize Lat - 90 to 90
%min(smosLat) %  -83.5171 0.3923-2.8199 degree spacing
%max(smosLat) % 81.9831

% Standardize Lon 0 to 360
%min(smosLon) %  -179.8703 0.5187 degree spacing
%max(smosLon) % 179.6109

% Adjust longitude values to be 0-360 degrees
tempIndex  = find(smosLon<0);
smosLon(tempIndex) = smosLon(tempIndex)+360;
clear tempIndex

% Sort adjusted longitude values and corresponding SMOS salinity data
[smosLon, tempIndex] = sort(smosLon);
smosS= smosS(tempIndex,:,:);
clear tempIndex

% % Quick plot to ensure looks right
% [x,y]=meshgrid(smosLat, smosLon);
% figure(2)
% contourf(y, x, smosS(:,:,1), 'LevelStep', 0.2,'LineStyle','none')

% Remove data outside month/year range of run
tempIndex = find(smosDateTime < datenum(startDate) | smosDateTime > datenum(endDate));
smosDateTime(tempIndex) = [];
smosS(:,:,tempIndex) = [];
clear tempIndex 

%% Load Second Dataset  - SMAP 70km

% Go to location of SMAP data
cd(smapDir);

% % Get list of all data files in directory
% One month does not have a file use method below instead
% tempFileList = dir(['RSS_smap_SSS_L3_monthly_','*', '.nc']);
% tempFileList=char(tempFileList.name);

% Create list of data files over processing time period
tempFileList = [repmat('RSS_smap_SSS_L3_monthly_',length(dateList),1),num2str(dateList(:,1)),...
    repmat('_',length(dateList),1),num2str(dateList(:,2),'%02i'),...
    repmat('_FNL_v04.0.nc',length(dateList),1)];

% Determine if all the data file names created exist in directory
tempFileExist = isfile(string(tempFileList));
if sum(tempFileExist)==0
    error('None of the filenames created exist in the current directory')
end

% Read in static variables from first file that exists
tempFileIndex = min(find(tempFileExist==1)); % Find index of first file that exists
fnTemp = tempFileList(tempFileIndex,:); % get first filename
smapLon = ncread(fnTemp,'lon'); % 1440x1 center longitude of grid cell degrees east 0 360 0.25degree
smapLat = ncread(fnTemp,'lat'); % 720x1 center latitude of grid cell degrees_north -90 90 0.25 degree
clear tempFileIndex fnTemp tempFileExist

% Initialize and preallocate dynamic variables with matrices of NaNs
smapTime = nan(numMonths,1); % 1x1 seconds since 2000-01-01T00:00:00Z reference time of analyzed variable field corresponding to center of the product time interval
smapS = nan(length(smapLon), length(smapLat), numMonths); % 1440x720 lon, lat sea_surface_salinity -9999 fill
smapGland = nan(length(smapLon), length(smapLat), numMonths); %1440x720 lon, lat Average land fraction weighted by antenna gain -9999 fill
smapFland = nan(length(smapLon), length(smapLat), numMonths); %1440x720 lon, lat  Average land fraction within 3dB contour -9999 fill
smapGice = nan(length(smapLon), length(smapLat), numMonths); %1440x720 lon, lat Average sea ice fraction (weighted by antenna gain) -9999 fill
smapSurtep = nan(length(smapLon), length(smapLat), numMonths);  %1440x720 lon, lat sea_surface_temperature -9999 fill

% Load and write data from monthly files to dynamic variables
for fileLoop = 1:size(tempFileList,1)
    if isfile(tempFileList(fileLoop,:)); % If the file exists read in data
        fnTemp = tempFileList(fileLoop,:); % Get name of next file
        smapTime(fileLoop,:) = ncread(fnTemp,'time'); %1x1 seconds since 2000-01-01T00:00:00Z reference time of analyzed variable field corresponding to center of the product time interval
        smapS(:,:,fileLoop) = ncread(fnTemp,'sss_smap'); % 1440x720 lon, lat sea_surface_salinity -9999 fill
        smapGland(:,:,fileLoop) = ncread(fnTemp, 'gland'); %1440x720 lon, lat Average land fraction weighted by antenna gain -9999 fill
        smapFland(:,:,fileLoop) = ncread(fnTemp,'fland'); %1440x720 lon, lat  Average land fraction within 3dB contour -9999 fill
        smapGice(:,:,fileLoop) = ncread(fnTemp,'gice'); %1440x720 lon, lat Average sea ice fraction (weighted by antenna gain) -9999 fill
        smapSurtep(:,:,fileLoop) = ncread(fnTemp,'surtep');  %1440x720 lon, lat sea_surface_temperature -9999 fill
        clear fnTemp
    else
        display(['Filename ', tempFileList(fileLoop,:), ' does not exist']);
    end % If file doesn't exist, do nothing. Already filled with nans from preallocation
end

clear tempFileList fileLoop

% Convert SMAP time to MATLAB datetime
% smapTime is seconds since seconds since 2000-01-01T 00:00:00Z
% smapTime is centered roughly on 15th of month
smapBaseTime = datevec('01-JAN-2000 00:00:00');
smapDateTime = datenum(datetime(smapBaseTime)+seconds(smapTime));
clear smapBaseTime clear smapTime

% Replace fill values with NaN
smapS(smapS == -9999) = NaN;
smapGland(smapGland == -9999) = NaN;
smapFland(smapFland == -9999) = NaN;
smapGice(smapGice == -9999) = NaN;
smapSurtep(smapSurtep == -9999) = NaN;

% Remove salinity data near land, near ice, cold SST region using thresholds provided with dataset
%  land_ice_exclusions = 'discard observations if land or seaice fraction exceeds threshold'
%  fland_fraction_threshold = 0.001
%  gland_fraction_threshold = 0.04
%  seaice_fraction_threshold = 0.003
smapS(smapFland>0.001 | smapGland>0.04 | smapGice>0.003 | smapSurtep<278.15)=NaN;
clear smapGland smapFland smapGice smapSurtep 

% Standardize Lat - 90 to 90
%min(smapLat) % -89.8750 0.25 degree spacing
%max(smapLat) % 89.8750

% Standardize Lon 0 to 360
%min(smapLon) % 0.1250 0.25 degree spacing
%max(smapLon) % 359.8750

% Remove data outside month/year range of run
tempIndex = find(smapDateTime < datenum(startDate) | smapDateTime > datenum(endDate));
smapDateTime(tempIndex) = [];
smapS(:,:,tempIndex) = [];
clear tempIndex

% Quick plot to ensure looks right
% [x,y]=meshgrid(smapLat, smapLon);
% figure(1)
% contourf(y, x, smapS(:,:,1),'LevelStep', 0.2,'LineStyle','none');

%% Load Third Dataset - Argo

% Go to location of Argo data
cd(argoDir);

fnArgoMainClimatology = 'RG_ArgoClim_Salinity_2019.nc'; % Covers 2004-2018 must be loaded regardless of times chosen

argoLon=ncread(fnArgoMainClimatology,'LONGITUDE'); % 360x1 degrees east
argoLat = ncread(fnArgoMainClimatology,'LATITUDE'); % 145x1 degrees north
argoPres = ncread(fnArgoMainClimatology,'PRESSURE'); % 58x1 dbar +down
argoMonths = ncread(fnArgoMainClimatology,'TIME'); % 180x1 'months since 2004-01-01 00:00:00' '01-JAN-2004 00:00:00'
argoSmean = ncread(fnArgoMainClimatology,'ARGO_SALINITY_MEAN'); %360x145x58 Lon, Lat, Press PSS78 -999 fill 'ARGO SALINITY MEAN Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY'
argoSanom = ncread(fnArgoMainClimatology,'ARGO_SALINITY_ANOMALY');% 360x145x58x180 Lon, Lat, Press PSS78 -999 fill 'ARGO SALINITY ANOMALY defined by Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY'
%argoBathMask = ncread(fnArgoMainClimatology,'BATHYMETRY_MASK');% 360x145x58 Lon, Lat, Press PSS78 -9
%argoMapMask = ncread(fnArgoMainClimatology,'MAPPING_MASK');% 360x145x58 Lon, Lat, Press PSS78 -9 'MAPPING MASK: pressure limits of mapping can be shallower than 2000dbar in marginal seas '

% If run year is larger than 2018, also load in data from monthly climatology extensions
if endYear > 2018
    % Get list of all unzipped files in directory that are not main climatology
    % One file for each month since end of 2018
    tempFileList = dir(['RG_ArgoClim_2','*', '.nc']);
    tempFileList=char(tempFileList.name);
    
    for fileLoop = 1:size(tempFileList,1) 
        fnTemp = tempFileList(fileLoop,:);
        argoMonths= vertcat(argoMonths, ncread(fnTemp,'TIME')); % 1x1 'months since 2004-01-01 00:00:00' '01-JAN-2004 00:00:00'
        argoSanom = cat(4,argoSanom,ncread(fnTemp,'ARGO_SALINITY_ANOMALY'));% 360x145x58x1 Lon, Lat, Press PSS78 -999 fill 'ARGO SALINITY ANOMALY defined by Jan 2004 - Dec 2018 (15.0 year) RG CLIMATOLOGY'
        clear fnTemp
    end
    clear tempFileList fileLoop
end

% Replace fill values -999 with NaN
argoSmean(argoSmean == -999) = NaN;
argoSanom(argoSanom == -999) = NaN;

% Convert climatology time to MATLAB datetime
% argoMonths is months since 1-Jan-2004
% argoMonths is centered on 15th of month so is given as(0.5, 1.5, 2.5, 3.5 etc.)
% Easiest way to create matlab datenum is to 
% First, add 14 days to 1-Jan-2004 so base is 15-Jan-2004
% Second, round down argoTime to closest integer (0,1,2,3, etc.) for use with calmonths function
argoBaseTime = datenum('01-JAN-2004 00:00:00')+14; 
argoDateTime = datenum(datetime(datevec(argoBaseTime)) + calmonths(floor(argoMonths)));
clear argoBaseTime

% Keep only shallowest pressure level data (2.5db)
tempIndex = find(argoPres==min(argoPres));
argoPres = argoPres(tempIndex);
argoSmean = argoSmean(:,:,tempIndex); %360x145 Lon, Lat
argoSanom = squeeze(argoSanom(:,:,tempIndex,:));% 360x145xN Lon, Lat, Time 
%argoBathMask = argoBathMask(:,:,tempIndex);% 360x145 Lon, Lat
%argoMapMask = argoMapMask(:,:,tempIndex);% 360x145 Lon, Lat 
clear tempIndex

% Calculate SSS from mean and anomalies
sss_argo = argoSmean+argoSanom;

% Standardize Lat - 90 to 90
%min(argoLat) %-64.5 1 degree spacing
%max(argoLat) % 79.5

% Standardize Lon 0 to 360
%min(argoLon) % 20.5 1 degree spacing
%max(argoLon) % 379.5

% Adjust Longitude so 0-360
tempIndex  = find(argoLon>360);
argoLon(tempIndex) = argoLon(tempIndex)-360;
clear tempIndex

[argoLon, tempIndex] = sort(argoLon);
argoSmean = argoSmean(tempIndex,:);
argoSanom = argoSanom(tempIndex,:,:);
%argoBathMask = argoBathMask(tempIndex,:);
%argoMapMask = argoMapMask(tempIndex,:);
sss_argo = sss_argo(tempIndex,:,:);
clear tempIndex

% Quick plot to ensure looks right
%[x,y]=meshgrid(argoLat, argoLon);
%figure(1)
%contourf(y, x, argoSmean(:,:,1))
%figure(2)
%contourf(y, x, sss_argo(:,:,1))

% Remove data outside month/year range of run
tempIndex = find(argoDateTime < datenum(startDate) | argoDateTime > datenum(endDate));
argoDateTime(tempIndex) = [];
argoSanom(:,:,tempIndex) = [];
sss_argo(:,:,tempIndex) = [];

clear fnArgoMainClimatology

%% Put all data on uniform space dimensions
% Argo is coarsest, so use argo lat lon

% Interpolate SMOS and SMAP data onto same spatial resolution as Argo
[X,Y]=meshgrid(argoLat,argoLon);
for monthLoop=1:numMonths
    sss_smos(:,:,monthLoop)= interp2(smosLat,smosLon,smosS(:,:,monthLoop),X,Y);
    sss_smap(:,:,monthLoop)=interp2(smapLat,smapLon,smapS(:,:,monthLoop),X,Y);
end

sss_argo=sss_argo;

% Quick plot to ensure looks right
% [x,y]=meshgrid(argoLat, argoLon);
% figure(1)
% contourf(Y,X,sss_smos(:,:,1),'Levelstep', 0.5)
% caxis([30 39])
% figure(2)
% contourf(Y,X,sss_smap(:,:,1),'Levelstep', 0.5)
% caxis([30 39])
% figure(3)
% contourf(Y,X,sss_argo(:,:,1),'Levelstep', 0.5)
% caxis([30 39])

%% Determine relative weight of each lat/lon box for global averaging
%  Weight is the area of each lat,lon bounded box
%  Area of each lat, lon bounded box is determined with function ctdarea from the Climate Data Toolbox
weightLatLonBox = cdtarea(X,Y, 'km2');
%weightLatLonBox = weightLatLonBox/sum(sum(weightLatLonBox,1));

%% Call triple point collocation code for each time step

% Preallocate RMSD variables for speed
rmsd_smos_temporal = NaN(numMonths,1);
rmsd_smap_temporal = NaN(numMonths,1);
rmsd_argo_temporal = NaN(numMonths,1);

% For each month, call basic_triple_point_collocation.m which will return
% domain average RMSD for each dataset for each month
for monthLoop=1:numMonths
    [loop_rmsd_smos, loop_rmsd_smap, loop_rmsd_argo,]=...
        covariance_triple_point_collocation(sss_smos(:,:,monthLoop),sss_smap(:,:,monthLoop), ...
        sss_argo(:,:,monthLoop), excludeLargeDiff);
    
    rmsd_smos_temporal(monthLoop,1) = loop_rmsd_smos;
    rmsd_smap_temporal(monthLoop,1) = loop_rmsd_smap;
    rmsd_argo_temporal(monthLoop,1) = loop_rmsd_argo;
    
    clear loop_rmsd_smos loop_rmsd_smap loop_rmsd_argo
end
clear monthLoop
cd(outputDir)

%% Plot RMSD for each dataset
figure(1)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [5 5 16 8]);
set(gcf,'PaperPositionMode','auto'); 

plot(argoDateTime,rmsd_smos_temporal, 'DisplayName','RMSD SMOS','color',[0.4660 0.6740 0.1880], 'LineWidth', 2);
hold on
plot(argoDateTime,rmsd_smap_temporal, 'DisplayName','RMSD SMAP','color',[0 0.4470 0.7410], 'LineWidth', 2);
plot(argoDateTime,rmsd_argo_temporal, 'DisplayName','RMSD ARGO','color',[0.9290 0.6940 0.1250], 'LineWidth', 2);
grid on
set(gca, 'Xlim', [datenum(startDate)-1 datenum(endDate)+1])
set(gca, 'FontSize', 14)
xticks(datenum(dateList(3:6:end,:)))
datetick('x', 'mmm-yyyy', 'keepticks', 'keeplimits')
xlabel('Date','FontSize',14)
ylabel('RMSD','FontSize',14)
legend

print(figure(1), '-dtiff', 'temporal_rmsd_smos_smap_rgargo_L3monthly_v1')

%% Call triple point collocation code for each lat, lon pair for all time steps

% Preallocate RMSD variables for speed
rmsd_smos_spatial = NaN(length(argoLon),length(argoLat));
rmsd_smap_spatial  = NaN(length(argoLon),length(argoLat));
rmsd_argo_spatial  = NaN(length(argoLon),length(argoLat));

% For each month, call covariance_triple_point_collocation.m which will return
% temporal average RMSD for each dataset for each domain location
for lonLoop=1:length(argoLon)
    for latLoop = 1:length(argoLat)
    [loop_rmsd_smos, loop_rmsd_smap, loop_rmsd_argo,]=...
        covariance_triple_point_collocation(squeeze(sss_smos(lonLoop,latLoop,:)),squeeze(sss_smap(lonLoop,latLoop,:)), ...
        squeeze(sss_argo(lonLoop,latLoop,:)), excludeLargeDiff);
    
    rmsd_smos_spatial(lonLoop,latLoop) = loop_rmsd_smos;
    rmsd_smap_spatial(lonLoop,latLoop) = loop_rmsd_smap;
    rmsd_argo_spatial(lonLoop,latLoop) = loop_rmsd_argo;
    
    clear loop_rmsd_smos loop_rmsd_smap loop_rmsd_argo
    end
end
clear lonLoop latLoop
cd(outputDir)

%% Plot rmsd for each dataset
% This plotting option relies on the m_map package

figure(2)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [5 5 8 12]);
set(gcf,'PaperPositionMode','auto');
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'color', 'w')

[x,y]=meshgrid(argoLat, argoLon);

% Set subplot common variables
lonLim = [0 360]; % plot longitude limits 
latLim = [-65 80]; % plot latitude limits
fs = 14; % plot font size
cLim = [0 0.5]; % colorbar limits
landColor = [.7 .7 .7];% Land color: dark grey: [.7 .7 .7], black:[0 0 0] 

subplot(3,1,1)
m_proj('mercator','lon',lonLim,'lat',latLim);
m_coast('patch',landColor,'edgecolor','k'); 
m_grid('box', 'on','tickdir','out','linewidth',3, 'fontsize', fs, 'FontWeight','bold');
hold on
m_pcolor(y, x, rmsd_smos_spatial(:,:))
cm =(cmocean('haline'));
colormap(gca, cm);
caxis(cLim)
colorbar
set(gca, 'FontSize', fs)
title('RMSD SMOS')

subplot(3,1,2)
m_proj('mercator','lon',lonLim,'lat',latLim);
m_coast('patch',landColor,'edgecolor','k'); 
m_grid('box', 'on','tickdir','out','linewidth',3, 'fontsize', fs, 'FontWeight','bold');
hold on
m_pcolor(y, x, rmsd_smap_spatial(:,:))
cm =(cmocean('haline'));
colormap(gca, cm);
caxis(cLim)
colorbar
set(gca, 'FontSize', fs)
title('RMSD SMAP')

subplot(3,1,3)
m_proj('mercator','lon',lonLim,'lat',latLim);
m_coast('patch',landColor,'edgecolor','k'); 
m_grid('box', 'on','tickdir','out','linewidth',3, 'fontsize', fs, 'FontWeight','bold');
hold on
m_pcolor(y, x, rmsd_argo_spatial(:,:))
cm =(cmocean('haline'));
colormap(gca, cm);
caxis(cLim)
colorbar
set(gca, 'FontSize', fs)
title('RMSD Argo')

print(figure(2), '-dtiff', 'spatial_rmsd_smos_smap_rgargo_L3monthly_v1')

%% Plot RMSD for each dataset
% This plotting option does not rely on the m_map package or cmocean
figure(3)
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', [5 5 8 12]);
set(gcf,'PaperPositionMode','auto'); 
[x,y]=meshgrid(argoLat, argoLon);

% Set subplot common variables
lonLim = [0 360];
latLim = [-65 80];
fs = 14; 
cLim = [0 0.5];

subplot(3,1,1)
%contourf(y, x, rmsd_smos_spatial(:,:),'Levelstep', 0.05)
t = pcolor(y, x, rmsd_smos_spatial(:,:));
set(t, 'edgecolor','none');
clear t
hold on
caxis(cLim)
colorbar
xlim(lonLim)
ylim(latLim)
set(gca, 'FontSize', fs)
title('RMSD SMOS')

subplot(3,1,2)
%contourf(y, x, rmsd_smap_spatial(:,:),'Levelstep', 0.05)
t = pcolor(y, x, rmsd_smap_spatial(:,:));
set(t, 'edgecolor','none');
clear t
hold on
hold on
caxis(cLim)
colorbar
xlim(lonLim)
ylim(latLim)
set(gca, 'FontSize', fs)
title('RMSD SMAP')

subplot(3,1,3)
%contourf(y, x, rmsd_argo_spatial(:,:),'Levelstep', 0.05)
t = pcolor(y, x, rmsd_argo_spatial(:,:));
set(t, 'edgecolor','none')
clear t
hold on
hold on
caxis(cLim)
colorbar
xlim(lonLim)
ylim(latLim)
set(gca, 'FontSize', fs)
title('RMSD Argo')
print(figure(3), '-dtiff', 'spatial_rmsd_smos_smap_rgargo_L3monthly_noland_v1')