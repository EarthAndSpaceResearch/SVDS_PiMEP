function [rmse_smos,rmse_smap,rmse_argo,std_smap_argo]=...
    basic_triple_point_collocation(start_date,months)

% Triple Point Collocation Analysis for satellite (SMAP,SMOS) Level 3 and 
% gridded Argo data
% 1. SMOS data is from CATDS CEC-Locean L3 Debiased V4, temporal sampling
% is 4 days, spatial resolution of 25 km with 9 day FWHM gaussian smoothing
% 2. SMAP data is from RSS V4.0 Level 3 monthly (RF) rain-filtered averages
% 70-km smoothed product.
% 3. Argo data is using the ancillary data from SMAP data generated from
% monthly 1-degree gridded interpolated Argo SSS field provided by Scripps.
% 
% Written by Hsun-Ying Kao (hkao@esr.org)

%start_date=datevec('4-1-2015');
%months=45;

year_start=start_date(1);
month_start=start_date(2);

org_dir=pwd; % save the original directory

yr0=year_start; mon0=month_start; p=1;

for mon=1:months
    
    % SET THE DIRECTORY WHERE SMOS L3 DATA ARE LOCATED
    cd ('../../data_aq/SMOS/debiasedSSS_09days_v4/') 
    fn_smos=dir(['SMOS_L3_DEBIAS_LOCEAN_AD_',num2str(yr0),num2str(mon0,'%02i'),'*']);
    fn_smos=char(fn_smos.name);
    
    for n=1:size(fn_smos)
        ncload(fn_smos(n,:),'SSS')
        smos_s(n,:,:)=SSS;
    end

    ncload(fn_smos(1,:),'lat','lon')
    smos_lat=lat; [smos_lon,I]=sort(deg360(lon));

    % calculate the monthly SMOS SSS
    sss_smos=smean(smos_s);

    clear smos_s

    cd(org_dir) % move to original directory

    % load monthly SMAP data
    
    % SET THE DIRECTORY WHERE SMAP L3 DATA ARE LOCATED
    cd ../data/V4/L3/monthly_RF/
    fn_smap=['RSS_smap_SSS_L3_monthly_RF_',num2str(yr0),'_', ...
        num2str(mon0,'%02i'),'_FNL_v04.0.nc'];
    ncload([num2str(yr0),'/',fn_smap],'sss_smap','sss_argo','lat','lon', ...
        'gland','gice','surtep')
    
    % set flags to remove data near land, near ice, cold SST region
    sss_smap(gland>0.001 | gice>0.001 | surtep<278.15)=nan;
    sss_smap(sss_smap<0) = nan;
    sss_argo(gland>0.001 | gice>0.001 | surtep<278.15)=nan;
    sss_argo(sss_argo<0) = nan;
    
    % interpolate SMOS data into same spatial resolution as SMAP
    [X,Y]=meshgrid(lon,lat);
    smos_interp=interp2(smos_lon,smos_lat,sss_smos(:,I),X,Y);

    sss_smap(isnan(sss_argo))=nan;
    smos_interp(isnan(sss_smap))=nan;
    
        
    %'SMOS - Argo'
    ds_smos_argo=smos_interp-sss_argo;
    bias_smos_argo=median(ds_smos_argo(:),'omitnan');
    std_smos_argo=std(ds_smos_argo(:),'omitnan');
    ms_smos_argo=bias_smos_argo.^2+std_smos_argo.^2;
    
    %'SMAP - Argo'
    ds_smap_argo=sss_smap-sss_argo;
    bias_smap_argo=median(ds_smap_argo(:),'omitnan');
    std_smap_argo=std(ds_smap_argo(:),'omitnan');
    ms_smap_argo=bias_smap_argo.^2+std_smap_argo.^2;
    
    %'SMOS - SMAP'
    ds_smos_smap=smos_interp-sss_smap;
    bias_smos_smap=median(ds_smos_smap(:),'omitnan');
    std_smos_smap=std(ds_smos_smap(:),'omitnan');
    ms_smos_smap=bias_smos_smap.^2+std_smos_smap.^2;
    
    % Root Mean Square Error calculated from MS
    rmse_smos(mon)=sqrt((ms_smos_argo+ms_smos_smap-ms_smap_argo)/2);
    rmse_smap(mon)=sqrt((ms_smap_argo+ms_smos_smap-ms_smos_argo)/2);
    rmse_argo(mon)=sqrt((ms_smos_argo+ms_smap_argo-ms_smos_smap)/2);
    
    month(p)=mon0; year(p)=yr0;
    
    p=p+1;
    mon0=mon0+1;
    if mon0>12, mon0=1; yr0=yr0+1; end
        
    cd (org_dir)
end

figure
plot(1:mon,rmse_smos,1:mon,rmse_smap,1:mon,rmse_argo,'linewidth',2)
title('RMSE')
legend('SMOS,','SMAP','Argo')


