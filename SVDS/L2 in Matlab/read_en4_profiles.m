% This is the script to read EN4 individual Argo float data and save  it
% locally with only the variables needed with quality control

e = referenceEllipsoid('GRS80','km'); % too slow
s=referenceSphere('Earth');
s.LengthUnit = 'kilometer';
currdir = pwd;

for theyear=2016:2022  
    for themonth=1:12
    
        yearname=sprintf('%04.4i',theyear);
        monname=sprintf('%02.2i',themonth);
        fname=['/salsci/data/ARGO_EN4/EN.4.2.2.f.profiles.g10.' yearname monname '.nc'];
        insitu_lat=nc_varget(fname,'LATITUDE');
        insitu_lon=nc_varget(fname,'LONGITUDE');
        insitu_juld=nc_varget(fname,'JULD');
        insitu_pnum=nc_varget(fname,'PLATFORM_NUMBER');

        % units = "days since 1950-01-01 00:00:00 utc" 
        % conventions = "relative julian days with decimal part (as parts of day)" 
        dayoffset=datetime('1950-1-1'); 
        dateoffset=convertTo(dayoffset,'datenum');
        insitu_datamode=nc_varget(fname,'DATA_MODE');
        insitu_depth=nc_varget(fname,'DEPH_CORRECTED');
        insitu_psal=nc_varget(fname,'PSAL_CORRECTED');
        insitu_psal_qc=nc_varget(fname,'PSAL_CORRECTED_QC'); 
        insitu_psal_surface=insitu_psal(:,1);
        insitu_psal_2nd=insitu_psal(:,2);
        insitu_profile_psal_qc=nc_varget(fname,'PROFILE_PSAL_QC');
        insitu_qc_flags_prof=nc_varget(fname,'QC_FLAGS_PROFILES');
        insitu_qc_flags_levels=nc_varget(fname,'QC_FLAGS_LEVELS');
        insitu_temp=nc_varget(fname,'TEMP');
        insitu_depth_surface=insitu_depth(:,1);
            
        % quality control
        good_sss=insitu_depth(:,1)<6;
        good_qc=str2num(insitu_psal_qc(:,1))<2; %#ok<*ST2NM>
        good_profile=str2num(insitu_profile_psal_qc(:,1))<2;
        good_qc_prof=insitu_qc_flags_prof<1;
        good_qc_level=insitu_qc_flags_levels(:,1)<1;
        good_data=good_sss & good_qc & good_profile & good_qc_level & ...
            good_qc_prof;
           
        in_situ_psal_surface=insitu_psal_surface(good_data);
        in_situ_psal_2nd=insitu_psal_2nd(good_data);
        in_situ_temp=insitu_temp(good_data);
        in_situ_latitude=insitu_lat(good_data);
        in_situ_longitude=insitu_lon(good_data);
        in_situ_juld=insitu_juld(good_data);
        in_situ_pnum=insitu_pnum(good_data,:);
        in_situ_depth_surface=insitu_depth_surface(good_data);
        in_situ_datamodes=insitu_datamode(good_data);
        in_situ_dates=in_situ_juld+dateoffset;
            
        clear insitu_juld good* insitu_depth insitu_psal ...
            insitu_psal_qc insitu_temp insitu_psal_2nd ...
            insitu_depth_surface insitu_psal_surface insitu_pnum ...
            insitu_qc* insitu_profile_psal_qc
            
        for theday=floor(min(in_situ_dates)):floor(max(in_situ_dates))
            dd=datevec(theday);
            dayname=sprintf('%02.2i',dd(3));
            fday=findrange(in_situ_dates,theday,theday+1);
                            
            insitu_sss=in_situ_psal_surface(fday);
            insitu_sl2=in_situ_psal_2nd(fday);
            insitu_sst=in_situ_temp(fday);
            insitu_lat=in_situ_latitude(fday);
            insitu_lon=in_situ_longitude(fday);
            insitu_obsdepth=in_situ_depth_surface(fday);
            insitu_datamode=in_situ_datamodes(fday);
            insitu_date=in_situ_dates(fday);
            insitu_platform_num=in_situ_pnum(fday,:);
                
            fnameout=['../data/argo_en4/EN422f_sss_moreqc_' ...
                yearname monname dayname '.mat'];

            save(fnameout,'insitu*','-v7.3','-nocompression');       
        end
    end
end


%tic;[obsdist,~] = distance(10, 125,latout,lonout,s);toc % 0.3 seconds with e, 0.1 with s
%tic;[len] = m_idist(125,10,lon(:),lat(:));toc % 1.8 seconds
