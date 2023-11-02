% This is the script to do matchup with RSS SMAP SSS data and individual
% EN4 Argo floats

clear

% set the Earth radius for calculating the distance for match-up
ref_earth=referenceSphere('Earth');
ref_earth.LengthUnit = 'kilometer';

% save the original directory
org_dir=pwd;
 
% set the date for analysis
day1=datetime('2015-4-1'); % start date for analysis
day2=datetime('2023-7-31'); % end date for analysis

dnum1=convertTo(day1,'datenum');
dnum2=convertTo(day2,'datenum');

for day=dnum1:dnum2

    dd2=datevec(day);
    yr=num2str(dd2(1)); mn=num2str(dd2(2),'%02i'); d3=num2str(dd2(3),'%02i');
    
    % read ARGO EN4 .mat files 
    fn=['../data/argo_en4/EN422f_sss_moreqc_' yr,mn,d3 '.mat'];
    
    if isempty(fn)==0 
        
        load(fn)
        ln=0;
        % Calculate how many SMAP data will be loaded.
        for sat_day=day-4:day+4
        
            dd=datevec(sat_day);
    
            %   move to the directory where the SMAP files are
            dir_name=['/salsci/data/SMAP_RSS/V05.3/FINAL/L2C/', ...
                num2str(dd(1)),'/',num2str(dd(2),'%02i')];
            cd(dir_name)
    
            doy=yearday(dd(1),dd(2),dd(3));
            fname=dir(['RSS_SMAP_SSS_L2C_r*',num2str(doy,'%03i'), ...
                '_FNL_V05.3.nc']);
       
            if isempty(fname)==0 % if the file exists
                ln=ln+length(fname);
            end
        
        end
    
    csmap_lat_asc=cell(1,ln);
    csmap_lat_des=cell(1,ln);
    csmap_lon_asc=cell(1,ln);
    csmap_lon_des=cell(1,ln);
    csmap_sss_asc=cell(1,ln); 
    csmap_sss_des=cell(1,ln);       
    csmap_time_asc=cell(1,ln);       
    csmap_time_des=cell(1,ln);  
    ciqc_flag_asc=cell(1,ln);   % iqc_flag
    ciqc_flag_des=cell(1,ln);   
    clandf_asc=cell(1,ln);      
    clandf_des=cell(1,ln);      
    csmap_sst=cell(1,ln);               
    chycom_sss=cell(1,ln);            
    cwindspd=cell(1,ln);
    cprecip=cell(1,ln);
    cicef=cell(1,ln);       % gice_est        
   
    % load SMAP data from 4 days before to 4 days after ARGO observation
    % date
    m=1;
    for sat_day=day-4:day+4
    
    dd=datevec(sat_day);
    
    % LOAD SMAP DATA
    %   move to the directory where the SMAP files are
    dir_name=['/salsci/data/SMAP_RSS/V05.3/FINAL/L2C/', ...
        num2str(dd(1)),'/',num2str(dd(2),'%02i')];
    cd(dir_name)
    
    doy=yearday(dd(1),dd(2),dd(3));
    fname=dir(['RSS_SMAP_SSS_L2C_r*',num2str(doy,'%03i'), ...
        '_FNL_V05.3.nc']);
   
    if isempty(fname)==0 % if the file exists
        lname=length(fname);
        a=m; b=m+lname-1;    
       
        for fn_num=a:b
            
            cd(dir_name) 
            filename_l2c=fname(fn_num-m+1).name;
            
            ncid=netcdf.open(filename_l2c);
            varid=netcdf.inqVarID(ncid,'sss_smap');
            sss_smap= netcdf.getVar(ncid,varid);              
            varid=netcdf.inqVarID(ncid,'gland');
            gland= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'gice_est');
            gice= netcdf.getVar(ncid,varid);       
            varid=netcdf.inqVarID(ncid,'iqc_flag');
            iqc_flag= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'cellon');
            cellon= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'cellat');
            cellat= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'surtep');
            surtep= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'time');
            time= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'winspd');
            winspd= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'sss_ref');
            sss_ref= netcdf.getVar(ncid,varid);  
            varid=netcdf.inqVarID(ncid,'rain');
            rain=netcdf.getVar(ncid,varid); 
            
            netcdf.close(ncid);
            
            cellat(cellat==-9999)=nan;
            cellon(cellon==-9999)=nan;
            cellon(cellon>180)=cellon(cellon>180)-360;
 
            % variables doesn't change with looking angle (winspd)
            surtep(surtep==-9999)=nan;   
            surtep=surtep-273.15;                 
            rain(rain==-9999)=nan;            
            sss_ref(sss_ref==-9999)=nan;                
            gice(gice==-9999)=nan;

            sss_smap(sss_smap==-9999)=nan;

            %%%%% use Q/C flags
            cd(org_dir)
            f1=read_smap_flags(iqc_flag,[1:8 11 17]);
            flag1=find(f1==1);

            sss_smap(flag1)=nan;

            cellon(isnan(sss_smap))=nan; cellat(isnan(sss_smap))=nan;

            % change the SMAP time to MATLAB time
            
            D=datetime(2000,1,1);
            dnum2000=convertTo(D,'datenum');
            smap_tt=time/86400+dnum2000;
            smap_tt(smap_tt<dnum2000)=nan;

            %separate LAT,LON,SSS,TIME,LANDF,FLAG into ascending and descending
        
            cellat_asc=squeeze(cellat(1,:,:));
            cellat_des=squeeze(cellat(2,:,:));    
            cellon_asc=squeeze(cellon(1,:,:));
            cellon_des=squeeze(cellon(2,:,:));
            sss_smap_asc=squeeze(sss_smap(1,:,:));
            sss_smap_des=squeeze(sss_smap(2,:,:));
            time_asc=squeeze(smap_tt(1,:,:));
            time_des=squeeze(smap_tt(2,:,:));
            landf_asc=squeeze(gland(1,:,:));
            landf_des=squeeze(gland(2,:,:));
            flag_asc=squeeze(iqc_flag(1,:,:));
            flag_des=squeeze(iqc_flag(2,:,:));


            fnnan=find(~isnan(sss_smap_asc) & ~isnan(sss_smap_des));
   
            csmap_lat_asc{fn_num}=cellat_asc(fnnan)';
            csmap_lat_des{fn_num}=cellat_des(fnnan)';
            csmap_lon_asc{fn_num}=cellon_asc(fnnan)';
            csmap_lon_des{fn_num}=cellon_des(fnnan)';
            csmap_sss_asc{fn_num}=sss_smap_asc(fnnan)';
            csmap_sss_des{fn_num}=sss_smap_des(fnnan)';
            csmap_time_asc{fn_num}=time_asc(fnnan)';
            csmap_time_des{fn_num}=time_des(fnnan)';
            ciqc_flag_asc{fn_num}=flag_asc(fnnan)';
            ciqc_flag_des{fn_num}=flag_des(fnnan)';
            clandf_asc{fn_num}=landf_asc(fnnan)';
            clandf_des{fn_num}=landf_des(fnnan)';
            csmap_sst{fn_num}=surtep(fnnan)';
            chycom_sss{fn_num}=sss_ref(fnnan)';
            cwindspd{fn_num}=winspd(fnnan)';
            cprecip{fn_num}=rain(fnnan)';
            cicef{fn_num}=gice(fnnan)';
            
    
        end
    
        m=m+lname;
    end
    end

    % transform data from cell to mat data format

    smap_lat_asc=cell2mat(csmap_lat_asc);   
    smap_lat_des=cell2mat(csmap_lat_des);          
    smap_lon_asc=cell2mat(csmap_lon_asc);         
    smap_lon_des=cell2mat(csmap_lon_des); 
    smap_sss_asc=cell2mat(csmap_sss_asc);  
    smap_sss_des=cell2mat(csmap_sss_des);  
    smap_time_asc=cell2mat(csmap_time_asc);
    smap_time_des=cell2mat(csmap_time_des); 
    landf_asc=cell2mat(clandf_asc);         
    landf_des=cell2mat(clandf_des);          
    iqc_flag_asc=cell2mat(ciqc_flag_asc);        
    iqc_flag_des=cell2mat(ciqc_flag_des);

    smap_sst=cell2mat(csmap_sst);              
    hycom_sss=cell2mat(chycom_sss);            
    windspd=cell2mat(cwindspd);            
    precip=cell2mat(cprecip);           
    icef=cell2mat(cicef);
    
    p=1;
    
  
    clear match_insitulon match_insitulat match_insituS match_insituT 
    clear match_smapcpa_asc match_smapcpa_des match_smapS50_asc match_smapS50_des
    clear match_insitutime filename
    clear match_smaplat_asc match_smaplat_des match_smaplon_asc 
    clear match_smaplon_des match_hycom_S match_sst
    clear match_ws match_precip match_landf_asc match_landf_des 
    clear match_icef match_insituP match_insitudepth     
    clear match_insitu_platformnum match_flag_asc match_flag_des 
    
    length(insitu_sss)
 
    smap_x_asc=[smap_lon_asc(:) smap_lat_asc(:)]; % put SMAP orbit location into a variable (x)
    smap_x_des=[smap_lon_des(:) smap_lat_des(:)]; 
    smap_Mdl_asc=KDTreeSearcher(smap_x_asc); % build a kd-tree model
    smap_Mdl_des=KDTreeSearcher(smap_x_des); 

    % n3: location of all points within distance 3
    % d3: distance between SMAP (x) and in situ 
    [smap_n3_asc,smap_d3_asc]=rangesearch(smap_Mdl_asc,[insitu_lon insitu_lat],3); 
    [smap_n3_des,smap_d3_des]=rangesearch(smap_Mdl_des,[insitu_lon insitu_lat],3); 
    
    for n = 1:length(insitu_sss)
            
        if size(smap_n3_asc{n},2)>0  % check if SMAP data nearby found
            
            nlat=insitu_lat(n);
            nlon=insitu_lon(n);

            % calculate the actual distance on a sphere (EARTH)
            dist3_smap=distance(nlat,nlon, ...
                smap_lat_asc(smap_n3_asc{n}),smap_lon_asc(smap_n3_asc{n}),ref_earth);

            ns = insitu_sss(n); ntime = insitu_date(n); nt = insitu_sst(n); 
            nd = insitu_obsdepth(n);
           
            % keep only data within time-space window (50 km & +-3.5 days)
            f50_smap=find(dist3_smap < 50 & ...
                abs(ntime-smap_time_asc(smap_n3_asc{n}))<3.5); 
            find_near_smap=smap_n3_asc{n};                        
            
            % find the SMAP footprints that are close to the in situ
            % data using 50 km criteria
            
            [shortest_dist_smap,Ismap] = min(dist3_smap(f50_smap));
            
            if shortest_dist_smap<50 
                
            match_insitulat(p) = nlat; match_insitulon(p) = nlon;
            match_insituS(p) = ns; match_insituT(p) = nt; 
            match_insitutime(p) = ntime; match_insitudepth(p) = nd;
            match_smaplat_asc(p) = smap_lat_asc(find_near_smap(f50_smap(Ismap))); 
            match_smaplat_des(p) = smap_lat_des(find_near_smap(f50_smap(Ismap))); 
            match_smaplon_asc(p) = smap_lon_asc(find_near_smap(f50_smap(Ismap)));
            match_smaplon_des(p) = smap_lon_des(find_near_smap(f50_smap(Ismap)));
            match_smapS50_asc(p) = ...
                median(smap_sss_asc(find_near_smap(f50_smap)),'omitnan');
            match_smapS50_des(p) = ...
                median(smap_sss_des(find_near_smap(f50_smap)),'omitnan');
            match_smapcpa_asc(p) = smap_sss_asc(find_near_smap(Ismap));
            match_smapcpa_des(p) = smap_sss_des(find_near_smap(Ismap));
            match_hycom_S(p) = ...
                median(hycom_sss(find_near_smap(f50_smap)),'omitnan');
            match_sst(p) = smap_sst(find_near_smap(f50_smap(Ismap)));
            match_ws(p) = windspd(find_near_smap(f50_smap(Ismap)));
            match_precip(p) = precip(find_near_smap(f50_smap(Ismap)));
            match_flag_asc(p) = iqc_flag_asc(find_near_smap(f50_smap(Ismap)));
            match_flag_des(p) = iqc_flag_des(find_near_smap(f50_smap(Ismap)));
            match_landf_asc(p) = landf_asc(find_near_smap(f50_smap(Ismap)));
            match_landf_des(p) = landf_des(find_near_smap(f50_smap(Ismap)));
            match_icef(p) = icef(find_near_smap(f50_smap(Ismap)));
            match_insitu_platformnum(p,:) = insitu_platform_num(n,:);
  

            p = p + 1;
            
        
            end
        end     
    end 

    if p==1
        disp('No Match')
    else

        ds_smapcpa_asc_insitu = match_smapcpa_asc - match_insituS;
        ds_smapcpa_des_insitu = match_smapcpa_des - match_insituS;
        ds_smapcpa_asc_insitu(abs(ds_smapcpa_asc_insitu)>10)=nan;
        ds_smapcpa_des_insitu(abs(ds_smapcpa_des_insitu)>10)=nan;
        std_ds_smapcpa_asc_insitu = std(ds_smapcpa_asc_insitu,'omitnan');
        std_ds_smapcpa_des_insitu = std(ds_smapcpa_des_insitu,'omitnan');
        med_ds_smapcpa_asc_insitu = median(ds_smapcpa_asc_insitu,'omitnan');
        med_ds_smapcpa_des_insitu = median(ds_smapcpa_des_insitu,'omitnan');
        
        ds_smapS50_asc_insitu = match_smapS50_asc - match_insituS;
        ds_smapS50_des_insitu = match_smapS50_des - match_insituS;
        ds_smapS50_asc_insitu(abs(ds_smapS50_asc_insitu)>10)=nan;
        ds_smapS50_des_insitu(abs(ds_smapS50_des_insitu)>10)=nan;
        std_ds_smapS50_asc_insitu = std(ds_smapS50_asc_insitu,'omitnan');
        std_ds_smapS50_des_insitu = std(ds_smapS50_des_insitu,'omitnan');
        med_ds_smapS50_asc_insitu = median(ds_smapS50_asc_insitu,'omitnan');
        med_ds_smapS50_des_insitu = median(ds_smapS50_des_insitu,'omitnan');
        
        
        cd(org_dir)        

        s_fname=['../smap_val/v53/svds_smapv53_en422qc_',num2str(day)];
        save(s_fname,'match_insitulat','match_insitulon','match_insitudepth', ...
            'match_insituS','match_insituT','match_insitutime', ...
            'match_smaplat_asc','match_smaplon_asc','match_smapS50_asc', ...
            'match_smaplat_des','match_smaplon_des','match_smapS50_des', ...
            'match_hycom_S','match_sst','match_ws','match_precip', ...
            'match_landf_asc','match_landf_des','match_icef', ...
            'match_smapcpa_asc','match_smapcpa_des','match_insitu_platformnum', ...
            'std_ds_smapcpa_asc_insitu','med_ds_smapcpa_asc_insitu', ...
            'std_ds_smapcpa_des_insitu','med_ds_smapcpa_des_insitu', ...
            'std_ds_smapS50_asc_insitu','med_ds_smapS50_asc_insitu', ...
            'std_ds_smapS50_des_insitu','med_ds_smapS50_des_insitu', ...
            'match_flag_asc','match_flag_des')
       
    end
    
        cd(org_dir)
   
    end
end


    
