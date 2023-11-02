% This is the script to generate the SMAP SSS validation results

clear

% load colormap data
load tideheight_cmap_v2

% save the original directory
org_dir=pwd;
 
% set the date for analysis
day1=datetime('2015-4-1'); % start date for analysis
day2=datetime('2023-7-31'); % end date for analysis

dnum1=convertTo(day1,'datenum');
dnum2=convertTo(day2,'datenum');
                        
p=1;

% declare the variable size to save computation time
ds=nan(1,1700000);
ancS=ds; ancT=ds; 
ice_frac=ds; insituS=ds; insitulat=ds; insitulon=ds;
landf_asc=ds; landf_des=ds; pcp=ds; smapS_asc=ds; smapS_des=ds; ws=ds; 

nsample=0; 
for day=dnum1:dnum2
    
    % move to the directory where validation results are saved
    cd ../smap_val/v53/

    fname=['svds_smapv53_en422qc_',num2str(day),'.mat'];
    
    if exist(fname,'file')==2
        
        % load the daily validation files and combine them together
        load(fname)
        
        insituS(nsample+1:nsample+length(match_insituS)) = match_insituS;
        smapS_asc(nsample+1:nsample+length(match_smapS50_asc)) = match_smapS50_asc;
        smapS_des(nsample+1:nsample+length(match_smapS50_des)) = match_smapS50_des;
        ancT(nsample+1:nsample+length(match_sst)) = match_sst;
        insitulat(nsample+1:nsample+length(match_insitulat)) = match_insitulat;
        insitulon(nsample+1:nsample+length(match_insitulon)) = deg20(match_insitulon);
        landf_asc(nsample+1:nsample+length(match_landf_asc)) = match_landf_asc;
        landf_des(nsample+1:nsample+length(match_landf_des)) = match_landf_des;
        ice_frac(nsample+1:nsample+length(match_icef)) = match_icef;
        %ancS(nsample+1:nsample+length(match_hycom_S)) = match_hycom_S;
        %ws(nsample+1:nsample+length(match_ws)) = match_ws;
        %pcp(nsample+1:nsample+length(match_precip)) = match_precip;
        
        nsample=nsample+length(match_insituS);
        
        % calculate the daily averaged and standard deviation of ascending
        % and descending salinity differences
        ds_day_asc=match_smapS50_asc-match_insituS;
        ds_day_des=match_smapS50_des-match_insituS;

        med_asc(p)=median(ds_day_asc,'omitnan');
        std_d_asc(p)=std(ds_day_asc,'omitnan');
        med_des(p)=median(ds_day_des,'omitnan');
        std_d_des(p)=std(ds_day_des,'omitnan');

    else
        % set the value as NAN when there's no data
        med_asc(p)=nan;
        std_d_asc(p)=nan;
        med_des(p)=nan;
        std_d_des(p)=nan;
    end

    p=p+1;

    % move back to the original directory
    cd(org_dir)
end

% calculate the salinity difference for the whole time period
ds_asc=smapS_asc-insituS;
ds_des=smapS_des-insituS;

% plot the global maps of dSSS (SMAP - ARGO SSS)

figure
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_scatter(insitulon,insitulat,5,(ds_asc+ds_des)/2,'filled')
m_coast('patch',[.5 .5 .5],'edgecolor','k');
m_grid('linewi',2,'tickdir','out');
clim([-1 1])
colorbar('southoutside')
axis equal tight, colormap(cmap)
title('SMAP SSS V5.3 - Argo SSS, 2015/05-2023/07')
set(gca,'tickdir','out')

print -dpng fig/SMAPV53_ArgoSSS_all


% CALCULATE STD and MEDIAN around each grid cell

lat=-75:80; lon=20:380; 
[clon,clat]=meshgrid(lon,lat);

x=[insitulon(:), insitulat(:)]; % put in situ location into a variable (x)
Mdl=KDTreeSearcher(x); % build a kd-tree model
  
% n3: location of all points within distance 3
% d3: distance between SMAP (x) and in situ 
[n3,d3]=rangesearch(Mdl,[clon(:) clat(:)],0.7); 
        
med_ds_asc=nan(size(clon)); std_ds_asc=nan(size(clon));
med_ds_des=nan(size(clon)); std_ds_des=nan(size(clon));
for n = 1:length(clon(:))
    med_ds_asc(n)=median(ds_asc(n3{n}),'omitnan');
    std_ds_asc(n)=std(ds_asc(n3{n}),'omitnan');
    med_ds_des(n)=median(ds_des(n3{n}),'omitnan');
    std_ds_des(n)=std(ds_des(n3{n}),'omitnan');
        
end

figure
subplot(3,1,1)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_pcolor(clon,clat,med_ds_asc), shading flat
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([-1 1])
colorbar, colormap(cmap)
title('median of dSSS (SMAP V53 ASC - ARGO)')
set(gca,'tickdir','out')

subplot(3,1,2)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_pcolor(clon,clat,med_ds_des), shading flat
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([-1 1])
colorbar, colormap(cmap)
title('median of dSSS (SMAP V53 DESC -ARGO)')
set(gca,'tickdir','out')

subplot(3,1,3)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_pcolor(clon,clat,med_ds_asc-med_ds_des), shading flat
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
    
clim([-1 1])
colorbar, colormap(cmap)
title('median of dSSS (SMAP ASC - DESC)')
set(gca,'tickdir','out')

print -dpng fig/region_dsss_ascdes_v53

figure
subplot(2,1,1)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_pcolor(clon,clat,std_ds_asc), shading flat
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([0 1])
colorbar
cptcmap('WhViBlGrYeOrRe')
title('STD of dSSS (SMAP V53 ASC -ARGO)')
set(gca,'tickdir','out')

subplot(2,1,2)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_pcolor(clon,clat,std_ds_des), shading flat
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([0 1])
colorbar
cptcmap('WhViBlGrYeOrRe')
title('STD of dSSS (SMAP V53 DESC -ARGO)')
set(gca,'tickdir','out')

print -dpng fig/region_dsss_std_v5


%sensitivity tests

% dSSS vs SST

figure
xbin=0:0.25:35;
[mat,xmid,ymid]=twodhist(ancT,ds_asc,xbin,-2:0.1:2);
subplot(2,1,1)
pcolor(xmid,ymid,mat), shading flat        
hold on
cptcmap('WhViBlGrYeOrRe')
title('SST (degree C)'),   ylabel('SMAP V53 ASC - ARGO')
set(gca,'tickdir','out')
ylim([-1.5 1.5])

med_ds_asc=nan(1,length(xmid)); std_ds_asc=nan(1,length(xmid));
for nbin=1:length(xmid)
    fsst=findrange(ancT,xbin(nbin),xbin(nbin+1));
    if isempty(fsst)==0
        med_ds_asc(nbin)=mean(ds_asc(fsst),'omitnan');
        std_ds_asc(nbin)=std(ds_asc(fsst),'omitnan');
    end
end

line([0 35],[0 0],'color','k','linewidth',2), colorbar
plot(xmid,med_ds_asc,'w','linewidth',2)
plot(xmid,std_ds_asc,'r','linewidth',2)

[mat,xmid,ymid]=twodhist(ancT,ds_des,xbin,-2:0.1:2);
subplot(2,1,2)
pcolor(xmid,ymid,mat), shading flat        
hold on
cptcmap('WhViBlGrYeOrRe')
title('SST (degree C)'),   ylabel('SMAP V53 DESC - ARGO')
set(gca,'tickdir','out')
ylim([-1.5 1.5])

med_ds_des=nan(1,length(xmid)); std_ds_des=nan(1,length(xmid));
for nbin=1:length(xmid)
    fsst=findrange(ancT,xbin(nbin),xbin(nbin+1));
    if isempty(fsst)==0
        med_ds_des(nbin)=mean(ds_des(fsst),'omitnan');
        std_ds_des(nbin)=std(ds_des(fsst),'omitnan');
    end
end

line([0 35],[0 0],'color','k','linewidth',2), colorbar
plot(xmid,med_ds_des,'w','linewidth',2)
plot(xmid,std_ds_des,'r','linewidth',2)

print -dpng fig/ds_vs_sst_v53


% dSSS vs landfrac

figure

xbin=0:0.00005:0.101; % open ocean
[mat,xmid,ymid]=twodhist(landf_asc,ds_asc,xbin,-2:0.1:2);
subplot(2,1,1)
pcolor(xmid,ymid,mat), shading flat   
set(gca,'XScale','log')
set(gca,'xtick',[0.0001 0.001 0.005])
xlim([0 0.02]), ylim([-1.5 1.5])
hold on
cptcmap('WhViBlGrYeOrRe'), colorbar
set(gca,'tickdir','out')
title('land fraction asc'),   ylabel('SMAP V53 ASC - ARGO')

med_ds_asc=nan(1,length(xmid)); std_ds_asc=nan(1,length(xmid));
for nbin=1:length(xmid)
    fwd=findrange(landf_asc,xbin(nbin),xbin(nbin+1));
    if isempty(fwd)==0
        med_ds_asc(nbin)=mean(ds_asc(fwd),'omitnan');
        std_ds_asc(nbin)=std(ds_asc(fwd),'omitnan');
    end
end

plot(xmid,zeros(size(xmid)),'k','linewidth',2)
plot(xmid,med_ds_asc,'w','linewidth',2)
plot(xmid,std_ds_asc,'r','linewidth',2)

[mat,xmid,ymid]=twodhist(landf_des,ds_des,xbin,-2:0.1:2);

subplot(2,1,2)
pcolor(xmid,ymid,mat), shading flat   
set(gca,'XScale','log')
set(gca,'xtick',[0.0001 0.001 0.005])
xlim([0 0.02]), ylim([-1.5 1.5])
hold on
cptcmap('WhViBlGrYeOrRe'), colorbar
set(gca,'tickdir','out')
title('land fraction desc'),   ylabel('SMAP V53 DESC - ARGO')

med_ds_des=nan(1,length(xmid)); std_ds_des=nan(1,length(xmid));
for nbin=1:length(xmid)
    fwd=findrange(landf_des,xbin(nbin),xbin(nbin+1));
    if isempty(fwd)==0
        med_ds_des(nbin)=mean(ds_des(fwd),'omitnan');
        std_ds_des(nbin)=std(ds_des(fwd),'omitnan');
    end
end

plot(xmid,zeros(size(xmid)),'k','linewidth',2)
plot(xmid,med_ds_des,'w','linewidth',2)
plot(xmid,std_ds_des,'r','linewidth',2)

print -dpng fig/ds_vs_landf_v53


% plot the global sea surface salinity maps

figure
subplot(2,1,1)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_scatter(insitulon,insitulat,3,smapS_asc,'filled')
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([32 38])
colorbar
title('SMAP V5.3 SSS ASC')

subplot(2,1,2)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_scatter(insitulon,insitulat,3,smapS_des,'filled')
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.5 .5 .5],'edgecolor','k');
clim([32 38])
colorbar
title('SMAP V5.3 SSS DESC')

print -dpng fig/svds_global_smap_v53_sss


% plot the global maps of the SMAP-in situ SSS

figure
subplot(2,1,1)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_scatter(insitulon,insitulat,3,ds_asc,'filled')
m_coast('patch',[.5 .5 .5],'edgecolor','k');
m_grid('linewi',2,'tickdir','out');
clim([-2 2])
colorbar, axis equal tight, colormap(cmap)
title('SMAP V5.3 SSS asc - Argo SSS')

subplot(2,1,2)
m_proj('miller','lat',[-80 80],'lon',[20 380]);
m_scatter(insitulon,insitulat,3,ds_des,'filled')
m_coast('patch',[.5 .5 .5],'edgecolor','k');
m_grid('linewi',2,'tickdir','out');
clim([-2 2])
colorbar, axis equal tight, colormap(cmap)
title('SMAP V5.3 SSS des - Argo SSS')

print -dpng fig/svds_global_smap_v53


% scattering plot (SMAP SSS vs ARGO SSS)

fnan=find(~isnan(smapS_asc) & insituS<50);
[r,p]=corrcoef(smapS_asc(fnan),insituS(fnan));

figure
subplot(2,1,1)
[mat,xmid,ymid]=twodhist(smapS_asc,insituS,30:0.05:40,30:0.05:40);
pcolor(xmid,ymid,mat), shading flat; outticks     
hold on
plot(25:46,25:46,'k'), colormap(jet), colorbar
axis([32 38 32 38]), axis('square'), clim([0 2000])
t=text(.05,.45,['cor = ',num2str(r(1,2),'%5.2f')],'verticalalignment','top', ...
    'units','normalized','fontsize',14);
t(1).Color='w';
xlabel('SMAP SSS V5.3 ASC'),   ylabel('Argo SSS')
fontsize 14 14 14 14

subplot(2,1,2)
[mat,xmid,ymid]=twodhist(smapS_des,insituS,30:0.05:40,30:0.05:40);
pcolor(xmid,ymid,mat), shading flat; outticks     
hold on
plot(25:46,25:46,'k'), colormap(jet), colorbar
axis([32 38 32 38]), axis('square'), clim([0 2000])
t=text(.05,.45,['cor = ',num2str(r(1,2),'%5.2f')],'verticalalignment','top', ...
    'units','normalized','fontsize',14);
t(1).Color='w';
xlabel('SMAP SSS V5.3 DESC'),   ylabel('Argo SSS')
fontsize 14 14 14 14

print -dpng fig/scat_svds_ascdes_v53


% time series of dSSS (SMAP - ARGO SSS)
% calculate the daily average and plot the median of daily data

dnum=day1:day2;

figure
subplot(3,1,1)          
plot(dnum,med_asc,'linewidth',2.5)
xtick(dnum(1):366:dnum(end))
datetick('x','mmmyy','keepticks')
hold on
line([dnum(1) dnum(end)],[0 0],'color','k','linewidth',2)
grid on
xlim([dnum(1) dnum(end)]), ylim([-0.4 0.4])
title('median of SMAP V5.3 ASC - ARGO SSS')

subplot(3,1,2)          
plot(dnum,med_des,'linewidth',2.5)
xtick(dnum(1):366:dnum(end))
datetick('x','mmmyy','keepticks')
hold on
line([dnum(1) dnum(end)],[0 0],'color','k','linewidth',2)
grid on
xlim([dnum(1) dnum(end)]), ylim([-0.4 0.4])
title('median of SMAP V5.3 DESC - ARGO SSS')

subplot(3,1,3)          
plot(dnum,med_asc-med_des,'linewidth',2.5)
xtick(dnum(1):366:dnum(end))
datetick('x','mmmyy','keepticks')
hold on
line([dnum(1) dnum(end)],[0 0],'color','k','linewidth',2)
grid on
xlim([dnum(1) dnum(end)]), ylim([-0.4 0.4])
title('median of SMAP V5.3 ASC - DESC')

print -dpng fig/ts_dsss_ascdes_v53

% plot the histograms of dSSS (SMAP - ARGO SSS)
figure
subplot(2,1,1)

x = -3:0.1:3;
histogram(ds_asc,x)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','k'), outticks
txtstd = ['STD = ',num2str(std(ds_asc,'omitnan'),'%4.2f')];
txtmedian = ['median = ',num2str(median(ds_asc,'omitnan'),'%4.2f')];
txtmean = ['mean = ',num2str(mean(ds_asc,'omitnan'),'%4.2f')];
    
text(.05,.45,[txtmedian,' \newline',txtmean,' \newline',txtstd] ...
    ,'verticalalignment','top','units','normalized','fontsize',12)
xlim([-2 2])
xlabel('SMAP V5.3 ASC - ARGO')
ylabel('count')

subplot(2,1,2)

x = -3:0.1:3;
histogram(ds_des,x)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','k'), outticks
txtstd = ['STD = ',num2str(std(ds_des,'omitnan'),'%4.2f')];
txtmedian = ['median = ',num2str(median(ds_des,'omitnan'),'%4.2f')];
txtmean = ['mean = ',num2str(mean(ds_des,'omitnan'),'%4.2f')];
    
text(.05,.45,[txtmedian,' \newline',txtmean,' \newline',txtstd] ...
    ,'verticalalignment','top','units','normalized','fontsize',12)
xlim([-2 2])
xlabel('SMAP V5.3 DESC - ARGO')
ylabel('count')

print -dpng fig/svds_hist_smap_ascdes_v53

%%%% calculate dSSS at different latitude

a=1;
lat_bin=-70:3:70;
for nlat=lat_bin
    flat=findrange(insitulat,nlat-1.5,nlat+1.5);
    ds_med_asc(a)=median(ds_asc(flat),'omitnan');
    ds_std_asc(a)=std(ds_asc(flat),'omitnan');
    ds_med_des(a)=median(ds_des(flat),'omitnan');
    ds_std_des(a)=std(ds_des(flat),'omitnan');
    a=a+1;
end

figure
subplot(2,2,3)
line(ds_med_asc,lat_bin,'color','b','linewidth',2)

ax1=gca;
ax1.XColor='b';
ax1.YColor='k';

ax1_pos=ax1.Position;
ax2=axes('Position',ax1_pos, ...
    'XAxisLocation','top', ...
    'YAxisLocation','right', ...
    'Color','none');

line(ds_std_asc,lat_bin,'Parent',ax2,'Color','r','linewidth',2)
ax2.XColor='r';
ax2.YColor='k';
axis(ax1,[-1 1 -70 70])
axis(ax2,[0 2 -70 70])
ax1.XLabel.String='median';
ax2.XLabel.String='STD';   
ax1.YTick=-70:10:70;  
ax2.YTick=-70:10:70;
ax1.YTickLabel={'','60\circS','','40\circS','','20\circS','','EQ','', ...
    '20\circN','','40\circN','','60\circN',''};
ax2.YTickLabel={'','60\circS','','40\circS','','20\circS','','EQ','', ...
    '20\circN','','40\circN','','60\circN',''};
hold on
line([0 0],[-70 70],'Parent',ax1,'Color','k')
grid on
title('SMAP V5.3 asc')

subplot(2,2,4)
line(ds_med_des,lat_bin,'color','b','linewidth',2)

ax1=gca;
ax1.XColor='b';
ax1.YColor='k';

ax1_pos=ax1.Position;
ax2=axes('Position',ax1_pos, ...
    'XAxisLocation','top', ...
    'YAxisLocation','right', ...
    'Color','none');

line(ds_std_des,lat_bin,'Parent',ax2,'Color','r','linewidth',2)
ax2.XColor='r';
ax2.YColor='k';
axis(ax1,[-1 1 -70 70])
axis(ax2,[0 2 -70 70])
ax1.XLabel.String='median';
ax2.XLabel.String='STD';   
ax1.YTick=-70:10:70;  
ax2.YTick=-70:10:70;
ax1.YTickLabel={'','60\circS','','40\circS','','20\circS','','EQ','', ...
    '20\circN','','40\circN','','60\circN',''};
ax2.YTickLabel={'','60\circS','','40\circS','','20\circS','','EQ','', ...
    '20\circN','','40\circN','','60\circN',''};
hold on
line([0 0],[-70 70],'Parent',ax1,'Color','k')
grid on
title('SMAP V5.3 des')

print -dpng fig/dsss_latband_smap_v53

