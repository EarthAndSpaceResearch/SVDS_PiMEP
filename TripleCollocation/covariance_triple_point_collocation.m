function [rmsd_d1,rmsd_d2,rmsd_d3]=...
    covariance_triple_point_collocation(data1, data2, data3, flag)

% Basic Triple Point Collocation Analysis for any combination of datasets
% Input datasets must be 2D with dimensions 1xN, Mx1, or MXN 
% All datasets must have uniform dimensions 
% Flag: 1 = exclude large dataset differences abs(>5PSU) from analysis, 0=include all dataset differences regardless of size
%
% J. Anderson (janderson@esr.org) based on Stoffelen 1998 & Gruber et al 2016 Covariance Notation 

%% Error check input data for uniform dimensions and valid flag value
narginchk(4,4); %Check that 3 datasets were provided

if nargin == 4
    if (sum(size(data1) == size(data2)) + sum(size(data1) == size(data3))) ~=4
        error('Input datasets must have uniform dimensions')
    end
    if ismember(flag, [0 1]) == 0
        error('Flag input must be 1 or 0')
    end
end

%% Prepare data for triple point collocation

% If flag = 1, remove data where absolute value of differences > 5 PSU

if flag == 1
    %'Data1 - Data3'
    ds_d1_d3=data1-data3; % Calculate the difference between Data1 & Data3
    tempLarge = find(abs(ds_d1_d3)>5);
    data1(tempLarge)=NaN;
    data3(tempLarge)=NaN;
    clear tempLarge ds_d1_d3

    %'Data2 - Data3'
    ds_d2_d3=data2-data3; % Calculate the difference between Data2 & Data3
    tempLarge = find(abs(ds_d2_d3)>5);
    data2(tempLarge)=NaN;
    data3(tempLarge)=NaN;
    clear tempLarge ds_d2_d3

    %'Data1 - Data2'
    ds_d1_d2=data1-data2; % Calculate the difference between Data1 & Data2
    tempLarge = find(abs(ds_d1_d2)>5);
    data1(tempLarge)=NaN;
    data2(tempLarge)=NaN;
    clear tempLarge ds_d1_d2
end

% Triple point collocation cannot be calculated if one dataset has a
% NaN. Make location of NaNs uniform for all datasets. 
data2(isnan(data1))=nan;
data3(isnan(data1))=nan;

data1(isnan(data2))=nan;
data3(isnan(data2))=nan;

data1(isnan(data3))=nan;
data2(isnan(data3))=nan;

%% Start triple point collocation

% Calculate variance of each dataset
var_d1 = var(data1(:),'omitnan'); % Collapse 2-D matrix into 1-D by calling vector notation (:)
var_d2 = var(data2(:),'omitnan');
var_d3 = var(data3(:),'omitnan');

% Calculate the covariance of Data1 and Data3
covar_d1_d3 = diag(cov(data1(:), data3(:),'omitrows'),1); % covariance of data1 and data3 is off diagonal element

% Calculate the covariance of Data2 and Data3
covar_d2_d3 = diag(cov(data2(:), data3(:),'omitrows'),1); % covariance of data2 and data3 is off diagonal element

% Calculate the covariance of Data1 and Data2
covar_d1_d2 = diag(cov(data1(:), data2(:),'omitrows'),1); % covariance of data1 and data2 is off diagonal element

% Calculate the unscaled error variances for each dataset
% Covariance is symmetrical, so cov(d1,d2)==cov(d2,d1)
errorvar_d1 = (var_d1-((covar_d1_d2*covar_d1_d3)/covar_d2_d3));
errorvar_d2 = (var_d2-((covar_d1_d2*covar_d2_d3)/covar_d1_d3));
errorvar_d3 = (var_d3-((covar_d1_d3*covar_d2_d3)/covar_d1_d2));

% % Determine rescaling parameters if desired
% % Rescale to d1
% rescale_d2 = (covar_d1_d3/covar_d2_d3);
% rescale_d3 = (covar_d1_d2/covar_d2_d3);
% 
% errorvar_d1 = errorvar_d1;
% errorvar_d2 = errorvar_d2*(rescale_d2^2);
% errorvar_d3 = errorvar_d3*(rescale_d3^2);

% Calculate RMSD from the unscaled error variances
rmsd_d1 = sqrt(errorvar_d1);
rmsd_d2 = sqrt(errorvar_d2);
rmsd_d3 = sqrt(errorvar_d3);

% It is possible to have negative unscaled variances which means imaginary
% rmsd. Only keep positive variances, real rmsd.
% Replace imaginary numbers with NaN
rmsd_d1(imag(rmsd_d1)~=0) = NaN; 
rmsd_d2(imag(rmsd_d2)~=0) = NaN; 
rmsd_d3(imag(rmsd_d3)~=0) = NaN;

rmsd_d1 = real(rmsd_d1);
rmsd_d2 = real(rmsd_d2);
rmsd_d3 = real(rmsd_d3);



