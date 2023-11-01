# SVDS and PiMEP
<i><b>Software and documents developed in support of SVDS and the PiMEP NASA-ESA Collaboration</i></b>

<br>
The Salinity Validation Data System (<a href="https://www.esr.org/data-products/svds/">SVDS</a>) was developed at <a href="https://www.esr.org/">Earth and Space Research (ESR)</a> to provide systematic estimation and assessment of satellite sea surface salinity over the global ocean, including but not limited to match-ups with in situ data and triple point collocation analysis.
<br>

<br>
ESA’s Salinity Pilot Mission Exploitation Platform (<a href="https://www.salinity-pimep.org/">Salinity Pi-MEP</a>) is an innovative platform that allows users to evaluate the quality of satellite salinity data using matchups within situ data, statistical tools, mapping and visualization tools, and automated report generation.  
<br><br>
Scientists at ESR continue to develop tools for calibration-validation of NASA Aquarius/SAC-D and NASA Soil Moisture Active Passive (SMAP) data and are providing them to the Pi-MEP platform as part of this NASA-ESA collaboration. 
<br><br>

<hr>
In this repository, we provide the software and documents developed in support of SVDS and the PiMEP project. Contents of each folder are described below. 

<h3>Documents:</h3>   
<ul>
<li><b>Comparing Satellite Salinity Retrievals with In Situ Measurements: A Recommendation for Aquarius and SMAP </b> (Version 1). Schanze, Julian J., David M. LeVine, Emmanuel P. Dinnat, Hsun‐Ying Kao (2020). <a href="https://doi.org/10.5281/zenodo.4769713"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4769713.svg" alt="DOI"></a> </li>
<li><b>NASA/RSS SMAP Salinity Version 5.0 Validation Analysis </b> Kao, Hsun-Ying, Anderson, Jesse E., Schanze, Julian J., Howard, Susan (2023). <a href="https://doi.org/10.5281/zenodo.8368125"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.8368125.svg" alt="DOI"></a> </li> 
</ul>

<h3>SVDS:   </h3>  
<ul>
<li><b>Coming Soon!</b> Salinity Data Validation System Code</li> 
</ul>

<h3>TripleCollocation:   </h3>  
<ul>
<li><b>covariance_triple_point_collocation.m:</b>   Basic Triple Point Collocation Analysis for any combination of datasets using covariance notation</li> 
<li><b>triple_point_usecase_monthlyL3.m:</b>   Example usecase for triple collocation analysis with Monthly, L3 SMAP, SMOS, and RG Argo </li> 
<li><b>Doumentation_TriplePointCollocationCodeandExampleUseCase_v1.pdf:</b>   Documentation for implementation of triple_point_usecase_monthlyL3.m and covariance_triple_point_collocation.m with monthly, L3 SMAP, SMOS, and RG Argo. Documentation includes figures of the temporal and spatial RMSD for each dataset. Links to access each of the publicly available datasets are also provided. For a zipped repository of the usecase data (2.62GB), email janderson@esr.org </li> 
<li><b>covariance_triple_point_collocation.ipynb:</b>   Basic Triple Point Collocation Analysis for any combination of datasets using covariance notation 
<li><b>triple_point_usecase_monthlyL3.ipynb:</b>   Example usecase for triple collocation analysis with Monthly, L3 SMAP, SMOS, and RG Argo </li> 
</ul>

![spatial_readmeplot](https://user-images.githubusercontent.com/40212307/181094552-69cf8161-fa10-4e9a-807c-4a70603b69b6.jpg)

<h3>Scripts:   </h3>  
<ul>
  <li><b>triple_point_collocation.m:</b>   Triple Point Collocation Analysis for satellite (SMAP, SMOS) Level 3 and gridded Argo data. </li> 
</ul>
 


