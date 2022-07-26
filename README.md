# PiMEP
<i><b>Software and documents developed in support of the PiMEP NASA-ESA Collaboration and Salinity Validation Data System</i></b>

<br>
ESA’s Salinity Pilot Mission Exploitation Platform (<a href="https://www.salinity-pimep.org/">Salinity Pi-MEP</a>) is an innovative platform that allows users to evaluate the quality of satellite salinity data using matchups within situ data, statistical tools, mapping and visualization tools, and automated report generation.  
<br><br>
Scientists at <a href="https://www.esr.org/">Earth and Space Research (ESR)</a> continue to develop tools for calibration-validation of NASA Aquarius/SAC-D and NASA Soil Moisture Active Passive (SMAP) data and are providing them to the Pi-MEP platform as part of this NASA-ESA collaboration. 
<br><br>

<hr>
In this repository, we provide the software and documents developed in support of this project: 

<h3>Documents:</h3>   
<ul>
<li><b><i>Comparing Satellite Salinity Retrievals with In Situ Measurements: A Recommendation for Aquarius and SMAP </b></i></li>
</ul>

<h3>TripleCollocation:   </h3>  
<ul>
<li><b>covariance_triple_point_collocation.m:</b>   Basic Triple Point Collocation Analysis for any combination of datasets using covariance notation</li> 
<li><b>triple_point_usecase_monthlyL3.m:</b>   Example usecase for triple collocation analysis with Monthly, L3 SMAP, SMOS, and RG Argo </li> 
<li><b>Doumentation_TriplePointCollocationCodeandExampleUseCase_v1.pdf:</b>   Documentation for implementation of triple_point_usecase_monthlyL3.m and covariance_triple_point_collocation.m with monthly, L3 SMAP, SMOS, and RG Argo. Documentation includes figures of the temporal and spatial RMSD for each dataset. Links to access each of the publicly available datasets are also provided. For a zipped repository of the usecase data (2.62GB), email janderson@esr.org </li> 
<li><b>covariance_triple_point_collocation.ipynb:</b>   Basic Triple Point Collocation Analysis for any combination of datasets using covariance notation </li> 
<li><b>triple_point_usecase_monthlyL3.ipynb:</b>   Example usecase for triple collocation analysis with Monthly, L3 SMAP, SMOS, and RG Argo <i>coming july 2022</i></li> 
</ul>

<h3>Scripts:   </h3>  
<ul>
  <li><b>triple_point_collocation.m:</b>   Triple Point Collocation Analysis for satellite (SMAP, SMOS) Level 3 and gridded Argo data. </li> 
</ul>
 

