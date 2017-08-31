### Normalization

Normalization refers to adjusting values measured on different scales to a more common standardized scale.  A proper normalization is important in the analysis, firstly to avoid any enrichment toward any probe type in the differential methylated analysis and secondly to reduce non-biological (technical) variability between samples.

Normalization method used is $method$.

The effect of normalization can be seen by generating a density bean plot using the normalized data and comparing it with the pre-normalization density plot.

![**_Density Bean Plot (pre-normalization)_**](images/densityBeanPlot.png)

The Post-Normalization Density Bean Plot looks like the folllowing:

![**_Density Bean Plot(post-normalization)_**](images/densityBeanPlot_norm.png) 

### Probe Filter

The probe filter step identifies poor quality probes and/or CpG sites in the dataset that should be excluded from downstream analysis.


