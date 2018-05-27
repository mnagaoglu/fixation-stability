## Computing fixation stability using ISOA and BCEA methods

ISOA stands for isoline area and BCEA stands for bivariate contour ellipse area.

[example.m](example.m) contains two examples using synthetic data, (I) normally distributed data, (2) non-normal data.

[ComputeFixationStability.m](ComputeFixationStability.m) takes eye/gaze position data and desired cumulative probability,  and computes corresponding ISOA and BCEA values. It also computes preferred retinal locus (PRL) (although the term "retinal" can be debated if the data were obtained using a tracker without retinal imaging). PRL is computed in two ways: the location corresponding to the maximum probability density, or the mean eye/gaze position location.

[kde2d_trimmed.m](kde2d_trimmed.m) is a "hacked" version of [kde2d.m](https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation?focused=5829342&tab=function) by Zdravko Botev. "Hacking" was done by [Bosco Tjan](https://en.wikipedia.org/wiki/Bosco_Tjan).


## Examples

Normally distributed data
![exampleI](exampleI.tif)

Multimodal data
![exampleII](exampleII.tif)
