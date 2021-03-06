# Ohia Fungal Endophytes

Fungal endophyte (FE) communities may be influenced not only by abiotic environmental conditions, but also by varying degrees of affinity to their host plants.  The landscape-dominant woody genus Metrosideros (Myrtaceae) contains many morphologically distinct taxa that occupy different habitats throughout the Hawaiian Islands.  This study used Metrosideros along an elevation gradient on O‘ahu to test the relative importance of environment versus host taxon on FE composition and diversity.  Variation in FE communities due to geographic distance (across and within sites) was also examined.  The fungal ITS1 region was sequenced from surface-sterilized leaf samples (N=176). Analyses revealed that variation in FE diversity was significantly explained by Metrosideros taxon, site, and geographic distance.  Considerable overlap in FE communities among host taxa and among sites was detected, however, and evidence for host-specificity of leaf FEs was weak and restricted to 700-1,000 m above sea level.  FE communities did not vary with elevation (environment); however, the elevation ranges examined may be too narrow for the detection of elevation/environmental effects.  Lastly, a significant pattern of isolation by distance on FE community composition was detected across the island as well as within each of the four sites.  These results suggest that within O‘ahu Metrosideros, leaf FE communities vary in diversity and composition across space, some of this variation is associated with host taxonomic effects and distance, and very little is associated with environmental variation across the island's elevation gradients.

___


### This repository contains all analysis code and output files for the project

|  File/Directory       	|  Contents                                                                                                                       	|
|-----------------------	|---------------------------------------------------------------------------------------------------------------------------------	|
| metadata.csv          	| Sample metadata                                                                                                                 	|
| 01_Process_Raw_Data.R 	|  Process raw Ion-Torrent reads, error correction, contaminant removal, exact sequence variant calling, generate sequence table 	|
| 02_Ohia_Analyses.R    	| Data analyses, statistical tests, figure generation                                                                             	|
| Ohia_Fungi.Rproj      	| R-Project file                                                                                                                  	|
| taxonomy/             	| Custom ITS1 taxonomic database, UNITE + Outgroups                                                                               	|
| output/               	| All statistical tables and figures                                                                                              	|

### Raw sequence data can be found on the SRA at Accession: [PRJNA606574](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA606574)


___


MIT License

Copyright (c) 2020 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
