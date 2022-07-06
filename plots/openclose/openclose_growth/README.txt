There are 3 groups of data (with corresponding subfolders):

- mus_all: all museums
- mus_indep: independent museums
- mus_lauth: local authority museums

Meaning of files:

xxxx-1var = museums are grouped by a single variable

xxxx-2vars = museums are grouped by 2 variables


The subfolders contain growth data grouped by variables. Each map shows only one value,
e.g. countries/regions by size=small.

Fields:

* diffN: difference in N of museums between 2017 and 1960 (e.g. 10 = 10 new museums in the area)
* open1960: N of museums open in 1960
* open2017: N of museums open in 2017
* close2017: N of museums closed cumulatively until 2017
* growth: growth rate (%) of musemums in the area. If N museums = 0 in 1960, growth will be NA.
* closedpc: percentage of closed musemums in the area (%). If N museums = 0 in 2017, growth will be NA.

"jenks" refers to the clustering technique used to identify the groups in the maps.

years_open_stats_1960_2017_pc_change-1var.xlsx:
	This file contains complete growth stats for single variables.

years_open_stats_1960_2017_pc_change-2vars.xlsx:
	This file contains complete growth stats for several combinations of 2 variables.
	The column VARS shows which variables are combined (var1 and var2).

To identify the LADs, a reference map is available in uk_lad_2016_reference_map.pdf
(with English regional boundaries too).