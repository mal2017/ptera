This workflow performs basic WGS analyis.

In its current iteration, the primary output is `results/copies/copies.tsv`.

```
Strain	| sample_name	sequence	length	bases	median.cov	est.copies
DGRP_313	DGRP_313_SAMN00014255	R1A1-element	5356	852468	159.16	64.96326530612244
DGRP_313	DGRP_313_SAMN00014255	roo	9092	1420089	156.19	63.75102040816326
...
...
```

Any modifications should honor or extend this format to fit into other workflows for the project.
