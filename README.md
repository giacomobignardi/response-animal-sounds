# Shared and unique taste in acoustic preferences

Materials and code to reproduce Bignardi., G., Vessel, E. A., Mosing, M. A., & Ullèn, F. (2026). Shared and unique taste in acoustic preferences [electronic response to James et al. (2026)]. This response was published in Science as an eLetter (<https://www.science.org/doi/10.1126/science.aea1202#elettersSection>).

Link to data, code, and materials for "Humans share acoustic preferences with other animals" (James et al., 2026): <https://github.com/themusiclab/animal-sounds>[\
\
](https://github.com/themusiclab/animal-sounds)

**Recommended citation**: Bignardi, G., Vessel, E.A, Mosing, M.A, & Ullén, F. Shared and unique taste in acoustic preferences [electronic response to James et al. (2026), Humans share acoustic preferences with other animals. Science, 391(6791), 1246-1249. <https://doi.org/10.1126/science.aea1202>].

## Directory structure

```         
R/
│ ├── reanalysis.R                       # Pulls data from github.com/themusiclab/animal-sounds and
│ │                                      # estimates partitioned repeatable variance into shared and unique via GLMM
│ └── vca.fun.R                          # Function to partition variance components into variance
│                                        # partitioning coefficients and beholder index
└── Bignardi_response_James_2026.pdf     # Letter including supplementary materials
```
