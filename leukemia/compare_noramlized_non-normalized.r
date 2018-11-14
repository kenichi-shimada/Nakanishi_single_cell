setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/non-normalization")
load("non-normalization.rda")
mat.non <- mat
stat.non <- stats
rt.non <- rt

setwd("/n/groups/mitchison/Kenichi/collaboration/Nakanishi_single_cell/Shared_Shimada_harvard/normalized")
load("normalized.rda")
mat.norm <- mat
stat.norm <- stats
rt.norm <- rt

identical(mat.non,mat.norm)
identical(stat.non,stat.norm)
identical(rt.non,rt.norm)

