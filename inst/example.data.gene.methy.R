### save example data with gene methylation
devtools::load_all();
data(example.data);
example.data.gene.methy <- gene.methylation(example.data);

usethis::use_data(example.data.gene.methy, overwrite = TRUE, compress = 'xz');
