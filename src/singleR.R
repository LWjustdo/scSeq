
singler = CreateSinglerObject(sc@raw.data, annot = NULL, 's106', min.genes = 0,
  technology = "DropSeq", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = sc@ident, do.main.types = T,
  reduce.file.size = T, numCores = 64)
singler$meta.data$orig.ident = sc@meta.data$orig.ident
singler$meta.data$xy = sc@dr$tsne@cell.embeddings
singler$meta.data$clusters = sc@ident
singler.new = convertSingleR2Browser(singler)