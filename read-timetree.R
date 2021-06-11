library(readr)
library(ape)
tree <- read.tree("~/aDNA/taxonomy/TimetreeOfLife2015.nwk.txt")
names <- read_delim("~/aDNA/taxonomy/names.dmp", "\t", escape_double = FALSE, 
		    col_names = FALSE, trim_ws = TRUE, 
		    col_types = cols(X2=col_skip(), X4=col_skip(), X5=col_skip(),
		                     X6=col_skip(), X7=col_skip(), X8=col_skip()))
colnames(names) <- c("taxid", "txt")
names$txt <- gsub(" ", "_", names$txt)
names$taxid <- paste("t", names$taxid, sep="")

asmbl <- read_delim("~/aDNA/taxonomy/current/db-assembly-summary.txt", "\t", escape_double = FALSE,
  col_names = c("assembly_accession", "bioproject", "biosample", "wgs_master",
                "refseq_category", "taxid", "species_taxid", "organism_name",
                "infraspecific_name", "isolate", "version_status",
                "assembly_level", "release_type", "genome_rep",
                "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm",
                "paired_asm_comp", "ftp_path", "excluded_from_refseq"),
col_types = cols(bioproject = col_skip(), biosample = col_skip(), wgs_master = col_skip(),
                 infraspecific_name = col_skip(), isolate = col_skip(), version_status = col_skip(),
                 assembly_level = col_skip(), release_type = col_skip(), genome_rep = col_skip(), 
                 seq_rel_date = col_skip(), asm_name = col_skip(), submitter = col_skip(), 
                 excluded_from_refseq = col_skip()), trim_ws = TRUE)
asmbl$taxid <- paste("t", asmbl$taxid, sep="")
asmbl$species_taxid <- paste("t", asmbl$species_taxid, sep="")
taxid <- unique(c(asmbl$taxid, asmbl$species_taxid))

names_in_ours <- names$taxid %in% taxid
names_in_tree <- names$txt %in% tree$tip.label
keep <- subset(names, names_in_ours & names_in_tree)

t_name <- keep$txt
names(t_name) <- keep$taxid

asmbl <- subset(asmbl, taxid %in% keep$taxid | species_taxid %in% keep$taxid)
asmbl$sname <- t_name[asmbl$species_taxid]

write.table(asmbl, "asmbl-timetree.txt", sep="\t", quote = FALSE, row.names = FALSE)

tree_in_ours <- tree$tip.label %in% asmbl$sname

tt <- drop.tip(tree, tree$tip.label[!tree_in_ours])
dist <- cophenetic.phylo(tt)

save(dist, file="timetree-dist.Rdata")
