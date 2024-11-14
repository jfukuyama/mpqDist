## code to prepare `gentry`
library(ape)
library(phyloseq)
library(phytools)
library(tidyverse)
library(taxonlookup)

load (url ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/gentry197.RData'))
gentry.coord <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/gentry.coord.txt', row.names = 1)

sp_counts = list()
for(i in 1:length(gentry197)) {
    sp_counts[[i]] = gentry197[[i]] |> melt() |> group_by(variable) |> summarise(count = sum(value)) |> head()
    sp_counts[[i]]$site = rownames(gentry.coord)[i]
}

full_sp_counts = dplyr::bind_rows(sp_counts)
gentry_wide = full_sp_counts |> pivot_wider(names_from = variable, values_from = count, values_fill = 0)
cleaner_names = colnames(gentry_wide[-1])
cleaner_names = str_replace(cleaner_names, " *$", replacement = "")
## remove cf, so we're going to be assigning the tentative species to what they look like
cleaner_names = str_replace(cleaner_names, pattern = "\\bcf|\\bcf\\.|\\bCF|\\bCF\\.", replacement = "")
## doing the same with aff
cleaner_names = str_replace(cleaner_names, pattern = "aff |aff\\.|AFF | AFF\\.", replacement = "")
cleaner_names = str_replace(cleaner_names, "\\(nn\\)", replacement = "")
cleaner_names = str_replace(cleaner_names, "  ", replacement = " ")
cleaner_names = str_replace(cleaner_names, " *$", replacement = "")
cleaner_names = str_replace(cleaner_names, " \\.", replacement = "")
cleaner_names = str_replace(cleaner_names, " var\\. .*$", replacement = "")
cleaner_names = str_replace_all(cleaner_names, "\"", replacement = "")
split_names = str_split(cleaner_names, pattern = " ")
cleaner_names[which(sapply(split_names, length) > 3)]

family = sapply(cleaner_names, function(name) str_split(name, pattern = " ")[[1]][1])
family = family |> tolower() |> tools::toTitleCase()
genus = sapply(cleaner_names, function(name) str_split(name, pattern = " ")[[1]][2])
genus = genus |> tolower() |> tools::toTitleCase()
species = sapply(cleaner_names, function(name) str_split(name, pattern = " ")[[1]][3])
species = species |> tolower()
species_df = data.frame(family = family, genus = genus, species = paste(genus, species, sep = "_"), species_trinomial = paste(family, genus, species, sep = "_"))

## fix the things identified by S.PhyloMaker
species_df$family[which(genus == "Azara" & family == "Flacourtiaceae")] = "Salicaceae"
species_df$family[which(genus == "Calophyllum" & family == "Clusiaceae")] = "Calophyllaceae"
species_df$family[which(genus == "Dracaena" & family == "Agavaceae")] = "Asparagaceae"
species_df$family[which(genus == "Eucryphia" & family == "Eucryphiaceae")] = "Cunoniaceae"
species_df$family[which(genus == "Laurelia" & family == "Monimiaceae")] = "Atherospermataceae"
species_df$family[which(genus == "Mangifera" & family == "Anacar\\diaceae")] = "Anacardiaceae"
species_df$family[which(genus == "Nothofagus" & family == "Fagaceae")] = "Nothofagaceae"
species_df$family[which(genus == "Thomandersia" & family == "Acanthaceae")] = "Thomandersiaceae"
species_df$family[which(genus == "Viburnum" & family == "Caprifoliaceae")] = "Adoxaceae"
species_df$family[which(genus == "Acer" & family == "Aceraceae")] = "Sapindaceae"
species_df$family[which(genus == "Aegiphila" & family == "Amaranthaceae/Verbenaceae")] = "Lamiaceae"
species_df$family[which(genus == "Agave" & family == "Agavaceae")] = "Asparagaceae"
species_df$family[which(genus == "Alangium" & family == "Alangiaceae")] = "Cornaceae"
species_df$family[which(genus == "Asclepias" & family == "Asclepiacaceae")] = "Apocynaceae"
species_df$family[which(genus == "Bauhinia" & family == "Clusiaceae")] = "Fabaceae"
species_df$family[which(genus == "Bridelia" & family == "Euphorbiaceae")] = "Phyllanthaceae"
species_df$family[which(genus == "Capparis" & family == "Capparidaceae")] = "Capparaceae"
species_df$family[which(genus == "Crateva" & family == "Capparidaceae")] = "Capparaceae"
species_df$family[which(genus == "Dalbergia" & family == "Clusiaceae")] = "Fabaceae"
species_df$family[which(genus == "Marsdenia" & family == "Asclepiadaceae")] = "Apocynaceae"
species_df$family[which(genus == "Phyllanthus" & family == "Euphorbiaceae")] = "Phyllanthaceae"
species_df$family[which(genus == "Sambucus" & family == "Caprifoliaceae")] = "Adoxaceae"
species_df$family[which(genus == "Schlegelia" & family == "Bignoniaceae")] = "Schlegeliaceae"




gentry_for_squashing = data.frame(t(gentry_wide[,-1]))
gentry_for_squashing$trinomial = species_df$species_trinomial
gentry_squashed = gentry_for_squashing |> group_by(trinomial) |> summarise_all(sum)
gentry_otu_squashed = gentry_squashed[,-1] |> t() |> data.frame()
rownames(gentry_otu_squashed) = gentry_wide$site
## every column of this data frame should be a unique trinomial name
colnames(gentry_otu_squashed) = gentry_squashed$trinomial

phylo <- read.tree("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/80497c52c8795cd1b432f51dd870fe16d5308733/PhytoPhylo")
nodes <- read.delim("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/80497c52c8795cd1b432f51dd870fe16d5308733/nodes",header=TRUE)
source("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/80497c52c8795cd1b432f51dd870fe16d5308733/R_codes%20for%20S.PhyloMaker")
result <- S.PhyloMaker(spList = species_df, tree = phylo, nodes = nodes, scenarios = "S1")
tr = result$Scenario.1

not_in_tree = which(!(colnames(gentry_otu_squashed) %in% subset(result$Species.list, status != "unmatch")$species_trinomial))
## <2% of counts in unmapped species
sum(gentry_otu_squashed[,not_in_tree]) / sum(gentry_otu_squashed)
gentry_otu_squashed_binomial_names = gentry_otu_squashed[,-not_in_tree]
rownames(result$Species.list) = result$Species.list$species_trinomial
colnames(gentry_otu_squashed_binomial_names) = result$Species.list[colnames(gentry_otu_squashed_binomial_names),"species"]
out = lookup_table(colnames(gentry_otu_squashed_binomial_names), by_species=TRUE, missing_action = "NA")
out = as.matrix(out)

gentry = phyloseq(otu_table(gentry_otu_squashed_binomial_names, taxa_are_rows = FALSE),
                  phy_tree(tr),
                  sample_data(gentry.coord),
                  tax_table(out))

usethis::use_data(gentry, overwrite = TRUE)
