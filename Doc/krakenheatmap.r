
library(tidyverse)

#Set your directory
setwd("")
# List files
files <- dir(pattern = "*.kraken2.species.tsv")
# Extract sample IDs from filenames
Names <- tibble(FileName = files) %>%
  mutate(sampleID = str_remove_all(FileName, ".kraken2..*"))

# Read in all tables
Tables <- map(files, ~ read_tsv(.,
                                col_names = c("Percentage", "TaxID", "SpeciesName"),
                                col_types = cols(.default = "c"))) %>%
  set_names(Names$sampleID)

Tables2Abundance <- map(Tables, ~ select(., SpeciesName, Percentage) %>%
                          mutate(Percentage = as.numeric(Percentage)))
AbundanceTable <- reduce(Tables2Abundance, full_join, by = "SpeciesName") %>%
  rename_with(~ Names$sampleID, -SpeciesName) %>%
  mutate(across(-SpeciesName, ~ replace_na(.x, 0)))

library(pheatmap)

AbundanceTable %>% 
  column_to_rownames("SpeciesName") %>%
  pheatmap()

MtrixForPH <- AbundanceTable %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
  column_to_rownames("SpeciesName") %>%
  as.matrix()


pheatmap(MtrixForPH)

pheatmap(MtrixForPH, cluster_row=F, cellwidth = 16)


library(viridis)
Color <- rev(inferno(500))
pheatmap(MtrixForPH, 
         cluster_row=F,
         cellwidth = 16,
         color=Color)



ColAnnot <- AbundanceTable %>%
  select(-SpeciesName) %>%
  colnames() %>%
  enframe(value = "Sample") %>%
  mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
  select(-name) %>% 
  column_to_rownames("Sample")


pheatmap(MtrixForPH, 
         cluster_row=F,
         cellwidth = 16,
         color=Color,
         annotation_col = ColAnnot)


library(RColorBrewer)

AnnotColor <- list(Method=AbundanceTable %>%
                     select(-SpeciesName) %>%
                     colnames() %>%
                     enframe(value = "Sample") %>%
                     mutate(Method=str_remove(Sample,"_\\d+.*")) %>%
                     select(-name) %>% 
                     column_to_rownames("Sample") %>%
                     distinct() %>%
                     mutate(Color=brewer.pal(3,"Dark2")) %>%
                     deframe())



KrakenPH <- pheatmap(MtrixForPH, 
                     cluster_row=F,
                     cellwidth = 16,
                     color=Color,
                     annotation_col = ColAnnot,
                     annotation_colors = AnnotColor)

KrakenPH


library(ggplotify)
KrakenPH <- as.ggplot(KrakenPH)

library(gridExtra)

QCK2HM <- grid.arrange(arrangeGrob(ReadLengthVP,N50VP),KrakenPH,ncol=2)
