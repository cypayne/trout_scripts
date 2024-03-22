## read quality checks

## Look for strangeness in basepair contents
# samtools stats is run in standard snakemake pipeline
# check the ACGT read content in samtools stats output
# gawk '/^GCT/ {print FILENAME, $0}' *.txt  > ACGT_read-content.txt
# if you see weird bp oscillation patterns, this could be
# due to low seq depth/low mapping rate
# check this by looking at the multiqc output
# M mapped reads and samtools stats alignment scores

library(tidyverse)

dir <- "qc_summary_stats/ccgp_trout/qc/samtools_stats/"
x <- read_table(paste0(dir,"ACGT_read-content.txt"), col_names = c("file", "tag", "cycle", "A", "C", "G", "T"))
y <- x %>% pivot_longer(A:T, names_to = "base",  values_to = "ppn")
ggplot(y, aes(x = cycle, y = ppn, colour = base)) + geom_line() + facet_wrap(~ file)
#ggsave(filename = "base-comps.pdf", width = 20, height = 20)


## ripping jsons from fastp
library(jsonlite)
dir <- "qc_summary_stats/test_jsons/fastp/"
lapply(dir, read_json)

files <- list.files(path = dir, pattern = ".json")
jsons <- paste0(dir, files) 
qc_list <- lapply(jsons, read_json)
names(qc_list) <- sub(".json", "", basename(jsons))

seq_content <- 
  data.frame(
  "A" = unlist(qc_list$`s001---1`$read1_before_filtering$content_curves$A),
  "C" = unlist(qc_list$`s001---1`$read1_before_filtering$content_curves$C),
  "G" = unlist(qc_list$`s001---1`$read1_before_filtering$content_curves$G),
  "T" = unlist(qc_list$`s001---1`$read1_before_filtering$content_curves$T)
  )

content_curve <- unlist(qc_list$`s001---1`$read1_before_filtering$content_curves$A)
ggplot(aes(x = c(1:length(content_curve)), y = content_curve)) + geom_line() 



# alternatively
jsons <- paste0(dir, files) 
files <- list.files(path = dir, pattern = ".json")
names <- sub(".json", "", files)

# load files
library(rjson)
library(tidyjson)
for(sample in names){
  sample_name <- str_replace_all(sample, "-","_")
  filepath <- file.path(dir,paste0(sample,".json"))
  assign(sample_name, fromJSON(file=filepath))
  
  seq_content_x <- eval(parse(text = sample_name))$read1_before_filtering$content_curves %>% 
    as.tibble() %>%
    mutate(file = sample_name, 
           cycle = c(1:nrow(seq_content_x)))
  #seq_content_x$file <- sample_name
  #seq_content_x$cycle <- c(1:nrow(seq_content_x))
    
  seq_content_y <- seq_content_x %>% pivot_longer(A:G, names_to = "base",  values_to = "ppn")
  ggplot(y, aes(x = cycle, y = ppn, colour = base)) + geom_line() + facet_wrap(~ file)
}

s001___1$read1_before_filtering$content_curves %>% unlist(A)

sample <- `s001---1`
cycle <- c(1:length(content_curve))
content_curve <- unlist(sample$read1_before_filtering$content_curves)
curve_df <- cbind(cycle, content_curve)
ggplot(content_curve, aes(x = c(1:length(content_curve)), y = content_curve)) + geom_line() 



library(purrr)
library(tidyverse)
library(jsonlite)

dir <- "qc_summary_stats/test_jsons/fastp/"
files <- dir(dir, pattern = "*.json")

data <- files %>%
  map_df(~fromJSON(file.path(dir, .), flatten = TRUE))

temp <- fromJSON(file.path(dir,files[1]))
temp_df <- as.data.frame(temp)


