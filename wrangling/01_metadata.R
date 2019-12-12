library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wrangled_syn <- "syn21432189"

# RNA-seq metadata -------------------------------------------------------------
###############################################################################T

meta_raw <- syn("syn21411559") %>%
  readxl::read_excel()

meta <- meta_raw %>%
  mutate(
    condition = paste0("time", Time, "dox", DOX, "conc", ERKi),
    Sample_ID = paste0("S", as.integer(str_extract_all(Sample_ID, "[0-9]+")))
  )
write_csv(
  meta,
  here("wrangling", "meta.csv")
)


# Store to synapse -------------------------------------------------------------
###############################################################################T


meta_wrangling_activity <- Activity(
  "Wranglinig metadata",
  used = c("syn21411559"),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/wrangling/01_metadata.R"
)

c(
  here("wrangling", "meta.csv")
) %>%
  synStoreMany(wrangled_syn, activity = meta_wrangling_activity)
