library(tidyverse)
library(ssh)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("genewalk")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

surface_fit <- syn("syn21444486") %>%
  read_csv()

# submit Genewalk job to cluster -----------------------------------------------
###############################################################################T

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- "/n/scratch2/ch305/jy/genewalk"
ssh_exec_wait(session, paste("mkdir", "-p", remote_wd, sep = " "))

write_csv(
  surface_fit %>%
    select(gene_id),
  file.path(wd, "surface_fit_genewalk_input.txt"),
  col_names = FALSE
)

scp_upload(session, here("scripts", "run_genewalk.sh"), remote_wd)
scp_upload(session, file.path(wd, "surface_fit_genewalk_input.txt"), remote_wd)

# submit job
ssh_exec_wait(
  session,
  paste(
    "sbatch",
    file.path(remote_wd, "run_genewalk.sh"),
    file.path(remote_wd, "surface_fit_genewalk_input.txt"),
    sep = " "
  )
)

ssh_disconnect(session)
