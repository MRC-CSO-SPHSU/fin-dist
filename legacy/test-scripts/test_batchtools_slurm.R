## Minimal reproducer: can future.batchtools launch ONE SLURM job and return ONE result?
## Run from inside an srun shell on a compute node, after `conda activate quarto`:
##   Rscript test_batchtools_slurm.R
##
## A working run prints a list with worker_host/worker_pid/worker_time.
## A failing run will either error at plan() (template not found / packages missing)
## or hang at value() (sbatch not submitting / worker not returning).

suppressPackageStartupMessages({
  library(future)
  library(future.batchtools)
})

message("=== ENVIRONMENT ===")
message("hostname:           ", Sys.info()[["nodename"]])
message("SLURM_JOB_ID:       '", Sys.getenv("SLURM_JOB_ID"), "'")
message("Sys.which sbatch:   '", Sys.which("sbatch"), "'")
message("future ver:         ", as.character(packageVersion("future")))
message("future.batchtools:  ", as.character(packageVersion("future.batchtools")))
message("batchtools ver:     ", as.character(packageVersion("batchtools")))
message("===================")

plan(
  batchtools_slurm,
  template  = "slurm.tmpl",
  resources = list(
    ncpus    = 1,
    memory   = 2048,    # MB per CPU
    walltime = 300,     # seconds
    account  = "none"
  )
)

message("--- submitting one future via SLURM ---")
t0 <- Sys.time()
f  <- future({
  list(
    worker_host = Sys.info()[["nodename"]],
    worker_pid  = Sys.getpid(),
    worker_time = format(Sys.time())
  )
})

message("--- waiting for the result (this blocks until the worker completes) ---")
res <- value(f)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message("got result after ", elapsed, "s")
print(res)
message("--- done ---")
