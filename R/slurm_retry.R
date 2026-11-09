# Retry batchtools' squeue/sacct status polls on transient slurmctld timeouts.
with_squeue_retry <- function(cf, tries = 6L, wait = 30) {
  wrap <- function(f) {
    force(f)
    function(reg) {
      for (i in seq_len(tries)) {
        out <- tryCatch(f(reg), error = identity)
        if (!inherits(out, "error")) return(out)
        # shortcut: matches Slurm's error wording; extend if MARS words timeouts differently
        transient <- grepl("Socket timed out|Unable to contact slurm controller|Connection refused",
                           conditionMessage(out))
        if (!transient || i == tries) stop(out)
        Sys.sleep(wait * i)   # 30, 60, 90, 120, 150 s: ~7.5 min before giving up
      }
    }
  }
  cf$listJobsQueued  <- wrap(cf$listJobsQueued)
  cf$listJobsRunning <- wrap(cf$listJobsRunning)
  cf
}
