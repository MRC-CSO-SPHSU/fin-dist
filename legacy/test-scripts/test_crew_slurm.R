## Minimal reproducer: can crew.cluster launch ONE worker and return ONE result?
## If succesful, the worker will print explicit information about the worker (manual launch, mode, etc.)

suppressPackageStartupMessages(library(crew.cluster))

message("--- ENVIRONMENT CHECK ---")
message("hostname:        ", Sys.info()[["nodename"]])
message("SLURM_JOB_ID:    '", Sys.getenv("SLURM_JOB_ID"), "'")
message("http_proxy:      '", Sys.getenv("http_proxy"), "'")
message("https_proxy:    '", Sys.getenv("https_proxy"), "'")
message("Sys.which sbatch: '", Sys.which("sbatch"), "'")
message("crew.cluster ver: ", as.character(packageVersion("crew.cluster")))
message("crew ver:        ", as.character(packageVersion("crew")))
message("mirai ver:       ", as.character(packageVersion("mirai")))
message("-------------------------")

ctl <- crew_controller_slurm(
  name            = "probe",
  workers         = 1,
  tls = crew::crew_tls(mode = "automatic"),
  seconds_idle    = 60,
  options_cluster = crew_options_slurm(
    verbose = TRUE,
    script_lines = c(
      "#SBATCH --account=none", # always add this option
      "module purge",
      "module load apps/miniforge",
      'source "$(conda info --base)/etc/profile.d/conda.sh"',
      "unset http_proxy",
      "unset https_proxy",
      "conda activate quarto"
    ),
    cpus_per_task = 1,
    time_minutes  = 5,
    log_output    = "probe-%A.out",
    log_error     = "probe-%A.err"
  )
)

message("--- starting controller ---")
ctl$start()

message("--- explicit launch of 1 worker ---")
ctl$launch(n = 1L)
Sys.sleep(2)
message("immediately after launch, sacct shows:")
system("sacct -u $USER --starttime now-1minutes  --format=JobID%15,JobName%15,NodeList,State,Submit,Start -n | tail -2")

message("--- pushing simple string task ---")
ctl$push(
  name    = "probe_task",
  command = paste0("worker on ", Sys.info()[["nodename"]],
                   " pid=", Sys.getpid(),
                   " time=", format(Sys.time()))
)

message("--- summary immediately after push (should show 1 pending task) ---")
print(ctl$summary())

message("--- waiting up to 240s for ONE task to complete ---")
t0 <- Sys.time()
ok <- ctl$wait(mode = "one", seconds_timeout = 240)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message("wait() returned: ", ok, "  (elapsed: ", elapsed, "s)")

message("--- explicit collect() to pull mirai results into the schedule ---")
ctl$collect()

message("--- summary after collect ---")
print(ctl$summary())

message("--- popping the completed task ---")
res <- ctl$pop()
print(res)

if (!is.null(res)) {
  if (length(res$result) >= 1L && !is.null(res$result[[1]])) {
    message("PAYLOAD: ", res$result[[1]])
  }
  if (length(res$error) >= 1L && !is.na(res$error[[1]])) {
    message("WORKER ERROR: ", res$error[[1]])
  }
}

message("--- final SLURM state for any probe worker jobs ---")
system("sacct -u $USER --starttime now-5minutes --format=JobID%15,JobName%15,NodeList,State,ExitCode,Elapsed,Reason%30 -n | tail -2")

message("--- worker log files crew created, if they exist ---")
system("ls -la probe-*.out probe-*.err 2>/dev/null")

message("--- shutting down ---")
ctl$terminate()
message("done.")
