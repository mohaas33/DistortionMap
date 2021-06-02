#!/usr/bin/env python

# coding=utf-8

introduction = [
    "# All local jobs are part of the vanilla universe.",
    "Universe        = vanilla",
    "",
    "# The requirement line specifies which machines we want to",
    "# run this job on.  Any arbitrary classad expression can",
    "# be used.",
    "#Requirements    = (CPU_Speed >= 1)",
    "",
    "# Rank is an expression that states how to rank machines which ",
    "# have already met the requirements expression.  Essentially, ",
    "# rank expresses preference.  A higher numeric value equals better ",
    "# rank.  Condor will give the job the machine with the highest rank.",
    "#    Rank = CPU_Speed",
    "",
    "# Jobs by default get 1.4Gb of RAM allocated, ask for more if needed",
    "# but if a job needs more than 2Gb it will not be able to run on the",
    "# older nodes",
    "request_memory = 7.1GB",
    "",
    "# If you need multiple cores you can ask for them, but the scheduling",
    "# may take longer the \"larger\" a job you ask for",
    "request_cpus = 1",
    "",
    "# This flag is used to order only one's own submitted jobs ",
    "# The jobs with the highest numbers get considered for ",
    "# scheduling first.",
    "#Priority        = 4",
    "",
    "# Copy all of the user's current shell environment variables ",
    "# at the time of job submission.",
    "#GetEnv          = True",
    "",
    "# Used to give jobs a directory with respect to file input ",
    "# and output.",
    "Initialdir      = /sphenix/user/shulga/Work/IBF/DistortionMap/",
    "",
    "# Input file given to the job.",
    "#Input           = /dev/null",
    "",
    "",
    "# This should be the last command and tells condor to queue the",
    "# job.  If a number is placed after the command (i.e. Queue 15)",
    "# then the job will be submitted N times.  Use the $(Process)",
    "# macro to make your input/output and log files unique.",
    "Queue"
]

ff= open("./run_all_300evts_pp_jobs.sh","w+")
ff.write("#!/usr/bin/bash"+"\n"),


evt_start = [0,  480,959, 1438,1917,2396,2875,3355,3834,4313,4792,5272]
evt_end = [480,959,1438,1917,2396,2875,3355,3834,4313,4792,5272,5751]

evt_bX = [1508006.0, 3016012.0, 4524021.0, 6032017.0, 7540029.0, 9048024.0, 10556035.0,12064032.0,13572039.0,15080048.0,16588047.0,18096053.0]
for j, (start,end) in enumerate(zip(evt_start,evt_end)):
    for i in range(start,end+1):
        filename = "./condor_macros/run_files_300evts_pp_{}_{}.sh".format(j,i)
        f= open(filename,"w+")
        f.write("#!/usr/bin/bash"+"\n")
        f.write("source macros/run_files_300evts_pp.sh {} {} {}".format(i,i+1,evt_bX[j])+"\n")
        f.close
        filename_job = "./condor_macros/condor_run_files_300evts_pp_{}_{}.job".format(j,i)
        ff.write("condor_submit {}".format(filename_job)+"\n")
        f_job= open(filename_job,"w+")
        n_line = 0
        for lines in introduction:
            f_job.write(lines+"\n")
            if n_line==3:
                f_job.write("# The executable we want to run."+"\n")
                f_job.write("Executable      = condor_macros/run_files_300evts_pp_{}_{}.sh".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write(""+"\n")
                f_job.write("# The argument to pass to the executable."+"\n")
                f_job.write("Arguments       = \"run job 300 evts AA {} {}\"".format(j,i)+"\n")
            if n_line==38:
                f_job.write("# The job's stdout is sent to this file."+"\n")
                f_job.write("Output          = /sphenix/user/shulga/Work/IBF/DistortionMap/Out/myjob_300evts_pp_{}_{}.out".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The job's stderr is sent to this file."+"\n")
                f_job.write("Error           = /sphenix/user/shulga/Work/IBF/DistortionMap/Out/myjob_300evts_pp_{}_{}.err".format(j,i)+"\n")
                f_job.write(""+"\n")
                f_job.write("# The condor log file for this job, useful when debugging."+"\n")
                f_job.write("Log             = /sphenix/user/shulga/Work/IBF/DistortionMap/Out/condor_300evts_pp_{}_{}.log".format(j,i)+"\n")
                f_job.write(""+"\n")
    
            n_line+=1
        f_job.close
ff.close
