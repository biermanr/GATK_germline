import re
import os

def collect_jobids_workdirs(slurm_log_path:str) -> list[dict]:
    """
    Take a slurm_log_path and return a list of jobs like:
    jobs = [{"job_id": 12345, "wordir": "/scratch/nf/a5/asdf50"}, {"job_id"...}]
    """
    job_workdir_re = re.compile("jobId: (?P<job_id>.*?);.*workDir: (?P<workdir>.*?) ")
    job_ids = set()
    jobs = []

    with open(slurm_log_path,'r') as slurm_log:
        for line in slurm_log:
            if "workDir" not in line:
                continue

            try:
                match = re.search(job_workdir_re,line)
                job_id = int(match.group("job_id"))
                workdir = match.group("workdir")

                if job_id in job_ids:
                    continue

                job_ids.add(job_id)
                job = {"job_id": job_id, "workdir": workdir}
                jobs.append(job)
            except:
                continue

    return jobs

def get_command(jobs: list[dict]) -> list[dict]:
    """
    Get the contents of .command.sh prompts for each job
    """
    new_jobs = []
    for job in jobs:
        path = os.path.join(job["workdir"],".command.sh")
        with open(path,"r") as f_in:
            for line in f_in:
                if not line.startswith("gatk"):
                    continue

                cmd = " ".join(line.strip().split()) #cleans whitespace
                job["cmd"] = cmd
                new_jobs.append(job)

    return new_jobs


def get_trace(jobs: list[dict]) -> list[dict]:
    """
    Get the contents of .command.trace for each job
    """
    new_jobs = []
    for job in jobs:
        path = os.path.join(job["workdir"],".command.trace")
        with open(path,"r") as f_in:
            for line in f_in:
                if line.startswith("realtime"):
                    miliseconds = int(line.strip().split("=")[1])
                    hours = miliseconds*(1/1000)*(1/60)*(1/60)
                    job["hours"] = hours

                if line.startswith("%cpu"):
                    cpu_eff = int(line.strip().split("=")[1])
                    while cpu_eff > 1:
                        cpu_eff /= 10
                    job["cpu_eff"] = cpu_eff

                if line.startswith("%mem"):
                    mem_eff = int(line.strip().split("=")[1])
                    while mem_eff > 1:
                        mem_eff /= 10
                    job["mem_eff"] = mem_eff

                if line.startswith("peak_rss"):
                    peak_rss = int(line.strip().split("=")[1])
                    peak_rss_gb = peak_rss/1e6
                    job["peak_rss_gb"] = peak_rss_gb

        new_jobs.append(job)

    return new_jobs

if __name__ == '__main__':
    slurm_log_path = ".nextflow.log"
    jobs = collect_jobids_workdirs(slurm_log_path)
    jobs = get_command(jobs)
    jobs = get_trace(jobs)

    with open("parsed_job_info.tsv","w") as f_out:
        fields = ["job_id", "hours", "cpu_eff", "mem_eff", "peak_rss_gb", "workdir", "cmd"]
        f_out.write("\t".join(fields)+"\n")

        for job in jobs:
            j_info = "\t".join([str(job.get(field, None)) for field in fields])
            f_out.write(j_info+"\n")

