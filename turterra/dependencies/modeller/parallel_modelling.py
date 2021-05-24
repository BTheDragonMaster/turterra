import subprocess
import sys
import typing
from pathlib import Path

"""
This file should be copied to a scratch directory 
as all the temp files for modelling will be made in
the same folder. These will be automatically deleted
after modelling is complete, leaving only the PDB files
and a log file listing the pdb files and their NDOPE scores.
"""

path_type = typing.Union[str, Path]


def run_modelling(script, template_dir, template_file, ali_dir, model_dir, key_file, num_threads):
    """
    Runs homology modelling in parallel for all the keys in key_file
    Parameters
    ----------
    script
    template_dir
    template_file
    ali_dir
    model_dir
    key_file
    num_threads
    """
    keys = []
    with open(key_file) as f:
        for line in f:
            keys.append(line.strip())
    commands_list = []
    # key, index, template_dir, ali_dir, model_dir
    for i, key in enumerate(keys):
        commands_list.append(f"python {script} {key} {i} {template_file} {template_dir} {ali_dir} {model_dir}")
    num_commands = len(commands_list)
    num_waits = num_commands // num_threads + 1
    print("Num Waits", num_waits)
    start_job = 0
    for i in range(num_waits):
        processes = [subprocess.Popen(commands_list[j], shell=True) for j in range(start_job, min(start_job + num_threads, num_commands))]
        print(len(processes))
        for process in processes:
            process.wait()
        start_job += num_threads


def main():
    arguments = sys.argv
    assert len(arguments) == 8
    script, template_dir, template_file, ali_dir, model_dir, key_file, num_threads = arguments[1:]
    run_modelling(script, template_dir, template_file, ali_dir, model_dir, key_file, int(num_threads))


if __name__ == "__main__":
    main()
