import os
import subprocess
import json
import math


SEEDS = ["324234", "234234", "66756", "85946", "87543095", "239576", "59456"]


# {EXP { SEED {ALGO [TIME]}}}
JSON_FILE = "serial.json"
complete_data = json.load(open(JSON_FILE, "r"))

TEMPLATE_SCRIPT = "./../build/test_serial --sizeExp EXP_SIZE --seed SEED"


def run_test(exp_size, seed):
    cmd = TEMPLATE_SCRIPT.replace("EXP_SIZE", str(exp_size)).replace("SEED", str(seed))
    print(f"Running: {cmd}")
    output = subprocess.check_output(cmd, shell=True).decode("utf-8")
    return output


def decode_output_and_save(output):
    lines = output.split("\n")

    # get the seed and exp_size
    size = str(math.log2(int(lines[0].split(": ")[1])))
    seed = str(int(lines[1].split(": ")[1]))

    # get the run times for each of the algorithms
    fftw_time = str(float(lines[3].split(" ")[1]))
    recursive_time = str(float(lines[4].split(" ")[1]))
    iterative_time = str(float(lines[5].split(" ")[1]))
    iteerativeIcp_time = str(float(lines[6].split(" ")[1]))

    # read the rest of the data and check if any of the autovalidates failed
    failed = False
    for line in lines[7:]:
        if "failed" in line:
            print(f"Failed: {line}")
            failed = True

    if failed:
        print("One of the tests failed. This data will not be saved.")
        return None

    # Ensure correct data structure
    if size not in complete_data:
        complete_data[size] = {}
    if seed not in complete_data[size]:
        complete_data[size][seed] = {}
    if "fftw" not in complete_data[size][seed]:
        complete_data[size][seed]["fftw"] = []
    if "recursive" not in complete_data[size][seed]:
        complete_data[size][seed]["recursive"] = []
    if "iterative" not in complete_data[size][seed]:
        complete_data[size][seed]["iterative"] = []
    if "iterativeIcp" not in complete_data[size][seed]:
        complete_data[size][seed]["iterativeIcp"] = []

    # Save the data
    complete_data[size][seed]["fftw"].append(fftw_time)
    complete_data[size][seed]["recursive"].append(recursive_time)
    complete_data[size][seed]["iterative"].append(iterative_time)
    complete_data[size][seed]["iterativeIcp"].append(iteerativeIcp_time)

    # Save the data to the json file incase of failure later we still have the data until this point
    with open(JSON_FILE, "w") as f:
        json.dump(complete_data, f, indent=4)


def main(total_reapeated_tests=5):
    """Runs the test_serial program with the different sizes and seeds and saves the data to the 
    json file and runs each version total_reapeated_tests times to be able to take the mean."""
    for run in range(total_reapeated_tests):
        for exp_size in range(10, 24):
            for seed in SEEDS:
                output = run_test(exp_size, seed)
                decode_output_and_save(output)


if __name__ == "__main__":
    main()
