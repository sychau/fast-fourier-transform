import os
import subprocess
import json
import math


SEEDS = ["324234", "234234", "66756"]


# {EXP { SEED {ALGO {THREADS [time]}}}}
JSON_FILE = "parallel.json"
if (not os.path.exists(JSON_FILE)):
    with open(JSON_FILE, "w") as f:
        json.dump({}, f, indent=4)
complete_data = json.load(open(JSON_FILE, "r"))

TEMPLATE_SCRIPT = (
    "./../build/test_parallel --nThreads THREADS --sizeExp EXP_SIZE --seed SEED"
)


def run_test(exp_size, seed, threads):
    cmd = (
        TEMPLATE_SCRIPT.replace("EXP_SIZE", str(exp_size))
        .replace("SEED", str(seed))
        .replace("THREADS", str(threads))
    )
    print(f"Running: {cmd}")
    output = subprocess.check_output(cmd, shell=True).decode("utf-8")
    return output


def decode_output_and_save(output):
    lines = output.split("\n")

    # get the seed and exp_size
    threads = str(lines[0].split(": ")[1])
    size = str(math.log2(int(lines[1].split(": ")[1])))
    seed = str(int(lines[2].split(": ")[1]))

    # get the run times for each of the algorithms
    iterativeIcpSerial_time = str(float(lines[4].split(" ")[1]))
    iterativeIcp_time = str(lines[5].split(": ")[1].split(" ")[0])

    # read the rest of the data and check if any of the autovalidates failed
    failed = False
    for line in lines[6:]:
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
    if "iterativeIcp_multithreaded" not in complete_data[size][seed]:
        complete_data[size][seed]["iterativeIcp_multithreaded"] = {}
    if "iterativeIcp_Serial" not in complete_data[size][seed]:
        complete_data[size][seed]["iterativeIcp_Serial"] = {}
    if threads not in complete_data[size][seed]["iterativeIcp_multithreaded"]:
        complete_data[size][seed]["iterativeIcp_multithreaded"][threads] = []
    if threads not in complete_data[size][seed]["iterativeIcp_Serial"]:
        complete_data[size][seed]["iterativeIcp_Serial"][threads] = []

    # Save the data
    complete_data[size][seed]["iterativeIcp_multithreaded"][threads].append(iterativeIcp_time)
    complete_data[size][seed]["iterativeIcp_Serial"][threads].append(iterativeIcpSerial_time)

    # Save the data to the json file incase of failure later we still have the data until this point
    with open(JSON_FILE, "w") as f:
        json.dump(complete_data, f, indent=4)


def main(total_reapeated_tests=3):
    """Runs the test_serial program with the different sizes and seeds and saves the data to the
    json file and runs each version total_reapeated_tests times to be able to take the mean.
    """
    for run in range(total_reapeated_tests):
        for threads in range(1, 9): 
            for exp_size in range(4, 24):
                for seed in SEEDS:
                    output = run_test(exp_size, seed, threads)
                    decode_output_and_save(output)


if __name__ == "__main__":
    main(3)
