import json
import matplotlib.pyplot as plt


PLOTS_DIR = "actual_plots"


def compute_means_parallel():
    with open("parallel.json", "r") as f:
        data = json.load(f)
    # {EXP { SEED {ALGO {THREADS [time]}}}} -> {EXP {ALGO {THREADS [time]}}}

    # Create the new data structure
    means = {}
    for exp in data:
        means[exp] = {}
        for seed in data[exp]:
            for algo in data[exp][seed]:
                if algo not in means[exp]:
                    means[exp][algo] = {}
                for threads in data[exp][seed][algo]:
                    if threads not in means[exp][algo]:
                        means[exp][algo][threads] = []

                    for value in data[exp][seed][algo][threads]:
                        means[exp][algo][threads].append(float(value))

    # Compute the means
    for exp in means:
        for algo in means[exp]:
            for threads in means[exp][algo]:
                means[exp][algo][threads] = sum(means[exp][algo][threads]) / len(
                    means[exp][algo][threads]
                )

    return means


def create_parallel_plots():
    means = compute_means_parallel()
    # {EXP {ALGO {THREADS [mean_time]}}}

    # Create time by threads for exp_size plots
    for exp in means:
        for algo in means[exp]:
            plt.plot(
                list(map(int, means[exp][algo].keys())),
                list(means[exp][algo].values()),
                ".-",
                label=algo,
            )
        plt.xlabel("Threads")
        plt.ylabel("Time")
        plt.title(f"Execution Time by Threads for exp_size = {exp}")
        plt.legend()
        plt.savefig(f"{PLOTS_DIR}/parallel_exp_size_{exp}.png")
        plt.close()

    # Create Speedup by exp_size for 1 to 8 threads
    for threads in range(2, 9):

        plt.plot(
            list(map(float, means.keys())),
            [
                means[exp]["iterativeIcp"][str(1)]
                / means[exp]["iterativeIcp"][str(threads)]
                for exp in means.keys()
            ],
            ".-",
            label=f"{threads} threads",
        )
        plt.xlabel("exp_size")
        plt.ylabel("Speedup")
        plt.title(f"Speedup by exp_size for {threads} threads")
        plt.legend()
        plt.savefig(f"{PLOTS_DIR}/parallel_speedup_{threads}.png")
        plt.close()


def compute_means_distributed():
    with open("distributed.json", "r") as f:
        data = json.load(f)
    # {EXP { SEED {ALGO {PROCESSES [time]}}}} -> {EXP {ALGO {PROCESSES [time]}}}

    # Create the new data structure
    means = {}
    for exp in data:
        means[exp] = {}
        for seed in data[exp]:
            for algo in data[exp][seed]:
                if algo not in means[exp]:
                    means[exp][algo] = {}
                for processes in data[exp][seed][algo]:
                    if processes not in means[exp][algo]:
                        means[exp][algo][processes] = []

                    for value in data[exp][seed][algo][processes]:
                        means[exp][algo][processes].append(float(value))

    # Compute the means
    for exp in means:
        for algo in means[exp]:
            for processes in means[exp][algo]:
                means[exp][algo][processes] = sum(means[exp][algo][processes]) / len(
                    means[exp][algo][processes]
                )

    return means


def create_distributed_plots():
    means = compute_means_parallel()
    # {EXP {ALGO {processes [mean_time]}}}

    # Create time by processes for exp_size plots
    for exp in means:
        for algo in means[exp]:
            plt.plot(
                list(map(int, means[exp][algo].keys())),
                list(means[exp][algo].values()),
                ".-",
                label=algo,
            )
        plt.xlabel("processes")
        plt.ylabel("Time")
        plt.title(f"Execution Time by processes for exp_size = {exp}")
        plt.legend()
        plt.savefig(f"{PLOTS_DIR}/parallel_exp_size_{exp}.png")
        plt.close()

    # Create Speedup by exp_size for 1 to 8 processes
    for processes in range(2, 9):

        plt.plot(
            list(map(float, means.keys())),
            [
                means[exp]["iterativeIcp"][str(1)]
                / means[exp]["iterativeIcp"][str(processes)]
                for exp in means.keys()
            ],
            ".-",
            label=f"{processes} processes",
        )
        plt.xlabel("exp_size")
        plt.ylabel("Speedup")
        plt.title(f"Speedup by exp_size for {processes} processes")
        plt.legend()
        plt.savefig(f"{PLOTS_DIR}/parallel_speedup_{processes}.png")
        plt.close()


if __name__ == "__main__":
    create_parallel_plots()
