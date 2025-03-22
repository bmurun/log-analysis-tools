import re
import matplotlib.pyplot as plt
from collections import defaultdict
import os
import numpy as np

def generate_measurement_delay_plots():

    log_list = []
    log_list.append("../output_eval/seq1.txt")
    log_list.append("../output_eval/seq2.txt")
    log_list.append("../output_eval/seq3.txt")

    topic_pattern = re.compile(r"Last message time for (\w+) is now")
    delay_pattern = re.compile(r"Received a measurement that was ([\d.eE+-]+) seconds in the past")

    # Nested dict: topic -> file -> list of delays
    delays_by_topic_and_file = defaultdict(lambda: defaultdict(list))

    for log_file in log_list:
        current_topic = None

        basename = os.path.basename(log_file)
        name = os.path.splitext(basename)[0]    
        with open(log_file, 'r') as file:
            for line in file:
                topic_match = topic_pattern.search(line)
                if topic_match:
                    current_topic = topic_match.group(1)

                delay_match = delay_pattern.search(line)
                if delay_match and current_topic:
                    delay_sec = float(delay_match.group(1))
                    delay_ms = delay_sec * 1000
                    delays_by_topic_and_file[current_topic][name].append(delay_ms)


    output_dir = "delay_histograms"
    os.makedirs(output_dir, exist_ok=True)
    
    for topic, file_delays in delays_by_topic_and_file.items():
        # Flatten all delays for this topic to find a common bin range
        all_delays = [delay for delays in file_delays.values() for delay in delays]
        min_delay, max_delay = min(all_delays), max(all_delays)
        bins = np.linspace(min_delay, max_delay, 100)

        plt.figure()
        for filename, delays in file_delays.items():
            plt.hist(delays, bins=bins, alpha=0.5, label=filename, edgecolor='black', density=True)
        plt.title(f"Measurement Delay for Topic: {topic}")
        plt.xlabel("Delay (ms)")
        plt.ylabel("Relative Frequency")
        plt.legend(title="Log File")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{topic}_overlay_histogram.png"), dpi=300)


if __name__ == "__main__":
    generate_measurement_delay_plots()