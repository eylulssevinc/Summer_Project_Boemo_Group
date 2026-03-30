import sys
import numpy as np
import matplotlib.pyplot as plt


def parse_iod_output(filepath):
    headers = {}
    sections = {}
    current_section = None
    current_data = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                parts = line[1:].split(' ', 1)
                if len(parts) == 2:
                    headers[parts[0]] = parts[1]
            elif line.startswith('>'):
                if current_section is not None:
                    sections[current_section] = current_data
                current_section = line[1:].rstrip(':')
                current_data = []
            elif line and current_section is not None:
                current_data.append(line)

    if current_section is not None:
        sections[current_section] = current_data

    return headers, sections


def main():
    if len(sys.argv) != 2:
        print("Usage: python plotIOD.py <iod_output_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    headers, sections = parse_iod_output(filepath)

    median_iod = float(headers.get('MedianIOD', 0))
    ci = headers.get('95ConfidenceInterval', None)

    data_single = np.array([float(x) for x in sections.get('DataSingleForkDistances', [])])
    sim_single = np.array([float(x) for x in sections.get('SimSingleForkDistances', [])])
    data_pair = np.array([float(x) for x in sections.get('DataForkPairDistances', [])])
    sim_pair = np.array([float(x) for x in sections.get('SimForkPairDistances', [])])

    landscape_raw = sections.get('Landscape', [])
    landscape_fr, landscape_iod, landscape_w = [], [], []
    for row in landscape_raw:
        parts = row.split('\t')
        landscape_fr.append(float(parts[0]))
        landscape_iod.append(float(parts[1]))
        landscape_w.append(float(parts[2]))
    landscape_fr = np.array(landscape_fr)
    landscape_iod = np.array(landscape_iod)
    landscape_w = np.array(landscape_w)

    fig, axes = plt.subplots(2, 3, figsize=(20, 9))

    # Panel 1: Objective landscape
    ax = axes[0,0]
    ax.plot(landscape_fr, landscape_w, 'o-', color='#2166ac', markersize=3, linewidth=1)
    ax.set_xscale('log')
    ax.set_xlabel('Firing rate (fr)')
    ax.set_ylabel('Wasserstein distance')
    ax.set_title('Objective landscape')

    # Panel 2: Objective vs IOD
    ax = axes[0,1]
    ax.plot(landscape_iod, landscape_w, 'o-', color='#b2182b', markersize=3, linewidth=1)
    ax.set_xlabel('Median IOD (kb)')
    ax.set_ylabel('Wasserstein distance')
    ax.set_title('Objective vs IOD')
    ax.axvline(median_iod, color='grey', linestyle='--', linewidth=0.8, label=f'Best IOD = {median_iod:.0f} kb')
    if ci:
        lo, hi = [float(x) for x in ci.split()]
        ax.axvspan(lo, hi, alpha=0.15, color='grey', label=f'95% CI [{lo:.0f}, {hi:.0f}] kb')
    ax.legend(fontsize=8)

    # Panel 3: Histogram — single-fork distances
    ax = axes[0,2]
    if len(data_single) > 0 or len(sim_single) > 0:
        all_single = np.concatenate([d for d in [data_single, sim_single] if len(d) > 0])
        bins_s = np.linspace(0, np.percentile(all_single, 99), 40)
        if len(data_single) > 0:
            ax.hist(data_single, bins=bins_s, density=True, alpha=0.5, color='#2166ac', label='Data')
        if len(sim_single) > 0:
            ax.hist(sim_single, bins=bins_s, density=True, alpha=0.5, color='#b2182b', label='Simulated')
    ax.set_xlabel('Single-fork behind-distance (kb)')
    ax.set_ylabel('Density')
    ax.set_title('Single-fork histograms')
    ax.legend(fontsize=8)

    # Panel 4: CDF comparison — single-fork behind-distances
    ax = axes[1,0]
    if len(data_single) > 0:
        s = np.sort(data_single)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#2166ac', linewidth=1.5, label='Data')
    if len(sim_single) > 0:
        s = np.sort(sim_single)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#b2182b', linewidth=1.5, label='Simulated')
    ax.set_xlabel('Single-fork behind-distance (kb)')
    ax.set_ylabel('Cumulative probability')
    ax.set_title('Single-fork distances')
    ax.legend(fontsize=8)

    # Panel 5: CDF comparison — fork-pair EdU start distances
    ax = axes[1,1]
    if len(data_pair) > 0:
        s = np.sort(data_pair)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#2166ac', linewidth=1.5, label='Data')
    if len(sim_pair) > 0:
        s = np.sort(sim_pair)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#b2182b', linewidth=1.5, label='Simulated')
    ax.set_xlabel('Fork-pair EdU start distance (kb)')
    ax.set_ylabel('Cumulative probability')
    ax.set_title('Fork-pair distances')
    ax.legend(fontsize=8)

    # Panel 6: Histogram — fork-pair distances
    ax = axes[1,2]
    if len(data_pair) > 0 or len(sim_pair) > 0:
        all_pair = np.concatenate([d for d in [data_pair, sim_pair] if len(d) > 0])
        bins_p = np.linspace(0, np.percentile(all_pair, 99), 40)
        if len(data_pair) > 0:
            ax.hist(data_pair, bins=bins_p, density=True, alpha=0.5, color='#2166ac', label='Data')
        if len(sim_pair) > 0:
            ax.hist(sim_pair, bins=bins_p, density=True, alpha=0.5, color='#b2182b', label='Simulated')
    ax.set_xlabel('Fork-pair EdU start distance (kb)')
    ax.set_ylabel('Density')
    ax.set_yscale('log')
    ax.set_title('Fork-pair histograms')
    ax.legend(fontsize=8)

    plt.tight_layout()
    outpath = filepath.rsplit('.', 1)[0] + '_IODplots.pdf'
    plt.savefig(outpath)
    print(f"Saved to {outpath}")
    plt.show()


if __name__ == '__main__':
    main()

