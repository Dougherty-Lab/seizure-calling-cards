#!/usr/bin/env python3
# calculate_profiles_parallel.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig
from pybedtools import BedTool
import scipy.stats as stats
import multiprocessing as mp
from functools import partial
import time

def process_peak_chunk_with_individual_data(peak_centers_chunk, bigwig_file, bin_size, bins_per_side):
    """
    Process a chunk of peaks in parallel, returning both profile data AND individual peak data
    
    Parameters:
    - peak_centers_chunk: List of (chrom, center) tuples for peaks
    - bigwig_file: BigWig file with signal values
    - bin_size: Size of each bin in bp
    - bins_per_side: Number of bins on each side of center
    
    Returns:
    - chunk_profile: Sum of signal values for this chunk (for profile plotting)
    - chunk_count: Count of valid values for this chunk (for profile plotting)
    - individual_peak_data: List of per-peak center/flank signals (for statistics)
    """
    # Open bigWig file
    bw = pyBigWig.open(bigwig_file)
    
    # Initialize arrays for profile data
    total_bins = bins_per_side * 2
    chunk_profile = np.zeros(total_bins)
    chunk_count = np.zeros(total_bins)
    
    # Initialize list for individual peak data
    individual_peak_data = []
    
    # Process each peak in this chunk
    for i, (chrom, center) in enumerate(peak_centers_chunk):
        # Arrays to store this peak's profile
        peak_profile = np.full(total_bins, np.nan)
        
        # Get profile data for plotting
        for j in range(total_bins):
            bin_start = center - (bin_size * bins_per_side) + (j * bin_size)
            bin_end = bin_start + bin_size
            
            # Skip if bin would have negative coordinates
            if bin_start < 0:
                continue
                
            try:
                # Get signal values for this bin
                values = bw.values(chrom, bin_start, bin_end)
                valid_values = [v for v in values if v is not None and not np.isnan(v)]
                
                if valid_values:
                    bin_signal = sum(valid_values) / len(valid_values)
                    chunk_profile[j] += bin_signal
                    chunk_count[j] += 1
                    peak_profile[j] = bin_signal
            except:
                # Handle potential errors (e.g., chromosome not in bigWig)
                continue
        
        # Calculate center and flank signals for this individual peak
        try:
            # Define center region (±500bp = ±2.5 bins for 200bp bins)
            center_bin_idx = bins_per_side  # Middle bin
            center_bins_start = max(0, center_bin_idx - 2)
            center_bins_end = min(total_bins, center_bin_idx + 3)
            
            # Define flank regions
            # Left flank: positions -7400 to -6400 bp
            left_flank_start_bin = max(0, int((center - 7400 - (center - bin_size * bins_per_side)) / bin_size))
            left_flank_end_bin = min(total_bins, int((center - 6400 - (center - bin_size * bins_per_side)) / bin_size))
            
            # Right flank: positions +6400 to +7400 bp  
            right_flank_start_bin = max(0, int((center + 6400 - (center - bin_size * bins_per_side)) / bin_size))
            right_flank_end_bin = min(total_bins, int((center + 7400 - (center - bin_size * bins_per_side)) / bin_size))
            
            # Extract signals
            center_signals = peak_profile[center_bins_start:center_bins_end]
            left_flank_signals = peak_profile[left_flank_start_bin:left_flank_end_bin]
            right_flank_signals = peak_profile[right_flank_start_bin:right_flank_end_bin]
            
            # Calculate means, excluding NaN values
            center_signal = np.nanmean(center_signals) if len(center_signals) > 0 else np.nan
            flank_signal = np.nanmean(np.concatenate([left_flank_signals, right_flank_signals])) if (len(left_flank_signals) + len(right_flank_signals)) > 0 else np.nan
            
            # Only include if both center and flank have valid signals
            if not (np.isnan(center_signal) or np.isnan(flank_signal)):
                individual_peak_data.append({
                    'center_signal': center_signal,
                    'flank_signal': flank_signal,
                    'chrom': chrom,
                    'center_pos': center
                })
                
        except Exception as e:
            # Skip this peak if there's an error
            continue
    
    # Close bigWig file
    bw.close()
    
    return chunk_profile, chunk_count, individual_peak_data

def calculate_profile_and_individual_data_parallel(peaks_file, bigwig_file, bin_size=200, bins_per_side=50, num_processes=None):
    """
    Calculate both the average signal profile AND individual peak data (parallel version)
    
    Parameters:
    - peaks_file: BED file with peak locations
    - bigwig_file: BigWig file with signal values
    - bin_size: Size of each bin in bp
    - bins_per_side: Number of bins on each side of center
    - num_processes: Number of processes to use (default: CPU count)
    
    Returns:
    - profile: Average signal profile (for plotting)
    - individual_data: List of per-peak center/flank signals (for statistics)
    """
    # Use all available CPUs if not specified
    if num_processes is None:
        num_processes = mp.cpu_count()
    
    # Load peaks
    peaks = BedTool(peaks_file)
    
    # Calculate center positions for each peak
    peak_centers = []
    for peak in peaks:
        center = int((int(peak.start) + int(peak.end)) / 2)
        peak_centers.append((peak.chrom, center))
    
    print(f"Processing {len(peak_centers)} peaks using {num_processes} processes")
    
    # Split peaks into chunks for parallel processing
    chunk_size = len(peak_centers) // num_processes
    if chunk_size == 0:
        chunk_size = 1
    peak_chunks = [peak_centers[i:i+chunk_size] for i in range(0, len(peak_centers), chunk_size)]
    
    # Create a partial function with fixed parameters
    process_chunk_partial = partial(
        process_peak_chunk_with_individual_data,
        bigwig_file=bigwig_file,
        bin_size=bin_size,
        bins_per_side=bins_per_side
    )
    
    # Process chunks in parallel
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_chunk_partial, peak_chunks)
    
    # Combine results from all processes
    total_bins = bins_per_side * 2
    profile = np.zeros(total_bins)
    count = np.zeros(total_bins)
    all_individual_data = []
    
    for chunk_profile, chunk_count, individual_data in results:
        profile += chunk_profile
        count += chunk_count
        all_individual_data.extend(individual_data)
    
    # Avoid division by zero for profile
    count[count == 0] = 1
    profile = profile / count
    
    return profile, all_individual_data

# Main analysis
def main():
    start_time = time.time()
    
    # File paths
    peaks_file = "peaks.bed"
    h3k27ac_bw = "ENCFF422JMO.bigWig"
    h3k27me3_bw = "ENCFF742UNH.bigWig"
    h3k4me1_bw = "ENCFF326LFP.bigWig"
    
    # Analysis parameters
    bin_size = 200  # 200bp bins
    bins_per_side = 50  # 10kb on each side
    
    # Use all available CPUs by default, or adjust as needed
    num_processes = mp.cpu_count()
    print(f"Using {num_processes} CPU cores for parallel processing")
    
    # Calculate profiles AND individual peak data for each histone mark
    print("Processing H3K27ac...")
    h3k27ac_profile, h3k27ac_individual = calculate_profile_and_individual_data_parallel(
        peaks_file, h3k27ac_bw, bin_size, bins_per_side, num_processes)
    
    print("Processing H3K27me3...")
    h3k27me3_profile, h3k27me3_individual = calculate_profile_and_individual_data_parallel(
        peaks_file, h3k27me3_bw, bin_size, bins_per_side, num_processes)
    
    print("Processing H3K4me1...")
    h3k4me1_profile, h3k4me1_individual = calculate_profile_and_individual_data_parallel(
        peaks_file, h3k4me1_bw, bin_size, bins_per_side, num_processes)
    
    # Save profile results to dataframe (for plotting)
    positions = np.arange(-10000, 10000, bin_size) + (bin_size / 2)
    results = pd.DataFrame({
        'position': positions,
        'H3K27ac': h3k27ac_profile,
        'H3K27me3': h3k27me3_profile,
        'H3K4me1': h3k4me1_profile
    })
    results.to_csv("histone_profiles.txt", sep="\t", index=False)
    
    # STATISTICAL ANALYSIS - Using peaks as n
    individual_data = {
        'H3K27ac': h3k27ac_individual,
        'H3K27me3': h3k27me3_individual,
        'H3K4me1': h3k4me1_individual
    }
    
    ratios = []
    for mark in ['H3K27ac', 'H3K27me3', 'H3K4me1']:
        peak_data = individual_data[mark]
        
        if len(peak_data) == 0:
            print(f"Warning: No valid data for {mark}")
            continue
            
        # Extract center and flank signals for each peak
        center_signals = [p['center_signal'] for p in peak_data]
        flank_signals = [p['flank_signal'] for p in peak_data]
        
        # Calculate summary statistics
        center_mean = np.mean(center_signals)
        flank_mean = np.mean(flank_signals)
        ratio = center_mean / flank_mean if flank_mean != 0 else np.inf
        
        # Perform t-test using individual peaks as observations
        # Use paired t-test since we're comparing center vs flank for the same peaks
        t_stat, p_value = stats.ttest_rel(center_signals, flank_signals)
        
        # Alternative: unpaired t-test if you prefer
        # t_stat, p_value = stats.ttest_ind(center_signals, flank_signals)
        
        # Calculate fold change
        fold_change = np.log2(ratio) if ratio > 0 and ratio != np.inf else np.nan
        
        ratios.append({
            'Mark': mark,
            'N_peaks': len(peak_data),
            'Center_mean': center_mean,
            'Flank_mean': flank_mean,
            'Center/Flank_ratio': ratio,
            'Log2_fold_change': fold_change,
            'T_statistic': t_stat,
            'P_value': p_value,
            'Significant': p_value < 0.05
        })
        
        print(f"{mark}: n={len(peak_data)} peaks, p={p_value:.2e}")
    
    # Save ratios with statistics
    pd.DataFrame(ratios).to_csv("histone_ratios.txt", sep="\t", index=False)
    
    # Save individual peak data for further analysis if needed
    for mark in ['H3K27ac', 'H3K27me3', 'H3K4me1']:
        if mark in individual_data:
            df = pd.DataFrame(individual_data[mark])
            df.to_csv(f"{mark}_individual_peak_data.txt", sep="\t", index=False)
    
    # Create simple plot for visualization
    plt.figure(figsize=(10, 6))
    plt.plot(positions, h3k27ac_profile, label='H3K27ac', color='blue')
    plt.plot(positions, h3k27me3_profile, label='H3K27me3', color='red')
    plt.plot(positions, h3k4me1_profile, label='H3K4me1', color='green')
    plt.axvline(-2000, linestyle='--', color='gray')
    plt.axvline(2000, linestyle='--', color='gray')
    plt.xlabel('Distance from peak center (bp)')
    plt.ylabel('Average signal')
    plt.title('Histone marks around CC peaks')
    plt.legend()
    plt.savefig("histone_profiles.png", dpi=300)
    
    # Create comparison of center vs flanking regions with statistics
    plt.figure(figsize=(10, 6))
    marks = ['H3K27ac', 'H3K27me3', 'H3K4me1']
    x = np.arange(len(marks))
    
    # Get center and flank means
    center_values = [next(item for item in ratios if item["Mark"] == mark)["Center_mean"] for mark in marks]
    flank_values = [next(item for item in ratios if item["Mark"] == mark)["Flank_mean"] for mark in marks]
    
    width = 0.35
    plt.bar(x - width/2, center_values, width, label='Peak Centers (±0.5kb)')
    plt.bar(x + width/2, flank_values, width, label='Flanking Regions (±6.4-7.4kb)')
    
    plt.ylabel('Average Signal')
    plt.title('Center vs Flanking Signal for Histone Marks (Statistics)')
    plt.xticks(x, marks)
    plt.legend()
    
    # Add p-values from t-tests
    for i, mark in enumerate(marks):
        p_value = next(item for item in ratios if item["Mark"] == mark)["P_value"]
        n_peaks = next(item for item in ratios if item["Mark"] == mark)["N_peaks"]
        asterisks = '*' * sum([p_value < 0.05, p_value < 0.01, p_value < 0.001])
        plt.text(i, max(center_values[i], flank_values[i]) + 0.1, 
                 f'{asterisks}\nn={n_peaks}', ha='center', va='bottom', fontweight='bold')
    
    plt.savefig("center_vs_flank.png", dpi=300)
    
    # Calculate execution time
    end_time = time.time()
    execution_time = end_time - start_time
    
    print("\nAnalysis complete. Results saved to:")
    print("  - histone_profiles.txt: Raw signal values for plotting")
    print("  - histone_ratios.txt: center/flank ratios with proper statistics")
    print("  - *_individual_peak_data.txt: Per-peak center/flank signals")
    print("  - histone_profiles.png: Line plot of signals")
    print("  - center_vs_flank.png: Bar chart with statistics")
    print(f"Total execution time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    main()