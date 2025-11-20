# Calcium-Fiber-Photometry-Analysis

These Matlab scripts and functions were created for Michaelides' Lab to analyze fiber photometry from a TDT system across multiple animals to visualize raw data, normalize that data, add it to structures, analyze AUC in 5-minute buckets, and analyze spikes and events in 5-minute buckets for saline vs. drug conditions.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FILES INCLUDED

- GCaMP_Drug_Script.m          : main script to read, normalize, plot, and analyze AUC for photometry data
- Spikes_Events_Drug_Script.m    : script to find, plot, and bucket spikes and events before and after drug injection
- rawTDT1.m                      : function to read and plot raw data from TDT photometry
- packageFP3.m                   : function to add individual data structure outputs to a larger data structure
- downsampleFP1.m                : function to downsample data using MATLAB resample function
- normalizeFP3.m                 : function to normalize raw fiber photometry data to ΔF/F or z-score
- phaucFP1.m                     : function to calculate AUC before and after an injection and in specified time buckets
- combineFP3.m                   : function to combine numerical data between entries with the same identifiers to generate average values and SEM

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT DATA STRUCTURE

Input folder should be the filepath for a folder of raw data from the TDT system

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HOW TO RUN

1. Paste filepath as the first input into rawTDT1 in the first line of GCaMP_Drug_Script, run to visualize raw data
2. Use packageFP3 to create/add to a structure of raw data, add identifiers for each animal and experiment
3. Use Normalize GCaMP Data section to smooth, downsample, and normalize raw data to the isosbestic signal, specifying the time range of the baseline signal
4. If desired, use Manually Adjust Baseline to 0 section to set the average value of the baseline recording to 0
5. Use packageFP3 to create/add to a structure of normalized data
6. Use AUC for GCaMPs section to analyze pre/post-injection AUCs and AUCs of specified time buckets
7. Use packageFP3 to create/add to a structure of AUC data
8. If desired, run Combine/Average recordings section on AUC structure to combine trace data from multiple animals of the same indentifiers
9. If combined, run/modify Plot Average Traces Together section to show mean trace with SEM shaded of combined recordings
10. To find spike and event data, use Spikes_Events_Drug_Script
11. Index your raw data structure to provide individual entries for spike and event analysis
12. Define window length, baseline length, recording length
13. Use packageFP3 to create/add to a structure of spike and event data

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
OUTPUTS

- raw data structure          : structure containing raw sensor and isosbestic signal data, time vector, sampling rate, and TTL information
- normalized data structure   : structure containing normalized signal data, time vector, sampling rate, and trace
- AUC data structure          : structure containing bucket size, post-injection window specified, signal vector and AUC pre and post-injection, bucket labels for plotting, and time vector, signal vector, and AUC data for each bucket
- Spike/Event data structure  : structure containing baseline and post-injection signal, bucket labels, baseline and post-injection spike and event rates, signals in buckets, buckets of spike amplitudes, widths and rates, and buckets of event amplitudes, widths, areas, and rates, as well as average values for these measures and values as a percent of that measure for the total baseline

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CONTACT

Maintainer: Will Dunne

Organization: NIH/NIDA

Email: will.dunne@nih.gov
