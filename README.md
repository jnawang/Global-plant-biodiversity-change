This folder contains R scripts used to predict global plant biodiversity redistribution and extinction risks under future (2081-2100) climate change using a series of methods. Below I document the objectives and sub-folder structure of this project. The number in sub-folder names generally indicate dependency relationships.

# Obectivies

-   Model the current and future suitable habitat distributions for \~ 68000 plant species under climate change.
-   Estimate species-specific range shift velocities to delineate their realized future species ranges.
-   Estimate plant species extinction rates and the effects of range shifts on extinction.
-   Identify global patterns of change in local plant species richness.

# Folder Structure

-   00_download_process_data
    -   00_01_process_species_ocurrence_data
    -   00_02_current_future_predictor_data_for_SDM
-   01_species_distribution_model
    -   01_00_SDM_history_validation
    -   01_01_Maxent
    -   01_02_RFDS
    -   01_03_Maxent_RFDS
-   02_predict_shift_rates
    -   02_01_process_data_for_shift_rate_prediction
    -   02_02_select_best_shift_rate_prediction_model
    -   02_03_perform_shift_rate_prediction
-   03_overlay_colonizable_suitable_ranges
    -   03_01_grid_distance_earth
    -   03_02_overlay_colonizable_suitable
-   04_merge_species_level_result
    -   04_01_merge_stat_files
    -   04_02_merge_area_shrink_files
-   05_get_global_extinction_pattern
-   06_get_diversity_change_pattern
    -   06_01_cal_richness
    -   06_02_grid_diversity_uncertainty

# Associated input data

-   Steps 00-03 require input data, which have been uploaded to Zenodo (DOI: ).

-   Steps 04-06 use the outputs from the previous steps. So we only provide key results files to Zenodo needed to reproduce our analyses, including data on species range size changes and richness changes under all modeled scenarios.

# Important note

Some of these scripts are designed to run on a computing cluster. Please make minor adjustments if you plan to run them on a personal computer.
