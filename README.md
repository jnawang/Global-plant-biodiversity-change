This folder contains data and R scripts used to predict global plant biodiversity redistribution and extinction under future (2081-2100) climate change using a series of methods. Below I document the objectives and sub-folder structure of this project. The number in sub-folder names generally indicate dependency relationships.

# Obectives

-   Estimate plant species extinction rates and global distribution of high-extinction-risk species
-   Obtain global patterns of changes in local plant species richness
-   Obtain global patterns of winner and loser species
-   Investigate what features make a species either winner or loser species

# Folder Structure

-   00_download_process_data
    -   00_01_process_species_ocurrence_data
    -   00_02_current_future_predictor_data_for_SDM
    -   00_03_climate_and_climate_change_in_species_ranges
-   01_species_distribution_model
    -   01_01_Maxent
    -   01_02_RFDS
    -   01_03_Maxent_RFDS
-   02_predict_shift_rates
    -   02_01_process_data_for_shift_rate_prediction
    -   02_02_select_best_shift_rate_prediction_model
    -   02_03_perform_shift_rate_prediction
-   03_overlay_colonizable_suitable_ranges
-   04_merge_species_level_result
-   05_get_extinction_result
-   06_get_diversity_result
