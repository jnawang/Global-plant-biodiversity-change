# Merge predicted bioshift velocities
# Junna Wang, 2/15/2025
setwd("/Users/junnawang/UCDLab/Biodiversity/pre_post_processing/02_predict_shift_rate/02_03_perform_shift_rate_prediction/")
slow <- read.csv('Final_prediction_shift_rate_slow_2025.csv')
medi <- read.csv('Final_prediction_shift_rate_median_2025.csv')
fast <- read.csv('Final_prediction_shift_rate_fast_2025.csv')
#
all <- data.frame(Species=medi$Species,
                  disp=medi$disp,
                  scenario=medi$scenario,
                  Lat50=medi$Lat50,
                  Lat975=fast$Lat975,
                  Lat025=slow$Lat025,	
                  Ele50=medi$Ele50,	
                  Ele975=fast$Ele975,	
                  Ele025=slow$Ele025)
write.csv(all, 'predict_shift_rate2025.csv', row.names = F)