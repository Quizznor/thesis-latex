awk '{printf "%.4f\t%.1f\n", ($3+1970.*31556952.)/(31556952.), $9}' AugerApj2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo > time_exposure.dat
