# Magnetic_storms_on_Earth-like_exoplanets_closer_to_their_host_stars
program codes and plot codes

You need to have the following files in your current directory to start the program:
#omni2_all_years.dat
#omni3save_july2016.p
#HELCATS_ICMECAT_v10_SCEQ.sav
they originate from the OMNI2 catalogue (ftp://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/extended) and 
the HELCATS ICMECAT catalogue (https://www.helcats-fp7.eu/catalogues/data/HELCATS_ICMECAT_v10_SCEQ.txt)

1. start the file "exoprogram_code.py"
it generates:
#exoprogram.txt
#dst_factors.tex
#scale_factors.txt
#factor_table.tex

2. start the plot programs "plotcode_single_ICME_event.py" and "plotcode_timeseries_histograms.py":
they generate the plots for the chosen single ICME event & the timeseries of ICME events as well as the Dst-index histograms for this timeseries 

