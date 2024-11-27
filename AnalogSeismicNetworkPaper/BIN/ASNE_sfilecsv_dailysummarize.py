#!/Users/thompsong/miniconda3/bin/python
# #!/usr/bin/env python
import os
import glob
import Seisan_Catalog

# Load the ASNE sfiles CSVs
SEISAN_DATA = os.environ['SEISAN_DATA']
#CSVpath = os.path.dirname(SEISAN_DATA)
#os.chdir(CSVpath)
list_of_csv_files = sorted(glob.glob('ASNE_sfiles??????.csv'))
eventtype_df, analyst_df = Seisan_Catalog.sfilecsv_daycount(list_of_csv_files)
print(eventtype_df)
print(analyst_df)
eventtype_df.to_csv('ASNE_dailyeventcounts.csv')
analyst_df.to_csv('ASNE_dailyanalystcounts.csv')
