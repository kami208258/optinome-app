import pandas as pd
import os

files = ["online_data_compilation.xlsx"]

for f in files:
    if os.path.exists(f):
        print(f"--- Inspecting {f} ---")
        try:
            xls = pd.ExcelFile(f)
            print("Sheet names:", xls.sheet_names)
            
            # check headers of the first sheet to see if it looks like data
            if xls.sheet_names:
                df = xls.parse(xls.sheet_names[0], header=None, nrows=5)
                print("First 5 rows of first sheet:")
                print(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    else:
        print(f"File {f} not found.")
