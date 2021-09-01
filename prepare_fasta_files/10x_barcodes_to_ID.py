import pandas as pd
import numpy as np

ID10x = pd.read_csv("./B10x/B10x_sequences.csv")
barcode = pd.read_csv("./B10x/B10x_sequences_with_barcodes_UMI.csv", index_col=0)

ID10x = ID10x.rename(columns={"Unnamed: 0":"10x_ID"})

print(ID10x.dtypes)
print(barcode.dtypes)

complete = pd.merge(ID10x, barcode, left_on=["alpha", "beta", "peptide", "mhc"], right_on=["alpha", "beta", "peptide", "mhc"])

print(complete)
print(np.shape(ID10x), np.shape(barcode))

complete.to_csv("./B10x/B10x_all_sequences_with_10xID_barcode_binding_UMI.csv")
