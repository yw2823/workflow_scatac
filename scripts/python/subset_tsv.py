import pandas as pd
import os
import sys
import argparse

##Read in cluster-specific barcode files

def read_barcode(barcode_file):
    try:

        df = pd.read_csv(barcode_file,header=None)
        print(f"Column names in barcode file: {df.columns}")  # Debug print
        
        if len(df) > 0 and df.iloc[0,0] == 'x':
            df = df.iloc[1:]
        
        print(f"The first row is now: {df.iloc[0,0]}")
        barcodes = set(df.iloc[:,0])
        print(f"The first few rows are: {list(barcodes)[:5]}")
        
        return barcodes
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        sys.exit(1)
    
def process_fragment_file(fragment_file,barcodes,output_file):
    try: 
        df = pd.read_csv(fragment_file, sep= "\t", compression="gzip",
            names=['chr','start','end','barcode','count'])

        print(f"first few lines:\n{df.head(5)}")
        print("Filtering")
        filtered_df = df[df['barcode'].isin(barcodes)].iloc[:,:4]

        print(f"Number of reads for cluster: {len(filtered_df)}" )
        filtered_df.to_csv(output_file,
                    sep="\t", index=False, header=False)
        print(f"Output written to {output_file}")
    except Exception as e:
        print(f"Error processing fragment file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Subset fragment file by barcodes')    
    parser.add_argument('--input', required = True, help='Fragment file path')
    parser.add_argument('--barcode', required = True,help='Barcode file path')
    parser.add_argument('--output', required = True,help='Output file path')
    
    args = parser.parse_args()
    
    # Prioritize long-form arguments


    print(f"Fragment file is {args.input}")
    print(f"Reading barcodes from {args.barcode}")
    barcodes = read_barcode(args.barcode)
    print(f"Found {len(barcodes)} barcodes")

    process_fragment_file(args.input,barcodes,args.output)
    print("Done")

if __name__ == "__main__":
    main()