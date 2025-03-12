import pandas as pd
import os
import sys
import traceback
import argparse
import gzip

def safe_read_file(file_path, read_mode='r'):
    """Safely read a file with comprehensive error checking."""
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File does not exist: {file_path}")
        
        # Check file size
        file_size = os.path.getsize(file_path)
        if file_size == 0:
            raise ValueError(f"File is empty: {file_path}")
        
        # Determine if it's a gzipped file
        is_gzipped = file_path.endswith('.gz')
        
        # Open file with appropriate method
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else read_mode
        
        with open_func(file_path, mode) as f:
            # Read first few lines to verify readability
            first_lines = [next(f) for _ in range(5)]
        
        print(f"File {file_path} checks passed.")
        print(f"First lines:\n{''.join(first_lines)}")
        
        return True
    
    except Exception as e:
        print(f"Error reading file {file_path}:")
        print(traceback.format_exc())
        return False

def read_barcode(barcode_file):
    try:
        # Comprehensive file checking
        if not safe_read_file(barcode_file):
            raise ValueError(f"Barcode file validation failed: {barcode_file}")
        
        # Read CSV with no header, handle 'x' case
        df = pd.read_csv(barcode_file, header=None)
        
        # Log dataframe details
        print(f"Barcode dataframe shape: {df.shape}")
        print(f"Barcode dataframe columns: {df.columns}")
        print(f"First few rows:\n{df.head()}")
        
        # Remove first row if it's 'x'
        if len(df) > 0 and df.iloc[0,0] == 'x':
            df = df.iloc[1:]
        
        # Extract barcodes, convert to string to handle potential type issues
        barcodes = set(df.iloc[:,0].astype(str))
        
        print(f"Total unique barcodes: {len(barcodes)}")
        print(f"First few barcodes: {list(barcodes)[:5]}")
        
        return barcodes
    
    except Exception as e:
        print(f"Comprehensive error in reading barcode file {barcode_file}:")
        print(traceback.format_exc())
        raise

def process_fragment_file(fragment_file, barcodes, output_file):
    try:
        # Comprehensive file checking
        if not safe_read_file(fragment_file):
            raise ValueError(f"Fragment file validation failed: {fragment_file}")
        
        # Read fragment file with error handling
        df = pd.read_csv(fragment_file, sep="\t", compression="gzip", 
                         names=['chr','start','end','barcode','count'])
        
        # Log dataframe details
        print(f"Fragment dataframe shape: {df.shape}")
        print(f"Fragment dataframe columns: {df.columns}")
        print(f"First few fragment rows:\n{df.head()}")
        
        # Validate barcode column
        if 'barcode' not in df.columns:
            raise ValueError("No 'barcode' column found in fragment file")
        
        # Filtering with type conversion to handle potential issues
        filtered_df = df[df['barcode'].astype(str).isin(barcodes)].iloc[:,:4]
        
        # Print filtering details
        print(f"Total fragments: {len(df)}")
        print(f"Filtered fragments: {len(filtered_df)}")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write output
        filtered_df.to_csv(output_file, sep="\t", index=False, header=False)
        print(f"Output written to {output_file}")
    
    except Exception as e:
        print(f"Comprehensive error processing fragment file {fragment_file}:")
        print(traceback.format_exc())
        raise

def main():
    parser = argparse.ArgumentParser(description='Subset fragment file by barcodes')
    parser.add_argument('--input', required=True, help='Input fragment file')
    parser.add_argument('--barcode', required=True, help='Barcode file')
    parser.add_argument('--output', required=True, help='Output fragment file')
    
    try:
        args = parser.parse_args()
        
        print("Subset Operation Details:")
        print(f"Input fragment file: {args.input}")
        print(f"Barcode file: {args.barcode}")
        print(f"Output file: {args.output}")
        
        # Read barcodes
        barcodes = read_barcode(args.barcode)
        
        # Process fragment file
        process_fragment_file(args.input, barcodes, args.output)
        
        print("Subset operation completed successfully.")
        return 0
    
    except Exception as e:
        print("Fatal error in subset operation:")
        print(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())