import argparse
import sys
import os
import pandas as pd
from engine.pipeline import run_pipeline

def main():
    parser = argparse.ArgumentParser(description="Medicinal Chemistry ADME Triage Engine")
    parser.add_argument("input_sdf", help="Path to input SDF file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    
    args = parser.parse_args()
    
    input_path = args.input_sdf
    output_path = args.output_csv
    
    if not os.path.exists(input_path):
        print(f"Error: Input file '{input_path}' does not exist.")
        sys.exit(1)
        
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Run Pipeline
    df = run_pipeline(input_path)
    
    if not df.empty:
        df.to_csv(output_path, index=False)
        print(f"Results saved to {output_path}")
    else:
        print("Processing failed or generated empty results.")
        sys.exit(1)

if __name__ == "__main__":
    main()
