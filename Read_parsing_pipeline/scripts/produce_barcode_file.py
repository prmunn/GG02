# Produce barcodes.fa file from xlsx spreadsheet

# Author: Paul Munn, Genomics Innovation Hub, Cornell University

# Version history:
# Created: 03/24/2021


import os
import argparse
# import sys
# import datetime
from colorama import Fore
import pandas as pd
import numpy as np


# Set up globals


def main(args):
    print_header()

    data_folder = ''
    data_file_name = args.input
    outputFile = args.output
    testMode = args.testmode

    sheet_name = args.sheet
    index_column = args.index
    data_column = args.data
    bcLength = args.bclength
    duplicateBCType = args.duplicate
    skiprows = range(0, 1)  # Skip over the top row
    verifyIntegrityFlag = True
    prefix_chars = '^'
    suffix_chars = ''

    engine = 'openpyxl'  # Support for xlxs file format
    if data_file_name.split('.')[1] == 'xls':
        engine = 'xlrd'  # Support for xls file format
    df = pd.read_excel(data_folder + data_file_name,
                       sheet_name=sheet_name, skiprows=skiprows, engine=engine, keep_default_na=False)

    index_names = df[df[index_column] == ''].index
    df.drop(index_names, inplace=True)
    # df.dropna(axis=0, subset=[index_column], inplace=True)  # Remove nulls from index

    df.set_index(index_column, drop=False, inplace=True, verify_integrity=verifyIntegrityFlag)

    # Write to output file
    with open(outputFile, 'w') as outf:
        for index, row in df.iterrows():
            if row['In use'] == 1:
                barcode = row[data_column]
                if duplicateBCType == 'i5' and bcLength < 10:
                    barcode = barcode[:bcLength]
                if duplicateBCType == 'i7i5' and bcLength < 10:
                    barcode = barcode[:10] + barcode[10:(10+bcLength)]
                if testMode:
                    print('>' + index)
                    print(prefix_chars + barcode + suffix_chars)
                outf.write('>' + index + '\n')
                outf.write(prefix_chars + barcode + suffix_chars + '\n')

    success_msg('Done')


def success_msg(text):
    print(Fore.LIGHTGREEN_EX + text + Fore.WHITE)


def error_msg(text):
    print(Fore.LIGHTRED_EX + text + Fore.WHITE)


def print_header():
    print(Fore.WHITE + '*********************************************')
    print(Fore.GREEN + '        Produce barcodes.fa whitelist')
    print(Fore.WHITE + '*********************************************')
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Produce barcode.fa whitelist for use with cutadapt")
    parser.add_argument("input", help="Input spreadsheet containing barcodes and headers (required)")
    parser.add_argument("-s", "--sheet", help="Sheet name within input spreadsheet", default="barcode whitelist")
    parser.add_argument("-n", "--index", help="Index column within input spreadsheet", default="barcode name")
    parser.add_argument("-d", "--data", help="Sheet name within input spreadsheet", default="barcode seq")
    parser.add_argument("-b", "--bclength", help="Length of i5tag barcode (default: 8)", type=int, default=8)
    parser.add_argument("-u", "--duplicate", help="Type of barcodes to duplicate (i5 or i7i5)", choices=["i5", "i7i5"], default="i7i5")
    parser.add_argument("-o", "--output", help="Output file name", choices=["i5tagBC_barcodes.fa", "i7tagBC_i5tagBC_barcodes.fa"], default="i7tagBC_i5tagBC_barcodes.fa")
    parser.add_argument("-t", "--testmode", help="Turn test mode on", action="store_true")
    args = parser.parse_args()

    main(args)
