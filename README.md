# ssGSEA
python script to run ssGSEA
```
usage: calcGSEA_mp.py [-h] [--veracyte] [--v] [--n_procs N_PROCS] input_reads input_sigs_json output

CalcGSEA

positional arguments:
  input_reads        path to input data - expect tab-sep raw RNA reads, samples in columns, genes in rows. Genes need to have ENSG in the name otherwise row will be ignored
  input_sigs_json    signature json file
  output             folder for output files

options:
  -h, --help         show this help message and exit
  --alt_input_reads  use this argument if input data is organized as samples in rows, genes in columns of reads or microarray data
  --v                output a lot! info on screen
  --n_procs N_PROCS  number of processes to be used during ssGSEA calculation
```

## input_sigs_json format

Provide signature genes as a json of a python dict:

keys -> Name of signature (append _up & _down in case you have a splitt up and down part and want to auto gen the average: [up - down]/2 for output)
values -> List of strings: GENE_ID - ENSG_ID Example: "KRT13 - ENSG00000171401"

```
{
  "SIGNATURE1_NAME": [
    "GENE_ID - ENSG_ID",
],
  "SIGNATURE2_NAME_up": [
    "GENE_ID - ENSG_ID",
],
  "SIGNATURE2_NAME_down": [
    "GENE_ID - ENSG_ID",
]
}
```

## Requirements
Python >= 3.13
```
numpy>=2.1.2
pandas>=2.2.3
python-dateutil>=2.9.0.post0
pytz>=2024.2
scipy>=1.14.1
six>=1.16.0
tzdata>=2024.2
```



  
