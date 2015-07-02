EvalAligments
-------------


```
usage: evaluate.py [-h] [-s SEQFILE] [-n NATIVE] [-c CONTACT] [-t THRESHOLD]
                   [-r REFORMAT] [-o OUTPUT]
                   alignment

Run alignment quality evaluation workflow. For given alignment it outputs PPV,
(normalized) number of contacts above score threshold, and maximum contact
score.

positional arguments:
  alignment             Input aligment file

optional arguments:
  -h, --help            show this help message and exit
  -s SEQFILE, --seqfile SEQFILE
                        Sequence file
  -n NATIVE, --native NATIVE
                        Reference pdb file to compare with
  -c CONTACT, --contact CONTACT
                        Path to contact predictor executable
  -t THRESHOLD, --threshold THRESHOLD
                        Contact score threshold
  -r REFORMAT, --reformat REFORMAT
                        Path to reformat.pl script from HHsuite
  -o OUTPUT, --output OUTPUT
                        Save output in csv format
```
