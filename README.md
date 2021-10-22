# FEPS-CFS
FEPS Complete Feature Set (FEPS-CFS) is a wrapper for the feature extraction tool previously developed by the KC Lab, also known as FEPS (Feature Extraction from Protein Sequences).  Currently, FEPS can extract a wide range of well published protein/DNA features.  However, if one wants to extract multiple features into one feature set, this must be done manually. The proposed solution is a wrapper that runs the FEPS processes on a loop through various user-defined parameters. The result is the ability to extract all features programmed into FEPS and place in a single feature vector file.
A comprehensive Feature Extraction from Protein Sequences (FEPS) web server was recently developed. FEPS uses published feature extraction methods of proteins from single or multiple-FASTA formatted files. In addition, FEPS also provides users with the ability to redefine some of the features by choosing one of the 544 physicochemical properties or to enter any set of user-defined amino acid indices, thereby increasing feature choices. The FEPS server includes 48 published feature extraction methods, six of which can use any of the 544 physicochemical properties. The total number of features calculated by FEPS is 2765, which exceeds the number of features computed by any other peer application. This exhaustive list of feature extraction methods enables researchers to develop machine learning-based approaches for various classification problems in bioinformatics. FEPS has been successfully applied for the prediction and classification of nuclear receptors, prediction of phosphorylation sites, and prediction of hydroxylation sites. 
FEPS is designed to extract various well published features and return them as a CSV file. However, in the original FEPS server, if researchers wanted to use multiple feature types to make predictions, they would need to manually combine the separate output files into one after extraction is finished.  To make it easier for researchers to use multiple feature types during the development of prediction algorithms, we developed a wrapper that allows concurrent runs of FEPS and places the output into a single CSV file. This gives users the ability to perform various analyses over the entire FEPS feature set. Importantly, this allows users to determine feature importance using the entire features set, which eliminates the need to try different feature sets to address a particular problem. 
## Installation
FEPS-CFS can be retrieved from https://github.com/dukkakc/FEPS_CFS.  Once downloaded, open a terminal and navigate to the repos root directory and run the following command:
```bash
python setup.py bdist_wheel
```
Then run the following command:
```bash
pip install /path/to/wheelfile.whl
```
To test that everything install properly run the following command from the FEPS_CFS root directory
```bash
python run_feps.py
```
## Usage
Once the FEPS library is installed, open a Python terminal or file and import it using:
```python
import FEPS.FEPS_CFS as feps
```
Next, we want to create two new directories (in the same directory as the Python program): One for input (ex: ‘infile’) and another for output (ex: ‘outfile’). Place the input fasta files inside the ‘infile’ folder (outfile will contain our output). After this, go back to the Python program and enter the following command:
```python
feps.feps('infile/','outfile/')#update in and out file here
```
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
