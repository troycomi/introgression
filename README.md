# introgression
> Some sort of short, description

## Background
Things about science

## Installation
All required packages are specified in the conda environment located in 
`code/environment.yml`.  The introgression environment can be generated with
```bash
conda env create -f environment.yml
```
To access the command line bindings of the main analyze class,
install the setup file using pip with
```bash
conda activate introgression
pip install --editable .
```
while in the code directory.

## Usage

### Configuration
A set of initial parameters are provided in `code/config.yaml` which need to
be set specifically for your system and dataset.

Strings of the form \_\_KEY\_\_
are substituted during execution and are used as a shortcut.  For example,
with 'output\_root' set to `/data/results`, the value `__OUTPUT_ROOT__/genes/`
becomes `/data/results/genes/`.

Strings of the form {state} are used for wildcards within the code.  Their
location and surrounding characters can change, but the wildcard must be the
same.  For example, `blocks_{state}.txt` can be changed to
`{state}_with-block.txt` but not `blocks_{st}.txt`.

### Command Line
With the package installed and the conda environment activated, main methods
are accessed with the `introgression` command. Some documentation is provided
by adding the argument `--help` to introgression or any of its subcommands.

#### introgression
Options include:
- --config: specify one or more configuration files.  Files are evaluated in
order.  Conflicting values are overwritten by the newest file. This allows a
base configuration for the system and analysis-specific configurations added
as needed.
- verbosity: set by varying the number of v's attached to the option, with 
`-v` indicating a log level of critical and `-vvvvv` indicating debug logging.

Most subcommand options will overwrite corresponding values in the config
file.  Leaving options unset without supplying a value in the config file
will raise an error.  Some values are only set through the config file 
including the list of chromosomes and the known states.

Available subcommands are:
##### predict
The predict subcommand uses an HMM to predict regions of introgression from
alignment files.  Several outputs are used in subsequent steps which refine
the predicted introgressed regions.

Test strains to predict introgression on can be supplied in the config file
under the name 'strains' or pulled from the directory structure of
test\_strains.

Available options are:
- --alignment: input alignment file location with wildcards for 
{prefix} (optional), {strain} and {chrom}.
- --prefix: An optional wildcard value for alignment files.  If left blank,
will default to the known states joined with an underscore.  Leaving the
{prefix} wildcard out of the alignment file will prevent its use as well.
- --blocks: An output file containing the regions predicted to belong to the
given state.  Must contain the {state} wildcard which will be populated with a
known state during analysis.  Columns are the strain, chromosomes, the 
predicted state, start position, end position, and the number of sites 
supporting the assignment.
- --test-strains. If strains are not provided in the config, this file with
{strain} and {chrom} wildcards will be used to populate the strains for 
prediction.
- --hmm-initial: Output file with the initial parameters of the HMM for each
strain.
- --hmm-trained: Output file with HMM parameters following Baum-Welch training.
- --positions: Output file with indices of sites which are non-gapped sequences
which differ between reference alignments.
- --probabilities: Output file with the probability of each position belonging
to the master reference strain.
- --threshold: The threshold value to apply when filtering the predicted HMM
path through the test sequence.  Either a float, indicating cutoff probability,
or 'viterbi' to indicate the Viterbi algorithm should be used to find the most
likely sequence of states.
- --only-poly-sites/--all-sites: A switch to indicate if all non-gapped,
sequenced sites should be considered during HMM training, or only polymorphic
sites.  Default is only polymorphic sites.

##### id-regions
id-regions prepends a column to block files with a unique region id, of the
form 'r#'.  Regions are sorted by the start position of the region.  Changing
the states to label will affect the region numbers as a different set of
regions will be considered.

Available options are:
- --blocks: The input file to label with the wildcard {state}.  This is the 
file produced by predict in the previous step.
- --labeled: The output file, also containing {state} wildcard.
- --state: May be specified multiple times to indicate which states to add
labels to.  Leaving unset will use the states in the config file (recommended).

##### summarize-regions
Analyzes the regions predicted to be introgressed.  Several columns are added
to the block file containing information about the region including the number
of matching sites to each state.

Available options are:
- --state: May be set multiple times for each state to summarize.  Leaving
unset will default to all states in the config file.
- --labeled: the labeled block file with {state} wildcard created in the
previous step.
- --masks: Sequence mask files with {strain} and {chrom} wildcards.
- --alignment: The input alignment file similar to the predict option.
- --positions: The position file created during predict.
- --quality: The output file with a {state} wildcard.
- --region: The alignment for each region in the labeled file with {state}
wildcard.  Each state file contains all regions for the state.
- --region-index: A pickled dictionary used for random access into the region
file.  Must have {state} wildcard.

##### filter-regions
From the quality files produced in `summarize-regions`, filter regions based
on several criteria including those with weak support for the alternative
hypothesis and those which can be assigned to multiple alternative states.

Regions passing the 'introgressed filter' satisfy all of the following:
- fraction of gaps masked in reference > 0.5
- fraction of gaps masked in predicted state > 0.5
- number of matches to predicted > 7
- number of matches to predicted > number of matches to reference
- sequence identity with predicted state is higher than reference
- sequence identity with reference is > 0.7

Regions passing the 'ambiguous filter' match only the predicted state.  No
other state has:
- sequence identity >= sequence identity with predicted state * threshold
- matching bases >= matching bases with predicted state * threshold

Available options are:
- --region: The region file from summarize-regions.
- --region-index: The region index file from summarize-regions.
- --quality: The quality file produced by summarize-regions with {state}
wildcard.
- --introgress-filter: The output file with only regions passing introgression
filter.  Must contain {state} wildcard.
- --introgress-inter: An output file with all regions.  Includes the reason
for filtering by the introgression filter or '' if passes. Must contain {state}
wildcard.
- --ambiguous-filter: Output file containing only regions which pass the 
ambiguous filter after passing the introgression filter.
Must contain {state} wildcard.
- --ambiguous-inter: Contains all regions from introgression filter with a 
column for the reason the region failed ambiguous filtering.  Must contain
{state} wildcard.
- --thresh: The threshold to apply to the ambiguous filter.
- --filter-sweep: If set and threshold values are supplied as arguments,
will output summary information for applying the ambiguous filter with various
threshold values.

filter-regions also accepts multiple threshold values as arguments to test
and output to the filter-sweep file.  Sample usage would be
```bash
introgression --config config.yml \
    filter-regions --threshold 0.995 \
    --filter-sweep sweep.txt 0.99 0.98 0.8
```
where 0.99, 0.98 and 0.8 are used as test threshold values as summarized in
sweep.txt.  Note that the ambiguous filter will only use the threshold 0.995
in this example.

##### summarize-strains
Summarize strains produces summary information for each test strain including
the number of regions and bases assigned to each hidden state, filtered at
each stage, and ambiguous between states.

Available options are:
- --introgress-inter: The introgressed filter file as used in `filter-regions`.
- --ambiguous-inter: The ambiguous filter file as used in `filter-regions`.
- --strain-info: Tab separate table with information on the strain to include
with the summary output.  Columns should be the strain name, alternate name,
location, environment, and population.
- --state-counts: The summary output file.

## License
TBD
