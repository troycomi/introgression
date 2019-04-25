# introgression
> Some sort of short, description

## Background
Things about science

## Installation
All required packages are specified in the conda environment located in 
`code/environment.yml`.  The introgression environment can be generated with
```
conda env create -f environment.yml
```
To access the command line bindings of the main class, install with pip with
```
conda activate introgression
pip install --editable .
```
while in the code directory.

## Usage

### Configuration
A set of initial parameters are provided in `code/config.yaml` which need to
be set specific for your system and dataset.

Strings of the form \_\_KEY\_\_
are substituted during execution and are used as a shortcut.  For example,
with 'output\_root' set to `/data/results`, the value `__OUTPUT_ROOT__/genes/`
becomes `/data/results/genes/`

Strings of the form {state} are used for wildcards within the code.  Their
location and surrounding characters can change, but the wildcard must be the
same.  For example, `blocks_{state}.txt` can be changed to
`{state}_with-block.txt` but not `blocks_{st}.txt`.

### Command Line
With the package installed and the conda environment activated, main methods
are accessed with the `introgression` command. Some documentation is provided
by adding the argument `--help` to introgression or any of its subcommands.

### introgression
Options include:
- --config: specify one or more configuration files.  Files are evaluated in
order.  Conflicting values are overwritten by the newest file. This allows a
base configuration for the system and analysis-specific configurations added
as needed.
- verbosity: set by varying the number of v's attached to the option, with 
`-v` indicating a log level of critical and `-vvvvv` indicating debug logging.
Available subcommands are:
- predict

## License
TBD
