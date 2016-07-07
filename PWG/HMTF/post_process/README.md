# MC generator comparisons

This readme discribes how to run the post analysis on the the `AnalysisResults.root` file produced by `/PWG/HMTF/macros/AddTaskHMTFMCMultEst.C`

## Running the post analysis on lxplus

### Prepare the environment
The post processing is implemented in python. Unfortunately, the python installation of alienv on lxplus is broken at the time of writing. On the other hand, the code downloads the AnalysisResults.root files with alien, so using the default root version is not an option either. Hence, it is necessary to prepare the shell environment before running the post processing. This is done by running the following while __not__ being in an alienv:

	$ cd ~
	$ source /afs/cern.ch/user/c/cbourjau/poor_man_alien.sh
	$ # Create a new alien token:
	$ alien-token-init

This makes the executable `hmtfmc` avialable on the command line. It is self documented:

``` shell
$ hmtfmc --help
usage: hmtfmc [-h] [-v] {prepare_plots,summarize,compare} ...

positional arguments:
  {prepare_plots,summarize,compare}

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose
```

## Preparing plot for later steps

This steps produces all plots and saves them in the downloaded AnalysisResults.root file. This step needs to be taken before creating the summary or comparison slides. For example:

``` shell
$ mkdir ~/hmtfmc_results
$ cd ~/hmtfmc_results
$ hmtfmc prepare_plots --help
usage: hmtfmc prepare_plots [-h] input_file

This will download AnalysisResults.root from the given alien path into the
current directory. The file name is composed of the train id and generator
name as given in the trains `env.sh` file. Subsequently, "all" plots are
created and saved in to the AnalysisResults.root file. This step needs to be
taken before creating summary or comparison PDFs. This operation requrires a
valid alien-token!

positional arguments:
  input_file  Path to file containing the input data

optional arguments:
  -h, --help  show this help message and exit

$ hmtfmc prepare_plots /alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/393_20160210-1718/merge/AnalysisResults.root
```
	
Now, there is a file called 393_AnalysisResults.root (393 being the train number). The file contains all the plots at this points

## Making the summary slides for one generator

Assuming that the previously downloaded file was for "Dipsy 13 TeV Ropes", the following will create the summary slides. Again, `--help` is available.

``` shell
$ hmtfmc summarize --help
usage: hmtfmc summarize [-h] [--gen_name GEN_NAME] input_file

Create summary slides for all triggers of a given generator. Remember to use
quotes if the generator name contains spaces and be careful with special
characters.

positional arguments:
  input_file           Path to file containing the input data

optional arguments:
  -h, --help           show this help message and exit
  --gen_name GEN_NAME  Name of the generator and tune. If not given, deduced
                       from filename

$ hmtfmc summary 393_AnalysisResults.root "Dipsy 13TeV Ropes"
```

the summary slides can now be found in the aptly called folders created in the current directory

## Making comparison slides
Two generators and/or tuning can be compared as follows

```shell
$ hmtfmc compare --help
usage: hmtfmc compare [-h] [--generator_name1 GENERATOR_NAME1]
                      [--generator_name2 GENERATOR_NAME2]
                      input_file1 {Inel,InelGt0,V0AND} input_file2
                      {Inel,InelGt0,V0AND}

Compare the 'highlight plots' for of two estimators for two given triggers.
Requires the plots to have been previously prepared by running
`prepare_plots`.

positional arguments:
  input_file1           Path to the first file to be compared
  {Inel,InelGt0,V0AND}  Trigger of interest in first file
  input_file2           Path to the second file to be compared
  {Inel,InelGt0,V0AND}  Trigger of interest in second file

optional arguments:
  -h, --help            show this help message and exit
  --generator_name1 GENERATOR_NAME1
                        Overwrite name and tune of first generator.
  --generator_name2 GENERATOR_NAME2
                        Overwrite name and tune of second generator

$ hmtfmc compare 393_AnalysisResults.root Inel "Dipsy INEL" 393_AnalysisResults.root InelGt0 "Dipsy V0AND"
```

## Customization
Customization of the created plots can be achieved as follows:

- Percentile bins, considered triggers (INEL, V0AND, INEL>0) and estimators:
  This is centralized in the `settings.py` file. If an estimator is not available in a given `AnalysisResults.root` file, it is simply ignored.
  
- Plots that are created with `$ hmtfmc prepare_plots`:
  Customization can be found in `plotting.py` and to a small extend in the `hmtfmc` file itself.
  
- Plots that are included in the summary slides:
  This can be customized the file `summarize.py`
  
- Plots that are included in the comparison slides:
  This can be customized the file `compare.py`
  
  
