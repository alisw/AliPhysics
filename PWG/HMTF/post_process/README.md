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

	$ hmtfmc --help
	usage: hmtfmc [-h] [-v] {compare,prepare_plots,summary} ...

	positional arguments:
	{compare,prepare_plots,summary}

	optional arguments:
	-h, --help            show this help message and exit
	-v, --verbose

## Preparing plot for later steps

This steps produces all plots and saves them in the downloaded AnalysisResults.root file. This step needs to be taken before creating the summary or comparison slides. For example:

	$ mkdir ~/hmtfmc_results
	$ cd ~/hmtfmc_results
	$ hmtfmc prepare_plots /alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/393_20160210-1718/merge/AnalysisResults.root
	
Now, there is a file called 393_AnalysisResults.root (393 being the train number). The file contains all the plots at this points

## Making the summary slides for one generator

Assuming that the previously downloaded file was for "Dipsy 13 TeV Ropes", the following will create the summary slides. Again, `--help` is available.

``` shell
$ hmtfmc summary --help
	usage: hmtfmc summary [-h] input_file gen_name

Create summary slides. Configurations can be made via a custom json file

positional arguments:
  input_file  Path to file containing the input data
  gen_name    Name of the generator used. Used in the title of the PDF

optional arguments:
  -h, --help  show this help message and exit
  
$ hmtfmc summary 393_AnalysisResults.root "Dipsy 13TeV Ropes"
```

the summary slides can now be found in the aptly called folders created in the current directory

## Making comparison slides
Two generators and/or tuning can be compared as follows

```shell
$ hmtfmc compare --help
usage: hmtfmc compare [-h]
                      input_file1 {Inel,InelGt0,V0AND} generator_name1
                      input_file2 {Inel,InelGt0,V0AND} generator_name2

Compare the 'highlight plots' for of two estimators for two given triggers.
Requires the plots to have been previously prepared by running
`prepare_plots`.

positional arguments:
  input_file1           Path to the first file to be compared
  {Inel,InelGt0,V0AND}  Trigger of interest in first file
  generator_name1       Name and tune of first generator
  input_file2           Path to the second file to be compared
  {Inel,InelGt0,V0AND}  Trigger of interest in second file
  generator_name2       Name and tune of second generator

optional arguments:
  -h, --help            show this help message and exit
$ hmtfmc compare 393_AnalysisResults.root Inel "Dipsy INEL" 393_AnalysisResults.root InelGt0 "Dipsy INEL>0"
```

