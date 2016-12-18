# Simulation Files for PWGLF/FORWARD Production {#pwglf_fwd_sim}

Christian Holm Christensen
<cholm@nbi.dk>

## Introduction

The files in this archive represents a generic approach to running
AliROOT based simulations.  The idea is that the user/committer only
specifies the run number for which an anchor run is needed, and
possibly an event generator. 

## Considerations

The shell and ROOT scripts, JDLs, and so on are based on previous
"official" production files, but with some differences.  The main
difference comes about because we want to _automatically_ configure
the simulation based on the relevant GRP OCDB entry as possible.

To that end, I created the ROOT script `GRP.C` which will query the
OCDB for GRP for a specified run.  In this way, we can get the beam
species, energy, and center of mass energy simply by passing a run
number.  We can also select a default event generator based on these
settings.

* For pp, the default generator is *Pythia6*
* For PbPb, the default generator is *Hijing*
* For pPb/Pbp, the default generator is *DpmJet*

However, since all steps need to have access to the GRP data, it means
that all steering scripts must load `GRP.C`, and hence we can at no
point use any of previously used scripts. 

## Setting up a simulation

The idea of these scripts are that we only change things in one
well-defined place, and that the scripts react to JDL parameters so
that we can re-use the scripts (and settings) for many kinds of
simulation passes.

### Configuring which detectors to include

Create a class called `DetCfg` deriving from `VirtualDetCfg` and
have it return true/false for the detectors you want to have on/off.
The class _must_ be declared in `DetConfig.C` and the function
`DetConfig()` _must_ set the global object pointer `detCfg` to point
to a new instance of the class `DetCfg`.

The script `DetConfig.C` is executed by both `Simulate.C` and
`Reconstruct.C` to get the list of enabled detectors.  The script is
also executed by `AOD.C` and `QA.C` to ensure that we do not add tasks
for which we have no data because the needed detectors was turned
off. 

### Configuring OCDB specific storage locations

Create a class called `OCDBCfg` deriving from `VirtualOCDBCfg` and put
it in the script `OCDBConfig.C`.  The function `OCDBConfig` _must_ set
the global object pointer `ocdbCfg` to point to a new instance of the
class `OCDBCfg`. 

The script `OCDBConfig.C` is executed by both `Simulate.C` and
`Reconstruct.C` to set the list of specific storage locations.

One can override `VirtualOCDBCfg::Prefix()` to return the default
prefix for specific storage locations.  The member function
`VirtualOCDBCfg::Init(bool forSim)` _must_ declare all specific storage
locations.  The parameter `forSim` is true when executed from
`Simulate.C` and false when executed from `Reconstruct.C`

### Configuring the QA tasks

Create the class `QACfg` deriving from `VirtualQACfg`, and override
member functions from `VirtualQACfg` to enable/disable specific QA
tasks and options.  Put the class `QACfg` in the script `QAConfig.C`.
The function `QAConfig` _must_ set the global object pointer
`qaCfg` to point to a new instance of the class `QACfg`.

The script `QAConfig.C` is executed by `QA.C` to set which tasks and
features to include. 

### Configuring the AOD tasks

Create the class `AODCfg` deriving from `VirtualAODCfg`, and override
member functions from `VirtualAODCfg` to enable/disable specific AOD
tasks and options.  Put the class `AODCfg` in the script `AODConfig.C`.
The function `AODConfig` _must_ set the global object pointer
`aodCfg` to point to a new instance of the class `AODCfg`.

The script `AODConfig.C` is executed by `AOD.C` to set which tasks and
features to include. 


## The files

* `run.sh`: main steering executable (same as
  `/alice/bin/aliroot_new`)
* `AOD.C`: AOD Filter train set-up.  Derived from `AODtrainsim.C` but
  modified to automatically get GRP parameters from OCDB using the
  script `GRP.C` and set-up the train accordingly. Also, selection of
  which tasks to include is done through the script `AODConfig.C`.
  That means this script should never be modified by the user. 
* `AODConfig.C`: Read by `AOD.C` to configure which tasks to include
  in the train. 
* `Check.C`: Perform ESD check. Derived from `CheckESD.C`. 
* `Config.C`: Simulation configuration script.  Uses `GRP.C` to
  automatically load the proper parameters for the Anchor run.
* `DetConfig.C`: Configuration script that sets which detectors to
  turn on. This is used by all passes to ensure consistency. 
* `Final.jdl.in`: Skeleton for final merging JDL.  This is used for
  both QA and AOD filtering.
* `GRP.C`: Script that defines the global variable `grp` which is
  filled with values from the appropriate OCDB GRP entry for the
  selected run.  This is used in `AOD.C`, `Config.C`, `QA.C`,
  `Reconstruct.C`, and `Simulate.C`. 
* `LHC10hPS.C`: Special Physics Selection set-up for LHC10h due to
  missing ZDC information to the trigger.  Not needed as far as I can
  tell. 
* `Merge.jdl.in`: Skeleton for intermediate merging JDL.  This is used for
  both QA and AOD filtering.
* `OCDBConfig.C`: Set-up specific storage locations for simulation
  and reconstruction.  It defines the global object ocdbCfg which is
  used in the `Simulate.C` and `Reconstruct.C` scripts. 
* `QA.C`: QA train set-up. Derived from `QAtrainsim.C` but
  modified to automatically get GRP parameters from OCDB using the
  script `GRP.C` and set-up the train accordingly. Also, selection of
  which tasks to include is done through the script `QAConfig.C`.
  That means this script should never be modified by the user. 
* `QAConfig.C`: Read by `QA.C` to configure which tasks to include
  in the train. 
* `Reconstruct.C`: Reconstruction steering script.  Derived from
  `rec.C` but modified to automatically get GRP parameters from OCDB
  using the script `GRP.C` and set-up the job accordingly. 
* `Run.jdl.in`: Skeleton for initial stage running.
* `simrun.sh`: Simulation shell script. Derived from basic `simrun.sh`
  but adapted to the files here, allows passing the number of events
  per job, pass the run number to relevant scripts, ignores
  meaningless options. 
* `Simulate.C`: Simulation steering script.  Derived from `sim.C` but
  modified to automatically get GRP parameters from OCDB using the
  script `GRP.C` and set-up the job accordingly.
* `Tag.C`: Tag ESD files.
* `merge.sh`: Shell script to steer merging jobs.  Used by both
  QA and AOD filtering. This is based on `/alice/bin/train_merge.sh`
  but modified to call either the QA or AOD train merging.

The shell script `doit.sh` is the omnipotent script that steers all
parts of the simulation. 

## Usage of `doit.sh`

    Usage: ./doit.sh [OPTIONS]
    
    Options:
    	-h,--help		     This help
    	-t,--tag     TAG	 Job tag (pp)
    	-i,--id      NAME	 Name of production 
    	-R,--run     RUN_NO  Run number 
    	-c,--copy		     Copy files to AliEn
    	-n,--jobs    JOBS	 Set number of jobs[**] (2)
    	-m,--events  EVENTS	 Set events/job[**] (1)
    	-s,--stage   STAGE	 Set the stage[***] (0)
    	-o,--output  DIR	 Set base output directory[*] ($ALIEN_HOME/test)
    	-d,--data    DIR     Set data directory ($ALIEN_HOME/mc)
    	-b,--bin     DIR	 Set base bin directory[*] ($ALIEN_HOME/bin)
    	-a,--aliroot RELEASE Set AliROOT release [*] (v5-04-Rev-20)
    	-r,--root    RELEASE Set ROOT release [*] ()
    	-g,--geant   RELEASE Set GEANT3 release [*] ()
    
    [*] Only make sense with option -c 
    [**] Only make sense for stage 0
    [***] If stage is set to 6, try to deduce the stage
          automatically. If less than zero, do not do anything but copy

Note, if the ROOT and GEANT release is not specified, it will be
deduced from the AliROOT version.

* Make JDLs and copy files to `/alice/cern.ch/user/a/aliprod/PWGLFforwardMC`,
  specifying that output should be based in `/alice/sim/2014`, and we
  should use AliROOT v5-05-Rev-11, but not submitting any jobs

    ./doit.sh -a v5-05-Rev-11 \
		      -d /alice/cern.ch/user/a/aliprod/PWGLFforwardMC \
			  -o /alice/sim/2014 \
			  -c \
			  -s -1

* Submit 100 jobs each with 10 events for run 118506 with the same
  parameters as above, using PhoJet

    ./doit.sh -d /alice/cern.ch/user/a/aliprod/PWGLFforwardMC \
	          -i LHC14z9a -r 118506 -s 0 -n 100 -m 10 -g PhoJet

  Output will end up in `/alice/sim/2014/LHC14z9a/118506`

* Submit merging jobs for the above production

    ./doit.sh -d /alice/cern.ch/user/a/aliprod/PWGLFforwardMC \
	          -i LHC14z9a -r 118506 -s 6

  Output will end up in `/alice/sim/2014/LHC14z9a/118506/`.  Repeat
  this until we hit the final merging. 

## `Config.C`

The `Config.C` configuration script has been trimmed down
considerably.  The script is now very generic since most settings are
derived from the GRP data.  

## `EGConfig.C`

When selecting the event generator, a more versatile and configurable
approach has been taken.  `Config.C` defines the class `Setup` with
the constructor `Setup::Setup(const char* genName)`.  Here `genName`
is passed verbatim from the `Run.jdl` JDL file. The string is then
passed to method `VirtualEGCfg::MakeGenerator` which in turn calls the
virtual `VirtualEGCfg::CreateGenerator`.  This can be overwritten in a
derived class to make particular event generators.  The provided
implementation in `EGConfig.C` (loaded by `Config.C`) internally has a
switch on this string to find the chosen event generator.  If no
suitable generator can be found, the script will fail with a `Fatal`
signal.

Generators are specified as

> _name_[_sub-type_[_sub-sub-type_]]

where _sub-type_ (and _sub-sub-type_) are optional.  Currently defined
event generators are (case insensitive)

* `pythia` Pythia6 Min.Bias. Optional sub-types:
    * `perugia0` Perugia0 tune. Optional sub-sub-types:
        * `chadr` Heavy flavour charm w/hadronic decay signals added
		* `bchadr` Heavy flavour beauty/charm w/hadronic decay signals added
		* `cele` Heavy flavour charm signals added
		* `bele` Heavy flavour beauty/charm signals added
		* `jspi2e` J/Psi to electrons signals added
		* `btojspi2e` Beauty to J/Psi to electrons signals added
	* `D6T` Rick Field's CDF Tune D6T
    * `atlas` Arthur Moraes' (new) ATLAS tune, possibly with
          sub-sub-type
    	* `_flat` for flat multiplicity probability from 0 to 200
	* `jets` Jets in central barrel
* `hijing` Hijing Min.Bias. For pPb and Pbp a
  cocktail with slow neutrons is made.  Optional sub-types:
    * `2000` No quenching and hard pT cut at 2.3GeV. Possible
      sub-sub-types are:
	    * `hf` Random heavy flavour signal added (see for pythia
          above).
* `ampt` AMPT min bias. Possible sub types are:
    * `hf` Random heavy flavour signal added (see for pythia
          above).
* `dpmjet`
* `phojet` Same as `dpmjet`
* `hydjet` 

More generators can easily be added.  The idea is to have a standard
`EGConfig.C` which can be expanded upon, but if a user has very
special requirements it is possible to provide ones own `EGConfig.C`
script. 


## JDL Parameters

### `Run.jdl`

* _run number_:  The run number
* _n jobs_: Number of sub-jobs to commit
* _n events_: Number of events per sub-job
* _tag_: The tag to put in produced files under (e.g., `LHC14z9a`)
* _other_: A colon (:) separated list of options and arguments for
  `simrun.sh`. These can include
   * `bmin`:_LEAST_B_ The smallest impact parameter to make
   * `bmax`:_LARGEST_B_ The largest impact parameter to make
   * `process`:_EG_STRING_ The event generator to use.
  that is, to specify using DpmJet with an impact parameter range of 0
  to 5fm, one must pass `process:dpmjet:bmin:0:bmax:5`
 
### `Merge.jdl`

* _run number_:  The run number
* _stage_: Stage number (1-4)
* _tag_: The tag to put in produced files under (e.g., `LHC14z9a`)
* _what_: What to merge: either `AOD` or `QA`

### `Final.jdl`

* _run number_:  The run number
* _stage_: Stage number (1-4)
* _tag_: The tag to put in produced files under (e.g., `LHC14z9a`)
* _what_: What to merge: either `AOD` or `QA`

## JDL Skeleton Substitutions

* _@out@_: The base of the output directory e.g., `/alice/sim/2014/`
* _@data@_: The directory where scripts, etc. are stored e.g.,
  `/alice/cern.sh/user/a/aliprod/PWGLFforwardMC`
* _@aliroot@_: AliROOT version to use e.g., `v5-04-Rev-20`
* _@root@_: ROOT version to use e.g., `v5-34-08-6`
* _@geant@_: GEANT3 version to use e.g., `v1-15a-1`

## ROOT Script Prototypes

### `AOD.C(UInt_t run, const char* xmlfile, Int_t stage)`

* _run_: The run number
* _xmlfile_: Not used, compatibility with `QA.C`
* _stage_: Running stage (0: production, otherwise merge)

### `QA.C(UInt_t run, const char* xmlfile, Int_t stage)`

* _run_: The run number
* _xmlfile_: Current list of input files.  Must not contain `Stage`
  for final merging. 
* _stage_: Running stage (0: production, otherwise merge)

### `GRP.C(UInt_t run)`

* _run_: The run number

### `Reconstruct.C(UInt_t run)`

* _run_: The run number

### `Simulate.C(UInt_t run, Int_t nev)`

* _run_: The run number
* _nev_: Number of events per sub-job

## Shell Script Arguments

### `run.sh`

All arguments are passed on to other script e.g., `simrun.sh`

### `simrun.sh`

* _--run_ *RUN*: The Run number
* _--event_ *INTEGER*: Number of events per sub-job
* _--process_ *LABEL*: The Event Generator to use.  If not specified
  or equal to `default` then select default EG based on collision
  system. 
* _--field_ *VALUE*: Ignored - obtained from OCDB
* _--energy_ *VALUE*: Ignored - obtained from OCDB
* _--physicslist_ *VALUE*: Ignored
* _--bmin_ FERMI: Least imapact parameter
* _--bmax_ FERMI: Largest imapact parameter
* _--pthardbin_ *VALUE*: Ignored
* _--quench_ *VALUE*: Quenching type to use
* _--sdd_: Ignored
* _--number-- *INTEGER*: Sub-job number
* _--qa__: Also run QA train
* _--aod_: Also run AOD train

### `merge.sh`

* _which_: What to do merging of (either `QA` or `AOD`)
* _run_:   The run number
* _dir_:   Directory where intermediate data can be found
* _stage_: Stage to execute

Local Variables:
  mode: markdown
End:
