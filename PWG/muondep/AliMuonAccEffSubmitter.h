#ifndef ALIMUONACCEFFSUBMITTER_H
#define ALIMUONACCEFFSUBMITTER_H

/**

  \ingroup pwg_muondep_submitter
  \class AliMuonAccEffSubmitter
  \brief Ease the handling of private muon-related MC productions. 
  
@todo Note that this class was developed prior to the existence of AliDPG. One day we may consider (or not ;-) depending on
how doable that is) changing this class to make things the DPG way. 

The class is dealing with three different kind of directories (template, local, remote), which are explained in the base class AliMuonGridSubmitter.

Various usages are described below :

- \ref basicusage
- \ref moreadvancedusage
- \ref pseudoidealsimulations

# Basic usage {#basicusage}

Let's assume you want to make simulations to compute the Acc x Eff of J/psi for pp at 13 TeV, period LHC16p.

You first locate in the template directory the generator name to be used, in this case 
`GenJPsi13TeV`

You then create an object of the AliMuonAccEffSubmitter class : 

```
AliMuonAccEffSubmitter submitter("GenJPsi13TeV");
```

This will compile the `GenJPsi13TeV.C` macro.

Specify the remote directory to be used. The directory will be created if it does not exist yet. 

```
submitter.SetRemoteDir("/alice/cern.ch/user/l/laphecet/simulations/tutorial/lhc16p");
```

Ask for as many events in the simulation as in real CMUL7 triggers (or any multiplier you see fit, e.g. 1.1 or 10, depending on the output statistics you want, depending on the real data statistics)

```
submitter.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUFAST",1);
```

Specify the run to be used as anchor : 

```
submitter.SetRunList(264197); 
```

To ease testing, we allow to overwrite output files (normally set to false to protect the innocent) : 

```
submitter.ShouldOverwriteFiles(true);
```

At this point the object is configured. It can be inspected : 

```
submitter.Print();
```

The Run method can then be used to perform some operations. 
If you are feeling lucky, try the full package to start with :

```
submitter.Run("full");
```

If you want to be more cautious, try first to generate the local files so you can inspect them :

```
submitter.Run("local");
```

You can also try a local test before going to the full-blown production on the Grid

```
submitter.MakeNofEventsFixed(10);
submitter.Run("localtest");
submitter.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUFAST",1); // important in case you do the submission in the same Root session afterwards...
```

At the end of the local test you can inspect the generated `*.log` files to check everything is as you expected.

Alternatively, you can use a two-stage operation, where you first do everything but the submit (you can then check the local files and also see how many jobs would be submitted, depending on your chosen number of triggers and maximum number of events per job) and then the submission itself.

```
submitter.Run("test");
submitter.Submit(false);
```

The submission will give you the masterjob identifier for each run as well as the number of events that will be generated.

```
run     firstChunk      lastChunk       eventsPerJob
----------------------
264197  1       43      4714
submit JDL 264197 1 43 4714 ...DONE
  --> the job Id is: 801549938

total number of jobs = 43
total number of generated events = 202702
```

# More advanced usage {#moreadvancedusage}

There are a number of options that can be changed to tune the behavior of the submitter. 

In particular, you can alter the runlist using the AliMuonGridSubmitter::SetRunList methods, the maximum number of events per job (\ref SetMaxEventsPerChunk), the AliPhysics version to be used with AliMuonGridSubmitter::SetAliPhysicsVersion, the kind of output you want to keep with \ref SetCompactMode. 

Variables used in the template files are changed with AliMuonGridSubmitter::SetVar (and displayed with either PrintVariables or Print or GetVar).

# Specific use for pseudo-ideal simulations {#pseudoidealsimulations}

In order to compute a quick acc x eff there is the possibility to generate _pseudo-ideal_ simulations. By pseudo-ideal simulation we mean : 

- use ideal pedestals (mean 0, sigma 1)
- complete configuration (i.e. all manus are there)
- raw/full alignment
- do _not_ cut on the pad status, i.e. disregard occupancy, HV, LV or RejectList  completely
 
For instance, making a simulation of 800k J/psi for pp 13 TeV :

```
AliMuonAccEffSubmitter a("GenJPsi13TeV",kFALSE,"4_25",800000,4000);
a->SetRemoteDir("/alice/cern.ch/user/l/laphecet/IdealSimulations5");
a->SetRunList(228936);
a->SetAliPhysicsVersion("VO_ALICE@AliPhysics::v5-08-13o-01-1");
a->Run("full");
```

Note that this will create (if not already done otherwise) OCDB snapshots in the local directory `./OCDB/228936`. 
Those snapshots are automatically generated (\ref AliMuonOCDBSnapshotGenerator) from a local OCDB create by the submitter in `./OCDB` with default objects for Config, Pedestals, OccupancyMap, HV, LV and RejectList. 
The other objects are taken directly from the `raw://` OCDB (or any OCDB the submitter is configured with), except the RecoParam which are patched to remove the cut on any of the defaults objects in the local OCDB (the idea being to produce an ideal but realistic simulation, as far as geometry, alignment, tracking and trigger parameters are concerned).

Those simulations are then to be dealt with using the _compact_ family of classes (@ref pwg_muondep_compact) 

\author Laurent Aphecetche, Subatech
 
*/

#include "AliMuonGridSubmitter.h"

class AliMuonAccEffSubmitter : public AliMuonGridSubmitter
{
public:

  AliMuonAccEffSubmitter(Bool_t localOnly=kFALSE);
  AliMuonAccEffSubmitter(const char* generator,
                         Bool_t localOnly=kFALSE,
                         const char* generatorVersion="8125");
  AliMuonAccEffSubmitter(const char* generator,
          Bool_t localOnly,
          const char* pythia6version,
          Int_t numberOfEventsForPseudoIdealSimulation,
          Int_t maxEventsPerChunk);

  virtual Bool_t Generate(const char* jdlname) const;
  virtual Bool_t Run(const char* mode);

  virtual ~AliMuonAccEffSubmitter();

  void MakeNofEventsPropToTriggerCount(const char* trigger="CMUL7-B-NOPF-MUON", Float_t ratio=1.0);
  
  void MakeNofEventsFixed(Int_t nevents);
 
  Bool_t UseOCDBSnapshots() const;

  void UseOCDBSnapshots(Bool_t flag);
  
  void UseExternalConfig(const char* externalConfigFullFilePath);
  
  void UseAODMerging(Bool_t flag);

  Bool_t SetupCollision ( Double_t cmsEnergy, Int_t lhapdf = 10800, const char *nucleons = "pp", const char *collSystem = "pp", Int_t npdf = 0, Int_t npdfErr = 1 );
  Bool_t SetupPowheg ( const char *particle, const char *version = "r3178-5" );
  void SetupPythia6 ( const char *version );
  void SetupPythia8 ( const char *version, const char* configStrings = "" );
  
  Bool_t Merge(Int_t stage, Bool_t dryRun=kTRUE);

  Int_t Submit(Bool_t dryRun=kTRUE);
  
  Int_t LocalTest();
  
  /// Return the name of the JDL file that will be produced
  TString RunJDLName() const { return "JDL"; }

  /// Return the name of the JDL file that will be produced for the merging phase
  TString MergeJDLName(Bool_t final) const { return (final ? "AOD_merge_final.jdl" : "AOD_merge.jdl"); }

  virtual void Print(Option_t* opt="") const;

  /// Set the input splitting level
  void SetSplitMaxInputFileNumber(Int_t n) { fSplitMaxInputFileNumber = n; }
  
  /// Get the input splitting level
  Int_t GetSplitMaxInputFileNumber() const { return fSplitMaxInputFileNumber; }
  
  void SetCompactMode(Int_t mode);
  
  /// Set the output files to be kept
  void SetCustomOutFiles(const char* logFiles, const char* rootFiles) { fLogOutToKeep = logFiles; fRootOutToKeep = rootFiles; }
  
  Bool_t MakeOCDBSnapshots();
  
  void SetOCDBPath(const char* ocdbPath);

  void SetOCDBSnapshotDir(const char* dir);

  Bool_t SetGenerator(const char* generator);

  /// Get the max number of events to be simulated by job
  Int_t MaxEventsPerChunk() const { return fMaxEventsPerChunk; }

  /// Set the max number of events to be simulated by job
  void SetMaxEventsPerChunk(Int_t n) { fMaxEventsPerChunk = n; SetVar("VAR_EVENTS_PER_JOB", Form("%i",n)); }

  /// Get the OCDB path used
  TString OCDBPath() const { return GetMapValue("OCDBPath"); }
  
  /// Get the directory of the remote directory used for OCDB snapshots
  TString RemoteSnapshotDir() const { return GetMapValue("OCDBsnapshot"); }
  /// Get the directory of the local directory used for OCDB snapshots
  TString LocalSnapshotDir() const { return GetMapValue("Snapshot"); }

  Int_t SplitRunList(const char* inputList, int maxJobs=1500);
  
private:

  Bool_t GenerateRunJDL(const char* name) const;
  
  Bool_t GenerateMergeJDL(const char* name) const;

  void SetDefaultVariables();

  void SetLocalOnly();

  void SetupCommon(Bool_t localOnly);

  /// Get the reference trigger (used to get number of events per job)
  TString ReferenceTrigger() const { return GetMapValue("ReferenceTrigger"); }
  
  void UpdateLocalFileList(Bool_t clearSnapshots=kFALSE);

private:
  /// not implemented
  AliMuonAccEffSubmitter(const AliMuonAccEffSubmitter& rhs);
  /// not implemented
  AliMuonAccEffSubmitter& operator=(const AliMuonAccEffSubmitter& rhs); //< not implemented
  
private:
  Float_t fRatio; ///< ratio simulated events vs real events
  Int_t fFixedNofEvents; ///< fixed number of events to be used per run
  Int_t fMaxEventsPerChunk; ///< max events to generate per subjob
  Int_t fSplitMaxInputFileNumber; ///< used for merging jdl
  TString fLogOutToKeep; ///< specify the log files to be kept
  TString fRootOutToKeep; ///< specify the output files to be kept
  TString fExternalConfig; ///< path to an (optional) external config file
  TString fSnapshotDir; ///< directory for OCDB snapshots
  Bool_t fUseAODMerging; ///< whether or not to perform (aod) merging
  
/// \cond CLASSIMP
  ClassDef(AliMuonAccEffSubmitter,5); // Helper class to submit AccxEff single particle simulations
/// \endcond
};

#endif

