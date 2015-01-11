/*
MTR DA for online

Contact: Franck Manso <manso@clermont.in2p3.fr>
Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mtrda.html
Reference run: 61898
Run Type:  PHYSICS
DA Type: MON
Number of events needed: 1000 events
Input files: MtgGlobalCrate.dat MtgRegionalCrate.dat MtgLocalMask.dat MtgLocalLut.dat MtgCurrent.dat DAConfig.dat
Output Files: ExportedFiles.dat MtgGlobalCrate.dat
Trigger types used: PHYSICS_EVENT CALIBRATION_EVENT
*/

//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for the MUON trigger configuration.                   //
//                                                                          //
// Calculates masks for the global trigger input, by looking at dead        //
// channels in calibration events and at noisy channels in physics events.  //
// Transfers trigger configuration files to the File Exchange Server and    //
// writes them (if modified) into the trigger configuration data base.      //
//                                                                          //
// Authors:                                                                 //
// Christian Fink (formerly at Subatech, Nantes)                            //
// Franck Manso (LPC Clermont-Ferrand, manso@clermont.in2p3.fr)             //
// Bogdan Vulpescu (LPC Clermont-Ferrand, vulpescu@clermont.in2p3.fr)       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//#define OFFLINE

extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

#include "AliRawReaderDate.h"

#include "AliMpConstants.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONVStore.h"
#include "AliMUON1DArray.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONLocalStruct.h"

#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TArrayS.h"

/// class for DA run parameters and DA working space
class AliDAConfig : TObject {

public:

  AliDAConfig() : 
    fDAConfigFileName("DAConfig.dat"),
    fCurrentFileName("MtgCurrent.dat"), 
    fLastCurrentFileName("MtgLastCurrent.dat"), 
    fSodName(""),
    fSodFlag(0),
    fDAName(""),
    fDAFlag(0),
    fDAMode(1),
    fGlobalFileName(""),
    fRegionalFileName(""),
    fLocalMaskFileName(""),
    fLocalLutFileName(""),
    fSignatureFileName(""),
    fGlobalFileVersion(0),
    fRegionalFileVersion(0),
    fLocalMaskFileVersion(0),
    fLocalLutFileVersion(0),
    fSignatureFileVersion(0),
    fGlobalFileLastVersion(0),
    fRegionalFileLastVersion(0),
    fLocalMaskFileLastVersion(0),
    fLocalLutFileLastVersion(0),
    fEventsN(0),
    fEventsD(0),
    fPrintLevel(0),
    fLocalMasks(0x0),
    fLocalMasksDA(0x0),
    fRegionalMasks(0x0),
    fGlobalMasks(0x0),
    fTriggerIO(new AliMUONTriggerIO),
    fAlgoNoisyInput(kFALSE),
    fAlgoDeadcInput(kFALSE),
    fThrN(0.1),
    fThrD(0.9),
    fThrLocN(0.1),
    fThrLocD(0.9),
    fMinEvents(10),
    fSkipEvents(0),
    fMaxEvents(65535),
    fWithWarnings(kFALSE),
    fUseFastDecoder(kFALSE),
    fNLocalBoard(AliMpConstants::TotalNofLocalBoards()+1),
    fPatternStoreN(0x0),
    fPatternStoreD(0x0)
  {
    /// default constructor
    for (Int_t ii = 0; ii < kGlobalInputs; ii++) {
      for (Int_t il = 0; il < kGlobalInputLength; il++) {
	fAccGlobalInputN[ii][il] = 0;
	fAccGlobalInputD[ii][il] = 0;
      }
    }
    fLocalMasks    = new AliMUON1DArray(fNLocalBoard);
    fLocalMasksDA  = new AliMUON1DArray(fNLocalBoard);
    fRegionalMasks = new AliMUONRegionalTriggerConfig();
    fGlobalMasks   = new AliMUONGlobalCrateConfig();
    fPatternStoreN = new AliMUON1DArray(fNLocalBoard);
    fPatternStoreD = new AliMUON1DArray(fNLocalBoard);

    // Generate local trigger masks store. All masks are set to FFFF
    for (Int_t i = 1; i <= AliMpConstants::TotalNofLocalBoards(); i++) {
      AliMUONVCalibParam* localBoard = new AliMUONCalibParamNI(1,8,i,0,0);
      for (Int_t x = 0; x < 2; x++) {
	for (Int_t y = 0; y < 4; y++) {
	  Int_t index = x*4+y;
	  localBoard->SetValueAsInt(index,0,0xFFFF);
	}
      }
      fLocalMasksDA->Add(localBoard);
    }
    
    // Create a default pattern store
    for (Int_t i = 1; i <= AliMpConstants::TotalNofLocalBoards(); i++) {
      AliMUONVCalibParam *patN = new AliMUONCalibParamND(2, 64, i, 0, 0.);
      AliMUONVCalibParam *patD = new AliMUONCalibParamND(2, 64, i, 0, 0.);
      fPatternStoreN->Add(patN);	
      fPatternStoreD->Add(patD);	
    }

  }

  virtual ~AliDAConfig()
  {
    /// destructor
    delete fLocalMasks;
    delete fLocalMasksDA;
    delete fRegionalMasks;
    delete fGlobalMasks; 
    delete fPatternStoreN;
    delete fPatternStoreD;
    delete fTriggerIO;
  }
  void PrintConfig()
  {
    /// print DA parameters
    printf("DA config file name: %s \n",GetDAConfigFileName());
    printf("Current file name: %s \n",GetCurrentFileName());
    printf("Last current file name: %s \n",GetLastCurrentFileName());
    printf("Global file name: %s (%d %d)\n",GetGlobalFileName(),GetGlobalFileVersion(),GetGlobalFileLastVersion());
    printf("Regional file name: %s (%d %d)\n",GetRegionalFileName(),GetRegionalFileVersion(),GetRegionalFileLastVersion());
    printf("Local mask file name: %s (%d %d)\n",GetLocalMaskFileName(),GetLocalMaskFileVersion(),GetLocalMaskFileLastVersion());
    printf("Local LUT file name: %s (%d %d)\n",GetLocalLutFileName(),GetLocalLutFileVersion(),GetLocalLutFileLastVersion());
    printf("Signature file name: %s (%d)\n",GetSignatureFileName(),GetSignatureFileVersion());
  }

     /// name of the DA configuration file from detector DB
  const Char_t* GetDAConfigFileName()    { return fDAConfigFileName.Data(); }
     /// file with current versions of the configuration files, usually MtgCurrent.dat
  const Char_t* GetCurrentFileName()     { return fCurrentFileName.Data(); }
     /// last known versions of the configuration files, usually MtgLastCurrent.dat
  const Char_t* GetLastCurrentFileName() { return fLastCurrentFileName.Data(); }

     /// name of the Start-of-data field in MtgCurrent.dat
  const Char_t* GetSodName() { return fSodName.Data(); }
     /// flag value of the Start-of-data field in MtgCurrent.dat
  Int_t GetSodFlag() const { return fSodFlag; }

     /// name of the Detector Algorithm field in MtgCurrent.dat
  const Char_t* GetDAName() { return fDAName.Data(); }
     /// flag value of the Detector Algorithm field in MtgCurrent.dat
  Int_t GetDAFlag() const { return fDAFlag; }
     /// DA active mode (GLOBAL or GLOBAL+LOCAL)
  Int_t GetDAMode() const { return fDAMode; }

     /// global crate configuration file name
  const Char_t* GetGlobalFileName()    { return fGlobalFileName.Data(); }
     /// regional crate configuration file name
  const Char_t* GetRegionalFileName()  { return fRegionalFileName.Data(); }
     /// masks for the local cards, file name
  const Char_t* GetLocalMaskFileName() { return fLocalMaskFileName.Data(); }
     /// transverse momentum Look-Up-Table, file name
  const Char_t* GetLocalLutFileName()  { return fLocalLutFileName.Data(); }
     /// signature file name
  const Char_t* GetSignatureFileName() { return fSignatureFileName.Data(); }

     /// version of the global crate configuration in the detector DB
  Int_t GetGlobalFileVersion()    const { return fGlobalFileVersion; }
     /// version of the regional crate configuration in the detector DB
  Int_t GetRegionalFileVersion()  const { return fRegionalFileVersion; }
     /// version of the masks for the local cards in the detector DB
  Int_t GetLocalMaskFileVersion() const { return fLocalMaskFileVersion; }
     /// version of the transverse momentum Look-Up-Table in the detector DB
  Int_t GetLocalLutFileVersion()  const { return fLocalLutFileVersion; }
     /// version of the signature file in the detector DB
  Int_t GetSignatureFileVersion() const { return fSignatureFileVersion; }

     /// last known version of the global crate configuration
  Int_t GetGlobalFileLastVersion()    const { return fGlobalFileLastVersion; }
     /// last known version of the regional crate configuration
  Int_t GetRegionalFileLastVersion()  const { return fRegionalFileLastVersion; }
     /// last known version of the masks for the local cards
  Int_t GetLocalMaskFileLastVersion() const { return fLocalMaskFileLastVersion; }
     /// last known version of the transverse momentum Look-Up-Table
  Int_t GetLocalLutFileLastVersion()  const { return fLocalLutFileLastVersion; }

     /// store for the masks for the local cards (own)
  AliMUONVStore*                GetLocalMasks()    const { return fLocalMasks; }
     /// store for the DA-calculated masks for the local cards (own)
  AliMUONVStore*                GetLocalMasksDA()  const { return fLocalMasksDA; }
     /// configuration object for the regional crate (own)
  AliMUONRegionalTriggerConfig* GetRegionalMasks() const { return fRegionalMasks; }
     /// configuration object for the global crate (own)
  AliMUONGlobalCrateConfig*     GetGlobalMasks()   const { return fGlobalMasks; }

     /// read/write configurations, masks and LUT to/from online files (own)
  AliMUONTriggerIO* GetTriggerIO() const { return fTriggerIO; }

     /// store for local strips patterns (noisy strips)
  AliMUONVStore* GetPatternStoreN() const { return fPatternStoreN; }
     /// store for local strips patterns (dead strips)
  AliMUONVStore* GetPatternStoreD() const { return fPatternStoreD; }

     /// number of accumulated PHYSICS events for noisy channels analysis
  Int_t GetEventsN() const { return fEventsN; }
     /// number of accumulated CALIBRATION events for dead channels analysis
  Int_t GetEventsD() const { return fEventsD; }

     /// print verbosity of the DA
  Int_t GetPrintLevel() const { return fPrintLevel; }

     /// select PHYSICS events for noisy channels analysis
  Bool_t GetAlgoNoisyInput() const { return fAlgoNoisyInput; }
     /// select CALIBRATION events for dead channels analysis
  Bool_t GetAlgoDeadcInput() const { return fAlgoDeadcInput; }

     /// threshold for noisy global inputs (fraction of events)
  Float_t GetThrN() const { return fThrN; }
     /// threshold for dead global inputs (fraction of events)
  Float_t GetThrD() const { return fThrD; }
     /// threshold for noisy local strips (fraction of events)
  Float_t GetThrLocN() const { return fThrLocN; }
     /// threshold for dead local strips (fraction of events)
  Float_t GetThrLocD() const { return fThrLocD; }

     /// minumum nr of events for rate calculation
  Int_t GetMinEvents()  const { return fMinEvents; }
     /// maximum number of events to analyze
  Int_t GetMaxEvents()  const { return fMaxEvents; }
     /// number of events to skip from start
  Int_t GetSkipEvents() const { return fSkipEvents; }

     /// show warnings from the raw data decoder
  Bool_t WithWarnings() const { return fWithWarnings; }
     /// use the high-performance (HP) decoder
  Bool_t UseFastDecoder() const { return fUseFastDecoder; }

     /// number of global input words
  Int_t GetGlobalInputs()      const { return kGlobalInputs; }
    /// length in bits of a global input word
  Int_t GetGlobalInputLength() const { return kGlobalInputLength; }

     /// get accumulated values for bit "ib" from global input word "ii", PHYSICS events
  Int_t GetAccGlobalInputN(Int_t ii, Int_t ib) const { return fAccGlobalInputN[ii][ib]; }
     /// get accumulated values for bit "ib" from global input word "ii", CALIBRATION events
  Int_t GetAccGlobalInputD(Int_t ii, Int_t ib) const { return fAccGlobalInputD[ii][ib]; }

     /// set the name of the Start-of-data field in MtgCurrent.dat
  void SetSodName(Char_t *name) { fSodName = TString(name); }
     /// set the flag value of the Start-of-data field in MtgCurrent.dat
  void SetSodFlag(Int_t flag)   { fSodFlag = flag; }

     /// set the name of the Detector Algorithm field in MtgCurrent.dat
  void SetDAName(Char_t *name) { fDAName = TString(name); }
     /// set the flag value of the Detector Algorithm field in MtgCurrent.dat
  void SetDAFlag(Int_t flag)   { fDAFlag = flag; }
     /// set DA active mode, 1 = GLOBAL (default), 2 = GLOBAL and LOCAL
  void SetDAMode(Int_t mode)   { fDAMode = mode; }

     /// set the global crate configuration file name
  void SetGlobalFileName(const Char_t *name)    { fGlobalFileName = TString(name); }
     /// set the regional crate configuration file name
  void SetRegionalFileName(const Char_t *name)  { fRegionalFileName = TString(name); }
     /// set the masks for the local cards, file name
  void SetLocalMaskFileName(const Char_t *name) { fLocalMaskFileName = TString(name); }
     /// set the transverse momentum Look-Up-Table, file name
  void SetLocalLutFileName(const Char_t *name)  { fLocalLutFileName = TString(name); }
     /// set the signature file name
  void SetSignatureFileName(const Char_t *name) { fSignatureFileName = TString(name); }

     /// set the version of the global crate configuration in the detector DB
  void SetGlobalFileVersion(Int_t ver)    { fGlobalFileVersion = ver; }
     /// set the version of the regional crate configuration in the detector DB
  void SetRegionalFileVersion(Int_t ver)  { fRegionalFileVersion = ver; }
     /// set the version of the masks for the local cards in the detector DB
  void SetLocalMaskFileVersion(Int_t ver) { fLocalMaskFileVersion = ver; }
     /// set the version of the transverse momentum Look-Up-Table in the detector DB
  void SetLocalLutFileVersion(Int_t ver)  { fLocalLutFileVersion = ver; }
     /// set the version of the signature file in the detector DB
  void SetSignatureFileVersion(Int_t ver) { fSignatureFileVersion = ver; }

     /// set the last known version of the global crate configuration
  void SetGlobalFileLastVersion(Int_t ver)    { fGlobalFileLastVersion = ver; }
     /// set the last known version of the regional crate configuration
  void SetRegionalFileLastVersion(Int_t ver)  { fRegionalFileLastVersion = ver; }
     /// set the last known version of the masks for the local cards
  void SetLocalMaskFileLastVersion(Int_t ver) { fLocalMaskFileLastVersion = ver; }
     /// set the last known version of the transverse momentum Look-Up-Table
  void SetLocalLutFileLastVersion(Int_t ver)  { fLocalLutFileLastVersion = ver; }

     /// increment the number of selected PHYSICS events
  void IncNoiseEvent() { fEventsN++; }
     /// increment the number of selected CALIBRATION events
  void IncDeadcEvent() { fEventsD++; }

     /// count the value of the bit "ib" of global input word "ii" (PHYSICS events)
  void AddAccGlobalInputN(Int_t ii, Int_t ib, Int_t val) { fAccGlobalInputN[ii][ib] += val; }
     /// count the value of the bit "ib" of global input word "ii" (CALIBRATION events)
  void AddAccGlobalInputD(Int_t ii, Int_t ib, Int_t val) { fAccGlobalInputD[ii][ib] += val; }

     /// set the print verbosity level of the DA
  void SetPrintLevel(Int_t level) { fPrintLevel = level; }

     /// select PHYSICS events for noisy channels analysis
  void SetAlgoNoisyInput(Bool_t val) { fAlgoNoisyInput = val; }
     /// select CALIBRATION events for dead channels analysis
  void SetAlgoDeadcInput(Bool_t val) { fAlgoDeadcInput = val; }

     /// set the threshold for noisy global inputs (fraction of events)
  void SetThrN(Float_t val) { fThrN = val; }
     /// set the threshold for dead global inputs (fraction of events)
  void SetThrD(Float_t val) { fThrD = val; }
     /// set the threshold for noisy local strips (fraction of events)
  void SetThrLocN(Float_t val) { fThrLocN = val; }
     /// set the threshold for dead local strips (fraction of events)
  void SetThrLocD(Float_t val) { fThrLocD = val; }

     /// set the minumum nr of events for rate calculation
  void SetMinEvents(Int_t val)  { fMinEvents = val; }
     /// set the maximum number of events to analyze
  void SetMaxEvents(Int_t val)  { fMaxEvents = val; }
     /// set the number of events to skip from start
  void SetSkipEvents(Int_t val) { fSkipEvents = val; }

     /// set/unset to show warnings from the raw data decoder
  void SetWithWarnings() { fWithWarnings = kTRUE; }
     /// set/unset the use of the high-performance (HP) decoder
  void SetUseFastDecoder() { fUseFastDecoder = kTRUE; }

    /// increment version of the global crate configuration file
  void IncGlobalFileVersion() { fGlobalFileVersion++; }
    /// increment version of the local mask configuration file
  void IncLocalMaskFileVersion() { fLocalMaskFileVersion++; }
    /// count skipped events
  void DecSkipEvents() { fSkipEvents--; }

private:

  /// copy constructor, not implemented
  AliDAConfig (const AliDAConfig& cfg);
  /// assignment operator, not implemented
  AliDAConfig& operator=(const AliDAConfig& cfg);

  const TString fDAConfigFileName;    //!< name of the DA configuration file from detector DB
  const TString fCurrentFileName;     //!< usually MtgCurrent.dat
  const TString fLastCurrentFileName; //!< usually MtgLastCurrent.dat

  TString fSodName; //!< name of the Start-of-data field in MtgCurrent.dat
  Int_t   fSodFlag; //!< flag value of the Start-of-data field in MtgCurrent.dat

  TString fDAName;  //!< name of the Detector Algorithm field in MtgCurrent.dat 
  Int_t   fDAFlag;  //!< flag value of the Detector Algorithm field in MtgCurrent.dat (enabled/disabled)
  Int_t   fDAMode;  //!< DA active mode, GLOBAL or GLOBAL+LOCAL

  TString fGlobalFileName;      //!< global crate configuration, file name
  TString fRegionalFileName;    //!< regional crate configuration, file name
  TString fLocalMaskFileName;   //!< masks for the local cards, file name
  TString fLocalLutFileName;    //!< transverse momentum Look-Up-Table, file name
  TString fSignatureFileName;   //!< signature file name

  Int_t   fGlobalFileVersion;    //!< version of the global crate configuration in the detector DB
  Int_t   fRegionalFileVersion;  //!< version of the regional crate configuration in the detector DB
  Int_t   fLocalMaskFileVersion; //!< version of the masks for the local cards in the detector DB
  Int_t   fLocalLutFileVersion;  //!< version of the transverse momentum Look-Up-Table in the detector DB
  Int_t   fSignatureFileVersion; //!< version of the signature file in the detector DB
  
  Int_t   fGlobalFileLastVersion;    //!< last known version of the global crate configuration
  Int_t   fRegionalFileLastVersion;  //!< last known version of the regional crate configuration
  Int_t   fLocalMaskFileLastVersion; //!< last known version of the masks for the local cards
  Int_t   fLocalLutFileLastVersion;  //!< last known version of the transverse momentum Look-Up-Table

  Int_t   fEventsN; //!< number of accumulated PHYSICS events
  Int_t   fEventsD; //!< number of accumulated CALIBRATION events

  Int_t   fPrintLevel;  //!< print verbosity of the DA

  AliMUONVStore*                fLocalMasks;    //!< store for the masks for the local cards
  AliMUONVStore*                fLocalMasksDA;  //!< store for the DA-calculated masks for the local cards
  AliMUONRegionalTriggerConfig* fRegionalMasks; //!< configuration object for the regional crate
  AliMUONGlobalCrateConfig*     fGlobalMasks;   //!< configuration object for the global crate

  AliMUONTriggerIO *fTriggerIO;  //!< read/write masks and LUT to/from online files

  Bool_t fAlgoNoisyInput; //!< select PHYSICS events for noisy channels analysis
  Bool_t fAlgoDeadcInput; //!< select CALIBRATION events for dead channels analysis

  Float_t fThrN;           //!< threshold for noisy global inputs (fraction of events)
  Float_t fThrD;           //!< threshold for dead global inputs (fraction of events)
  Float_t fThrLocN;        //!< threshold for noisy local strips (fraction of events)
  Float_t fThrLocD;        //!< threshold for dead local strips (fraction of events)
  Int_t   fMinEvents;      //!< minumum nr of events for rate calculation
  Int_t   fSkipEvents;     //!< number of events to skip from start
  Int_t   fMaxEvents;      //!< maximum number of events to analyze
  Bool_t  fWithWarnings;   //!< show warnings from the raw data decoder
  Bool_t  fUseFastDecoder; //!< use the high-performance (HP) decoder

  const Int_t fNLocalBoard; //!< number of local boards

  enum { kGlobalInputs = 4,         //!< number of global input words
	 kGlobalInputLength = 32    //!< length in bits of a global input word
  };

  AliMUONVStore *fPatternStoreN; //! store for local strips patterns
  AliMUONVStore *fPatternStoreD; //! store for local strips patterns

  Int_t fAccGlobalInputN[kGlobalInputs][kGlobalInputLength]; //!< storage for global input (PHYSICS events)
  Int_t fAccGlobalInputD[kGlobalInputs][kGlobalInputLength]; //!< storage for global input (CALIBRATION events)

};

//__________________________________________________________________
Bool_t ReadDAConfig(AliDAConfig& cfg)
{
    /// read run parameters for the DA

    char line[80];

    TString file;
    file = cfg.GetDAConfigFileName();
    std::ifstream in(gSystem->ExpandPathName(file.Data()));
    if (!in.good()) {
      printf("Cannot open DA configuration file %s ; use default values.\n",file.Data());
      return kTRUE;
    }

    TString tmp;
    Int_t pos;
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetPrintLevel(tmp.Atoi());
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetThrN(tmp.Atof());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetThrD(tmp.Atof());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetMinEvents(tmp.Atoi());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetSkipEvents(tmp.Atoi());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetMaxEvents(tmp.Atoi());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    if (tmp.Atoi() != 0) cfg.SetWithWarnings();
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    if (tmp.Atoi() != 0) cfg.SetUseFastDecoder();
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetThrLocN(tmp.Atof());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetThrLocD(tmp.Atof());

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.SetDAMode(tmp.Atoi());

    return kTRUE;

}

//__________________________________________________________________
void WriteLastCurrentFile(AliDAConfig& cfg, TString currentFile)
{
    /// write last current file

    ofstream out;
    TString file;
    file = currentFile;
    out.open(file.Data());
    out << cfg.GetSodName() << " " << cfg.GetSodFlag() << endl;
    out << cfg.GetDAName()  << " " << cfg.GetDAFlag()  << endl;

    out << cfg.GetGlobalFileName()    << " " << cfg.GetGlobalFileVersion()    << endl;
    out << cfg.GetRegionalFileName()  << " " << cfg.GetRegionalFileVersion()  << endl;
    out << cfg.GetLocalMaskFileName() << " " << cfg.GetLocalMaskFileVersion() << endl;
    out << cfg.GetLocalLutFileName()  << " " << cfg.GetLocalLutFileVersion()  << endl;
    out << cfg.GetSignatureFileName() << " " << cfg.GetSignatureFileVersion() << endl;

    out.close();
}

//___________________________________________________________________________________________
Bool_t ReadCurrentFile(AliDAConfig& cfg, TString currentFile, Bool_t lastCurrentFlag = kFALSE)
{
    /// read last current file name and version

    char line[80];
    char name[80];
    Int_t flag;

    TString file;
    file = currentFile;
    std::ifstream in(gSystem->ExpandPathName(file.Data()));
    if (!in.good()) {
      printf("Cannot open last current file %s\n",currentFile.Data());
      return kFALSE;
    }
    
    // read SOD 
    in.getline(line,80);  
    sscanf(line, "%s %d", name, &flag);
    cfg.SetSodName(name);
    cfg.SetSodFlag(flag);
    if (cfg.GetPrintLevel()) printf("Sod Flag %d\n", cfg.GetSodFlag());

    //read DA
    in.getline(line,80);    
    sscanf(line, "%s %d", name, &flag);
    cfg.SetDAName(name);
    cfg.SetDAFlag(flag);
    if (cfg.GetPrintLevel()) printf("DA Flag: %d\n", cfg.GetDAFlag());

    // read global
    in.getline(line,80);    
    TString tmp(line);
    Int_t pos =  tmp.First(" ");
    TString tmp1 = tmp(0, pos);
    cfg.SetGlobalFileName(tmp1.Data());
    
    if (!lastCurrentFlag) {
      cfg.SetGlobalFileVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Global File Name: %s version: %d\n", 
				      cfg.GetGlobalFileName(), cfg.GetGlobalFileVersion());
    } else {
      cfg.SetGlobalFileLastVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Global File Name: %s last version: %d\n", 
				      cfg.GetGlobalFileName(), cfg.GetGlobalFileLastVersion());
    }

    // read regional
    in.getline(line,80);
    tmp = line;
    pos = tmp.First(" ");
    tmp1 = tmp(0, pos);
    cfg.SetRegionalFileName(tmp1.Data());

    if (!lastCurrentFlag) {
      cfg.SetRegionalFileVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Regional File Name: %s version: %d\n", 
				      cfg.GetRegionalFileName(), cfg.GetRegionalFileVersion());

    } else {
      cfg.SetRegionalFileLastVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Regional File Name: %s last version: %d\n", 
				      cfg.GetRegionalFileName(), cfg.GetRegionalFileLastVersion());
    }

    // read mask
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    tmp1 = tmp(0, pos);
    cfg.SetLocalMaskFileName(tmp1.Data());

    if (!lastCurrentFlag) {
      cfg.SetLocalMaskFileVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Mask File Name: %s version: %d\n", 
				      cfg.GetLocalMaskFileName(), cfg.GetLocalMaskFileVersion());
    } else {
      cfg.SetLocalMaskFileLastVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Mask File Name: %s last version: %d\n", 
				      cfg.GetLocalMaskFileName(), cfg.GetLocalMaskFileLastVersion());
    }
    // read Lut
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    tmp1 = tmp(0, pos);
    cfg.SetLocalLutFileName(tmp1.Data());

    if (!lastCurrentFlag) {
      cfg.SetLocalLutFileVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Lut File Name: %s version: %d\n", 
				      cfg.GetLocalLutFileName(), cfg.GetLocalLutFileVersion());
    } else {
      cfg.SetLocalLutFileLastVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
      if (cfg.GetPrintLevel()) printf("Lut File Name: %s last version: %d\n", 
				      cfg.GetLocalLutFileName(), cfg.GetLocalLutFileLastVersion());
    }

    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    tmp1 = tmp(0, pos);
    cfg.SetSignatureFileName(tmp1.Data());
    cfg.SetSignatureFileVersion(atoi(tmp(pos+1, tmp.Length()-pos).Data()));
    if (cfg.GetPrintLevel()) printf("Lut File Name: %s version: %d\n", 
				    cfg.GetSignatureFileName(), cfg.GetSignatureFileVersion());

    return kTRUE;
}

//_____________
void ReadFileNames(AliDAConfig& cfg)
{
  /// if last current file does not exist than read current file

  if (!ReadCurrentFile(cfg,cfg.GetLastCurrentFileName(), kTRUE)) 
    {
      ReadCurrentFile(cfg,cfg.GetCurrentFileName(), kTRUE);
    } 
  
  // any case read current file
  ReadCurrentFile(cfg,cfg.GetCurrentFileName());
  
}

//__________________
Bool_t ExportFiles(AliDAConfig& cfg)
{
    /// Export files to FES

    // env variables have to be set (suppose by ECS ?)
    // setenv DATE_FES_PATH
    // setenv DATE_RUN_NUMBER
    // setenv DATE_ROLE_NAME
    // setenv DATE_DETECTOR_CODE

#ifdef OFFLINE
    gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");
    gSystem->Setenv("DAQDA_TEST_DIR", "/alisoft/FES");
#endif

    // update files
    Int_t status = 0;

    Bool_t modified = kFALSE;
    Bool_t globalExported = kFALSE;

    ofstream out;
    TString fileExp("ExportedFiles.dat");
    TString file;

    out.open(fileExp.Data());
    if (!out.good()) {
	printf("Failed to create file: %s\n",file.Data());
	return kFALSE;
    }      

    // check if MtgLastCurrent.dat exists
    // if not, do initial export of all files
    Bool_t initFES = kFALSE;
    if (gSystem->AccessPathName("MtgLastCurrent.dat"))
      initFES = kTRUE;
    if (initFES) printf("Copy all configuration files to the FES.\n");

    file = cfg.GetLocalMaskFileName();  
    if ((cfg.GetLocalMaskFileLastVersion() != cfg.GetLocalMaskFileVersion()) || initFES) {
      modified = kTRUE;
      status = daqDA_FES_storeFile(file.Data(), "LOCAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetLocalMaskFileName());
	return kFALSE;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetLocalMaskFileName());
      out << cfg.GetLocalMaskFileName() << endl;
    }

    file = cfg.GetLocalLutFileName();
    if ((cfg.GetLocalLutFileLastVersion() != cfg.GetLocalLutFileVersion()) || initFES) {
      modified = kTRUE;
      status = daqDA_FES_storeFile(file.Data(), "LUT");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetLocalLutFileName());
	return kFALSE;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetLocalLutFileName());
      out << cfg.GetLocalLutFileName() << endl;

    }

    file = cfg.GetGlobalFileName();
    if ((cfg.GetGlobalFileLastVersion() != cfg.GetGlobalFileVersion()) || modified || initFES) {
      modified = kTRUE;
      globalExported = kTRUE;
      status = daqDA_FES_storeFile(file.Data(), "GLOBAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetGlobalFileName());
	return kFALSE;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetGlobalFileName());
      out << cfg.GetGlobalFileName() << endl;
    }

    file = cfg.GetRegionalFileName();
    if ( (cfg.GetRegionalFileLastVersion() != cfg.GetRegionalFileVersion()) || modified || initFES) {
      status = daqDA_FES_storeFile(file.Data(), "REGIONAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetRegionalFileName());
	return kFALSE;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetRegionalFileName());
      out << cfg.GetRegionalFileName() << endl;

      // needed for the initialisation of the mapping
      if (!globalExported) {
	file = cfg.GetGlobalFileName();
	status = daqDA_FES_storeFile(file.Data(), "GLOBAL");
	if (status) {
	  printf("Failed to export file: %s\n",cfg.GetGlobalFileName());
	  return kFALSE;
	}
	if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetGlobalFileName());
	out << cfg.GetGlobalFileName() << endl;
      }

    }

    out.close();

    // export Exported file to FES anyway
    status = daqDA_FES_storeFile(fileExp.Data(), "EXPORTED");
    if (status) {
      printf("Failed to export file: %s\n", fileExp.Data());
      return kFALSE;
    }
    if(cfg.GetPrintLevel()) printf("Export file: %s\n",fileExp.Data());

    // write last current file
    WriteLastCurrentFile(cfg,cfg.GetLastCurrentFileName());

    return kTRUE;
}

//__________________
Bool_t ImportFiles(AliDAConfig& cfg)
{
    /// copy locally a file from daq detector config db 
    /// The current detector is identified by detector code in variable
    /// DATE_DETECTOR_CODE. It must be defined.
    /// If environment variable DAQDA_TEST_DIR is defined, files are copied from 
    /// DAQDA_TEST_DIR instead of the database. 
    /// The usual environment variables are not needed.

    Int_t status = 0;

#ifdef OFFLINE
    gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");
#endif

    status = daqDA_DB_getFile(cfg.GetDAConfigFileName(), cfg.GetDAConfigFileName());
    if (status) {
      printf("Failed to get DA config file from DB: %s\n",cfg.GetDAConfigFileName());
      return kFALSE;
    }
 
    ReadDAConfig(cfg);
    
    status = daqDA_DB_getFile(cfg.GetCurrentFileName(), cfg.GetCurrentFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetCurrentFileName());
      return kFALSE;
    }
    
    ReadFileNames(cfg);

    status = daqDA_DB_getFile(cfg.GetGlobalFileName(), cfg.GetGlobalFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n", cfg.GetGlobalFileName());
      return kFALSE;
    }

    status = daqDA_DB_getFile(cfg.GetRegionalFileName(), cfg.GetRegionalFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetRegionalFileName());
      return kFALSE;
    }

    status = daqDA_DB_getFile(cfg.GetLocalMaskFileName(), cfg.GetLocalMaskFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetLocalMaskFileName());
      return kFALSE;
    }

    status = daqDA_DB_getFile(cfg.GetLocalLutFileName(), cfg.GetLocalLutFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetLocalLutFileName());
      return kFALSE;
    }
 
    return kTRUE;
}

//_____________
void ReadMaskFiles(AliDAConfig& cfg)
{
    /// read mask files

    const Char_t* localFile    = cfg.GetLocalMaskFileName();
    const Char_t* regionalFile = cfg.GetRegionalFileName();
    const Char_t* globalFile   = cfg.GetGlobalFileName();

    cfg.GetTriggerIO()->ReadConfig(localFile, regionalFile, globalFile, cfg.GetLocalMasks(), cfg.GetRegionalMasks(), cfg.GetGlobalMasks());			
}

//______________________________________________________________
UInt_t GetFetMode(const AliDAConfig & cfg)
{
  /// FET mode = 3 to run algorithm for dead global inputs
  /// 0x3 prepulse
  /// 0x0 internal

  return cfg.GetGlobalMasks()->GetFetRegister(3);

}

//______________________________________________________________
void StoreGlobalInput(AliDAConfig& cfg, const UInt_t * const globalInput) 
{
  /// accumulate and build statistics of global input values
  
  for (Int_t ii = 0; ii < cfg.GetGlobalInputs(); ii++) {
    for (Int_t ib = 0; ib < cfg.GetGlobalInputLength(); ib++) {
      // lsb -> msb
      if (cfg.GetAlgoNoisyInput())
	cfg.AddAccGlobalInputN(ii,ib,((globalInput[ii] >> ib) & 0x1));
      if (cfg.GetAlgoDeadcInput())
	cfg.AddAccGlobalInputD(ii,ib,((globalInput[ii] >> ib) & 0x1));
    }
  }

}

//______________________________________________________________
void UpdateGlobalMasks(AliDAConfig& cfg) 
{
  /// update the global masks
  
#ifdef OFFLINE
  gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");
#endif

  Float_t rateN = 0.0, rateD = 0.0;
  UInt_t gmask[4], omask;
  Bool_t noise, deadc, withEvN, withEvD, updated = kFALSE;

  for (Int_t ii = 0; ii < cfg.GetGlobalInputs(); ii++) {
    gmask[ii] = 0;

    for (Int_t ib = 0; ib < cfg.GetGlobalInputLength(); ib++) {
      // lsb -> msb
      noise = kFALSE;
      deadc = kFALSE;
      withEvN = kFALSE;
      withEvD = kFALSE;
      if (cfg.GetEventsN() > cfg.GetMinEvents()) {
	rateN = (Float_t)cfg.GetAccGlobalInputN(ii,ib)/(Float_t)cfg.GetEventsN();
	noise = (rateN > cfg.GetThrN());	
	withEvN = kTRUE;
      }
      if (cfg.GetEventsD() > cfg.GetMinEvents()) {
	rateD = (Float_t)cfg.GetAccGlobalInputD(ii,ib)/(Float_t)cfg.GetEventsD();
	deadc = (rateD < cfg.GetThrD());
	withEvD = kTRUE;
      }
      if (!withEvN && !withEvD) {
	// - copy the bit from the old mask
	gmask[ii] |= ((cfg.GetGlobalMasks()->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	if (cfg.GetPrintLevel()) 
	  printf("Mask not changed (just copy the old values)\n");
      }
      if (!withEvN && withEvD) {
	if (!deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.GetGlobalMasks()->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  if (cfg.GetPrintLevel()) 
	    printf("Found dead  channel %1d:%02d (%4.2f) \n",ii,ib,rateD);
	}
      }
      if (withEvN && !withEvD) {
	if (!noise) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.GetGlobalMasks()->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  if (cfg.GetPrintLevel()) 
	    printf("Found noisy channel %1d:%02d (%4.2f) \n",ii,ib,rateN);
	}
      }
      if (withEvN && withEvD) {
	if (!noise && !deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.GetGlobalMasks()->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  if (cfg.GetPrintLevel()) {
	    if (noise)
	      printf("Found noisy channel %1d:%02d (%4.2f) \n",ii,ib,rateN);
	    if (deadc)
	      printf("Found dead  channel %1d:%02d (%4.2f) \n",ii,ib,rateD);
	  }
	}
      }
    }
  }

  // check if at least one mask value has been changed from previous version
  for (Int_t ii = 0; ii < cfg.GetGlobalInputs(); ii++) {
    printf("Global mask [%1d] %08x \n",ii,gmask[ii]);
    omask = cfg.GetGlobalMasks()->GetGlobalMask(ii);
    if (gmask[ii] != omask) {
      updated = kTRUE;
      cfg.GetGlobalMasks()->SetGlobalMask(ii,gmask[ii]);
    }
  }

  Int_t status = 0;
  if (updated) {
    
    // update version
    cfg.IncGlobalFileVersion();
    
    // don't change the file version ("-x.dat")
    
    cfg.GetTriggerIO()->WriteGlobalConfig(cfg.GetGlobalFileName(),cfg.GetGlobalMasks());
    
    // write last current file
    WriteLastCurrentFile(cfg,cfg.GetCurrentFileName());

    status = daqDA_DB_storeFile(cfg.GetGlobalFileName(), cfg.GetGlobalFileName());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.GetGlobalFileName());
      return;
    }
    
    status = daqDA_DB_storeFile(cfg.GetCurrentFileName(), cfg.GetCurrentFileName());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.GetCurrentFileName());
      return;
    }

  }
  
}

//______________________________________________________________
void UpdateLocalMask(AliDAConfig& cfg, Int_t localBoardId, Int_t connector, Int_t strip) 
{

  /// update local strip mask

  AliMUONVCalibParam* localMask = 
    static_cast<AliMUONVCalibParam*>(cfg.GetLocalMasksDA()->FindObject(localBoardId));
  
  UShort_t mask = static_cast<UShort_t>(localMask->ValueAsInt(connector,0)); 
  
  mask &= ~(0x1 << strip); // set strip mask to zero
  
  localMask->SetValueAsInt(connector, 0, mask);  

}

//______________________________________________________________
void MakePattern(AliDAConfig& cfg, Int_t localBoardId, const TArrayS& xPattern, const TArrayS& yPattern) 
{
  /// calculate the hit map for each strip in x and y direction
  
  AliMUONVCalibParam* pat = 0x0;

  if (cfg.GetAlgoNoisyInput())
    pat = static_cast<AliMUONVCalibParam*>(cfg.GetPatternStoreN()->FindObject(localBoardId));
  if (cfg.GetAlgoDeadcInput())
    pat = static_cast<AliMUONVCalibParam*>(cfg.GetPatternStoreD()->FindObject(localBoardId));

  if (!pat) {
    pat = new AliMUONCalibParamND(2, 64, localBoardId, 0, 0.);
    if (cfg.GetAlgoNoisyInput())
      cfg.GetPatternStoreN()->Add(pat);	
    if (cfg.GetAlgoDeadcInput())
      cfg.GetPatternStoreD()->Add(pat);	
  }

  for (Int_t i = 0; i < 4; ++i) {
    for (Int_t j = 0; j < 16; ++j) {
	
      Int_t xMask = xPattern[i];
      Int_t yMask = yPattern[i];
      
      Int_t index = 16*i + j;
      Double_t patOcc = 0.;
      
      if ( (xMask >> j ) & 0x1 ) {
	patOcc  = pat->ValueAsDouble(index, 0) + 1.;
	pat->SetValueAsDouble(index, 0, patOcc);
      }
      if ( (yMask >> j ) & 0x1 ) {
	patOcc  = pat->ValueAsDouble(index, 1) + 1.;
	pat->SetValueAsDouble(index, 1, patOcc);
      }
    }
  }

}

//______________________________________________________________
void MakePatternStore(AliDAConfig& cfg) 
{
  /// analyse patterns for local strips (calculate occupancy)
  
#ifdef OFFLINE
  gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");
#endif

  AliMUONVCalibParam* pat;
  Int_t localBoardId = 0;
  Int_t nEventsN = 0, nEventsD = 0;  
  UShort_t strip = 0;
  Int_t connector = 0;
  Bool_t updated = kFALSE;
  
  if (cfg.GetEventsN() > cfg.GetMinEvents()) {
    nEventsN = cfg.GetEventsN();
  }
  if (cfg.GetEventsD() > cfg.GetMinEvents()) {
    nEventsD = cfg.GetEventsD();
  }

  // noisy strips
  if (nEventsN > 0) {

    TIter next(cfg.GetPatternStoreN()->CreateIterator());

    while ( ( pat = dynamic_cast<AliMUONVCalibParam*>(next() ) ) ) {
      
      localBoardId  = pat->ID0();
      
      for (Int_t index = 0; index < pat->Size(); index++) {
	
	Double_t patXOcc  = pat->ValueAsDouble(index, 0)/(Double_t)nEventsN;
	Double_t patYOcc  = pat->ValueAsDouble(index, 1)/(Double_t)nEventsN;
	
	pat->SetValueAsDouble(index, 0, patXOcc);
	pat->SetValueAsDouble(index, 1, patYOcc);
	
	// check for x strip
	if (patXOcc > cfg.GetThrLocN()) {
	  strip  = index % 16;
	  connector = index/16;
	  UpdateLocalMask(cfg, localBoardId, connector, strip);
	}
	// check for y strip
	if (patYOcc > cfg.GetThrLocN()) {
	  strip  = index % 16;
	  connector = index/16 + 4;
	  UpdateLocalMask(cfg, localBoardId, connector, strip);
	}
	
      }
    }
    
  }
  
  // dead strips
  if (nEventsD > 0) {

    TIter next(cfg.GetPatternStoreD()->CreateIterator());

    while ( ( pat = dynamic_cast<AliMUONVCalibParam*>(next() ) ) ) {
      
      localBoardId  = pat->ID0();
      
      for (Int_t index = 0; index < pat->Size(); index++) {
	
	Double_t patXOcc  = pat->ValueAsDouble(index, 0)/(Double_t)nEventsD;
	Double_t patYOcc  = pat->ValueAsDouble(index, 1)/(Double_t)nEventsD;
	
	pat->SetValueAsDouble(index, 0, patXOcc);
	pat->SetValueAsDouble(index, 1, patYOcc);
	
	// check for x strip
	if (patXOcc < cfg.GetThrLocD()) {
	  strip  = index % 16;
	  connector = index/16;
	  UpdateLocalMask(cfg, localBoardId, connector, strip);
	}
	// check for y strip
	if (patYOcc < cfg.GetThrLocD()) {
	  strip  = index % 16;
	  connector = index/16 + 4;
	  UpdateLocalMask(cfg, localBoardId, connector, strip);
	}
	
      }
    }
    
  }
  
  // make and AND with the previous version of the mask and
  // check if the mask has changed
  UShort_t maskDA, mask;
  Int_t nMaskBits = AliMpConstants::TotalNofLocalBoards()*8*16;
  Int_t nMaskBitsChanged = 0;
  for (localBoardId = 1; localBoardId <= AliMpConstants::TotalNofLocalBoards(); localBoardId++) {
    AliMUONVCalibParam* localMaskDA = static_cast<AliMUONVCalibParam*>(cfg.GetLocalMasksDA()->FindObject(localBoardId));
    AliMUONVCalibParam* localMask = static_cast<AliMUONVCalibParam*>(cfg.GetLocalMasks()->FindObject(localBoardId));
    for (connector = 0; connector < 8; connector++) {
      maskDA = static_cast<UShort_t>(localMaskDA->ValueAsInt(connector,0)); 
      mask = static_cast<UShort_t>(localMask->ValueAsInt(connector,0));
      maskDA &= mask; 
      localMaskDA->SetValueAsInt(connector, 0, maskDA);  
      if (maskDA != mask) {
	updated = kTRUE;
	// calculated percentage of mask bits changed
	for (Int_t iBit = 0; iBit < 16; iBit++) {
	  if (((maskDA >> iBit) & 0x1) != ((mask >> iBit) &0x1)) {
	    nMaskBitsChanged++;
	  }
	}
      }
    }
  }

  printf("LOCAL mask bits changed = %5d (%7.3f %%) \n",nMaskBitsChanged,100*(Float_t)nMaskBitsChanged/(Float_t)nMaskBits);

  Int_t status = 0;
  if (updated) {

    // update version
    cfg.IncLocalMaskFileVersion();
    
    // don't change the file version ("-x.dat")
    
    cfg.GetTriggerIO()->WriteLocalMasks(cfg.GetLocalMaskFileName(),*cfg.GetLocalMasksDA());

    // write last current file
    WriteLastCurrentFile(cfg,cfg.GetCurrentFileName());

    status = daqDA_DB_storeFile(cfg.GetLocalMaskFileName(), cfg.GetLocalMaskFileName());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.GetLocalMaskFileName());
      return;
    }
    
    status = daqDA_DB_storeFile(cfg.GetCurrentFileName(), cfg.GetCurrentFileName());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.GetCurrentFileName());
      return;
    }

  }

}

//*************************************************************//
int main(Int_t argc, Char_t **argv) 
{
    /// main routine
  
    // needed for streamer application
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo", "*", "TStreamerInfo", "RIO", "TStreamerInfo()"); 

    /* check that we got some arguments = list of files */
    if (argc<2) {
      printf("Wrong number of arguments\n");
      return -1;
    }

    AliDAConfig cfg;

    Char_t inputFile[256] = "";
    inputFile[0] = 0;
    if (argc > 1)
      if (argv[1] != NULL)
        strncpy(inputFile, argv[1], 256);
      else {
        printf("MUONTRGda : No input File !\n");
        return -1;
      }

    // decoding the events
  
    Int_t status = 0;
    Int_t nDateEvents = 0;

    void* event;

    // containers
    // old decoder
    AliMUONDDLTrigger*       ddlTrigger     = 0x0;
    AliMUONDarcHeader*       darcHeader     = 0x0;
    AliMUONRegHeader*        regHeader   = 0x0;
    AliMUONLocalStruct*      localStruct = 0x0;
    // new (fast) decoder
    const AliMUONRawStreamTriggerHP::AliHeader*      darcHeaderHP  = 0x0;
    const AliMUONRawStreamTriggerHP::AliLocalStruct* localStructHP = 0x0;

    TArrayS xPattern(4);
    TArrayS yPattern(4);
    Int_t localBoardId = 0;

    TStopwatch timers;

    timers.Start(kTRUE); 

    // comment out, since we do not retrieve files from database
    if (!ImportFiles(cfg)) {
      printf("Import from DB failed\n");
      printf("For local test set DAQDA_TEST_DIR to the local directory where the Mtg files are located \n");
      return -1;
    }
    
    ReadMaskFiles(cfg);

#ifdef OFFLINE
    // the run number extracted from the file name
    TString tmp(inputFile);
    Int_t pos1 = tmp.First('d');
    Int_t pos2 = tmp.Last('.');
    Int_t len  = pos2 - (pos1+3);
    tmp = tmp(pos1+3,len);
    gSystem->Setenv("DATE_RUN_NUMBER",tmp.Data());
    gSystem->Exec("echo \"DATE_RUN_NUMBER = \" $DATE_RUN_NUMBER");
#endif

    if(!ExportFiles(cfg)) {
      printf("ExportFiles failed\n");
      return -1;
    }

    if (!cfg.GetDAFlag()) {

      cout << "MUONTRGda: DA enable: " << cfg.GetDAFlag() << endl;
      cout << "MUONTRGda: Print level: " << cfg.GetPrintLevel() << endl;
      
      printf("MUONTRGda: Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

      return status;

    }

    // FET is triggered by CTP
    Bool_t modeFET3 = kTRUE;
    if (GetFetMode(cfg) != 3) {
      printf("FET is not in mode 3. Only PHYSICS events will be analysed (noisy channels)\n");
      modeFET3 = kFALSE;
    }

    // All 5 global cards are controlled by the Mts proxy
    if (cfg.GetGlobalMasks()->GetGlobalCrateEnable() != 0x1F) {
      printf("The MTS proxy does not control all global cards\n");
      return -1;
    }

    // The global cards are ON (active on the global inputs)
    if (!cfg.GetGlobalMasks()->GetMasksOn()) {
      printf("Global masks are not ON\n");
      return -1;
    }
  
    // make sure to catch the "rare" calib events (1 every 50s in physics)
    const Char_t* tableSOD[]  = {"ALL", "yes", "CAL", "all", NULL, NULL};
    monitorDeclareTable(const_cast<char**>(tableSOD));

    status = monitorSetDataSource(inputFile);
    if (status) {
      cerr << "ERROR : monitorSetDataSource status (hex) = " << hex << status
	   << " " << monitorDecodeError(status) << endl;
      return -1;
    }
    status = monitorDeclareMp("MUON Trigger monitoring");
    if (status) {
      cerr << "ERROR : monitorDeclareMp status (hex) = " << hex << status
	   << " " << monitorDecodeError(status) << endl;
      return -1;
    }

    /* define wait event timeout - 1s max */
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);

    cout << "MUONTRGda : Reading data from file " << inputFile <<endl;

    UInt_t *globalInput = new UInt_t[4];
    Bool_t doUpdate = kFALSE;
    Int_t runNumber = 0;
    Int_t nEvents = 0;

    while(1) 
    {
      if (nEvents >= cfg.GetMaxEvents()) break;
      if (cfg.GetPrintLevel()) {
	if (nEvents && nEvents % 1000 == 0) 	
	  cout<<"Cumulated events " << nEvents << endl;
      }
      // check shutdown condition 
      if (daqDA_checkShutdown()) 
	  break;

      // Skip Events if needed
      while (cfg.GetSkipEvents()) {
	status = monitorGetEventDynamic(&event);
	cfg.DecSkipEvents();
      }

      // starts reading
      status = monitorGetEventDynamic(&event);
      if (status < 0)  {
	cout << "MUONTRGda : EOF found" << endl;
	break;
      }

      nDateEvents++;

      // decoding rawdata headers
      AliRawReader *rawReader = new AliRawReaderDate(event);
 
      Int_t eventType = rawReader->GetType();
      runNumber = rawReader->GetRunNumber();
    
      // L1Swc1
      // CALIBRATION_EVENT 
      // SYSTEM_SOFTWARE_TRIGGER_EVENT
      // DETECTOR_SOFTWARE_TRIGGER_EVENT
      cfg.SetAlgoNoisyInput(kFALSE);
      cfg.SetAlgoDeadcInput(kFALSE);
      if (eventType == PHYSICS_EVENT) {
	cfg.SetAlgoNoisyInput(kTRUE);
	doUpdate = kTRUE;
	cfg.IncNoiseEvent();
      } else if (modeFET3 && eventType == CALIBRATION_EVENT) {
	cfg.SetAlgoDeadcInput(kTRUE);
	doUpdate = kTRUE;
	cfg.IncDeadcEvent();
      } else {
	continue;
      }
      
      nEvents++;
      if (cfg.GetPrintLevel() == 2) printf("\nEvent # %d\n",nEvents);

      // decoding MUON payload
      AliMUONVRawStreamTrigger   *rawStream   = 0x0;
      if (cfg.UseFastDecoder()) {
	rawStream = new AliMUONRawStreamTriggerHP(rawReader);
      } else {
	rawStream = new AliMUONRawStreamTrigger(rawReader);
      }

      // ... without warnings from the decoder !!!
      if (!cfg.WithWarnings())
	rawStream->DisableWarnings();

      // loops over DDL 
      while((status = rawStream->NextDDL())) {

	if (cfg.GetPrintLevel() == 2) printf("iDDL %d\n", rawStream->GetDDL());

	if (cfg.UseFastDecoder()) {
	  darcHeaderHP = static_cast<AliMUONRawStreamTriggerHP*>(rawStream)->GetHeaders();
	  if (cfg.GetPrintLevel() == 2) printf("Global output (fast decoder) %x\n", (Int_t)darcHeaderHP->GetGlobalOutput());
	  for (Int_t ig = 0; ig < cfg.GetGlobalInputs(); ig++) {
	    globalInput[ig] = darcHeaderHP->GetGlobalInput(ig);
	  }
	  // loop over regional structure
	  Int_t nReg = (Int_t)static_cast<AliMUONRawStreamTriggerHP*>(rawStream)->GetRegionalHeaderCount();
	  for(Int_t iReg = 0; iReg < nReg; iReg++) {
	    // loop over local structures
	    Int_t nLoc = (Int_t)static_cast<AliMUONRawStreamTriggerHP*>(rawStream)->GetLocalStructCount(iReg);
	    for(Int_t iLoc = 0; iLoc < nLoc; iLoc++) {
	      localStructHP = static_cast<AliMUONRawStreamTriggerHP*>(rawStream)->GetLocalStruct(iReg, iLoc);
	      localBoardId = cfg.GetTriggerIO()->LocalBoardId(rawStream->GetDDL(),iReg,localStructHP->GetId());
	      if (localBoardId > 0) {
		localStructHP->GetXPattern(xPattern);
		localStructHP->GetYPattern(yPattern);
		MakePattern(cfg,localBoardId,xPattern,yPattern);
	      }
 	    }
	  }
	} else {
	  ddlTrigger = rawStream->GetDDLTrigger();
	  darcHeader = ddlTrigger->GetDarcHeader();
	  if (cfg.GetPrintLevel() == 2) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());
	  globalInput = darcHeader->GetGlobalInput();
	  // loop over regional structure
	  Int_t nReg = darcHeader->GetRegHeaderEntries();
	  for(Int_t iReg = 0; iReg < nReg; iReg++) {
	    regHeader = darcHeader->GetRegHeaderEntry(iReg);
	    // loop over local structures
	    Int_t nLoc = regHeader->GetLocalEntries();
	    for(Int_t iLoc = 0; iLoc < nLoc; iLoc++) {  
	      localStruct = regHeader->GetLocalEntry(iLoc);
	      localBoardId = cfg.GetTriggerIO()->LocalBoardId(rawStream->GetDDL(),iReg,localStruct->GetId());
	      if (localBoardId > 0) {
		localStruct->GetXPattern(xPattern);
		localStruct->GetYPattern(yPattern);
		MakePattern(cfg,localBoardId,xPattern,yPattern);
	      }
	    }
	  }
	}
	if (rawStream->GetDDL() == 0) {
	  StoreGlobalInput(cfg,globalInput);
	}
	
      } // NextDDL

      delete rawReader;
      delete rawStream;

    } // while (1)

    // update configuration files ifrequested event types were found
    if (doUpdate) {
      if (cfg.GetDAMode() > 0) UpdateGlobalMasks(cfg);
      if (cfg.GetDAMode() > 1) MakePatternStore(cfg);
    }

    timers.Stop();

    cout << "MUONTRGda: DA enable: " << cfg.GetDAFlag() << endl;
    cout << "MUONTRGda: Run number: " << runNumber << endl;
    cout << "MUONTRGda: Nb of DATE events: " << nDateEvents << endl;
    cout << "MUONTRGda: Nb of events used: " << nEvents << endl;
    cout << "MUONTRGda: Nb of events used (noise): " << cfg.GetEventsN() << endl;
    cout << "MUONTRGda: Nb of events used (deadc): " << cfg.GetEventsD() << endl;
    cout << "MUONTRGda: Minumum nr of events for rate calculation: " << cfg.GetMinEvents() << endl;
    cout << "MUONTRGda: Maximum nr of analyzed events: " << cfg.GetMaxEvents() << endl;
    cout << "MUONTRGda: Skip events from start: " << cfg.GetSkipEvents() << endl;
    cout << "MUONTRGda: Threshold for global noisy inputs: " << 100*cfg.GetThrN() << "%" << endl;
    cout << "MUONTRGda: Threshold for global dead inputs: " << 100*cfg.GetThrD() << "%" << endl;
    cout << "MUONTRGda: Threshold for local noisy inputs: " << 100*cfg.GetThrLocN() << "%" << endl;
    cout << "MUONTRGda: Threshold for local dead inputs: " << 100*cfg.GetThrLocD() << "%" << endl;
    cout << "MUONTRGda: Print level: " << cfg.GetPrintLevel() << endl;
    cout << "MUONTRGda: Show decoder warnings: " << cfg.WithWarnings() << endl;
    cout << "MUONTRGda: Use the fast decoder: " << cfg.UseFastDecoder() << endl;
    cout << "MUONTRGda: DA mode (1=GLOBAL, 2=GLOBAL+LOCAL): " << cfg.GetDAMode() << endl;

    printf("MUONTRGda: Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

    return status;

}

