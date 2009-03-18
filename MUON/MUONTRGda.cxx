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

//AliRoot
#include "AliRawReaderDate.h"

#include "AliMpConstants.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONVStore.h"
#include "AliMUON1DArray.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONTriggerCrateConfig.h"

//ROOT
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TPluginManager.h"

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
    fRegionalMasks(0x0),
    fGlobalMasks(0x0),
    fTriggerIO(new AliMUONTriggerIO),
    fAlgoNoisyInput(false),
    fAlgoDeadcInput(false),
    fThrN(0.1),
    fThrD(0.9),
    fMinEvents(10),
    fSkipEvents(0),
    fMaxEvents(65535),
    fWithWarnings(false),
    fUseFastDecoder(false),
    fNLocalBoard(AliMpConstants::TotalNofLocalBoards()+1)
  {
    /// default constructor
    for (Int_t ii = 0; ii < kGlobalInputs; ii++) {
      for (Int_t il = 0; il < kGlobalInputLength; il++) {
	fAccGlobalInputN[ii][il] = 0;
	fAccGlobalInputD[ii][il] = 0;
      }
    }
    fLocalMasks    = new AliMUON1DArray(fNLocalBoard);
    fRegionalMasks = new AliMUONRegionalTriggerConfig();
    fGlobalMasks   = new AliMUONGlobalCrateConfig();
  }

  AliDAConfig (const AliDAConfig& cfg); // copy constructor
  AliDAConfig& operator=(const AliDAConfig& cfg);  // assignment operator
  virtual ~AliDAConfig()
  {
    /// destructor
    delete fLocalMasks;
    delete fRegionalMasks;
    delete fGlobalMasks; 
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

  const Char_t* GetDAConfigFileName()    { return fDAConfigFileName.Data(); }
  const Char_t* GetCurrentFileName()     { return fCurrentFileName.Data(); }
  const Char_t* GetLastCurrentFileName() { return fLastCurrentFileName.Data(); }

  const Char_t* GetSodName() { return fSodName.Data(); }
  Int_t GetSodFlag() const { return fSodFlag; }

  const Char_t* GetDAName() { return fDAName.Data(); }
  Int_t GetDAFlag() const { return fDAFlag; }

  const Char_t* GetGlobalFileName()    { return fGlobalFileName.Data(); }
  const Char_t* GetRegionalFileName()  { return fRegionalFileName.Data(); }
  const Char_t* GetLocalMaskFileName() { return fLocalMaskFileName.Data(); }
  const Char_t* GetLocalLutFileName()  { return fLocalLutFileName.Data(); }
  const Char_t* GetSignatureFileName() { return fSignatureFileName.Data(); }

  Int_t GetGlobalFileVersion()    const { return fGlobalFileVersion; }
  Int_t GetRegionalFileVersion()  const { return fRegionalFileVersion; }
  Int_t GetLocalMaskFileVersion() const { return fLocalMaskFileVersion; }
  Int_t GetLocalLutFileVersion()  const { return fLocalLutFileVersion; }
  Int_t GetSignatureFileVersion() const { return fSignatureFileVersion; }

  Int_t GetGlobalFileLastVersion()    const { return fGlobalFileLastVersion; }
  Int_t GetRegionalFileLastVersion()  const { return fRegionalFileLastVersion; }
  Int_t GetLocalMaskFileLastVersion() const { return fLocalMaskFileLastVersion; }
  Int_t GetLocalLutFileLastVersion()  const { return fLocalLutFileLastVersion; }

  AliMUONVStore*                GetLocalMasks()    const { return fLocalMasks; }
  AliMUONRegionalTriggerConfig* GetRegionalMasks() const { return fRegionalMasks; }
  AliMUONGlobalCrateConfig*     GetGlobalMasks()   const { return fGlobalMasks; }

  AliMUONTriggerIO* GetTriggerIO() const { return fTriggerIO; }

  Int_t GetEventsN() const { return fEventsN; }
  Int_t GetEventsD() const { return fEventsD; }

  Int_t GetPrintLevel() const { return fPrintLevel; }

  Bool_t GetAlgoNoisyInput() const { return fAlgoNoisyInput; }
  Bool_t GetAlgoDeadcInput() const { return fAlgoDeadcInput; }

  Float_t GetThrN() const { return fThrN; }
  Float_t GetThrD() const { return fThrD; }

  Int_t GetMinEvents()  const { return fMinEvents; }
  Int_t GetMaxEvents()  const { return fMaxEvents; }
  Int_t GetSkipEvents() const { return fSkipEvents; }

  Bool_t WithWarnings() const { return fWithWarnings; }
  Bool_t UseFastDecoder() const { return fUseFastDecoder; }

  Int_t GetGlobalInputs()      const { return kGlobalInputs; }
  Int_t GetGlobalInputLength() const { return kGlobalInputLength; }

  Int_t GetAccGlobalInputN(Int_t ii, Int_t ib) const { return fAccGlobalInputN[ii][ib]; }
  Int_t GetAccGlobalInputD(Int_t ii, Int_t ib) const { return fAccGlobalInputD[ii][ib]; }

  void SetSodName(Char_t *name) { fSodName = TString(name); }
  void SetSodFlag(Int_t flag)   { fSodFlag = flag; }

  void SetDAName(Char_t *name) { fDAName = TString(name); }
  void SetDAFlag(Int_t flag)   { fDAFlag = flag; }

  void SetGlobalFileName(const Char_t *name)    { fGlobalFileName = TString(name); }
  void SetRegionalFileName(const Char_t *name)  { fRegionalFileName = TString(name); }
  void SetLocalMaskFileName(const Char_t *name) { fLocalMaskFileName = TString(name); }
  void SetLocalLutFileName(const Char_t *name)  { fLocalLutFileName = TString(name); }
  void SetSignatureFileName(const Char_t *name) { fSignatureFileName = TString(name); }

  void SetGlobalFileVersion(Int_t ver)    { fGlobalFileVersion = ver; }
  void SetRegionalFileVersion(Int_t ver)  { fRegionalFileVersion = ver; }
  void SetLocalMaskFileVersion(Int_t ver) { fLocalMaskFileVersion = ver; }
  void SetLocalLutFileVersion(Int_t ver)  { fLocalLutFileVersion = ver; }
  void SetSignatureFileVersion(Int_t ver) { fSignatureFileVersion = ver; }

  void SetGlobalFileLastVersion(Int_t ver)    { fGlobalFileLastVersion = ver; }
  void SetRegionalFileLastVersion(Int_t ver)  { fRegionalFileLastVersion = ver; }
  void SetLocalMaskFileLastVersion(Int_t ver) { fLocalMaskFileLastVersion = ver; }
  void SetLocalLutFileLastVersion(Int_t ver)  { fLocalLutFileLastVersion = ver; }

  void IncNoiseEvent() { fEventsN++; }
  void IncDeadcEvent() { fEventsD++; }

  void AddAccGlobalInputN(Int_t ii, Int_t ib, Int_t val) { fAccGlobalInputN[ii][ib] += val; }
  void AddAccGlobalInputD(Int_t ii, Int_t ib, Int_t val) { fAccGlobalInputD[ii][ib] += val; }

  void SetPrintLevel(Int_t level) { fPrintLevel = level; }

  void SetAlgoNoisyInput(Bool_t val) { fAlgoNoisyInput = val; }
  void SetAlgoDeadcInput(Bool_t val) { fAlgoDeadcInput = val; }

  void SetThrN(Float_t val) { fThrN = val; }
  void SetThrD(Float_t val) { fThrD = val; }

  void SetMinEvents(Int_t val)  { fMinEvents = val; }
  void SetMaxEvents(Int_t val)  { fMaxEvents = val; }
  void SetSkipEvents(Int_t val) { fSkipEvents = val; }

  void SetWithWarnings() { fWithWarnings = true; }
  void SetUseFastDecoder() { fUseFastDecoder = true; }

  void IncGlobalFileVersion() { fGlobalFileVersion++; }
  void DecSkipEvents() { fSkipEvents--; }

private:

  const TString fDAConfigFileName;
  const TString fCurrentFileName;
  const TString fLastCurrentFileName;

  TString fSodName;
  Int_t   fSodFlag;

  TString fDAName;
  Int_t   fDAFlag;

  TString fGlobalFileName;
  TString fRegionalFileName;
  TString fLocalMaskFileName;
  TString fLocalLutFileName;
  TString fSignatureFileName;

  Int_t   fGlobalFileVersion;
  Int_t   fRegionalFileVersion;
  Int_t   fLocalMaskFileVersion;
  Int_t   fLocalLutFileVersion;
  Int_t   fSignatureFileVersion;
  
  Int_t   fGlobalFileLastVersion;
  Int_t   fRegionalFileLastVersion;
  Int_t   fLocalMaskFileLastVersion;
  Int_t   fLocalLutFileLastVersion;

  Int_t   fEventsN;
  Int_t   fEventsD;

  Int_t   fPrintLevel;

  AliMUONVStore*                fLocalMasks;
  AliMUONRegionalTriggerConfig* fRegionalMasks;
  AliMUONGlobalCrateConfig*     fGlobalMasks;    

  AliMUONTriggerIO *fTriggerIO;

  Bool_t fAlgoNoisyInput;
  Bool_t fAlgoDeadcInput;

  Float_t fThrN;
  Float_t fThrD;
  Int_t   fMinEvents;
  Int_t   fSkipEvents;
  Int_t   fMaxEvents;
  Bool_t  fWithWarnings;
  Bool_t  fUseFastDecoder;

  const Int_t fNLocalBoard;

  enum { kGlobalInputs = 4, kGlobalInputLength = 32 };

  Int_t fAccGlobalInputN[kGlobalInputs][kGlobalInputLength];
  Int_t fAccGlobalInputD[kGlobalInputs][kGlobalInputLength];

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
      return true;
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
    
    return true;

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
Bool_t ReadCurrentFile(AliDAConfig& cfg, TString currentFile, Bool_t lastCurrentFlag = false)
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
      return false;
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

    return true;
}

//_____________
void ReadFileNames(AliDAConfig& cfg)
{
  /// if last current file does not exist than read current file

  if (!ReadCurrentFile(cfg,cfg.GetLastCurrentFileName(), true)) 
    {
      ReadCurrentFile(cfg,cfg.GetCurrentFileName(), true);
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

    Bool_t modified = false;

    ofstream out;
    TString fileExp("ExportedFiles.dat");
    TString file;

    out.open(fileExp.Data());
    if (!out.good()) {
	printf("Failed to create file: %s\n",file.Data());
	return false;
    }      

    // check if MtgLastCurrent.dat exists
    // if not, do initial export of all files
    Bool_t initFES = false;
    if (gSystem->AccessPathName("MtgLastCurrent.dat"))
      initFES = true;
    if (initFES) printf("Copy all configuration files to the FES.\n");

    file = cfg.GetGlobalFileName();
    if ((cfg.GetGlobalFileLastVersion() != cfg.GetGlobalFileVersion()) || initFES) {
      status = daqDA_FES_storeFile(file.Data(), "GLOBAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetGlobalFileName());
	return false;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetGlobalFileName());
      out << cfg.GetGlobalFileName() << endl;
    }

    file = cfg.GetLocalMaskFileName();  
    if ((cfg.GetLocalMaskFileLastVersion() != cfg.GetLocalMaskFileVersion()) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), "LOCAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetLocalMaskFileName());
	return false;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetLocalMaskFileName());
      out << cfg.GetLocalMaskFileName() << endl;
    }

    file = cfg.GetLocalLutFileName();
    if ((cfg.GetLocalLutFileLastVersion() != cfg.GetLocalLutFileVersion()) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), "LUT");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetLocalLutFileName());
	return false;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetLocalLutFileName());
      out << cfg.GetLocalLutFileName() << endl;

    }

    // exported regional file whenever mask or/and Lut are modified
    file = cfg.GetRegionalFileName();
    if ( (cfg.GetRegionalFileLastVersion() != cfg.GetRegionalFileVersion()) || modified || initFES) {
      status = daqDA_FES_storeFile(file.Data(), "REGIONAL");
      if (status) {
	printf("Failed to export file: %s\n",cfg.GetRegionalFileName());
	return false;
      }
      if(cfg.GetPrintLevel()) printf("Export file: %s\n",cfg.GetRegionalFileName());
      out << cfg.GetRegionalFileName() << endl;
    }

    out.close();

    // export Exported file to FES anyway
    status = daqDA_FES_storeFile(fileExp.Data(), "EXPORTED");
    if (status) {
      printf("Failed to export file: %s\n", fileExp.Data());
      return false;
    }
    if(cfg.GetPrintLevel()) printf("Export file: %s\n",fileExp.Data());

    // write last current file
    WriteLastCurrentFile(cfg,cfg.GetLastCurrentFileName());

    return true;
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
      return false;
    }
 
    ReadDAConfig(cfg);
    
    status = daqDA_DB_getFile(cfg.GetCurrentFileName(), cfg.GetCurrentFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetCurrentFileName());
      return false;
    }
    
    ReadFileNames(cfg);

    status = daqDA_DB_getFile(cfg.GetGlobalFileName(), cfg.GetGlobalFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n", cfg.GetGlobalFileName());
      return false;
    }

    status = daqDA_DB_getFile(cfg.GetRegionalFileName(), cfg.GetRegionalFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetRegionalFileName());
      return false;
    }

    status = daqDA_DB_getFile(cfg.GetLocalMaskFileName(), cfg.GetLocalMaskFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetLocalMaskFileName());
      return false;
    }

    status = daqDA_DB_getFile(cfg.GetLocalLutFileName(), cfg.GetLocalLutFileName());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.GetLocalLutFileName());
      return false;
    }
 
    return true;
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
  Bool_t noise, deadc, withEvN, withEvD, updated = false;

  for (Int_t ii = 0; ii < cfg.GetGlobalInputs(); ii++) {
    gmask[ii] = 0;

    for (Int_t ib = 0; ib < cfg.GetGlobalInputLength(); ib++) {
      // lsb -> msb
      noise = false;
      deadc = false;
      withEvN = false;
      withEvD = false;
      if (cfg.GetEventsN() > cfg.GetMinEvents()) {
	rateN = (Float_t)cfg.GetAccGlobalInputN(ii,ib)/(Float_t)cfg.GetEventsN();
	noise = (rateN > cfg.GetThrN());	
	withEvN = true;
      }
      if (cfg.GetEventsD() > cfg.GetMinEvents()) {
	rateD = (Float_t)cfg.GetAccGlobalInputD(ii,ib)/(Float_t)cfg.GetEventsD();
	deadc = (rateD < cfg.GetThrD());
	withEvD = true;
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
      updated = true;
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
  
    Int_t status;
    Int_t nDateEvents = 0;

    void* event;

    // containers
    // old decoder
    AliMUONDDLTrigger*       ddlTrigger  = 0x0;
    AliMUONDarcHeader*       darcHeader  = 0x0;
    // new (fast) decoder
    const AliMUONRawStreamTriggerHP::AliHeader* darcHeaderHP = 0x0;

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
    Int_t pos = tmp.First("daq");
    tmp = tmp(pos+3,5);
    gSystem->Setenv("DATE_RUN_NUMBER",tmp.Data());
    gSystem->Exec("echo \"DATE_RUN_NUMBER = \" $DATE_RUN_NUMBER");
#endif

    if(!ExportFiles(cfg)) {
      printf("ExportFiles failed\n");
      return -1;
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
    Bool_t doUpdate = false;
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
      cfg.SetAlgoNoisyInput(false);
      cfg.SetAlgoDeadcInput(false);
      if (eventType == PHYSICS_EVENT) {
	cfg.SetAlgoNoisyInput(true);
	doUpdate = true;
	cfg.IncNoiseEvent();
      } else if (modeFET3 && eventType == CALIBRATION_EVENT) {
	cfg.SetAlgoDeadcInput(true);
	doUpdate = true;
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

	if (rawStream->GetDDL() == 0) {
	  if (cfg.UseFastDecoder()) {
	    darcHeaderHP = static_cast<AliMUONRawStreamTriggerHP*>(rawStream)->GetHeaders();
	    if (cfg.GetPrintLevel() == 2) printf("Global output (fast decoder) %x\n", (Int_t)darcHeaderHP->GetGlobalOutput());
	    for (Int_t ig = 0; ig < cfg.GetGlobalInputs(); ig++) 
	      globalInput[ig] = darcHeaderHP->GetGlobalInput(ig);
	  } else {
	    ddlTrigger = rawStream->GetDDLTrigger();
	    darcHeader = ddlTrigger->GetDarcHeader();
	    if (cfg.GetPrintLevel() == 2) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());
	    globalInput = darcHeader->GetGlobalInput();
	  }
	  StoreGlobalInput(cfg,globalInput);
	}

      } // NextDDL

      delete rawReader;
      delete rawStream;

    } // while (1)

    // update configuration files ifrequested event types were found
    if (doUpdate && cfg.GetDAFlag()) 
      UpdateGlobalMasks(cfg);

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
    cout << "MUONTRGda: Threshold for noisy inputs: " << 100*cfg.GetThrN() << "%" << endl;
    cout << "MUONTRGda: Threshold for dead inputs: " << 100*cfg.GetThrD() << "%" << endl;
    cout << "MUONTRGda: Print level: " << cfg.GetPrintLevel() << endl;
    cout << "MUONTRGda: Show decoder warnings: " << cfg.WithWarnings() << endl;
    cout << "MUONTRGda: Use the fast decoder: " << cfg.UseFastDecoder() << endl;

    printf("MUONTRGda: Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());


    return status;

}

