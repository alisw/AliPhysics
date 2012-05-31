#ifndef ALIITSONLINESPDFO_H
#define ALIITSONLINESPDFO_H  
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                     // 
// This class is used within the detector algorithm framework //
// to write and read FO scan data.                            //
////////////////////////////////////////////////////////////////

#include <TString.h>
#include <THashList.h>

class TFile;
class TArrayS;
class AliITSOnlineSPDfoInfo;

class AliITSOnlineSPDfo {

 public:
  AliITSOnlineSPDfo();//ctor
  AliITSOnlineSPDfo(TString inputfile, Int_t runNr, Int_t eqId); 
  AliITSOnlineSPDfo(const AliITSOnlineSPDfo &c);
  
  virtual ~AliITSOnlineSPDfo(){ delete fArray; delete fDACnames;} //dctor
 
  enum {kFOPOL=0, kCONVPOL=1, kCOMPREF=2, kCGPOL =3, kPreVTH=4, kIdFOPOL=20, kIdCONVPOL=17, kIdCOMPREF=16, kIdCGPOL=14, kIdPreVTH=39};
  
  // GENERAL METHODS
  void    AddMeasurement(const TArrayS dac, Short_t measure[3], Int_t hs, Int_t chipId);
  Int_t   CheckDACEntry(const TArrayS dac);
  TString CreateDACEntry(const TArrayS dacs) const;
  TArrayS CreateDACArray(const TArrayS dacs, const TArrayS dacId) const;
  void    CreateOutputFile(); 
  void    WriteToFile();
  
  // SETTER
 void SetNdacs(UInt_t ndacs) {fNdacs=ndacs;}
 void SetFOscanParams(AliITSOnlineSPDfoInfo *info) {fInfo=info;}
 void SetFile(TString inputfile);
 void SetDACArray(TObjArray *obj) {if(!fArray) fArray = obj; else printf("The fArray is alreay available, no need to set it again.\n");}
 
 // GETTER 
  TFile* GetFile()   const           {return fFile;}
  UInt_t GetNdacs()  const           {return fNdacs;}
  
  TObjArray *GetDACArray() const     {return fArray;}
  THashList * GetDACnameList() const {return fDACnames;}
  
  TArrayI GetDACscanParams() const;                               // retrieves per each DAC the range and the step used in the scan 
  AliITSOnlineSPDfoInfo * GetFOscanInfo() const {return fInfo;}
  Int_t *GetDACvalues(TString s, const Int_t ndacs) const;        // translates the string in the corresponding DACS                                                      // the user has to delete the pointer after use!
  Double_t *GetDACvaluesD(TString s, const Int_t ndacs) const;    // translates the string in the corresponding DACS
  

 protected:
 Int_t fRunNr;
 Int_t fNdacs;
 TString fFileName;
 TFile *fFile; 
 AliITSOnlineSPDfoInfo *fInfo;
 THashList *fDACnames;
 TObjArray *fArray;               // array of the 10 chips in the 6 HS per DAC set
 Double_t fCheckIndex;            // check array index (to speed up)
 Int_t fIndex; 
 TString fInitialConfiguration;
  
  private:
  AliITSOnlineSPDfo& operator= (const AliITSOnlineSPDfo& c);
  
    ClassDef(AliITSOnlineSPDfo,1)
  };
    
#endif
