#ifndef ALI_ITS_ONLINESPDPHYS_H
#define ALI_ITS_ONLINESPDPHYS_H

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online         //
// physics run.                                           //
// Directly connected to a TFile with all containers.     //
// Handles reading and writing of this TFile. Hitmaps are //
// stored in this file (AliITSOnlineSPDHitArray).         //
// Also some general information is stored                //
// (AliITSOnlineSPDphysInfo).                             //
////////////////////////////////////////////////////////////

#include <TString.h>

class TFile;
class AliITSOnlineSPDphysInfo;
class AliITSOnlineSPDHitArray;

class AliITSOnlineSPDphys {

 public:
  AliITSOnlineSPDphys():fFile(NULL),fWrite(kFALSE),fModified(kFALSE),fInfoModified(kFALSE),fPhysInfo(NULL),fFileName("."){}
  AliITSOnlineSPDphys(const Char_t *fileName, Bool_t readFromGridFile=kFALSE);
  AliITSOnlineSPDphys(const AliITSOnlineSPDphys& phys);
  virtual ~AliITSOnlineSPDphys();
  AliITSOnlineSPDphys& operator=(const AliITSOnlineSPDphys& phys);

  virtual void       AddPhys(AliITSOnlineSPDphys* phys2);
  virtual void       ClearThis();
  void               InitializeHitMap() {InitHitmap();} // online monitoring
  // SET METHODS ***********************************
  void     AddRunNr(UInt_t val);
  void     SetEqNr(UInt_t val);

  void     SetNrEvents(UInt_t val);
  void     AddNrEvents(Int_t val);
  void     IncrementNrEvents();

  void     SetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val);
  void     AddHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, Int_t val);
  void     IncrementHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  // GET METHODS ***********************************
  UInt_t   GetNrRuns() const;
  UInt_t   GetRunNr(UInt_t posi) const;
  UInt_t   GetEqNr() const;
  UInt_t   GetNrEvents() const;

  UInt_t   GetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  Float_t  GetHitsEfficiency(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  Float_t  GetHitsEfficiencyError(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  Float_t  GetAverageMultiplicity(UInt_t hs, UInt_t chipi);
  Float_t  GetAverageMultiplicityTot(UInt_t hs);

 protected:
  TFile    *fFile;                  // file to read and write from
  Bool_t   fWrite;                  // is file opened for writing?
  Bool_t   fModified;               // is the hitmap modified (needs saving)?
  Bool_t   fInfoModified;           // is the overall phys information modified (needs saving)?
  AliITSOnlineSPDphysInfo *fPhysInfo;    // overall phys information
  AliITSOnlineSPDHitArray *fHitArray[6]; // hit array, one for each halfstave
  TString  fFileName;                    // filename of file to read write

  void     InitHitmap();
  void     ReadHitmap();
  void     SaveHitmap();
  
};

#endif
