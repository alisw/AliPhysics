#ifndef ALI_ITS_ONLINESPDSCAN_H
#define ALI_ITS_ONLINESPDSCAN_H

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan.   //
// Directly connected to a TFile with all containers.     //
// Handles reading and writing of this TFile.             //
// Hitmaps and information on nr of events with hits      //
// is stored in this file (AliITSOnlineSPDHitArray and    //
// AliITSOnlineSPDHitEvent). Also some general            //
// information is stored (AliITSOnlineSPDscanInfo).       //
////////////////////////////////////////////////////////////

#include <Rtypes.h> 

class TFile;
class AliITSOnlineSPDscanInfo;
class AliITSOnlineSPDHitArray;
class AliITSOnlineSPDHitEvent;

class AliITSOnlineSPDscan {

 public:
  AliITSOnlineSPDscan():fFile(NULL),fWrite(kFALSE),fCurrentStep(-1),fModified(kFALSE),fInfoModified(kFALSE),fScanInfo(NULL){}
  AliITSOnlineSPDscan(Char_t *fileName);
  AliITSOnlineSPDscan(const AliITSOnlineSPDscan& scan);
  virtual ~AliITSOnlineSPDscan();
  AliITSOnlineSPDscan& operator=(const AliITSOnlineSPDscan& scan);

  virtual UInt_t     AddScanStep(); // returns the index (nsi) of the added step
  virtual void       ClearThis();       // clear all steps
  // SET METHODS ***********************************
  void     SetType(UInt_t val);
  void     SetRunNr(UInt_t val);
  void     SetRouterNr(UInt_t val);
  void     SetTriggers(UInt_t nsi, UInt_t val);
  void     SetChipPresent(UInt_t hs, UInt_t chipi, Bool_t val);
  void     SetRowStart(UInt_t val);
  void     SetRowEnd(UInt_t val);
  void     SetDacStart(UInt_t val);
  void     SetDacEnd(UInt_t val);
  void     SetDacStep(UInt_t val);

  void     SetHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val);
  void     IncrementTriggers(UInt_t nsi);
  void     IncrementHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  void     SetHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi, Int_t val);
  void     SetHitEventsTot(UInt_t nsi, UInt_t hs, Int_t val);
  void     IncrementHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi);
  void     IncrementHitEventsTot(UInt_t nsi, UInt_t hs);
  // GET METHODS ***********************************
  UInt_t   GetNSteps() const;
  UInt_t   GetType() const;
  UInt_t   GetRunNr() const;
  UInt_t   GetRouterNr() const;
  UInt_t   GetTriggers(UInt_t nsi) const;
  Bool_t   GetChipPresent(UInt_t hs, UInt_t chipi) const;
  UInt_t   GetRowStart() const;
  UInt_t   GetRowEnd() const;
  UInt_t   GetDacStart() const;
  UInt_t   GetDacEnd() const;
  UInt_t   GetDacStep() const;

  UInt_t   GetHits(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  Float_t  GetHitsEfficiency(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  Float_t  GetHitsEfficiencyError(UInt_t nsi, UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi);
  UInt_t   GetHitEvents(UInt_t nsi, UInt_t hs, UInt_t chipi);
  UInt_t   GetHitEventsTot(UInt_t nsi, UInt_t hs);
  Float_t  GetHitEventsEfficiency(UInt_t nsi, UInt_t hs, UInt_t chipi);
  Float_t  GetHitEventsTotEfficiency(UInt_t nsi, UInt_t hs);
  Float_t  GetHitEventsEfficiencyError(UInt_t nsi, UInt_t hs, UInt_t chipi);
  Float_t  GetHitEventsTotEfficiencyError(UInt_t nsi, UInt_t hs);
  Float_t  GetAverageMultiplicity(UInt_t nsi, UInt_t hs, UInt_t chipi);
  Float_t  GetAverageMultiplicityTot(UInt_t nsi, UInt_t hs);

 protected:
  TFile    *fFile;                  // file to read and write from
  Bool_t   fWrite;                  // is file opened for writing?
  Int_t    fCurrentStep;            // index of current step (kept in memory)
  Bool_t   fModified;               // is the current step modified (needs saving)?
  Bool_t   fInfoModified;           // is the overall scan information modified (needs saving)?
  AliITSOnlineSPDscanInfo *fScanInfo;           // overall scan information
  AliITSOnlineSPDHitArray *fCurrentHitArray[6]; // hit array, one for each halfstave
  AliITSOnlineSPDHitEvent *fCurrentHitEvent[6]; // hit events, one for each halfstave
  Char_t   fFileName[200];                      // filename of file to read write

  void     Init();
  void     CreateNewStep();
  void     SwitchToStep(UInt_t nsi);
  void     FillGap(UInt_t nsi);
  void     ReadCurrentStep();
  void     SaveCurrentStep();
  
};

#endif
