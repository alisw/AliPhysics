#ifndef ALI_ITS_ONLINESPDSCANINFO_H
#define ALI_ITS_ONLINESPDSCANINFO_H

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds information needed for a scan.                     //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>

class AliITSOnlineSPDscanInfo :  public TObject {

 public:
  AliITSOnlineSPDscanInfo();
  virtual ~AliITSOnlineSPDscanInfo();

  virtual UInt_t AddScanStep(); // returns the index (nsi) of the added step
  virtual void   ClearThis();
  // SET METHODS ***********************************
  void     SetType(UInt_t val) {fType=val;}
  void     SetRunNr(UInt_t val) {fRunNr=val;}
  void     SetRouterNr(UInt_t val) {fRouterNr=val;}
  void     SetTriggers(UInt_t nsi, UInt_t val);
  void     SetChipPresent(UInt_t hs, UInt_t chipi, Bool_t val)
    {fChipPresent[hs*10+chipi]=val;}
  void     SetRowStart(UInt_t val){fRowStart=val;}
  void     SetRowEnd(UInt_t val){fRowEnd=val;}
  void     SetDacStart(UInt_t val){fDacStart=val;}
  void     SetDacEnd(UInt_t val){fDacEnd=val;}  
  void     SetDacStep(UInt_t val){fDacStep=val;}

  void     IncrementTriggers(UInt_t nsi);

  // GET METHODS ***********************************
  UInt_t   GetNSteps() const {return fNSteps;}
  UInt_t   GetType() const {return fType;}
  UInt_t   GetRunNr() const {return fRunNr;}
  UInt_t   GetRouterNr() const {return fRouterNr;}
  UInt_t   GetTriggers(UInt_t nsi) const ;
  Bool_t   GetChipPresent(UInt_t hs, UInt_t chipi) const {return fChipPresent[hs*10+chipi];}
  UInt_t   GetRowStart() const {return fRowStart;}
  UInt_t   GetRowEnd() const {return fRowEnd;}
  UInt_t   GetDacStart() const {return fDacStart;}
  UInt_t   GetDacEnd() const {return fDacEnd;}
  UInt_t   GetDacStep() const {return fDacStep;}

 protected:
  UInt_t   fType;              // type of calibration scan
  UInt_t   fRunNr;             // run nr
  UInt_t   fRouterNr;          // router nr
  UInt_t   fNSteps;            // nr of s-curve steps
  TArrayI  fTriggers;          // number of triggers for the different steps of the scan
  Bool_t   fChipPresent[60];   // which chips are present
  UInt_t   fRowStart;          // row start
  UInt_t   fRowEnd;            // row end
  UInt_t   fDacStep;           // dac step
  UInt_t   fDacStart;          // dac start
  UInt_t   fDacEnd;            // dac end

  ClassDef(AliITSOnlineSPDscanInfo,1)
    };
    
#endif
