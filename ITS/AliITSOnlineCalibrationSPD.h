#ifndef ALI_ITS_ONLINECALIBRATIONSPD_H
#define ALI_ITS_ONLINECALIBRATIONSPD_H

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// This class is used as a container to keep the dead,noisy,active   //
// pixels online. Each object corresponds to one DDL (eq).           //
// Note: This class should not be used directly.                     //
// Use it via AliITSOnlineCalibrationSPDhandler instead.             //
///////////////////////////////////////////////////////////////////////

/* $Id$  */

#include <TObject.h>
#include <TArrayI.h>

class AliITSOnlineCalibrationSPD : public TObject {

 public:
    AliITSOnlineCalibrationSPD();
    virtual ~AliITSOnlineCalibrationSPD() {}

    void   SetEqNr(UInt_t eq) {fEqNr=eq;}
    UInt_t GetEqNr() const {return fEqNr;}
    void   SetBadList(TArrayI badlist) {fBadChannels=badlist;}
    void   SetNrBad(UInt_t nr) {fNrBad=nr;}

    UInt_t GetNrBad() const {return fNrBad;}
    Int_t  GetKeyAt(UInt_t index) const; //returns -1 if out of bounds

    void   ClearBad() {fBadChannels.Reset(); fNrBad=0;}

    void   ActivateALL();
    void   ActivateEq(Bool_t setval = kTRUE);
    void   ActivateHS(UInt_t hs, Bool_t setval = kTRUE);
    void   ActivateChip(UInt_t hs, UInt_t chip, Bool_t setval = kTRUE);

    Bool_t IsActiveEq() const;
    Bool_t IsActiveHS(UInt_t hs) const;
    Bool_t IsActiveChip(UInt_t hs, UInt_t chip) const;

    void   UnSetDeadALL();
    void   SetDeadEq(Bool_t setval = kTRUE);
    void   SetDeadHS(UInt_t hs, Bool_t setval = kTRUE);
    void   SetDeadChip(UInt_t hs, UInt_t chip, Bool_t setval = kTRUE);
    
    Bool_t IsDeadEq() const;
    Bool_t IsDeadHS(UInt_t hs) const;
    Bool_t IsDeadChip(UInt_t hs, UInt_t chip) const;

 private:
    UInt_t   fEqNr;                // eq nr
    UInt_t   fNrBad;               // nr of bad (single) pixels
    TArrayI  fBadChannels;         // array of keys for the bad (single) pixels
    Bool_t   fActiveEq;            // active bit for each equipment
    Bool_t   fActiveHS[6];         // active bit for each half-stave
    Bool_t   fActiveChip[60];      // active bit for each chip
    Bool_t   fDeadEq;              // dead bit for each equipment
    Bool_t   fDeadHS[6];           // dead bit for each half-stave
    Bool_t   fDeadChip[60];        // dead bit for each chip

    ClassDef(AliITSOnlineCalibrationSPD,3)
};

#endif
