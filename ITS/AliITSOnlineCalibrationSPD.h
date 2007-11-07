#ifndef ALI_ITS_ONLINECALIBRATIONSPD_H
#define ALI_ITS_ONLINECALIBRATIONSPD_H

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// This class is used as a container to keep the dead and noisy      //
// pixels online. Each object corresponds to one DDL (eq).           //
// Note: This class should not be used directly.                     //
// Use it via AliITSOnlineCalibrationSPDhandler instead.             //
///////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>

class AliITSOnlineCalibrationSPD : public TObject {

 public:
    AliITSOnlineCalibrationSPD();
    virtual ~AliITSOnlineCalibrationSPD() {}

    void   SetEqNr(UInt_t mod) {fEqNr=mod;}
    UInt_t GetEqNr() const {return fEqNr;}
    void   SetBadList(TArrayI badlist) {fBadChannels=badlist;}
    void   SetNrBad(UInt_t nr) {fNrBad=nr;}

    UInt_t GetNrBad() const {return fNrBad;}
    Int_t  GetKeyAt(UInt_t index) const; //returns -1 if out of bounds

    void   ClearBad() {fBadChannels.Reset(); fNrBad=0;}

 private:
    UInt_t   fEqNr;         // eq nr
    UInt_t   fNrBad;        // nr of bad pixels
    TArrayI  fBadChannels;  // array of keys for the bad

    ClassDef(AliITSOnlineCalibrationSPD,1)
};

#endif
