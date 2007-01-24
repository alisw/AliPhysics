#ifndef ALI_ITS_ONLINECALIBRATIONSPD_H
#define ALI_ITS_ONLINECALIBRATIONSPD_H

///////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                            //
// This class is used as a container to keep the dead and noisy      //
// pixels online. Each object corresponds to one module.             //
// Note: This class should not be used directly.                     //
// Use it via AliITSOnlineCalibrationSPDhandler instead.             //
///////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>

class AliITSOnlineCalibrationSPD : public TObject {

 public:
    AliITSOnlineCalibrationSPD();
    virtual ~AliITSOnlineCalibrationSPD() {}

    void   SetModuleNr(UInt_t mod) {fModuleNr=mod;}
    UInt_t GetModuleNr() const {return fModuleNr;}
    void   SetDeadList(TArrayI deadlist) {fDeadChannels=deadlist;}
    void   SetNoisyList(TArrayI noisylist) {fNoisyChannels=noisylist;}
    void   SetNrDead(UInt_t nr) {fNrDead=nr;}
    void   SetNrNoisy(UInt_t nr) {fNrNoisy=nr;}

    void   AddDead(UInt_t col, UInt_t row);
    Int_t  GetNrDead() const {return fNrDead;}
    Int_t  GetDeadColAt(UInt_t index) const; //returns -1 if out of bounds
    Int_t  GetDeadRowAt(UInt_t index) const; //returns -1 if out of bounds
    void   ClearDead() {fDeadChannels.Reset(); fNrDead=0;}
    Bool_t IsPixelDead(Int_t col, Int_t row) const ;
    void   AddNoisy(UInt_t col, UInt_t row);
    Int_t  GetNrNoisy() const {return fNrNoisy;}
    Int_t  GetNoisyColAt(UInt_t index) const; //returns -1 if out of bounds
    Int_t  GetNoisyRowAt(UInt_t index) const; //returns -1 if out of bounds
    void   ClearNoisy() {fNoisyChannels.Reset(); fNrNoisy=0;}
    Bool_t IsPixelNoisy(Int_t col, Int_t row) const ;

 private:
    UInt_t   fModuleNr;       // module nr
    UInt_t   fNrDead;         // nr of dead pixels
    TArrayI  fDeadChannels;   // index 0 is nr of dead, then pairs of col and row
    UInt_t   fNrNoisy;        // nr of noisy pixels
    TArrayI  fNoisyChannels;  // index 0 is nr of noisy, then pairs of col and row

    ClassDef(AliITSOnlineCalibrationSPD,1)
};

#endif
