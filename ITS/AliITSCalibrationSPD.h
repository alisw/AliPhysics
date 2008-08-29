#ifndef ALIITSCALIBRATIONSPD_H
#define ALIITSCALIBRATIONSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "TRandom.h"
#include "AliITSCalibration.h"
#include "TArrayS.h"

////////////////////////////////////////////////////
//                                                //
// ITS response class for SPD                     //
////////////////////////////////////////////////////
class AliITSCalibrationSPD :  public AliITSCalibration {
 public:
    AliITSCalibrationSPD(); // default constructor
    virtual ~AliITSCalibrationSPD() {;} // destructror

    virtual void   ClearBad();

    virtual Int_t  GetNrBad() const;
    virtual Int_t  GetNrBadInChip(Int_t chip) const;
    virtual Int_t  GetNrBadInColumn(Int_t col) const;

    virtual Int_t  GetBadColAt(UInt_t index) const;
    virtual Int_t  GetBadRowAt(UInt_t index) const;
    virtual void   GetBadPixel(Int_t index, Int_t &row, Int_t &col) const;

    virtual Int_t  GetNrBadSingle() const {return fNrBad;}
    virtual void   SetNrBadSingle(UInt_t nr) {fNrBad=nr;} // used to be called SetNrBad, but misleading
    virtual void   SetBadList(TArrayS badlist) {fBadChannels=badlist;}
    virtual void   SetNrBad(UInt_t /*nr*/); // Use SetNrBadSingle!!!

    virtual Bool_t IsBad() const;
    virtual Bool_t IsChipBad(Int_t chip) const;
    virtual Bool_t IsColumnBad(Int_t col) const;
    virtual Bool_t IsPixelBad(Int_t col, Int_t row) const;

    virtual void   SetChipBad(UInt_t chip);
    virtual void   UnSetChipBad(UInt_t chip);

    virtual void   AddBad(UInt_t col, UInt_t row);

    virtual Int_t  GetChipIndexFromCol(UInt_t col) const;
    //    virtual Int_t  GetChipFromChipIndex(UInt_t index) const;

    virtual void    GiveCompressParam(Int_t *) const
      {NotImplemented("GiveCompressParam");}
    virtual  void   SetDetParam(Double_t *)
      {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Double_t *) const 
      {NotImplemented("GetDetParam");}
    virtual  void   SetNDetParam(Int_t /* n */)
      {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam() const
      {NotImplemented("NDetParam"); return 0;}
    virtual void    SetSigmaSpread(Double_t, Double_t) 
      {NotImplemented("SetSigmaSpread");}
    virtual void    SigmaSpread(Double_t & /* p1 */,Double_t & /* p2 */) const 
      {NotImplemented("SigmaSpread");}

 protected:
    UInt_t   fNrBad;           // Nr of SINGLE bad pixels
    TArrayS  fBadChannels;     // Array with bad channels info (col0,row0,col1...rowN) N = fNrBad
    Bool_t   fBadChip[5];      // Is chip completely dead?

    ClassDef(AliITSCalibrationSPD,8) // SPD response
};

#endif
