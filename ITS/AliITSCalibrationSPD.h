#ifndef ALIITSCALIBRATIONSPD_H
#define ALIITSCALIBRATIONSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "TRandom.h"
#include "AliITSCalibration.h"
#include "TArrayS.h"
#include "AliITSresponseSPD.h"


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


    // Set Threshold and noise + threshold fluctuations parameter values
    virtual  void   SetThresholds(Double_t thresh, Double_t sigma)
	{fThresh=thresh; fSigma=sigma;}
    // Get Threshold and noise + threshold fluctuations parameter values
    virtual  void   Thresholds(Double_t &thresh, Double_t &sigma) const
	{thresh=fThresh; sigma=fSigma;}
    // Set Bias Voltage parameter
    virtual void    SetBiasVoltage(Double_t bias=18.182) {fBiasVoltage=bias;}
    Double_t  GetBiasVoltage() const {return fBiasVoltage;}
    //Returns just baseline value
    Double_t GetBaseline() const {return fBaseline;}
    // Set noise and baseline in one (abstract method of AliITSCalibration)
    virtual void SetNoiseParam(Double_t n,Double_t b)
        {fNoise = n;fBaseline = b;}
    // Get noise and baseline in one (abstract method of AliITSCalibration)
    virtual void GetNoiseParam(Double_t &n,Double_t &b) const
        {n =fNoise;b = fBaseline;}
    // Returns just noise value
    Double_t GetNoise() const {return fNoise;} 
    //Declaration of member functions peculiar to this class
    // Applies a random noise and addes the baseline
    Double_t ApplyBaselineAndNoise() const {return fBaseline+
                                               fNoise*gRandom->Gaus();}
    // Set coupling parameters 
    virtual void SetCouplingParam(Double_t col, Double_t row)
        {fCouplCol = col; fCouplRow = row;}
    // Get coupling parameters 
    virtual void GetCouplingParam(Double_t &col, Double_t &row) const
        {col = fCouplCol; row = fCouplRow;}


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
    // static const Double_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const Double_t fgkThreshDefault; //default for fThresh
    static const Double_t fgkSigmaDefault; //default for fSigma
    static const Double_t fgkCouplColDefault; //default for fCouplCol
    static const Double_t fgkCouplRowDefault; //default for fCouplRow
    static const Double_t fgkBiasVoltageDefault; //default for fBiasVoltage
    Double_t fBaseline;        // Base-line value
    Double_t fNoise;           // Gaussian noise scale
    Double_t fThresh;          // Threshold value
    Double_t fSigma;           // Noise + threshold fluctuations value
    Double_t fCouplCol;        // Coupling parameter along the cols
    Double_t fCouplRow;        // Coupling parameter along the rows
    Double_t fBiasVoltage;     // Bias Voltage for the SPD (used to compute DistanceOverVoltage)
    UInt_t   fNrBad;           // Nr of SINGLE bad pixels
    TArrayS  fBadChannels;     // Array with bad channels info (col0,row0,col1...rowN) N = fNrBad
    Bool_t   fBadChip[5];     // Is chip completely dead?

    ClassDef(AliITSCalibrationSPD,7) // SPD response
};

#endif
