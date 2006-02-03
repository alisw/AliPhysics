#ifndef ALIITSCALIBRATIONSPD_H
#define ALIITSCALIBRATIONSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TRandom.h"
#include "AliITSCalibration.h"
#include "TArrayI.h"
#include "AliITSresponseSPD.h"


////////////////////////////////////////////////////
//                                                //
// ITS response class for SPD                     //
////////////////////////////////////////////////////
class AliITSCalibrationSPD :  public AliITSCalibration {
 public:
    AliITSCalibrationSPD(); // default constructor
    virtual ~AliITSCalibrationSPD() {;} // destructror

    // Set Threshold and noise + threshold fluctuations parameter values
    virtual  void   SetThresholds(Double_t thresh, Double_t sigma)
	{fThresh=thresh; fSigma=sigma;}
    // Get Threshold and noise + threshold fluctuations parameter values
    virtual  void   Thresholds(Double_t &thresh, Double_t &sigma) const
	{thresh=fThresh; sigma=fSigma;}
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

    //Declaration of member functions peculiar to this class
    // Sets the fraction of Dead SPD Pixels
    void SetFractionDead(Double_t d=0.01){ fDeadPixels = d;}
    // Retruns the fraction of Dead SPD Pixels
    Double_t GetFractionDead() const {return fDeadPixels;}
    // Returns a logical if a pixels is dead or not
    Bool_t IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const ;

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
 
    void SetNBadChannels(Int_t n) {fBadChannels.Set(n);}
    void AddBadChannel(Int_t bad,Int_t i) {fBadChannels.AddAt(bad,i);}
    TArrayI GetBadChannels() const {return fBadChannels;}
    virtual void SetCouplingParam(Double_t col, Double_t row) {((AliITSresponseSPD*)fResponse)->SetCouplingParam(col,row);}
    virtual void GetCouplingParam(Double_t &col, Double_t &row) const {((AliITSresponseSPD*)fResponse)->GetCouplingParam(col, row);} 
    virtual void SetParamOptions(const char* a,const char* b) {((AliITSresponseSPD*)fResponse)->SetParamOptions(a,b);}
    virtual void GetParamOptions(char *a,char* b) const {((AliITSresponseSPD*)fResponse)->ParamOptions(a,b);}

 protected:
    // static const Double_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const Double_t fgkThreshDefault; //default for fThresh
    static const Double_t fgkSigmaDefault; //default for fSigma
    Double_t fBaseline;        // Base-line value
    Double_t fNoise;           // Gaussian noise scale
    Double_t fThresh;          // Threshold value
    Double_t fSigma;           // Noise + threshold fluctuations value
    Double_t fDeadPixels;      // the fraction of dead pixels
    TArrayI fBadChannels;     // Array with bad channels info (N, row0,col0,row1,col1...rowN,colN)

    ClassDef(AliITSCalibrationSPD,2) // SPD response
};

#endif
