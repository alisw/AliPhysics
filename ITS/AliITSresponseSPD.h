#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$
*/

#include "AliITSresponse.h"

//----------------------------------------------
//
// ITS response class for SPD
//
class AliITSresponseSPD :  public AliITSresponse {
 public:
    AliITSresponseSPD(); // default constructor
    virtual ~AliITSresponseSPD() {} // destructror

    // Implementation of virtual member functions declared in AliITSresponse 
    // sets the diffusion coeffecient.
    virtual  void   SetDiffCoeff(Float_t p1, Float_t /*dummy */) {fDiffCoeff=p1;}
    // returns the diffusion coeffeciant
    virtual  void   DiffCoeff(Float_t &p1, Float_t & /*p2 */) const 
                             {p1 = fDiffCoeff;}
    // Set Threshold and noise + threshold fluctuations parameter values
    virtual  void   SetThresholds(Float_t thresh, Float_t sigma)
	{fThresh=thresh; fSigma=sigma;}
    // Get Threshold and noise + threshold fluctuations parameter values
    virtual  void   Thresholds(Float_t &thresh, Float_t &sigma) const
	{thresh=fThresh; sigma=fSigma;}
    // set coupling parameters
    virtual  void   SetNoiseParam(Float_t col, Float_t row)
	{fCouplCol=col; fCouplRow=row;}   
    // get coupling parameters
    virtual  void   GetNoiseParam(Float_t &col, Float_t &row) const 
	{col=fCouplCol; row=fCouplRow;}

//Declaration of member functions peculiar to this class
    // Sets the fraction of Dead SPD Pixels
    void SetFractionDead(Float_t d=0.01){ fDeadPixels = d;}
    // Retruns the fraction of Dead SPD Pixels
    Float_t GetFractionDead() const {return fDeadPixels;}

    //abstract methods in AliITSresponse not implemented in this class
    virtual void    SetDriftSpeed(Float_t /* p1 */)
      {NotImplemented("SetDriftSpeed");}
    virtual Float_t DriftSpeed() const 
      {NotImplemented("DrifSpeed"); return 0.;}
    virtual void    GiveCompressParam(Int_t *) const
      {NotImplemented("GiveCompressParam");}
    virtual void   SetElectronics(Int_t /* i */) 
                    {NotImplemented("SetElectronics");}
    virtual Int_t Electronics() const {NotImplemented("Electronics"); return 0;}
    virtual void SetMaxAdc(Float_t /* adc */) {NotImplemented("SetMaxAdc");}
    virtual Float_t MaxAdc() const {NotImplemented("MaxAdc"); return 0.;}
    virtual void    SetDynamicRange(Float_t /*dr */) 
      {NotImplemented("SetDynamicRange");}
    virtual Float_t DynamicRange() const 
      {NotImplemented("DynamicRange"); return 0.;}
    virtual void    SetChargeLoss(Float_t /* cl */)
      {NotImplemented("SetChargeLoss"); }
    virtual Float_t ChargeLoss() const 
      {NotImplemented("ChargeLoss"); return 0.;}
    virtual  void   SetDetParam(Float_t *)
      {NotImplemented("SetDetParam");}
    virtual void   GetDetParam(Float_t *) const 
      {NotImplemented("GetDetParam");}
    virtual  void   SetNDetParam(Int_t /* n */)
      {NotImplemented("SetNDetParam");}
    virtual Int_t  NDetParam() const
      {NotImplemented("NDetParam"); return 0;}
    virtual void   SetParamOptions(const char* /* a */,const char* /* b */)
      {NotImplemented("SetParamOptions");}
    virtual void   ParamOptions(char *,char*) const
      {NotImplemented("ParamOptions");} 
    virtual void   SetZeroSupp(const char*)
      {NotImplemented("SetZeroSupp");}
    virtual const char *ZeroSuppOption() const 
      {NotImplemented("ZeroSuppression"); return "";}
    virtual void    SetNSigmaIntegration(Float_t)
      {NotImplemented("SetNSigmaIntegration");}
    virtual Float_t NSigmaIntegration() const
      {NotImplemented("NSigmaIntegration"); return 0.;}
    virtual void    SetNLookUp(Int_t) 
      {NotImplemented("SetNLookUp");}
    virtual void    SetSigmaSpread(Float_t, Float_t) 
      {NotImplemented("SetSigmaSpread");}
    virtual void    SigmaSpread(Float_t & /* p1 */,Float_t & /* p2 */) const 
      {NotImplemented("SigmaSpread");}

 protected:
    static const Float_t fgkDiffCoeffDefault; //default for fDiffCoeff
    static const Float_t fgkThreshDefault; //default for fThresh
    static const Float_t fgkSigmaDefault; //default for fSigma
    Float_t fDiffCoeff;       // Sigma diffusion coefficient (not used) 
    Float_t fThresh;          // Threshold value
    Float_t fSigma;           // Noise + threshold fluctuations value
    Float_t fCouplCol;        // Coupling probability along a column
    Float_t fCouplRow;        // Coupling probability along a row
    Float_t fDeadPixels;      // the fraction of dead pixels


    ClassDef(AliITSresponseSPD,2) // SPD response
};

#endif
