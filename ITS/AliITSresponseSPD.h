#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$
*/

#include "AliITSresponse.h"
#include <TString.h>

//----------------------------------------------
//
// ITS response class for SPD
//
class AliITSresponseSPD :  public AliITSresponse {
 public:
    AliITSresponseSPD(); // default constructor
    virtual ~AliITSresponseSPD() {} // destructror
    // Configuration methods
    // sets the diffusion coeffeciant.
    virtual  void   SetDiffCoeff(Float_t p1=0) {fDiffCoeff=p1;}
    // returns the diffusion coeffeciant
    virtual  Float_t   DiffCoeff() {return fDiffCoeff;}
    // Set Threshold and noise + threshold fluctuations parameter values
    virtual  void   SetThresholds(Float_t thresh=2000, Float_t sigma=280)
	{fThresh=thresh; fSigma=sigma;}
    // Get Threshold and noise + threshold fluctuations parameter values
    virtual  void   Thresholds(Float_t &thresh, Float_t &sigma)
	{thresh=fThresh; sigma=fSigma;}
    // set coupling parameters
    virtual  void   SetNoiseParam(Float_t col=0., Float_t row=0.)
	{fCouplCol=col; fCouplRow=row;}   
    // get coupling parameters
    virtual  void   GetNoiseParam(Float_t &col, Float_t &row)
	{col=fCouplCol; row=fCouplRow;}
    // Sets the fraction of Dead SPD Pixels
    virtual void SetFractionDead(Float_t d=0.01){ fDeadPixels = d;}
    // Retruns the fraction of Dead SPD Pixels
    virtual Float_t GetFractionDead(){return fDeadPixels;}
    // Type of data - real or simulated
    virtual void    SetDataType(char *data="simulated") {fDataType=data;}
    // Get data typer
    virtual const char  *DataType() const {return fDataType.Data();}

 protected:
    Float_t fDiffCoeff;       // Sigma diffusion coefficient (not used) 
    Float_t fThresh;          // Threshold value
    Float_t fSigma;           // Noise + threshold fluctuations value
    Float_t fCouplCol;        // Coupling probability along a column
    Float_t fCouplRow;        // Coupling probability along a row
    Float_t fDeadPixels;      // the fraction of dead pixels

    TString fDataType;        // Type of data - real or simulated

    ClassDef(AliITSresponseSPD,1) // SPD response
};

#endif
