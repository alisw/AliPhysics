#ifndef ALIVAODHEADER_H
#define ALIVAODHEADER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Virtual event header class
//     We need a virtual class to abstract the AOD and NanoAOD header classes
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include <TVector2.h>

#include "AliVHeader.h"
//#include "AliAODVertex.h"
#include <TString.h>
#include <TBits.h>
#include "AliCentrality.h"
#include "AliEventplane.h"

class TGeoHMatrix;
class TString;


class AliVAODHeader : public AliVHeader {

 public :
  AliVAODHeader() : AliVHeader() {};
 
  
  virtual ~AliVAODHeader() =0;

  virtual void     SetMagneticField(Double_t magFld)        = 0;
  virtual void     SetMuonMagFieldScale(Double_t magFldScl) = 0;
  virtual void     SetDiamond(Float_t xy[2],Float_t cov[3]) = 0; 
  virtual void     SetDiamondZ(Float_t z, Float_t sig2z)    = 0;
  virtual Int_t    GetRunNumber()                  const    = 0;
  virtual Double_t GetMagneticField()              const    = 0;
  virtual Double_t GetMuonMagFieldScale()          const    = 0;
  virtual Double_t GetDiamondX()                   const    = 0;
  virtual Double_t GetDiamondY()                   const    = 0;
  virtual Double_t GetDiamondZ()                   const    = 0;
  virtual void     GetDiamondCovXY(Float_t cov[3]) const    = 0;
  virtual Double_t GetSigma2DiamondX()             const    = 0;
  virtual Double_t GetSigma2DiamondY()             const    = 0;
  virtual Double_t GetSigma2DiamondZ()             const    = 0;

};

#endif
