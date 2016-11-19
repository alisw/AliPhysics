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
 
  
  virtual ~AliVAODHeader() {};

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

  virtual Bool_t   InitMagneticField()      const       = 0;
  virtual void     SetRunNumber(Int_t nRun)             = 0;
  virtual void     SetOrbitNumber(UInt_t nOr)           = 0;
  virtual void     SetPeriodNumber(UInt_t nPer)         = 0;
  virtual void     SetBunchCrossNumber(UShort_t nBx)    = 0;
  virtual void     SetTimeStamp(UInt_t t)               = 0;
  virtual void     SetEventType(UInt_t evttype)         = 0;
  virtual UInt_t   GetEventType()           const       = 0;
  virtual void     SetTriggerMask(ULong64_t trigMsk)    = 0;
  virtual void     SetTriggerMaskNext50(ULong64_t trigMsk) = 0;
  virtual ULong64_t GetTriggerMask()        const       = 0;
  virtual ULong64_t GetTriggerMaskNext50()  const       = 0;
  virtual void     SetTriggerCluster(UChar_t trigClus)  = 0;
  virtual void     SetFiredTriggerClasses(TString trig) = 0;
  virtual TString  GetFiredTriggerClasses() const       = 0;
  virtual Double_t GetZDCN1Energy()         const       = 0;
  virtual Double_t GetZDCP1Energy()         const       = 0;
  virtual Double_t GetZDCN2Energy()         const       = 0;
  virtual Double_t GetZDCP2Energy()         const       = 0;

  virtual Double_t GetZDCEMEnergy(Int_t /* i */)            const  = 0;
  virtual Int_t    GetNumberOfESDTracks()                   const  = 0;
  virtual UInt_t   GetNumberOfITSClusters(Int_t /* ilay */) const  = 0;
  virtual Float_t  GetT0spread(Int_t /* i */)               const  = 0;
  // FIXME: THIS IS UGLY!!!!                                       
  // FIXME: use dynamic cast in AliAODEVent?                       
  virtual AliCentrality* GetCentralityP()                   const  = 0;
  virtual AliEventplane* GetEventplaneP()                   const  = 0;
  virtual Double_t       GetCentrality () const = 0;
  virtual const Float_t* GetVZEROEqFactors()                const  = 0;
  virtual Float_t        GetVZEROEqFactors(Int_t /* i */)   const  = 0;
  virtual void           SetVZEROEqFactors(const Float_t* /*factors*/)  = 0;
  virtual UInt_t         GetOfflineTrigger()  = 0;
  virtual Int_t          GetRefMultiplicity()    const  =0;
  virtual Double_t       GetEventplane()         const =0;
  //
  virtual UInt_t  GetDAQAttributes()             const;
  virtual void    SetDAQAttributes(UInt_t);


};

#endif
