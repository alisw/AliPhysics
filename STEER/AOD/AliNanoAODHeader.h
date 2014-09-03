#ifndef _ALINANOAODHEADER_H_
#define _ALINANOAODHEADER_H_

#include "AliVAODHeader.h"
#include "AliNanoAODStorage.h"



class AliNanoAODHeader : public AliVAODHeader, public AliNanoAODStorage
{
public:
  using AliVHeader::ClassName;
  AliNanoAODHeader()  {;}
  AliNanoAODHeader(Int_t size){ AllocateInternalStorage(size);}
  virtual ~AliNanoAODHeader(){;}


  // Interface methods
  // AliNanoAODHeader(const AliVHeader& evt); 
  AliNanoAODHeader& operator=(const AliNanoAODHeader& evt);
  
  virtual UShort_t  GetBunchCrossNumber()   const { AliError("Not implemented");return 0;}
  virtual UInt_t    GetOrbitNumber()        const { AliError("Not implemented");return 0;}
  virtual UInt_t    GetPeriodNumber()       const { AliError("Not implemented");return 0;}
  virtual ULong64_t GetTriggerMask()        const { AliError("Not implemented");return 0;}
  virtual UChar_t   GetTriggerCluster()     const { AliError("Not implemented");return 0;}
  virtual UInt_t    GetEventType()          const { AliError("Not implemented");return 0;}

  virtual Bool_t   InitMagneticField()             const    {AliError("Not Implemented"); return 0;};
  virtual void     SetRunNumber(Int_t /*n*/)                    {AliError("Not Implemented");};
  virtual void     SetMagneticField(Double_t /*magFld*/)        {AliError("Not Implemented");};
  virtual void     SetMuonMagFieldScale(Double_t /*magFldScl*/) {AliError("Not Implemented");};
  virtual void     SetDiamond(Float_t */*xy[2]*/,Float_t */*cov[3]*/) {AliError("Not Implemented");}; 
  virtual void     SetDiamondZ(Float_t /*z*/, Float_t /*sig2z*/)    {AliError("Not Implemented");};
  virtual Int_t    GetRunNumber()                  const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetMuonMagFieldScale()          const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetDiamondX()                   const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetDiamondY()                   const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetDiamondZ()                   const    {AliError("Not Implemented"); return 0;};
  virtual void     GetDiamondCovXY(Float_t */*cov[3]*/) const    {AliError("Not Implemented");};
  virtual Double_t GetSigma2DiamondX()             const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetSigma2DiamondY()             const    {AliError("Not Implemented"); return 0;};
  virtual Double_t GetSigma2DiamondZ()             const    {AliError("Not Implemented"); return 0;};

  virtual void     SetOrbitNumber(UInt_t /* nOr */)           {AliError("Not Implemented"); };
  virtual void     SetPeriodNumber(UInt_t /* nPer */)         {AliError("Not Implemented"); };
  virtual void     SetBunchCrossNumber(UShort_t /* nBx */)    {AliError("Not Implemented"); };
  virtual void     SetEventType(UInt_t /* evttype */)         {AliError("Not Implemented"); };
  virtual void     SetTriggerMask(ULong64_t /* trigMsk */)    {AliError("Not Implemented"); };
  virtual void     SetTriggerCluster(UChar_t /* trigClus */)  {AliError("Not Implemented"); };
  virtual void     SetFiredTriggerClasses(TString /* trig */) {AliError("Not Implemented"); };
  virtual TString  GetFiredTriggerClasses() const             {AliError("Not Implemented"); return "";};
  virtual Double_t GetZDCN1Energy()         const             {AliError("Not Implemented"); return 0;};
  virtual Double_t GetZDCP1Energy()         const             {AliError("Not Implemented"); return 0;};
  virtual Double_t GetZDCN2Energy()         const             {AliError("Not Implemented"); return 0;};
  virtual Double_t GetZDCP2Energy()         const             {AliError("Not Implemented"); return 0;};

  virtual Double_t GetZDCEMEnergy(Int_t /* i */) const            {AliError("Not Implemented"); return 0;};
  virtual Int_t    GetNumberOfESDTracks()  const            {AliError("Not Implemented"); return 0;};
  virtual UInt_t   GetNumberOfITSClusters(Int_t /* ilay */) const {AliError("Not Implemented"); return 0;};
  virtual Float_t  GetT0spread(Int_t /* i */)               const {AliError("Not Implemented"); return 0;};
  // FIXME: THIS IS UGLY!!!!
  // FIXME: use dynamic cast in AliAODEVent?
  virtual AliCentrality* GetCentralityP()  const {AliError("Not Implemented"); return 0;};
  virtual AliEventplane* GetEventplaneP()  const {AliError("Not Implemented"); return 0;};
  virtual Double_t       GetEventplane()     const {AliError("Not Implemented"); return 0;};
  virtual const Float_t* GetVZEROEqFactors() const {AliError("Not Implemented"); return 0;};
  virtual Float_t        GetVZEROEqFactors(Int_t /* i */) const {AliError("Not Implemented"); return 0;};
  virtual void           SetVZEROEqFactors(const Float_t* /*factors*/) {AliError("Not Implemented"); } 

  virtual UInt_t GetOfflineTrigger()  {AliError("Not Implemented"); return 0;};


  virtual void Print(Option_t* /*option = ""*/) const  {Printf("I'm a special header!");}
 
  virtual void Clear(Option_t * opt) ;

  virtual Int_t  GetIRInt2ClosestInteractionMap()                  const {AliError("Not Implemented"); return 0;};
  virtual Int_t  GetIRInt1ClosestInteractionMap(Int_t /*gap = 3*/) const {AliError("Not Implemented"); return 0;};


  virtual Int_t     GetRefMultiplicity()    const { AliError("Not Impletented"); return 0; }

  Double_t  GetMagneticField()      const { return GetVar(1); }
  Double_t  GetCentrality (/*estimator = "V0M"*/) const { return GetVar(0);}
  
  ClassDef(AliNanoAODHeader, 1)

};

#endif /* _ALINANOAODHEADER_H_ */
