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
  
  virtual UShort_t  GetBunchCrossNumber()   const { NotImplemented();return 0;}
  virtual UInt_t    GetOrbitNumber()        const { NotImplemented();return 0;}
  virtual UInt_t    GetPeriodNumber()       const { NotImplemented();return 0;}
  virtual UInt_t    GetTimeStamp()          const { NotImplemented();return 0;}
  virtual ULong64_t GetTriggerMask()        const { NotImplemented();return 0;}
  virtual ULong64_t GetTriggerMaskNext50()  const { NotImplemented();return 0;}
  virtual UChar_t   GetTriggerCluster()     const { NotImplemented();return 0;}
  virtual UInt_t    GetEventType()          const { NotImplemented();return 0;}

  virtual Bool_t   InitMagneticField()             const    {NotImplemented(); return 0;};
  virtual void     SetRunNumber(Int_t /*n*/)                    {NotImplemented();};
  virtual void     SetMagneticField(Double_t /*magFld*/)        {NotImplemented();};
  virtual void     SetMuonMagFieldScale(Double_t /*magFldScl*/) {NotImplemented();};
  virtual void     SetDiamond(Float_t */*xy[2]*/,Float_t */*cov[3]*/) {NotImplemented();}; 
  virtual void     SetDiamondZ(Float_t /*z*/, Float_t /*sig2z*/)    {NotImplemented();};
  virtual Int_t    GetRunNumber()                  const    {NotImplemented(); return 0;};
  virtual Double_t GetMuonMagFieldScale()          const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondX()                   const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondY()                   const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondZ()                   const    {NotImplemented(); return 0;};
  virtual void     GetDiamondCovXY(Float_t */*cov[3]*/) const    {NotImplemented();};
  virtual Double_t GetSigma2DiamondX()             const    {NotImplemented(); return 0;};
  virtual Double_t GetSigma2DiamondY()             const    {NotImplemented(); return 0;};
  virtual Double_t GetSigma2DiamondZ()             const    {NotImplemented(); return 0;};

  virtual void     SetOrbitNumber(UInt_t /* nOr */)           {NotImplemented(); };
  virtual void     SetPeriodNumber(UInt_t /* nPer */)         {NotImplemented(); };
  virtual void     SetBunchCrossNumber(UShort_t /* nBx */)    {NotImplemented(); };
  virtual void     SetTimeStamp(UInt_t /* t */)               {NotImplemented(); };
  virtual void     SetEventType(UInt_t /* evttype */)         {NotImplemented(); };
  virtual void     SetTriggerMask(ULong64_t /* trigMsk */)    {NotImplemented(); };
  virtual void     SetTriggerMaskNext50(ULong64_t /* trigMsk */) {NotImplemented(); };
  virtual void     SetTriggerCluster(UChar_t /* trigClus */)  {NotImplemented(); };
  virtual void     SetFiredTriggerClasses(TString /* trig */) {NotImplemented(); };
  virtual TString  GetFiredTriggerClasses() const             {NotImplemented(); return "";};
  virtual Double_t GetZDCN1Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCP1Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCN2Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCP2Energy()         const             {NotImplemented(); return 0;};

  virtual Double_t GetZDCEMEnergy(Int_t /* i */) const            {NotImplemented(); return 0;};
  virtual Int_t    GetNumberOfESDTracks()  const            {NotImplemented(); return 0;};
  virtual UInt_t   GetNumberOfITSClusters(Int_t /* ilay */) const {NotImplemented(); return 0;};
  virtual Float_t  GetT0spread(Int_t /* i */)               const {NotImplemented(); return 0;};
  // FIXME: THIS IS UGLY!!!!
  // FIXME: use dynamic cast in AliAODEVent?
  virtual AliCentrality* GetCentralityP()  const {NotImplemented(); return 0;};
  virtual AliEventplane* GetEventplaneP()  const {NotImplemented(); return 0;};
  virtual Double_t       GetEventplane()     const {NotImplemented(); return 0;};
  virtual const Float_t* GetVZEROEqFactors() const {NotImplemented(); return 0;};
  virtual Float_t        GetVZEROEqFactors(Int_t /* i */) const {NotImplemented(); return 0;};
  virtual void           SetVZEROEqFactors(const Float_t* /*factors*/) {NotImplemented(); } 

  virtual UInt_t GetOfflineTrigger()  {NotImplemented(); return 0;};


  virtual void Print(Option_t* /*option = ""*/) const  {Printf("I'm a special header!");}
 
  virtual void Clear(Option_t * opt) ;

  virtual Int_t  GetIRInt2ClosestInteractionMap()                  const {NotImplemented(); return 0;};
  virtual Int_t  GetIRInt1ClosestInteractionMap(Int_t /*gap = 3*/) const {NotImplemented(); return 0;};


  virtual Int_t     GetRefMultiplicity()    const { NotImplemented(); return 0; }

  Double_t  GetMagneticField()      const { return GetVar(1); }
  Double_t  GetCentrality (/*estimator = "V0M"*/) const { return GetVar(0);}
  
  ClassDef(AliNanoAODHeader, 1)
private:
  void NotImplemented() const;
};

#endif /* _ALINANOAODHEADER_H_ */
