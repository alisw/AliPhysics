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



  virtual void Print(Option_t* /*option = ""*/) const  {Printf("I'm a special header!");}
 
  virtual void Clear(Option_t * opt) ;


  Double_t  GetMagneticField()      const { return GetVar(1); }
  Double_t  GetCentrality (const char *estimator = "V0M") { return GetVar(0);}
  
  ClassDef(AliNanoAODHeader, 1)

};

#endif /* _ALINANOAODHEADER_H_ */
