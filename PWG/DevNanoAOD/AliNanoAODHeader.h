#ifndef _ALINANOAODHEADER_H_
#define _ALINANOAODHEADER_H_

#include "AliVHeader.h"
#include "AliNanoAODStorage.h"



class AliNanoAODHeader : public AliVHeader, public AliNanoAODStorage
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
  virtual void Print(Option_t* /*option = ""*/) const  {Printf("I'm a special header!");}
 
  virtual void Clear(Option_t * opt) ;


  Double_t  GetMagneticField()      const { return GetVar(1); }
  Double_t  GetCentrality (const char *estimator = "V0M") { return GetVar(0);}
  
  ClassDef(AliNanoAODHeader, 1)

};

#endif /* _ALINANOAODHEADER_H_ */
