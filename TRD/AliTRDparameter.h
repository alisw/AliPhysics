#ifndef ALITRDPARAMETER_H
#define ALITRDPARAMETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD parameter class                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TObjArray;
class AliTRDgeometry;
class AliTRDpadPlane;

class AliTRDparameter : public TNamed {

 public:

  enum { kNplan = 6, kNcham = 5, kNsect = 18, kNdet = 540 };
 
  AliTRDparameter();
  AliTRDparameter(const Text_t* name, const Text_t* title);
  AliTRDparameter(const AliTRDparameter &p);   
  virtual ~AliTRDparameter();
  AliTRDparameter &operator=(const AliTRDparameter &p); 

  virtual void     Copy(TObject &p) const;
  virtual void     Init();
  virtual void     ReInit();
 
  virtual void     SetExpandTimeBin(Int_t nbefore, Int_t nafter)
                                                                  { fTimeBefore = nbefore;
                                                                    fTimeAfter  = nafter;       };



          Int_t    GetTimeMax()                             const { return fTimeMax; };
          Int_t    GetTimeBefore()                          const { return fTimeBefore;        }; 
          Int_t    GetTimeAfter()                           const { return fTimeAfter;         }; 

          Float_t  GetDriftVelocity()                       const { return fDriftVelocity;     };

  


          void     PrintDriftVelocity();


 protected:

  Int_t                fTimeMax;
  Int_t                fTimeBefore;                         //  Number of timebins before the drift region
  Int_t                fTimeAfter;                          //  Number of timebins after the drift region


  Float_t              fDriftVelocity;                      //  Drift velocity (cm / mus)



 private:


  ClassDef(AliTRDparameter,7)                               //  TRD parameter class

};

#endif
