#ifndef ALI_TPC_INVERSE_CORRECTION_H
#define ALI_TPC_INVERSE_CORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliTPCInverseCorrection class                                              //
//                                                                            //
// This is a wrapper that inverts an AliTPCCorrection. This is done by        //
// swapping the CalculateCorrection and CalculateInverseCorrection functions. //
// The wrapped correction is supplied as a pointer and the class relies       //
// on the fact, that this pointer keeps pointing to the right object.         //
// However, the ownership is not changed, i.e. the wrapped correction         //
// will not be deleted when this correction is destructed.                    //
//                                                                            //
// date: 27/04/2010                                                           //
// Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                       //
////////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"


class AliTPCInverseCorrection : public AliTPCCorrection {
public:
  AliTPCInverseCorrection();
  AliTPCInverseCorrection(AliTPCCorrection *correction);
  virtual ~AliTPCInverseCorrection();

  void SetCorrection(const AliTPCCorrection *correction) {fCorrection=(AliTPCCorrection*) correction;}
  const AliTPCCorrection* GetCorrection() const {return fCorrection;}
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
  virtual void GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]);

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2);
 
  // convenience functions
  virtual void Print(Option_t* option="") const;
 
private:
  AliTPCCorrection *fCorrection; // The correction to be inverted.

  AliTPCInverseCorrection & operator = (const AliTPCInverseCorrection);
  AliTPCInverseCorrection(const AliTPCInverseCorrection&); //dummy copy contructor

  ClassDef(AliTPCInverseCorrection,1);
};

#endif
