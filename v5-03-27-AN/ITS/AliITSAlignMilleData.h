#ifndef ALIITSALIGNMILLEDATA_H
#define ALIITSALIGNMILLEDATA_H
/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


/// \ingroup rec
/// \class AliITSAlignMilleData
/// \brief Class for alignment of ITS
//
// Authors: Marcello Lunardon
#include <TObject.h> 

#define ITSMILLENPARCH         6
#define ITSMILLENLOCAL         5

class AliITSAlignMilleData : public TObject
{
 public:
  AliITSAlignMilleData();
  virtual ~AliITSAlignMilleData();
  Double_t  GetMeasX() const {return fMeasX;}
  Double_t  GetSigmaX() const {return fSigmaX;}
  void      SetMeasX(Double_t meas) {fMeasX=meas;}
  void      SetSigmaX(Double_t meas) {fSigmaX=meas;}

  Double_t  GetMeasZ() const {return fMeasZ;}
  Double_t  GetSigmaZ() const {return fSigmaZ;}
  void      SetMeasZ(Double_t meas) {fMeasZ=meas;}
  void      SetSigmaZ(Double_t meas) {fSigmaZ=meas;}

  Int_t    *GetIdxlocX() const {return (Int_t*)fIdxlocX;}
  Int_t    *GetIdxgloX() const {return (Int_t*)fIdxgloX;}
  Double_t *GetDerlocX() const {return (Double_t*)fDerlocX;}
  Double_t *GetDergloX() const {return (Double_t*)fDergloX;}    

  Int_t    *GetIdxlocZ() const {return (Int_t*)fIdxlocZ;}
  Int_t    *GetIdxgloZ() const {return (Int_t*)fIdxgloZ;}
  Double_t *GetDerlocZ() const {return (Double_t*)fDerlocZ;}
  Double_t *GetDergloZ() const {return (Double_t*)fDergloZ;}    

 private:  
  /// structure to store data for 2 LocalEquations (X and Z)
  Double_t fMeasX;  ///
  Double_t fSigmaX; ///
  Int_t    fIdxlocX[ITSMILLENLOCAL]; ///
  Double_t fDerlocX[ITSMILLENLOCAL]; ///
  Int_t    fIdxgloX[ITSMILLENPARCH]; ///
  Double_t fDergloX[ITSMILLENPARCH]; ///

  Double_t fMeasZ;  ///
  Double_t fSigmaZ; ///
  Int_t    fIdxlocZ[ITSMILLENLOCAL]; ///
  Double_t fDerlocZ[ITSMILLENLOCAL]; ///
  Int_t    fIdxgloZ[ITSMILLENPARCH]; ///
  Double_t fDergloZ[ITSMILLENPARCH]; ///

  //AliITSAlignMilleData(const AliITSAlignMilleData& rhs);
  //AliITSAlignMilleData& operator=(const AliITSAlignMilleData& rhs);

  ClassDef(AliITSAlignMilleData, 0)

};

#endif
