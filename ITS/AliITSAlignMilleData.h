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

  Int_t    *GetIdxlocX()  {return fIdxlocX;}
  Int_t    *GetIdxgloX()  {return fIdxgloX;}
  Double_t *GetDerlocX()  {return fDerlocX;}
  Double_t *GetDergloX()  {return fDergloX;}    

  Int_t    *GetIdxlocZ()  {return fIdxlocZ;}
  Int_t    *GetIdxgloZ()  {return fIdxgloZ;}
  Double_t *GetDerlocZ()  {return fDerlocZ;}
  Double_t *GetDergloZ()  {return fDergloZ;}    

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
