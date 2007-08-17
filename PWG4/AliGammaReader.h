#ifndef ALIGAMMAREADER_H
#define ALIGAMMAREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for reading data in order to do prompt gamma correlations 
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "AliStack.h"
#include "TObject.h" 
 
class AliESD ; 

class TH2F ; 

class AliGammaReader : public TObject {

public: 

  AliGammaReader() ; // ctor
  AliGammaReader(const AliGammaReader & g) ; // cpy ctor
  AliGammaReader & operator = (const AliGammaReader & g) ;//cpy assignment
  virtual ~AliGammaReader() {;} //virtual dtor

  enum datatype_t {kData, kMC, kMCData};
  
  void InitParameters();

  Int_t GetDataType(){ return fDataType ; }
  void SetDataType(Int_t data ){fDataType = data ; }

  virtual Float_t   GetCTSEtaCut() const {return fCTSEtaCut ; }
  virtual Float_t   GetEMCALEtaCut() const {return fEMCALEtaCut ; }
  virtual Float_t   GetPHOSEtaCut() const {return fPHOSEtaCut ; }
  virtual Float_t  GetPhiEMCALCut(Int_t i) { return  fPhiEMCALCut[i]  ; }
  virtual Float_t  GetPhiPHOSCut(Int_t i) { return  fPhiPHOSCut[i]  ; }
  virtual Float_t  GetNeutralPtCut()    {  return fNeutralPtCut  ; }
  virtual Float_t  GetChargedPtCut()  {  return fChargedPtCut  ; }

  virtual void Print(const Option_t * opt)const;
  
  virtual void SetCTSEtaCut(Float_t eta){ fCTSEtaCut= eta ; }
  virtual void SetEMCALEtaCut(Float_t eta){ fEMCALEtaCut= eta ; }
  virtual void SetPHOSEtaCut(Float_t eta){ fPHOSEtaCut= eta ; }
  virtual void SetPhiEMCALCut(Float_t  phi0, Float_t  phi1)
  { fPhiEMCALCut[0]= phi0 ; fPhiEMCALCut[1]= phi1 ;}
  virtual void SetPhiPHOSCut(Float_t  phi0, Float_t  phi1)
  { fPhiPHOSCut[0]= phi0 ; fPhiPHOSCut[1]= phi1 ;}
  virtual void SetNeutralPtCut(Float_t  pt){  fNeutralPtCut = pt ; }
  virtual void SetChargedPtCut(Float_t  pt){  fChargedPtCut = pt ; }

  virtual void CreateParticleList(TObject* data, TObject * data2, 
				  TClonesArray * plCh, TClonesArray * plEMCAL, 
				  TClonesArray * plPHOS, TClonesArray * parton) {;}
 protected:
  Int_t        fDataType ;
  Float_t      fCTSEtaCut ;//CTS  pseudorapidity acceptance
  Float_t      fEMCALEtaCut ;//EMCAL pseudorapidity acceptance
  Float_t      fPHOSEtaCut ;//PHOS pseudorapidity acceptance
  Float_t      fPhiEMCALCut[2]; //EMCAL phi acceptance 
  Float_t      fPhiPHOSCut[2];  //PHOS phi acceptance
  Float_t      fNeutralPtCut; //
  Float_t      fChargedPtCut;  // 

   ClassDef(AliGammaReader,0)
} ;
 

#endif //ALIGAMMAREADER_H



