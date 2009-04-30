#ifndef ALIISOLATIONCUT_H
#define ALIISOLATIONCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  */

//_________________________________________________________________________

// Class containing methods for the isolation cut. 
//
//
// 

// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system --- 
#include <TObject.h>
class TRefArray ;

// --- ANALYSIS system ---
class AliAODPWG4ParticleCorrelation ;

class AliIsolationCut : public TObject {
  
 public: 
  
  AliIsolationCut() ; // default ctor
  AliIsolationCut(const AliIsolationCut & g) ; // cpy ctor
  AliIsolationCut & operator = (const AliIsolationCut & g) ;//cpy assignment
  virtual ~AliIsolationCut() {;} //virtual dtalr
  
  enum type {kPtThresIC, kSumPtIC, kPtFracIC, kSumPtFracIC};
  
  Float_t    GetConeSize()        const {return fConeSize ; }
  Float_t    GetPtThreshold()     const {return fPtThreshold ; }
  Float_t    GetPtFraction()      const {return fPtFraction ; }
  Int_t      GetICMethod()        const {return fICMethod ; }
  TString    GetICParametersList() ; 
  
  void MakeIsolationCut(TRefArray * plCTS, TRefArray * plNe, Double_t * vertex, 
			const Bool_t fillAOD, AliAODPWG4ParticleCorrelation  * pCandidate, const TString aodRefArrayName,
			Int_t &n, Int_t & nfrac, Float_t &ptsum, Bool_t & isolated) const ;  
  
  void Print(const Option_t * opt)const;
  
  void SetConeSize(Float_t r)       {fConeSize = r ; }
  void SetPtThreshold(Float_t pt)   {fPtThreshold = pt; }
  void SetPtFraction(Float_t pt)    {fPtFraction = pt; }
  void SetICMethod(Int_t i )        {fICMethod = i ; }
  
  void InitParameters();
  
  
 private:
  
  Float_t      fConeSize ;    //Size of the isolation cone 
  Float_t      fPtThreshold ; //Mimium pt of the particles in the cone or sum in cone 
  Float_t      fPtFraction ;  //Fraction of the momentum of particles in cone or sum in cone
  Int_t        fICMethod ;    //Isolation cut method to be used
                              // kPtIC: Pt threshold method
                              // kSumPtIC: Cone pt sum method
                              // kPtFracIC:   Pt threshold, fraction of candidate pt, method
                              // kSumPtFracIC:   Cone pt sum , fraction of cone sum, method
  
  ClassDef(AliIsolationCut,1)
} ;


#endif //ALIISOLATIONCUT_H



