#ifndef ALIANAGAMMADIRECT_H
#define ALIANAGAMMADIRECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.5  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.4.4.3  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________

// Class for the analysis of prompt gamma, isolation cut. 
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "TObject.h" 
#include <TH2F.h>
#include <TNtuple.h>

class TList ;

class AliAnaGammaDirect : public TObject {

public: 

  AliAnaGammaDirect() ; // default ctor
  AliAnaGammaDirect(const AliAnaGammaDirect & g) ; // cpy ctor
  AliAnaGammaDirect & operator = (const AliAnaGammaDirect & g) ;//cpy assignment
  virtual ~AliAnaGammaDirect() ; //virtual dtor
  
  enum anatype_t {kNoIC, kPtIC, kSumPtIC, kSeveralIC};
  
  Double_t  GetMinGammaPt()    const {return fMinGammaPt ; }
  Float_t     GetConeSize()          const {return fConeSize ; }
  Float_t     GetPtThreshold()      const {return fPtThreshold ; }
  Float_t     GetPtSumThres()     const {return fPtSumThreshold ; }
  Int_t        GetICMethod()          const {return fICMethod ; }
  Bool_t     IsMC() const {return fAnaMC ; };

  TList *  GetCreateOutputObjects();
  void GetPromptGamma(TClonesArray * plNe, TClonesArray * plCTS, TClonesArray * plPrimNe,  TParticle * pGamma, Bool_t &Is)  const;
  
  void MakeSeveralICAnalysis(TClonesArray * plCalo, TClonesArray * plCTS); 
  void MakeIsolationCut(TClonesArray * plCTS, TClonesArray * plNe, 
			TParticle *pCandidate, Int_t index, Int_t &n,
			Bool_t &imcpt, Bool_t &icms, Float_t &ptsum) const ;  
  
  void Print(const Option_t * opt)const;
  
  void SetMinGammaPt(Double_t ptcut){fMinGammaPt =ptcut;}
  void SetConeSize(Float_t r)              {fConeSize = r ; }
  void SetPtThreshold(Float_t pt)        {fPtThreshold = pt; };
  void SetPtSumThreshold(Float_t pt) {fPtSumThreshold = pt; };
  void SetICMethod(Int_t i )          {fICMethod = i ; }
  void SetMC()    {fAnaMC = kTRUE ; }

  Int_t    GetNCones()                  const {return fNCones ; }
  Int_t    GetNPtThresholds()                const {return fNPtThres ; }
  Float_t GetConeSizes(Int_t i)      const {return fConeSizes[i] ; }
  Float_t GetPtThresholds(Int_t i)  const {return fPtThresholds[i] ; }
  
  void InitParameters();
 
  void SetNCones(Int_t ncs)              {fNCones = ncs ; }
  void SetNPtThresholds(Int_t npt)        {fNPtThres = npt; }
  void SetConeSizes(Int_t i, Float_t r)         {fConeSizes[i] = r ; }
  void SetPtThresholds(Int_t i, Float_t pt)   {fPtThresholds[i] = pt; }

  
  private:
     
  Double_t    fMinGammaPt ;  // Min pt in Calorimeter
  Float_t      fConeSize ; //Size of the isolation cone 
  Float_t      fPtThreshold ; //Mimium pt of the particles in the cone to set isolation
  Float_t      fPtSumThreshold ; //Mimium pt sum of the particles in the cone to set isolation  
  Int_t        fICMethod ; //Isolation cut method to be used
                                           // kNoIC: No isolation
                                           // kPtIC: Pt threshold method
                                           // kSumPtIC: Cone pt sum method
                                           // kSeveralIC: Analysis for several cuts
  Bool_t fAnaMC ; //Set in case of using MCData reader 

  //Histograms  
  TH1F * fhNGamma    ;  //Number of (isolated) gamma identified
  TH2F * fhPhiGamma  ; // Phi of identified gamma
  TH2F * fhEtaGamma  ; // eta of identified gamma
  TH2F * fhConeSumPt ; // Sum Pt in the cone

  TNtuple *    fntuplePrompt ; //List of found prompt photons, pt, eta and phi. Also primary information.

  //Prompt photon analysis data members for multiple cones and pt thresholds kIsolationCut
  Int_t         fNCones   ; //Number of cone sizes to test
  Int_t         fNPtThres ; //Number of ptThres to test
  Float_t     fConeSizes[10] ; // Array with cones to test
  Float_t     fPtThresholds[10] ; // Array with pt thresholds to test
  
  TH1F* fhPtThresIsolated[20][20]; // Isolated gamma with pt threshold 
  TH2F* fhPtSumIsolated[20] ;  //  Isolated gamma with threshold on cone pt sume
  TNtuple *    fntSeveralIC[20] ; //ntuple 

  ClassDef(AliAnaGammaDirect,1)
} ;
 

#endif //ALIANAGAMMADIRECT_H



