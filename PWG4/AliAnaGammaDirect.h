#ifndef ALIANAGAMMADIRECT_H
#define ALIANAGAMMADIRECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.3  2007/03/08 10:24:32  schutz
 * Coding convention
 *
 * Revision 1.2  2007/02/09 18:40:40  schutz
 * BNew version from Gustavo
 *
 * Revision 1.1  2007/01/23 17:17:29  schutz
 * New Gamma package
 *
 *
 */

//_________________________________________________________________________

// Class for the analysis of gamma 
// This class only contains 3 methods: one to fill lists of particles (ESDs) comming 
//  from the CTS (ITS+TPC) and the calorimeters; the other search in the 
//  corresponing calorimeter for the highest energy cluster, identify it as 
//  prompt photon(Shower Shape and Isolation Cut), the last method does the 
//  isolation selection
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include <TTree.h> 
#include "AliAnalysisTask.h" 
#include <TH2F.h>

class AliESD ; 
 
class AliAnaGammaDirect : public AliAnalysisTask {

public: 

  AliAnaGammaDirect(const char *name) ; // default ctor
  AliAnaGammaDirect(const AliAnaGammaDirect & g) ; // cpy ctor
  AliAnaGammaDirect & operator = (const AliAnaGammaDirect & g) ;//cpy assignment
  virtual ~AliAnaGammaDirect() ; //virtual dtor

  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");
  
  void InitParameters();
  TTree *     GetChain()                const {return fChain ; }
  AliESD *    GetESD()                  const {return fESD ; }
  TObjArray * GetOutputContainer()      const {return fOutputContainer ; }
  Double_t  GetMinGammaPt()    const {return fMinGammaPt ; }
  TString    GetCalorimeter()       const {return fCalorimeter ; }
  Bool_t      GetPrintInfo()           const {return fPrintInfo ; }
  Float_t     GetConeSize()          const {return fConeSize ; }
  Float_t     GetPtThreshold()      const {return fPtThreshold ; }
  Float_t     GetPtSumThres()     const {return fPtSumThreshold ; }
  Int_t        GetICMethod()          const {return fMakeICMethod ; }

  Bool_t   IsEMCALPIDOn() const {return fEMCALPID ; }
  Bool_t   IsPHOSPIDOn() const {return fPHOSPID ; }
  Float_t  GetEMCALPhotonWeight() { return  fEMCALPhotonWeight  ; }
  Float_t  GetEMCALPi0Weight()    {  return fEMCALPi0Weight  ; }
  Float_t  GetPHOSPhotonWeight()  {  return fPHOSPhotonWeight  ; }

  void Print(const Option_t * opt)const;

  void SetMinGammaPt(Double_t ptcut){fMinGammaPt =ptcut;}
  void SetCalorimeter(TString calo){ fCalorimeter= calo ; }
  void SetPrintInfo(Bool_t print){ fPrintInfo = print ; }
  void SetConeSize(Float_t r)              {fConeSize = r ; }
  void SetPtThreshold(Float_t pt)        {fPtThreshold = pt; };
  void SetPtSumThreshold(Float_t pt) {fPtSumThreshold = pt; };
  void SetICMethod(Int_t i )          {fMakeICMethod = i ; }
  
  void SetEMCALPIDOn(Bool_t pid){ fEMCALPID= pid ; }
  void SetPHOSPIDOn(Bool_t pid){ fPHOSPID= pid ; }
  void SetEMCALPhotonWeight(Float_t  w){  fEMCALPhotonWeight = w ; }
  void SetEMCALPi0Weight(Float_t  w){  fEMCALPi0Weight = w ; }
  void SetPHOSPhotonWeight(Float_t  w){  fPHOSPhotonWeight = w ; }

  void CreateParticleList(TClonesArray * particleList, 
			  TClonesArray * plCh, TClonesArray * plNe, 
			  TClonesArray * plNePHOS);
  
  
  void GetPromptGamma(TClonesArray * plNe, TClonesArray * plCTS, TParticle * pGamma, Bool_t &Is)  const;
  
  void MakeIsolationCut(TClonesArray * plCTS, TClonesArray * plNe, 
			TParticle *pCandidate, Int_t index, 
			Bool_t &imcpt, Bool_t &icms, Float_t &ptsum) const ;  
  
 private:

  TTree       *fChain ;   //!pointer to the analyzed TTree or TChain
  AliESD       *fESD ;     //! Declaration of leave types
  TObjArray  *fOutputContainer ; //! output data container
  Bool_t        fPrintInfo ;      //Print most interesting information on screen
  Double_t    fMinGammaPt ;  // Min pt in Calorimeter
  TString      fCalorimeter ; //PHOS or EMCAL detects Gamma
  Bool_t       fEMCALPID ;//Fill EMCAL particle lists with particles with corresponding pid
  Bool_t       fPHOSPID;  //Fill PHOS particle lists with particles with corresponding pid
  Float_t      fEMCALPhotonWeight; //Bayesian PID weight for photons in EMCAL 
  Float_t      fEMCALPi0Weight;  //Bayesian PID weight for pi0 in EMCAL 
  Float_t      fPHOSPhotonWeight; //Bayesian PID weight for photons in PHOS 
  Float_t      fConeSize ; //Size of the isolation cone 
  Float_t      fPtThreshold ; //Mimium pt of the particles in the cone to set isolation
  Float_t      fPtSumThreshold ; //Mimium pt sum of the particles in the cone to set isolation  
  Int_t        fMakeICMethod ; //Isolation cut method to be used
                                           // 0: No isolation
                                           // 1: Pt threshold method
                                           // 2: Cone pt sum method
  //Histograms  
  TH1F * fhNGamma    ; 
  TH2F * fhPhiGamma    ; 
  TH2F * fhEtaGamma    ; 

  ClassDef(AliAnaGammaDirect,0)
} ;
 

#endif //ALIANAGAMMADIRECT_H



