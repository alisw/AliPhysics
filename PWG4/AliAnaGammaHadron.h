#ifndef ALIANAGAMMAHADRON_H
#define ALIANAGAMMAHADRON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/02/09 18:40:40  schutz
 * BNew version from Gustavo
 *
 * Revision 1.1  2007/01/23 17:17:29  schutz
 * New Gamma package
 *
 *
 */

//_________________________________________________________________________
//  Class for the analysis of gamma-jet correlations.     
//  Basically it seaches for a prompt photon in the Calorimeters acceptance, 
//  if so we construct a jet around the highest pt particle in the opposite 
//  side in azimuth. This jet has to fullfill several conditions to be 
//  accepted. Then the fragmentation function of this jet is constructed 
//  Class created from old AliPHOSGammaPion

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TROOT.h>
#include <TChain.h>
#include "TTask.h"
#include "TArrayD.h"
#include "TChain.h"
#include <TH2F.h>
#include <TTree.h> 
#include "AliAnaGammaDirect.h" 

class AliESD ; 
 
class AliAnaGammaHadron : public AliAnaGammaDirect {

public: 

  AliAnaGammaHadron(const char *name) ; // default ctor
  AliAnaGammaHadron(const AliAnaGammaHadron & gj) ; // cpy ctor
  AliAnaGammaHadron & operator = (const AliAnaGammaHadron & g) ;//cpy assignment
  virtual ~AliAnaGammaHadron() ; //virtual dtor

  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

  void InitParameters(); 
  Double_t GetAngleMaxParam(Int_t i) const {return fAngleMaxParam.At(i) ; }
  Double_t GetInvMassMaxCut() const {return fInvMassMaxCut ; }
  Double_t GetInvMassMinCut() const {return fInvMassMinCut ; }
  Double_t GetPhiMaxCut() const {return fPhiMaxCut ; }
  Double_t GetPhiMinCut() const {return fPhiMinCut ; }
  Float_t    GetMinPtPion() const {return fMinPtPion ; }

  void Print(const Option_t * opt)const;
  
  void SetAngleMaxParam(Int_t i, Double_t par)
  {fAngleMaxParam.AddAt(par,i) ; }
  void SetMinPtPion(Float_t pt){fMinPtPion = pt; };
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
    {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	
  void SetPhiCutRange(Double_t phimin, Double_t phimax)
  {fPhiMaxCut =phimax;  fPhiMinCut =phimin;}

  private:

  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e);
  void MakeGammaChargedCorrelation(TClonesArray * pl, TParticle *pGamma) const ;
  void MakeGammaNeutralCorrelation(TClonesArray * pl, TParticle *pGamma)  ;

 private:

  Double_t   fPhiMaxCut ;      // 
  Double_t   fPhiMinCut ;      // 
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
 
  Double_t   fMinPtPion;       // Minimum pt of pion
  TObjArray  *fOutputContainer ; //! output data container
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters

  //Histograms

  TH2F * fhPhiCharged  ; 
  TH2F * fhPhiNeutral   ; 
  TH2F * fhEtaCharged  ; 
  TH2F * fhEtaNeutral   ; 
  TH2F * fhDeltaPhiGammaCharged  ;  
  TH2F * fhDeltaPhiGammaNeutral   ; 
  TH2F * fhDeltaEtaGammaCharged  ; 
  TH2F * fhDeltaEtaGammaNeutral  ; 

  TH2F * fhCorrelationGammaNeutral  ; 
  TH2F * fhCorrelationGammaCharged  ; 

  TH2F * fhAnglePairAccepted  ; 
  TH2F * fhAnglePairNoCut  ; 
  TH2F * fhAnglePairAzimuthCut  ; 
  TH2F * fhAnglePairOpeningAngleCut   ; 
  TH2F * fhAnglePairAllCut   ; 
  TH2F * fhInvMassPairNoCut    ; 
  TH2F * fhInvMassPairAzimuthCut  ; 
  TH2F * fhInvMassPairOpeningAngleCut  ; 
  TH2F * fhInvMassPairAllCut   ; 
  
  ClassDef(AliAnaGammaHadron,0)
} ;
 

#endif //ALIANAGAMMAHADRON_H



