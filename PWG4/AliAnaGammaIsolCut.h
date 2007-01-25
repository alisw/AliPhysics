#ifndef ALIANAGAMMAISOLCUT_H
#define ALIANAGAMMAISOLCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1  2007/01/23 17:17:29  schutz
 * New Gamma package
 *
 *
 */

//_________________________________________________________________________

// Class for the analysis of gamma  (gamma-jet, 
// gamma-hadron(Arleo, TODO))
// This class make the Isolation Cut analysis for 2 methods, pt cone sum 
// and pt particle threshold and with different  cone sizes and pt thresholds.

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TROOT.h>
#include <TChain.h>
#include "TTask.h"
#include "TArrayD.h"
#include "TChain.h"
#include <TH2F.h>
#include <TTree.h> 
#include <TParticle.h> 
#include "AliAnaGammaDirect.h" 

class AliESD ; 
 
class AliAnaGammaIsolCut : public AliAnaGammaDirect {

public: 

  AliAnaGammaIsolCut(const char *name) ; // default ctor
  AliAnaGammaIsolCut(const AliAnaGammaIsolCut & g) ; // cpy ctor
  virtual ~AliAnaGammaIsolCut() ; //virtual dtor
  virtual void Exec(Option_t * opt = "") ;
  virtual void Init(Option_t * opt = "");
  virtual void Terminate(Option_t * opt = "");

  void Print(const Option_t * opt)const;

  Int_t    GetNCones()                  const {return fNCones ; }
  Int_t    GetNPtThres()                const {return fNPtThres ; }
  Float_t GetConeSizes(Int_t i)      const {return fConeSizes[i] ; }
  Float_t GetPtThresholds(Int_t i)  const {return fPtThresholds[i] ; }
  
  void SetNCone(Int_t ncs)              {fNCones = ncs ; }
  void SetNPtThres(Int_t npt)        {fNPtThres = npt; }
  void SetConeSizes(Int_t i, Float_t r)         {fConeSizes[i] = r ; }
  void SetPtThresholds(Int_t i, Float_t pt)   {fPtThresholds[i] = pt; }
  
  void MakeHistos() ;
  
  
 private:

  TObjArray  *fOutputContainer ; //! output data container

  Int_t         fNCones   ; //Number of cone sizes to test
  Int_t         fNPtThres ; //Number of ptThres to test
  Float_t     fConeSizes[10] ; // Arrat with cones to test
  Float_t     fPtThresholds[10] ; // Array with pt thresholds to test

  //Histograms  
 
  TH1F * fhPtCandidate ;
  TH1F* fhPtThresIsolated[10][10] ;
  TH2F* fhPtSumIsolated[10] ;

  ClassDef(AliAnaGammaIsolCut,0)
} ;
 

#endif //ALIANAGAMMAISOLCUT_H



