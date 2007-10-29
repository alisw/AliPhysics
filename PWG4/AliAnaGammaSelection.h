#ifndef ALIANAGAMMASELECTION_H
#define ALIANAGAMMASELECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 *
 *
 */

//_________________________________________________________________________
// Class for plotting particle/cluster/track distributions without cuts 
// and select clusters/tracks/particles needed in the analysis
// depending on PID criteria or other. 
//
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "TObject.h" 
#include <TH2F.h>
#include <TNtuple.h>
#include <TString.h>

class TList ;

class AliAnaGammaSelection : public TObject {

public: 

  AliAnaGammaSelection() ; // default ctor
  AliAnaGammaSelection(const AliAnaGammaSelection & g) ; // cpy ctor
  AliAnaGammaSelection & operator = (const AliAnaGammaSelection & g) ;//cpy assignment
  virtual ~AliAnaGammaSelection() ; //virtual dtor
  
  enum Det {kEMCAL, kPHOS, kCTS};

  void InitParameters();
  
  Bool_t     IsMC() const {return fAnaMC ; };
  Bool_t      FillCTS() const {return fFillCTS ;};

  TList *  GetCreateOutputObjects();
  void Selection(TString det, TClonesArray * pl, TClonesArray * plPrim)  const;
  
  void Print(const Option_t * opt)const;

  void SetFillCTS(Bool_t f) {fFillCTS = f;}
  void SetMC()    {fAnaMC = kTRUE ; }
  
  
  private:
  Bool_t fAnaMC ; //Set in case of using MCData reader 
  Bool_t fFillCTS ; //Keep CTS info in fntCTS
  
  //Histograms  
  TNtuple *    fntEMCAL ; //ntuple of EMCAL particles before selection.
  TNtuple *    fntPHOS ; //ntuple of PHOS particles before selection.
  TNtuple *    fntCTS ; //ntuple of CTS particles before selection.

  ClassDef(AliAnaGammaSelection,1)
} ;
 

#endif //ALIANAGAMMASELECTION_H



