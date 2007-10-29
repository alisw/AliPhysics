/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
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
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TParticle.h>
#include <TH2.h>
#include <TList.h>
#include "Riostream.h"
#include "TROOT.h"

// --- AliRoot system --- 
#include "AliAnaGammaSelection.h" 
#include "AliLog.h"

ClassImp(AliAnaGammaSelection)
  
//____________________________________________________________________________
  AliAnaGammaSelection::AliAnaGammaSelection() : 
    TObject(), fAnaMC(0), fFillCTS(0), 
    fntEMCAL(0), fntPHOS(0), fntCTS(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaGammaSelection::AliAnaGammaSelection(const AliAnaGammaSelection & g) : 
  TObject(g), fAnaMC(g.fAnaMC), fFillCTS(g.fFillCTS), 
  fntEMCAL(g.fntEMCAL),  fntPHOS(g.fntPHOS),  fntCTS(g.fntCTS)
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaGammaSelection & AliAnaGammaSelection::operator = (const AliAnaGammaSelection & source)
{
  // assignment operator
  
  if(&source == this) return *this;

  fAnaMC = source.fAnaMC ;
  fFillCTS = source.fFillCTS ;
  fntEMCAL = source.fntEMCAL ;  
  fntPHOS = source.fntPHOS ;
  fntCTS = source.fntCTS ;

  return *this;
  
}

//____________________________________________________________________________
AliAnaGammaSelection::~AliAnaGammaSelection() 
{
  // Remove all pointers
  
  delete fntEMCAL    ;  
  delete fntPHOS    ;  
  delete fntCTS    ;  

}

//________________________________________________________________________
TList *  AliAnaGammaSelection::GetCreateOutputObjects()
{  

  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("SelectionHistos") ; 
  
//   //Histograms of highest gamma identified in Event
//   fhNGamma  = new TH1F("NGamma","Number of #gamma over calorimeter",240,0,120); 
//   fhNGamma->SetYTitle("N");
//   fhNGamma->SetXTitle("p_{T #gamma}(GeV/c)");
//   outputContainer->Add(fhNGamma) ; 
  
 
  //NTUPLE
  fntEMCAL = new TNtuple("ntEMCAL", "Tree of EMCAL particles before selection", "pt:phi:eta:pdg:status:ptprimary:phiprimary:etaprimary:pdgprimary:statusprimary");
  outputContainer->Add(fntEMCAL) ;
  fntPHOS = new TNtuple("ntPHOS", "Tree of PHOS particles before selection", "pt:phi:eta:pdg:status:ptprimary:phiprimary:etaprimary:pdgprimary:statusprimary");
  //outputContainer->Add(fntPHOS) ;
  fntCTS = new TNtuple("ntCTS", "Tree of CTS particles before selection", "pt:phi:eta:pdg:status:ptprimary:phiprimary:etaprimary:pdgprimary:statusprimary");
  //outputContainer->Add(fntCTS) ;
  
  gROOT->cd();

  return outputContainer ;

}

//____________________________________________________________________________
void AliAnaGammaSelection::Selection(TString det, TClonesArray * pl, TClonesArray * plPrim) const 
{
  //Select particles from detector "det"

  //Keep some particle info before selection
  Float_t pt = 0, ptprimary = 0 ;
  Float_t phi = 0, phiprimary = 0 ;
  Float_t eta = 0, etaprimary = 0 ;
  Int_t pdg = 0, pdgprimary = 0 ;
  Int_t status = 0, statusprimary = 0 ;

  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;
    pt = particle->Pt();
    phi = particle->Phi();
    eta = particle->Eta();
    pdg =  TMath::Abs(particle->GetPdgCode()) ;
    status = particle->GetStatusCode(); // = 2 means decay gamma!!!
    
    if(fAnaMC){
      TParticle * primary = dynamic_cast<TParticle *>(plPrim->At(ipr)) ;
      ptprimary = primary->Pt();
      phiprimary = primary->Phi();
      etaprimary = primary->Eta();
      pdgprimary =  TMath::Abs(primary->GetPdgCode()) ;
      statusprimary = primary->GetStatusCode(); // = 2 means decay gamma!!!
    }

    gROOT->cd();
    if(det == "EMCAL")
    fntEMCAL->Fill(pt,phi,eta,pdg,status,ptprimary,phiprimary, etaprimary,pdgprimary,statusprimary);
    else if(det == "PHOS")
      fntPHOS->Fill(pt,phi,eta,pdg,status,ptprimary,phiprimary, etaprimary,pdgprimary,statusprimary);
    else
      fntCTS->Fill(pt,phi,eta,pdg,status,ptprimary,phiprimary, etaprimary,pdgprimary,statusprimary);
  }
  
}

  //____________________________________________________________________________
void AliAnaGammaSelection::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  fAnaMC = kFALSE ;
 
  fFillCTS = kFALSE ;
}


//__________________________________________________________________
void AliAnaGammaSelection::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ; 
  printf("Is MC option selected     %d\n",  fAnaMC) ;
  printf("    \n") ;
  
} 
