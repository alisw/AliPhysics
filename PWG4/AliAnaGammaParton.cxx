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
 * Revision 1.3  2007/09/26 11:07:19  schutz
 * Update classes for the new analysis framwork
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */
//_________________________________________________________________________
// Class for the analysis of gamma-parton correlations
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "Riostream.h"

//---- AliRoot system ----
#include "AliAnaGammaParton.h" 
#include "AliLog.h"
  
  ClassImp(AliAnaGammaParton)
  

//____________________________________________________________________________
  AliAnaGammaParton::AliAnaGammaParton() : 
    AliAnaGammaCorrelation(),   
    fhDeltaEtaParton(0), fhDeltaPhiParton(0), 
    fhDeltaPtParton(0), fhPtRatParton(0)
{
  //Default Ctor

  SetCorrelationType(kParton);
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaGammaParton::AliAnaGammaParton(const AliAnaGammaParton & g) :   
  AliAnaGammaCorrelation(g),   
  fhDeltaEtaParton(g.fhDeltaEtaParton), fhDeltaPhiParton(g.fhDeltaPhiParton), 
  fhDeltaPtParton(g.fhDeltaPtParton), fhPtRatParton(g.fhPtRatParton)
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaGammaParton & AliAnaGammaParton::operator = (const AliAnaGammaParton & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((AliAnaGammaCorrelation *)this)->operator=(source);
  fhDeltaEtaParton = source.fhDeltaEtaParton;
  fhDeltaPhiParton = source.fhDeltaPhiParton;
  fhDeltaPtParton = source.fhDeltaPtParton;
  fhPtRatParton = source.fhPtRatParton;
  
  return *this;

}

//____________________________________________________________________________
AliAnaGammaParton::~AliAnaGammaParton() 
{
  // Remove all pointers except analysis output pointers.

 
}


//________________________________________________________________________
TList *  AliAnaGammaParton::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer

  AliDebug(1,"Init parton histograms");

  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GammaPartonHistos") ; 

  //---kParton---
  fhDeltaPhiParton  = new TH2F
    ("DeltaPhiParton","#phi_{#gamma} - #phi_{parton} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiParton->SetYTitle("#Delta #phi");
  fhDeltaPhiParton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhDeltaPhiParton);

  fhDeltaEtaParton  = new TH2F
    ("DeltaEtaParton","#eta_{#gamma} - #eta_{parton} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaParton->SetYTitle("#Delta #eta");
  fhDeltaEtaParton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhDeltaEtaParton);

  fhDeltaPtParton  = new TH2F
    ("DeltaPtParton","#p_{T #gamma} - #p_{T parton} vs p_{T #gamma}",
     200,0,120,100,-10,10); 
  fhDeltaPtParton->SetYTitle("#Delta #p_{T}");
  fhDeltaPtParton->SetXTitle("p_{T #gamma} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtParton);

  fhPtRatParton  = new TH2F
    ("PtRatParton","#p_{T parton} / #p_{T #gamma} vs p_{T #gamma}",
     200,0,120,200,0,5); 
  fhPtRatParton->SetYTitle("ratio");
  fhPtRatParton->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPtRatParton);

  SetOutputContainer(outputContainer);

  return outputContainer;
}

 //____________________________________________________________________________
void AliAnaGammaParton::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  ;

}

//__________________________________________________________________
void AliAnaGammaParton::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
} 

//__________________________________________________________________
void  AliAnaGammaParton::MakeGammaCorrelation(TParticle * pGamma, TClonesArray *pl, TClonesArray *) 
{
  //Gamma Parton Correlation Analysis
  AliDebug(2, "Begin parton analysis");
  
  Double_t ptg  = pGamma->Pt();
  Double_t phig = pGamma->Phi();
  Double_t etag = pGamma->Eta();
  
  Double_t pt    = -100.;
  Double_t eta   = -100.; 
  Double_t phi   = -100. ;
  for(Int_t ipr = 0;ipr < pl->GetEntries() ; ipr ++ ){//fCaloList==parton list
    TParticle * particle = dynamic_cast<TParticle *>(pl->At(ipr)) ;
    
    if(particle->GetPdgCode() !=22 && (ipr ==4 || ipr == 5)){// 6 or 7 in list.
      //Only good for gamma-jet events
      pt    = particle->Pt();
      phi   = particle->Phi() ;
      eta   = particle->Eta() ;
    }
  }
  
  fhDeltaEtaParton->Fill(ptg,etag-eta);
  fhDeltaPhiParton->Fill(ptg,phig-phi);
  fhDeltaPtParton->Fill(ptg,ptg-pt);
  fhPtRatParton->Fill(ptg,pt/ptg);
  
  AliDebug(2, "End of parton analysis");
  
} 
