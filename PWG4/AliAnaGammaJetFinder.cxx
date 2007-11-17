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
 * Revision 1.3  2007/10/29 13:48:42  gustavo
 * Corrected coding violations
 *
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for the analysis of gamma-jet (standard jet finder) correlations
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "Riostream.h"

//---- AliRoot system ----
#include "AliAnaGammaJetFinder.h" 
#include "AliLog.h"
  
  ClassImp(AliAnaGammaJetFinder)
  

//____________________________________________________________________________
  AliAnaGammaJetFinder::AliAnaGammaJetFinder() : 
    AliAnaGammaCorrelation(),   
    fhDeltaEtaJet(0), fhDeltaPhiJet(0), 
    fhDeltaPtJet(0), fhPtRatJet(0)
{
  //Default Ctor

  SetCorrelationType(kJetFinder);
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaGammaJetFinder::AliAnaGammaJetFinder(const AliAnaGammaJetFinder & g) :   
  AliAnaGammaCorrelation(g),   
  fhDeltaEtaJet(g.fhDeltaEtaJet), fhDeltaPhiJet(g.fhDeltaPhiJet), 
  fhDeltaPtJet(g.fhDeltaPtJet), fhPtRatJet(g.fhPtRatJet)
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaGammaJetFinder & AliAnaGammaJetFinder::operator = (const AliAnaGammaJetFinder & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((AliAnaGammaCorrelation *)this)->operator=(source);
  fhDeltaEtaJet = source.fhDeltaEtaJet;
  fhDeltaPhiJet = source.fhDeltaPhiJet;
  fhDeltaPtJet = source.fhDeltaPtJet;
  fhPtRatJet = source.fhPtRatJet;
  
  return *this;

}

//____________________________________________________________________________
AliAnaGammaJetFinder::~AliAnaGammaJetFinder() 
{
   // Remove all pointers except analysis output pointers.
 
}


//________________________________________________________________________
TList *  AliAnaGammaJetFinder::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer

  AliDebug(1,"Init jet histograms");

  TList * outputContainer = new TList() ; 
  outputContainer->SetName("GammaJetHistos") ; 

  //---kJet---
  fhDeltaPhiJet  = new TH2F
    ("DeltaPhiJet","#phi_{#gamma} - #phi_{jet} vs p_{T #gamma}",
     200,0,120,200,0,6.4); 
  fhDeltaPhiJet->SetYTitle("#Delta #phi");
  fhDeltaPhiJet->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhDeltaPhiJet);

  fhDeltaEtaJet  = new TH2F
    ("DeltaEtaJet","#eta_{#gamma} - #eta_{jet} vs p_{T #gamma}",
     200,0,120,200,-2,2); 
  fhDeltaEtaJet->SetYTitle("#Delta #eta");
  fhDeltaEtaJet->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhDeltaEtaJet);

  fhDeltaPtJet  = new TH2F
    ("DeltaPtJet","#p_{T #gamma} - #p_{T jet} vs p_{T #gamma}",
     200,0,120,100,-10,10); 
  fhDeltaPtJet->SetYTitle("#Delta #p_{T}");
  fhDeltaPtJet->SetXTitle("p_{T #gamma} (GeV/c)"); 
  outputContainer->Add(fhDeltaPtJet);

  fhPtRatJet  = new TH2F
    ("PtRatJet","#p_{T jet} / #p_{T #gamma} vs p_{T #gamma}",
     200,0,120,200,0,5); 
  fhPtRatJet->SetYTitle("ratio");
  fhPtRatJet->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPtRatJet);

  SetOutputContainer(outputContainer);

  return outputContainer;
}

 //____________________________________________________________________________
void AliAnaGammaJetFinder::InitParameters()
{
 
  //Initialize the parameters of the analysis.

  ;

}

//__________________________________________________________________
void AliAnaGammaJetFinder::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
} 

//__________________________________________________________________
void  AliAnaGammaJetFinder::MakeGammaCorrelation(TParticle * pGamma, TClonesArray *pl, TClonesArray *) 
{
  //Gamma -Jet  Correlation Analysis
  AliDebug(2, "Begin jet analysis");
  cout<<pGamma<<" "<<pl<<endl;
  AliInfo("Not implemented");  
} 
