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
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for the analysis of gamma correlations  
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TParticle.h>
#include <TH2.h>

//---- AliRoot system ----
#include "AliAnaGammaCorrelation.h" 
#include "AliNeutralMesonSelection.h" 
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliAnaGammaCorrelation)


//____________________________________________________________________________
  AliAnaGammaCorrelation::AliAnaGammaCorrelation() : 
    TObject(), fOutputContainer(0x0),   
    fNeutralMesonSelection(0x0), fCorrelationType(0),
    fJetsOnlyInCTS(0),
    fMinPtHadron(0),
    fDeltaPhiMaxCut(0.), fDeltaPhiMinCut(0.), 
    fRatioMaxCut(0.), fRatioMinCut(0.)
{
  //Default Ctor

  //Initialize parameters

  if(!fNeutralMesonSelection)
    fNeutralMesonSelection = new AliNeutralMesonSelection();
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliAnaGammaCorrelation::AliAnaGammaCorrelation(const AliAnaGammaCorrelation & g) :   
  TObject(), fOutputContainer(g.fOutputContainer),  
  fNeutralMesonSelection(g.fNeutralMesonSelection),
  fCorrelationType(g.fCorrelationType),
  fJetsOnlyInCTS(g.fJetsOnlyInCTS), 
  fMinPtHadron(g.fMinPtHadron),
  fDeltaPhiMaxCut(g.fDeltaPhiMaxCut), fDeltaPhiMinCut(g.fDeltaPhiMinCut), 
  fRatioMaxCut(g.fRatioMaxCut), fRatioMinCut(g.fRatioMinCut) 
{
  // cpy ctor

}

//_________________________________________________________________________
AliAnaGammaCorrelation & AliAnaGammaCorrelation::operator = (const AliAnaGammaCorrelation & source)
{
  // assignment operator

  if(this == &source)return *this;
  ((TObject *)this)->operator=(source);

  fOutputContainer = source.fOutputContainer;  
  fNeutralMesonSelection = source.fNeutralMesonSelection ;
  fCorrelationType = source.fCorrelationType;
  fJetsOnlyInCTS = source.fJetsOnlyInCTS ;

  fMinPtHadron = source.fMinPtHadron ;
  fDeltaPhiMaxCut = source.fDeltaPhiMaxCut ; fDeltaPhiMinCut = source.fDeltaPhiMinCut ; 
  fRatioMaxCut = source.fRatioMaxCut ; fRatioMinCut = source.fRatioMinCut ; 



  return *this;

}

//____________________________________________________________________________
AliAnaGammaCorrelation::~AliAnaGammaCorrelation() 
{
   // Remove all pointers except analysis output pointers.
 
}

 //____________________________________________________________________________
void AliAnaGammaCorrelation::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  fCorrelationType = kHadron ;
  //-----------kHadron----------------
  fMinPtHadron = 0.   ;
  //-----------kHadron & kJetLeadCone----------------
  fJetsOnlyInCTS = kFALSE ;
  fDeltaPhiMaxCut      = 4.5;
  fDeltaPhiMinCut      = 1.5 ;
  //-----------kJetLeadCone----------------
  fRatioMaxCut = 1.0 ; 
  fRatioMinCut = 0.1 ; 

}

//__________________________________________________________________
void AliAnaGammaCorrelation::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Correlation           =     %d\n", fCorrelationType) ;  
  
  printf("pT Hadron       >    %f\n", fMinPtHadron) ; 
  printf("Phi gamma-Hadron      <     %f\n", fDeltaPhiMaxCut) ; 
  printf("Phi gamma-Hadron      >     %f\n", fDeltaPhiMinCut) ;
  printf("Ratio pt hadron/gamma      <     %f\n", fRatioMaxCut) ; 
  printf("Ratio pt hadron/gamma      >     %f\n", fRatioMinCut) ;

} 
