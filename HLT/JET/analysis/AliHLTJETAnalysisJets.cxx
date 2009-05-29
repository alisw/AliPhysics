//-*- Mode: C++ -*-
// $Id: AliHLTJETAnalysisJets.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTJETAnalysisJets.cxx
    @author Jochen Thaeder
    @date   
    @brief  Container holding analysis objects
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TH2F.h"

#include "AliHLTJETAnalysisJets.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAnalysisJets)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTJETAnalysisJets::AliHLTJETAnalysisJets() :
  fJets(NULL ),
  fJetsMC(NULL ),
  fMatchedJets(NULL),
  fMatchedJetsMC(NULL),
  fDeltaEt(NULL),
  fDeltaEta(NULL),
  fDeltaPhi(NULL),
  fDeltaEtaDeltaPhi(NULL),
  fSpectraEt(NULL),
  fSpectraEta(NULL),
  fSpectraPhi(NULL),
  fCorrelationsJetEt(NULL),
  fResolutionsJetEt(NULL),
  fResolutionsDiJetEt(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTJETAnalysisJets::~AliHLTJETAnalysisJets() {
  // see header file for class documentation

  if ( fJetsMC )
    delete fJetsMC;
  fJetsMC = NULL;

  if ( fMatchedJets )
    delete fMatchedJets;
  fMatchedJets = NULL;

  if ( fMatchedJetsMC )
    delete fMatchedJetsMC;
  fMatchedJetsMC = NULL;

  if ( fDeltaEt ) {
    fDeltaEt->Clear();
    delete fDeltaEt;
  }
  fDeltaEt = NULL;

  if ( fDeltaEta ) {
    fDeltaEta->Clear();
    delete fDeltaEta;
  }
  fDeltaEta = NULL;

  if ( fDeltaPhi ) {
    fDeltaPhi->Clear();
    delete fDeltaPhi;
  }
  fDeltaPhi = NULL;
  
  if ( fDeltaEtaDeltaPhi ) {
    fDeltaEtaDeltaPhi->Clear();
    delete fDeltaEtaDeltaPhi;
  }
  fDeltaEtaDeltaPhi = NULL;

  if ( fSpectraEt ) {
    fSpectraEt->Clear();
    delete fSpectraEt;
  }
  fSpectraEt = NULL;

  if ( fSpectraEta ) {
    fSpectraEta->Clear();
    delete fSpectraEta;
  }
  fSpectraEta = NULL;

  if ( fSpectraPhi ) {
    fSpectraPhi->Clear();
    delete fSpectraPhi;
  }
  fSpectraPhi = NULL;

  if ( fCorrelationsJetEt ) {
    fCorrelationsJetEt->Clear();
    delete fCorrelationsJetEt;   
  }
  fCorrelationsJetEt = NULL;

  if ( fResolutionsJetEt ) {
    fResolutionsJetEt->Clear();
    delete fResolutionsJetEt;   
  }
  fResolutionsJetEt = NULL;

  if ( fResolutionsDiJetEt ) {
    fResolutionsDiJetEt->Clear();
    delete fResolutionsDiJetEt;   
  }
  fResolutionsDiJetEt = NULL;

}

/*
 * ---------------------------------------------------------------------------------
 *                                   Initialize / Reset
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETAnalysisJets::Initialize() {
  // see header file for class documentation  

  Int_t iResult = 0;

  // -- Setup match arrays
  fMatchedJets = new TArrayI(100);
  fMatchedJetsMC = new TArrayI(100);

  // -- Setup Delta histograms  
  SetupDeltaHistograms();

  // -- Setup Spectra histograms
  SetupSpectraHistograms();

  // -- Setup Matched histograms
  SetupMatchedHistograms();

  return iResult;
}

//##################################################################################
void AliHLTJETAnalysisJets::ResetEvent() {
  // see header file for class documentation  

  if ( fJetsMC )
    delete fJetsMC;
  fJetsMC = NULL;

  fJets = NULL;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Setter - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisJets::SetHLTMC( AliHLTMCEvent* mcEvent ) {
  // see header file for class documentation

  // -- New MC jets 
  if ( fJetsMC )
    delete fJetsMC;
  fJetsMC = new AliHLTJETJets();
  
  AliAODJet* jet = NULL;

  while ( (jet = mcEvent->NextGenJet()) )
    fJetsMC->AddJet(jet);

  // -- Sort jets
  fJetsMC->Sort();

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::SetMC( AliMCEvent* /*mcEvent*/ ) {
  // see header file for class documentation

  HLTFatal("No implemented!");
  
  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::SetJets( AliHLTJETJets* jets ) {
  // see header file for class documentation

  fJets = jets;

  // -- Sort jets
  fJets->Sort();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Getter - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
TH1* AliHLTJETAnalysisJets::GetHistogram( Int_t histIdx, Int_t plotIdx ) {
  // see header file for class documentation
  
  if ( histIdx >= AliHLTJETAnalysisBase::kHistMax ) {
    HLTError("Histogram index out of bound : %d - max = %d", histIdx, AliHLTJETAnalysisBase::kHistMax );
    return NULL; 
  }
  
  TClonesArray *hist = NULL;
  Int_t maxIdx = -1;
  
  switch (histIdx ) {
  case AliHLTJETAnalysisBase::kHistDeltaEt : 
    hist = fDeltaEt;
    maxIdx = AliHLTJETAnalysisBase::kDeltaMax;
    break;
  case AliHLTJETAnalysisBase::kHistDeltaEta : 
    hist = fDeltaEta;
    maxIdx = AliHLTJETAnalysisBase::kDeltaMax;
    break;
  case AliHLTJETAnalysisBase::kHistDeltaPhi : 
    hist = fDeltaPhi;
    maxIdx = AliHLTJETAnalysisBase::kDeltaMax;
    break;
  case AliHLTJETAnalysisBase::kHistDeltaEtaDeltaPhi : 
    hist = fDeltaEtaDeltaPhi;
    maxIdx = AliHLTJETAnalysisBase::kDeltaMax;
    break;
  case AliHLTJETAnalysisBase::kHistSpectraEt :
    hist = fSpectraEt;
    maxIdx = AliHLTJETAnalysisBase::kSpectraMax;
    break;
  case AliHLTJETAnalysisBase::kHistSpectraEta :
    hist = fSpectraEta;
    maxIdx = AliHLTJETAnalysisBase::kSpectraMax;
    break;
  case AliHLTJETAnalysisBase::kHistSpectraPhi :
    hist = fSpectraPhi;
    maxIdx = AliHLTJETAnalysisBase::kSpectraMax;
    break;
  case AliHLTJETAnalysisBase::kHistCorrelationsJetEt :
    hist = fCorrelationsJetEt;
    maxIdx = AliHLTJETAnalysisBase::kPlotMax;
    break;
  case AliHLTJETAnalysisBase::kHistResolutionsJetEt :
    hist = fResolutionsJetEt;
    maxIdx = AliHLTJETAnalysisBase::kPlotMax;
    break;
  case AliHLTJETAnalysisBase::kHistResolutionsDiJetEt :
    hist = fResolutionsDiJetEt;
    maxIdx = AliHLTJETAnalysisBase::kPlotMax;
    break;
  }
 
  // -- Check Boundaries
  if ( plotIdx >= maxIdx ) {
    HLTError("Index out of bound : %d - max = %d", plotIdx, maxIdx );
    return NULL; 
  }
  else if ( maxIdx == -1 ) {
    HLTError("Histogram index not found." );
    return NULL; 
  }

  // -- Retrieve histogram
  return reinterpret_cast<TH1*>((*hist)[plotIdx]);  
}

/*
 * ---------------------------------------------------------------------------------
 *                              Analysis - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETAnalysisJets::Analyze() {
  // see header file for class documentation
  
  Int_t iResult = 0;

  if ( !fJets ) {
    HLTError("No input jets set.");
    iResult = -1;
  }
  if ( !fJetsMC ) {
    HLTError("No input MC jets set.");
    iResult = -1;
  }

  if ( !iResult) {    
    
    // -- Fill unmatched jets into histograms
    FillBasicSpectraHistograms();
    FillUnmatchedDeltaHistograms();
    
    // -- Match jets
    MatchJets();
   
    // -- Fill matched jets into histograms
    FillMatchedDeltaHistograms();  
    FillMatchedSpectraHistograms();  
    FillMatchedHistograms();  
  }

  return iResult;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                              ///
//////                                               PRIVATE                                                  //////
///                                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * ---------------------------------------------------------------------------------
 *                             Setup / Reset - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisJets::SetupDeltaHistograms() {
  // see header file for class documentation
  
  // ---------------------------------------------------
  // -- Difference in reconstruction
  //    All / Matched jets
  // ---------------------------------------------------

  fDeltaEt                 = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kDeltaMax );
  fDeltaEta                = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kDeltaMax );
  fDeltaPhi                = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kDeltaMax );
  fDeltaEtaDeltaPhi        = new TClonesArray( "TH2F", AliHLTJETAnalysisBase::kDeltaMax );

  for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kDeltaMax; ++idx ) {

    // -- Delta Et -------------------------------------------------
    new ((*fDeltaEt)[idx]) TH1F(Form("delta E_{t} : %s", 
				     AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
				Form("#DeltaE_{t} : %s;#DeltaE_{t};dN/d#DeltaE_{t}", 
				     AliHLTJETAnalysisBase::fgkDeltaType[idx] ),
				100, -200., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaEt)[idx]));
        
    // -- Delta Eta ------------------------------------------------
    new ((*fDeltaEta)[idx]) TH1F(Form("delta Eta : %s", 
				      AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
				 Form("#Delta#eta : %s;#Delta#eta;dN/d#Delta#eta", 
				      AliHLTJETAnalysisBase::fgkDeltaType[idx] ),
				 100, -1.2, 1.2);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaEta)[idx]));
    
    // -- Delta Phi ------------------------------------------------
    new ((*fDeltaPhi )[idx]) TH1F(Form("delta Phi : %s", 
				       AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
				  Form("#Delta#phi : %s;#Delta #phi;dN/d#Delta#phi", 
				       AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
				  100, -7., 7.);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaPhi)[idx]));
    
    // -- Delta Eta Delta Phi --------------------------------------
    new ((*fDeltaEtaDeltaPhi) [idx] ) TH2F(Form("delta Eta delta Phi : %s", 
						AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
					   Form("#Delta#eta #Delta#phi : %s;#Delta#eta;#Delta#phi", 
						AliHLTJETAnalysisBase::fgkDeltaType[idx] ), 
					   100, -1.2, 1.2, 100, -7., 7.);
    SetupHist(reinterpret_cast<TH2F*>((*fDeltaEtaDeltaPhi)[idx]));

  } //  for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kDeltaMax; ++idx ) {

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::SetupSpectraHistograms() {
  // see header file for class documentation
    
  // ---------------------------------------------------
  // -- Jet spectra
  // ---------------------------------------------------

  fSpectraEt  = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kSpectraMax );
  fSpectraEta = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kSpectraMax );
  fSpectraPhi = new TClonesArray( "TH1F", AliHLTJETAnalysisBase::kSpectraMax );

  for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kSpectraMax; ++idx ) {
    
    // -- Spectra Et -----------------------------------------------
    new ((*fSpectraEt)[idx]) TH1F(Form("E_{t} : %s", 
				       AliHLTJETAnalysisBase::fgkSpectraType[idx] ), 
				  Form("E_{t} : %s;E_{t} (GeV/c);dN/dE_{t}", 
				       AliHLTJETAnalysisBase::fgkSpectraType[idx] ),
				  100, 0., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fSpectraEt)[idx]));
    
    // -- Spectra Eta ----------------------------------------------
    new ((*fSpectraEta)[idx]) TH1F(Form("#eta : %s", 
					AliHLTJETAnalysisBase::fgkSpectraType[idx] ), 
				   Form("#eta : %s;#eta;dN/d#eta",
					AliHLTJETAnalysisBase::fgkSpectraType[idx] ),
				   80, -0.9, 0.9);
    SetupHist(reinterpret_cast<TH1F*>((*fSpectraEta)[idx]));
    
    // -- Spectra Phi ----------------------------------------------
    new ((*fSpectraPhi)[idx]) TH1F(Form("#phi : %s", 
					AliHLTJETAnalysisBase::fgkSpectraType[idx] ), 
				   Form("#phi : %s;#phi;dN/d#phi",
					AliHLTJETAnalysisBase::fgkSpectraType[idx] ),
				   50, 0., 7.);
    SetupHist(reinterpret_cast<TH1F*>((*fSpectraPhi)[idx]));
    
  } // for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kSpectraMax; ++idx ) {

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::SetupMatchedHistograms() {
  // see header file for class documentation
  
  // ---------------------------------------------------
  // -- Correlations
  // -- Resolutions
  // ---------------------------------------------------
  
  fCorrelationsJetEt  = new TClonesArray( "TH2F", AliHLTJETAnalysisBase::kPlotMax );
  fResolutionsJetEt   = new TClonesArray( "TH2F", AliHLTJETAnalysisBase::kPlotMax );
  fResolutionsDiJetEt = new TClonesArray( "TH2F", AliHLTJETAnalysisBase::kPlotMax );

  for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kPlotMax; ++idx ) {

    // -- Correlations ---------------------------------------------
    new ((*fCorrelationsJetEt)[idx]) TH2F(Form("Correlations : %s", 
					       AliHLTJETAnalysisBase::fgkPlotType[idx] ), 
					  Form("Correlations : %s; E_{t,Pythia} (GeV/c); E_{t,Rec} (GeV/c)", 
					       AliHLTJETAnalysisBase::fgkPlotType[idx] ),
					  100, 0., 200.,100, 0., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fCorrelationsJetEt)[idx]));
    
    // -- Jet Resolutions ------------------------------------------
    new ((*fResolutionsJetEt)[idx]) TH2F(Form("Resolutions : %s", 
					      AliHLTJETAnalysisBase::fgkPlotType[idx] ), 
					 Form("Resolutions : %s; E_{t,Pythia} (GeV/c); f", 
					      AliHLTJETAnalysisBase::fgkPlotType[idx] ),
					 200, 0., 200.,200, -2., 2.);
					
    SetupHist(reinterpret_cast<TH1F*>((*fResolutionsJetEt)[idx]));

    // -- Di-Jet Resolutions ---------------------------------------
    new ((*fResolutionsDiJetEt)[idx]) TH2F(Form("Di Jet Resolutions : %s", 
						AliHLTJETAnalysisBase::fgkPlotType[idx] ), 
					   Form("Di Jet Resolutions : %s; E_{t,Pythia} (GeV/c); f", 
						AliHLTJETAnalysisBase::fgkPlotType[idx] ),
					   100, -200., 200.,100, -200., 200.);
					
    SetupHist(reinterpret_cast<TH1F*>((*fResolutionsDiJetEt)[idx]));

  } // for ( Int_t idx = 0; idx < kJetPlotMax; ++idx ) {

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Analysis - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETAnalysisJets::MatchJets() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Reset match arrays
  Int_t maxNJets;
  if ( fJets->GetNAODJets() > fJetsMC->GetNAODJets() )
    maxNJets = fJets->GetNAODJets();
  else
    maxNJets = fJetsMC->GetNAODJets();

  for ( Int_t idx = 0; idx < maxNJets; ++idx ) {
    (*fMatchedJets)[idx] = -1;
    (*fMatchedJetsMC)[idx] = -1;
  }

  // -- Match Jets - pythia fixed
  for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {

    Float_t minDistance2 = 100.;
    Int_t idxClosest = -1;

    for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {

      // -- very Jet matched only once
      if ( (*fMatchedJets)[jetIter] != -1 )
	continue;
                          
      Float_t distance2 = GetDistance2( fJetsMC->GetJet(jetIterMC), 
					fJets->GetJet(jetIter) );
      
      // -- Get Closest
      if ( distance2 < minDistance2 ) {
	minDistance2 = distance2;
	idxClosest = jetIter;
      }

    } //     for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {
 
    // -- Match found
    if  ( minDistance2 < 10. ) {
      (*fMatchedJetsMC)[jetIterMC] = idxClosest;
      (*fMatchedJets)[idxClosest] = jetIterMC;
    }

  } //  for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Fill - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisJets::FillBasicSpectraHistograms() {
  // see header file for class documentation

  // -- Fill basic MC spectras
  // ---------------------------
  for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {
    FillHist(fSpectraEt, AliHLTJETAnalysisBase::kSpectraPythiaAll, fJetsMC->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraPythiaAll, fJetsMC->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraPythiaAll, fJetsMC->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {

  // -- Fill basic Reco spectras
  // -----------------------------
  for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {
    FillHist(fSpectraEt, AliHLTJETAnalysisBase::kSpectraRecoAll, fJets->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraRecoAll, fJets->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraRecoAll, fJets->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {

  // -- Fill basic Reco leading spectras
  // -------------------------------------
  if ( fJets->GetNAODJets() > 0 ) {
    FillHist(fSpectraEt, AliHLTJETAnalysisBase::kSpectraRecoLeadAll, fJets->GetJet(0)->Pt());
    FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraRecoLeadAll, fJets->GetJet(0)->Eta());
    FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraRecoLeadAll, fJets->GetJet(0)->Phi());
  } 

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillUnmatchedDeltaHistograms() {
  // see header file for class documentation

  for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {
    AliAODJet *jetMC = fJetsMC->GetJet(jetIterMC);
    
    for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {
      AliAODJet *jet = fJets->GetJet(jetIter);

      FillHist(fDeltaEt, AliHLTJETAnalysisBase::kDeltaAll, jetMC->Pt()-jet->Pt());
      FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaAll, jetMC->Eta()-jet->Eta());
      FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaAll, jetMC->Phi()-jet->Phi());
      FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaAll,
	       jetMC->Eta()-jet->Eta(), jetMC->Phi()-jet->Phi());
      
    } // for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {
  } // for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {

  // -- Leading Jets
  // -----------------
  if ( fJets->GetNAODJets() > 0 && fJetsMC->GetNAODJets() > 0 ) {
    AliAODJet *jetMC = fJetsMC->GetJet(0);
    AliAODJet *jet = fJets->GetJet(0);

    FillHist(fDeltaEt, AliHLTJETAnalysisBase::kDeltaLead, jetMC->Pt()-jet->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaLead, jetMC->Eta()-jet->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaLead, jetMC->Phi()-jet->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaLead,
	     jetMC->Eta()-jet->Eta(), jetMC->Phi()-jet->Phi());
  }

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedDeltaHistograms() {
  // see header file for class documentation

  for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {

    // -- Not matched jet
    if ( (*fMatchedJetsMC)[jetIterMC] == -1 )
      continue;
    
    AliAODJet *jetMC = fJetsMC->GetJet(jetIterMC);
    AliAODJet *jet = fJets->GetJet((*fMatchedJetsMC)[jetIterMC]);

    FillHist(fDeltaEt, AliHLTJETAnalysisBase::kDeltaMatchedAll, jetMC->Pt()-jet->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaMatchedAll, jetMC->Eta()-jet->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedAll, jetMC->Phi()-jet->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedAll,
	     jetMC->Eta()-jet->Eta(), jetMC->Phi()-jet->Phi());
    
  } // for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {
  
  // -- Leading Jets
  // -----------------
  if ((*fMatchedJetsMC)[0] != -1 ) {

    AliAODJet *jetMC = fJetsMC->GetJet(0);
    AliAODJet *jet = fJets->GetJet((*fMatchedJetsMC)[0]);
    
    FillHist(fDeltaEt, AliHLTJETAnalysisBase::kDeltaMatchedLead, jetMC->Pt()-jet->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaMatchedLead, jetMC->Eta()-jet->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedLead, jetMC->Phi()-jet->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedLead,
	     jetMC->Eta()-jet->Eta(), jetMC->Phi()-jet->Phi());
  }
  
  return;
}
 
//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedSpectraHistograms() {
  // see header file for class documentation

  Int_t idx = 0;

  // -- Fill matched MC spectras
  // ---------------------------
  for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {
    if ( (*fMatchedJetsMC)[jetIter] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraPythiaUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraPythiaMatched;

    FillHist(fSpectraEt, idx, fJetsMC->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, idx, fJetsMC->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, idx, fJetsMC->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {

  // -- Fill matched Reco spectras
  // -----------------------------
  for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {
    if ( (*fMatchedJets)[jetIter] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraRecoUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraRecoMatched;
    
    FillHist(fSpectraEt, idx, fJets->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, idx, fJets->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, idx, fJets->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJets->GetNAODJets(); ++jetIter ) {

  // -- Fill matched Reco leading spectras
  // -------------------------------------
  if ( fJets->GetNAODJets() > 0 ) {
    if ( (*fMatchedJets)[0] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraRecoLeadUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraRecoLeadMatched;

    FillHist(fSpectraEt, idx, fJets->GetJet(0)->Pt());
    FillHist(fSpectraEta, idx, fJets->GetJet(0)->Eta());
    FillHist(fSpectraPhi, idx, fJets->GetJet(0)->Phi());
  } 
  
  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedHistograms() {
  // see header file for class documentation

  for ( Int_t jetIterMC = 0; jetIterMC < fJetsMC->GetNAODJets(); ++jetIterMC ) {

    // -- Not matched jet
    if ( (*fMatchedJetsMC)[jetIterMC] == -1 )
      continue;
    
    AliAODJet *jetMC = fJetsMC->GetJet(jetIterMC);
    AliAODJet *jet = fJets->GetJet((*fMatchedJetsMC)[jetIterMC]);

    // -- Correlations    
    FillHist(fCorrelationsJetEt, AliHLTJETAnalysisBase::kPlotAll, jetMC->Pt(), jet->Pt());

    // -- Resolutions
    FillHist(fResolutionsJetEt, AliHLTJETAnalysisBase::kPlotAll, jetMC->Pt(),
	     (jetMC->Pt()-jet->Pt())/(jetMC->Pt()+jet->Pt()));

  } // for ( Int_t jetIter = 0; jetIter < fJetsMC->GetNAODJets(); ++jetIter ) {
  
  // -- Leading Jets
  // -----------------
  if ((*fMatchedJetsMC)[0] != -1 ) {

    AliAODJet *jetMC = fJetsMC->GetJet(0);
    AliAODJet *jet = fJets->GetJet((*fMatchedJetsMC)[0]);
    
    // -- Correlations    
    FillHist(fCorrelationsJetEt, AliHLTJETAnalysisBase::kPlotLead, jetMC->Pt(), jet->Pt());
    
    // -- Resolutions
    FillHist(fResolutionsJetEt, AliHLTJETAnalysisBase::kPlotLead, jetMC->Pt(),
	     (jetMC->Pt()-jet->Pt())/(jetMC->Pt()+jet->Pt()));
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                               Helper - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisJets::SetupHist( TH1* hist ) {
  // see header file for class documentation
  
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->SetMarkerStyle(kFullCircle);
  //  hist->SetStats(kFALSE);

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillHist( TClonesArray* array, Int_t idx, Float_t valueX ) {
    // see header file for class documentation
  
  reinterpret_cast<TH1F*>((*array)[idx])->Fill(valueX);
  
  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillHist( TClonesArray* array, Int_t idx, Float_t valueX, Float_t valueY ) {
  // see header file for class documentation
  
  reinterpret_cast<TH2F*>((*array)[idx])->Fill(valueX, valueY);
  
  return;
}

//################################################################################## 
Float_t AliHLTJETAnalysisJets::GetDistance2( AliAODJet *jet1, AliAODJet *jet2 ) {
  // see header file for class documentation
  
  return ( (jet1->Eta()-jet2->Eta())*(jet1->Eta()-jet2->Eta()) ) + 
    ( (jet1->Phi()-jet2->Phi())*(jet1->Phi()-jet2->Phi()) );
}
