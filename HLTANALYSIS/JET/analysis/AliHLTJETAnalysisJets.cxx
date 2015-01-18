//-*- Mode: C++ -*-
// $Id: AliHLTJETAnalysisJets.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
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
    @author Jochen Thaeder <jochen@thaeder.de>
    @brief  Container holding analysis objects
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TH2F.h"

#include "AliHLTJETAnalysisJets.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAnalysisJets)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTJETAnalysisJets::AliHLTJETAnalysisJets() :
  fHasMC(kFALSE),
  fJetsRec(NULL), fJetsCmp(NULL),
  fMatchedJetsRec(NULL), fMatchedJetsCmp(NULL),
  fMatchingThreshold(5.),
  fDeltaEt(NULL), fDeltaEta(NULL), fDeltaPhi(NULL), fDeltaEtaDeltaPhi(NULL),
  fSpectraEt(NULL), fSpectraEta(NULL), fSpectraPhi(NULL),
  fCorrelationsJetEt(NULL),
  fResolutionsJetEt(NULL), fResolutionsDiJetEt(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTJETAnalysisJets::~AliHLTJETAnalysisJets() {
  // see header file for class documentation

  if ( fHasMC ) {
    if ( fJetsCmp )
      delete fJetsCmp;
    fJetsCmp = NULL;
  }
  
  if ( fMatchedJetsRec )
    delete fMatchedJetsRec;
  fMatchedJetsRec = NULL;

  if ( fMatchedJetsCmp )
    delete fMatchedJetsCmp;
  fMatchedJetsCmp = NULL;

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
  fMatchedJetsRec = new TArrayI(100);
  fMatchedJetsCmp = new TArrayI(100);

  // -- Setup Delta histograms  
  SetupDeltaHistograms();

  // -- Setup Spectra histograms
  SetupSpectraHistograms();

  // -- Setup Matched histograms
  SetupMatchedHistograms();

  fMatchingThreshold = 10.;

  return iResult;
}

//##################################################################################
void AliHLTJETAnalysisJets::ResetEvent() {
  // see header file for class documentation  

  fJetsRec = NULL;

  if ( fJetsCmp && fHasMC )
    delete fJetsCmp;
  fJetsCmp = NULL;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Setter - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisJets::SetJetsRec( AliHLTJets* jets ) {
  // see header file for class documentation

  fJetsRec = jets;

  // -- Sort jets
  fJetsRec->Sort();

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::SetJetsCmp( AliHLTMCEvent* hltMcEvent = NULL, 
					AliMCEvent* mcEvent = NULL, 
					AliHLTJets* jets = NULL ) {
  // see header file for class documentation

  // -- Fill HLT MC event
  if ( hltMcEvent ) {   
    fHasMC = kTRUE;

    // -- New MC jets 
    if ( fJetsCmp )
      delete fJetsCmp;
    fJetsCmp = new AliHLTJets();
  
    AliAODJet* jet = NULL;

    while ( (jet = hltMcEvent->NextGenJet()) )
      fJetsCmp->AddJet(jet);
  }

  // -- Fill Off-Line MC event
  else if ( mcEvent ) {
    fHasMC = kTRUE;
    HLTFatal("No implemented!");
  }

  // -- Fill reconstructed jets
  else  
    fJetsCmp = jets;

  // -- Sort jets
  fJetsCmp->Sort();

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

  if ( !fJetsRec ) {
    HLTError("No input jets set.");
    iResult = -1;
  }

  if ( !iResult) {    
    
    // -- Fill unmatched jets into histograms
    FillBasicSpectraHistograms();
    FillUnmatchedDeltaHistograms();
    
    // -- Match jets and fill matched jets into histograms
    if ( ! MatchJets() ) {
      FillMatchedDeltaHistograms();  
      FillMatchedSpectraHistograms();  
      FillMatchedHistograms();  
    }
  }
  return iResult;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///                                                                              ///
//////                             PRIVATE                                    //////
///                                                                              ///
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

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

    const Char_t *type = AliHLTJETAnalysisBase::fgkDeltaType[idx];

    // -- Delta Et -------------------------------------------------
    new ((*fDeltaEt)[idx]) TH1F(Form("delta E_{t} : %s", type), 
				Form("#DeltaE_{t} : %s;#DeltaE_{t};dN/d#DeltaE_{t}", type),
				100, -200., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaEt)[idx]));
        
    // -- Delta Eta ------------------------------------------------
    new ((*fDeltaEta)[idx]) TH1F(Form("delta Eta : %s", type), 
				 Form("#Delta#eta : %s;#Delta#eta;dN/d#Delta#eta", type), 
				 100, -1.2, 1.2);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaEta)[idx]));
    
    // -- Delta Phi ------------------------------------------------
    new ((*fDeltaPhi )[idx]) TH1F(Form("delta Phi : %s", type), 
				  Form("#Delta#phi : %s;#Delta #phi;dN/d#Delta#phi", type), 
				  100, -7., 7.);
    SetupHist(reinterpret_cast<TH1F*>((*fDeltaPhi)[idx]));
    
    // -- Delta Eta Delta Phi --------------------------------------
    new ((*fDeltaEtaDeltaPhi) [idx] ) TH2F(Form("delta Eta delta Phi : %s", type), 
					   Form("#Delta#eta #Delta#phi : %s;#Delta#eta;#Delta#phi", type), 
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

  const Char_t *type = NULL;

  for ( Int_t idx = 0; idx < AliHLTJETAnalysisBase::kSpectraMax; ++idx ) {
    
    if (fHasMC)
      type = AliHLTJETAnalysisBase::fgkSpectraTypeMC[idx];
    else
      type = AliHLTJETAnalysisBase::fgkSpectraType[idx];

    // -- Spectra Et -----------------------------------------------
    new ((*fSpectraEt)[idx]) TH1F(Form("E_{t} : %s", type), 
				  Form("E_{t} : %s;E_{t} (GeV/c);dN/dE_{t}", type), 
				  100, 0., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fSpectraEt)[idx]));
    
    // -- Spectra Eta ----------------------------------------------
    new ((*fSpectraEta)[idx]) TH1F(Form("#eta : %s", type), 
				   Form("#eta : %s;#eta;dN/d#eta", type), 
				   80, -0.9, 0.9);
    SetupHist(reinterpret_cast<TH1F*>((*fSpectraEta)[idx]));
    
    // -- Spectra Phi ----------------------------------------------
    new ((*fSpectraPhi)[idx]) TH1F(Form("#phi : %s", type), 
				   Form("#phi : %s;#phi;dN/d#phi", type), 
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

    const Char_t *type = AliHLTJETAnalysisBase::fgkPlotType[idx];

    // -- Correlations ---------------------------------------------
    new ((*fCorrelationsJetEt)[idx]) TH2F(Form("Correlations : %s", type), 
					  Form("Correlations : %s; E_{t,Pythia} (GeV/c); E_{t,Rec} (GeV/c)", type), 
					  100, 0., 200.,100, 0., 200.);
    SetupHist(reinterpret_cast<TH1F*>((*fCorrelationsJetEt)[idx]));
    
    // -- Jet Resolutions ------------------------------------------
    new ((*fResolutionsJetEt)[idx]) TH2F(Form("Resolutions : %s", type), 
					 Form("Resolutions : %s; E_{t,Pythia} (GeV/c); f", type), 
					 200, 0., 200.,200, -2., 2.);
					
    SetupHist(reinterpret_cast<TH1F*>((*fResolutionsJetEt)[idx]));

    // -- Di-Jet Resolutions ---------------------------------------
    new ((*fResolutionsDiJetEt)[idx]) TH2F(Form("Di Jet Resolutions : %s", type), 
					   Form("Di Jet Resolutions : %s; E_{t,Pythia} (GeV/c); f", type),  
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

  // -- Check for both inputs
  if ( !fJetsRec || !fJetsCmp )
    return -1;

  // -- Check for jets in both inputs
  if ( fJetsRec->GetNAODJets() == 0 || fJetsCmp->GetNAODJets() == 0 )
    return -1;

  // -- Reset match arrays
  Int_t maxNJets;
  if ( fJetsRec->GetNAODJets() > fJetsCmp->GetNAODJets() )
    maxNJets = fJetsRec->GetNAODJets();
  else
    maxNJets = fJetsCmp->GetNAODJets();

  for ( Int_t idx = 0; idx < maxNJets; ++idx ) {
    (*fMatchedJetsRec)[idx] = -1;
    (*fMatchedJetsCmp)[idx] = -1;
  }

  // -- Match Jets - compare fixed
  for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {

    Float_t minDistance2 = 100.;
    Int_t idxClosest = -1;

    for ( Int_t jetIterRec = 0; jetIterRec < fJetsRec->GetNAODJets(); ++jetIterRec ) {

      // -- every Jet matched only once
      if ( (*fMatchedJetsRec)[jetIterRec] != -1 )
	continue;
                          
      Float_t distance2 = GetDistance2( fJetsCmp->GetJet(jetIterCmp), 
					fJetsRec->GetJet(jetIterRec) );
      
      // -- Get Closest
      if ( distance2 < minDistance2 ) {
	minDistance2 = distance2;
	idxClosest = jetIterRec;
      }

    } //     for ( Int_t jetIterRec = 0; jetIterRec < fJetsRec->GetNAODJets(); ++jetIterRec ) {
 
    // -- Match found
    if  ( minDistance2 < fMatchingThreshold ) {
      (*fMatchedJetsCmp)[jetIterCmp] = idxClosest;
      (*fMatchedJetsRec)[idxClosest] = jetIterCmp;
    }

  } //  for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {

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

  if (fJetsRec) {
    // -- Fill basic Reco spectras
    // -----------------------------
    for ( Int_t jetIter = 0; jetIter < fJetsRec->GetNAODJets(); ++jetIter ) {
      FillHist(fSpectraEt,  AliHLTJETAnalysisBase::kSpectraRecAll, fJetsRec->GetJet(jetIter)->Pt());
      FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraRecAll, fJetsRec->GetJet(jetIter)->Eta());
      FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraRecAll, fJetsRec->GetJet(jetIter)->Phi());
    } // for ( Int_t jetIter = 0; jetIter < fJetsRec->GetNAODJets(); ++jetIter ) {
    
    // -- Fill basic Reco leading spectras
    // -------------------------------------
    if ( fJetsRec->GetNAODJets() > 0 ) {
      FillHist(fSpectraEt,  AliHLTJETAnalysisBase::kSpectraRecLeadAll, fJetsRec->GetJet(0)->Pt());
      FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraRecLeadAll, fJetsRec->GetJet(0)->Eta());
      FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraRecLeadAll, fJetsRec->GetJet(0)->Phi());
    } 
  }

  if (fJetsCmp) {
    // -- Fill basic Compare spectras
    // --------------------------------
    for ( Int_t jetIter = 0; jetIter < fJetsCmp->GetNAODJets(); ++jetIter ) {
      FillHist(fSpectraEt,  AliHLTJETAnalysisBase::kSpectraCmpAll, fJetsCmp->GetJet(jetIter)->Pt());
      FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraCmpAll, fJetsCmp->GetJet(jetIter)->Eta());
      FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraCmpAll, fJetsCmp->GetJet(jetIter)->Phi());
    } // for ( Int_t jetIter = 0; jetIter < fJetsCmp->GetNAODJets(); ++jetIter ) {

    // -- Fill basic Reco leading spectras
    // -------------------------------------
    if ( fJetsCmp->GetNAODJets() > 0 ) {
      FillHist(fSpectraEt,  AliHLTJETAnalysisBase::kSpectraCmpLeadAll, fJetsCmp->GetJet(0)->Pt());
      FillHist(fSpectraEta, AliHLTJETAnalysisBase::kSpectraCmpLeadAll, fJetsCmp->GetJet(0)->Eta());
      FillHist(fSpectraPhi, AliHLTJETAnalysisBase::kSpectraCmpLeadAll, fJetsCmp->GetJet(0)->Phi());
    }
  }

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillUnmatchedDeltaHistograms() {
  // see header file for class documentation

  if ( !fJetsRec || !fJetsCmp )
    return;

  for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {
    AliAODJet *jetCmp = fJetsCmp->GetJet(jetIterCmp);
 
    for ( Int_t jetIterRec = 0; jetIterRec < fJetsRec->GetNAODJets(); ++jetIterRec ) {
      AliAODJet *jetRec = fJetsRec->GetJet(jetIterRec);
 
      FillHist(fDeltaEt,  AliHLTJETAnalysisBase::kDeltaAll, jetCmp->Pt()  - jetRec->Pt());
      FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaAll, jetCmp->Eta() - jetRec->Eta());
      FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaAll, jetCmp->Phi() - jetRec->Phi());
      FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaAll,
	       jetCmp->Eta() - jetRec->Eta(), 
	       jetCmp->Phi() - jetRec->Phi());
      
    } // for ( Int_t jetIterRec = 0; jetIterRec < fJetsRec->GetNAODJets(); ++jetIterRec ) {
  } // for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {

  // -- Leading Jets
  // -----------------
  if ( fJetsRec->GetNAODJets() > 0 && fJetsCmp->GetNAODJets() > 0 ) {
    AliAODJet *jetCmp = fJetsCmp->GetJet(0);
    AliAODJet *jetRec = fJetsRec->GetJet(0);

    FillHist(fDeltaEt,  AliHLTJETAnalysisBase::kDeltaLead, jetCmp->Pt()  - jetRec->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaLead, jetCmp->Eta() - jetRec->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaLead, jetCmp->Phi() - jetRec->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaLead,
	     jetCmp->Eta() - jetRec->Eta(), 
	     jetCmp->Phi() - jetRec->Phi());
  }

  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedDeltaHistograms() {
  // see header file for class documentation

  for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {

    // -- Not matched jet
    if ( (*fMatchedJetsCmp)[jetIterCmp] == -1 )
      continue;
    
    AliAODJet *jetCmp = fJetsCmp->GetJet(jetIterCmp);
    AliAODJet *jetRec = fJetsRec->GetJet((*fMatchedJetsCmp)[jetIterCmp]);

    FillHist(fDeltaEt,  AliHLTJETAnalysisBase::kDeltaMatchedAll, jetCmp->Pt()  - jetRec->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaMatchedAll, jetCmp->Eta() - jetRec->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedAll, jetCmp->Phi() - jetRec->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedAll,
	     jetCmp->Eta() - jetRec->Eta(), 
	     jetCmp->Phi() - jetRec->Phi());
    
  } // for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {
  
  // -- Leading Jets
  // -----------------
  if ((*fMatchedJetsCmp)[0] != -1 ) {

    AliAODJet *jetCmp = fJetsCmp->GetJet(0);
    AliAODJet *jetRec = fJetsRec->GetJet((*fMatchedJetsCmp)[0]);
    
    FillHist(fDeltaEt,  AliHLTJETAnalysisBase::kDeltaMatchedLead, jetCmp->Pt()  - jetRec->Pt());
    FillHist(fDeltaEta, AliHLTJETAnalysisBase::kDeltaMatchedLead, jetCmp->Eta() - jetRec->Eta());
    FillHist(fDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedLead, jetCmp->Phi() - jetRec->Phi());
    FillHist(fDeltaEtaDeltaPhi, AliHLTJETAnalysisBase::kDeltaMatchedLead,
 	     jetCmp->Eta() - jetRec->Eta(), 
	     jetCmp->Phi() - jetRec->Phi());
  }
  
  return;
}
 
//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedSpectraHistograms() {
  // see header file for class documentation

  Int_t idx = 0;

  // -- Fill matched compare spectras
  // --------------------------------
  for ( Int_t jetIter = 0; jetIter < fJetsCmp->GetNAODJets(); ++jetIter ) {
    if ( (*fMatchedJetsCmp)[jetIter] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraCmpUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraCmpMatched;

    FillHist(fSpectraEt,  idx, fJetsCmp->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, idx, fJetsCmp->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, idx, fJetsCmp->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJetsCmp->GetNAODJets(); ++jetIter ) {

  // -- Fill matched compare leading spectras
  // ----------------------------------------
  if ( fJetsCmp->GetNAODJets() > 0 ) {
    if ( (*fMatchedJetsCmp)[0] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraCmpLeadUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraCmpLeadMatched; 

    FillHist(fSpectraEt,  idx, fJetsCmp->GetJet(0)->Pt());
    FillHist(fSpectraEta, idx, fJetsCmp->GetJet(0)->Eta());
    FillHist(fSpectraPhi, idx, fJetsCmp->GetJet(0)->Phi());
  }

  // -- Fill matched Reco spectras
  // -----------------------------
  for ( Int_t jetIter = 0; jetIter < fJetsRec->GetNAODJets(); ++jetIter ) {
    if ( (*fMatchedJetsRec)[jetIter] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraRecUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraRecMatched;
    
    FillHist(fSpectraEt,  idx, fJetsRec->GetJet(jetIter)->Pt());
    FillHist(fSpectraEta, idx, fJetsRec->GetJet(jetIter)->Eta());
    FillHist(fSpectraPhi, idx, fJetsRec->GetJet(jetIter)->Phi());
  } // for ( Int_t jetIter = 0; jetIter < fJetsRec->GetNAODJets(); ++jetIter ) {

  // -- Fill matched Reco leading spectras
  // -------------------------------------
  if ( fJetsRec->GetNAODJets() > 0 ) {
    if ( (*fMatchedJetsRec)[0] == -1 )
      idx = AliHLTJETAnalysisBase::kSpectraRecLeadUnmatched;
    else
      idx = AliHLTJETAnalysisBase::kSpectraRecLeadMatched;

    FillHist(fSpectraEt,  idx, fJetsRec->GetJet(0)->Pt());
    FillHist(fSpectraEta, idx, fJetsRec->GetJet(0)->Eta());
    FillHist(fSpectraPhi, idx, fJetsRec->GetJet(0)->Phi());
  } 
  
  return;
}

//##################################################################################
void AliHLTJETAnalysisJets::FillMatchedHistograms() {
  // see header file for class documentation

  for ( Int_t jetIterCmp = 0; jetIterCmp < fJetsCmp->GetNAODJets(); ++jetIterCmp ) {

    // -- Not matched jet
    if ( (*fMatchedJetsCmp)[jetIterCmp] == -1 )
      continue;
    
    AliAODJet *jetCmp = fJetsCmp->GetJet(jetIterCmp);
    AliAODJet *jetRec = fJetsRec->GetJet((*fMatchedJetsCmp)[jetIterCmp]);

    // -- Correlations    
    FillHist(fCorrelationsJetEt, AliHLTJETAnalysisBase::kPlotAll, jetCmp->Pt(), jetRec->Pt());

    // -- Resolutions
    FillHist(fResolutionsJetEt, AliHLTJETAnalysisBase::kPlotAll, jetCmp->Pt(),
	     (jetCmp->Pt()-jetRec->Pt())/(jetCmp->Pt()+jetRec->Pt()));

  } // for ( Int_t jetIter = 0; jetIter < fJetsCmp->GetNAODJets(); ++jetIter ) {
  
  // -- Leading Jets
  // -----------------
  if ((*fMatchedJetsCmp)[0] != -1 ) {

    AliAODJet *jetCmp = fJetsCmp->GetJet(0);
    AliAODJet *jetRec = fJetsRec->GetJet((*fMatchedJetsCmp)[0]);
    
    // -- Correlations    
    FillHist(fCorrelationsJetEt, AliHLTJETAnalysisBase::kPlotLead, jetCmp->Pt(), jetRec->Pt());
    
    // -- Resolutions
    FillHist(fResolutionsJetEt, AliHLTJETAnalysisBase::kPlotLead, jetCmp->Pt(),
	     (jetCmp->Pt()-jetRec->Pt())/(jetCmp->Pt()+jetRec->Pt()));
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
