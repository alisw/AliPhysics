//-*- Mode: C++ -*-
// $Id: AliHLTJETAnalysisMerge.cxx  $
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

/** @file   AliHLTJETAnalysisMerge.cxx
    @author Jochen Thaeder
    @date   
    @brief  Container merging analysis objects
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TCanvas.h"
#include "TH2F.h"

#include "AliHLTJETJets.h"

#include "AliHLTJETAnalysisMerge.h"
#include "AliHLTJETAnalysisJets.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAnalysisMerge)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTJETAnalysisMerge::AliHLTJETAnalysisMerge() :
  fCanvasArray(NULL),
  fAnalysisJetsArray(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTJETAnalysisMerge::~AliHLTJETAnalysisMerge() {
  // see header file for class documentation

  if ( fCanvasArray ) {
    fCanvasArray->Clear();
    delete fCanvasArray;
  }
  fCanvasArray = NULL;

  if ( fAnalysisJetsArray ) {
    fAnalysisJetsArray->Clear();
    delete fAnalysisJetsArray;
  }
  fAnalysisJetsArray = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                   Initialize / Reset
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETAnalysisMerge::Initialize() {
  // see header file for class documentation  

  Int_t iResult = 0;

  fCanvasArray = new TObjArray();
  fCanvasArray->SetOwner();

  fAnalysisJetsArray = new TObjArray();

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Setter - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisMerge::AddJets( AliHLTJETAnalysisJets* jets ) {
  // see header file for class documentation  

  fAnalysisJetsArray->Add( jets );

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Output - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisMerge::CreateCanvas() {
  // see header file for class documentation  

  CreateCanvasSpectra();

  CreateCanvasDelta();
  
  CreateCanvasMatched();

  return;
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
 *                             Output - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETAnalysisMerge::CreateCanvasSpectra() {
  // see header file for class documentation  

  TString canvasName = "JetSpectra";
  
  AliHLTJETAnalysisJets* jets = 
    reinterpret_cast<AliHLTJETAnalysisJets*>((*fAnalysisJetsArray)[fAnalysisJetsArray->GetLast()]);     
  
  // -- Spectra E_t
  // ----------------
  TCanvas* canvas = AddCanvas( Form("E_{t} ")+canvasName, 3, 3 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 5, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 6, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 7, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraPythiaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 8, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraPythiaMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 9, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEt, 
					      AliHLTJETAnalysisBase::kSpectraPythiaUnmatched),kTRUE, kFALSE);
  
  // -- Spectra eta
  // ----------------
  canvas = AddCanvas( Form("#eta ")+canvasName, 3, 3 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 5, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 6, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 7, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraPythiaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 8, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraPythiaMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 9, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraEta, 
					      AliHLTJETAnalysisBase::kSpectraPythiaUnmatched),kTRUE, kFALSE);

  // -- Spectra phi
  // ----------------
  canvas = AddCanvas( Form("#phi ")+canvasName, 3, 3);
  
  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 5, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 6, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraRecoLeadUnmatched),kTRUE, kFALSE);

  DrawHistogram(canvas, 7, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraPythiaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 8, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraPythiaMatched),kTRUE, kFALSE);
  DrawHistogram(canvas, 9, jets->GetHistogram(AliHLTJETAnalysisBase::kHistSpectraPhi, 
					      AliHLTJETAnalysisBase::kSpectraPythiaUnmatched),kTRUE, kFALSE);

  return;
}

//##################################################################################
void AliHLTJETAnalysisMerge::CreateCanvasDelta() {
  // see header file for class documentation  

  TString canvasName = "JetDelta";
  
  AliHLTJETAnalysisJets* jets = 
    reinterpret_cast<AliHLTJETAnalysisJets*>((*fAnalysisJetsArray)[fAnalysisJetsArray->GetLast()]);     
  
  // -- Delta E_t
  // --------------
  TCanvas* canvas = AddCanvas( Form("#DeltaE_{t} ")+canvasName, 2, 2 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEt, 
					      AliHLTJETAnalysisBase::kDeltaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEt, 
					      AliHLTJETAnalysisBase::kDeltaMatchedAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEt, 
					      AliHLTJETAnalysisBase::kDeltaLead),kTRUE, kFALSE);
  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEt, 
					      AliHLTJETAnalysisBase::kDeltaMatchedLead),kTRUE, kFALSE);

  // -- Delta Eta
  // --------------
  canvas = AddCanvas( Form("#Delta#eta} ")+canvasName, 2, 2 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEta, 
					      AliHLTJETAnalysisBase::kDeltaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEta, 
					      AliHLTJETAnalysisBase::kDeltaMatchedAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEta, 
					      AliHLTJETAnalysisBase::kDeltaLead),kTRUE, kFALSE);
  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaEta, 
					      AliHLTJETAnalysisBase::kDeltaMatchedLead),kTRUE, kFALSE);

  // -- Delta Phi
  // --------------
  canvas = AddCanvas( Form("#Delta#phi ")+canvasName, 2, 2 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaPhi, 
					      AliHLTJETAnalysisBase::kDeltaAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaPhi, 
					      AliHLTJETAnalysisBase::kDeltaMatchedAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaPhi, 
					      AliHLTJETAnalysisBase::kDeltaLead),kTRUE, kFALSE);
  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistDeltaPhi, 
					      AliHLTJETAnalysisBase::kDeltaMatchedLead),kTRUE, kFALSE);

  // -- Delta Eta Delta Phi
  // ------------------------
  canvas = AddCanvas( Form("#Delta#eta#Delta#phi ")+canvasName, 2, 2 );

  DrawHistogram(canvas, 1, jets->GetHistogram( AliHLTJETAnalysisBase::kHistDeltaEtaDeltaPhi, 
					       AliHLTJETAnalysisBase::kDeltaAll ),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram( AliHLTJETAnalysisBase::kHistDeltaEtaDeltaPhi, 
					       AliHLTJETAnalysisBase::kDeltaMatchedAll ),kTRUE, kFALSE);
  DrawHistogram(canvas, 3, jets->GetHistogram( AliHLTJETAnalysisBase::kHistDeltaEtaDeltaPhi, 
					       AliHLTJETAnalysisBase::kDeltaLead ),kTRUE, kFALSE);
  DrawHistogram(canvas, 4, jets->GetHistogram( AliHLTJETAnalysisBase::kHistDeltaEtaDeltaPhi, 
					       AliHLTJETAnalysisBase::kDeltaMatchedLead ),kTRUE, kFALSE);

  return;
}

//##################################################################################
void AliHLTJETAnalysisMerge::CreateCanvasMatched() {
  // see header file for class documentation  

  TString canvasName = "JetMatched";
  
  AliHLTJETAnalysisJets* jets = 
    reinterpret_cast<AliHLTJETAnalysisJets*>((*fAnalysisJetsArray)[fAnalysisJetsArray->GetLast()]);     
  
  // -- Resolutions
  // ----------------
  TCanvas* canvas = AddCanvas( Form("Correlations + Resolutions")+canvasName, 2, 3 );

  DrawHistogram(canvas, 1, jets->GetHistogram(AliHLTJETAnalysisBase::kHistCorrelationsJetEt, 
					      AliHLTJETAnalysisBase::kPlotAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 2, jets->GetHistogram(AliHLTJETAnalysisBase::kHistCorrelationsJetEt, 
					      AliHLTJETAnalysisBase::kPlotLead),kTRUE, kFALSE);

  DrawHistogram(canvas, 3, jets->GetHistogram(AliHLTJETAnalysisBase::kHistResolutionsJetEt,
					      AliHLTJETAnalysisBase::kPlotAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 4, jets->GetHistogram(AliHLTJETAnalysisBase::kHistResolutionsJetEt, 
					      AliHLTJETAnalysisBase::kPlotLead),kTRUE, kFALSE);

  DrawHistogram(canvas, 5, jets->GetHistogram(AliHLTJETAnalysisBase::kHistResolutionsDiJetEt, 
					      AliHLTJETAnalysisBase::kPlotAll),kTRUE, kFALSE);
  DrawHistogram(canvas, 6, jets->GetHistogram(AliHLTJETAnalysisBase::kHistResolutionsDiJetEt, 
					      AliHLTJETAnalysisBase::kPlotLead),kTRUE, kFALSE);

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                               Helper - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
TCanvas* AliHLTJETAnalysisMerge::AddCanvas( TString name, Int_t divideX=1, Int_t divideY=1 ) {
  // see header file for class documentation

  fCanvasArray->Add (new TCanvas(name, name, 10, 10, 1400, 800));
  reinterpret_cast<TCanvas*>((*fCanvasArray)[fCanvasArray->GetLast()])->Divide(divideX,divideY);
  
  return reinterpret_cast<TCanvas*>((*fCanvasArray)[fCanvasArray->GetLast()]);
}

//##################################################################################
void AliHLTJETAnalysisMerge::DrawHistogram( TCanvas* canvas, Int_t idx, TH1* hist, 
					    Bool_t bScale=kFALSE, Bool_t bLogY=kFALSE) {
  // see header file for class documentation

  if ( hist == NULL )
    return;

  TVirtualPad* pad = canvas->cd(idx);
  
  if ( bScale ) 
    hist->Scale( 1./hist->GetBinWidth(0) );
  
  if ( bLogY && hist->GetEntries() != 0 )
    pad->SetLogy();
  
  pad->SetGridy();
  pad->SetGridx();
  
  if ( !strcmp(hist->ClassName(),"TH2F") )
    hist->Draw("COLZ");
  else 
    hist->Draw("");
  
  return;
}
