/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskFlowSingleMu.cxx 55545 2012-04-04 07:16:39Z pcrochet $ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskFlowSingleMu
/// Analysis task for flow of single muons in the spectrometer.
/// The output is a list of histograms and CF containers.
/// The macro class can run on AODs or ESDs.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#define AliAnalysisTaskFlowSingleMu_cxx

#include "AliAnalysisTaskFlowSingleMu.h"

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TStyle.h"
//#include "TMCProcess.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TLatex.h"

// STEER includes
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliVHeader.h"
#include "AliAODMCHeader.h"
#include "AliStack.h"
#include "AliEventplane.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"

// CORRFW includes
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

// PWG includes
#include "AliVAnalysisMuon.h"
#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"
#include "AliMuonTrackCuts.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskFlowSingleMu) // Class implementation in ROOT context
/// \endcond

using std::ifstream;

//________________________________________________________________________
AliAnalysisTaskFlowSingleMu::AliAnalysisTaskFlowSingleMu() :
  AliVAnalysisMuon(),
  fEPKeys(0x0),
  fRandom(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskFlowSingleMu::AliAnalysisTaskFlowSingleMu(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fEPKeys(0x0),
  fRandom(0x0)
{
  //
  /// Constructor.
  //
  TString epMethods = "V0A V0C TPC V0Cr1 V0Cr4";
  fEPKeys = epMethods.Tokenize(" ");
}


//________________________________________________________________________
AliAnalysisTaskFlowSingleMu::~AliAnalysisTaskFlowSingleMu()
{
  //
  /// Destructor
  //
  delete fEPKeys;
  delete fRandom;
}

//___________________________________________________________________________
void AliAnalysisTaskFlowSingleMu::MyUserCreateOutputObjects()
{
  //
  /// Create outputs
  //

  if ( ! fRandom ) fRandom = new TRandom3();
  // Set seed here, to avoid having the same seed for all jobs in grid
  fRandom->SetSeed(0);

  Int_t nPtBins = 160;
  Double_t ptMin = 0., ptMax = 80.;
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");
  
  Int_t nEtaBins = 25;
  Double_t etaMin = -4.5, etaMax = -2.;
  TString etaName("Eta"), etaTitle("#eta"), etaUnits("");
  
  Int_t nPhiBins = 36;
  Double_t phiMin = 0.; Double_t phiMax = 2.*TMath::Pi();
  TString phiName("Phi"), phiTitle("#phi"), phiUnits("rad");

  Int_t nDeltaPhiBins = 18;
  Double_t deltaPhiMin = 0.; Double_t deltaPhiMax = TMath::Pi();
  TString deltaPhiName("DeltaPhi"), deltaPhiTitle("#phi_{#mu} - #psi_{plane}"), deltaPhiUnits("rad");
    
  Int_t nChargeBins = 2;
  Double_t chargeMin = -2, chargeMax = 2.;
  TString chargeName("Charge"), chargeTitle("charge"), chargeUnits("e");
  
  Int_t nMotherTypeBins = kNtrackSources;
  Double_t motherTypeMin = -0.5, motherTypeMax = (Double_t)kNtrackSources - 0.5;
  TString motherType("MotherType"), motherTypeTitle("motherType"), motherTypeUnits("");

  Int_t nbins[kNvars] = {nPtBins, nEtaBins, nPhiBins, nDeltaPhiBins, nChargeBins, nMotherTypeBins};
  Double_t xmin[kNvars] = {ptMin, etaMin, phiMin, deltaPhiMin, chargeMin, motherTypeMin};
  Double_t xmax[kNvars] = {ptMax, etaMax, phiMax, deltaPhiMax, chargeMax, motherTypeMax};
  TString axisTitle[kNvars] = {ptTitle, etaTitle, phiTitle, deltaPhiTitle, chargeTitle, motherTypeTitle};
  TString axisUnits[kNvars] = {ptUnits, etaUnits, phiUnits, deltaPhiUnits, chargeUnits, motherTypeUnits};

  Int_t iobj=0;
  TH1* histo = 0x0;
  TString currTitle = "";
  for ( Int_t iep=0; iep<fEPKeys->GetEntries(); iep++ ) {
    TString currEP = fEPKeys->At(iep)->GetName();
    AliCFContainer* cfContainer = new AliCFContainer(Form("FlowSingleMuContainer%s", currEP.Data()),"Container for tracks",kNsteps,kNvars,nbins);
  
    TString currName = "";
    for ( Int_t idim = 0; idim<kNvars; idim++){
      currName = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
      currName.ReplaceAll("()","");

      cfContainer->SetVarTitle(idim, currName.Data());
      cfContainer->SetBinLimits(idim, xmin[idim], xmax[idim]);
    }
    
    TString stepTitle[kNsteps] = {"reconstructed", "generated"};

    TAxis* currAxis = 0x0;
    for (Int_t istep=0; istep<kNsteps; istep++){
      cfContainer->SetStepTitle(istep, stepTitle[istep].Data());
      AliCFGridSparse* gridSparse = cfContainer->GetGrid(istep);

      currAxis = gridSparse->GetAxis(kHvarMotherType);
      for ( Int_t ibin=0; ibin<fSrcKeys->GetEntries(); ibin++ ) {
        currAxis->SetBinLabel(ibin+1, fSrcKeys->At(ibin)->GetName());
      }
    }

    AddObjectToCollection(cfContainer, iobj++);

    Bool_t isTPC = currEP.Contains("TPC");
    Double_t minPsiPlane = ( isTPC ) ? 0. : -TMath::Pi()/2.;
    Double_t maxPsiPlane = ( isTPC ) ? TMath::Pi() : TMath::Pi()/2.;
    
    histo = new TH1D(Form("hEventPlane%s", currEP.Data()), Form("Event plane distribution: %s", currEP.Data()), 18, minPsiPlane, maxPsiPlane);
    histo->GetXaxis()->SetTitle(Form("#psi_{%s} (rad)", currEP.Data()));
    AddObjectToCollection(histo, iobj++);

    for ( Int_t jep=iep+1; jep<fEPKeys->GetEntries(); jep++ ) {
      TString auxEP = fEPKeys->At(jep)->GetName();
      currTitle = Form("cos(#psi_{%s} - #psi_{%s})",currEP.Data(),auxEP.Data());
      histo = new TH1D(Form("hResoSub3%s%s",currEP.Data(),auxEP.Data()), currTitle.Data(), 100, -1.,1.);
      histo->GetXaxis()->SetTitle(currTitle.Data());
      AddObjectToCollection(histo, iobj++);
    }

    currTitle = Form("cos(#psi_{1} - #psi_{2})");
    histo = new TH1D(Form("hResoSub2%s",currEP.Data()), Form("%s: %s", currEP.Data(), currTitle.Data()), 100, -1.,1.);
    histo->GetXaxis()->SetTitle(currTitle.Data());

    AddObjectToCollection(histo, iobj++);
  } // loop on EPs

  fMuonTrackCuts->Print("mask");
}

//________________________________________________________________________
void AliAnalysisTaskFlowSingleMu::ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality)
{
  //
  /// Fill output objects
  //

  //
  // Base event cuts
  //
  if ( GetVertexSPD()->GetNContributors() < fMinNvtxContirbutors ) return;
  //Bool_t isPileupFromSPD = ( fAODEvent && ! fAODEvent->GetTracklets() ) ? InputEvent()->IsPileupFromSPD(3, 0.8, 3., 2., 5.) : InputEvent()->IsPileupFromSPDInMultBins(); // Avoid break when reading Muon AODs (tracklet info is not present and IsPileupFromSPDInMultBins crashes
  //if ( isPileupFromSPD ) return;
  if ( TMath::Abs(GetVertexSPD()->GetZ()) > 10. ) return;

  TArrayD psiPlane(fEPKeys->GetEntriesFast());
  psiPlane.Reset(-999.);
  const Double_t kLowestVal = -10.;
  Int_t harmonic = 2;

  //
  // Fill EP info
  //
  AliEventplane* eventPlane = InputEvent()->GetEventplane();
  psiPlane[0] = eventPlane->GetEventplane("V0A",InputEvent(),harmonic);
  psiPlane[1] = eventPlane->GetEventplane("V0C",InputEvent(),harmonic);
  psiPlane[2] = eventPlane->GetEventplane("Q");
  if ( psiPlane[2] < 0.) psiPlane[2] = -999.;
  Double_t qx = 0., qy = 0.;
  psiPlane[3] = eventPlane->CalculateVZEROEventPlane(InputEvent(),0,harmonic,qx,qy); // V0C ring 1
  psiPlane[4] = eventPlane->CalculateVZEROEventPlane(InputEvent(),3,harmonic,qx,qy); // V0C ring 3
  //psiPlane[3] = fRandom->Uniform(-TMath::Pi()/2., TMath::Pi()/2.);

  Bool_t noValidEP = kTRUE;
  TArrayD psiPlaneCalc(fEPKeys->GetEntriesFast());
  for ( Int_t iep=0; iep<fEPKeys->GetEntriesFast(); iep++ ) {
    psiPlaneCalc[iep] = psiPlane[iep];
    if ( psiPlane[iep] < kLowestVal ) continue;
    if ( psiPlane[iep] < 0. ) psiPlaneCalc[iep] += TMath::Pi();
    noValidEP = kFALSE;
  }
  if ( noValidEP ) return;

  for ( Int_t iep=0; iep<fEPKeys->GetEntriesFast(); iep++ ) {
    if ( psiPlane[iep] < kLowestVal ) continue; 
    for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
      TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
      TString currEP = fEPKeys->At(iep)->GetName();
      ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, Form("hEventPlane%s",currEP.Data())))->Fill(psiPlane[iep]);
      for ( Int_t jep=iep+1; jep<fEPKeys->GetEntriesFast(); jep++ ) {
        if ( psiPlane[jep] < kLowestVal ) continue; 
        ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, Form("hResoSub3%s%s",currEP.Data(), fEPKeys->At(jep)->GetName())))->Fill(TMath::Cos(2.*(psiPlaneCalc[iep]-psiPlaneCalc[jep])));
      } // loop on auxialiary EP
      if ( iep == 2 ) {
        TVector2* qsub1 = eventPlane->GetQsub1();
        TVector2* qsub2 = eventPlane->GetQsub2();
        if ( qsub1 && qsub2 )
          ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, Form("hResoSub2%s",currEP.Data())))->Fill(TMath::Cos(qsub1->Phi()/2.-qsub2->Phi()/2.));
      }
    } // loop on trigger
  } // loop on EP

  //
  // Fill track info
  //
  Double_t containerInput[kNvars];
  AliVParticle* track = 0x0;
  for ( Int_t istep = 0; istep<2; ++istep ) {
    Int_t nTracks = ( istep == kStepReconstructed ) ? GetNTracks() : GetNMCTracks();
    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      track = ( istep == kStepReconstructed ) ? GetTrack(itrack) : GetMCTrack(itrack);
      
      Bool_t isSelected = ( istep == kStepReconstructed ) ? fMuonTrackCuts->IsSelected(track) : ( TMath::Abs(track->PdgCode()) == 13 );
      if ( ! isSelected ) continue;
      
      Int_t trackSrc = ( istep == kStepReconstructed ) ? GetParticleType(track) : RecoTrackMother(track);
 

      containerInput[kHvarPt]         = track->Pt();
      containerInput[kHvarEta]        = track->Eta();
      containerInput[kHvarPhi]        = track->Phi();
      containerInput[kHvarCharge]     = track->Charge()/3.;
      containerInput[kHvarMotherType] = (Double_t)trackSrc;

      for ( Int_t iep=0; iep<fEPKeys->GetEntriesFast(); iep++ ) {
        if ( psiPlane[iep] < kLowestVal ) continue; 
        TString currEP = fEPKeys->At(iep)->GetName();
        Double_t deltaPhi = containerInput[kHvarPhi] - psiPlane[iep];

        // The difference between phi - psiPlane must be between (0, pi)
        // Since the reaction plain is symmetric in pi,
        // we can do in such a way that:
        // -pi < delta_phi + N*pi < pi
        // We choose N accrodingly
        Double_t nPi = TMath::Floor( 1. - deltaPhi / TMath::Pi() );
        deltaPhi += nPi * TMath::Pi();
        containerInput[kHvarDeltaPhi] = deltaPhi;

        for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
          TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
          if ( istep == kStepReconstructed && ! TrackPtCutMatchTrigClass(track, trigClassName) ) continue;
          ((AliCFContainer*)GetMergeableObject(physSel, trigClassName, centrality, Form("FlowSingleMuContainer%s",currEP.Data())))->Fill(containerInput,istep);
        } // loop on selected trigger classes
      } // loop on event plane
    } // loop on tracks
  } // loop on container steps
}


//________________________________________________________________________
TArrayD AliAnalysisTaskFlowSingleMu::GetCentralityRange(TString sRange)
{
  //
  /// Get the range from string
  //
  TArrayD centralityRange(2);
  centralityRange.Reset(-999.);
  TObjArray* array = sRange.Tokenize("_");
  if ( array->GetEntries() == 2 ) {
    for ( Int_t ientry=0; ientry<2; ientry++ ) {
      centralityRange[ientry] = ((TObjString*)array->At(ientry))->GetString().Atof();
    }
  }
  delete array;
  return centralityRange;
}

//________________________________________________________________________
void AliAnalysisTaskFlowSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histograms at the end.
  //

  AliVAnalysisMuon::Terminate("");

  if ( ! fMergeableCollection ) return;
  
  TString physSel = fTerminateOptions->At(0)->GetName();
  TString trigClassName = fTerminateOptions->At(1)->GetName();
  TString centralityRange = fTerminateOptions->At(2)->GetName();
  TString furtherOpt = fTerminateOptions->At(3)->GetName();

  TObjArray* centralityClasses = centralityRange.Tokenize(",");
  TArrayD centralityArray(centralityClasses->GetEntries()+1);
  Int_t ib = 0;
  for ( Int_t icent=0; icent<centralityClasses->GetEntries(); icent++ ) {
    TArrayD range = GetCentralityRange(centralityClasses->At(icent)->GetName());
    if ( icent == 0 ) centralityArray.AddAt(range[0],ib++);
    centralityArray.AddAt(range[1],ib++);
  }
  
  TString selectEP = "";
  TObjArray resoEParray;
  resoEParray.SetOwner();
  
  TString outFilename = "";
  TObjArray* optArr = furtherOpt.Tokenize(" ");
  TString currName = "";
  TString resoFilename = "";
  for ( Int_t iopt=0; iopt<optArr->GetEntries(); iopt++ ) {
    currName = optArr->At(iopt)->GetName();
    if ( currName.Contains(".root") ) outFilename = currName;
    else if ( currName.Contains(".txt") ) resoFilename = currName;
    else {
      for ( Int_t iep=0; iep<fEPKeys->GetEntries(); iep++ ) {
        if ( ! currName.CompareTo(fEPKeys->At(iep)->GetName()) ) resoEParray.Add(new TObjString(currName.Data()));
      }
    }
  }
  delete optArr;

  if ( resoEParray.At(0) ) selectEP = resoEParray.At(0)->GetName();

  furtherOpt.ToUpper();

  Bool_t v2only = furtherOpt.Contains("V2ONLY");

  TAxis* customBins = 0x0;
  
  //Double_t myBins[] = {0., 0.5, 1., 1.5, 2., 3., 4., 5., 6., 8., 10.};
  Double_t myBins[] = {4., 5., 6., 8., 10.}; // REMEMBER TO CHOOSE
  //Double_t myBins[] = {2., 3., 4., 5., 6., 8., 10.};
  //Double_t myBins[] = {0., 0.5, 1., 1.5, 2., 3., 4., 10.};
  //Double_t myBins[] = {2., 3., 4., 5., 6., 8., 10., 80.};
  Int_t nMyBins = sizeof(myBins)/sizeof(myBins[0])-1;
  customBins = new TAxis(nMyBins, myBins); // REMEMBER TO CHOOSE

  TCanvas *can = 0x0, *showCan = 0x0;
  Int_t xshift = 100;
  Int_t yshift = 20;
  Int_t igroup1 = 0;
  Int_t igroup2 = 0;

  //
  /// Read plane resolution file (if present)
  //
  TObjArray resoCentrality;
  resoCentrality.SetOwner();
  TArrayD resolution[3];
  for ( Int_t ires=0; ires<3; ires++ ) {
    resolution[ires].Set(50);
  }

  Bool_t externalReso = ( ! resoFilename.IsNull() && ! gSystem->AccessPathName(resoFilename.Data()) );
  Bool_t internalReso = kTRUE;
  if ( externalReso ) {
    // If plane resolution file exists, fill resolution
    TString currLine = "";
    ifstream inFile(resoFilename.Data());
    if (inFile.is_open()) {
      ib = 0;
      while (! inFile.eof() ) {
        currLine.ReadLine(inFile,kFALSE);
        if ( currLine.BeginsWith("#") || currLine.Length() < 3 ) continue;
        TObjArray* array = currLine.Tokenize(" ");
        resoCentrality.Add(new TObjString(array->At(0)->GetName()));
        for ( Int_t ires=0; ires<3; ires++ ) {
          resolution[ires].AddAt(((TObjString*)array->At(ires+1))->GetString().Atof(),ib);
        }
        delete array;
        ib++;
      } // ! EOF
    } // file is open
    for ( Int_t ires=0; ires<3; ires++ ) {
      resolution[ires].Set(ib);
    }
  } // file exists
  else {
    // If no plane resolution file exist,
    // use the centralities in terminateOption
    // and set the resolution to 1 and errors to 0.

    for ( Int_t ires=0; ires<3; ires++ ) {
      resolution[ires].Set(centralityClasses->GetEntries());
      Double_t resetVal = ( ires == 0 ) ? 1. : 0.;
      resolution[ires].Reset(resetVal);
    }
    TArrayD currCos(3);
    for ( Int_t icent=0; icent<centralityClasses->GetEntries(); icent++ ) {
      resoCentrality.Add(new TObjString(centralityClasses->At(icent)->GetName()));
      Int_t icos = 0;
      Double_t sumRelErrSquare = 0.;
      for ( Int_t isub=0; isub<resoEParray.GetEntries(); isub++ ) {
        for ( Int_t jsub=0; jsub<resoEParray.GetEntries(); jsub++ ) {
          TH1* histo = static_cast<TH1*>(GetSum(physSel,trigClassName,resoCentrality.At(icent)->GetName(),Form("hResoSub3%s%s", resoEParray.At(isub)->GetName(),resoEParray.At(jsub)->GetName())));
          if ( ! histo || histo->Integral() == 0 ) continue;
          currCos[icos] = histo->GetMean();
          if ( currCos[icos] <= 0. ) continue;
          Double_t relErr = histo->GetMeanError() / currCos[icos];
          sumRelErrSquare += relErr * relErr;
          icos++;
        } // loop on sub events
      } // loop on sub events
      if ( icos != 3  ) {
        printf("Warning: resolution cannot be estimated for %s\n", resoCentrality.At(icent)->GetName());
        internalReso = kFALSE;
        continue;
      }
      resolution[0][icent] = TMath::Sqrt((currCos[0]*currCos[1]/currCos[2]));
      resolution[1][icent] = resolution[0][icent]/2.*TMath::Sqrt(sumRelErrSquare);
      printf("Resolution %s  %g +- %g\n", centralityClasses->At(icent)->GetName(), resolution[0][icent], resolution[1][icent]);
    }
  }

  Bool_t hasResolution = ( externalReso || internalReso );


  TArrayD resoCentrArray(resoCentrality.GetEntries()+1);
  ib = 0;
  for ( Int_t icent=0; icent<resoCentrality.GetEntries(); icent++ ) {
    TArrayD range = GetCentralityRange(resoCentrality.At(icent)->GetName());
    if ( icent == 0 ) resoCentrArray.AddAt(range[0],ib++);
    resoCentrArray.AddAt(range[1], ib++);
  }

  if ( hasResolution ) {
    currName = externalReso ? Form("externalReso") : Form("reso_%s_%s_%s", resoEParray.At(0)->GetName(), resoEParray.At(1)->GetName(), resoEParray.At(2)->GetName());
    can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
    TGraphErrors* resoGraph = new TGraphErrors(centralityClasses->GetEntries());
    for ( Int_t icent=0; icent<resoCentrality.GetEntries(); icent++ ) {

      resoGraph->SetPoint(icent, 0.5*(resoCentrArray[icent+1] + resoCentrArray[icent]), resolution[0][icent]);
      resoGraph->SetPointError(icent, 0.5*(resoCentrArray[icent+1] - resoCentrArray[icent]), resolution[1][icent]);
    }
    resoGraph->SetTitle(Form("Resolution %s", selectEP.Data()));
    resoGraph->GetXaxis()->SetTitle("Centrality");
    resoGraph->GetYaxis()->SetTitle("Resolution");
    resoGraph->Draw("apz");
  }

  const Int_t kNresoCentr = resoCentrality.GetEntries();

  //
  // Create gridSparse array
  //
  AliCFGridSparse* gridSparseArray[kNsteps*kNresoCentr];
  TString stepName[2];
  for ( Int_t icent=0; icent<kNresoCentr; icent++ ) {
    for ( Int_t istep=0; istep<kNsteps; istep++ ) {
      gridSparseArray[istep*kNresoCentr+icent] = 0x0;
    }

    AliCFContainer* cfContainer = static_cast<AliCFContainer*> ( GetSum(physSel,trigClassName,resoCentrality.At(icent)->GetName(),Form("FlowSingleMuContainer%s", selectEP.Data())) );
    if ( ! cfContainer ) continue;
    for ( Int_t istep=0; istep<kNsteps; istep++ ) {
      stepName[istep] = cfContainer->GetStepTitle(istep);
      gridSparseArray[istep*kNresoCentr+icent] = cfContainer->GetGrid(istep);
    } // loop on step
  } // loop on centrality

  Int_t srcColors[kNtrackSources] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange};
  
  Bool_t isMC = furtherOpt.Contains("MC");
  Int_t firstSrc = ( isMC ) ? 0 : kUnidentified;
  Int_t lastSrc  = ( isMC ) ? kNtrackSources - 1 : kUnidentified;
  if ( ! isMC ) srcColors[kUnidentified] = 1;

  TList outList;

  TString graphTypeName[3] = {"raw","correctStat","correctSyst"};
  Int_t drawOrder[3] = {0, 2, 1};
  TString drawOpt[3] = {"apz","pz","a2"};
  Int_t nTypes = ( hasResolution ) ? 3 : 1;

  TString histoName = "", histoPattern = "";
  ///////////
  // Flow  //
  ///////////
  gStyle->SetOptFit(1111);
  TString funcFormula = ( v2only ) ? "[0] * (1. + 2.*[1]*TMath::Cos(2*x))" : "[0] * (1. + 2.*([1]*TMath::Cos(x) + [2]*TMath::Cos(2*x) + [3]*TMath::Cos(3*x)))";
  //
  TF1* func = new TF1("funcFlow",funcFormula.Data(), 0., TMath::Pi());
  if  ( v2only ) func->SetParNames("scale", "v2");
  else func->SetParNames("scale","v1", "v2", "v3");
  Int_t v2par = ( v2only ) ? 1 : 2;

  for ( Int_t istep=0; istep<kNsteps; istep++ ) {
    for ( Int_t isrc = firstSrc; isrc <= lastSrc; ++isrc ) {
      TString baseName = Form("%s_%s", stepName[istep].Data(), fSrcKeys->At(isrc)->GetName());
      TGraphErrors* flowVsCentrality[3];
      for ( Int_t itype=0; itype<nTypes; itype++ ) {
        histoName = Form("v2VsCentrality_%s_%s", baseName.Data(), graphTypeName[itype].Data());
        flowVsCentrality[itype] = new TGraphErrors();
        flowVsCentrality[itype]->SetName(histoName.Data());
        flowVsCentrality[itype]->SetTitle(histoName.Data());
      }
      for ( Int_t icent=0; icent<centralityClasses->GetEntries(); icent++ ) {
        TH2* hDeltaPhi2D = 0x0;
        TString baseNameCent = Form("%s_%s", baseName.Data(), centralityClasses->At(icent)->GetName());
        TArrayD weightedReso(3);
        weightedReso.Reset(0.);
        Double_t nMuons = 0.;
        TString debugString = Form("\n%s", centralityClasses->At(icent)->GetName());
        for ( Int_t rcent=0; rcent<kNresoCentr; rcent++ ) {
          if ( resoCentrArray[rcent] < centralityArray[icent] ||
               resoCentrArray[rcent+1] > centralityArray[icent+1] ) continue;
          AliCFGridSparse* gridSparse = gridSparseArray[istep*kNresoCentr+rcent];
          if ( ! gridSparse ||  gridSparse->GetEntries() == 0. ) continue;

          SetSparseRange(gridSparse, kHvarEta, "", -3.999, -2.501);
          SetSparseRange(gridSparse, kHvarMotherType, "", isrc+1, isrc+1, "USEBIN");
          TH1* auxHisto = gridSparse->Project(kHvarDeltaPhi,kHvarPt);

          Double_t currYield = auxHisto->Integral();
          if ( currYield == 0. ) continue;
          nMuons += currYield;
          debugString += Form("\nAdding %s  yield %g  reso", resoCentrality.At(rcent)->GetName(), currYield);
          for ( Int_t ires=0; ires<3; ires++ ) {
            weightedReso[ires] += resolution[ires][rcent] * currYield;
            debugString += Form("  %g", resolution[ires][rcent]);
          }
          histoName = Form("ptVsDeltaPhi_%s", baseNameCent.Data());
          if ( ! hDeltaPhi2D ) hDeltaPhi2D = static_cast<TH2*>(auxHisto->Clone(histoName.Data()));
          else hDeltaPhi2D->Add(auxHisto);
          delete auxHisto;
        }

        if ( ! hDeltaPhi2D ) continue;

        debugString += Form("\nWeighted reso ");
        for ( Int_t ires=0; ires<3; ires++ ) {
          weightedReso[ires] = weightedReso[ires]/nMuons;
          debugString += Form("  %g", weightedReso[ires]);
        }
        AliDebug(1,debugString.Data());

        TAxis* ptAxis = hDeltaPhi2D->GetYaxis();

        Int_t ipad = 0;
        TGraphErrors* flowVsPt[3];
        for ( Int_t itype=0; itype<nTypes; itype++ ) {
          histoName = Form("v2VsPt_%s_%s", baseNameCent.Data(), graphTypeName[itype].Data());
          flowVsPt[itype] = new TGraphErrors();
          flowVsPt[itype]->SetName(histoName.Data());
          flowVsPt[itype]->SetTitle(histoName.Data());
        }

        if ( ! customBins ) customBins = static_cast<TAxis*>(ptAxis->Clone());


        for ( Int_t ipt=0; ipt<=customBins->GetNbins(); ipt++ ) {
          Int_t minCustomBin = ( ipt == 0 ) ? 1 : ipt;
          Int_t maxCustomBin = ( ipt == 0 ) ? customBins->GetNbins() : ipt;
          Int_t ptMinBin = ptAxis->FindBin(customBins->GetBinLowEdge(minCustomBin)+1.e-4);
          Int_t ptMaxBin = ptAxis->FindBin(customBins->GetBinUpEdge(maxCustomBin)-1.e-4);
          Double_t ptMin = ptAxis->GetBinLowEdge(ptMinBin);
          Double_t ptMax = ptAxis->GetBinUpEdge(ptMaxBin);
          TH1* projHisto = hDeltaPhi2D->ProjectionX(Form("checkHisto_%s_ptBin%i", baseNameCent.Data(), ipt), ptMinBin, ptMaxBin);
          projHisto->SetTitle(Form("%g < p_{t} (GeV/c) < %g", ptMin, ptMax));
          if ( projHisto->Integral() < 50 ) break;
          if ( ipad%4 == 0 ) {
            currName = projHisto->GetName();
            currName.Append("_can");
            showCan = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
            showCan->Divide(2,2);
            ipad = 0;
          }
          showCan->cd(++ipad);
          if ( v2only ) func->SetParameters(projHisto->Integral(), 0.01);
          else {
            func->SetParameters(projHisto->Integral(), 0., 0.01, 0.);
            func->FixParameter(1,0.);
          }
          projHisto->Fit(func,"R");
          for ( Int_t itype=0; itype<nTypes; itype++ ) {
            TGraphErrors* currGraph = ( ipt == 0 ) ? flowVsCentrality[itype] : flowVsPt[itype];
            Double_t resoVal = ( itype == 0 ) ? 1. : weightedReso[0];
            Double_t resoErr = ( itype == 0 ) ? 0. : weightedReso[itype];
            Double_t rawVal = func->GetParameter(v2par);
            Double_t rawErr = (itype == 2 ) ? 0. : func->GetParError(v2par);
            Double_t yVal = rawVal/resoVal;
            Double_t xVal = ( ipt == 0 ) ? 0.5*(centralityArray[icent+1] + centralityArray[icent]) : 0.5*(ptMax + ptMin);
            Double_t xErr = ( ipt == 0 ) ? 0.5*(centralityArray[icent+1] - centralityArray[icent]) : 0.5*(ptMax - ptMin);
            Int_t ipoint = currGraph->GetN();
            currGraph->SetPoint(ipoint,xVal,yVal);
            if ( itype > 0 && ipt != 0 ) {
              // For the plot vs. pt the resolution error is fully correlated.
              // Hence we do not use add it bin by bin, but rather write it
              // on the title
              resoErr = 0.;
              TString graphTitle = currGraph->GetTitle();
              if ( ! graphTitle.Contains("|") ) {
                graphTitle += Form(" | Rel. error on EP resolution: %.2g (stat.) %.2g (syst.)", weightedReso[1]/weightedReso[0], weightedReso[2]/weightedReso[0]);
                currGraph->SetTitle(graphTitle.Data());
              }
            }
            Double_t rawErrRel = ( rawVal != 0. ) ? rawErr/rawVal : 0.;
            Double_t resoErrRel = ( resoVal != 0. ) ? resoErr/resoVal : 0.;
            Double_t yErr = TMath::Abs(yVal) * TMath::Sqrt(rawErrRel*rawErrRel+resoErrRel*resoErrRel);
            currGraph->SetPointError(ipoint,xErr,yErr);
          } // loop on type
        } // loop on pt bins
        for ( Int_t itype=0; itype<nTypes; itype++ ) {
          Int_t currType = drawOrder[itype];
          if ( itype < 2 ) {
            currName = flowVsPt[currType]->GetName();
            currName.Append("Can");
            can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
            igroup2++;
            can->cd();
          }
          flowVsPt[currType]->GetXaxis()->SetTitle(ptAxis->GetTitle());
          flowVsPt[currType]->GetYaxis()->SetTitle("v2");
          flowVsPt[currType]->GetYaxis()->SetRangeUser(0., 0.25);
          flowVsPt[currType]->SetFillStyle(0);
          flowVsPt[currType]->Draw(drawOpt[currType].Data());
          outList.Add(flowVsPt[currType]);
        } // loop on types
      } // loop on centrality
      for ( Int_t itype=0; itype<nTypes; itype++ ) {
        Int_t currType = drawOrder[itype];
        if ( flowVsCentrality[currType]->GetN() == 0 ) continue;
        if ( itype < 2 ) {
          currName = flowVsCentrality[currType]->GetName();
          currName.Append("Can");
          can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
          can->cd();
        }
        flowVsCentrality[currType]->GetXaxis()->SetTitle("Centrality");
        flowVsCentrality[currType]->GetYaxis()->SetTitle("v2");
        flowVsCentrality[currType]->SetFillStyle(0);
        flowVsCentrality[currType]->GetYaxis()->SetRangeUser(0., 0.25);
        flowVsCentrality[currType]->Draw(drawOpt[currType].Data());
        outList.Add(flowVsCentrality[currType]);
        igroup2++;
      } // loop on type
      igroup1++;
      igroup2 = 0;
    } // loop on track sources      
  } // loop on steps

  delete customBins;


  ///////////////////////
  // Event plane check //
  ///////////////////////
  igroup1++;
  igroup2 = 0;
  Int_t icolor=0;
  currName ="eventPlaneDistrib";
  can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
  can->SetGridy();
  igroup2++;

  Bool_t addOffsetToCentrality = kFALSE;

  TLegend* leg = new TLegend(0.7,0.7,0.92,0.92);
  leg->SetBorderSize(1);
  for ( Int_t icent=0; icent<centralityClasses->GetEntries(); icent++ ) {
    TH1* histo = static_cast<TH1*> ( GetSum(physSel,trigClassName,centralityClasses->At(icent)->GetName(),Form("hEventPlane%s", selectEP.Data())) );
    if ( ! histo ) continue;
    histo->SetName(Form("%s_%s", histo->GetName(), centralityClasses->At(icent)->GetName()));
    if ( histo->Integral() < 50. ) continue;
    if ( histo->GetSumw2N() == 0 ) histo->Sumw2();
    histo->Fit("pol0","R0");
    TF1* pol0 = histo->GetFunction("pol0");
    histo->Scale(1./pol0->GetParameter(0));
    Int_t offset = ( addOffsetToCentrality ) ? icolor : 0;
    TGraphErrors* graph = new TGraphErrors();
    graph->SetTitle(Form("%s %s", histo->GetTitle(), trigClassName.Data()));
    for ( Int_t ipoint=0; ipoint<histo->GetXaxis()->GetNbins(); ipoint++ ) {
      Int_t ibin = ipoint+1;
      graph->SetPoint(ipoint, histo->GetXaxis()->GetBinCenter(ibin), histo->GetBinContent(ibin) + offset);
      graph->SetPointError(ipoint, histo->GetXaxis()->GetBinWidth(ibin)/2., histo->GetBinError(ibin));
    }
    graph->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    graph->GetYaxis()->SetTitle("Data/Fit (pol0)");
    TString legTitle = Form("%s %%", centralityClasses->At(icent)->GetName());
    if ( addOffsetToCentrality ) {
      graph->GetYaxis()->SetRangeUser(0, 1.5*(Double_t)centralityClasses->GetEntries());
      graph->GetYaxis()->SetNdivisions(2*centralityClasses->GetEntries(),0,0);
      legTitle += Form(" ( + %i)", offset);
    }
    Int_t currColor = ( icolor < kNtrackSources ) ? srcColors[icolor] : 20+icolor;
    graph->SetLineColor(currColor);
    graph->SetMarkerColor(currColor);
    graph->SetMarkerStyle(20+icolor);
    TString currDrawOpt = "pz";
    if ( icolor == 0 ) currDrawOpt += "a";
    graph->Draw(currDrawOpt.Data());
    legTitle.ReplaceAll("_"," - ");
    leg->AddEntry(graph,legTitle.Data(),"lp");
    if ( addOffsetToCentrality ) {
      TLatex* latex = new TLatex(histo->GetXaxis()->GetXmax()+0.5*histo->GetXaxis()->GetBinWidth(1), 1.+(Double_t)offset, Form("#chi^{2}/NDF = %.3g", pol0->GetChisquare()/(Double_t)pol0->GetNDF()));
      latex->SetTextColor(currColor);
      latex->SetTextSize(0.025);
      latex->Draw("same");
    }
    icolor++;
  }
  leg->Draw("same");


  //////////////////////
  // Event statistics //
  //////////////////////
  printf("\nTotal analyzed events:\n");
  TString evtSel = Form("trigger:%s", trigClassName.Data());
  fEventCounters->PrintSum(evtSel.Data());
  printf("Physics selected analyzed events:\n");
  evtSel = Form("trigger:%s/selected:yes", trigClassName.Data());
  fEventCounters->PrintSum(evtSel.Data());
  
  TString countPhysSel = "any";
  if ( physSel.Contains(fPhysSelKeys->At(kPhysSelPass)->GetName()) ) countPhysSel = "yes";
  else if ( physSel.Contains(fPhysSelKeys->At(kPhysSelReject)->GetName()) ) countPhysSel="no";
  countPhysSel.Prepend("selected:");
  printf("Analyzed events vs. centrality:\n");
  evtSel = Form("trigger:%s/%s", trigClassName.Data(), countPhysSel.Data());
  fEventCounters->Print("centrality",evtSel.Data(),kTRUE);

  if ( ! outFilename.IsNull() ) {
    printf("\nWriting output file %s\n", outFilename.Data());
    TFile* file = TFile::Open(outFilename.Data(), "RECREATE");
    outList.Write();
    file->Close();
  }
}
