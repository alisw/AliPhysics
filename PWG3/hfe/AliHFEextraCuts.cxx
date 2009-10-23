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
//
// Extra cuts implemented by the ALICE Heavy Flavour Electron Group
// Cuts stored here:
// - ITS pixels
// - TPC cluster ratio
// - TRD tracklets
//
// Authors:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TClass.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>

#include "AliESDtrack.h"
#include "AliMCParticle.h"

#include "AliHFEextraCuts.h"

ClassImp(AliHFEextraCuts)

//______________________________________________________
AliHFEextraCuts::AliHFEextraCuts(const Char_t *name, const Char_t *title):
  AliCFCutBase(name, title),
  fCutCorrelation(0),
  fRequirements(0),
  fClusterRatioTPC(0.),
  fMinTrackletsTRD(0),
  fPixelITS(0),
  fCheck(kTRUE),
  fQAlist(0x0),
  fDebugLevel(0)
{
  //
  // Default Constructor
  //
  memset(fImpactParamCut, 0, sizeof(Float_t) * 4);
}

//______________________________________________________
AliHFEextraCuts::AliHFEextraCuts(const AliHFEextraCuts &c):
  AliCFCutBase(c),
  fCutCorrelation(c.fCutCorrelation),
  fRequirements(c.fRequirements),
  fClusterRatioTPC(c.fClusterRatioTPC),
  fMinTrackletsTRD(c.fMinTrackletsTRD),
  fPixelITS(c.fPixelITS),
  fCheck(c.fCheck),
  fQAlist(0x0),
  fDebugLevel(c.fDebugLevel)
{
  //
  // Copy constructor
  // Performs a deep copy
  //
  memcpy(fImpactParamCut, c.fImpactParamCut, sizeof(Float_t) * 4);
  if(IsQAOn()){
    fIsQAOn = kTRUE;
    fQAlist = dynamic_cast<TList *>(c.fQAlist->Clone());
    fQAlist->SetOwner();
  }
}

//______________________________________________________
AliHFEextraCuts &AliHFEextraCuts::operator=(const AliHFEextraCuts &c){
  //
  // Assignment operator
  //
  if(this != &c){
    AliCFCutBase::operator=(c);
    fCutCorrelation = c.fCutCorrelation;
    fRequirements = c.fRequirements;
    fClusterRatioTPC = c.fClusterRatioTPC;
    fMinTrackletsTRD = c.fMinTrackletsTRD;
    fPixelITS = c.fPixelITS;
    fCheck = c.fCheck;
    fDebugLevel = c.fDebugLevel;

    memcpy(fImpactParamCut, c.fImpactParamCut, sizeof(Float_t) * 4);
    if(IsQAOn()){
      fIsQAOn = kTRUE;
      fQAlist = dynamic_cast<TList *>(c.fQAlist->Clone());
      fQAlist->SetOwner();
    }else fQAlist = 0x0;
  }
  return *this;
}

//______________________________________________________
AliHFEextraCuts::~AliHFEextraCuts(){
  //
  // Destructor
  //
  if(fQAlist){
    fQAlist->Delete();
    delete fQAlist;
  }
}

//______________________________________________________
Bool_t AliHFEextraCuts::IsSelected(TObject *o){
  //
  // Steering function for the track selection
  //
  if(TString(o->IsA()->GetName()).CompareTo("AliESDtrack") == 0){
    return CheckESDCuts(dynamic_cast<AliESDtrack *>(o));
  }
  return CheckMCCuts(dynamic_cast<AliMCParticle *>(o));
}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckESDCuts(AliESDtrack *track){
  //
  // Checks cuts on reconstructed tracks
  // returns true if track is selected
  // QA histograms are filled before track selection and for
  // selected tracks after track selection
  //
  ULong64_t survivedCut = 0;	// Bitmap for cuts which are passed by the track, later to be compared with fRequirements
  if(IsQAOn()) FillQAhistosESD(track, kBeforeCuts);
  // Apply cuts
  Float_t impactR, impactZ, ratioTPC;
  track->GetImpactParameters(impactR, impactZ);
  // printf("Check TPC findable clusters: %d, found Clusters: %d\n", track->GetTPCNclsF(), track->GetTPCNcls());
  ratioTPC = track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t>(track->GetTPCNclsF()) : 1.;
  UChar_t trdTracklets;
  trdTracklets = track->GetTRDntrackletsPID();
  UChar_t itsPixel = track->GetITSClusterMap();
  Int_t det, status1, status2;
  Float_t xloc, zloc;
  track->GetITSModuleIndexInfo(0, det, status1, xloc, zloc);
  track->GetITSModuleIndexInfo(1, det, status2, xloc, zloc);
  if(TESTBIT(fRequirements, kMinImpactParamR)){
    // cut on min. Impact Parameter in Radial direction
    if(TMath::Abs(impactR) >= fImpactParamCut[0]) SETBIT(survivedCut, kMinImpactParamR);
  }
  if(TESTBIT(fRequirements, kMinImpactParamZ)){
    // cut on min. Impact Parameter in Z direction
    if(TMath::Abs(impactZ) >= fImpactParamCut[1]) SETBIT(survivedCut, kMinImpactParamZ);
  }
  if(TESTBIT(fRequirements, kMaxImpactParamR)){
    // cut on max. Impact Parameter in Radial direction
    if(TMath::Abs(impactR) <= fImpactParamCut[2]) SETBIT(survivedCut, kMaxImpactParamR);
  }
  if(TESTBIT(fRequirements, kMaxImpactParamZ)){
    // cut on max. Impact Parameter in Z direction
    if(TMath::Abs(impactZ) <= fImpactParamCut[3]) SETBIT(survivedCut, kMaxImpactParamZ);
  }
  if(TESTBIT(fRequirements, kClusterRatioTPC)){
    // cut on min ratio of found TPC clusters vs findable TPC clusters
    if(ratioTPC >= fClusterRatioTPC) SETBIT(survivedCut, kClusterRatioTPC);
  }
  if(TESTBIT(fRequirements, kMinTrackletsTRD)){
    // cut on minimum number of TRD tracklets
    if(fDebugLevel > 0){
      printf("Min TRD cut: [%d|%d]\n", fMinTrackletsTRD, trdTracklets);
    }
    if(trdTracklets >= fMinTrackletsTRD) SETBIT(survivedCut, kMinTrackletsTRD);
  }
  if(TESTBIT(fRequirements, kPixelITS)){
    // cut on ITS pixel layers
    if(fDebugLevel > 0){
      printf("ITS cluster Map: ");
      PrintBitMap(itsPixel);
    }
    switch(fPixelITS){
      case kFirst: 
	if(!TESTBIT(itsPixel, 0)) {
	  if(fCheck){
	    if(!CheckITSstatus(status1)) {
	      SETBIT(survivedCut, kPixelITS);
	    }
	  }
	}
	else {
	  SETBIT(survivedCut, kPixelITS);
	}
		    break;
      case kSecond: 
	if(!TESTBIT(itsPixel, 1)) {
	  if(fCheck) {
	    if(!CheckITSstatus(status2)) {
	      SETBIT(survivedCut, kPixelITS);
	    }
	  }
	}
	else { 
	  SETBIT(survivedCut, kPixelITS);
	}
		    break;
      case kBoth: 
	if(!(TESTBIT(track->GetITSClusterMap(),0))) {
	  if(fCheck) {  
	    if(!CheckITSstatus(status1)) {
	      if(!(TESTBIT(track->GetITSClusterMap(),1))) {
		if(!CheckITSstatus(status2)) {
		  SETBIT(survivedCut, kPixelITS);
		}
	      }
	      else SETBIT(survivedCut, kPixelITS);
	    }
	  }
	}
	else {
	  
	  if(!(TESTBIT(track->GetITSClusterMap(),1))) {
	    if(fCheck) {
	      if(!CheckITSstatus(status2)) {
		SETBIT(survivedCut, kPixelITS);
	      }
	    }
	  }
	  else SETBIT(survivedCut, kPixelITS);
	
	}
	           break;
      case kAny: 
	if((!TESTBIT(itsPixel, 0)) && (!TESTBIT(itsPixel, 1))){
	  if(fCheck){
	    if(!CheckITSstatus(status1) || (!CheckITSstatus(status2))) {
	      SETBIT(survivedCut, kPixelITS);
	    }
	  }
	}
	else { 
	  SETBIT(survivedCut, kPixelITS);
	}
		    break;
      default: break;
    }
    if(fDebugLevel > 0)printf("Survived Cut: %s\n", TESTBIT(survivedCut, kPixelITS) ? "YES" : "NO");
  }
  if(fRequirements == survivedCut){
    //
    // Track selected
    //
    if(IsQAOn()) FillQAhistosESD(track, kAfterCuts);
    return kTRUE;
  }
  if(IsQAOn()) FillCutCorrelation(survivedCut);
  return kFALSE;
}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckMCCuts(AliMCParticle */*track*/) const {
  //
  // Checks cuts on Monte Carlo tracks
  // returns true if track is selected
  // QA histograms are filled before track selection and for
  // selected tracks after track selection
  //
  return kTRUE;	// not yet implemented
}

//______________________________________________________
void AliHFEextraCuts::FillQAhistosESD(AliESDtrack *track, UInt_t when){
  //
  // Fill the QA histograms for ESD tracks
  // Function can be called before cuts or after cut application (second argument)
  //
  TList *container = dynamic_cast<TList *>(fQAlist->At(when));
  Float_t impactR, impactZ;
  track->GetImpactParameters(impactR, impactZ);
  (dynamic_cast<TH1F *>(container->At(0)))->Fill(impactR);
  (dynamic_cast<TH1F *>(container->At(1)))->Fill(impactZ);
  // printf("TPC findable clusters: %d, found Clusters: %d\n", track->GetTPCNclsF(), track->GetTPCNcls());
  (dynamic_cast<TH1F *>(container->At(2)))->Fill(track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t>(track->GetTPCNclsF()) : 1.);
  (dynamic_cast<TH1F *>(container->At(3)))->Fill(track->GetTRDntrackletsPID());
  UChar_t itsPixel = track->GetITSClusterMap();
  TH1 *pixelHist = dynamic_cast<TH1F *>(container->At(4));
  Int_t firstEntry = pixelHist->GetXaxis()->GetFirst();
  if(!((itsPixel & BIT(0)) || (itsPixel & BIT(1))))
    pixelHist->Fill(firstEntry + 3);
  else{
    if(itsPixel & BIT(0)){
      pixelHist->Fill(firstEntry);
      if(itsPixel & BIT(1)) pixelHist->Fill(firstEntry + 2);
      else pixelHist->Fill(firstEntry + 4);
    }
    if(itsPixel & BIT(1)){
      pixelHist->Fill(firstEntry + 1);
      if(!(itsPixel & BIT(0))) pixelHist->Fill(firstEntry + 5);
    }
  }
}

// //______________________________________________________
// void AliHFEextraCuts::FillQAhistosMC(AliMCParticle *track, UInt_t when){
//   //
//   // Fill the QA histograms for MC tracks
//   // Function can be called before cuts or after cut application (second argument)
//   // Not yet implemented
//   //
// }

//______________________________________________________
void AliHFEextraCuts::FillCutCorrelation(ULong64_t survivedCut){
  //
  // Fill cut correlation histograms for tracks that didn't pass cuts
  //
  TH2 *correlation = dynamic_cast<TH2F *>(fQAlist->At(2));
  for(Int_t icut = 0; icut < kNcuts; icut++){
    if(!TESTBIT(fRequirements, icut)) continue;
    for(Int_t jcut = icut; jcut < kNcuts; jcut++){
      if(!TESTBIT(fRequirements, jcut)) continue;
      if(TESTBIT(survivedCut, icut) && TESTBIT(survivedCut, jcut))
	correlation->Fill(icut, jcut);
    }
  }
}

//______________________________________________________
void AliHFEextraCuts::AddQAHistograms(TList *qaList){
  //
  // Add QA histograms
  // For each cut a histogram before and after track cut is created
  // Histos before respectively after cut are stored in different lists
  // Additionally a histogram with the cut correlation is created and stored
  // in the top directory
  //
  TList *histos[2];
  TH1 *histo1D = 0x0;
  TH2 *histo2D = 0x0;
  histos[0] = new TList();
  histos[0]->SetName("BeforeCut");
  histos[0]->SetOwner();
  histos[1] = new TList();
  histos[1]->SetName("AfterCut");
  histos[1]->SetOwner();
  TString cutstr[2] = {"before", "after"};
  for(Int_t icond = 0; icond < 2; icond++){
    histos[icond]->AddAt((histo1D = new TH1F(Form("impactParamR%s", cutstr[icond].Data()), "Radial Impact Parameter", 100, 0, 10)), 0);
    histo1D->GetXaxis()->SetTitle("Impact Parameter");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    histos[icond]->AddAt((histo1D = new TH1F(Form("impactParamZ%s", cutstr[icond].Data()), "Z Impact Parameter", 200, 0, 20)), 1);
    histo1D->GetXaxis()->SetTitle("Impact Parameter");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    histos[icond]->AddAt((histo1D = new TH1F(Form("tpcClr%s", cutstr[icond].Data()), "Cluster Ratio TPC", 10, 0, 1)), 2);
    histo1D->GetXaxis()->SetTitle("Cluster Ratio TPC");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    histos[icond]->AddAt((histo1D = new TH1F(Form("trdTracklets%s", cutstr[icond].Data()), "Number of TRD tracklets", 7, 0, 7)), 3);
    histo1D->GetXaxis()->SetTitle("Number of TRD Tracklets");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    histos[icond]->AddAt((histo1D = new TH1F(Form("itsPixel%s", cutstr[icond].Data()), "ITS Pixel Hits", 6, 0, 6)), 4);
    histo1D->GetXaxis()->SetTitle("ITS Pixel");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    Int_t first = histo1D->GetXaxis()->GetFirst();
    TString binNames[6] = { "First", "Second", "Both", "None", "Exclusive First", "Exclusive Second"};
    for(Int_t ilabel = 0; ilabel < 6; ilabel++)
      histo1D->GetXaxis()->SetBinLabel(first + ilabel, binNames[ilabel].Data());
  }
  fQAlist = new TList();
  fQAlist->SetOwner();
  fQAlist->SetName("HFelectronExtraCuts");
  fQAlist->AddAt(histos[0], 0);
  fQAlist->AddAt(histos[1], 1);
  // Add cut correlation
  fQAlist->AddAt((histo2D = new TH2F("cutcorrellation", "Cut Correlation", kNcuts, 0, kNcuts - 1, kNcuts, 0, kNcuts -1)), 2);
  TString labels[kNcuts] = {"MinImpactParamR", "MaxImpactParamR", "MinImpactParamZ", "MaxImpactParamZ", "ClusterRatioTPC", "MinTrackletsTRD", "ITSpixel"};
  Int_t firstx = histo2D->GetXaxis()->GetFirst(), firsty = histo2D->GetYaxis()->GetFirst();
  for(Int_t icut = 0; icut < kNcuts; icut++){
    histo2D->GetXaxis()->SetBinLabel(firstx + icut, labels[icut].Data());
    histo2D->GetYaxis()->SetBinLabel(firsty + icut, labels[icut].Data());
  }
  qaList->AddLast(fQAlist);
}

//______________________________________________________
void AliHFEextraCuts::PrintBitMap(Int_t bitmap){
  for(Int_t ibit = 32; ibit--; )
    printf("%d", bitmap & BIT(ibit) ? 1 : 0);
  printf("\n");
}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckITSstatus(Int_t itsStatus) const {
  //
  // Check whether ITS area is dead
  //
  Bool_t status;
  switch(itsStatus){
    case 2: status = kFALSE; break;
    case 3: status = kFALSE; break;
    case 7: status = kFALSE; break;
    default: status = kTRUE;
  }
  return status;
}
