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
#include <TBits.h>
#include <TClass.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>

#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliVertexerTracks.h"
#include "AliVVertex.h"

#include "AliHFEextraCuts.h"

ClassImp(AliHFEextraCuts)

const Int_t AliHFEextraCuts::fgkNQAhistos = 7;

//______________________________________________________
AliHFEextraCuts::AliHFEextraCuts(const Char_t *name, const Char_t *title):
  AliCFCutBase(name, title),
  fEvent(NULL),
  fCutCorrelation(0),
  fRequirements(0),
  fMinNClustersTPC(0),
  fClusterRatioTPC(0.),
  fMinTrackletsTRD(0),
  fPixelITS(0),
  fTPCclusterDef(0),
  fTPCclusterRatioDef(0),
  fFractionTPCShared(-1.0),
  fCheck(kFALSE),
  fQAlist(0x0) ,
  fDebugLevel(0)
{
  //
  // Default Constructor
  //
  memset(fImpactParamCut, 0, sizeof(Float_t) * 4);
  memset(fIPcutParam, 0, sizeof(Float_t) * 4);
}

//______________________________________________________
AliHFEextraCuts::AliHFEextraCuts(const AliHFEextraCuts &c):
  AliCFCutBase(c),
  fEvent(c.fEvent),
  fCutCorrelation(c.fCutCorrelation),
  fRequirements(c.fRequirements),
  fMinNClustersTPC(c.fMinNClustersTPC),
  fClusterRatioTPC(c.fClusterRatioTPC),
  fMinTrackletsTRD(c.fMinTrackletsTRD),
  fPixelITS(c.fPixelITS),
  fTPCclusterDef(c.fTPCclusterDef),
  fTPCclusterRatioDef(c.fTPCclusterRatioDef),
  fFractionTPCShared(c.fFractionTPCShared),
  fCheck(c.fCheck),
  fQAlist(0x0),
  fDebugLevel(0)
{
  //
  // Copy constructor
  // Performs a deep copy
  //
  memcpy(fImpactParamCut, c.fImpactParamCut, sizeof(Float_t) * 4);
  memcpy(fIPcutParam, c.fIPcutParam, sizeof(Float_t) * 4);
  if(IsQAOn()){
    fIsQAOn = kTRUE;
    fQAlist = dynamic_cast<TList *>(c.fQAlist->Clone());
    if(fQAlist) fQAlist->SetOwner();
  }
}

//______________________________________________________
AliHFEextraCuts &AliHFEextraCuts::operator=(const AliHFEextraCuts &c){
  //
  // Assignment operator
  //
  if(this != &c){
    AliCFCutBase::operator=(c);
    fEvent = c.fEvent;
    fCutCorrelation = c.fCutCorrelation;
    fRequirements = c.fRequirements;
    fClusterRatioTPC = c.fClusterRatioTPC;
    fMinNClustersTPC = c.fMinNClustersTPC;
    fMinTrackletsTRD = c.fMinTrackletsTRD;
    fPixelITS = c.fPixelITS;
    fTPCclusterDef = c.fTPCclusterDef;
    fTPCclusterRatioDef = c.fTPCclusterRatioDef;
    fFractionTPCShared = c.fFractionTPCShared;
    fCheck = c.fCheck;
    fDebugLevel = c.fDebugLevel;
    memcpy(fImpactParamCut, c.fImpactParamCut, sizeof(Float_t) * 4);
    memcpy(fIPcutParam, c.fIPcutParam, sizeof(Float_t) * 4);
    if(IsQAOn()){
      fIsQAOn = kTRUE;
      fQAlist = dynamic_cast<TList *>(c.fQAlist->Clone());
      if(fQAlist) fQAlist->SetOwner();
    }else fQAlist = 0x0;
  }
  return *this;
}

//______________________________________________________
AliHFEextraCuts::~AliHFEextraCuts(){
  //
  // Destructor
  //
  if(fQAlist) delete fQAlist;
}

//______________________________________________________
void AliHFEextraCuts::SetRecEventInfo(const TObject *event){
  //
  // Set Virtual event an make a copy
  //
  if (!event) {
    AliError("Pointer to AliVEvent !");
    return;
  }
  TString className(event->ClassName());
  if (! (className.CompareTo("AliESDEvent")==0 || className.CompareTo("AliAODEvent")==0)) {
    AliError("argument must point to an AliESDEvent or AliAODEvent !");
    return ;
  }
  fEvent = (AliVEvent*) event;

}

//______________________________________________________
Bool_t AliHFEextraCuts::IsSelected(TObject *o){
  //
  // Steering function for the track selection
  //
  TString type = o->IsA()->GetName();
  AliDebug(2, Form("Object type %s", type.Data()));
  if(!type.CompareTo("AliESDtrack") || !type.CompareTo("AliAODTrack")){
    return CheckRecCuts(dynamic_cast<AliVTrack *>(o));
  }
  return CheckMCCuts(dynamic_cast<AliVParticle *>(o));
}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckRecCuts(AliVTrack *track){
  //
  // Checks cuts on reconstructed tracks
  // returns true if track is selected
  // QA histograms are filled before track selection and for
  // selected tracks after track selection
  //
  AliDebug(1, "Called");
  ULong64_t survivedCut = 0;	// Bitmap for cuts which are passed by the track, later to be compared with fRequirements
  if(IsQAOn()) FillQAhistosRec(track, kBeforeCuts);
  // Apply cuts
  Float_t impactR, impactZ;
  Double_t hfeimpactR, hfeimpactnsigmaR;
  Double_t hfeimpactRcut, hfeimpactnsigmaRcut;
  Double_t maximpactRcut; 
  GetImpactParameters(track, impactR, impactZ);
  if(TESTBIT(fRequirements, kMinHFEImpactParamR) || TESTBIT(fRequirements, kMinHFEImpactParamNsigmaR)){
    // Protection for PbPb
    GetHFEImpactParameterCuts(track, hfeimpactRcut, hfeimpactnsigmaRcut);
    GetHFEImpactParameters(track, hfeimpactR, hfeimpactnsigmaR);
  }
  UInt_t nclsTPC = GetTPCncls(track);
  // printf("Check TPC findable clusters: %d, found Clusters: %d\n", track->GetTPCNclsF(), track->GetTPCNcls());
  Float_t fractionSharedClustersTPC = GetTPCsharedClustersRatio(track);
  Double_t ratioTPC = GetTPCclusterRatio(track);
  UChar_t trdTracklets;
  trdTracklets = track->GetTRDntrackletsPID();
  UChar_t itsPixel = track->GetITSClusterMap();
  Int_t status1 = GetITSstatus(track, 0);
  Int_t status2 = GetITSstatus(track, 1);
  Bool_t statusL0 = CheckITSstatus(status1);
  Bool_t statusL1 = CheckITSstatus(status2);
  if(TESTBIT(fRequirements, kTPCfractionShared)) {
    if(TMath::Abs(fractionSharedClustersTPC) >= fFractionTPCShared) SETBIT(survivedCut, kTPCfractionShared);    
  }
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
  if(TESTBIT(fRequirements, kMaxImpactParameterRpar)) {
    // HFE Impact parameter cut
    GetMaxImpactParameterCutR(track,maximpactRcut);
    if(TMath::Abs(impactR) >= maximpactRcut) SETBIT(fRequirements, kMaxImpactParameterRpar);
  }
  if(TESTBIT(fRequirements, kMaxImpactParamZ)){
    // cut on max. Impact Parameter in Z direction
    if(TMath::Abs(impactZ) <= fImpactParamCut[3]) SETBIT(survivedCut, kMaxImpactParamZ);
  }
  if(TESTBIT(fRequirements, kMinHFEImpactParamR)){
    // cut on min. HFE Impact Parameter in Radial direction
    if(TMath::Abs(hfeimpactR) >= hfeimpactRcut) SETBIT(survivedCut, kMinHFEImpactParamR);
  }
  if(TESTBIT(fRequirements, kMinHFEImpactParamNsigmaR)){
    // cut on max. HFE Impact Parameter n sigma in Radial direction
    if(TMath::Abs(hfeimpactnsigmaR) >= hfeimpactnsigmaRcut) SETBIT(survivedCut, kMinHFEImpactParamNsigmaR);
  }
  if(TESTBIT(fRequirements, kClusterRatioTPC)){
    // cut on min ratio of found TPC clusters vs findable TPC clusters
    if(ratioTPC >= fClusterRatioTPC) SETBIT(survivedCut, kClusterRatioTPC);
  }
  if(TESTBIT(fRequirements, kMinTrackletsTRD)){
    // cut on minimum number of TRD tracklets
    AliDebug(1, Form("Min TRD cut: [%d|%d]\n", fMinTrackletsTRD, trdTracklets));
    if(trdTracklets >= fMinTrackletsTRD) SETBIT(survivedCut, kMinTrackletsTRD);
  }
  if(TESTBIT(fRequirements, kMinNClustersTPC)){
    // cut on minimum number of TRD tracklets
    AliDebug(1, Form("Min TPC cut: [%d|%d]\n", fMinNClustersTPC, nclsTPC));
    if(nclsTPC >= fMinNClustersTPC) SETBIT(survivedCut, kMinNClustersTPC);
  }
  if(TESTBIT(fRequirements, kPixelITS)){
    // cut on ITS pixel layers
    AliDebug(1, "ITS cluster Map: ");
    //PrintBitMap(itsPixel);
    switch(fPixelITS){
      case kFirst:
        AliDebug(2, "First");
	      if(TESTBIT(itsPixel, 0)) 
	        SETBIT(survivedCut, kPixelITS);
        else
	        if(fCheck && !statusL0)
	          SETBIT(survivedCut, kPixelITS);
		    break;
      case kSecond: 
        AliDebug(2, "Second");
	      if(TESTBIT(itsPixel, 1))
	        SETBIT(survivedCut, kPixelITS);
        else
	        if(fCheck && !statusL1)
	          SETBIT(survivedCut, kPixelITS);
		    break;
      case kBoth: 
        AliDebug(2, "Both");
	      if(TESTBIT(itsPixel,0) && TESTBIT(itsPixel,1))
		      SETBIT(survivedCut, kPixelITS);
	      else
          if(fCheck && !(statusL0 && statusL1)) 
		        SETBIT(survivedCut, kPixelITS);
	      break;
      case kAny: 
        AliDebug(2, "Any");
	      if(TESTBIT(itsPixel, 0) || TESTBIT(itsPixel, 1))
	        SETBIT(survivedCut, kPixelITS);
        else
	        if(fCheck && !(statusL0 || statusL1))
	            SETBIT(survivedCut, kPixelITS);
		    break;
      case kExclusiveSecond:
        AliDebug(2, "Exlusive second");
        if(fCheck){ // Cut out tracks which pass a dead ITS Layer 0
          if(TESTBIT(itsPixel,1) && !TESTBIT(itsPixel,0) && statusL0)
            SETBIT(survivedCut, kPixelITS);
        } else {
          if(TESTBIT(itsPixel,1) && !TESTBIT(itsPixel,0))
            SETBIT(survivedCut, kPixelITS);
        }
        break;
      default: 
        AliDebug(2, "None");
        break;
    }
    AliDebug(1, Form("Survived Cut: %s\n", TESTBIT(survivedCut, kPixelITS) ? "YES" : "NO"));
  }

  if(TESTBIT(fRequirements, kTOFPID)){
    // cut on TOF pid
    if(track->GetStatus() & AliESDtrack::kTOFpid) SETBIT(survivedCut, kTOFPID);
  }

  if(TESTBIT(fRequirements, kTOFmismatch)){
    // cut on TOF mismatch
    if(!(track->GetStatus() & AliESDtrack::kTOFmismatch)) SETBIT(survivedCut, kTOFmismatch);
  }

  if(TESTBIT(fRequirements, kTPCPIDCleanUp)){
      // cut on TPC PID cleanup
      Int_t fClusterdEdx=GetTPCnclusdEdx(track);
      Bool_t fBitsAboveThreshold=GetTPCCountSharedMapBitsAboveThreshold(track);
      if((fBitsAboveThreshold==0)&&(fClusterdEdx>80)) SETBIT(survivedCut, kTPCPIDCleanUp);
  }

  if(TESTBIT(fRequirements, kEMCALmatch)){
    if(track->GetStatus() && AliESDtrack::kEMCALmatch) SETBIT(survivedCut, kEMCALmatch);
  }

  if(fRequirements == survivedCut){
    //
    // Track selected
    //
    AliDebug(2, "Track Survived cuts\n");
    if(IsQAOn()) FillQAhistosRec(track, kAfterCuts);
    return kTRUE;
  }
  AliDebug(2, "Track cut");
  if(IsQAOn()) FillCutCorrelation(survivedCut);
  return kFALSE;
}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckMCCuts(AliVParticle */*track*/) const {
  //
  // Checks cuts on Monte Carlo tracks
  // returns true if track is selected
  // QA histograms are filled before track selection and for
  // selected tracks after track selection
  //
  return kTRUE;	// not yet implemented
}

//______________________________________________________
void AliHFEextraCuts::FillQAhistosRec(AliVTrack *track, UInt_t when){
  //
  // Fill the QA histograms for ESD tracks
  // Function can be called before cuts or after cut application (second argument)
  //
  Float_t impactR, impactZ;
  GetImpactParameters(track, impactR, impactZ);
  TH1 *htmp = NULL;
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(0 + when * fgkNQAhistos)))) htmp->Fill(impactR);
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(1 + when * fgkNQAhistos)))) htmp->Fill(impactZ);
  // printf("TPC findable clusters: %d, found Clusters: %d\n", track->GetTPCNclsF(), track->GetTPCNcls());
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(2 + when * fgkNQAhistos)))) htmp->Fill(GetTPCclusterRatio(track));
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(3 + when * fgkNQAhistos)))) htmp->Fill(track->GetTRDntrackletsPID());
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(4 + when * fgkNQAhistos)))) htmp->Fill(GetTPCncls(track));
  UChar_t itsPixel = track->GetITSClusterMap();
  TH1 *pixelHist = dynamic_cast<TH1F *>(fQAlist->At(5 + when * fgkNQAhistos));
  //Int_t firstEntry = pixelHist->GetXaxis()->GetFirst();
  if(pixelHist){
    Double_t firstEntry = 0.5;
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
  // Fill histogram with the status bits
  TH1 *hStatusBits = dynamic_cast<TH1 *>(fQAlist->At(6 + when * fgkNQAhistos));
  hStatusBits->Fill(0);  // Fill first bin with all tracks
  if(track->GetStatus() && AliESDtrack::kTOFpid) hStatusBits->Fill(1);
  if(!(track->GetStatus() && AliESDtrack::kTOFmismatch)) hStatusBits->Fill(2);
  if(track->GetStatus() && AliESDtrack::kEMCALmatch) hStatusBits->Fill(3);
  if(GetTPCCountSharedMapBitsAboveThreshold(track)==0) hStatusBits->Fill(4);
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(7 + when * fgkNQAhistos)))) htmp->Fill(GetTPCnclusdEdx(track));
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
  TH2 *correlation = dynamic_cast<TH2F *>(fQAlist->At(2 * fgkNQAhistos));
  if(correlation){
    for(Int_t icut = 0; icut < kNcuts; icut++){
      if(!TESTBIT(fRequirements, icut)) continue;
      for(Int_t jcut = icut; jcut < kNcuts; jcut++){
        if(!TESTBIT(fRequirements, jcut)) continue;
        if(TESTBIT(survivedCut, icut) && TESTBIT(survivedCut, jcut))
	        correlation->Fill(icut, jcut);
      }
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

  TH1 *histo1D = 0x0;
  TH2 *histo2D = 0x0;
  TString cutstr[2] = {"before", "after"};

  if(!fQAlist) fQAlist = new TList;  // for internal representation, not owner
  for(Int_t icond = 0; icond < 2; icond++){
    qaList->AddAt((histo1D = new TH1F(Form("%s_impactParamR%s",GetName(),cutstr[icond].Data()), "Radial Impact Parameter", 100, 0, 10)), 0 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 0 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Impact Parameter");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_impactParamZ%s",GetName(),cutstr[icond].Data()), "Z Impact Parameter", 200, 0, 20)), 1 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 1 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Impact Parameter");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_tpcClr%s",GetName(),cutstr[icond].Data()), "Cluster Ratio TPC", 100, 0, 1)), 2 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 2 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Cluster Ratio TPC");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_trdTracklets%s",GetName(),cutstr[icond].Data()), "Number of TRD tracklets", 7, 0, 7)), 3 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 3 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Number of TRD Tracklets");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_tpcClusters%s",GetName(),cutstr[icond].Data()), "Number of TPC clusters", 161, 0, 160)), 4 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 4 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Number of TPC clusters");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_itsPixel%s",GetName(),cutstr[icond].Data()), "ITS Pixel Hits", 6, 0, 6)), 5 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 5 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("ITS Pixel");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    Int_t first = histo1D->GetXaxis()->GetFirst();
    TString binNames[6] = { "First", "Second", "Both", "None", "Exclusive First", "Exclusive Second"};
    for(Int_t ilabel = 0; ilabel < 6; ilabel++)
      histo1D->GetXaxis()->SetBinLabel(first + ilabel, binNames[ilabel].Data());
    qaList->AddAt((histo1D = new TH1F(Form("%s_trackStatus%s",GetName(),cutstr[icond].Data()), "Track Status", 5, 0, 5)), 6 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 6 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Track Status Bit");
    histo1D->GetYaxis()->SetTitle("Number of tracks");
    TString tsnames[5] = {"All", "TOFPID", "No TOFmismatch", "EMCALmatch","No TPC shared bit"};
    for(Int_t istat = 0; istat < 5; istat++) histo1D->GetXaxis()->SetBinLabel(istat + 1, tsnames[istat].Data());
    qaList->AddAt((histo1D = new TH1F(Form("%s_tpcdEdxClusters%s",GetName(),cutstr[icond].Data()), "Number of TPC clusters for dEdx calculation", 161, 0, 160)), 7 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 7 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Number of TPC clusters for dEdx calculation");
    histo1D->GetYaxis()->SetTitle("counts");
  }
  // Add cut correlation
  qaList->AddAt((histo2D = new TH2F(Form("%s_cutcorrelation",GetName()), "Cut Correlation", kNcuts, 0, kNcuts - 1, kNcuts, 0, kNcuts -1)), 2 * fgkNQAhistos);
  fQAlist->AddAt(histo2D, 2 * fgkNQAhistos);
  TString labels[kNcuts] = {"MinImpactParamR", "MaxImpactParamR", "MinImpactParamZ", "MaxImpactParamZ", "ClusterRatioTPC", "MinTrackletsTRD", "ITSpixel", "kMinHFEImpactParamR", "kMinHFEImpactParamNsigmaR", "TPC Number of clusters" "Fraction Shared TPC clusters", "TOFPID", "No TOFmismatch", "EMCALmatch", "ImpactParam"};
  Int_t firstx = histo2D->GetXaxis()->GetFirst(), firsty = histo2D->GetYaxis()->GetFirst();
  for(Int_t icut = 0; icut < kNcuts; icut++){
    histo2D->GetXaxis()->SetBinLabel(firstx + icut, labels[icut].Data());
    histo2D->GetYaxis()->SetBinLabel(firsty + icut, labels[icut].Data());
  }
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

//______________________________________________________
Int_t AliHFEextraCuts::GetITSstatus(AliVTrack *track, Int_t layer){
	//
	// Check ITS layer status
	//
	Int_t status = 0;
	if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
		Int_t det;
		Float_t xloc, zloc;
		AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
		if(esdtrack) esdtrack->GetITSModuleIndexInfo(layer, det, status, xloc, zloc);
	}
	return status;
}


//______________________________________________________
Bool_t AliHFEextraCuts::GetTPCCountSharedMapBitsAboveThreshold(AliVTrack *track){
  //
  // Checks if number of shared bits is above 1
  //
  Int_t nsharebit = 1;
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack){ // coverity
	TBits shared = esdtrack->GetTPCSharedMap();
	if(shared.CountBits() >= 2) nsharebit=1;
        else nsharebit=0;
    }
  }
  return nsharebit;
}



//______________________________________________________
UInt_t AliHFEextraCuts::GetTPCncls(AliVTrack *track){
	//
	// Get Number of findable clusters in the TPC
	//
	Int_t nClusters = 0; // in case no Information available consider all clusters findable
	TString type = track->IsA()->GetName();
	if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack){ // coverity
      if(TESTBIT(fTPCclusterDef, kFoundIter1)){
		    nClusters = esdtrack->GetTPCNclsIter1();
        AliDebug(2, ("Using def kFoundIter1"));
      } else if(TESTBIT(fTPCclusterDef, kCrossedRows)){
        AliDebug(2, ("Using def kCrossedRows"));
        nClusters = static_cast<UInt_t>(esdtrack->GetTPCClusterInfo(2,1));
      } else{
        AliDebug(2, ("Using def kFound"));
		    nClusters = esdtrack->GetTPCNcls();
      }
    }
  }
	else if(!type.CompareTo("AliAODTrack")){
		AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if(aodtrack){
	    const TBits &tpcmap = aodtrack->GetTPCClusterMap();
	    for(UInt_t ibit = 0; ibit < tpcmap.GetNbits(); ibit++)
		    if(tpcmap.TestBitNumber(ibit)) nClusters++;
    }
	}
	return nClusters;
}

//______________________________________________________
UInt_t AliHFEextraCuts::GetTPCnclusdEdx(AliVTrack *track){
  //
  // Get number of TPC cluster used to calculate dEdx
  //
  Int_t nClustersdEdx = 0;
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack){ // coverity
      nClustersdEdx = esdtrack->GetTPCsignalN();
    }
  }
  return nClustersdEdx;
}



//______________________________________________________
Float_t AliHFEextraCuts::GetTPCsharedClustersRatio(AliVTrack *track){
  //
  // Get fraction of shared TPC clusters
  //
  Float_t fracClustersTPCShared = 0.0; 
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack){ // coverity
      Float_t nClustersTPC = esdtrack->GetTPCclusters(0);
      Int_t nClustersTPCShared = esdtrack->GetTPCnclsS();
      if (nClustersTPC!=0) {
	fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
      }
    }
  }
  return fracClustersTPCShared;
}

//______________________________________________________
Double_t AliHFEextraCuts::GetTPCclusterRatio(AliVTrack *track){
  //
  // Get Ratio of found / findable clusters for different definitions
  // Only implemented for ESD tracks
  //
  Double_t clusterRatio = 1.; // in case no Information available consider all clusters findable
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack){ // coverity
      if(TESTBIT(fTPCclusterRatioDef, kCROverFindable)){
        AliDebug(2, "Using ratio def kCROverFindable");
        clusterRatio = esdtrack->GetTPCNclsF() ? esdtrack->GetTPCClusterInfo(2,1)/static_cast<Double_t>(esdtrack->GetTPCNclsF()) : 1.; // crossed rows/findable
      } else if(TESTBIT(fTPCclusterRatioDef, kFoundOverCR)){
        AliDebug(2, "Using ratio def kFoundOverCR");
        clusterRatio = esdtrack->GetTPCClusterInfo(2,0); // found/crossed rows
      } else if(TESTBIT(fTPCclusterRatioDef, kFoundOverFindableIter1)){
        AliDebug(2, "Using ratio def kFoundOverFindableIter1");
        clusterRatio = esdtrack->GetTPCNclsFIter1() ? static_cast<Double_t>(esdtrack->GetTPCNclsIter1())/static_cast<Double_t>(esdtrack->GetTPCNclsFIter1()) : 1.; // crossed
      } else {
        AliDebug(2, "Using ratio def kFoundOverFindable");
        clusterRatio = esdtrack->GetTPCNclsF() ? static_cast<Double_t>(esdtrack->GetTPCNcls())/static_cast<Double_t>(esdtrack->GetTPCNclsF()) : 1.; // found/findable
      }
    }
  }
  return clusterRatio;
}

//______________________________________________________
void AliHFEextraCuts::GetImpactParameters(AliVTrack *track, Float_t &radial, Float_t &z){
	//
	// Get impact parameter
	//
	TString type = track->IsA()->GetName();
	if(!type.CompareTo("AliESDtrack")){
		AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(esdtrack) esdtrack->GetImpactParameters(radial, z);
	}
	else if(!type.CompareTo("AliAODTrack")){
		AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if(aodtrack){
		  Double_t xyz[3];
		  aodtrack->XYZAtDCA(xyz);
		  z = xyz[2];
		  radial = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]+xyz[1]);
    }
	}
}

//______________________________________________________
void AliHFEextraCuts::GetHFEImpactParameters(AliVTrack *track, Double_t &dcaxy, Double_t &dcansigmaxy){
	//
	// Get HFE impact parameter (with recalculated primary vertex)
	//
	dcaxy=0;
	dcansigmaxy=0;
  if(!fEvent){
    AliDebug(1, "No Input event available\n");
    return;
  }
  const Double_t kBeampiperadius=3.;
  TString type = track->IsA()->GetName();
  Double_t dcaD[2]={-999.,-999.};
  Double_t covD[3]={-999.,-999.,-999.};
  Float_t dcaF[2]={-999.,-999.};
  Float_t covF[3]={-999.,-999.,-999.};
  
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fEvent);
  const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
  if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
   // recalculate primary vertex for peri. and pp
   AliVertexerTracks vertexer(fEvent->GetMagneticField());
   vertexer.SetITSMode();
   vertexer.SetMinClusters(4);
	 Int_t skipped[2];
   skipped[0] = track->GetID();
   vertexer.SetSkipTracks(1,skipped);
   
 //----diamond constraint-----------------------------
   vertexer.SetConstraintOn();
   Float_t diamondcovxy[3];
   esdevent->GetDiamondCovXY(diamondcovxy);
   Double_t pos[3]={esdevent->GetDiamondX(),esdevent->GetDiamondY(),0.};
   Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
   AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
   vertexer.SetVtxStart(diamond);
   delete diamond; diamond=NULL;
 //----------------------------------------------------
    
   vtxESDSkip = vertexer.FindPrimaryVertex(fEvent);

   // Getting the DCA
   // Propagation always done on a working copy to not disturb the track params of the original track
   AliESDtrack *esdtrack = NULL;
   if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
     // Case ESD track: take copy constructor
     AliESDtrack *tmptrack = dynamic_cast<AliESDtrack *>(track);
     if(tmptrack) esdtrack = new AliESDtrack(*tmptrack);
   } else {
    // Case AOD track: take different constructor
    esdtrack = new AliESDtrack(track);
   }
   if(esdtrack && esdtrack->PropagateToDCA(vtxESDSkip, fEvent->GetMagneticField(), kBeampiperadius, dcaD, covD)){
    // protection
    dcaxy = dcaD[0];
    if(covD[0]) dcansigmaxy = dcaxy/TMath::Sqrt(covD[0]);
    if(!covD[0]) dcansigmaxy = -49.;
   }
   delete esdtrack;
   delete vtxESDSkip;
  } 
  else{
   AliESDtrack *esdtrack = NULL;
   if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
    // Case ESD track: take copy constructor
    AliESDtrack *tmptrack = dynamic_cast<AliESDtrack *>(track);
    if(tmptrack) esdtrack = new AliESDtrack(*tmptrack);
   } else {
    // Case AOD track: take different constructor
    esdtrack = new AliESDtrack(track);
   }
   if(esdtrack) esdtrack->GetImpactParameters(dcaF, covF); 
   dcaxy = dcaF[0];
   if(covF[0]) dcansigmaxy = dcaxy/TMath::Sqrt(covF[0]);
   if(!covF[0]) dcansigmaxy = -49.;
   delete esdtrack;
  }
}


//______________________________________________________
void AliHFEextraCuts::GetHFEImpactParameterCuts(AliVTrack *track, Double_t &hfeimpactRcut, Double_t &hfeimpactnsigmaRcut){
	//
	// Get HFE impact parameter cut(pt dependent)
	//
  
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(!esdtrack) return;
    Double_t pt = esdtrack->Pt();	
    hfeimpactRcut = fIPcutParam[0]+fIPcutParam[1]*exp(fIPcutParam[2]*pt);  // abs R cut
    hfeimpactnsigmaRcut = fIPcutParam[3];                                  // sigma cut
  }
}
//______________________________________________________
void AliHFEextraCuts::GetMaxImpactParameterCutR(AliVTrack *track, Double_t &maximpactRcut){
	//
	// Get max impact parameter cut r (pt dependent)
	//
  
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if(!esdtrack) return;
    Double_t pt = esdtrack->Pt();	
    if(pt > 0.15) {
      maximpactRcut = 0.0182 + 0.035/TMath::Power(pt,1.01);  // abs R cut
    }
    else maximpactRcut = 9999999999.0;
  }
}
