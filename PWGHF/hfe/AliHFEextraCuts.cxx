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

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
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
#include "AliExternalTrackParam.h"

#include "AliHFEextraCuts.h"

ClassImp(AliHFEextraCuts)

const Int_t AliHFEextraCuts::fgkNQAhistos = 10;

//______________________________________________________
AliHFEextraCuts::AliHFEextraCuts(const Char_t *name, const Char_t *title):
  AliCFCutBase(name, title),
  fEvent(NULL),
  fCutCorrelation(0),
  fRequirements(0),
  fMinNClustersTPC(0),
  fMinNClustersTPCPID(0),
  fClusterRatioTPC(0.),
  fMinTrackletsTRD(0),
  fMaxChi2TRD(5.0),
  fMinNbITScls(0),
  fTRDtrackletsExact(0),
  fPixelITS(0),
  fDriftITS(0),
  fTPCclusterDef(0),
  fTPCclusterRatioDef(0),
  fFractionTPCShared(-1.0),
  fOppSideIPcut(kFALSE),
  fTOFsignalDx(1.0),
  fTOFsignalDz(1.0),
  fMagField(-10000),
  fCheck(kFALSE),
  fQAlist(0x0) ,
  fDebugLevel(0)
{
  //
  // Default Constructor
  //
  //printf("Set the number of min ITS clusters %d\n",(Int_t)fMinNbITScls);
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
  fMinNClustersTPCPID(c.fMinNClustersTPCPID),
  fClusterRatioTPC(c.fClusterRatioTPC),
  fMinTrackletsTRD(c.fMinTrackletsTRD),
  fMaxChi2TRD(c.fMaxChi2TRD),
  fMinNbITScls(c.fMinNbITScls),
  fTRDtrackletsExact(c.fTRDtrackletsExact),
  fPixelITS(c.fPixelITS),
  fDriftITS(c.fDriftITS),
  fTPCclusterDef(c.fTPCclusterDef),
  fTPCclusterRatioDef(c.fTPCclusterRatioDef),
  fFractionTPCShared(c.fFractionTPCShared),
  fOppSideIPcut(c.fOppSideIPcut),
  fTOFsignalDx(c.fTOFsignalDx),
  fTOFsignalDz(c.fTOFsignalDz),
  fMagField(c.fMagField),
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
    fMinNClustersTPCPID = c.fMinNClustersTPCPID;
    fMinTrackletsTRD = c.fMinTrackletsTRD;
    fMaxChi2TRD      = c.fMaxChi2TRD;
    fMinNbITScls = c.fMinNbITScls;
    fTRDtrackletsExact = c.fTRDtrackletsExact;
    fPixelITS = c.fPixelITS;
    fDriftITS = c.fDriftITS;
    fTPCclusterDef = c.fTPCclusterDef;
    fTPCclusterRatioDef = c.fTPCclusterRatioDef;
    fFractionTPCShared = c.fFractionTPCShared;
    fOppSideIPcut = c.fOppSideIPcut;
    fTOFsignalDx = c.fTOFsignalDx;
    fTOFsignalDz = c.fTOFsignalDz;
    fMagField = c.fMagField;
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
  TClass *type = o->IsA();
  AliDebug(2, Form("Object type %s", type->GetName()));
  if((type == AliESDtrack::Class()) || (type == AliAODTrack::Class())){
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
  AliDebug(1, Form("%s: Called", GetName()));
  ULong64_t survivedCut = 0;	// Bitmap for cuts which are passed by the track, later to be compared with fRequirements
  if(IsQAOn()) FillQAhistosRec(track, kBeforeCuts);
  // Apply cuts
  Float_t impactR = -999.;
  Float_t impactZ = -999.;
  Double_t hfeimpactR, hfeimpactnsigmaR;
  Double_t hfeimpactRcut, hfeimpactnsigmaRcut;
  Double_t maximpactRcut; 
  GetImpactParameters(track, impactR, impactZ);
  if(TESTBIT(fRequirements, kMinHFEImpactParamR) || TESTBIT(fRequirements, kMinHFEImpactParamNsigmaR)||TESTBIT(fRequirements, kMinHFEImpactParamRcharge)){
    // Protection for PbPb
    GetHFEImpactParameterCuts(track, hfeimpactRcut, hfeimpactnsigmaRcut);
    GetHFEImpactParameters(track, hfeimpactR, hfeimpactnsigmaR);
  }
  Int_t  nclsITS = GetITSNbOfcls(track);
  UInt_t nclsTPC = GetTPCncls(track);
  UInt_t nclsTPCPID = track->GetTPCsignalN();
  // printf("Check TPC findable clusters: %d, found Clusters: %d\n", track->GetTPCNclsF(), track->GetTPCNcls());
  Float_t fractionSharedClustersTPC = GetTPCsharedClustersRatio(track);
  Double_t ratioTPC = GetTPCclusterRatio(track);
  UChar_t trdTracklets;
  trdTracklets = track->GetTRDntrackletsPID();
  Float_t trdchi2=-999.;
  trdchi2=GetTRDchi(track);
  UChar_t itsPixel = track->GetITSClusterMap();
  Int_t status1 = GetITSstatus(track, 0);
  Int_t status2 = GetITSstatus(track, 1);
  Bool_t statusL0 = CheckITSstatus(status1);
  Bool_t statusL1 = CheckITSstatus(status2);
  Double_t tofsignalDx = 0.0;
  Double_t tofsignalDz = 0.0;
  GetTOFsignalDxDz(track,tofsignalDx,tofsignalDz);

  if(TESTBIT(fRequirements, kTPCfractionShared)) {
    // cut on max fraction of shared TPC clusters
    if(TMath::Abs(fractionSharedClustersTPC) < fFractionTPCShared) SETBIT(survivedCut, kTPCfractionShared);    
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
  if(TESTBIT(fRequirements, kMinHFEImpactParamRcharge)){
    // cut on min. HFE Impact Parameter in Radial direction, multiplied by particle charge
    Double_t charge = (Double_t)track->Charge();
    
    if(fMagField < 0)charge *= -1.;//the IP distribution side to be chosen depends on magnetic field polarization
    if(fOppSideIPcut)charge*=-1.;//in case we choose to select electrons from the side of the photon peak
    Double_t hfeimpactRcharge = hfeimpactR*charge;//switch selected side of the peak
    if(hfeimpactRcharge >= hfeimpactRcut) SETBIT(survivedCut, kMinHFEImpactParamRcharge);
  }
  if(TESTBIT(fRequirements, kMinHFEImpactParamNsigmaR)){
    // cut on max. HFE Impact Parameter n sigma in Radial direction
    // if(fAbsHFEImpactParamNsigmaR) {
    //   //if((TMath::Abs(hfeimpactnsigmaR) >= hfeimpactnsigmaRcut) && (TMath::Abs(hfeimpactnsigmaR) < 8)) { //mj debug
    //     if(TMath::Abs(hfeimpactnsigmaR) >= hfeimpactnsigmaRcut) {
    //       SETBIT(survivedCut, kMinHFEImpactParamNsigmaR);
    //       //  printf("0: hfeimpactnsigmaR= %lf hfeimpactnsigmaRcut= %lf\n",hfeimpactnsigmaR,hfeimpactnsigmaRcut);
    //     }
    //   }
    //   else {
    if(hfeimpactnsigmaRcut>0){
      if(hfeimpactnsigmaR >= hfeimpactnsigmaRcut) {
        SETBIT(survivedCut, kMinHFEImpactParamNsigmaR);
        //printf("1: hfeimpactnsigmaR= %lf hfeimpactnsigmaRcut= %lf\n",hfeimpactnsigmaR,hfeimpactnsigmaRcut);
      }
    }
    else{
      if(hfeimpactnsigmaR <= hfeimpactnsigmaRcut) {
        SETBIT(survivedCut, kMinHFEImpactParamNsigmaR);
        //printf("2: hfeimpactnsigmaR= %lf hfeimpactnsigmaRcut= %lf\n",hfeimpactnsigmaR,hfeimpactnsigmaRcut);
      }
    }
    //}
  }
  if(TESTBIT(fRequirements, kClusterRatioTPC)){
    // cut on min ratio of found TPC clusters vs findable TPC clusters
    if(ratioTPC >= fClusterRatioTPC) SETBIT(survivedCut, kClusterRatioTPC);
  }
  if(TESTBIT(fRequirements, kMinTrackletsTRD)){
    // cut on minimum number of TRD tracklets
    AliDebug(1, Form("%s: Min TRD cut: [%d|%d], Require exact number [%s]\n", GetName(), fMinTrackletsTRD, trdTracklets, fTRDtrackletsExact ? "Yes" : "No"));
    if(fTRDtrackletsExact){
      if(trdTracklets == fMinTrackletsTRD) {
        SETBIT(survivedCut, kMinTrackletsTRD);
        AliDebug(1, Form("%s: Track Selected", GetName()));
      }
    }else{
      if(trdTracklets >= fMinTrackletsTRD){ 
        SETBIT(survivedCut, kMinTrackletsTRD);
        AliDebug(1, Form("%s: Track Selected", GetName()));
      }
      //printf("Min number of tracklets %d\n",fMinTrackletsTRD);
    }
  }

  if(TESTBIT(fRequirements, kMaxTRDChi2)){
    // cut on TRD chi2
    AliDebug(1, Form("%s: Cut on TRD chi2: [%f|%f]\n", GetName(),fMaxChi2TRD, trdchi2));
    if(trdchi2 < fMaxChi2TRD) {
	SETBIT(survivedCut, kMaxTRDChi2);
        AliDebug(1,Form("%s: survived %f\n",GetName(),trdchi2));
    }
  }

  if(TESTBIT(fRequirements, kMinNbITScls)){
    // cut on minimum number of ITS clusters
    //printf(Form("Min ITS clusters: [%d|%d]\n", (Int_t)fMinNbITScls, nclsITS));
    AliDebug(1, Form("%s: Min ITS clusters: [%d|%d]\n", GetName(), fMinNbITScls, nclsITS));
    if(nclsITS >= ((Int_t)fMinNbITScls)) SETBIT(survivedCut, kMinNbITScls);
  }
  
  if(TESTBIT(fRequirements, kMinNClustersTPC)){
    // cut on minimum number of TPC tracklets
    //printf(Form("Min TPC cut: [%d|%d]\n", fMinNClustersTPC, nclsTPC));
    AliDebug(1, Form("%s: Min TPC cut: [%d|%d]\n", GetName(), fMinNClustersTPC, nclsTPC));
    if(nclsTPC >= fMinNClustersTPC) SETBIT(survivedCut, kMinNClustersTPC);
  }
  if(TESTBIT(fRequirements, kMinNClustersTPCPID)){
    AliDebug(1, Form("%s: Min TPC PID cut: [%d|%d]\n", GetName(), fMinNClustersTPCPID, nclsTPCPID));
    if(nclsTPCPID >= fMinNClustersTPCPID) SETBIT(survivedCut, kMinNClustersTPCPID);
  }
  if(TESTBIT(fRequirements, kDriftITS)){
    //printf("Require drift\n");
    switch(fDriftITS){
    case kFirstD:
      if(TESTBIT(itsPixel, 2)) SETBIT(survivedCut, kDriftITS);
      break;
    default: 
      SETBIT(survivedCut, kDriftITS);
      break;
  }
  }
  if(TESTBIT(fRequirements, kPixelITS)){
    // cut on ITS pixel layers
    AliDebug(1, Form("%s: ITS cluster Map: ", GetName()));
    //PrintBitMap(itsPixel);
    switch(fPixelITS){
      case kFirst:
        AliDebug(2, "First");
	//printf("First\n");
	      if(TESTBIT(itsPixel, 0)) 
	        SETBIT(survivedCut, kPixelITS);
        else
	        if(fCheck && !statusL0)
	          SETBIT(survivedCut, kPixelITS);
		    break;
      case kSecond: 
	//printf("Second\n");
        AliDebug(2, "Second");
	      if(TESTBIT(itsPixel, 1))
	        SETBIT(survivedCut, kPixelITS);
        else
	        if(fCheck && !statusL1)
	          SETBIT(survivedCut, kPixelITS);
		    break;
      case kBoth: 
	//printf("Both\n");
        AliDebug(2, "Both");
	      if(TESTBIT(itsPixel,0) && TESTBIT(itsPixel,1))
		      SETBIT(survivedCut, kPixelITS);
	      else
          if(fCheck && !(statusL0 && statusL1)) 
		        SETBIT(survivedCut, kPixelITS);
	      break;
      case kAny: 
	//printf("Any\n");
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
      case kNone:
        // No cut applied, set as survived
        SETBIT(survivedCut, kPixelITS);
        break;
      default: 
        // default, defined as no cut applied 
        AliDebug(2, "None");
        SETBIT(survivedCut, kPixelITS);
        break;
    }
    AliDebug(1, Form("%s: Survived Cut: %s\n", GetName(), TESTBIT(survivedCut, kPixelITS) ? "YES" : "NO"));
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
      Bool_t fBitsAboveThreshold=GetTPCCountSharedMapBitsAboveThreshold(track);
      if((fBitsAboveThreshold==0)&&(nclsTPCPID>80)) SETBIT(survivedCut, kTPCPIDCleanUp);
  }

  if(TESTBIT(fRequirements, kEMCALmatch)){
    if(track->GetStatus() && AliESDtrack::kEMCALmatch) SETBIT(survivedCut, kEMCALmatch);
  }

  if(TESTBIT(fRequirements, kRejectKinkDaughter)){
    //printf("test daughter\n");
    if(!IsKinkDaughter(track)) SETBIT(survivedCut, kRejectKinkDaughter);
    //else printf("Found kink daughter\n");
  }

  if(TESTBIT(fRequirements, kRejectKinkMother)){
    //printf("test mother\n");
    if(!IsKinkMother(track)) SETBIT(survivedCut, kRejectKinkMother);
    //else printf("Found kink mother\n");
  }
  if(TESTBIT(fRequirements, kTOFsignalDxy)){
    // cut on TOF matching cluster
    if((TMath::Abs(tofsignalDx) <= fTOFsignalDx) && (TMath::Abs(tofsignalDz) <= fTOFsignalDz)) SETBIT(survivedCut, kTOFsignalDxy);
  }
  if(TESTBIT(fRequirements, kITSpattern)){
    // cut on ITS pattern (every layer with a working ITS module must have an ITS cluster)
    if(CheckITSpattern(track)) SETBIT(survivedCut, kITSpattern); 
  }
  
  if(fRequirements == survivedCut){
    // 
    // Track selected
    //
    AliDebug(2, Form("%s: Track Survived cuts\n", GetName()));
    if(IsQAOn()) FillQAhistosRec(track, kAfterCuts);
    return kTRUE;
  }
  AliDebug(2, Form("%s: Track cut", GetName()));
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
  if(hStatusBits) {
    hStatusBits->Fill(0);  // Fill first bin with all tracks
    if(track->GetStatus() && AliESDtrack::kTOFpid) hStatusBits->Fill(1);
    if(!(track->GetStatus() && AliESDtrack::kTOFmismatch)) hStatusBits->Fill(2);
    if(track->GetStatus() && AliESDtrack::kEMCALmatch) hStatusBits->Fill(3);
    if(GetTPCCountSharedMapBitsAboveThreshold(track)==0) hStatusBits->Fill(4);
  }
  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(7 + when * fgkNQAhistos)))) htmp->Fill(track->GetTPCsignalN());

  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(8 + when * fgkNQAhistos)))) htmp->Fill(GetTRDchi(track));

  if((htmp = dynamic_cast<TH1F *>(fQAlist->At(9 + when * fgkNQAhistos)))) htmp->Fill(GetITSNbOfcls(track));
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
    qaList->AddAt((histo1D = new TH1F(Form("%s_tpcClr%s",GetName(),cutstr[icond].Data()), "Cluster Ratio TPC", 100, 0, 1.3)), 2 + icond * fgkNQAhistos);
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
    qaList->AddAt((histo1D = new TH1F(Form("%s_trdchi2perTracklet%s",GetName(),cutstr[icond].Data()), "chi2 per TRD tracklet", 100, 0, 10)), 8 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 8 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("Chi2 per TRD Tracklet");
    histo1D->GetYaxis()->SetTitle("Number of Tracks");
    qaList->AddAt((histo1D = new TH1F(Form("%s_ITScls%s",GetName(),cutstr[icond].Data()), "ITScls", 6, 0, 6)), 9 + icond * fgkNQAhistos);
    fQAlist->AddAt(histo1D, 9 + icond * fgkNQAhistos);
    histo1D->GetXaxis()->SetTitle("ITScls");
    histo1D->GetYaxis()->SetTitle("Number of ITS cls");

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
Int_t AliHFEextraCuts::GetITSstatus(const AliVTrack * const track, Int_t layer) const {
	//
	// Check ITS layer status
	//
	Int_t status = 0;
	if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
		Int_t det;
		Float_t xloc, zloc;
		const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
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
  const TBits *shared;
  if(track->IsA() == AliESDtrack::Class())
    shared = &((static_cast<AliESDtrack *>(track))->GetTPCSharedMap());
  else if(track->IsA() == AliAODTrack::Class())
    shared = &((static_cast<AliAODTrack *>(track))->GetTPCSharedMap());
   else 
    return 0;

	if(shared->CountBits() >= 2) nsharebit=1;
  else nsharebit=0;
 
  return nsharebit;
}

//______________________________________________________
UInt_t AliHFEextraCuts::GetTPCncls(AliVTrack *track){
  //
  // Get Number of findable clusters in the TPC
  //
  Int_t nClusters = 0; // in case no Information available consider all clusters findable
  TClass *type = track->IsA();
  if(TESTBIT(fTPCclusterDef, kFoundIter1)){
    if(type == AliESDtrack::Class()){
      AliDebug(2, ("Using def kFoundIter1"));
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      nClusters = esdtrack->GetTPCNclsIter1();
    } else {
      AliDebug(2, ("Number of clusters in iteration 1 not available for AOD tracks"));
    }
  } else if(TESTBIT(fTPCclusterDef, kCrossedRows)){
    AliDebug(2, ("Using def kCrossedRows"));
    nClusters = static_cast<UInt_t>(track->GetTPCClusterInfo(2,1));
  } else if(TESTBIT(fTPCclusterDef, kFound)){
    AliDebug(2, ("Using def kFound"));
    nClusters = track->GetTPCNcls();
  }
  else {
    AliDebug(2, ("Using def kFoundAll"));
    if(type == AliESDtrack::Class()){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      const TBits &clusterTPC = esdtrack->GetTPCClusterMap();
      nClusters = clusterTPC.CountBits();
    } 
    else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      const TBits &clusterTPC = aodtrack->GetTPCClusterMap();
      nClusters = clusterTPC.CountBits();
    }  
  }
  return nClusters;
}

//______________________________________________________
Float_t AliHFEextraCuts::GetTPCsharedClustersRatio(AliVTrack *track){
  //
  // Get fraction of shared TPC clusters
  //
  Float_t fracClustersTPCShared = 0.0; 
  Int_t nClustersTPC = track->GetTPCNcls();
  Int_t nClustersTPCShared = 0;
  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    nClustersTPCShared = esdtrack->GetTPCnclsS();
  } else if(type == AliAODTrack::Class()){
    AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
    const TBits &shared = aodtrack->GetTPCSharedMap();
    nClustersTPCShared = shared.CountBits(0) - shared.CountBits(159);
  }
  if (nClustersTPC!=0) {
    fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
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
  TClass *type = track->IsA();
  if(TESTBIT(fTPCclusterRatioDef, kCROverFindable)){
    AliDebug(2, "Using ratio def kCROverFindable");
    clusterRatio = track->GetTPCNclsF() ? track->GetTPCClusterInfo(2,1)/static_cast<Double_t>(track->GetTPCNclsF()) : 1.; // crossed rows/findable
  } else if(TESTBIT(fTPCclusterRatioDef, kFoundOverCR)){
    AliDebug(2, "Using ratio def kFoundOverCR");
    clusterRatio = track->GetTPCClusterInfo(2,0); // found/crossed rows
  } else if(TESTBIT(fTPCclusterRatioDef, kFoundOverFindableIter1)){
    if(type == AliESDtrack::Class()){
      AliDebug(2, "Using ratio def kFoundOverFindableIter1");
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      clusterRatio = esdtrack->GetTPCNclsFIter1() ? static_cast<Double_t>(esdtrack->GetTPCNclsIter1())/static_cast<Double_t>(esdtrack->GetTPCNclsFIter1()) : 1.; // crossed
    } else {
      AliDebug(2, "Cluster ratio after iteration 1 not possible for AOD tracks");
      clusterRatio = 1.;
    }
  } else if(TESTBIT(fTPCclusterRatioDef, kFoundOverFindable)) {
    AliDebug(2, "Using ratio def kFoundOverFindable");
    clusterRatio = track->GetTPCNclsF() ? static_cast<Double_t>(track->GetTPCNcls())/static_cast<Double_t>(track->GetTPCNclsF()) : 1.; // found/findable
  }
  else {
    AliDebug(2, "Using ratio def kFoundAllOverFindable");
    Int_t allclusters = 0;
    if(type == AliESDtrack::Class()){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      const TBits &clusterTPC = esdtrack->GetTPCClusterMap();
      allclusters = clusterTPC.CountBits();
    }
    else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      const TBits &clusterTPC = aodtrack->GetTPCClusterMap();
      allclusters = clusterTPC.CountBits();
    }
    clusterRatio = track->GetTPCNclsF() ? static_cast<Double_t>(allclusters)/static_cast<Double_t>(track->GetTPCNclsF()) : 1.; // foundall/findable
  }
  return clusterRatio;

}
//___________________________________________________
void AliHFEextraCuts::GetImpactParameters(AliVTrack *track, Float_t &radial, Float_t &z){
  //
  // Get impact parameter
  //

  const Double_t kBeampiperadius=3.;
  Double_t dcaD[2]={-999.,-999.},
           covD[3]={-999.,-999.,-999.};
  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    esdtrack->GetImpactParameters(radial, z);
  }
  else if(type == AliAODTrack::Class()){

    //case of AOD tracks
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(fEvent);
    if(!aodevent) {
      AliDebug(1, "No aod event available\n");
      return;
    }

    //Case ESD track: take copy constructor
    AliAODTrack *aodtrack = NULL;
    AliAODTrack *tmptrack = dynamic_cast<AliAODTrack *>(track);
    if(tmptrack) aodtrack = new AliAODTrack(*tmptrack);

    AliAODVertex *vtxAODSkip  = aodevent->GetPrimaryVertex();
    if(!vtxAODSkip) return;
    AliExternalTrackParam etp; etp.CopyFromVTrack(aodtrack);
    if(etp.PropagateToDCA(vtxAODSkip, aodevent->GetMagneticField(), kBeampiperadius, dcaD, covD)) {
      radial = dcaD[0];
      z = dcaD[1];
    }
    //if(vtxAODSkip) delete vtxAODSkip;
    if(aodtrack) delete aodtrack;
  }
}
//______________________________________________________
Bool_t AliHFEextraCuts::IsKinkDaughter(AliVTrack *track){
  //
  // Is kink Daughter
  //
  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    if(!esdtrack) return kFALSE;
    if(esdtrack->GetKinkIndex(0)>0) return kTRUE;
    else return kFALSE;

  }
  else if(type == AliAODTrack::Class()){
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if(aodtrack){
      AliAODVertex *aodvertex = aodtrack->GetProdVertex();
      if(!aodvertex) return kFALSE;
      if(aodvertex->GetType()==AliAODVertex::kKink) return kTRUE;
      else return kFALSE;
    }
    else return kFALSE;
  }

  return kFALSE;
}
//______________________________________________________
Bool_t AliHFEextraCuts::IsKinkMother(AliVTrack *track){
  //
  // Is kink Mother: only for ESD since need to loop over vertices for AOD
  //
  //

  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    if(!esdtrack) return kFALSE;
    if(esdtrack->GetKinkIndex(0)!=0) return kTRUE;
    else return kFALSE;
  }

  return kFALSE;

}

//______________________________________________________
Float_t AliHFEextraCuts::GetTRDchi(AliVTrack *track){
  //
  // Get TRDchi2
  //
  Int_t ntls(0);
  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      ntls = esdtrack->GetTRDntracklets();
      return ntls ? esdtrack->GetTRDchi2()/ntls : -999;
  }
  else if(type == AliAODTrack::Class()){
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if(aodtrack){
      return  999.;
    }
  }

  return 999.;

}

//______________________________________________________
Int_t AliHFEextraCuts::GetITSNbOfcls(AliVTrack *track){
  //
  // Get ITS nb of clusters
  //
  TClass *type = track->IsA();
  if(type == AliESDtrack::Class()){
    AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
    return esdtrack->GetITSclusters(0);

  }
  else if(type == AliAODTrack::Class()){
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if(aodtrack){
      return  aodtrack->GetITSNcls();
    }
  }

  return 0;
}
//______________________________________________________
void AliHFEextraCuts::GetHFEImpactParameters(const AliVTrack * const track, Double_t &dcaxy, Double_t &dcansigmaxy){
  //
  // Get HFE impact parameter (with recalculated primary vertex)
  //
  dcaxy=0;
  dcansigmaxy=0;
  if(!fEvent){
    AliDebug(1, "No Input event available\n");
    return;
  }
  TString type = track->IsA()->GetName();
  const Double_t kBeampiperadius=3.;
  Double_t dcaD[2]={-999.,-999.},
           covD[3]={-999.,-999.,-999.};
  Bool_t isRecalcVertex(kFALSE);

  if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
    //case of ESD tracks
    AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fEvent);
    if(!esdevent) {
      AliDebug(1, "No esd event available\n");
      return;
    }

    const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
    if(!vtxESDSkip) return;

    //case ESD track: take copy constructor
    const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
    if(tmptrack){

      if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
        vtxESDSkip = RemoveDaughtersFromPrimaryVtx(esdevent, track);
        isRecalcVertex = kTRUE;
      }

      if(vtxESDSkip){
        AliESDtrack esdtrack(*tmptrack);
        fMagField = fEvent->GetMagneticField();
        if(esdtrack.PropagateToDCA(vtxESDSkip, fMagField, kBeampiperadius, dcaD, covD)){
          dcaxy = dcaD[0];
          if(covD[0]) dcansigmaxy = dcaxy/TMath::Sqrt(covD[0]);
        }
        if(isRecalcVertex) delete vtxESDSkip;
      }
    }
  }
  else {
    //case of AOD tracks
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(fEvent);
    if(!aodevent) {
      AliDebug(1, "No aod event available\n");
      return;
    }

    AliAODVertex *vtxAODSkip  = aodevent->GetPrimaryVertex();
    if(!vtxAODSkip) return;

    //Case ESD track: take copy constructor
    const AliAODTrack *tmptrack = dynamic_cast<const AliAODTrack *>(track);
    if(tmptrack){ 

      if(vtxAODSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
        vtxAODSkip = RemoveDaughtersFromPrimaryVtx(aodevent, track);
        isRecalcVertex = kTRUE;
      } 
      if(vtxAODSkip){
        AliAODTrack aodtrack(*tmptrack);
        AliExternalTrackParam etp; etp.CopyFromVTrack(&aodtrack);
        fMagField = aodevent->GetMagneticField();
        if(etp.PropagateToDCA(vtxAODSkip,fMagField, kBeampiperadius, dcaD, covD)) {
          dcaxy = dcaD[0];
          if(covD[0]) dcansigmaxy = dcaxy/TMath::Sqrt(covD[0]);
        }
        if(isRecalcVertex) delete vtxAODSkip;
      }
    }
  }
}

//______________________________________________________
void AliHFEextraCuts::GetHFEImpactParameters(const AliVTrack * const track, Double_t dcaD[2], Double_t covD[3]){
	//
	// Get HFE impact parameter (with recalculated primary vertex)
	//
  if(!fEvent){
    AliDebug(1, "No Input event available\n");
    return;
  }
  const Double_t kBeampiperadius=3.;
  TString type = track->IsA()->GetName();
  Bool_t isRecalcVertex(kFALSE);
  
  if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
    //case of ESD tracks
    AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fEvent);
    if(!esdevent) {
      AliDebug(1, "No esd event available\n");
      return;
    }

    // Check whether primary vertex is available
    const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
    if(!vtxESDSkip) return;

    const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
    if(tmptrack){
      //case ESD track: take copy constructor

      if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
        vtxESDSkip = RemoveDaughtersFromPrimaryVtx(esdevent, track);
        isRecalcVertex = kTRUE;
      }
      if(vtxESDSkip){
        AliESDtrack esdtrack(*tmptrack);
        fMagField = fEvent->GetMagneticField();
        esdtrack.PropagateToDCA(vtxESDSkip, fMagField, kBeampiperadius, dcaD, covD);

        if(isRecalcVertex) delete vtxESDSkip;
      }
    }
  }
  else {
    //case of AOD tracks
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(fEvent);
    if(!aodevent) {
      AliDebug(1, "No aod event available\n");
      return;
    }

    // Check whether primary vertex is available
    AliAODVertex *vtxAODSkip  = aodevent->GetPrimaryVertex();
    if(!vtxAODSkip) return;

    //Case ESD track: take copy constructor
    const AliAODTrack *tmptrack = dynamic_cast<const AliAODTrack *>(track);
    if(tmptrack){

      if(vtxAODSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
        vtxAODSkip = RemoveDaughtersFromPrimaryVtx(aodevent, track);
        isRecalcVertex = kTRUE;
      }
      if(vtxAODSkip){
        AliAODTrack aodtrack(*tmptrack);
        AliExternalTrackParam etp; etp.CopyFromVTrack(&aodtrack);
        fMagField = aodevent->GetMagneticField();
        etp.PropagateToDCA(vtxAODSkip, fMagField, kBeampiperadius, dcaD, covD);

        if(isRecalcVertex) delete vtxAODSkip;
      }
    }
  }
}

//______________________________________________________
void AliHFEextraCuts::GetHFEImpactParameterCuts(const AliVTrack * const track, Double_t &hfeimpactRcut, Double_t &hfeimpactnsigmaRcut){
	//
	// Get HFE impact parameter cut(pt dependent)
	//
  
  Double_t pt = track->Pt();	
  hfeimpactRcut = fIPcutParam[0]+fIPcutParam[1]*exp(fIPcutParam[2]*pt);  // abs R cut
  hfeimpactnsigmaRcut = fIPcutParam[3];                                  // sigma cut
}
//______________________________________________________
void AliHFEextraCuts::GetMaxImpactParameterCutR(const AliVTrack * const track, Double_t &maximpactRcut){
	//
	// Get max impact parameter cut r (pt dependent)
	//
  
  Double_t pt = track->Pt();	
  if(pt > 0.15) {
    maximpactRcut = 0.0182 + 0.035/TMath::Power(pt,1.01);  // abs R cut
  }
  else maximpactRcut = 9999999999.0;
}
//______________________________________________________
void AliHFEextraCuts::GetTOFsignalDxDz(const AliVTrack * const track, Double_t &tofsignalDx, Double_t &tofsignalDz){
  //
  // TOF matching 
  //
  
  TString type = track->IsA()->GetName();
  if(!type.CompareTo("AliESDtrack")){
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(!esdtrack) return;
    tofsignalDx = esdtrack->GetTOFsignalDx();
    tofsignalDz = esdtrack->GetTOFsignalDz();
  }

}

//______________________________________________________
Bool_t AliHFEextraCuts::CheckITSpattern(const AliVTrack *const track) const {
  //
  // Check if every ITS layer, which has a module which is alive, also
  // has an ITS cluster
  //
  Bool_t patternOK(kTRUE);
  Int_t status(0);
  for(Int_t ily = 0; ily < 6; ily++){
    status = GetITSstatus(track, ily);
    if(CheckITSstatus(status)){
      // pixel alive, check whether layer has a cluster
      if(!TESTBIT(track->GetITSClusterMap(),ily)){
        // No cluster even though pixel is alive - reject track
        patternOK = kFALSE;
        break;
      }
    }
  }
  return patternOK;
}

//---------------------------------------------------------------------------
const AliVVertex* AliHFEextraCuts::RemoveDaughtersFromPrimaryVtx(const AliESDEvent * const esdevent, const AliVTrack * const track) {
  //
  // This method returns a primary vertex without the daughter tracks of the 
  // candidate and it recalculates the impact parameters and errors for ESD tracks.
  // 
  // The output vertex is created with "new". 
  //

  //const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();

  AliVertexerTracks vertexer(fEvent->GetMagneticField());
  vertexer.SetITSMode();
  vertexer.SetMinClusters(4);
  Int_t skipped[2];
  skipped[0] = track->GetID();
  vertexer.SetSkipTracks(1,skipped);

  //diamond constraint
  vertexer.SetConstraintOn();
  Float_t diamondcovxy[3];
  esdevent->GetDiamondCovXY(diamondcovxy);
  Double_t pos[3]={esdevent->GetDiamondX(),esdevent->GetDiamondY(),0.};
  Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
  AliESDVertex diamond(pos,cov,1.,1);
  vertexer.SetVtxStart(&diamond);

  const AliVVertex *vtxESDSkip = vertexer.FindPrimaryVertex(fEvent);

  return vtxESDSkip;
}

//---------------------------------------------------------------------------
AliAODVertex* AliHFEextraCuts::RemoveDaughtersFromPrimaryVtx(const AliAODEvent * const aod, const AliVTrack * const track) {
  //
  // This method returns a primary vertex without the daughter tracks of the 
  // candidate and it recalculates the impact parameters and errors for AOD tracks.
  // The output vertex is created with "new". 
  //

  AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
  if(!vtxAOD) return 0;
  TString title=vtxAOD->GetTitle();
  if(!title.Contains("VertexerTracks")) return 0;

  AliVertexerTracks vertexer(aod->GetMagneticField());

  vertexer.SetITSMode();
  vertexer.SetMinClusters(3);
  vertexer.SetConstraintOff();

  if(title.Contains("WithConstraint")) {
    Float_t diamondcovxy[3];
    aod->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3]={aod->GetDiamondX(),aod->GetDiamondY(),0.};
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
    AliESDVertex diamond(pos,cov,1.,1);
    vertexer.SetVtxStart(&diamond);
  }
  Int_t skipped[2]; for(Int_t i=0;i<2;i++) skipped[i]=-1;
  Int_t id = (Int_t)track->GetID();
  if(!(id<0)) skipped[0] = id;

  /*// exclude tracks with global constrained parameters
  Int_t nTracks=aod->GetNumberOfTracks();
  for(Int_t i=0; i<nTracks; i++){
    t = aod->GetTrack(i);
    if(t->TestFilterMask(512)){
      id = (Int_t)t->GetID();
      skipped[nTrksToSkip++] = id;
    }
  }*/

  vertexer.SetSkipTracks(1,skipped);
  AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(aod);

  if(!vtxESDNew) return 0;
  if(vtxESDNew->GetNContributors()<=0) {
    delete vtxESDNew; vtxESDNew=NULL;
    return 0;
  }

  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vtxESDNew->GetXYZ(pos); // position
  vtxESDNew->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vtxESDNew->GetChi2toNDF();
  delete vtxESDNew; vtxESDNew=NULL;

  AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);

  return vtxAODNew;
}
