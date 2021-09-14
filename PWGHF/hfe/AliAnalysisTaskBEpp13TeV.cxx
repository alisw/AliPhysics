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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliHFEextraCuts.h"
#include "AliHFEtools.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskBEpp13TeV.h"

class AliAnalysisTaskBEpp13TeV;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskBEpp13TeV) // classimp: necessary for root

AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV()
:AliAnalysisTaskSE()
,fAOD(0)
,fOutputList(0)
,fPIDResponse(0)
,fExtraCuts(0)

,fIsMC(false)
,fMinTPCnCrossedRow(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(3)
,fITSlayer(2)
,fTPCnsigmaLow(-1)
,fTPCnsigmaHigh(3)
,fTOFnsigma(3)

,hNrEvents(0)
,fHistPt(0)

,hFilterMask(0)
,hTPCnCrossedRow(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,dcaTrack(0)
,dcaPion(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV(const char* name)
:AliAnalysisTaskSE(name)
,fAOD(0)
,fOutputList(0)
,fPIDResponse(0)
,fExtraCuts(0)

,fIsMC(false)
,fMinTPCnCrossedRow(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(3)
,fITSlayer(2)
,fTPCnsigmaLow(-1)
,fTPCnsigmaHigh(3)
,fTOFnsigma(3)

,hNrEvents(0)
,fHistPt(0)

,hFilterMask(0)
,hTPCnCrossedRow(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,dcaTrack(0)
,dcaPion(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::~AliAnalysisTaskBEpp13TeV()
{
  // destructor
  if(fOutputList) delete fOutputList;
  if(fExtraCuts) delete fExtraCuts;
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserCreateOutputObjects()
{
  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  int nPtBins = 11;
  double ptbinningX[12] = { 1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10. };

  int nBinsPID = 400;
  double minPID = -10;
  double maxPID = 10;
  double binLimPID[nBinsPID+1];
  for(int i=0; i<=nBinsPID; i++) binLimPID[i] = minPID + (maxPID-minPID)/nBinsPID*(double)i;

  int nBinsIP = 4000;
  double minIP = -0.2;
  double maxIP = 0.2;
  double binLimIP[nBinsIP+1];
  for(int i=0; i<=nBinsIP; i++) binLimIP[i] = minIP + (maxIP-minIP)/nBinsIP*(double)i;
	
  // example of a histogram
  fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);
  fOutputList->Add(fHistPt);
  hNrEvents = new TH1F("hNrEvents","number of events",1,0.,1.);
  fOutputList->Add(hNrEvents);
  hFilterMask = new TH1F("hFilterMask", "", 2, 0., 2.);
  fOutputList->Add(hFilterMask);
  hTPCnCrossedRow = new TH1F("hTPCnCrossedRow", "", 200, 0., 200.);
  fOutputList->Add(hTPCnCrossedRow);
  hTPCclsPID = new TH1F("hTPCclsPID", "", 200, 0., 200.);
  fOutputList->Add(hTPCclsPID);
  hTPCchi2 = new TH1F("hTPCchi2", "", 100, 0., 10.);
  fOutputList->Add(hTPCchi2);
  hTPCclsRatio = new TH1F("hTPCclsRatio", "", 15, 0., 1.5);
  fOutputList->Add(hTPCclsRatio);
  hITSNcls= new TH1F("hITSNcls", "", 10, 0., 10.);
  fOutputList->Add(hITSNcls);
  hITSlayer = new TH1F("hITSlayer", "", 3, 0.5, 3.5);
  fOutputList->Add(hITSlayer);
  hDCAxy = new TH1F("hDCAxy", "", 600, -3., 3.);
  fOutputList->Add(hDCAxy);
  hDCAz = new TH1F("hDCAz", "", 600, -3., 3.);
  fOutputList->Add(hDCAz);
  hPt = new TH1F("hPt", "pt; (GeV/c)", 300, 0., 30.);
  fOutputList->Add(hPt);
  hEta = new TH1F("hEta", "", 200, -1., 1.);
  fOutputList->Add(hEta);
  hPhi = new TH1F("hPhi", "", 700, -0.5, 6.5);
  fOutputList->Add(hPhi);
  hTPCnsigma = new TH2F("hTPCnsigma", "n#sigma_{TPC} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTPCnsigma);
  hTPCnsigmaTOFcut = new TH2F("hTPCnsigmaTOFcut", "n#sigma_{TPC |TOFcut} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTPCnsigmaTOFcut);
  hTPCnsigmaTOFcutPt = new TH2F("hTPCnsigmaTOFcutPt", "n#sigma_{TPC |TOFcut} vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hTPCnsigmaTOFcutPt);
  hTPCnsigmaQA = new TH2F("hTPCnsigmaQA", "n#sigma_{TPC} QA vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTPCnsigmaQA);
  hTPCnsigmaPiQA = new TH2F("hTPCnsigmaPiQA", "n#sigma_{TPC |TOFcut} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTPCnsigmaPiQA);
  hTOFnsigma = new TH2F("hTOFnsigma", "n#sigma_{TOF} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TOF}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTOFnsigma);
  hTOFnsigmaQA = new TH2F("hTOFnsigmaQA", "n#sigma_{TOF} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TOF}", 500, 0., 10., 400, -10., 10.);
  fOutputList->Add(hTOFnsigmaQA);
  dcaTrack = new TH2F("dcaTrack", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaTrack);
  dcaPion = new TH2F("dcaPion", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaPion);


    
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserExec(Option_t *)
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  
  if(!fExtraCuts)
	fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  fExtraCuts->SetRecEventInfo(fAOD);

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
	AliDebug(1,"Using default PID Response");
	fPIDResponse = AliHFEtools::GetDefaultPID(false, fInputEvent->IsA()==AliAODEvent::Class());
  }

  // Event selection
  if(!PassEventCuts(fAOD)) return;
  if(PassPileUpEvent(fAOD)) return;

  //Look for kink mother
  double *fListOfMotherKink = 0;
  int fNumberOfVertices = 0;
  int fNumberOfMotherKink = 0;

  fNumberOfVertices = fAOD->GetNumberOfVertices();
  fListOfMotherKink = new double[fNumberOfVertices];

  for(int iVertex=0; iVertex<fNumberOfVertices; iVertex++){
	AliAODVertex *aodvtx = fAOD->GetVertex(iVertex);
	if(!aodvtx) continue;
	if(aodvtx->GetType()==AliAODVertex::kKink){
	  AliAODTrack *mother = (AliAODTrack*)aodvtx->GetParent();
	  if(!mother) continue;
	  int idmother = mother->GetID();
	  fListOfMotherKink[fNumberOfMotherKink] = idmother;
	  fNumberOfMotherKink++;
	}
  }

  double fBz=-999.;
  if(fAOD->GetMagneticField()<0) fBz = -1.;
  else if(fAOD->GetMagneticField()>0) fBz = 1.;
  else return;
  
  hNrEvents->Fill(0);
  
  for(int iTracks = 0; iTracks<fAOD->GetNumberOfTracks(); iTracks++){
    AliAODTrack *aodTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
    if(!aodTrack) continue;

	// Track selection
	if(!PassTrackCuts(aodTrack)) continue;

	// Reject kink
	bool kinkmotherpass = true;
	for(int kinkmother=0; kinkmother<fNumberOfMotherKink; kinkmother++){
	  if(aodTrack->GetID()==fListOfMotherKink[kinkmother]){
		kinkmotherpass=false;
		continue;
	  }
	}
	if(!kinkmotherpass) continue;
    
    fHistPt->Fill(aodTrack->Pt());
	
	double pt = aodTrack->Pt();
	double hfeImpactParam = -999., hfeImpactParamResol = -999.;
	fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);
	
	// electron identification
	double fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
	double fTOFnSigma = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);
	
	hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
	hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);
	if(TMath::Abs(fTOFnSigma)>fTOFnsigma) continue;

	hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
	hTPCnsigmaTOFcutPt->Fill(pt, fTPCnSigma);

	if(fTPCnSigma>-5 && fTPCnSigma<-3){
	  hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
	  dcaPion->Fill(pt, hfeImpactParam*fBz*aodTrack->Charge());
	}

	if(fTPCnSigma<fTPCnsigmaLow || fTPCnSigma>fTPCnsigmaHigh) continue;
		
	hTPCnsigmaQA->Fill(aodTrack->P(), fTPCnSigma);
	hTOFnsigmaQA->Fill(aodTrack->P(), fTOFnSigma);

	dcaTrack->Fill(pt, hfeImpactParam*fBz*aodTrack->Charge());

  }
  
  PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                     // the output manager which will take care of writing
                                                        // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//__________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassEventCuts(AliAODEvent *event){
  
  //event selection cuts
  AliAODVertex *vtx = event->GetPrimaryVertex();
  if(vtx->GetNContributors()<2) return false;

  // need to confirm -----
  AliAODVertex *vtxSPD = event->GetPrimaryVertexSPD();
  double cov[6]={0};
  vtxSPD->GetCovarianceMatrix(cov);
  double zRes = TMath::Sqrt(cov[5]);
  if(vtxSPD->IsFromVertexerZ() && (zRes>0.25)) return false;
  //---------------------
  
  if(TMath::Abs(vtx->GetZ()) > 10) return false;
  
  // need to confirm -----
  double zvtx = vtx->GetZ();
  double zvtxSPD = vtxSPD->GetZ();
  if(TMath::Abs(zvtx-zvtxSPD)>0.5) return false;
  //----------------------

  return true;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassPileUpEvent(AliAODEvent *event){
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5); //Multi Vertex pileup selection
  utils.SetMaxPlpChi2MV(5); // max value of Chi2perNDF of the pileup and multi-vertex
  utils.SetMinWDistMV(15); // min of the sqrt of weighted distance between the primary and the pileup vertex, multi-vertex
  utils.SetCheckPlpFromDifferentBCMV(false);
  bool isPileupFromMV = utils.IsPileUpMV(event);
  return isPileupFromMV;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassTrackCuts(AliAODTrack *track){

	if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;
	if(TMath::Abs(track->Eta()) >= 0.8) return false;
	if(track->Pt()<0.5 || track->Pt()>15.) return false;

	// basic tracking
	ULong_t status = track->GetStatus();
	if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;
	
	// TPC cut
	unsigned short findableTPC = track->GetTPCNclsF();
	unsigned short TPCsignalN = track->GetTPCsignalN();
	double FoundOverFindable = (findableTPC ? static_cast<float>(TPCsignalN)/static_cast<float>(findableTPC) : 0);
	if(track->GetTPCNCrossedRows() < fMinTPCnCrossedRow || track->GetTPCsignalN() < fMinTPCNclsPID || track->Chi2perNDF() > fMaxTPCchi2 || FoundOverFindable < fMinTPCclsRatio) return false;

	// ITS cut
	//std::cout<<fMinITSNcls<<" : "<<track->GetITSNcls()<<std::endl;
	if(track->GetITSNcls() < fMinITSNcls) return false;
	//std::cout<<"kFirst: "<<AliHFEextraCuts::kFirst<<", kAny: "<<AliHFEextraCuts::kAny<<", kBoth: "<<AliHFEextraCuts::kBoth<<std::endl;
	//std::cout<<"fITSlayer: "<<fITSlayer<<std::endl;
	if(fITSlayer==AliHFEextraCuts::kFirst){
		if(!(track->HasPointOnITSLayer(0))) return false;
		hITSlayer->Fill(1);
	}else if(fITSlayer==AliHFEextraCuts::kAny){
		if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(2);
	}else if(fITSlayer==AliHFEextraCuts::kBoth){
		if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(3);
	}else return false;
	
	// dca cut
	float dcaxy = -999.; float dcaz = -999.;
	fExtraCuts->GetImpactParameters((AliVTrack *)track, dcaxy, dcaz);
	//std::cout<<dcaxy<<" : "<<dcaz<<std::endl;
	if(TMath::Abs(dcaxy)>1. || TMath::Abs(dcaz)>2.) return false;
	
	hFilterMask->Fill(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA));
	hEta->Fill(track->Eta());
	hPhi->Fill(track->Phi());
	hPt->Fill(track->Pt());
	hDCAxy->Fill(dcaxy);
	hDCAz->Fill(dcaz);
	hTPCnCrossedRow->Fill(track->GetTPCNCrossedRows());
	hTPCclsPID->Fill(track->GetTPCsignalN());
	hTPCchi2->Fill(track->Chi2perNDF());
	hTPCclsRatio->Fill(FoundOverFindable);
	hITSNcls->Fill(track->GetITSNcls());
	return true;
}

