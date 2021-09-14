/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Jonghan Park (jonghan@cern.ch)                                 *
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
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
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
//,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fAODMCParticle(0)

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

,hBhadronPt(0)
,hD0Pt(0)
,hLcPt(0)

,hGenBePt(0)
,hRecBePt_track(0)
,hRecBePt_tof(0)
,hRecBePt_tpc(0)

,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,dcaTrack(0)
,dcaPion(0)
,dcaBeauty(0)
,dcaCharm(0)
,dcaDalitz(0)
,dcaConv(0)
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
//,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fAODMCParticle(0)

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

,hBhadronPt(0)
,hD0Pt(0)
,hLcPt(0)

,hGenBePt(0)
,hRecBePt_track(0)
,hRecBePt_tof(0)
,hRecBePt_tpc(0)

,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,dcaTrack(0)
,dcaPion(0)
,dcaBeauty(0)
,dcaCharm(0)
,dcaDalitz(0)
,dcaConv(0)
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
  double ptbinningD0[13] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36. };
  double ptbinningLc[7] = { 1., 2., 4., 6., 8., 12., 24. };

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
  hBhadronPt = new TH1F("hBhadronPt", "", 100, 0, 100);
  fOutputList->Add(hBhadronPt);
  hD0Pt = new TH1F("hD0Pt", "", 12, ptbinningD0);
  fOutputList->Add(hD0Pt);
  hLcPt = new TH1F("hLcPt", "", 6, ptbinningLc);
  fOutputList->Add(hLcPt);
  hGenBePt = new TH1F("hGenBePt", "", nPtBins, ptbinningX);
  fOutputList->Add(hGenBePt);
  hRecBePt_track = new TH1F("hRecBePt_track", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_track);
  hRecBePt_tof = new TH1F("hRecBePt_tof", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_tof);
  hRecBePt_tpc = new TH1F("hRecBePt_tpc", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_tpc);
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
  dcaBeauty = new TH2F("dcaBeauty", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBeauty);
  dcaCharm = new TH2F("dcaCharm", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaCharm);
  dcaDalitz = new TH2F("dcaDalitz", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDalitz);
  dcaConv = new TH2F("dcaConv", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaConv);
    
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


  if(fIsMC){
		
	fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if(!fAODArrayMCInfo){
	  AliError("No AOD MC particles");
	  return;
	}
	
	for(int iMC = 0; iMC<fAODArrayMCInfo->GetEntries(); iMC++){
	
	  fAODMCParticle = (AliAODMCParticle*) fAODArrayMCInfo->At(iMC);
	  
	  int hf = -999;
	  double hfpt = -999., hfeta = -999.;
	  hf = GetHeavyFlavours(fAODMCParticle, hfpt, hfeta);
	  
	  if(TMath::Abs(hfeta)<0.5){
		if(hf==kPromptD0) hD0Pt->Fill(hfpt);
		if(hf==kPromptLc) hLcPt->Fill(hfpt);
	  }
	  if(TMath::Abs(hfeta<0.8)){
		if(hf==kPromptB || hf==kNonPromptD) hBhadronPt->Fill(hfpt);
	  }

	  int src = -999, srcPdg = -999;
	  double srcPt = -999.;
	  src = GetElecSource(fAODMCParticle, srcPt, srcPdg);

	  if(TMath::Abs(fAODMCParticle->Eta()) < 0.8){
		if(src==kDirectBeauty || src==kBeautyCharm) hGenBePt->Fill(fAODMCParticle->Pt());
	  }
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
    
	double pt = aodTrack->Pt();
	double hfeImpactParam = -999., hfeImpactParamResol = -999.;
	fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);
	double IP = hfeImpactParam*fBz*aodTrack->Charge();
  
	int mcelectronSource=-999, mcelectronSourcePDG=-999;
	double mcelectronSourcePt=-999.;
	if(fIsMC){
	  fAODMCParticle = NULL;
			
	  int label = TMath::Abs(aodTrack->GetLabel());
	  if(label < fAODArrayMCInfo->GetEntriesFast())
		fAODMCParticle = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
	  if(fAODMCParticle){
		AliDebug(2, "Associated MC particle found");
		mcelectronSource = GetElecSource(fAODMCParticle, mcelectronSourcePt, mcelectronSourcePDG);
	  }
	  
	  // Fill beauty dca information
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
		hRecBePt_track->Fill(pt);
		dcaBeauty->Fill(pt, IP);
	  }
	  // Fill charm dca information
	  if(mcelectronSource==kDirectCharm){
		dcaCharm->Fill(pt, IP);
	  }
	  // Fill Dalitz dca information
	  if(mcelectronSource>=5 && mcelectronSource<=15){
		dcaDalitz->Fill(pt, IP);
	  }
	  // Fill conversion dca information
	  if(mcelectronSource>=18 && mcelectronSource<=28){
		dcaConv->Fill(pt, IP);
	  }
	}

	// electron identification
	double fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
	double fTOFnSigma = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);
	
	hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
	hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);
	if(TMath::Abs(fTOFnSigma)>fTOFnsigma) continue;
	if(fIsMC){
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBePt_tof->Fill(pt);
	}


	hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
	hTPCnsigmaTOFcutPt->Fill(pt, fTPCnSigma);

	if(fTPCnSigma>-5 && fTPCnSigma<-3){
	  hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
	  dcaPion->Fill(pt, hfeImpactParam*fBz*aodTrack->Charge());
	}

	if(fTPCnSigma<fTPCnsigmaLow || fTPCnSigma>fTPCnsigmaHigh) continue;
	if(fIsMC){
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBePt_tpc->Fill(pt);
	}
		
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
  if(!vtx || vtx->GetNContributors()<2) return false;

  // need to confirm -----
  //AliAODVertex *vtxSPD = event->GetPrimaryVertexSPD();
  //double cov[6]={0};
  //vtxSPD->GetCovarianceMatrix(cov);
  //double zRes = TMath::Sqrt(cov[5]);
  //if(vtxSPD->IsFromVertexerZ() && (zRes>0.25)) return false;
  //---------------------
  
  if(TMath::Abs(vtx->GetZ()) > 10) return false;
  
  // need to confirm -----
  //double zvtx = vtx->GetZ();
  //double zvtxSPD = vtxSPD->GetZ();
  //if(TMath::Abs(zvtx-zvtxSPD)>0.5) return false;
  //----------------------

  return true;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassPileUpEvent(AliAODEvent *event){
  //This function checks if there was a pile up reconstructed with SPD
  bool isPileupfromSPDmulbins = event->IsPileupFromSPDInMultBins();
  if(isPileupfromSPDmulbins) return true;
  
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
	if(track->GetITSNcls() < fMinITSNcls) return false;
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
//_________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg){

  if(!mcpart) return kMisID;
  if(!fAODArrayMCInfo) return -1;
	
  if(TMath::Abs(mcpart->GetPdgCode()) != 11 ) return kElse;

  int origin = -1;
  bool isFinalOpenCharm = kFALSE;

  int iLabel = mcpart->GetMother();
  if((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
	AliDebug(1, "label is out of range, return\n");
	return -1;
  }
	
  AliAODMCParticle *mctrack = NULL; // will change all the time
  int tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1;
  AliAODMCParticle *partMother = mctrack;	//mtrack 
  AliAODMCParticle *partMotherCopy = mctrack;	//mtrack
  int maPdgcode = partMother->GetPdgCode();	//mpdg
  mpt = partMother->Pt();	//mpt
  mpdg = partMother->GetPdgCode();	//mpdg
  int gmaPdgcode, ggmaPdgcode;
  double gmpt, ggmpt;
  int gmpdg, ggmpdg;

  // if the mother is charmed hadron
  if((int(TMath::Abs(maPdgcode)/100.)%10)==4 || (int(TMath::Abs(maPdgcode)/1000.)%10)==4){
	if(TMath::Abs(maPdgcode)==411 || TMath::Abs(maPdgcode)==421 || TMath::Abs(maPdgcode)==431 || TMath::Abs(maPdgcode)==4122){
	  mpt = partMother->Pt();
	  mpdg = partMother->GetPdgCode();
	  isFinalOpenCharm = kTRUE;
	}
	if(!isFinalOpenCharm){
	  return -1;
	}
		
	// iterate until find B hadron as a  mother
	for(int i=1; i<100; i++){
	  int jLabel = partMother->GetMother();
	  if(jLabel == -1){
		return kDirectCharm;
	  }
	  if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
		AliDebug(1, "Stack label is negative, return\n");
		return -1;
	  }
			
	  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))){
		return -1;
	  }
	  int grandMaPDG = mctrack->GetPdgCode();
	  if(TMath::Abs(grandMaPDG)==511 || TMath::Abs(grandMaPDG)==521 || TMath::Abs(grandMaPDG)==531 || TMath::Abs(grandMaPDG)==5122){
		mpt = mctrack->Pt();
		mpdg = mctrack->GetPdgCode();
		return kBeautyCharm;
	  }
	  partMother = mctrack;
	} // end of iteration 
  }
  
  // if the mother is beauty hadron
  else if((int(TMath::Abs(maPdgcode)/100.)%10)==5 || (int(TMath::Abs(maPdgcode)/1000.)%10)==5){
	if(TMath::Abs(maPdgcode)==511 || TMath::Abs(maPdgcode)==521 || TMath::Abs(maPdgcode)==531 || TMath::Abs(maPdgcode)==5122){
	  mpt = partMotherCopy->Pt();
	  mpdg = partMotherCopy->GetPdgCode();
	  return kDirectBeauty;
	}
  }
	
  // if the mother is gamma
  else if(TMath::Abs(maPdgcode)==22){
	tmpMomLabel = partMotherCopy->GetMother();  // mother of photon
	mpt = partMotherCopy->Pt(); // pT of photon
	mpdg = partMotherCopy->GetPdgCode();
	if(tmpMomLabel==-1) return kGamma;  // no grandmother
	if((tmpMomLabel<0) || (tmpMomLabel>=fAODArrayMCInfo->GetEntriesFast())) {
	  return -1;
	}
	if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
	  return -1;
	}
	partMother = mctrack; // gmtrack
	partMotherCopy = mctrack; // gmtrack
	mpt = partMother->Pt(); // grand mother pT
	mpdg = partMother->GetPdgCode(); // grand mother PDG
	maPdgcode = partMother->GetPdgCode(); // grand mother PDG
		
	// check if the ligth meson is the decay product of heavy mesons
	tmpMomLabel = partMother->GetMother(); // grand grand mother of photon
	if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){//grand grand mother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
		partMother = mctrack; //ggmtrack
        gmaPdgcode = partMother->GetPdgCode(); //grand grand mother PDG
		mpt = partMother->Pt(); // grand grand mother pT
		mpdg = partMother->GetPdgCode(); // grand grand mother pT
		gmpt = partMother->Pt(); // grand grand mother pt
		gmpdg = partMother->GetPdgCode(); // grand grand mother pT

		if(TMath::Abs(maPdgcode)==111){
		  mpt = gmpt;
		  mpdg = gmpdg;
		  if(gmaPdgcode == 310) return kGammaK0s2P;
		  else if(gmaPdgcode == 130) return kGammaK0l2P;
		  else if(TMath::Abs(gmaPdgcode) == 321) return kGammaK2P;
		  else if(TMath::Abs(gmaPdgcode) == 3122) return kGammaLamda2P;
		  else if(gmaPdgcode == 3222) return kGammaSigma2P;
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaPi0;
		}
		else if(TMath::Abs(maPdgcode)==221){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaEta;
		}
		else if(TMath::Abs(maPdgcode)==223){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaOmega;
		}
		else if(TMath::Abs(maPdgcode)==333){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaPhi;
		}
		else if(TMath::Abs(maPdgcode)==331){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaEtaPrime;
		}
		else if(TMath::Abs(maPdgcode)==113){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaRho0;
		}
		else origin = kElse;//grand grand mother but nothing we identify
	  }//mctrack grandgrandmother
	}
	else{
	  // grandmother is primary
	  if(TMath::Abs(maPdgcode)==111){
		return kGammaPi0;
	  }
	  else if(TMath::Abs(maPdgcode)==221){
		return kGammaEta;
	  }
	  else if(TMath::Abs(maPdgcode)==223){
		return kGammaOmega;
	  }
	  else if(TMath::Abs(maPdgcode)==333){
		return kGammaPhi;
	  }
	  else if(TMath::Abs(maPdgcode)==331){
		return kGammaEtaPrime;
	  }
	  else if(TMath::Abs(maPdgcode)==113){
		return kGammaRho0;
	  }
	  else origin = kElse;//grandmother is primary but nothing we identify
	}
	return origin;
  }

  // if the mother is light meson
  else{
	
	tmpMomLabel = partMotherCopy->GetMother(); // grand mother
	mpt = partMotherCopy->Pt(); // mother pT
	mpdg = partMotherCopy->GetPdgCode(); // mother PDG
	if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){// grand mother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
		partMother = mctrack; // grand mother
		gmaPdgcode = partMother->GetPdgCode(); // grand mother PDG
		mpt = partMother->Pt(); // grand mother pT
		mpdg = partMother->GetPdgCode(); // grand mother PDG
		gmpt = partMother->Pt(); // grand mother pT
		gmpdg = partMother->GetPdgCode(); // grand mother PDG

		if(TMath::Abs(maPdgcode)==111){
		  mpt = gmpt;
		  mpdg = gmpdg;
		  if(gmaPdgcode == 310) return kK0s2P;
		  else if(gmaPdgcode == 130) return kK0l2P;
		  else if(TMath::Abs(gmaPdgcode) == 321) return kK2P;
		  else if(TMath::Abs(gmaPdgcode) == 3122) return kLamda2P;
		  else if(gmaPdgcode == 3222) return kSigma2P;
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kPi0;
		}
		else if(TMath::Abs(maPdgcode)==221){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kEta;
		}
		else if(TMath::Abs(maPdgcode)==223){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kOmega;
		}
		else if(TMath::Abs(maPdgcode)==333){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kPhi;
		}
		else if(TMath::Abs(maPdgcode)==331){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kEtaPrime;
		}
		else if(TMath::Abs(maPdgcode)==113){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kRho0;
		}
		else if(TMath::Abs(maPdgcode)==321){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kKe3;
		}
		else if(TMath::Abs(maPdgcode)==130){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kK0L;
		}
		else origin = kElse;//grandmother but nothing we identidy
	  }//mctrack grandmother
	}
	else {
	  // no grandmother
	  if(TMath::Abs(maPdgcode)==111) return kPi0;
	  else if(TMath::Abs(maPdgcode)==221) return kEta;
	  else if(TMath::Abs(maPdgcode)==223) return kOmega;
	  else if(TMath::Abs(maPdgcode)==333) return kPhi;
	  else if(TMath::Abs(maPdgcode)==331) return kEtaPrime;
	  else if(TMath::Abs(maPdgcode)==113) return kRho0;
	  else if(TMath::Abs(maPdgcode)==321) return kKe3;
	  else if(TMath::Abs(maPdgcode)==130) return kK0L;
	  else origin = kElse;//mother but nothing we identify
	}
  }//mother is something different from J/psi,charm,beauty or gamma
	
  return origin;
}
//_______________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &hfpt, double &hfeta){

  if(!mcpart) return -1;
  if(!fAODArrayMCInfo) return -1;
  
  int pdgHF = TMath::Abs(mcpart->GetPdgCode());
  hfpt = mcpart->Pt();
  hfeta = mcpart->Eta();
  if(!(pdgHF/100==4 || pdgHF/100==5 || pdgHF/1000==4 || pdgHF/1000==5)) return -1;

  AliAODMCParticle *mctrack = NULL;
  AliAODMCParticle *partMother = NULL;
  
  if(pdgHF==411 || pdgHF==421 || pdgHF==431 || pdgHF==4122 || pdgHF==4132 || pdgHF==4232 || pdgHF==4332){
    // iterate until find B hadron as a mother
    int jLabel = -999;
    int maPdgcode = -999;
    for(int i=1; i<100; i++){
      if(i==1) jLabel = mcpart->GetMother();
      if(i!=1) jLabel = partMother->GetMother();
     
      if(jLabel==-1){
        if(pdgHF==421) return kPromptD0;
        if(pdgHF==4122) return kPromptLc;
      }    
      if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }    
      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))) {
        return -1;
      }    
      maPdgcode = TMath::Abs(mctrack->GetPdgCode());
      if(maPdgcode==511 || maPdgcode==521 || maPdgcode==531 || maPdgcode==5122 || maPdgcode==5132 || maPdgcode==5232 || maPdgcode==5332){
        hfpt = mctrack->Pt();
        hfeta = mctrack->Eta();
        return kNonPromptD;
      }    
      partMother = mctrack;
    }// end of iteration 
  }
  
  // prompt B mesons
  else if(pdgHF==511 || pdgHF==521 || pdgHF==531 || pdgHF==5122 || pdgHF==5132 || pdgHF==5232 || pdgHF==5332){
    return kPromptB;
  }

  return -1;
}
