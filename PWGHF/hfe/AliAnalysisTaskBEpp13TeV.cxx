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

///////////////////////////////////////////////////////////////////////
//                                                                   //
//    Task for beauty electron analysis in pp collisions at 13TeV    //
//                                                                   //
//    Authors                                                        //
//    Jonghan Park (jonghan@cern.ch)                                 //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskBEpp13TeV.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "AliHFEextraCuts.h"
#include "AliHFEtools.h"
#include "AliMultSelection.h"
#include "AliGenEventHeader.h"
#include "TAxis.h"
#include "AliKFParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "AliAODv0.h"
#include "AliAODv0KineCuts.h"
#include "AliESDtrack.h"
#include "AliHFEV0taginfo.h"

//______________________________________________________________________
ClassImp(AliAnalysisTaskBEpp13TeV)

//______________________________________________________________________
AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV(const char *name) 
  : AliAnalysisTaskSE(name)
,fIsMC(0)
,fMinTPCNcls(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4.)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(4)
,fITSlayer("kBoth")
,fEta(0.8)
,fMinPt(0.5)
,fMaxPt(30.)
,fDCAxy(0.1)
,fDCAz(0.2)
,fTPCnsigmaLow(0.)
,fTPCnsigmaHigh(3.)
,fTOFnsigma(3.)

,fAOD(0)
,fOutputList(0)
,fPidResponse(0)
,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fExtraCuts(0)
,fAODMCParticle(0)
,fAODv0(0)
,fAODV0Cuts(0)
,fV0Tagger(0)

,hVertexZ(0)
,hNrEvents(0)

,hFilterMask(0)
,hTPCNcls(0)
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

,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)
,hV0ElecTPCnsigmaDeno(0)
,hV0ElecTPCnsigmaNume(0)

,dcaTrack(0)
,dcaPion(0)

,fRnd(0)
{
  //Named constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
	fV0Tagger = new AliHFEV0taginfo("Tagger");
}

//________________________________________________________________________
AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV() 
  : AliAnalysisTaskSE()

,fIsMC(0)
,fMinTPCNcls(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4.)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(4)
,fITSlayer("kBoth")
,fEta(0.8)
,fMinPt(0.5)
,fMaxPt(30.)
,fDCAxy(0.1)
,fDCAz(0.2)
,fTPCnsigmaLow(0.)
,fTPCnsigmaHigh(3.)
,fTOFnsigma(3.)

,fAOD(0)
,fOutputList(0)
,fPidResponse(0)
,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fExtraCuts(0)
,fAODMCParticle(0)
,fAODv0(0)
,fAODV0Cuts(0)
,fV0Tagger(0)

,hVertexZ(0)
,hNrEvents(0)

,hFilterMask(0)
,hTPCNcls(0)
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

,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)
,hV0ElecTPCnsigmaDeno(0)
,hV0ElecTPCnsigmaNume(0)

,dcaTrack(0)
,dcaPion(0)

,fRnd(0)
{
	// default constructor
}

//______________________________________________________________________
AliAnalysisTaskBEpp13TeV::~AliAnalysisTaskBEpp13TeV()
{
	//Destructor 
	if(fOutputList) delete fOutputList;
	if(fExtraCuts) delete fExtraCuts;
	if(fV0Tagger) delete fV0Tagger;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTaskBEpp13TeV::UserCreateOutputObjects()
{
  
	fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
	fAODV0Cuts = new AliAODv0KineCuts();

	fOutputList = new TList();
	fOutputList->SetOwner();	

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
	
	int nBinsInvMassK0s = 100;
	double minInvMassK0s = 0;
	double maxInvMassK0s = 1;
	double binLimInvMassK0s[nBinsInvMassK0s+1];
	for(int i=0; i<=nBinsInvMassK0s; i++) binLimInvMassK0s[i] = minInvMassK0s + double(i)*(maxInvMassK0s/nBinsInvMassK0s);
	
	int nBinsProdRadi = 200;
	double minProdRadi = 0.;
	double maxProdRadi = 20.;
	double binLimProdRadi[nBinsProdRadi+1];
	for(int i=0; i<=nBinsProdRadi; i++) binLimProdRadi[i] = minProdRadi + double(i)*(maxProdRadi/nBinsProdRadi);
	
	int nBinsCent = 20;
	double minCent = 0.;
	double maxCent = 100.;
	double binLimCent[nBinsCent+1];
	for(int i=0; i<=nBinsCent; i++) binLimCent[i] = minCent + double(i)*(maxCent/nBinsCent);
	
	int nBinsMult = 60;
	double minMult = 0.;
	double maxMult = 30000.;
	double binLimMult[nBinsMult+1];
	for(int i=0; i<=nBinsMult; i++) binLimMult[i] = minMult + double(i)*(maxMult/nBinsMult);
	
	// event qa
	hVertexZ = new TH1F("hVertexZ", "", 600, -30, 30);
	fOutputList->Add(hVertexZ);

	hNrEvents = new TH1F("hNrEvents", "Number of Events", 1, 0, 1);
	fOutputList->Add(hNrEvents);

	hFilterMask = new TH1F("hFilterMask", "", 2, 0., 2.);
	fOutputList->Add(hFilterMask);

	hTPCNcls = new TH1F("hTPCNcls", "", 200, 0., 200.);
	fOutputList->Add(hTPCNcls);

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

	hTPCnsigma = new TH2F("hTPCnsigma", "TPC n#sigma", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigma);

	hTPCnsigmaTOFcut = new TH2F("hTPCnsigmaTOFcut", "TPC n#sigma after TOF cut", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaTOFcut);

	hTPCnsigmaTOFcutPt = new TH2F("hTPCnsigmaTOFcutPt", "TPC n#sigma after TOF cut", nPtBins, ptbinningX, nBinsPID, binLimPID);
	fOutputList->Add(hTPCnsigmaTOFcutPt);

	hTPCnsigmaQA = new TH2F("hTPCnsigmaQA", "TPC pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaQA);

	hTPCnsigmaPiQA = new TH2F("hTPCnsigmaPiQA", "TPC pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaPiQA);

	hTOFnsigma = new TH2F("hTOFnsigma", "TOF n#sigma", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTOFnsigma);

	hTOFnsigmaQA = new TH2F("hTOFnsigmaQA", "TOF pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTOFnsigmaQA);

	hV0ElecTOFnsigmaDeno = new TH2F("hV0ElecTOFnsigmaDeno", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTOFnsigmaDeno);

	hV0ElecTOFnsigmaNume = new TH2F("hV0ElecTOFnsigmaNume", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTOFnsigmaNume);

	hV0ElecTPCnsigmaDeno = new TH2F("hV0ElecTPCnsigmaDeno", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTPCnsigmaDeno);

	hV0ElecTPCnsigmaNume = new TH2F("hV0ElecTPCnsigmaNume", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTPCnsigmaNume);

	dcaTrack = new TH2F("dcaTrack", "inclusive electron's dca", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaTrack);

	dcaPion = new TH2F("dcaPion", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaPion);

	fRnd = new TRandom3(0);

	PostData(1, fOutputList);
	
}

//Main loop
//Called for each event
void AliAnalysisTaskBEpp13TeV::UserExec(Option_t *) 
{
	//Check Event
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if(!fAOD)
	{
		printf("ERROR: fAOD not available\n");
		return;
	}

	if(fIsMC){
		fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
		if(!fAODMCHeader){
			AliError("No AliAODMCHeader");
			return;
		}
		
		fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
		if(!fAODArrayMCInfo){
			AliError("No AOD MC particles");
			return;
		}
		if(fAODArrayMCInfo->GetEntries() < 1) return;
	}
	
	//Check HFEextraCut
	if(!fExtraCuts)
	{
		fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
	}
	fExtraCuts->SetRecEventInfo(fAOD);

	//PID response
	fPidResponse = fInputHandler->GetPIDResponse();
	if(!fPidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
	}
	
	//Initialize V0 electron tagger
	if(fV0Tagger){
		fV0Tagger->Reset();
		fV0Tagger->TagV0Tracks(fAOD);
	}

	//=================
	// Event selection
	//=================
	if(!PassEventCuts(fAOD)) return;		
	
	//==============================
	// Look for kink mother for AOD
	//==============================
	double *fListOfmotherkink = 0;
	int fNumberOfVertices = 0; 
	int fNumberOfMotherkink = 0;

	fNumberOfVertices = fAOD->GetNumberOfVertices();

	fListOfmotherkink = new double[fNumberOfVertices];

	for(int ivertex=0; ivertex < fNumberOfVertices; ivertex++) 
	{
		AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
		if(!aodvertex) continue;
		if(aodvertex->GetType()==AliAODVertex::kKink) 
		{
			AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
			if(!mother) continue;
			int idmother = mother->GetID();
			fListOfmotherkink[fNumberOfMotherkink] = idmother;
			fNumberOfMotherkink++;
		}
	}

	//--------------------
	double fSignB = -999.;
	if(fAOD->GetMagneticField()<0) fSignB = -1;
	if(fAOD->GetMagneticField()>0) fSignB = 1;
	
	hNrEvents->Fill(0);

	//=======================================================================
	///Track loop
	for(int iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
	{

		AliAODTrack *aodTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
		
		double pt = aodTrack->Pt();
	
		//=========    
		// RecKink
		//=========
		bool kinkmotherpass = true;
		for(int kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++) 
		{
			if(aodTrack->GetID() == fListOfmotherkink[kinkmother]) 
			{
				kinkmotherpass = kFALSE;
				continue;
			}
		}
		if(!kinkmotherpass) continue;
		
		double hfeImpactParam = -999., hfeImpactParamResol = -999.;
		fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);

		if(!PassTrackCuts(aodTrack)) continue;
		
		//=======================================================================
		// QA plots after track selection
		//=======================================================================
	
		double fTPCnSigma = fPidResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
		double fTOFnSigma = fPidResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);
		double fTPCnSigmaPi = fPidResponse->NumberOfSigmasTPC(aodTrack, AliPID::kPion);
		double fTOFnSigmaPi = fPidResponse->NumberOfSigmasTOF(aodTrack, AliPID::kPion);
		
		hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
		hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);

		//V0 electrons from systematic studies of TOF PID cut
		AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(aodTrack->GetID()); /// enum EParticleType: kElectron = 0, kMuon = 1, kPion = 2, etc
		if(myv0pid == AliPID::kElectron){
			// TOF eID systematics
			if(fTPCnSigma >= fTPCnsigmaLow && fTPCnSigma <= fTPCnsigmaHigh){
				hV0ElecTOFnsigmaDeno->Fill(aodTrack->Pt(), fTOFnSigma);
				if(TMath::Abs(fTOFnSigma) <= fTOFnsigma) hV0ElecTOFnsigmaNume->Fill(aodTrack->Pt(), fTOFnSigma);
			}
			// TPC eID systematics
			if(TMath::Abs(fTOFnSigma) <= fTOFnsigma){
				hV0ElecTPCnsigmaDeno->Fill(aodTrack->Pt(), fTPCnSigma);
				if(fTPCnSigma >= fTPCnsigmaLow && fTPCnSigma <= fTPCnsigmaHigh) hV0ElecTPCnsigmaNume->Fill(aodTrack->Pt(), fTPCnSigma);
			}
		}


		if(TMath::Abs(fTOFnSigma) > fTOFnsigma) continue;

		hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
		hTPCnsigmaTOFcutPt->Fill(aodTrack->Pt(), fTPCnSigma);

		if(fTPCnSigma > -5 && fTPCnSigma < -3){
			hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
			dcaPion->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
		}

		if(fTPCnSigma < fTPCnsigmaLow || fTPCnSigma > fTPCnsigmaHigh) continue;
		
		hTPCnsigmaQA->Fill(aodTrack->P(), fTPCnSigma);
		hTOFnsigmaQA->Fill(aodTrack->P(), fTOFnSigma);

		
		dcaTrack->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());

	}//End of track loop
	
	//=======================================================================
	delete fListOfmotherkink;
	PostData(1, fOutputList);
}      

//=======================================================================
void AliAnalysisTaskBEpp13TeV::Terminate(Option_t *) 
{
//Draw result to the screen
//Called once at the end of the query

	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	
	if(!fOutputList) 
	{
		printf("ERROR: Output list not available\n");
		return;
	}
}

//=======================================================================

//_________________________________________
bool AliAnalysisTaskBEpp13TeV::PassEventCuts(AliAODEvent *event){

	//event selection cuts
	AliAODVertex *vtx = event->GetPrimaryVertex();
	if(vtx->GetNContributors()<2) return false;

	AliAODVertex *vtxSPD = event->GetPrimaryVertexSPD();
	double cov[6]={0};
	vtxSPD->GetCovarianceMatrix(cov);
	double zRes = TMath::Sqrt(cov[5]);
	if(vtxSPD->IsFromVertexerZ() && (zRes>0.25)) return false;
	
	// To be confirmed
	double zvtx = vtx->GetZ();
	if(TMath::Abs(zvtx) > 10) return false;
	
	double zvtxSPD = vtxSPD->GetZ();
	if(TMath::Abs(zvtx-zvtxSPD)>0.5) return false;
	
	return true;
}

bool AliAnalysisTaskBEpp13TeV::PassTrackCuts(AliAODTrack *track){

	if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;
	if(TMath::Abs(track->Eta()) > fEta) return false;

	// basic tracking
	ULong_t status = track->GetStatus();
	if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;
	
	// dca cut
	float dcaxy = -999.; float dcaz = -999.;
	fExtraCuts->GetImpactParameters((AliVTrack *)track, dcaxy, dcaz);
	if(TMath::Abs(dcaxy) > fDCAxy || TMath::Abs(dcaz) > fDCAz) return false;
	
	// TPC cut
	unsigned short findableTPC = track->GetTPCNclsF();
	unsigned short nclustersTPC = track->GetTPCNcls();
	double FoundOverFindable = (findableTPC ? static_cast<float>(nclustersTPC)/static_cast<float>(findableTPC) : 0);
	if(track->GetTPCNcls() < fMinTPCNcls || track->GetTPCsignalN() < fMinTPCNclsPID || track->Chi2perNDF() > fMaxTPCchi2 || FoundOverFindable < fMinTPCclsRatio) return false;
	
	// ITS cut
	if(track->GetITSNcls() < fMinITSNcls) return false;
	if(fITSlayer.Contains("kFirst")){
		if(!(track->HasPointOnITSLayer(0))) return false;
		hITSlayer->Fill(1);
	}else if(fITSlayer.Contains("kAny")){
		if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(2);
	}else if(fITSlayer.Contains("kBoth")){
		if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(3);
	}else return false;
	
	hFilterMask->Fill(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA));
	hEta->Fill(track->Eta());
	hPhi->Fill(track->Phi());
	hPt->Fill(track->Pt());
	hDCAxy->Fill(dcaxy);
	hDCAz->Fill(dcaz);
	hTPCNcls->Fill(track->GetTPCNcls());
	hTPCclsPID->Fill(track->GetTPCsignalN());
	hTPCchi2->Fill(track->Chi2perNDF());
	hTPCclsRatio->Fill(FoundOverFindable);
	hITSNcls->Fill(track->GetITSNcls());
	return true;
}


//_________________________________________

int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	if(TMath::Abs(mcpart->GetPdgCode()) != 11 ) return kMisID;

	int origin = -1;
	Bool_t isFinalOpenCharm = kFALSE;

	int iLabel = mcpart->GetMother();
	if ((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
		AliDebug(1, "label is out of range, return\n");
		return -1;
	}
	
	AliAODMCParticle *mctrack = NULL; // will change all the time
	int tmpMomLabel=0;
	if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1;
	AliAODMCParticle *partMother = mctrack;	//mtrack 
	AliAODMCParticle *partMotherCopy = mctrack;	//mtrack
	int maPdgcode = mctrack->GetPdgCode();	//mpdg
	mpt = partMother->Pt();	//mpt
	mpdg = partMother->GetPdgCode();	//mpdg
	int grmaPdgcode;
	int ggrmaPdgcode;
	double gmpt, ggmpt;
	int gmpdg, ggmpdg;

	// if the mother is charmed hadron
	if( (int(TMath::Abs(maPdgcode)/100.)%10) == 4 || (int(TMath::Abs(maPdgcode)/1000.)%10) == 4 ) {
		if(TMath::Abs(maPdgcode)==411 || TMath::Abs(maPdgcode)==421 || TMath::Abs(maPdgcode)==431 || TMath::Abs(maPdgcode)==4122){
			mpt = partMother->Pt();
			mpdg = partMother->GetPdgCode();
			isFinalOpenCharm = kTRUE;
		}
		if (!isFinalOpenCharm) {
			return -1;
		}
		
		// iterate until find B hadron as a  mother
		for (int i=1; i<100; i++){
			int jLabel = partMother->GetMother();
			if (jLabel == -1) {
				return kDirectCharm;
			}
			if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
				AliDebug(1, "Stack label is negative, return\n");
				return -1;
			}
			
			if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))) {
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
	else if( (int(TMath::Abs(maPdgcode)/100.)%10) == 5 || (int(TMath::Abs(maPdgcode)/1000.)%10) == 5 ) {
		if (TMath::Abs(maPdgcode)==511 || TMath::Abs(maPdgcode)==521 || TMath::Abs(maPdgcode)==531 || TMath::Abs(maPdgcode)==5122){
			mpt = partMotherCopy->Pt();
			mpdg = partMotherCopy->GetPdgCode();
			return kDirectBeauty;
		}
	}
	
	// if the mother is gamma
	else if ( TMath::Abs(maPdgcode) == 22 ) {
		
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
		partMother = mctrack; // grand mother
		partMotherCopy = mctrack; // grand mother
		mpt = partMother->Pt(); // grand mother pT
		mpdg = partMother->GetPdgCode();
		maPdgcode = partMother->GetPdgCode(); // grand mother PDG
		
		// check if the ligth meson is the decay product of heavy mesons
		tmpMomLabel = partMother->GetMother();
		if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack; //grand grand mother
        grmaPdgcode = partMother->GetPdgCode(); //grand grand mother PDG
				mpt = partMother->Pt(); // grand grand mother PDG
				mpdg = partMother->GetPdgCode();
				gmpt = partMother->Pt(); // grand grand mother PDG
				gmpdg = partMother->GetPdgCode();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 5 ) {
					return kGammaB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 4 ) {
					return kGammaD2M;
				}
				
				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack; // grand grand grand mother
						ggrmaPdgcode = partMother->GetPdgCode(); // grand grand grand mother PDG
						mpt = partMother->Pt(); // grand grand grand mother PDG
						mpdg = partMother->GetPdgCode();
						ggmpt = partMother->Pt(); // grand grand grand mother PDG
						ggmpdg = partMother->GetPdgCode();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 5 ) {
							return kGammaB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 4 ) {
							return kGammaD2M;
						}
					}
				}
				
				if ( TMath::Abs(maPdgcode) == 111 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					else if(grmaPdgcode == 310) return kGammaK0s2P;
					else if(grmaPdgcode == 130) return kGammaK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
					else if(grmaPdgcode == 3222) return kGammaSigma2P;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaRho0;
				}
				else origin = kElse;//grandgrandmother but nothing we identify
			}//mctrack grandgrandmother
		}
		else {
			// grandmother is primary
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kGammaPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kGammaEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kGammaOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kGammaPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kGammaEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kGammaRho0;
			}
			else origin = kElse;//grandmother is primary but nothing we identify
		}
		return origin;
	}

	// if the mother is light meson
	else {
		
		tmpMomLabel = partMotherCopy->GetMother();
		mpt = partMotherCopy->Pt(); // mother pT
		mpdg = partMotherCopy->GetPdgCode();
		if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack; // grand mother
				grmaPdgcode = partMother->GetPdgCode(); // grand mother PDG
				mpt = partMother->Pt(); // grand mother pT
				mpdg = partMother->GetPdgCode();
				gmpt = partMother->Pt(); // grand mother pT
				gmpdg = partMother->GetPdgCode();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 5 ) {
					return kB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 4 ) {
					return kD2M;
				}
				
				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack; // grand grand mother
						ggrmaPdgcode = partMother->GetPdgCode(); // grand grand mother PDG
						mpt = partMother->Pt(); // grand grand mother pT
						mpdg = partMother->GetPdgCode();
						ggmpt = partMother->Pt(); // grand grand mother pT
						ggmpdg = partMother->GetPdgCode();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 5 ) {
							return kB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 4 ) {
							return kD2M;
						}
					}
				}
				
				if ( TMath::Abs(maPdgcode) == 111 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					else if(grmaPdgcode == 310) return kK0s2P;
					else if(grmaPdgcode == 130) return kK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
					else if(grmaPdgcode == 3222) return kSigma2P;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kRho0;
				}
				else if ( TMath::Abs(maPdgcode) == 321 ) {
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kKe3;
				}
				else if ( TMath::Abs(maPdgcode) == 130 ) {
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kK0L;
				}
				else origin = kElse;//grandmother but nothing we identidy
			}//mctrack grandmother
		}
		else {
			// no grandmother
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kRho0;
			}
			else if ( TMath::Abs(maPdgcode) == 321 ) {
				return kKe3;
			}
			else if ( TMath::Abs(maPdgcode) == 130 ) {
				return kK0L;
			}
			else origin = kElse;//mother but nothing we identify
		}
	}//mother is something different from J/psi,charm,beauty or gamma
	
	return origin;

}

//_________________________________________

int AliAnalysisTaskBEpp13TeV::GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &mpt, double &meta){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	int pdgHF = TMath::Abs(mcpart->GetPdgCode());
	mpt = mcpart->Pt();
	meta = mcpart->Eta();
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
				mpt = mctrack->Pt();
				meta = mctrack->Eta();
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

int AliAnalysisTaskBEpp13TeV::GetGammaPt(const AliAODMCParticle * const mcpart, double &mpt){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	if(TMath::Abs(mcpart->GetPdgCode())!=11) return -1;

	int moLabel = mcpart->GetMother();
	if ((moLabel<0) || (moLabel>=fAODArrayMCInfo->GetEntriesFast())){
		AliDebug(1, "label is out of range, return\n");
		return -1;
	}
	
	AliAODMCParticle *mother = NULL;
	if(!(mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(moLabel))))) return -1;
	int moPdg = mother->GetPdgCode();	//mpdg
	
	// if the mother is gamma
	if(TMath::Abs(moPdg)==22){
		mpt = mother->Pt();	//mpt
		
		int gmoLabel=0;
		gmoLabel = mother->GetMother();  // mother of photon
		if(gmoLabel==-1) return kDirectGamma;  // no grandmother
		if((gmoLabel<0) || (gmoLabel>=fAODArrayMCInfo->GetEntriesFast())){
			return -1;
		}
		AliAODMCParticle *gmother = NULL;
		if(!(gmother=dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(gmoLabel))))){
			return -1;
		}
		int gmoPdg = gmother->GetPdgCode(); // grand mother PDG
		if(gmoPdg/100==4||gmoPdg/1000==4) return -1;
		if(gmoPdg/100==5||gmoPdg/1000==5) return -1;
		int ggmoLabel=0;
		ggmoLabel = gmother->GetMother();
		if((ggmoLabel>=0) && (ggmoLabel<fAODArrayMCInfo->GetEntriesFast())){//grandgrandmother
			AliAODMCParticle *ggmother = NULL;
			if((ggmother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(ggmoLabel))))){
        int ggmoPdg = ggmother->GetPdgCode(); //grand grand mother PDG
				if(ggmoPdg/100==4||ggmoPdg/1000==4) return -1;
				if(ggmoPdg/100==5||ggmoPdg/1000==5) return -1;

				if(TMath::Abs(gmoPdg)==111){
					if(ggmoPdg==221 || ggmoPdg==223 || ggmoPdg==333 || ggmoPdg==331 || ggmoPdg==113) return -1;
					else if(ggmoPdg == 310 || ggmoPdg == 130 || TMath::Abs(ggmoPdg) == 321 || TMath::Abs(ggmoPdg) == 3122 || ggmoPdg == 3222) return kBaryonGamma;
					return kDalitzGamma;
				}
				else if(TMath::Abs(gmoPdg)==221){
					if(ggmoPdg==111 || ggmoPdg==223 || ggmoPdg==333 || ggmoPdg==331 || ggmoPdg==113) return -1;
					return kDalitzGamma;
				}
			}
			int gggmoLabel=0;
			gggmoLabel = ggmother->GetMother();
			if((gggmoLabel>=0) && (gggmoLabel<fAODArrayMCInfo->GetEntriesFast())){
				AliAODMCParticle *gggmother = NULL;
				if((gggmother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(gggmoLabel))))){
					int gggmoPdg = gggmother->GetPdgCode();
					if(gggmoPdg/100==4||gggmoPdg/1000==4) return -1;
					if(gggmoPdg/100==5||gggmoPdg/1000==5) return -1;
				}
			}
		}else{
			if(TMath::Abs(gmoPdg)==111) return kDalitzGamma;
			else if(TMath::Abs(gmoPdg)==221) return kDalitzGamma;
		}
		return -1;
	}
	return -1;
}







