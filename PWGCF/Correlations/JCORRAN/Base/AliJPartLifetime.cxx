/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

//==================================================================
// Simple class for the pseudo rapidity distribution
// by Taesoo Kim, Beomkyu Kim
// updated by Jasper Parkkila
//==================================================================

#include <TRandom.h>
#include <TMath.h>
#include <TRegexp.h>
#include <TVector.h>

#include "Riostream.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "AliESD.h"
#include "AliAnalysisManager.h"
#include "AliJPartLifetime.h"
#include "AliESDtrackCuts.h"
#include "AliESDFMD.h"
#include "AliStack.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliPWG0Helper.h"

//______________________________________________________________________________
AliJPartLifetime::AliJPartLifetime():
     AliAnalysisTaskSE()
    , fOutput(NULL)
    , fcard(NULL)
    , fTrackCuts(NULL)
    , fEsd(NULL)
    , fPIDResponse(NULL)
    , fEventNumbers(NULL)
{
	//DefineOutput (1, TList::Class());
}

//______________________________________________________________________________
AliJPartLifetime::AliJPartLifetime(const char *name, TString inputformat):
        AliAnalysisTaskSE(name),
    fOutput(NULL),
    fcard(NULL),
    fTrackCuts(NULL),
    fEsd(NULL),
    fPIDResponse(NULL),
    fEventNumbers(NULL)
{
	// Constructor
	AliInfo("---- AliJPartLifetime Constructor ----");
	JUNUSED(inputformat);
	DefineOutput (1, TList::Class());
}

//____________________________________________________________________________
AliJPartLifetime::AliJPartLifetime(const AliJPartLifetime& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fOutput( ap.fOutput ),
	fcard( ap.fcard ),
	fTrackCuts(ap.fTrackCuts),
	fEsd(ap.fEsd),
    fPIDResponse(ap.fPIDResponse),
	fEventNumbers(ap.fEventNumbers)
{

	AliInfo("----DEBUG AliJPartLifetime COPY ----");

}

//_____________________________________________________________________________
AliJPartLifetime& AliJPartLifetime::operator = (const AliJPartLifetime& ap)
{
	AliInfo("----DEBUG AliJPartLifetime operator= ----");
	this->~AliJPartLifetime();
	new(this) AliJPartLifetime(ap);
	return *this;
}

//___________________________________________________________________
AliJPartLifetime::~AliJPartLifetime()
{
	delete fOutput;
	//delete fTrigger;
	//delete[] fTriggerList;
	delete fTrackCuts;
	delete fEventNumbers;
    //delete fPIDResponse;
	for (uint i = 0; i < kAll; i++){
		uint triggc = fcard->GetNoOfBins(kTriggType);
		for(uint j = 0; j < triggc; ++j){
			delete fV0Mass[i][j];
			delete fhTriggPtBin[i][j];
			delete fDecayLength[i][j];
			//delete fDCA[i][j];
			delete fProperTime[i][j];

            delete fV0MassPID[i][j];
            delete fProperTimePID[i][j];
		}

		delete fpTspectra[i];
	}
	//delete fcard;
}

//___________________________________________________________________
void AliJPartLifetime::UserCreateOutputObjects()
{
	// Histograms container
	fOutput = new TList();
	fOutput->SetOwner(kTRUE);

	fcard->PrintOut();
	//-----------------------------------------------------------------------

	// TrackCuts for strangeness measure-------------------------------------
	fTrackCuts = new AliESDtrackCuts();
	//fTrackCuts -> GetStandardITSTPCTrackCuts2010(1,0);
	// same to the dN/dPt note
	fTrackCuts -> SetMaxDCAToVertexXYPtDep("(0.0182+0.0350/pt^1.01)");
	fTrackCuts -> SetMinNCrossedRowsTPC(120);
	fTrackCuts -> SetMaxDCAToVertexZ(2);
	fTrackCuts -> SetEtaRange(-0.8,0.8);
	fTrackCuts -> SetMaxChi2PerClusterTPC(4);
	fTrackCuts -> SetRequireTPCRefit(kTRUE);
	fTrackCuts -> SetRequireITSRefit(kTRUE);
	fTrackCuts -> SetClusterRequirementITS
		(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fTrackCuts -> SetAcceptKinkDaughters(kFALSE);
	fTrackCuts -> SetMaxChi2PerClusterITS(36);
	fTrackCuts -> SetMaxChi2TPCConstrainedGlobal(36);
	fTrackCuts -> SetPtRange(0.15);
	fTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fTrackCuts -> SetMaxFractionSharedTPCClusters(0.4);

	// Vertex distributions--------------------------------------------------

	const char *lname[] = {"AllEvents","PassingPS","AfterPileUpRejection","GoodVertex","GoodVertexCut"};
	const uint binsize = sizeof(lname)/sizeof(lname[0])-1;
	fEventNumbers = new TH1D("EventNumbers","",binsize,0,binsize);
	for (uint i = 0; i < binsize; i++)
		fEventNumbers->GetXaxis()->SetBinLabel(i+1,lname[i]);
	fOutput->Add(fEventNumbers);

	uint triggc = fcard->GetNoOfBins(kTriggType);

	const char *tname[] = {"K0S","Lambda","ALambda"};
	for (uint i = 0; i < kAll; i++){
		//==================================
		//  trigger pt histos
		//==================================
		double ptbw = 10/100.0;  //see hPt histo below, let's make 10 bins per 1GeV/c
		for (uint hit = 0; hit < triggc; hit++){
			double pTt1 = fcard->GetBinBorder(kTriggType, hit);
			double pTt2 = fcard->GetBinBorder(kTriggType, hit + 1);

			uint bc = (uint)TMath::Ceil((pTt2-pTt1)/ptbw);

			fV0Mass[i][hit] = new TH1D(Form("mass_%s_%u",tname[i],hit),tname[i],180*5,0.2,2); //0.2 to 1.5 GeV
			fOutput->Add(fV0Mass[i][hit]);

			fDecayLength[i][hit] = new TH1D(Form("DL_%s_%u",tname[i],hit),"DecayLength",1000,0,100); //cm
			fOutput->Add(fDecayLength[i][hit]);

			/*fDCA[i][hit] = new TH1D(Form("DCA_%s_%u",tname[i],hit),"DCA",1000,0,100); //cm
			fOutput->Add(fDCA[i][hit]);*/

			fProperTime[i][hit] = new TH1D(Form("Tau_%s_%u",tname[i],hit),"Tau",300,0,30);
			fOutput->Add(fProperTime[i][hit]);

            fV0MassPID[i][hit] = new TH1D(Form("mass_PID_%s_%u",tname[i],hit),tname[i],180*5,0.2,2); //0.2 to 1.5 GeV
			fOutput->Add(fV0MassPID[i][hit]);

            fProperTimePID[i][hit] = new TH1D(Form("Tau_PID_%s_%u",tname[i],hit),"Tau PID",300,0,30);
			fOutput->Add(fProperTimePID[i][hit]);

			fhTriggPtBin[i][hit] = new TH1D(Form("hTriggPtBin%02d",hit),Form("pTt: %3.1f-%3.1f",pTt1,pTt2),bc,pTt1,pTt2);
			fhTriggPtBin[i][hit]->Sumw2();
			fOutput->Add(fhTriggPtBin[i][hit]);
		}

		fpTspectra[i] = new TH1D(Form("pTspectra_%s",tname[i]),"pTbins",10,0,10);
		//fOutput->Add(fV0pT[i]);
		//fOutput->Add(fV0p[i]);
		fOutput->Add(fpTspectra[i]);
	}

	OpenFile(1);
	TDirectory *dircard = gDirectory;
	fcard->WriteCard(dircard);

	PostData(1,fOutput);
}

//___________________________________________________________________
void AliJPartLifetime::UserExec(Option_t* )
{
	AliVEvent *event = InputEvent();
	if (!event) {
		Printf("ERROR: Could not retrieve event");
		return;
	}

	// connect to ESD tree --------------------------------------------------
	fEsd = dynamic_cast<AliESDEvent*>(event);
	if (!fEsd) {
		AliError("Cannot get the ESD event");
		return;
	}

	AliMCEvent *mcEvent = MCEvent();

	Bool_t IsMC = kFALSE;
	TArrayF vtxMC3d(3);

	AliStack *stack = NULL;
	if (mcEvent) {
		IsMC = kTRUE;
		stack = mcEvent -> Stack();
		mcEvent->GenEventHeader()->PrimaryVertex(vtxMC3d);
	}

	/*AliGenPythiaEventHeader* pythiaGenHeader = NULL;
	AliGenDPMjetEventHeader* dpmHeader = NULL;

	if(IsMC) { // get access to event generator headers
		pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
		dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(mcEvent->GenEventHeader());

		Double_t eta=-20;
		Int_t nPrim  = stack->GetNprimary();
		for (Int_t i = 0; i < nPrim; ++i){
			//TParticle *  part = (TParticle *) fParticles.At(i);
			TParticle* part = stack->Particle(i);
			eta = part->Eta();
		} // end of the particle loop

	}*/
	// -----------------------------------------------------------------------

	// Load InputHandler for each event---------------------------------------
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)
		AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	// -----------------------------------------------------------------------

    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(!fPIDResponse)
        printf("!fPIDResponse\n");

	// Decide TRIGGER information---------------------------------------------
	fEventNumbers->Fill("AllEvents",1);


	Bool_t IsEventSelected = kFALSE;
	if (inputHandler->IsEventSelected() & AliVEvent::kMB){
		IsEventSelected =  kTRUE;
		fEventNumbers->Fill("PassingPS",1);
	}

	// Reject pile-up events, ALICE official----------------------------------
	if (!IsMC && event->IsPileupFromSPD(3.,0.8,3.,2.,5.))
		return;
	// -----------------------------------------------------------------------


	if (IsEventSelected)
		fEventNumbers->Fill("AfterPileUpRejection",1);

	// SPD Vertices-----------------------------------------------------------
	Bool_t IsGoodVertex = kFALSE;
	Bool_t IsGoodVertexCut = kFALSE;
	const AliESDVertex* vertex_SPD = AliPWG0Helper::GetVertex(fEsd, AliPWG0Helper::kSPD);
	Double_t spdzvtx = -50.;
	Double_t position[3];
	if (vertex_SPD){
		spdzvtx = vertex_SPD->GetZ();
		if (AliPWG0Helper::TestVertex(vertex_SPD, AliPWG0Helper::kSPD)){
			IsGoodVertex = kTRUE;
			vertex_SPD->GetXYZ(position);
			if(TMath::Abs(spdzvtx) < 10.)
				IsGoodVertexCut = kTRUE;
		}
	}

	Bool_t IsGoodTPCVertex = kFALSE;
	Bool_t IsGoodTPCVertexCut = kFALSE;

	const AliESDVertex* vertex_TPC = AliPWG0Helper::GetVertex(fEsd, AliPWG0Helper::kTPCITS);
	Double_t tpczvtx = -50.;
	Double_t tpcposition[3];
	if(vertex_TPC) {
		tpczvtx = vertex_TPC->GetZ();
		if (AliPWG0Helper::TestVertex(vertex_TPC,AliPWG0Helper::kTPCITS)){
			IsGoodTPCVertex = kTRUE;
			vertex_TPC->GetXYZ(tpcposition);
			if(TMath::Abs(tpczvtx) < 10.)
				IsGoodTPCVertexCut = kTRUE;
		}

	}
	if (IsGoodVertexCut && IsGoodTPCVertexCut && TMath::Abs(tpczvtx-spdzvtx) > 0.5)
		IsGoodTPCVertexCut = kFALSE;

	if (IsEventSelected){
		if (IsGoodVertex)
			fEventNumbers->Fill("GoodVertex",1);
		if (IsGoodVertexCut) {
			fEventNumbers->Fill("GoodVertexCut",1);
			this->StrangenessMeasure();
		}
	}

	PostData(1, fOutput);
}
//___________________________________________________________________
void AliJPartLifetime::FinishTaskOutput()
{

}
//___________________________________________________________________
void AliJPartLifetime::Terminate(Option_t*)
{

}

void AliJPartLifetime::StrangenessMeasure(){
	Double_t tPrimaryVtxPosition[3];
	Int_t nv0s = 0;
	nv0s = fEsd->GetNumberOfV0s();
	Int_t    lOnFlyStatus = 0, nv0sOn = 0, nv0sOff = 0;
	Double_t lChi2V0 = 0;
	Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
	Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
	Double_t lV0CosineOfPointingAngle = 0;
	Double_t lV0Radius = 0;
	Double_t lV0DecayLength = 0;
	Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
	Double_t lPt       = 0, lRapK0s = 0, lRapLambda = 0;
	Double_t lAlphaV0  = 0, lPtArmV0 = 0;
	Double_t tV0Position[3];
	Double_t lMagneticField      = 999;
	Double_t lP        = 0;

	const AliESDVertex *primaryVtx = ((AliESDEvent*)fEsd)->GetPrimaryVertex();
	tPrimaryVtxPosition[0] = primaryVtx->GetX();
	tPrimaryVtxPosition[1] = primaryVtx->GetY();
	tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	if (TMath::Abs(tPrimaryVtxPosition[2]) > 10)
		return;

	for(Int_t iV0 = 0; iV0 < nv0s; iV0++){
		AliESDv0 *v0 = fEsd->GetV0(iV0);
		UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());
		lMagneticField = ((AliESDEvent*)fEsd)->GetMagneticField();

		AliESDtrack *pTrack = fEsd->GetTrack(lKeyPos);
		AliESDtrack *nTrack = fEsd->GetTrack(lKeyNeg);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retrieve one of the daughter tracks");
			continue;
		}

		// Remove like-sign
		if ( pTrack->GetSign() == nTrack->GetSign()){
			//cout<< "like sign, continue"<< endl;
			continue;
		}
		// Tracks quality cuts
		if (((pTrack->GetTPCNcls())<80)||((nTrack->GetTPCNcls())<80))
			continue;

		// TPC refit condition (done during reconstruction for
		// Offline but not for On-the-fly)
		if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit))
			continue;

		if (pTrack) lDcaPosToPrimVertex =
			TMath::Abs(pTrack->GetD(
						tPrimaryVtxPosition[0]
						, tPrimaryVtxPosition[1]
						, lMagneticField) );

		if (nTrack) lDcaNegToPrimVertex =
			TMath::Abs(nTrack->GetD(
						tPrimaryVtxPosition[0]
						, tPrimaryVtxPosition[1]
						, lMagneticField) );

		lOnFlyStatus             = v0->GetOnFlyStatus();
		lChi2V0                  = v0->GetChi2V0();
		lDcaV0Daughters          = v0->GetDcaV0Daughters();
		lDcaV0ToPrimVertex       = v0->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);

		lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(
				tPrimaryVtxPosition[0]
				, tPrimaryVtxPosition[1]
				, tPrimaryVtxPosition[2]);
		v0->GetXYZ(tV0Position[0], tV0Position[1], tV0Position[2]);

		lV0Radius = TMath::Sqrt(tV0Position[0]*tV0Position[0]+tV0Position[1]*tV0Position[1]);
		lV0DecayLength = TMath::Sqrt(
				TMath::Power(tV0Position[0] - tPrimaryVtxPosition[0],2) +
				TMath::Power(tV0Position[1] - tPrimaryVtxPosition[1],2) +
				TMath::Power(tV0Position[2] - tPrimaryVtxPosition[2],2));

		// Armenteros variables: !!
		lAlphaV0      = v0->AlphaV0();
		lPtArmV0      = v0->PtArmV0();

		// Selections:
		if ( (lDcaPosToPrimVertex      < 0.05 )||
				(lDcaNegToPrimVertex      < 0.05 )||
				(lDcaV0Daughters          > 0.5 )  ||
				(lV0CosineOfPointingAngle < 0.99)
		   ) continue;


		v0->ChangeMassHypothesis(310);
		lPt = v0->Pt();
		lP  = v0->P();

		int pttbin = fcard->GetBin(kTriggType,lPt);
		if(pttbin < 0)
			continue;

		fV0Mass[kK0s][pttbin]->Fill(v0->GetEffMass());
        int pidn = kUnknown, pidp = kUnknown;
        pidp = GetPID(fPIDResponse,pTrack);
        pidn = GetPID(fPIDResponse,nTrack);

		if (v0->GetEffMass() > 0.482 && v0->GetEffMass() < 0.509) {
			fhTriggPtBin[kK0s][pttbin]->Fill(lPt);
			fpTspectra[kK0s]->Fill(pttbin);
			fDecayLength[kK0s][pttbin]->Fill(lV0DecayLength);
			//fDCA[kK0s][pttbin]->Fill(lV0DecayLength);

            double tau = lV0DecayLength*0.497614/lP; //proper time = M_K0s * decay length / P
			fProperTime[kK0s][pttbin]->Fill(tau);

            //int pidp = fTrackCuts->AcceptTrack(pTrack)?GetPID(fPIDResponse,pTrack):kUnknown;
			//if (fTrackCuts->AcceptTrack(pTrack)) { // good ITSTPC track
			//}
			//if (fTrackCuts->AcceptTrack(nTrack)) { // good ITSTPC track
			//}

            if(pidp == kPion && pidn == kPion)
                fProperTimePID[kK0s][pttbin]->Fill(tau);
		}
        if(pidp == kPion && pidn == kPion)
            fV0MassPID[kK0s][pttbin]->Fill(v0->GetEffMass()); //move this up?

		v0->ChangeMassHypothesis(3122);
		fV0Mass[kLamb][pttbin]->Fill(v0->GetEffMass());

		if (v0->GetEffMass() > 1.11 && v0->GetEffMass() < 1.12) {
			fhTriggPtBin[kLamb][pttbin]->Fill(lPt);
			fpTspectra[kLamb]->Fill(pttbin);

            /*int pidp = fTrackCuts->AcceptTrack(pTrack)?GetPID(fPIDResponse,pTrack):kUnknown;
            int pidn = fTrackCuts->AcceptTrack(pTrack)?GetPID(fPIDResponse,pTrack):kUnknown;
            if(pidp == kProton && pidn == kPion){ //p+ + pi-

                //fProperTimePID[kLamb][pttbin]->Fill(tau);
            }*/
            /*if (fTrackCuts->AcceptTrack(pTrack)) { // good ITSTPC track
			}
			if (fTrackCuts->AcceptTrack(nTrack)) { // good ITSTPC track
			}*/
		}

		v0->ChangeMassHypothesis(-3122);
		fV0Mass[kALamb][pttbin]->Fill(v0->GetEffMass());

		if (v0->GetEffMass() > 1.11 && v0->GetEffMass() < 1.12) {
			fhTriggPtBin[kALamb][pttbin]->Fill(lPt);
			fpTspectra[kALamb]->Fill(pttbin);

            /*int pidp = fTrackCuts->AcceptTrack(pTrack)?GetPID(fPIDResponse,pTrack):kUnknown;
            int pidn = fTrackCuts->AcceptTrack(pTrack)?GetPID(fPIDResponse,pTrack):kUnknown;
            if(pidn == kProton && pidp == kPion) //p- + pi+
                fProperTimePID[kALamb][pttbin]->Fill(tau);*/
			/*if (fTrackCuts->AcceptTrack(pTrack)) { // good ITSTPC track
			}
			if (fTrackCuts->AcceptTrack(nTrack)) { // good ITSTPC track
			}*/
		}
	}
}


Int_t AliJPartLifetime::GetPID(AliPIDResponse *pid, const AliVTrack *trk){
    if(!pid)
        return kUnknown; // no pid available

    Double_t sigmas[] = {-999,-999,-999,-999};

    uint ipid = kUnknown;
    Double_t lsigma = 3.0;
    sigmas[kPion] = pid -> NumberOfSigmasTPC(trk,AliPID::kPion);
    sigmas[kKaon] = pid -> NumberOfSigmasTPC(trk,AliPID::kKaon);
    sigmas[kProton] = pid -> NumberOfSigmasTPC(trk,AliPID::kProton);
    sigmas[kElectron] = pid -> NumberOfSigmasTPC(trk,AliPID::kElectron);
    for (uint i = 0; i < kUnknown; i++){
        if (fabs(sigmas[i]) < lsigma) {
            lsigma = fabs(sigmas[i]);
            ipid = i;
        }
    }

    // derive information, whether tof pid is available
    /*if (0){
        const Bool_t ka = !(trk->GetStatus() & AliESDtrack::kTOFmismatch);
        const Bool_t kb =  (trk->GetStatus() & AliESDtrack::kTOFpid);
        const Bool_t ktof = ka && kb;
    }*/

    return ipid;//(lsigma > 3.0 || ipid == kElectron)?kPion:ipid;
}
