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

//______________________________________________________________________________
// Analysis task for high pt particle correlations
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects
//////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
//#include "TList.h"
//#include "TTree.h"
#include "TFile.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliGenHijingEventHeader.h"
#include "AliJFFlucTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliJTrack.h"
#include "AliJMCTrack.h"
//#include "AliJPhoton.h"
#include "AliJEventHeader.h"
#include "AliJHistManager.h"
#include "AliInputEventHandler.h"
#include "AliJEfficiency.h"
#include "AliMultSelection.h"
#include "AliJRunTable.h"
//#pragma GCC diagnostic warning "-Wall"
//______________________________________________________________________________
AliJFFlucTask::AliJFFlucTask():
	AliAnalysisTaskSE(),
	fInputList(0),
	fOutput(0),
	fFFlucAna(0),
	h_ratio(0)
{
	fEvtNum=0;
	fFilterBit = 0;
	fEta_min = 0;
	fEta_max = 0;
	fDebugLevel = 0;
	fEffMode =0;
	fEffFilterBit=0;
	fPt_min=0;
	fPt_max=0;
	fCentDetName="V0M";
	fInFileName="";
	IsMC = kFALSE;
	IsKineOnly = kFALSE;
	IsExcludeWeakDecay = kFALSE;
	IsCentFlat = kFALSE;
	IsPhiModule = kFALSE;
	IsSCptdep = kFALSE;
	IsSCwithQC= kFALSE;
	IsEbEWeighted = kFALSE;
	fQC_eta_min=-0.8;
	fQC_eta_max=0.8;

	for(int icent=0; icent<7; icent++){
		for(int isub=0; isub<2; isub++){
			h_ModuledPhi[icent][isub]=NULL;
		}
	}
	fzvtxCut = 10;  // default z vertex cut
	//  DefineOutput(1, TDirectory::Class());
}

//______________________________________________________________________________
AliJFFlucTask::AliJFFlucTask(const char *name,  Bool_t IsMC, Bool_t IsExcludeWeakDecay):
	AliAnalysisTaskSE(name),
	fInputList(0),
	fOutput(0),
	fFFlucAna(0x0),
	h_ratio(0)
{
	DefineOutput(1, TDirectory::Class());
	fTaskName = name;

	fEvtNum=0;
	fFilterBit = 0;
	fEta_min = 0;
	fEta_max = 0;
	fDebugLevel = 0;
	fEffMode =0;
	fEffFilterBit=0;
	fPt_min=0;
	fPt_max=0;
	fInFileName="";
	fCentDetName="V0M";
	IsMC = kFALSE;
	IsKineOnly = kFALSE;
	IsExcludeWeakDecay = kFALSE;
	IsCentFlat = kFALSE;
	IsPhiModule = kFALSE;
	IsSCptdep = kFALSE;
	IsSCwithQC = kFALSE;
	IsEbEWeighted = kFALSE;

	fQC_eta_min=-0.8;
	fQC_eta_max=0.8;

	for(int icent=0; icent<7; icent++){
		for(int isub=0; isub<2; isub++){
			h_ModuledPhi[icent][isub]=NULL;
		}
	}
	fzvtxCut = 10;
}

//____________________________________________________________________________
AliJFFlucTask::AliJFFlucTask(const AliJFFlucTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fInputList(ap.fInputList),
	fOutput(ap.fOutput),
	fFFlucAna(ap.fFFlucAna)
{
	AliInfo("----DEBUG AliJFFlucTask COPY ----");
}

//_____________________________________________________________________________
AliJFFlucTask& AliJFFlucTask::operator = (const AliJFFlucTask& ap)
{
	// assignment operator
	AliInfo("----DEBUG AliJFFlucTask operator= ----");
	this->~AliJFFlucTask();
	new(this) AliJFFlucTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJFFlucTask::~AliJFFlucTask()
{
	//delete pfOutlierLowCut;
	//delete pfOutlierHighCut;
	delete fFFlucAna;
	delete fInputList;
	delete fOutput;
	delete h_ratio;
}

//________________________________________________________________________
void AliJFFlucTask::UserCreateOutputObjects()
{
	fFFlucAna =  new AliJFFlucAnalysis( fTaskName );
	fFFlucAna->SetDebugLevel(fDebugLevel);
	fFFlucAna->SetIsPhiModule( IsPhiModule);
	fFFlucAna->SetIsSCptdep( IsSCptdep ) ;
	fFFlucAna->SetSCwithQC( IsSCwithQC );
	fFFlucAna->SetEbEWeight( IsEbEWeighted );
	fFFlucAna->SetQCEtaCut( fQC_eta_min, fQC_eta_max, 0.5 );

//	fFFlucAna->SetSCwithFineCentbin( IsSCwithFineCentBin );
	// setting histos for phi modulation
	if( IsPhiModule==kTRUE){
		for(int icent=0; icent<7; icent++){
			for(int isub=0; isub<2; isub++){
				fFFlucAna->SetPhiModuleHistos( icent, isub, h_ModuledPhi[icent][isub] );
			}
		}
	}

	fInputList = new TClonesArray("AliJBaseTrack" , 2500);
	fInputList->SetOwner(kTRUE);

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();
	fFFlucAna->SetEffConfig( fEffMode, fEffFilterBit );
	fFFlucAna->UserCreateOutputObjects();

	PostData(1, fOutput);
}

//______________________________________________________________________________
void AliJFFlucTask::UserExec(Option_t* /*option*/)
{
	// Processing of one event
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));
	// initializing variables from last event
	fInputList->Clear();

	float fCent = -1.0f;
	float fImpactParameter = -1.0f;
	double fvertex[3];

	fEvtNum++;
	if(fEvtNum % 100 == 0)
		cout << "evt : " << fEvtNum <<endl;

	// load current event and save track, event info
	if(IsKineOnly) {
		AliMCEvent*  mcEvent = MCEvent();
		if (!mcEvent) {
			AliError("ERROR: mcEvent not available");
			return;
		}

		AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
		if (headerH) {
			//Double_t gReactionPlane = headerH->ReactionPlaneAngle();
			Double_t gImpactParameter = headerH->ImpactParameter();
			fCent = GetCentralityFromImpactPar(gImpactParameter);
			if( fALICEIPinfo == kTRUE){
				//force to use ALICE impact parameter setting
				double ALICE_Cent[8] = {0, 5, 10, 20, 30, 40, 50, 60};
				double ALICE_IPinfo[8] = {0, 3.50, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09};
				for(int icent=0; icent<8; icent++){
					if(fImpactParameter >= ALICE_IPinfo[icent] && fImpactParameter < ALICE_IPinfo[icent+1])
						fCent = 0.5f*(ALICE_Cent[icent]+ALICE_Cent[icent+1]);
				}
			}
		}
		if( fEvtNum == 1 ){
			int runN = 1234;
			fFFlucAna->GetAliJEfficiency()->SetRunNumber ( runN );
			fFFlucAna->GetAliJEfficiency()->Load();
		}
		ReadKineTracks( mcEvent, fInputList ) ; // read tracklist
		AliGenEventHeader *header = mcEvent->GenEventHeader();
		if(!header)
			return;
		TArrayF gVertexArray;
		header->PrimaryVertex(gVertexArray);
		for(int i = 0; i < 3; i++)
			fvertex[i] = gVertexArray.At(i);
		fFFlucAna->Init();
		fFFlucAna->SetInputList( fInputList );
		fFFlucAna->SetEventCentrality( fCent );
		fFFlucAna->SetEventImpactParameter( fImpactParameter );
		fFFlucAna->SetEventVertex( fvertex );
		fFFlucAna->SetEtaRange( fEta_min, fEta_max ) ;
		fFFlucAna->UserExec(""); // doing some analysis here.

	} else { // Kine
		AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
		fCent = ReadCentrality(currentEvent,fCentDetName);
		//fCent = ReadAODCentrality( currentEvent, fCentDetName  ) ;
		//fCent = ReadMultSelectionCentrality(currentEvent,fCentDetName);
		if( fEvtNum == 1 ){
			int runN = currentEvent->GetRunNumber();
			fFFlucAna->GetAliJEfficiency()->SetRunNumber ( runN );
			fFFlucAna->GetAliJEfficiency()->Load();
		}

		if( IsGoodEvent( currentEvent )){
			ReadAODTracks( currentEvent, fInputList ) ; // read tracklist
			ReadVertexInfo( currentEvent, fvertex); // read vertex info
			// Analysis Part
			fFFlucAna->Init();
			fFFlucAna->SetInputList( fInputList );
			fFFlucAna->SetEventCentrality( fCent );
			fFFlucAna->SetEventImpactParameter( fImpactParameter); // need this??
			fFFlucAna->SetEventVertex( fvertex );
			fFFlucAna->SetEtaRange( fEta_min, fEta_max );
			fFFlucAna->SetEventTracksQA( TPCTracks, GlobTracks);
			fFFlucAna->SetEventFB32TracksQA( FB32Tracks, FB32TOFTracks );
			fFFlucAna->UserExec(""); // doing some analysis here.
			//
		}
	} // AOD
}

//______________________________________________________________________________
void AliJFFlucTask::ReadVertexInfo ( AliAODEvent *aod, double*  fvtx )
{
	fvtx[0] = aod->GetPrimaryVertex()->GetX();
	fvtx[1] = aod->GetPrimaryVertex()->GetY();
	fvtx[2] = aod->GetPrimaryVertex()->GetZ();
}
//______________________________________________________________________________
/*float AliJFFlucTask::ReadAODCentrality( AliAODEvent *aod, TString Trig )
{
	AliCentrality *cent = aod->GetCentrality();
	if(!cent){ Printf("No Cent event"); return -999; }
	return cent->GetCentralityPercentile(Trig.Data());
}
//______________________________________________________________________________
float AliJFFlucTask::ReadMultSelectionCentrality( AliAODEvent *aod, TString Trig )
{
	AliMultSelection *pms = (AliMultSelection*)aod->FindListObject("MultSelection");
	if(!pms){ Printf("No MultSelection task"); return -999.0f; } //NaN?
	return pms->GetMultiplicityPercentile(Trig.Data());
}*/
//______________________________________________________________________________

float AliJFFlucTask::ReadCentrality( AliAODEvent *aod, TString Trig ){
	AliMultSelection *pms = (AliMultSelection*)aod->FindListObject("MultSelection");
	return pms?pms->GetMultiplicityPercentile(Trig.Data()):aod->GetCentrality()->GetCentralityPercentile(Trig.Data());
}

//______________________________________________________________________________

void AliJFFlucTask::Init()
{
	AliInfo("Doing initialization") ;
}
//______________________________________________________________________________
void AliJFFlucTask::ReadAODTracks( AliAODEvent *aod , TClonesArray *TrackList)
{
	//aod->Print();
	if( IsMC == kTRUE ){  // how to get a flag to check  MC or not !
		TClonesArray *mcArray = (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray){ Printf("Error not a proper MC event"); };  // check mc array

		Int_t nt = mcArray->GetEntriesFast();
		Int_t ntrack =0;
		for( int it=0; it < nt ; it++){
			AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
			if(!track) { Error("ReadEventAODMC", "Could not receive particle %d",(int) it); continue; };
			if( track->IsPhysicalPrimary() ){
				// insert AMTP weak decay switch here
				if(IsExcludeWeakDecay == kTRUE){
					//cout << "finding weak decaying particle ... " << endl;
					Int_t gMotherIndex = track->GetMother(); // check and ask about this to DJ changed to mother from firstmother
					if(gMotherIndex != -1) {
						//cout << "this mother is " << gMotherIndex << endl;
						AliAODMCParticle *motherParticle = (AliAODMCParticle*)mcArray->At(gMotherIndex);
						//cout << "mother pdg code is " << motherParticle->GetPdgCode() << endl;
						if(motherParticle) {
							if(IsThisAWeakDecayingParticle(motherParticle)){
								//Exclude the decay products of weakly decaying particles
								continue;
							}
						}
					}
				} // weak decay particles are exclude

				if( fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}
				Int_t pdg = track->GetPdgCode();
				Char_t ch = (Char_t) track->Charge();
				// partile charge selection
				if( fPcharge != 0){ // fPcharge : 0 all particle
					if( fPcharge==1 && ch<0 )
						continue; // 1 for + charge
					if( fPcharge==-1 && ch>0)
						continue; // -1 for - charge
				}
				Int_t label = track->GetLabel();
				AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
				itrack->SetLabel( label );
				itrack->SetParticleType( pdg);
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetCharge(ch) ;
			}
		}
	}else{
		Int_t nt = aod->GetNumberOfTracks();
		Int_t ntrack =0;
		for( int it=0; it<nt ; it++){
			AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
			if(!track){
				Error("ReadEventAOD", "Could not receive partice %d", (int) it);
				continue;
			}
			if(track->TestFilterBit( fFilterBit )){ //
				Char_t ch = (Char_t)track->Charge();
				if( fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}
				if( fPcharge !=0){ // fPcharge 0 : all particle
					if( fPcharge==1 && ch<0)
						continue; // 1 for + particle
					if( fPcharge==-1 && ch>0)
						continue; // -1 for - particle
				}
				AliJBaseTrack *itrack = new( (*TrackList)[ntrack++]) AliJBaseTrack;
				//itrack->SetID( track->GetID() );
				itrack->SetID( TrackList->GetEntriesFast() );
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetParticleType(kJHadron);
				itrack->SetCharge(track->Charge() );
				itrack->SetStatus(track->GetStatus() );
			}
		}
	} //read aod reco track done.
}
//______________________________________________________________________________
Bool_t AliJFFlucTask::IsGoodEvent( AliAODEvent *event){
	//event selection here!
	//check vertex
	AliVVertex *vtx = event->GetPrimaryVertex();
	if(!vtx || vtx->GetNContributors() <= 0)
		return kFALSE;
	double zvert = vtx->GetZ();
	if(zvert < -fzvtxCut || zvert > fzvtxCut)
		return kFALSE;

	// event cent flatting  --- do it only when IsCentFlat is true
	if(IsCentFlat == kTRUE){
		//float centrality = ReadAODCentrality( event, fCentDetName); //"V0M"
		float centrality = ReadCentrality(event,fCentDetName);
		double cent_flat_ratio = h_ratio->GetBinContent( (h_ratio->GetXaxis()->FindBin(centrality))) ;
		if (gRandom->Uniform(0, cent_flat_ratio) > 1 )
			return kFALSE;
	}

	if(IsKineOnly == kTRUE)
		return kTRUE;

	int frunNumber = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEvent()->GetRunNumber();
	if(frunNumber < 0)
		cout << "ERROR: unknown run number" << endl;
	AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
	fRunTable->SetRunNumber( frunNumber );

	int fperiod = fRunTable->GetRunNumberToPeriod(frunNumber);
	if(fperiod == AliJRunTable::kLHC15o){
		const AliVVertex* vtTrc = event->GetPrimaryVertex();
		const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
		double covTrc[6],covSPD[6];
		vtTrc->GetCovarianceMatrix(covTrc);
		vtSPD->GetCovarianceMatrix(covSPD);
		double dz = vtTrc->GetZ()-vtSPD->GetZ();
		double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
		double errTrc = TMath::Sqrt(covTrc[5]);
		double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
		if(TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
			return kFALSE;

		AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
		if(!pms){
			AliError("MultSelection unavailable.");
			return kFALSE;
		}

		double v0mcent = pms->GetMultiplicityPercentile("V0M");
		double cl0cent = pms->GetMultiplicityPercentile("CL0");
		double center = 0.973488*cl0cent+0.0157497;
		double sigma = 0.673612+cl0cent*(0.0290718+cl0cent*(-0.000546728+cl0cent*5.82749e-06));
		if(v0mcent < center-5.0*sigma || v0mcent > center+5.5*sigma || v0mcent < 0.0 || v0mcent > 60.0)
			return kFALSE;
	}

	TPCTracks = 0;
	GlobTracks = 0;
	Int_t nTracks = event->GetNumberOfTracks();
	for(int it = 0; it < nTracks; it++){
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
		//AliAODTrack* trackAOD = event->GetTrack(itracks);
		if (!trackAOD || !trackAOD->TestFilterBit(1))
			continue;
		if ((trackAOD->Pt() < 0.2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > 0.8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2) )
			continue;
		TPCTracks++;
	}

	for(int it = 0; it < nTracks; it++){
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
		if (!trackAOD || !trackAOD->TestFilterBit(16))
			continue;
		if ((trackAOD->Pt() < 0.2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > 0.8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1) )
			continue;
		Double_t b[2];
		Double_t bCov[3];
		if (!(trackAOD->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov) ))
			continue;
		if ( (TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3) )
			continue;
		GlobTracks++;
	}

	FB32Tracks = 0;
	FB32TOFTracks = 0;
	for (int it = 0; it < nTracks; it++){
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
		if (!trackAOD || !trackAOD->TestFilterBit(32))
			continue;
		FB32Tracks++;
		if (TMath::Abs(trackAOD->GetTOFsignalDz()) <= 10 && trackAOD->GetTOFsignal() >= 12000 && trackAOD->GetTOFsignal() <= 25000)
			FB32TOFTracks++;
	}

	if(fCutOutliers == kTRUE){
		if(fperiod == AliJRunTable::kLHC15o){
			AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
			if(!pms){
				AliError("MultSelection unavailable.");
				return kFALSE;
			}

			double v0mcent = pms->GetMultiplicityPercentile("V0M");
			double tfbtpc = (double)TPCTracks;
			double lcut = 2.31181837e+03+v0mcent*(-7.79946952e+01+v0mcent*(8.45194500e-01+v0mcent*(-1.72787009e-03-1.86192490e-05*v0mcent)));
			if(tfbtpc < lcut)
				return kFALSE;
			double hcut = 3.15901050e+03+v0mcent*(-9.42636072e+01+v0mcent*(8.06432447e-01+v0mcent*(3.37574557e-03-6.14272547e-05*v0mcent)));
			if(tfbtpc > hcut)
				return kFALSE;

			double tfb32 = (double)FB32Tracks;
			double tfb32tof = (double)FB32TOFTracks;
			double mu32tof = -1.0178+tfb32*(0.333132+tfb32*(9.10282e-05-1.61861e-08*tfb32));
			double sigma32tof = 1.47848+tfb32*(0.0385923+tfb32*(-5.06153e-05+tfb32*(4.37641e-08+tfb32*(-1.69082e-11+tfb32*2.35085e-15))));
			double nsigma[] = {4.0,4.0};
			if(tfb32tof < mu32tof-nsigma[0]*sigma32tof || tfb32tof > mu32tof+nsigma[1]*sigma32tof)
				return kFALSE;

		}else{
			if(!((double)TPCTracks > (-40.3+1.22*GlobTracks) && (double)TPCTracks < (32.1+1.59*GlobTracks)))
				return kFALSE;
		}
			
	}

	return kTRUE;
}
//
//______________________________________________________________________________
void AliJFFlucTask::Terminate(Option_t *)
{
	//fFFlucAna->Terminate("");
	cout<<"AliJFFlucTask Analysis DONE !!"<<endl;
}


//______________________________________________________________________________
Bool_t AliJFFlucTask::IsThisAWeakDecayingParticle(AliMCParticle *thisGuy)
{
	// In order to prevent analyzing daughters from weak decays
	// - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
	Int_t pdgcode = TMath::Abs( thisGuy->PdgCode() );
	Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
		3122, 3112, // Lambda0 Sigma+-
		130, 310 // K_L0 K_S0
	};
	for(Int_t i=0; i!=7; ++i)
		if( myWeakParticles[i] == pdgcode ){
			return kTRUE;
		}

	return kFALSE;
}
//===============================================================================
Bool_t AliJFFlucTask::IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy)
{
	// In order to prevent analyzing daughters from weak decays
	// - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
	Int_t pdgcode = TMath::Abs( thisGuy->GetPdgCode() );
	Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
		3122, 3112, // Lambda0 Sigma+-
		130, 310 // K_L0 K_S0
	};
	for(Int_t i=0; i!=7; ++i)
		if( myWeakParticles[i] == pdgcode ) {
			return kTRUE;
		}
	return kFALSE;
}
//______________________________________________________________________________
void AliJFFlucTask::SetEffConfig( int effMode, int FilterBit)
{
	fEffMode = effMode;
	fEffFilterBit = 0; // as default
	if( FilterBit == 128 )
		fEffFilterBit = 0;
	if( FilterBit == 768 )
		fEffFilterBit = 5;
	cout << "setting to EffCorr Mode : " << effMode << endl;
	cout << "setting to EffCorr Filter bit : " << FilterBit  << " = " << fEffFilterBit << endl;
}
//______________________________________________________________________________
void AliJFFlucTask::SetIsCentFlat( Bool_t isCentFlat ){
	cout << "Setting to flatting Centrality with LHC11h data : " << isCentFlat  << endl;
	IsCentFlat = isCentFlat;

	if(IsCentFlat == kTRUE){
		if( IsMC == kTRUE )
			Printf("this is not LHC11h data!!!!");
		// ratio from ExtractCentRatio.C // do not change this manually
		h_ratio = new TH1D("h_ratio","",240,0,60);
		h_ratio->SetBinContent(1,1.04895);
		h_ratio->SetBinContent(2,1);
		h_ratio->SetBinContent(3,1.04325);
		h_ratio->SetBinContent(4,1.03867);
		h_ratio->SetBinContent(5,1.01948);
		h_ratio->SetBinContent(6,1.08264);
		h_ratio->SetBinContent(7,1.07225);
		h_ratio->SetBinContent(8,1.09455);
		h_ratio->SetBinContent(9,1.09914);
		h_ratio->SetBinContent(10,1.07226);
		h_ratio->SetBinContent(11,1.07857);
		h_ratio->SetBinContent(12,1.07964);
		h_ratio->SetBinContent(13,1.08201);
		h_ratio->SetBinContent(14,1.10538);
		h_ratio->SetBinContent(15,1.06169);
		h_ratio->SetBinContent(16,1.10798);
		h_ratio->SetBinContent(17,1.08845);
		h_ratio->SetBinContent(18,1.08419);
		h_ratio->SetBinContent(19,1.09068);
		h_ratio->SetBinContent(20,1.07633);
		h_ratio->SetBinContent(21,1.36282);
		h_ratio->SetBinContent(22,1.35824);
		h_ratio->SetBinContent(23,1.39488);
		h_ratio->SetBinContent(24,1.34785);
		h_ratio->SetBinContent(25,1.41281);
		h_ratio->SetBinContent(26,1.36358);
		h_ratio->SetBinContent(27,1.35688);
		h_ratio->SetBinContent(28,1.37658);
		h_ratio->SetBinContent(29,1.31628);
		h_ratio->SetBinContent(30,1.38905);
		h_ratio->SetBinContent(31,1.36245);
		h_ratio->SetBinContent(32,1.3199);
		h_ratio->SetBinContent(33,1.31438);
		h_ratio->SetBinContent(34,1.29465);
		h_ratio->SetBinContent(35,1.32072);
		h_ratio->SetBinContent(36,1.24303);
		h_ratio->SetBinContent(37,1.18857);
		h_ratio->SetBinContent(38,1.11082);
		h_ratio->SetBinContent(39,1.06849);
		h_ratio->SetBinContent(40,1);
		h_ratio->SetBinContent(41,3.33011);
		h_ratio->SetBinContent(42,2.9697);
		h_ratio->SetBinContent(43,2.66142);
		h_ratio->SetBinContent(44,2.24475);
		h_ratio->SetBinContent(45,2.00121);
		h_ratio->SetBinContent(46,1.67406);
		h_ratio->SetBinContent(47,1.43596);
		h_ratio->SetBinContent(48,1.29798);
		h_ratio->SetBinContent(49,1.19522);
		h_ratio->SetBinContent(50,1.11261);
		h_ratio->SetBinContent(51,1.04576);
		h_ratio->SetBinContent(52,1.04513);
		h_ratio->SetBinContent(53,1.07547);
		h_ratio->SetBinContent(54,1.04186);
		h_ratio->SetBinContent(55,1.06598);
		h_ratio->SetBinContent(56,1.03554);
		h_ratio->SetBinContent(57,1.06338);
		h_ratio->SetBinContent(58,1.051);
		h_ratio->SetBinContent(59,1.03166);
		h_ratio->SetBinContent(60,1.06029);
		h_ratio->SetBinContent(61,1.01288);
		h_ratio->SetBinContent(62,1);
		h_ratio->SetBinContent(63,1.0638);
		h_ratio->SetBinContent(64,1.04379);
		h_ratio->SetBinContent(65,1.04036);
		h_ratio->SetBinContent(66,1.06243);
		h_ratio->SetBinContent(67,1.0516);
		h_ratio->SetBinContent(68,1.04946);
		h_ratio->SetBinContent(69,1.03621);
		h_ratio->SetBinContent(70,1.06466);
		h_ratio->SetBinContent(71,1.01806);
		h_ratio->SetBinContent(72,1.04632);
		h_ratio->SetBinContent(73,1.01945);
		h_ratio->SetBinContent(74,1.03469);
		h_ratio->SetBinContent(75,1.04348);
		h_ratio->SetBinContent(76,1.02205);
		h_ratio->SetBinContent(77,1.06899);
		h_ratio->SetBinContent(78,1.04297);
		h_ratio->SetBinContent(79,1.05137);
		h_ratio->SetBinContent(80,1.07058);
		h_ratio->SetBinContent(81,1.05674);
		h_ratio->SetBinContent(82,1.01647);
		h_ratio->SetBinContent(83,1.03579);
		h_ratio->SetBinContent(84,1.0688);
		h_ratio->SetBinContent(85,1.07705);
		h_ratio->SetBinContent(86,1.05545);
		h_ratio->SetBinContent(87,1.05272);
		h_ratio->SetBinContent(88,1.06507);
		h_ratio->SetBinContent(89,1.027);
		h_ratio->SetBinContent(90,1.0703);
		h_ratio->SetBinContent(91,1.06285);
		h_ratio->SetBinContent(92,1.05242);
		h_ratio->SetBinContent(93,1.05984);
		h_ratio->SetBinContent(94,1.03455);
		h_ratio->SetBinContent(95,1.07065);
		h_ratio->SetBinContent(96,1.08061);
		h_ratio->SetBinContent(97,1.05791);
		h_ratio->SetBinContent(98,1.01024);
		h_ratio->SetBinContent(99,1.03655);
		h_ratio->SetBinContent(100,1.0229);
		h_ratio->SetBinContent(101,1.04371);
		h_ratio->SetBinContent(102,1.02322);
		h_ratio->SetBinContent(103,1.00899);
		h_ratio->SetBinContent(104,1.04831);
		h_ratio->SetBinContent(105,1.08);
		h_ratio->SetBinContent(106,1.04502);
		h_ratio->SetBinContent(107,1.03318);
		h_ratio->SetBinContent(108,1.04894);
		h_ratio->SetBinContent(109,1.01745);
		h_ratio->SetBinContent(110,1.05296);
		h_ratio->SetBinContent(111,1.07492);
		h_ratio->SetBinContent(112,1.08634);
		h_ratio->SetBinContent(113,1.02578);
		h_ratio->SetBinContent(114,1.00636);
		h_ratio->SetBinContent(115,1.00809);
		h_ratio->SetBinContent(116,1.03827);
		h_ratio->SetBinContent(117,1.01522);
		h_ratio->SetBinContent(118,1);
		h_ratio->SetBinContent(119,1.04044);
		h_ratio->SetBinContent(120,1.02615);
		h_ratio->SetBinContent(121,1.07061);
		h_ratio->SetBinContent(122,1.03836);
		h_ratio->SetBinContent(123,1.0654);
		h_ratio->SetBinContent(124,1.06016);
		h_ratio->SetBinContent(125,1.02215);
		h_ratio->SetBinContent(126,1.04882);
		h_ratio->SetBinContent(127,1.03332);
		h_ratio->SetBinContent(128,1.02424);
		h_ratio->SetBinContent(129,1);
		h_ratio->SetBinContent(130,1.02053);
		h_ratio->SetBinContent(131,1.04855);
		h_ratio->SetBinContent(132,1.034);
		h_ratio->SetBinContent(133,1.0276);
		h_ratio->SetBinContent(134,1.01399);
		h_ratio->SetBinContent(135,1.02446);
		h_ratio->SetBinContent(136,1.03738);
		h_ratio->SetBinContent(137,1.0156);
		h_ratio->SetBinContent(138,1.06815);
		h_ratio->SetBinContent(139,1.06946);
		h_ratio->SetBinContent(140,1.03094);
		h_ratio->SetBinContent(141,1.01122);
		h_ratio->SetBinContent(142,1.04319);
		h_ratio->SetBinContent(143,1.03308);
		h_ratio->SetBinContent(144,1.01178);
		h_ratio->SetBinContent(145,1.04965);
		h_ratio->SetBinContent(146,1.04205);
		h_ratio->SetBinContent(147,1.01552);
		h_ratio->SetBinContent(148,1.07672);
		h_ratio->SetBinContent(149,1.03083);
		h_ratio->SetBinContent(150,1.00266);
		h_ratio->SetBinContent(151,1.04081);
		h_ratio->SetBinContent(152,1.05451);
		h_ratio->SetBinContent(153,1.05444);
		h_ratio->SetBinContent(154,1.03539);
		h_ratio->SetBinContent(155,1.06241);
		h_ratio->SetBinContent(156,1.04532);
		h_ratio->SetBinContent(157,1.03856);
		h_ratio->SetBinContent(158,1.06928);
		h_ratio->SetBinContent(159,1.02747);
		h_ratio->SetBinContent(160,1.01488);
		h_ratio->SetBinContent(161,1.08372);
		h_ratio->SetBinContent(162,1.03469);
		h_ratio->SetBinContent(163,1.06914);
		h_ratio->SetBinContent(164,1.05105);
		h_ratio->SetBinContent(165,1.05741);
		h_ratio->SetBinContent(166,1.06204);
		h_ratio->SetBinContent(167,1.04757);
		h_ratio->SetBinContent(168,1.061);
		h_ratio->SetBinContent(169,1.03418);
		h_ratio->SetBinContent(170,1.07949);
		h_ratio->SetBinContent(171,1.01637);
		h_ratio->SetBinContent(172,1.06089);
		h_ratio->SetBinContent(173,1.04924);
		h_ratio->SetBinContent(174,1.04922);
		h_ratio->SetBinContent(175,1.0725);
		h_ratio->SetBinContent(176,1.09398);
		h_ratio->SetBinContent(177,1.03688);
		h_ratio->SetBinContent(178,1.08006);
		h_ratio->SetBinContent(179,1.01428);
		h_ratio->SetBinContent(180,1.07294);
		h_ratio->SetBinContent(181,1.06006);
		h_ratio->SetBinContent(182,1.06782);
		h_ratio->SetBinContent(183,1.09687);
		h_ratio->SetBinContent(184,1.07964);
		h_ratio->SetBinContent(185,1.07039);
		h_ratio->SetBinContent(186,1.09371);
		h_ratio->SetBinContent(187,1.01951);
		h_ratio->SetBinContent(188,1.03736);
		h_ratio->SetBinContent(189,1.08374);
		h_ratio->SetBinContent(190,1.04944);
		h_ratio->SetBinContent(191,1.03483);
		h_ratio->SetBinContent(192,1.09588);
		h_ratio->SetBinContent(193,1.04325);
		h_ratio->SetBinContent(194,1);
		h_ratio->SetBinContent(195,1.07037);
		h_ratio->SetBinContent(196,1.07313);
		h_ratio->SetBinContent(197,1.11307);
		h_ratio->SetBinContent(198,1.03543);
		h_ratio->SetBinContent(199,1.05078);
		h_ratio->SetBinContent(200,1.02593);
		h_ratio->SetBinContent(201,14.7486);
		h_ratio->SetBinContent(202,14.6304);
		h_ratio->SetBinContent(203,15.05);
		h_ratio->SetBinContent(204,14.4104);
		h_ratio->SetBinContent(205,13.8821);
		h_ratio->SetBinContent(206,14.2444);
		h_ratio->SetBinContent(207,13.7572);
		h_ratio->SetBinContent(208,13.4553);
		h_ratio->SetBinContent(209,13.1558);
		h_ratio->SetBinContent(210,12.4867);
		h_ratio->SetBinContent(211,12.123);
		h_ratio->SetBinContent(212,12.2009);
		h_ratio->SetBinContent(213,11.1698);
		h_ratio->SetBinContent(214,10.6532);
		h_ratio->SetBinContent(215,9.40398);
		h_ratio->SetBinContent(216,8.81708);
		h_ratio->SetBinContent(217,7.87832);
		h_ratio->SetBinContent(218,6.94789);
		h_ratio->SetBinContent(219,5.71367);
		h_ratio->SetBinContent(220,4.74537);
		h_ratio->SetBinContent(221,4.04486);
		h_ratio->SetBinContent(222,3.01101);
		h_ratio->SetBinContent(223,2.35428);
		h_ratio->SetBinContent(224,1.81681);
		h_ratio->SetBinContent(225,1.4577);
		h_ratio->SetBinContent(226,1.13349);
		h_ratio->SetBinContent(227,1.14531);
		h_ratio->SetBinContent(228,1.07198);
		h_ratio->SetBinContent(229,1.03062);
		h_ratio->SetBinContent(230,1.09616);
		h_ratio->SetBinContent(231,1.00376);
		h_ratio->SetBinContent(232,1.01719);
		h_ratio->SetBinContent(233,1.02229);
		h_ratio->SetBinContent(234,1);
		h_ratio->SetBinContent(235,1.02041);
		h_ratio->SetBinContent(236,1.05855);
		h_ratio->SetBinContent(237,1.04459);
		h_ratio->SetBinContent(238,1.03411);
		h_ratio->SetBinContent(239,1.04969);
		h_ratio->SetBinContent(240,1.03707);
		//
	}
}
//
//

void AliJFFlucTask::ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList)
{
	Int_t nt = mcEvent->GetNumberOfPrimaries();
	Int_t ntrack =0;
	for (Int_t it = 0; it < nt; it++) {
		AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(it));
		if(mcEvent->IsPhysicalPrimary(it)) {
			double Pt = track->Pt();
			if( Pt < fPt_min || Pt > fPt_max )
				continue ; // pt cut
			TParticle *particle = track->Particle();
			if(!particle)
				continue;

			if( IsExcludeWeakDecay == kTRUE){
				Int_t gMotherIndex = particle->GetFirstMother(); //
				if(gMotherIndex != -1){
					DEBUG( 4,  "this particle has a mother " );
					AliMCParticle* motherParticle= dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(gMotherIndex));
					if(motherParticle){
						if(IsThisAWeakDecayingParticle( motherParticle)){
							DEBUG ( 4, Form("this particle will be removed because it comes from : %d", motherParticle->PdgCode() ));
							continue;
						}
					}
				}
			}

			Int_t pdg = particle->GetPdgCode();
			Char_t ch = (Char_t) track->Charge();
			if(fPcharge != 0){ // fPcharge 0 : all particle
				if(fPcharge == 1 && ch < 0)
					continue; // 1 for + particle
				if(fPcharge == -1 && ch > 0)
					continue; // -1 for - particle
			}
			Int_t label = track->GetLabel();
			AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
			itrack->SetLabel( label );
			itrack->SetParticleType( pdg);
			itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
			itrack->SetCharge(ch) ;
		}
	}
}

double AliJFFlucTask::GetCentralityFromImpactPar(double ip) {
/*
\begin{tabular}{ |c|c c c c c c c c| }
 \hline
Centrality(\%)&0-5      &5-10     &10-20    &20-30    &30-40     &40-50      &50-60      &60-70\\
\hline
b(fm) AMPT    &0.00-3.72&3.72-5.23&5.23-7.31&7.31-8.88&8.88-10.20&10.20-11.38&11.38-12.47&12.47-13.50\\
b(fm) HIJING  &0.00-3.60&3.60-5.09&5.09-7.20&7.20-8.83&8.83-10.20&10.20-11.40&11.40-12.49&12.49-13.49\\
b(fm) ALICE   &0.00-3.50&3.50-4.94&4.94-6.98&6.98-    &    -9.88 &9.81-      &     -12.09&12.09-\\
 \hline
\end{tabular}
   \begin{tablenotes}
 \url{https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies}
    \end{tablenotes}
*/
	static double bmin[10]={0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
	static double centmean[10]={2.5,7.5,15,25,35,45,55,65,75,90};
	for(int i=0;i<9;i++){
		if(bmin[i] < ip && ip <= bmin[i+1])
			return centmean[i];
	}
	return 0.0;
}



void AliJFFlucTask::SetInFileName( TString inName ){
	if( IsPhiModule == kFALSE){
		cout << "Phi Modulation option is setted OFF : InFile will be ignored" << endl;
		return ;
	}
	cout << "Setting InFIle as " << inName.Data() << endl;
	fInFileName = inName;
	TGrid::Connect("alien:");
	TFile *inclusFile = TFile::Open( fInFileName.Data() , "read" );
	cout << "File connected to Alien" << endl;
	for(int icent=0; icent< 7; icent++){
		for(int isub=0; isub<2; isub++){
			h_ModuledPhi[icent][isub] = dynamic_cast<TH1D *>(inclusFile->Get(Form("h_phi_moduleC%02dS%02d", icent, isub)));
		}
	}
}
