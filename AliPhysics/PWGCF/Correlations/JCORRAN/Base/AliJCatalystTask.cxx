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
// author: Jasper Parkila, D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliGenHijingEventHeader.h"
#include "AliJCatalystTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJEventHeader.h"
#include "AliJHistManager.h"
#include "AliInputEventHandler.h"
#include "AliJEfficiency.h"
#include "AliMultSelection.h"
#include "AliJRunTable.h"
//#pragma GCC diagnostic warning "-Wall"
//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask():
	AliAnalysisTaskSE(),
	fRunNum(-1),
	fcent(-999),
	fZvert(-999),
	fnoCentBin(false),
	fInputList(0),
	fInputListALICE(0),
	fOutput(0)
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
	flags = 0;

	fzvtxCut = 10;  // default z vertex cut
	fJCatalystEntry = 0;
	//  DefineOutput(1, TDirectory::Class());
}

//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const char *name):
	AliAnalysisTaskSE(name),
	fRunNum(-1),
	fcent(-999),
	fZvert(-999),
	fnoCentBin(false),
	fInputList(0),
	fInputListALICE(0),
	fOutput(0)
{
	DefineOutput(1, TDirectory::Class());
	fTaskName = name;

	fEvtNum=0;
	fFilterBit = 0;
	fEta_min = 0;
	fEta_max = 0;
	fDebugLevel = 0;
	fEffFilterBit=0;
	fPt_min=0;
	fPt_max=0;
	fCentDetName="V0M";
	flags = 0;

	fzvtxCut = 10;
	fJCatalystEntry = 0;
}

//____________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const AliJCatalystTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fRunNum(ap.fRunNum),
	fcent(ap.fcent),
	fZvert(ap.fZvert),
	fnoCentBin(ap.fnoCentBin),
	fInputList(ap.fInputList),
	fInputListALICE(ap.fInputListALICE),
	fOutput(ap.fOutput)
{
	AliInfo("----DEBUG AliJCatalystTask COPY ----");
}

//_____________________________________________________________________________
AliJCatalystTask& AliJCatalystTask::operator = (const AliJCatalystTask& ap)
{
	// assignment operator
	AliInfo("----DEBUG AliJCatalystTask operator= ----");
	this->~AliJCatalystTask();
	new(this) AliJCatalystTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJCatalystTask::~AliJCatalystTask()
{
	delete fInputList;
	delete fInputListALICE;
	delete fOutput;
}

//________________________________________________________________________
void AliJCatalystTask::UserCreateOutputObjects()
{
	fInputList = new TClonesArray("AliJBaseTrack" , 2500);
	fInputList->SetOwner(kTRUE);
	fInputListALICE = new TClonesArray("AliJBaseTrack" , 2500);
	fInputListALICE->SetOwner(kTRUE);

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();
	PostData(1, fOutput);
}

//______________________________________________________________________________
void AliJCatalystTask::UserExec(Option_t* /*option*/)
{
	// Processing of one event
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));
	// initializing variables from last event
        fJCatalystEntry = fEntry;
	fInputList->Clear();
	fInputListALICE->Clear();

	float fImpactParameter = -1.0f;
	double fvertex[3];

	fEvtNum++;
	if(fEvtNum % 1000 == 0)
		cout << "evt : " << fEvtNum <<endl;

	// load current event and save track, event info
	if(flags & FLUC_KINEONLY) {
		AliMCEvent *mcEvent = MCEvent();
		if (!mcEvent) {
			AliError("ERROR: mcEvent not available");
			return;
		}

		if(!fnoCentBin) {
			AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
			if(!headerH)
				return;
			//Double_t gReactionPlane = headerH->ReactionPlaneAngle();
			Double_t gImpactParameter = headerH->ImpactParameter();
			fcent = GetCentralityFromImpactPar(gImpactParameter);
		}
		if(flags & FLUC_ALICE_IPINFO){
			//force to use ALICE impact parameter setting
			double ALICE_Cent[8] = {0, 5, 10, 20, 30, 40, 50, 60};
			double ALICE_IPinfo[8] = {0, 3.50, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09};
			for(int icent=0; icent<8; icent++){
				if(fImpactParameter >= ALICE_IPinfo[icent] && fImpactParameter < ALICE_IPinfo[icent+1])
					fcent = 0.5f*(ALICE_Cent[icent]+ALICE_Cent[icent+1]);
			}
		}
		if(fnoCentBin) fcent = 1.0; // forcing no centrality selection
		ReadKineTracks( mcEvent, fInputList, fInputListALICE, fcent ) ; // read tracklist
		AliGenEventHeader *header = mcEvent->GenEventHeader();
		if(!header)
			return;
		TArrayF gVertexArray;
		header->PrimaryVertex(gVertexArray);
		for(int i = 0; i < 3; i++)
			fvertex[i] = gVertexArray.At(i);
	} else { // Kine
		AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
		if(fnoCentBin) {
			fcent = 1.0;
		} else {
			fcent = ReadCentrality(currentEvent,fCentDetName);
		}
		fRunNum = currentEvent->GetRunNumber();

		if(!IsGoodEvent( currentEvent ))
			return;
		ReadAODTracks( currentEvent, fInputList, fcent ) ; // read tracklist
		ReadVertexInfo( currentEvent, fvertex); // read vertex info
	} // AOD
	fZvert = fvertex[2];
}

//______________________________________________________________________________
void AliJCatalystTask::ReadVertexInfo ( AliAODEvent *aod, double*  fvtx )
{
	fvtx[0] = aod->GetPrimaryVertex()->GetX();
	fvtx[1] = aod->GetPrimaryVertex()->GetY();
	fvtx[2] = aod->GetPrimaryVertex()->GetZ();
}
//______________________________________________________________________________
float AliJCatalystTask::ReadCentrality( AliAODEvent *aod, TString Trig ){
	AliMultSelection *pms = (AliMultSelection*)aod->FindListObject("MultSelection");
	return pms?pms->GetMultiplicityPercentile(Trig.Data()):aod->GetCentrality()->GetCentralityPercentile(Trig.Data());
}

//______________________________________________________________________________

void AliJCatalystTask::Init()
{
	AliInfo("Doing initialization") ;
}
//______________________________________________________________________________
void AliJCatalystTask::ReadAODTracks(AliAODEvent *aod, TClonesArray *TrackList, float fcent)
{
	//aod->Print();
	if(flags & FLUC_MC){  // how to get a flag to check  MC or not !
		TClonesArray *mcArray = (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray){ Printf("Error not a proper MC event"); };  // check mc array

		Int_t nt = mcArray->GetEntriesFast();
		Int_t ntrack = 0;
		for( int it=0; it < nt ; it++){
			AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
			if(!track) {
				Error("ReadEventAODMC","Could not read particle %d",it);
				continue;
			}
			if( track->IsPhysicalPrimary() ){
				// insert AMTP weak decay switch here
				if(flags & FLUC_EXCLUDEWDECAY){
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
				Char_t ch = (Char_t) track->Charge();
				if (ch == 0 ) continue; // Remove all others than charged particles

				if(fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}

				Int_t pdg = track->GetPdgCode();
				// partile charge selection
				if( fPcharge != 0){ // fPcharge : 0 all particle
					if( fPcharge==1 && ch<0 )
						continue; // 1 for + charge
					if( fPcharge==-1 && ch>0)
						continue; // -1 for - charge
				}
				if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
				Int_t label = track->GetLabel();
				AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
				itrack->SetLabel( label );
				itrack->SetParticleType(pdg);
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
				Error("ReadEventAOD", "Could not read particle %d", (int) it);
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
				if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
				AliJBaseTrack *itrack = new( (*TrackList)[ntrack++]) AliJBaseTrack;
				itrack->SetID( TrackList->GetEntriesFast() );
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetParticleType(kJHadron);
				itrack->SetCharge(track->Charge() );
				itrack->SetStatus(track->GetStatus() );
			}
		}
	} //read aod reco track done.
	if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
    
}
//______________________________________________________________________________
Bool_t AliJCatalystTask::IsGoodEvent( AliAODEvent *event){
	//event selection here!
	//check vertex
	AliVVertex *vtx = event->GetPrimaryVertex();
	if(!vtx || vtx->GetNContributors() <= 0)
		return kFALSE;
	double zvert = vtx->GetZ();
	if(zvert < -fzvtxCut || zvert > fzvtxCut)
		return kFALSE;

	if(flags & FLUC_KINEONLY)
		return kTRUE;

	//int frunNumber = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEvent()->GetRunNumber();
	AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
	fRunTable->SetRunNumber(fRunNum);

	int fperiod = fRunTable->GetRunNumberToPeriod(fRunNum);
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

	if(flags & FLUC_CUT_OUTLIERS){
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
void AliJCatalystTask::Terminate(Option_t *)
{
	//fFFlucAna->Terminate("");
	cout<<"AliJCatalystTask Analysis DONE !!"<<endl;
}


//______________________________________________________________________________
Bool_t AliJCatalystTask::IsThisAWeakDecayingParticle(AliMCParticle *thisGuy)
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
Bool_t AliJCatalystTask::IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy)
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
void AliJCatalystTask::SetEffConfig( int effMode, int FilterBit)
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
void AliJCatalystTask::ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fcent)
{
	Int_t nt = mcEvent->GetNumberOfPrimaries();
	Int_t ntrack = 0;
	for (Int_t it = 0; it < nt; it++) {
		AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(it));
		if(mcEvent->IsPhysicalPrimary(it)) {
			double Pt = track->Pt();
			if( Pt < fPt_min || Pt > fPt_max )
				continue ; // pt cut
			TParticle *particle = track->Particle();
			if(!particle)
				continue;

			if(flags & FLUC_EXCLUDEWDECAY){
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
			if(ch==0) continue; // Remove all others than charged particles
			if(fPcharge != 0){ // fPcharge 0 : all particle
				if(fPcharge == 1 && ch < 0)
					continue; // 1 for + particle
				if(fPcharge == -1 && ch > 0)
					continue; // -1 for - particle
			}
			if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
			Int_t label = track->GetLabel();
			AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
			itrack->SetLabel( label );
			itrack->SetParticleType(pdg);
			itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
			itrack->SetCharge(ch) ;
			if(TMath::Abs(track->Eta()) < 0.8) {
				AliJBaseTrack *jtrack =  new ((*TrackListALICE)[TrackListALICE->GetEntriesFast()])AliJBaseTrack;
				jtrack->SetLabel( label );
				jtrack->SetParticleType(pdg);
				jtrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				jtrack->SetCharge(ch) ;
			}
		}
	}
}

double AliJCatalystTask::GetCentralityFromImpactPar(double ip) {
	//https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
	static double bmin[12] = {0.0,1.60,2.27,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
	static double centmean[12] = {0.5,1.5,3.5,7.5,15,25,35,45,55,65,75,90};
	for(UInt_t i = 0; i < 11; i++){
		if(bmin[i+1] > ip)
			return centmean[i];
	}
	return 0.0;
}

