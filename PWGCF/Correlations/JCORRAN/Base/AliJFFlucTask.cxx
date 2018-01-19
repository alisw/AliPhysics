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
	flags = 0;

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
AliJFFlucTask::AliJFFlucTask(const char *name):
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
	fCentDetName="V0M";
	flags = 0;

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
	if(flags & FLUC_SCPT)
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_SCPT);
	if(flags & FLUC_EBE_WEIGHTING)
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING);
	if(flags & FLUC_PHI_MODULATION){
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_PHI_MODULATION);
		//setting histos for phi modulation
		for(int icent=0; icent<7; icent++){
			for(int isub=0; isub<2; isub++){
				fFFlucAna->SetPhiModuleHistos( icent, isub, h_ModuledPhi[icent][isub] );
			}
		}

		if(flags & FLUC_PHI_INVERSE)
			fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_PHI_INVERSE);
		//if(flags & FLUC_PHI_REJECTION);
			//fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_PHI_REJECTION);
	}
	//fFFlucAna->SetIsPhiModule( IsPhiModule);
	//fFFlucAna->SetIsSCptdep( IsSCptdep ) ;
	//fFFlucAna->SetSCwithQC( IsSCwithQC );
	//fFFlucAna->SetEbEWeight( IsEbEWeighted );

	fFFlucAna->SetQCEtaCut( fQC_eta_min, fQC_eta_max, 0.5 );

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
	if(flags & FLUC_KINEONLY) {
		AliMCEvent *mcEvent = MCEvent();
		if (!mcEvent) {
			AliError("ERROR: mcEvent not available");
			return;
		}

		AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
		if(!headerH)
			return;
		//Double_t gReactionPlane = headerH->ReactionPlaneAngle();
		Double_t gImpactParameter = headerH->ImpactParameter();
		fCent = GetCentralityFromImpactPar(gImpactParameter);
		if(flags & FLUC_ALICE_IPINFO){
			//force to use ALICE impact parameter setting
			double ALICE_Cent[8] = {0, 5, 10, 20, 30, 40, 50, 60};
			double ALICE_IPinfo[8] = {0, 3.50, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09};
			for(int icent=0; icent<8; icent++){
				if(fImpactParameter >= ALICE_IPinfo[icent] && fImpactParameter < ALICE_IPinfo[icent+1])
					fCent = 0.5f*(ALICE_Cent[icent]+ALICE_Cent[icent+1]);
			}
		}
		if( fEvtNum == 1 ){
			int runN = 1234;
			fFFlucAna->GetAliJEfficiency()->SetRunNumber ( runN );
			fFFlucAna->GetAliJEfficiency()->Load();
		}
		ReadKineTracks( mcEvent, fInputList, fCent ) ; // read tracklist
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

		if(!IsGoodEvent( currentEvent ))
			return;
		ReadAODTracks( currentEvent, fInputList, fCent ) ; // read tracklist
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
void AliJFFlucTask::ReadAODTracks(AliAODEvent *aod, TClonesArray *TrackList, float fCent)
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

				if(fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}

				if(flags & FLUC_PHI_REJECTION){
					int isub = (int)(track->Eta() > 0.0);
					int cbin = AliJFFlucAnalysis::GetCentralityClass(fCent);
					int pbin = h_ModuledPhi[cbin][isub]->GetXaxis()->FindBin(TMath::Pi()-track->Phi());
					if(gRandom->Uniform(0,1) > h_ModuledPhi[cbin][isub]->GetBinContent(pbin)/h_ModuledPhi[cbin][isub]->GetMaximum())
						continue;
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

				if(flags & FLUC_PHI_REJECTION){
					int isub = (int)(track->Eta() > 0.0);
					int cbin = AliJFFlucAnalysis::GetCentralityClass(fCent);
					int pbin = h_ModuledPhi[cbin][isub]->GetXaxis()->FindBin(TMath::Pi()-track->Phi());
					if(gRandom->Uniform(0,1) > h_ModuledPhi[cbin][isub]->GetBinContent(pbin)/h_ModuledPhi[cbin][isub]->GetMaximum())
						continue;
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
	if(flags & FLUC_CENT_FLATTENING){
		//float centrality = ReadAODCentrality( event, fCentDetName); //"V0M"
		float centrality = ReadCentrality(event,fCentDetName);
		double cent_flat_ratio = h_ratio->GetBinContent( (h_ratio->GetXaxis()->FindBin(centrality))) ;
		if (gRandom->Uniform(0, cent_flat_ratio) > 1 )
			return kFALSE;
	}

	if(flags & FLUC_KINEONLY)
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
void AliJFFlucTask::ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList, float fCent)
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
			
			if(flags & FLUC_PHI_REJECTION){
				int isub = (int)(track->Eta() > 0.0);
				int cbin = AliJFFlucAnalysis::GetCentralityClass(fCent);
				int pbin = h_ModuledPhi[cbin][isub]->GetXaxis()->FindBin(TMath::Pi()-track->Phi());
				if(gRandom->Uniform(0,1) > h_ModuledPhi[cbin][isub]->GetBinContent(pbin)/h_ModuledPhi[cbin][isub]->GetMaximum())
					continue;
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

void AliJFFlucTask::EnableCentFlat(const TString fname){
	//cout << "Setting to flatting Centrality with LHC11h data : " << isCentFlat  << endl;
	//flags |= FLUC_CENT_FLATTENING;
	h_ratio = 0; //TODO: rename, add loading if we ever need this again.
}

void AliJFFlucTask::EnablePhiModule(const TString fname){
	//flags |= FLUC_PHI_MODULATION;
	cout << "Setting InFIle as " << fname.Data() << endl;
	TGrid::Connect("alien:");
	TFile *inclusFile = TFile::Open( fname.Data() , "read" );
	for(int icent=0; icent< 7; icent++){
		for(int isub=0; isub<2; isub++){
			h_ModuledPhi[icent][isub] = dynamic_cast<TH1D *>(inclusFile->Get(Form("h_phi_moduleC%02dS%02d", icent, isub)));
		}
	}
}

