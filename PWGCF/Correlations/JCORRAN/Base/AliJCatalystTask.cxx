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

#include <TGrid.h>
#include <TRandom.h>
#include <AliAnalysisUtils.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODMCParticle.h>
#include <AliMCEvent.h>
#include <AliGenHijingEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliGenHepMCEventHeader.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAODEvent.h>
#include <AliPWG0Helper.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
//#include <AliInputEventHandler.h>
#include <AliMultSelection.h>
#include "AliJCatalystTask.h"
#include "AliJTrack.h"
#include "AliJHistManager.h"
#include "AliJRunTable.h"
#include "AliJFFlucAnalysis.h" // TEMP for getting bins

//#pragma GCC diagnostic warning "-Wall"
//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask():
	AliAnalysisTaskSE(),
	fInputList(0),
	fInputListALICE(0),
	//fOutput(0),
	fCentDetName("V0M"),
	paodEvent(0),
	fcent(-999),
	fZvert(-999),
	fnoCentBin(false),
	fDebugLevel(0),
	fEvtNum(0),
	fFilterBit(0),
	fNumTPCClusters(70),
	fEffMode(0),
	fEffFilterBit(0),
	fRunNum(-1),
	fPcharge(0),
	fEta_min(-0.8),
	fEta_max(0.8),
	fPt_min(0.2),
	fPt_max(5.0),
	fzvtxCut(10.0),
	fremovebadarea(kFALSE),
	flags(0),
	fJCatalystEntry(0),
	fIsGoodEvent(false),
	fJCorMapTask(NULL),
	fJCorMapTaskName("JCorrectionMapTask"),
	pPhiWeights(0),
	grEffCor(0),
	fCentBinEff(0),
	phiMapIndex(0),
// QA part.
	fMainList(NULL),
	bSaveAllQA(kFALSE),
	bSaveHMOhist(kFALSE),
	fCentralityBins(16),
	fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.)
{
	//
}

//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const char *name):
	AliAnalysisTaskSE(name),
	fInputList(0),
	fInputListALICE(0),
	//fOutput(0),
	fTaskName(name),
	fCentDetName("V0M"),
	paodEvent(0),
	fcent(-999),
	fZvert(-999),
	fnoCentBin(false),
	fDebugLevel(0),
	fEvtNum(0),
	fFilterBit(0),
	fNumTPCClusters(70),
	fEffMode(0),
	fEffFilterBit(0),
	fRunNum(-1),
	fPcharge(0),
	fEta_min(-0.8),
	fEta_max(0.8),
	fPt_min(0.2),
	fPt_max(5.0),
	fzvtxCut(10.0),
	fremovebadarea(kFALSE),
	flags(0),
	fJCatalystEntry(0),
	fIsGoodEvent(false),
	fJCorMapTask(NULL),
	fJCorMapTaskName("JCorrectionMapTask"),
	pPhiWeights(0),
	grEffCor(0),
	fCentBinEff(0),
	phiMapIndex(0),
// QA part.
	fMainList(NULL),
	bSaveAllQA(kFALSE),
	bSaveHMOhist(kFALSE),
	fCentralityBins(16),
	fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
 fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.)
{
// Main list to save the output of the QA.
  fMainList = new TList();
  fMainList->SetName("fJCatalystOutput");
  fMainList->SetOwner(kTRUE);

	InitializeArrays();

	//DefineOutput(1, TDirectory::Class());	// Uncommented in the current version on AliPhysics.
	DefineOutput(1, TList::Class());
}

//____________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const AliJCatalystTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fInputList(ap.fInputList),
	fInputListALICE(ap.fInputListALICE),
	//fOutput(ap.fOutput),
	fcent(ap.fcent),
	fZvert(ap.fZvert),
	fnoCentBin(ap.fnoCentBin),
	fRunNum(ap.fRunNum),
	fJCorMapTask(ap.fJCorMapTask),
	fJCorMapTaskName(ap.fJCorMapTaskName),
	pPhiWeights(ap.pPhiWeights),
	grEffCor(ap.grEffCor),
	fCentBinEff(ap.fCentBinEff),
	phiMapIndex(ap.phiMapIndex),
// QA part.
	fMainList(ap.fMainList),
	bSaveAllQA(ap.bSaveAllQA),
	bSaveHMOhist(ap.bSaveHMOhist),
	fCentralityBins(ap.fCentralityBins),
	fcent_0(ap.fcent_0), fcent_1(ap.fcent_1), fcent_2(ap.fcent_2), fcent_3(ap.fcent_3), fcent_4(ap.fcent_4), fcent_5(ap.fcent_5), fcent_6(ap.fcent_6), fcent_7(ap.fcent_7), fcent_8(ap.fcent_8), fcent_9(ap.fcent_9), fcent_10(ap.fcent_10), fcent_11(ap.fcent_11), fcent_12(ap.fcent_12), fcent_13(ap.fcent_13), fcent_14(ap.fcent_14), fcent_15(ap.fcent_15), fcent_16(ap.fcent_16)
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
	if (fMainList) {delete fMainList;}
}

//________________________________________________________________________
void AliJCatalystTask::UserCreateOutputObjects()
{
	fInputList = new TClonesArray("AliJBaseTrack" , 2500);
	fInputList->SetOwner(kTRUE);
	fInputListALICE = new TClonesArray("AliJBaseTrack" , 2500);
	fInputListALICE->SetOwner(kTRUE);

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	OpenFile(1);
	
	fJCorMapTask = (AliJCorrectionMapTask*)man->GetTask(fJCorMapTaskName);
	if(!fJCorMapTask ) AliInfo("----CHECK if AliJCorrectionMapTask Missing ----");
	if( fJCorMapTask ) {
		fCentBinEff = fJCorMapTask->GetCentBinEff();
		fCentBinEff->Print();
	}

	gRandom->SetSeed();

	BookControlHistograms();
	PostData(1, fMainList);

	//OpenFile(1);
	//fOutput = gDirectory;
	//fOutput->cd();
	//PostData(1, fOutput);
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

	float fImpactParameter = .0; // setting 0 for the generator which doesn't have this info. 
	double fvertex[3];

	fEvtNum++;
	if(fEvtNum % 1000 == 0)
		cout << "evt : " << fEvtNum << ", fEntry: " << fEntry <<  endl;

	// load current event and save track, event info
	if(flags & FLUC_KINEONLY) {
		AliMCEvent *mcEvent;
		if(flags & FLUC_KINEONLYEXT) {
			AliInputEventHandler*  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
			mcEvent = fMcHandler->MCEvent();
		}
		else {
			mcEvent = MCEvent();
		}

		if (!mcEvent) {
			AliError("ERROR: mcEvent not available");
			return;
		}

		if(!fnoCentBin) {
			AliGenHijingEventHeader* hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
			AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(mcEvent->GenEventHeader());
			AliGenHepMCEventHeader* hepHeader = dynamic_cast<AliGenHepMCEventHeader*>(mcEvent->GenEventHeader());
			if (hijingHeader) {
				fImpactParameter = hijingHeader->ImpactParameter();
    	}
    	else if (dpmHeader) {
				fImpactParameter = dpmHeader->ImpactParameter();
    	}
    	else if (hepHeader) {
				fImpactParameter = hepHeader->impact_parameter();
			}
			else {
				DEBUG( 4,  "KineOnly no header in event generator" );       			      
			}

			fcent = GetCentralityFromImpactPar(fImpactParameter);
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
		if(flags & FLUC_KINEONLYEXT) {
			ReadKineTracks( mcEvent->Stack(), fInputList, fInputListALICE, fcent ) ; // read tracklist
		}
		else {
			ReadKineTracks( mcEvent, fInputList, fInputListALICE, fcent ) ; // read tracklist
		}

		AliGenEventHeader *header = mcEvent->GenEventHeader();
		if(!header)
			return;
		TArrayF gVertexArray;
		header->PrimaryVertex(gVertexArray);
		for(int i = 0; i < 3; i++)
			fvertex[i] = gVertexArray.At(i);
	}

	else { // Kine
		paodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
		if(fnoCentBin) {
			fcent = 1.0;
		}
		else {
			fcent = ReadCentrality(paodEvent,fCentDetName);
		}
		fRunNum = paodEvent->GetRunNumber();

		// Get the centrality bin for the current centrality event.
		Int_t centBin = GetCentralityBin(fcent);
		if (centBin == -1) {return;}

		ReadVertexInfo(paodEvent, fvertex);	// Read vertex info before selection.
		if (bSaveAllQA) {
			fVertexXHistogram[centBin][0]->Fill(fvertex[0]);
			fVertexYHistogram[centBin][0]->Fill(fvertex[1]);
			fVertexZHistogram[centBin][0]->Fill(fvertex[2]);
		}

		fIsGoodEvent = IsGoodEvent(paodEvent, centBin);
		if(!fIsGoodEvent) {
			return;
		}

		if (bSaveAllQA) {	// Save the vertex and centrality if PV passes the selection.
			fVertexXHistogram[centBin][1]->Fill(fvertex[0]);
			fVertexYHistogram[centBin][1]->Fill(fvertex[1]);
			fVertexZHistogram[centBin][1]->Fill(fvertex[2]);
			fCentralityHistogram[centBin]->Fill(fcent);
		}

		// Load correction maps in the event loop
		//if( fEvtNum == 1 && fJCorMapTask ) {

		if(fJCorMapTask) {
			grEffCor = fJCorMapTask->GetEffCorrectionMap(fRunNum,fcent,fEffFilterBit); 
			int fcBin = 
				AliJFFlucAnalysis::GetBin(fcent,AliJFFlucAnalysis::BINNING_CENT_PbPb);
			pPhiWeights = fJCorMapTask->GetCorrectionMap(phiMapIndex,fRunNum,fcBin);
		}

		ReadAODTracks( paodEvent, fInputList, fcent ) ; // read tracklist
		ReadVertexInfo( paodEvent, fvertex); // read vertex info
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
	if(flags & FLUC_MC) {  // how to get a flag to check  MC or not !
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
				} // weak decay particles are excluded.

				if(fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}

				Int_t pdg = track->GetPdgCode();
				Char_t ch = (Char_t) track->Charge();
				if(ch < 0){
					if(fPcharge == 1)
						continue;
				}
				else if(ch > 0){
					if(fPcharge == -1)
						continue;
				}
				else continue;
        
        if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue; // Need to check this here also
				AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
				itrack->SetLabel(track->GetLabel());
				itrack->SetParticleType( pdg);
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetCharge(ch) ;
			}
		}	// for( int it=0; it < nt ; it++){
	} // end: if(flags & FLUC_MC)

	else {	// Read AOD reco tracks.
		Int_t nt = aod->GetNumberOfTracks();
		if (bSaveAllQA) {fMultHistogram[GetCentralityBin(fcent)][0]->Fill(nt);}

		Int_t ntrack =0;
		Double_t PV[3] = {0.};
		ReadVertexInfo(aod, PV);

		for( int it=0; it<nt ; it++){
			AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
			if(!track) {
				Error("ReadEventAOD", "Could not read particle %d", (int) it);
				continue;
			}

			if (bSaveAllQA) {FillControlHistograms(track, 0, fcent, PV);}	// Fill the QA histograms before the track selection.

			if(track->GetTPCNcls() < fNumTPCClusters)
				continue;

			if(track->TestFilterBit( fFilterBit )){ //
				if( fPt_min > 0){
					double Pt = track->Pt();
					if( Pt < fPt_min || Pt > fPt_max )
						continue ; // pt cut
				}

				Char_t ch = (Char_t)track->Charge();
				if(ch < 0){
					if(fPcharge == 1)
						continue;
				}
				else if(ch > 0){
					if(fPcharge == -1)
						continue;
				}
				else continue;

                if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue; // Need to check this here also
                // Removal of bad area, now only with eta symmetric
				Bool_t isBadArea = TMath::Abs(track->Eta()) > 0.6;
				if(fremovebadarea) {
					if(isBadArea) continue;
				} 
				
				if (bSaveAllQA) {FillControlHistograms(track, 1, fcent, PV);}	// Fill the QA histograms after the track selection.

				AliJBaseTrack *itrack = new( (*TrackList)[ntrack++]) AliJBaseTrack;
				itrack->SetID( TrackList->GetEntriesFast() );
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetParticleType(kJHadron);
				itrack->SetCharge(track->Charge() );
				itrack->SetStatus(track->GetStatus() );
				// Apply eff correction !!!
				Double_t effCorr = 1.0;
				if(grEffCor && fJCorMapTask) {
					effCorr = fJCorMapTask->GetEffCorrection(grEffCor,track->Pt());//fEfficiency->GetCorrection(pt,fEffFilterBit,fCent);
				} 
				itrack->SetTrackEff(effCorr);
				
				// Adding phi weight for a track
				Double_t phi_module_corr = 1.0;
				if(fJCorMapTask){
					Double_t w;
					if(pPhiWeights) {
						Double_t phi = itrack->Phi();
						Double_t eta = itrack->Eta();
						w = pPhiWeights->GetBinContent(pPhiWeights->FindBin(phi,eta,fZvert));
					}
					else {
						w = 1.0;
					}
					if(w > 1e-6) phi_module_corr = w;
				}
				itrack->SetWeight(phi_module_corr);
			}	// End: if(track->TestFilterBit( fFilterBit ))
		}

	} //read aod reco track done.
	if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
	if (bSaveAllQA) {fMultHistogram[GetCentralityBin(fcent)][1]->Fill(TrackList->GetEntriesFast());}
    
}
//______________________________________________________________________________
Bool_t AliJCatalystTask::IsGoodEvent( AliAODEvent *event, Int_t thisCent){
	//event selection here!
	//check vertex
	AliVVertex *vtx = event->GetPrimaryVertex();
	if(!vtx || vtx->GetNContributors() <= 0)
		return kFALSE;
	double zvert = vtx->GetZ();
	if(zvert < -fzvtxCut || zvert > fzvtxCut)
		return kFALSE;

	if(flags & FLUC_CENT_FLATTENING){
		float fCent = ReadCentrality(event,fCentDetName);
		TH1 *pweightMap = fJCorMapTask->GetCentCorrection();
		if(gRandom->Uniform(0,1) > pweightMap->GetBinContent(pweightMap->GetXaxis()->FindBin(fCent)))
			return kFALSE;
	}

	if(flags & FLUC_KINEONLY)
		return kTRUE;

	//int frunNumber = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEvent()->GetRunNumber();
	AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
	fRunTable->SetRunNumber(fRunNum);

	int fperiod = fRunTable->GetRunNumberToPeriod(fRunNum);
	if(fperiod == AliJRunTable::kLHC15o || fperiod == AliJRunTable::kLHC18q || fperiod == AliJRunTable::kLHC18r){
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
		if(fperiod == AliJRunTable::kLHC15o || fperiod == AliJRunTable::kLHC18q || fperiod == AliJRunTable::kLHC18r){
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

		}
		else if (fperiod == AliJRunTable::kLHC10h) {	// High multiplicity outlier cuts for LHC10h based on the SCklm analysis.
			UInt_t MTPC = 0;
			UInt_t Mglobal = 0;
			for(int it = 0; it < nTracks; it++){
				AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
				if (!trackAOD) {continue;}
				if (trackAOD->TestFilterBit(128)) {MTPC++;}
				if (trackAOD->TestFilterBit(256)) {Mglobal++;}
			}

			if (bSaveAllQA && bSaveHMOhist) {fHMOsHistogram[thisCent][0]->Fill(Mglobal, MTPC);}	// Fill the HMO histogram only if both bSave* are true.
			if(!((double)MTPC > (-65.0 + 1.54*Mglobal) && (double)MTPC < (90.0 + 2.30*Mglobal))) {
				return kFALSE;
			}
			if (bSaveAllQA && bSaveHMOhist) {fHMOsHistogram[thisCent][1]->Fill(Mglobal, MTPC);}
		}
		else {	// TBA: QA histo if chosen, possibility to choose between the two versions for LHC10h
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
        switch(FilterBit){
        case 128:
                fEffFilterBit = 0;
                break;
        case 96:
                fEffFilterBit = 4;
                break;
        case 768:
                fEffFilterBit = 5;
                break;
        default:
                fEffFilterBit = 0;
                break;
        }
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
			if(ch < 0){
				if(fPcharge == 1)
					continue;
			}else
			if(ch > 0){
				if(fPcharge == -1)
					continue;
			}else continue;
			if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
			AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
			Int_t label = track->GetLabel();
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
	if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
}
// To read the track generated from a external alievent generators
void AliJCatalystTask::ReadKineTracks( AliStack *stack, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fcent)
{
	Int_t nt = stack->GetNprimary();
	Int_t ntrack = 0;
	for (Int_t it = 0; it < nt; it++) {
		TParticle* track = stack->Particle(it); if(!track) continue;
		if (!AliPWG0Helper::IsPrimaryCharged(track, nt)) continue; 
			double Pt = track->Pt();
			if( Pt < fPt_min || Pt > fPt_max )
				continue ; // pt cut
			Int_t pdg = track->GetPdgCode();
			Char_t ch = (Char_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
			if(ch < 0){
				if(fPcharge == 1)
					continue;
			}else
			if(ch > 0){
				if(fPcharge == -1)
					continue;
			}else continue;
			if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
			AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
			Int_t label = 100;//track->GetLabel();
			itrack->SetLabel( label );
			itrack->SetParticleType(pdg);
			itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->Energy() );
			itrack->SetCharge(ch) ;
			if(TMath::Abs(track->Eta()) < 0.8) {
				AliJBaseTrack *jtrack =  new ((*TrackListALICE)[TrackListALICE->GetEntriesFast()])AliJBaseTrack;
				jtrack->SetLabel( label );
				jtrack->SetParticleType(pdg);
				jtrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->Energy() );
				jtrack->SetCharge(ch) ;
			}
	}
	if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
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

//______________________________________________________________________________
void AliJCatalystTask::InitializeArrays() {
	for(Int_t icent=0; icent<fCentralityBins; icent++){

		fControlHistogramsList[icent] = NULL;
		for(Int_t i=0; i<2; i++) {
		  fVertexXHistogram[icent][i] = NULL;
		  fVertexYHistogram[icent][i] = NULL;
		  fVertexZHistogram[icent][i] = NULL;
		  fTPCClustersHistogram[icent][i] = NULL;
		  fITSClustersHistogram[icent][i] = NULL;
		  fChiSquareTPCHistogram[icent][i] = NULL;
		  fDCAzHistogram[icent][i] = NULL;
		  fDCAxyHistogram[icent][i] = NULL;
	    fChargeHistogram[icent][i] = NULL;
	    fPTHistogram[icent][i] = NULL;
		  fPhiHistogram[icent][i] = NULL;
		  fEtaHistogram[icent][i] = NULL;
		  fMultHistogram[icent][i] = NULL;
		  fHMOsHistogram[icent][i] = NULL;
		}

		fCentralityHistogram[icent] = NULL;
	}
}

//______________________________________________________________________________
void AliJCatalystTask::BookControlHistograms(){

for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
{

	//Check if value of this centrality bin is negative -> if yes: break. We do not need anymore
	if(fcentralityArray[icent+1] < 0)
	{
	   break; //The next edge is a breaking point -> this bin does not exist anymore
	}

	fControlHistogramsList[icent] = new TList();
	fControlHistogramsList[icent]->SetName(Form("ControlHistograms_%.1f-%.1f", fcentralityArray[icent], fcentralityArray[icent+1]));
	fControlHistogramsList[icent]->SetOwner(kTRUE);
	if(bSaveAllQA){	fMainList->Add(fControlHistogramsList[icent]); }

	 // a) Book histogram to hold pt spectra:
	 fPTHistogram[icent][0] = new TH1F("fPTHist_BeforeTrackSelection","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][0]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][0]);

	 fPTHistogram[icent][1] = new TH1F("fPTHist_AfterTrackSelection","Pt Distribution",1000,0.,10.);
	 fPTHistogram[icent][1]->GetXaxis()->SetTitle("P_t");
	 fPTHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPTHistogram[icent][1]);
	 
	 // b) Book histogram to hold phi spectra
	 fPhiHistogram[icent][0] = new TH1F("fPhiHist_BeforeTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][0]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][0]);

	 fPhiHistogram[icent][1] = new TH1F("fPhiHist_AfterTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
	 fPhiHistogram[icent][1]->GetXaxis()->SetTitle("Phi");
	 fPhiHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fPhiHistogram[icent][1]);

	 // c) Book histogram to hold eta distribution before track selection:
	 fEtaHistogram[icent][0] = new TH1F("fEtaHist_BeforeTrackSelection","Eta Distribution",1000,-1.,1.); 
	 fEtaHistogram[icent][0]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][0]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][0]);

	 fEtaHistogram[icent][1] = new TH1F("fEtaHist_AfterTrackSelection","Eta Distribution",1000,-1.,1.);
	 fEtaHistogram[icent][1]->GetXaxis()->SetTitle("Eta");
	 fEtaHistogram[icent][1]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fEtaHistogram[icent][1]);

	 // d) Book histogam to hold multiplicty distributions 
	 fMultHistogram[icent][0] = new TH1F("fMultiHisto_BeforeTrackSelection","Multiplicity",30000,0.,30000.); 
	 fMultHistogram[icent][0]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][0]);
	 
	 fMultHistogram[icent][1] = new TH1F("fMultiHisto_AfterTrackSelection","Multiplicity",30000,0.,30000.); 
	 fMultHistogram[icent][1]->GetXaxis()->SetTitle("Multiplicity M");
	 fControlHistogramsList[icent]->Add(fMultHistogram[icent][1]);

	 // e) Book histogam for Vertex X 
	 fVertexXHistogram[icent][0] = new TH1F("fVertexX_BeforeEventSelection","VertexXBefore",1000,-20.,20.); 
	 fVertexXHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][0]);

	 fVertexXHistogram[icent][1] = new TH1F("fVertexX_AfterEventSelection","VertexXAfter",1000,-20.,20.); 
	 fVertexXHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][1]);

	 // f) Book histogam for Vertex Y 
	 fVertexYHistogram[icent][0] = new TH1F("fVertexY_BeforeEventSelection","VertexYBefore",1000,-20.,20.); 
	 fVertexYHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][0]);

	 fVertexYHistogram[icent][1] = new TH1F("fVertexY_AfterEventSelection","VertexYAfter",1000,-20.,20.); 
	 fVertexYHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][1]);

	 // g) Book histogam for Vertex Z 
	 fVertexZHistogram[icent][0] = new TH1F("fVertexZ_BeforeEventSelection","VertexZBefore",1000,-20.,20.); 
	 fVertexZHistogram[icent][0]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][0]);

	 fVertexZHistogram[icent][1] = new TH1F("fVertexZ_AfterEventSelection","VertexZAfter",1000,-20.,20.); 
	 fVertexZHistogram[icent][1]->GetXaxis()->SetTitle("");
	 fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][1]);

	 // i) Book histogram for number of TPC clustes 
	 fTPCClustersHistogram[icent][0] = new TH1F("fTPCClusters_BeforeCut","TPCClustersBeforeCut",170,0.,170.); 
	 fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][0]);

	 fTPCClustersHistogram[icent][1] = new TH1F("fTPCClusters_AfterCut","TPCClustersAfterCut",170,0.,170.); 
	 fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][1]);

	 //j) Book histogram for number of ITC clusters 
	 fITSClustersHistogram[icent][0] = new TH1F("fITSClusters_BeforeCut","ITSClustersBeforeCut",10,0.,10.); 
	 fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][0]);

	 fITSClustersHistogram[icent][1] = new TH1F("fITSClusters_AfterCut","ITSClustersAfterCut",10,0.,10.); 
	 fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][1]);

	 // k) Book histogram for chi square TPC 
	 fChiSquareTPCHistogram[icent][0] = new TH1F("fChiSquareTPC_BeforeCut","ChiSquareTPCBeforeCut",1000,0.,20.); 
	 fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][0]);

	 fChiSquareTPCHistogram[icent][1] = new TH1F("fChiSquareTPC_AfterCut","ChiSquareTPCAfterCut",1000,0.,20.); 
	 fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][1]);

	  // l) Book histogram for DCAz
	 fDCAzHistogram[icent][0] = new TH1F("fDCAz_BeforeCut","DCAzBeforeCut",1000,-10.,10.);  
	 fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][0]);

	 fDCAzHistogram[icent][1] = new TH1F("fDCAz_AfterCut","DCAzAfterCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][1]);
	 
	 // m) Book histogram for DCAxy
	 fDCAxyHistogram[icent][0] = new TH1F("fDCAxy_BeforeCut","DCAxyBeforeCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][0]);

	 fDCAxyHistogram[icent][1] = new TH1F("fDCAxy_AfterCut","DCAxyAfterCut",1000,-10.,10.); 
	 fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][1]); 

	 // n) Book histogram for Charge
	 fChargeHistogram[icent][0] = new TH1I("fCharge_BeforeCut","ChargeBeforeCut",11,-5.5,5.5); 
	 fControlHistogramsList[icent]->Add(fChargeHistogram[icent][0]);

	 fChargeHistogram[icent][1] = new TH1I("fCharge_AfterCut","ChargeAfterCut",11,-5.5,5.5); 
	 fControlHistogramsList[icent]->Add(fChargeHistogram[icent][1]);

	 // o) Book histogram Centrality 
	 fCentralityHistogram[icent]= new TH1F("fCentralityHistogram_After","CentralityHistogramAfter",22,0.,110.);
	 fCentralityHistogram[icent]->GetXaxis()->SetTitle("Centrality");
	 fCentralityHistogram[icent]->SetLineColor(4);
	 fControlHistogramsList[icent]->Add(fCentralityHistogram[icent]);

	 // p) Book the TH2D for the HMOs in LHC10h.
	 fHMOsHistogram[icent][0] = new TH2D("fHMOsHistogram_Before","Correlations before HMO cuts", 1000, 0.,5000., 1000, 0., 5000.);
	 fHMOsHistogram[icent][0]->GetXaxis()->SetTitle("M_{global}");
	 fHMOsHistogram[icent][0]->GetYaxis()->SetTitle("M_{TPC}");
	 if (bSaveHMOhist) {fControlHistogramsList[icent]->Add(fHMOsHistogram[icent][0]);}

	 fHMOsHistogram[icent][1] = new TH2D("fHMOsHistogram_After","Correlations after HMO cuts", 1000, 0.,5000., 1000, 0., 5000.);
	 fHMOsHistogram[icent][1]->GetXaxis()->SetTitle("M_{global}");
	 fHMOsHistogram[icent][1]->GetYaxis()->SetTitle("M_{TPC}");
	 if (bSaveHMOhist) {fControlHistogramsList[icent]->Add(fHMOsHistogram[icent][1]);}

  }//for(Int_t icent=0; icent<fCentralityBins; icent++)
}

//______________________________________________________________________________
void AliJCatalystTask::FillControlHistograms(AliAODTrack *thisTrack, Int_t whichHisto, Float_t cent, Double_t *v) {

// Get the corresponding centrality bin.
	Int_t CentralityBin = GetCentralityBin(cent);

	Float_t ValueDCAxy = 999.;   // DCA in the xy-plane.
	Float_t ValueDCAz = 999.;    // DCA along z.

	if (fFilterBit == 128)  // These methods work only for constrained TPConly tracks.
	{ //These two quantities are the DCA from global tracks but not what we will cut on.
	  ValueDCAxy = thisTrack->DCA();
	  ValueDCAz = thisTrack->ZAtDCA();
	}
	else  //For the unconstrained tracks.
	{
	  Double_t pos[3];  //Coordinates of the track closest to PV?

	  thisTrack->GetXYZ(pos);
	  ValueDCAxy = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
	  ValueDCAz = pos[2] - v[2];
	}

	fPTHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Pt());
	fPhiHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Phi());
  fEtaHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Eta());
  fTPCClustersHistogram[CentralityBin][whichHisto]->Fill(thisTrack->GetTPCNcls());
  fITSClustersHistogram[CentralityBin][whichHisto]->Fill(thisTrack->GetITSNcls());
  fChiSquareTPCHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Chi2perNDF());
  fDCAzHistogram[CentralityBin][whichHisto]->Fill(ValueDCAz);
  fDCAxyHistogram[CentralityBin][whichHisto]->Fill(ValueDCAxy);
  fChargeHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Charge());
}

//______________________________________________________________________________
Int_t AliJCatalystTask::GetCentralityBin(Float_t cent)
{
  //if this functions returns a negative value -> Error. No Centrality could be selected
	
  //Check for centrality bin
  for(Int_t icent=0; icent<fCentralityBins+1; icent++) //loop over all centrality bins
  {
		if(fcentralityArray[icent]<0) {return -1;}
		if(cent >= fcentralityArray[icent]) { continue; } 
		else { return icent-1; } 
  }	

  //We went through all centrality edges without returning. This means, that the measured value is bigger than the maximum centrality that we want for our analyis
  return -1;

}

//______________________________________________________________________________
void AliJCatalystTask::SetInitializeCentralityArray()
{

 TString sMethodName = "void AliJCatalystTask::SetInitializeCentralityArray()";

  Float_t ListCentralities[17] = { fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16 };

  for(Int_t i=0; i<17; i++) { fcentralityArray[i] = ListCentralities[i]; }

  //Protections
  if(fcentralityArray[0] < 0 || fcentralityArray[1] < 0) { Fatal(sMethodName.Data(),"First Centrality bin not defined"); } //They need at least one well defined centrality bin

  for(Int_t icent=0; icent<fCentralityBins; icent++)
  {
		//The next bin should be a valid boundery, i.e. it is > 0. but it is also smaller than the previous boundery -> Wrong ordering
		if( fcentralityArray[icent+1] > 0. && fcentralityArray[icent+1] < fcentralityArray[icent] ) { Fatal(sMethodName.Data(),"Wrong ordering of centrality bounderies"); }
  }
}

