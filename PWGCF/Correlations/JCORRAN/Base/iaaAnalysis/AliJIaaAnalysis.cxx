/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// iaaAnalysis main class
// used in local and grid execution

#include <TH1D.h>
#include "AliJIaaAnalysis.h"

#include "AliJTrackCounter.h"
#include <TClonesArray.h>

#include "../AliJCard.h"
#include "AliJIaaHistograms.h"
#include "AliJIaaCorrelations.h"
#include "../AliJEventPool.h"

#include "../AliJEventHeader.h"
#include "../AliJRunHeader.h"
#include "../AliJTrack.h"
#include "../AliJPhoton.h"
#include "../AliJMCTrack.h"
#include "../AliJAcceptanceCorrection.h"



#include "../AliJEfficiency.h"
#include <iostream>

ClassImp(AliJIaaAnalysis)

AliJIaaAnalysis::AliJIaaAnalysis() :
	TObject(),
	fExecLocal(kFALSE),
	fFirstEvent(kTRUE),
	fjtrigg((particleType)-100),
	fjassoc((particleType)-100),
	fcard(0),
	finputFile(0),
	fInclusiveFile(""),
	fevt(0),
	fhistos(0),
	fcorrelations(0),
	fAcceptanceCorrection(0x0),
	fassocPool(0),
	fchargedHadronList(0),
	ftriggList(0),
	fassocList(0),
	finputList(0),
	fdmg(0),
	feventHeader(0),
	frunHeader(0),
	fcent(0),
	fbTriggCorrel(0),
	fbLPCorrel(0),
	fMinimumPt(0),
    fMCTruthRun(false),
	fEventBC(0),
	fEfficiency(0),
	fRunTable(0),
	fHadronSelectionCut(0)
{
	// constructor
}

AliJIaaAnalysis::AliJIaaAnalysis(Bool_t execLocal) :
	TObject(),
	fExecLocal(execLocal),
	fFirstEvent(kTRUE),
	fjtrigg((particleType)-100),
	fjassoc((particleType)-100),
	fcard(0),
	finputFile(0),
	fInclusiveFile(""),
	fevt(0),
	fhistos(0),
	fcorrelations(0),
	fAcceptanceCorrection(0x0),
	fassocPool(0),
	fchargedHadronList(0),
	ftriggList(0),
	fassocList(0),
	finputList(0),
	fdmg(0),
	feventHeader(0),
	frunHeader(0),
	fcent(0),
	fbTriggCorrel(0),
	fbLPCorrel(0),
	fMinimumPt(0),
    fMCTruthRun(false),
	fEventBC(0),
	fEfficiency(0),
	fRunTable(0),
	fHadronSelectionCut(0)
{
	// constructor
}

AliJIaaAnalysis::~AliJIaaAnalysis(){
	// destructor

	delete fhistos;
	delete fAcceptanceCorrection;
	delete fcorrelations;

	delete fassocPool;

	delete fchargedHadronList;
	delete ftriggList;
	delete fassocList;

	delete fdmg;
	delete fEfficiency;
}

AliJIaaAnalysis::AliJIaaAnalysis(const AliJIaaAnalysis& obj) :
	TObject(),
	fExecLocal(obj.fExecLocal),
	fFirstEvent(obj.fFirstEvent),
	fjtrigg(obj.fjtrigg),
	fjassoc(obj.fjassoc),
	fcard(obj.fcard),
	finputFile(obj.finputFile),
	fInclusiveFile(obj.fInclusiveFile),
	fevt(obj.fevt),
	fhistos(obj.fhistos),
	fcorrelations(obj.fcorrelations),
	fAcceptanceCorrection(obj.fAcceptanceCorrection),
	fassocPool(obj.fassocPool),
	fchargedHadronList(obj.fchargedHadronList),
	ftriggList(obj.ftriggList),
	fassocList(obj.fassocList),
	finputList(obj.finputList),
	fdmg(obj.fdmg),
	feventHeader(obj.feventHeader),
	frunHeader(obj.frunHeader),
	fcent(obj.fcent),
	fbTriggCorrel(obj.fbTriggCorrel),
	fbLPCorrel(obj.fbLPCorrel),
	fMinimumPt(obj.fMinimumPt),
    fMCTruthRun(obj.fMCTruthRun),
	fEventBC(obj.fEventBC),
	fEfficiency(obj.fEfficiency),
	fRunTable(obj.fRunTable),
	fHadronSelectionCut(obj.fHadronSelectionCut)
{
	// copy constructor
	JUNUSED(obj);
}

AliJIaaAnalysis& AliJIaaAnalysis::operator=(const AliJIaaAnalysis& obj){
	// equal sign operator
	JUNUSED(obj);
	return *this;
}


void AliJIaaAnalysis::Initialize() const{
	// init

}

void AliJIaaAnalysis::UserCreateOutputObjects(){
	// local init


	cout << "jtAnalysis user create output objects ----------------" << endl;

	fHadronSelectionCut =int ( fcard->Get("HadronSelectionCut"));

    // Switch between regular analysis and MC truth analysis
    fMCTruthRun = false;
    if(fcard->Get("AnalyseMCTruth") == 1) fMCTruthRun = true;

	// Initialize the histograms needed to store the output
	fhistos = new AliJIaaHistograms( fcard );
	//if(fcard->Get("QualityControlLevel")>1) fhistos->Set2DHistoCreate(true);
	//if(fcard->Get("QualityControlLevel")>0) fhistos->SetAcceptanceCorrectionQA(true);
	fhistos->CreateEventTrackHistos();
	fhistos->CreateCorrelationHistograms();

	fhistos->fHMG->Print();

	fEventBC = (Int_t)(fcard->Get( "eventBC" ));

	// Create a class for acceptance correction
	fAcceptanceCorrection = new AliJAcceptanceCorrection(fcard);

	// Set the number of hits per bin required in the acceptance correction histograms
	int hitsPerBin = fcard->Get("HitsPerBinAcceptance");
	fAcceptanceCorrection->SetMinCountsPerBinInclusive(hitsPerBin);
	if(fcard->Get("AcceptanceTestMode") == 1){
		fAcceptanceCorrection->SetTestMode(true);
	}

	// Create the class doing correlation analysis
	fcorrelations = new AliJIaaCorrelations( fcard, fhistos);
	cout<<endl<< " -----" <<endl;

	// If inclusive file is specified, set inclusive sampling for correlation analysis and
	// read the inclusive histograms for the acceptance correction and histogram class
	if( fInclusiveFile.Length() ) {
		fhistos->ReadInclusiveHistos(fInclusiveFile);
		fcorrelations->SetSamplingInclusive(); //kperp background and triangle. Default is flat
		fAcceptanceCorrection->ReadMixedEventHistograms(fInclusiveFile);
		cout<<"Background and acceptance sampling from " << fInclusiveFile <<endl;
	} else {
		cout << "Background and acceptance sampled from flat distributions." <<endl;
	}
	cout<< " -----" <<endl <<endl;

	// Tell the correlation analysis to use the defined acceptance correction
	//fcorrelations->SetAcceptanceCorrection(fAcceptanceCorrection);
	if(fcard->Get("UseZVertexBins") == 1){
		fcorrelations->UseZVertexAcceptance(true);
	}


	// EventPool for Mixing
	fassocPool   = new AliJEventPool( fcard, fhistos, fcorrelations, fjassoc);

	fchargedHadronList  = new TClonesArray(kParticleProtoType[fMCTruthRun ? kJHadronMC : kJHadron],1500);
	ftriggList  = new TClonesArray(kParticleProtoType[fjtrigg],1500);
	fassocList  = new TClonesArray(kParticleProtoType[fjassoc],1500);
	finputList = NULL;

	fdmg = new AliJDataManager(fcard, fhistos, fcorrelations, fExecLocal);
	fdmg->SetExecLocal( fExecLocal );

	//==== Read the Data files =====
	if( fExecLocal ){
		fdmg->ChainInputStream(finputFile);

		// for grid running, numberEvents is filled by the encapsulating
		// grid task, which has access to the input handlers and can
		// extract event numbers out of it
		int nEvents = fdmg->GetNEvents();
		frunHeader = fdmg->GetRunHeader();
		cout<<"RunID = "<<frunHeader->GetRunNumber()<< " Looping over "<<nEvents<<" events"<<endl;

	}

	// RunHeader was set from the task macro

	//==== Efficiency ====
	fEfficiency = new AliJEfficiency;
	fEfficiency->SetMode( fcard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	if(fExecLocal) {
		fEfficiency->SetDataPath("/mnt/flustre/alice/taxi_jcorran/2013/EFF/data"); // Efficiency root file location local or alien
	} else {
		fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien
	}

	//-------------------------------------------------------------------
	fbTriggCorrel  = fcard->Get("CorrelationType")==0;
	fbLPCorrel     = fcard->Get("CorrelationType")==1;
	fMinimumPt = fcard->GetBinBorder(kAssocType, 0);

	// Initialize counters
	fcent = -1;
	fhistos->fHMG->WriteConfig();
	fFirstEvent = kTRUE;
	fevt = -1;

	cout << "end of jcorran user create output objects ----------------" << endl;

}

void AliJIaaAnalysis::UserExec(){

	// Variables needed inside loops
	AliJBaseTrack *triggerTrack;    // Track for the trigger particle
	AliJBaseTrack *associatedTrack; // Track for the associated particle
	double ptt;                     // pT of the trigger particle
	double effCorr;                 // Efficiency correction
	int iptt;                       // Index of trigger pT bin
	int ipta;                       // Index of associated pT bin

	// event loop
	fevt++;

	if( fFirstEvent ) {

		// ==== Set the RunTable only in the first event ====
		fRunTable = &AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( frunHeader->GetRunNumber() );
		double sqrtS = fRunTable->GetBeamEnergy(fRunTable->GetPeriod());
		cout << "sqrts = "<< sqrtS << ",runnumber = "<< frunHeader->GetRunNumber() << endl;
		fEfficiency->SetRunNumber( frunHeader->GetRunNumber() );
		fEfficiency->Load();
		double b_polarity = (frunHeader->GetL3MagnetFieldIntensity() < 0) ? -1 : 1;
		fcorrelations->SetMagneticFieldPolarity(b_polarity);
		fFirstEvent = kFALSE;
	}

	// Load the event to Data Manager and count events
	fdmg->LoadEvent(fevt);
	fhistos->fhEvents->Fill( 0 );

	if(!fdmg->IsGoodEvent()) return;  // Vertex cut applied in IsGoodEvent and histo saved there too

	// Read event header and z-vertex position
	feventHeader  = fdmg->GetEventHeader();
	double zVert    = feventHeader->GetZVertex();
	//----------------------------------------------------------

	// Read the trigger mask
	UInt_t triggermaskJCorran = feventHeader->GetTriggerMaskJCorran();

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran)){
		fhistos->fhEvents->Fill( 5 );
	}

	// select only some BC%4
	if( feventHeader->GetBunchCrossNumber() % 4 != fEventBC && fEventBC > -1 ){
		return;
	}

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran)){
		fhistos->fhEvents->Fill( 6 );
	}

	if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
		fcent = feventHeader->GetCentrality();
	}else  {
		fcent = 1; //ntracks;
	}

	int cBin        = fcard->GetBin(kCentrType, fcent);
	if(cBin<0) return;

	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran)){
		fhistos->fhEvents->Fill( 7 );
	}

	int zBin        = fcard->GetBin(kZVertType, zVert); //should be alway >0; checked in fdmg->IsGoodEvent()

	// do not fill MB in case of MB mixing
	if( fdmg->IsSelectedTrigger((int) triggermaskJCorran)){
		fhistos->fhZVert[cBin]->Fill(zVert);
	}

	//------------------------------------------------------------------
	// Triggers and associated
	//----------------------ooooo---------------------------------------

	if(fjtrigg==kJHadron || fjassoc==kJHadron || fjtrigg==kJHadronMC || fjassoc==kJHadronMC){
		fchargedHadronList->Clear();
		fdmg->RegisterList(fchargedHadronList, NULL, cBin, zBin, fMCTruthRun ? kJHadronMC : kJHadron);
		// apply efficiencies

		for( int i = 0; i < fchargedHadronList->GetEntries(); i++ ){

			triggerTrack = (AliJBaseTrack*)fchargedHadronList->At(i);
			ptt = triggerTrack->Pt();

            if(fMCTruthRun){  
                fhistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1.0 );
                triggerTrack->SetTrackEff(1.0);
            } else {
                effCorr = 1./fEfficiency->GetCorrection(ptt, fHadronSelectionCut, fcent);  // here you generate warning if ptt>30
                fhistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1./effCorr );
                triggerTrack->SetTrackEff( 1./effCorr );
            }
		}
	}

	//---- assign input list ----
	if(fjtrigg==kJHadron || fjtrigg==kJHadronMC) finputList = fchargedHadronList;
	int noAllTriggTracks = finputList->GetEntries();
	int noAllChargedTracks = fchargedHadronList->GetEntries();
	fhistos->fhChargedMult[cBin]->Fill(noAllChargedTracks);

	fhistos->fhCentr->Fill(fcent);
	fhistos->fhiCentr->Fill(cBin);

	//----------------------------------------------------
	//----- Generate trigg list and find LP             --
	//----------------------------------------------------
	AliJTrackCounter *lpTrackCounter = new AliJTrackCounter();
	AliJBaseTrack *lPTr = NULL;
	int noTriggs=0;
	ftriggList->Clear();
	for(int itrack=0; itrack<noAllTriggTracks; itrack++){
		triggerTrack = (AliJBaseTrack*)finputList->At(itrack);
		triggerTrack->SetTriggBin( fcard->GetBin(kTriggType, triggerTrack->Pt()) );

		ptt = triggerTrack->Pt();

		effCorr = 1.0/triggerTrack->GetTrackEff();

		if( ptt>fMinimumPt ){
			fhistos->fhChargedPt[cBin]->Fill(ptt, effCorr );
            fhistos->fhChargedPtPublished[cBin]->Fill(ptt, effCorr );
			fhistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
			fhistos->fhChargedEta->Fill(triggerTrack->Eta(), effCorr);
		}

		if( !triggerTrack->IsInTriggerBin() ) continue;
		iptt = triggerTrack->GetTriggBin();
		fhistos->fhIphiTrigg[cBin][iptt]->Fill( triggerTrack->Phi(), effCorr);
		fhistos->fhIetaTrigg[cBin][iptt]->Fill( triggerTrack->Eta(), effCorr);

		if( ptt > lpTrackCounter->Pt() ) {
			lpTrackCounter->Store(noTriggs, ptt, iptt);
			lPTr = triggerTrack;
		}

		new ((*ftriggList)[noTriggs++]) AliJBaseTrack(*triggerTrack);
	}

	//--------------------------------------------------
	//---   Generate assoc list and pool             ---
	//--------------------------------------------------
	fassocList->Clear();
	int noAssocs=0;
	if(fjassoc==kJHadron || fjassoc==kJHadronMC) finputList = fchargedHadronList;

	int noAllAssocTracks = finputList->GetEntries();


	for(int itrack=0;itrack<noAllAssocTracks; itrack++){

		associatedTrack = (AliJBaseTrack*)finputList->At(itrack);
		associatedTrack->SetAssocBin( fcard->GetBin(kAssocType, associatedTrack->Pt()) );

		if(associatedTrack->IsInAssocBin()){

			ipta  = associatedTrack->GetAssocBin();
			effCorr = 1.0/associatedTrack->GetTrackEff();
			fhistos->fhIphiAssoc[cBin][ipta]->Fill( associatedTrack->Phi(), effCorr);
			fhistos->fhIetaAssoc[cBin][ipta]->Fill( associatedTrack->Eta(), effCorr);
			new ((*fassocList)[noAssocs++]) AliJBaseTrack(*associatedTrack);
		}
	}

	//-----------------------------------------------
	// Leading particle pT and eta
	//-----------------------------------------------
	if( lpTrackCounter->Exists() ){
        if(fMCTruthRun){
            effCorr = 1.0;
        } else {
            effCorr = 1./fEfficiency->GetCorrection(lpTrackCounter->Pt(), fHadronSelectionCut, fcent );
        }
		fhistos->fhLPpt->Fill(lpTrackCounter->Pt(), effCorr);
		fhistos->fhLPeta->Fill(lPTr->Eta(), effCorr);
	}

	if(noAssocs>0 ) fassocPool->AcceptList(fassocList, fcent, zVert, noAllChargedTracks, fevt);

	//------------------------------------------------------------------
	// Do the Correlation
	//----------------------ooooo---------------------------------------
	int nTriggerTracks=-1;
	nTriggerTracks = fbTriggCorrel ? noTriggs : 1;
	if(fbLPCorrel && !lpTrackCounter->Exists())
	{
		delete lpTrackCounter;
        //delete fbLPCorrel;
		return;
	}
	triggerTrack = NULL;

	for(int ii=0;ii<nTriggerTracks;ii++){ // trigger loop
		if (fbTriggCorrel)  triggerTrack = (AliJBaseTrack*)ftriggList->At(ii);
		if (fbLPCorrel)     triggerTrack = (AliJBaseTrack*)ftriggList->At(lpTrackCounter->GetIndex());

		ptt = triggerTrack->Pt();
		iptt   = triggerTrack->GetTriggBin();
		if(iptt<0) {
			cout<<"Not registered trigger ! I better stop here." <<endl;
			exit(-1);
		}
		effCorr = 1.0/triggerTrack->GetTrackEff();
		fhistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(ptt, effCorr);//inclusive

		for(int jj=0;jj<noAssocs;jj++){ // assoc loop
			associatedTrack = (AliJBaseTrack*)fassocList->At(jj);

			//-------------------------------------------------------------
			fcorrelations->FillCorrelationHistograms(kReal, cBin, zBin, triggerTrack, associatedTrack);
			//-------------------------------------------------------------
		} // end assoc loop
	} // end trigg loop

	// ===== Event mixing =====
	fassocPool->Mix(ftriggList, kAzimuthFill, fcent, zVert, noAllChargedTracks, fevt, fbLPCorrel);

	//--------------------------------------------------------------
	// End of the Correlation
	//--------------------------------------------------------------

	delete lpTrackCounter;


}

void AliJIaaAnalysis::Terminate() {
	// termination

	fcorrelations->PrintOut();
	fassocPool->PrintOut();
	fEfficiency->Write();
	fhistos->fHMG->Print();

}

particleType  AliJIaaAnalysis::GetParticleType(const char *inchar){
	// part type
	for(int i=0;i<kNumberOfParticleTypes;i++) {
		if(strcmp(inchar,kParticleTypeStrName[i])==0) return (particleType)i;
	}
	std::cout<<"Particle type <"<<inchar<<"> not recognized"<<std::endl;
	exit(1);
}

double AliJIaaAnalysis::DeltaPhi(double phi1, double phi2) {
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	return res>-kJPi*0.5 ? res : kJTwoPi+res ;
}
