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
#include "AliJIaaAna.h"

#include <TClonesArray.h>

#include "../AliJCard.h"
#include "AliJIaaHistograms.h"
#include "AliJIaaCorrelations.h"
#include "../AliJEventPool.h"

#include "../AliJTrack.h"
//#include "../AliJAcceptanceCorrection.h"

#include "../AliJEfficiency.h"
//#include <iostream>

ClassImp(AliJIaaAna)

//----------------------------------------------------------------------------------------------------------------------------
AliJIaaAna::AliJIaaAna() :
	TObject(),
	fFirstEvent(kTRUE),
	fRunNumber(-999),
	fcent(-999),
	fZvert(-999),
	fenableEP(kFALSE),
 	fPsi2(-999),
	fEPmin(0.),
	fEPmax(90.),
	fjtrigg((particleType)-100),
	fjassoc((particleType)-100),
	fcard(0),
	finputFile(0),
	fInclusiveFile(""),
	fevt(0),
	fhistos(0),
	fcorrelations(0),
	//fAcceptanceCorrection(0x0),
	fassocPool(0),
	ftriggList(0),
	fassocList(0),
	finputList(0),
	fbTriggCorrel(0),
	fbLPCorrel(0),
	fMinimumPt(0),
	fMCTruthRun(false),
	fEfficiency(0),
	fRunTable(0),
	fHadronSelectionCut(0)
{
	// constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJIaaAna::AliJIaaAna(Bool_t execLocal) :
	TObject(),
	fFirstEvent(kTRUE),
	fRunNumber(-999),
	fcent(-999),
	fZvert(-999),
	fenableEP(kFALSE),
 	fPsi2(-999),
	fEPmin(0.),
	fEPmax(90.),
	fjtrigg((particleType)-100),
	fjassoc((particleType)-100),
	fcard(0),
	finputFile(0),
	fInclusiveFile(""),
	fevt(0),
	fhistos(0),
	fcorrelations(0),
	//fAcceptanceCorrection(0x0),
	fassocPool(0),
	ftriggList(0),
	fassocList(0),
	finputList(0),
	fbTriggCorrel(0),
	fbLPCorrel(0),
	fMinimumPt(0),
	fMCTruthRun(false),
	fEfficiency(0),
	fRunTable(0),
	fHadronSelectionCut(0)
{
	// constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJIaaAna::~AliJIaaAna(){
	// destructor

	delete fhistos; //
	//delete fAcceptanceCorrection; //
	delete fcorrelations; //

	delete fassocPool;//

	delete ftriggList; //
	delete fassocList; //

	delete fEfficiency; //
}

//----------------------------------------------------------------------------------------------------------------------------
AliJIaaAna::AliJIaaAna(const AliJIaaAna& obj) :
	TObject(),
	fFirstEvent(obj.fFirstEvent),
	fRunNumber(obj.fRunNumber),
	fcent(obj.fcent),
	fZvert(obj.fZvert),
	fenableEP(obj.fenableEP),
 	fPsi2(obj.fPsi2),
	fEPmin(obj.fEPmin),
	fEPmax(obj.fEPmax),
	fjtrigg(obj.fjtrigg),
	fjassoc(obj.fjassoc),
	fcard(obj.fcard),
	finputFile(obj.finputFile),
	fInclusiveFile(obj.fInclusiveFile),
	fevt(obj.fevt),
	fhistos(obj.fhistos),
	fcorrelations(obj.fcorrelations),
	//fAcceptanceCorrection(obj.fAcceptanceCorrection),
	fassocPool(obj.fassocPool),
	ftriggList(obj.ftriggList),
	fassocList(obj.fassocList),
	finputList(obj.finputList),
	fbTriggCorrel(obj.fbTriggCorrel),
	fbLPCorrel(obj.fbLPCorrel),
	fMinimumPt(obj.fMinimumPt),
	fMCTruthRun(obj.fMCTruthRun),
	fEfficiency(obj.fEfficiency),
	fRunTable(obj.fRunTable),
	fHadronSelectionCut(obj.fHadronSelectionCut)
{
	// copy constructor
	JUNUSED(obj);
}

//----------------------------------------------------------------------------------------------------------------------------
AliJIaaAna& AliJIaaAna::operator=(const AliJIaaAna& obj){
//----------------------------------------------------------------------------------------------------------------------------
	// equal sign operator
	JUNUSED(obj);
	return *this;
}


//----------------------------------------------------------------------------------------------------------------------------
void AliJIaaAna::Initialize() const{
//----------------------------------------------------------------------------------------------------------------------------
	// init

}

//----------------------------------------------------------------------------------------------------------------------------
void AliJIaaAna::UserCreateOutputObjects(){
//----------------------------------------------------------------------------------------------------------------------------
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

	// Create a class for acceptance correction
	//fAcceptanceCorrection = new AliJAcceptanceCorrection(fcard);

	// Set the number of hits per bin required in the acceptance correction histograms
	//int hitsPerBin = fcard->Get("HitsPerBinAcceptance");
	//fAcceptanceCorrection->SetMinCountsPerBinInclusive(hitsPerBin);
	//if(fcard->Get("AcceptanceTestMode") == 1){
	//	fAcceptanceCorrection->SetTestMode(true);
	//}

	// Create the class doing correlation analysis
	fcorrelations = new AliJIaaCorrelations( fcard, fhistos);
	cout<<endl<< " -----" <<endl;

	// If inclusive file is specified, set inclusive sampling for correlation analysis and
	// read the inclusive histograms for the acceptance correction and histogram class
	if( fInclusiveFile.Length() ) {
		fhistos->ReadInclusiveHistos(fInclusiveFile);
		fcorrelations->SetSamplingInclusive(); //kperp background and triangle. Default is flat
		//fAcceptanceCorrection->ReadMixedEventHistograms(fInclusiveFile);
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

	ftriggList  = new TClonesArray(kParticleProtoType[fjtrigg],1500);
	fassocList  = new TClonesArray(kParticleProtoType[fjassoc],1500);
	finputList = NULL;

	// RunHeader was set from the task macro

	//==== Efficiency ====
	fEfficiency = new AliJEfficiency;
	fEfficiency->SetMode( fcard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien
	

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

//----------------------------------------------------------------------------------------------------------------------------
void AliJIaaAna::UserExec(){
//----------------------------------------------------------------------------------------------------------------------------

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
		fRunTable->SetRunNumber( fRunNumber );
		double sqrtS = fRunTable->GetBeamEnergy(fRunTable->GetPeriod());
		cout << "sqrts = "<< sqrtS << ",runnumber = "<< fRunNumber << endl;
		fEfficiency->SetRunNumber( fRunNumber );
		fEfficiency->Load();
		double b_polarity = 1; // find put what ??
		fcorrelations->SetMagneticFieldPolarity(b_polarity);
		fFirstEvent = kFALSE;
	}

	// Load the event to Data Manager and count events
	fhistos->fhEvents->Fill( 0 );
	//----------------------------------------------------------

	int cBin        = fcard->GetBin(kCentrType, fcent);
	if(cBin<0) return;

	int zBin        = fcard->GetBin(kZVertType, fZvert); //should be alway >0; checked in fdmg->IsGoodEvent()
  if (zBin<0) return;

	fhistos->fhZVert[cBin]->Fill(fZvert);

	//------------------------------------------------------------------
	// Triggers and associated
	// finputList is from External Task no need to clear
	//----------------------ooooo---------------------------------------
	if(fjtrigg==kJHadron || fjassoc==kJHadron){
		for( int i = 0; i < finputList->GetEntries(); i++ ){
			triggerTrack = (AliJBaseTrack*)finputList->At(i);
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
	int noAllTriggTracks = finputList->GetEntries();
	fhistos->fhChargedMult[cBin]->Fill(noAllTriggTracks);

	// Filling external values..
	fhistos->fhCentr->Fill(fcent);
	fhistos->fhiCentr->Fill(cBin);
	fhistos->fhEP[cBin]->Fill(fPsi2);

	//----------------------------------------------------
	//----- Generate trigg list             --
	//----------------------------------------------------
	int noTriggs=0;
	ftriggList->Clear();
	for(int itrack=0; itrack<noAllTriggTracks; itrack++){
		triggerTrack = (AliJBaseTrack*)finputList->At(itrack);
		triggerTrack->SetTriggBin( fcard->GetBin(kTriggType, triggerTrack->Pt()) );

		ptt = triggerTrack->Pt();

		effCorr = 1.0/triggerTrack->GetTrackEff();

		if( ptt>fMinimumPt ){
			fhistos->fhChargedPt[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			fhistos->fhChargedPtPublished[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			fhistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
			fhistos->fhChargedEta->Fill(triggerTrack->Eta(), effCorr);
		}

		if( !triggerTrack->IsInTriggerBin() ) continue;
		iptt = triggerTrack->GetTriggBin();
    
		fhistos->fhIphiTrigg[cBin][iptt]->Fill( triggerTrack->Phi(), effCorr);
		fhistos->fhIetaTrigg[cBin][iptt]->Fill( triggerTrack->Eta(), effCorr);

		new ((*ftriggList)[noTriggs++]) AliJBaseTrack(*triggerTrack);
	}

	//--------------------------------------------------
	//---   Generate assoc list and pool             ---
	//--------------------------------------------------
	fassocList->Clear();
	int noAssocs=0;
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
	//cout << "NoTrigg = "<< noTriggs <<", noAssocs = "<< noAssocs << endl;
	if(noAssocs>0 ) fassocPool->AcceptList(fassocList, fcent, fZvert, noAllTriggTracks, fevt);
	//------------------------------------------------------------------
	// Do the Correlation cBin, zBin from the main LOOP.
	//----------------------ooooo---------------------------------------
	RunCorrelations( ftriggList, fassocList, noAllTriggTracks, cBin, zBin); 
	// Mixing is done in RunCorrelations...
	//--------------------------------------------------------------
	// End of the Correlation
	//--------------------------------------------------------------
}


/* These valuses are from EPTask ...
        double fcent;
        double fZvert;
        Bool_t fenableEP;
        double fPsi2;
        double fEPmin;
        double fEPmax;
*/
//----------------------------------------------------------------------------------------------------------------------------
void AliJIaaAna::RunCorrelations(TClonesArray *triggList, TClonesArray *assoList, int noAllTriggTracks, int cBin, int zBin) {
//----------------------------------------------------------------------------------------------------------------------------
	AliJBaseTrack *triggerTrack = NULL;
	AliJBaseTrack *associatedTrack = NULL;
	int noTrigg  = triggList->GetEntriesFast();
	int noAssocs = assoList->GetEntriesFast();
        for(int ii=0;ii<noTrigg;ii++){ // trigger loop
                triggerTrack = (AliJBaseTrack*)ftriggList->At(ii);
                double ptt = triggerTrack->Pt();
                int iptt   = triggerTrack->GetTriggBin();
                if(iptt<0) {
                        cout<<"Not registered trigger ! I better stop here." <<endl;
                        return;
                }
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Correlation w.r.t \psi_2
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if(fenableEP) {
                        double phis = GetPhiS2(triggerTrack->Phi(),fPsi2);
                        phis *= 180./TMath::Pi(); // card in degree 0-360
                        if (!AccecptEPBins( phis, fEPmin, fEPmax)) continue; // EP Enable && in phiS bin
                        fhistos->fhPhiS[cBin][iptt]->Fill(phis);
                }
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                double effCorr = 1.0/triggerTrack->GetTrackEff();
                fhistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(ptt, effCorr);//inclusive

                for(int jj=0;jj<noAssocs;jj++){ // assoc loop
                        associatedTrack = (AliJBaseTrack*)fassocList->At(jj);
                        //-------------------------------------------------------------
                        fcorrelations->FillCorrelationHistograms(kReal, cBin, zBin, triggerTrack, associatedTrack);
                        //-------------------------------------------------------------
                } // end assoc loop
        } // end trigg loop
	// ===== Event mixing =====
	fassocPool->Mix(ftriggList, kAzimuthFill, fcent, fZvert, noAllTriggTracks, fevt, fbLPCorrel);
	// fcent, fZvert from Task and fevt from Exec, fbLPCorrel is global
}

//----------------------------------------------------------------------------------------------------------------------------
void AliJIaaAna::Terminate() {
//----------------------------------------------------------------------------------------------------------------------------
	// termination

	fcorrelations->PrintOut();
	fassocPool->PrintOut();
	fEfficiency->Write();
	fhistos->fHMG->Print();

}

//----------------------------------------------------------------------------------------------------------------------------
particleType  AliJIaaAna::GetParticleType(char const *inchar){
//----------------------------------------------------------------------------------------------------------------------------
	// part type
	for(int i=0;i<kNumberOfParticleTypes;i++) {
		if(strcmp(inchar,kParticleTypeStrName[i])==0) return (particleType)i;
	}
	std::cout<<"Particle type <"<<inchar<<"> not recognized"<<std::endl;
	exit(1);
}

//----------------------------------------------------------------------------------------------------------------------------
double AliJIaaAna::DeltaPhi(double phi1, double phi2) {
//----------------------------------------------------------------------------------------------------------------------------
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	return res>-kJPi*0.5 ? res : kJTwoPi+res ;
}

//----------------------------------------------------------------------------------------------------------------------------
double AliJIaaAna::GetPhiS2(double phit, double ep) {
//----------------------------------------------------------------------------------------------------------------------------
	double dphi = (phit - ep)/TMath::Pi() ;
	while( dphi < 0 || dphi>=2 ){
		dphi += dphi<0?2:-2;
	} 
	return dphi*TMath::Pi(); // 0 -2pi
}

//----------------------------------------------------------------------------------------------------------------------------
Bool_t AliJIaaAna::AccecptEPBins(double phis, double min, double max) {
//----------------------------------------------------------------------------------------------------------------------------
	// Selecting all slices [min,max],[180-max,180-min],[180+min,180+max],[360-max,360-min]
	// give min max 0-90 degree intuitive... only works for max<180
	Bool_t accept = kFALSE;
	if( phis>min      && phis<max      ) accept = kTRUE; 
	if( phis>180.-max && phis<180.-min ) accept = kTRUE; 
	if( phis>180.+min && phis<180.+max ) accept = kTRUE; 
	if( phis>360.-max && phis<360.-min ) accept = kTRUE; 

	return accept;
}
//----------------------------------------------------------------------------------------------------------------------------

