/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
//===========================================================

#include <TRandom.h>
#include <TMath.h>

#include "AliJBaseTrack.h"
#include "AliJCorrelations.h"
#include "AliJEventPool.h"
#include "AliJHistos.h"
#include "AliJEfficiency.h"

#include "AliJJet.h"
#include "AliJJetCORRAnalysis.h"
#include "AliJHistManager.h"

AliJJetCORRAnalysis::AliJJetCORRAnalysis():
    fInputList(NULL),
    fJetList(NULL),
    fJetListOfList(),
    fIsFullAcceptance(0),
    fJetEtaRange(0),
    fTrackJetMap(NULL),
    fJetPtBins(NULL),
    fJetPtMinCut(0),
    fCard(0x0),
    fHistos(0x0),
    fEfficiency(0x0),
    ftriggList(0x0),
    fassocList(0x0),
    fcorrelations(0x0),
    fassocPool(0x0),
    cBin(-1),
	fcent(-999),
    zBin(-1),
	zVert(-999),
    fevt(0),
	fTargetJetIndex(0),
    fDebugMode(0)
{
    // Constaractor
}
   
AliJJetCORRAnalysis::AliJJetCORRAnalysis( AliJCard * card ):
    fInputList(NULL),
    fJetList(NULL),
    fJetListOfList(),
    fIsFullAcceptance(0),
    fJetEtaRange(0),
    fTrackJetMap(NULL),
    fJetPtBins(NULL),
    fJetPtMinCut(0),
    fCard(card),
    fHistos(0x0),
    fEfficiency(0x0),
    ftriggList(0x0),
    fassocList(0x0),
    fcorrelations(0x0),
    fassocPool(0x0),
    cBin(-1),
	fcent(-999),
    zBin(-1),
	zVert(-999),
    fevt(0),
	fTargetJetIndex(0),
    fDebugMode(0)
{
    // Constaractor
}

AliJJetCORRAnalysis::AliJJetCORRAnalysis( const AliJJetCORRAnalysis & obj ):
    fInputList(obj.fInputList),
    fJetList(obj.fJetList),
    fJetListOfList(obj.fJetListOfList),
    fIsFullAcceptance(obj.fIsFullAcceptance),
    fJetEtaRange(obj.fJetEtaRange),
    fTrackJetMap(obj.fTrackJetMap),
    fJetPtBins(obj.fJetPtBins),
	fJetPtMinCut(obj.fJetPtMinCut),
    fCard(obj.fCard),
    fHistos(obj.fHistos),
    fEfficiency(obj.fEfficiency),
    ftriggList(obj.ftriggList),
    fassocList(obj.fassocList),
    fcorrelations(obj.fcorrelations),
    fassocPool(obj.fassocPool),
    cBin(-1),
	fcent(-999),
    zBin(-1),
	zVert(-999),
    fevt(0),
	fTargetJetIndex(0),
    fDebugMode(0)
{
    // Constaractor
}

AliJJetCORRAnalysis& AliJJetCORRAnalysis::operator=(const AliJJetCORRAnalysis & obj){
    if( this != &obj ){
        this->~AliJJetCORRAnalysis();
        new(this) AliJJetCORRAnalysis(obj);
    }
    return *this;
}

AliJJetCORRAnalysis::~AliJJetCORRAnalysis(){
  // destructor
	delete fInputList;
	delete ftriggList;
	delete fassocList;
	delete fHistos;
	delete fcorrelations;
	delete fassocPool;
	delete fEfficiency;
	
}

void AliJJetCORRAnalysis::UserCreateOutputObjects(){
    // comment needed
    fJetPtBins = fCard->GetVector("Jet:PtBins");
    fJetPtMinCut = fCard->Get("Jet:PtMinCut");

    fInputList  = new TClonesArray("AliJBaseTrack",1500);
    fInputList->SetOwner(kTRUE);
    ftriggList  = new TClonesArray("AliJBaseTrack",1500);
    fassocList  = new TClonesArray("AliJBaseTrack",1500);

	fHistos = new AliJHistos(fCard);
	fHistos->CreateEventTrackHistos();
	fHistos->CreateAzimuthCorrHistos();
	fHistos->CreateIAAMoons();
	fHistos->CreateXEHistos();
	fHistos->CreateXtHistos();
	fHistos->CreatePairPtCosThetaStar();

	fHistos->fHMG->Print();

	fcorrelations = new AliJCorrelations(fCard, fHistos);
	fassocPool   = new AliJEventPool( fCard, fHistos, fcorrelations, kJHadron);

	fEfficiency = new AliJEfficiency();
	fEfficiency->SetMode( fCard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien

	fHistos->fHMG->WriteConfig();

}

void AliJJetCORRAnalysis::ClearBeforeEvent(){

}

void AliJJetCORRAnalysis::UserExec(){
	fevt++;
	TObjArray *jets = (TObjArray*) fJetListOfList[fTargetJetIndex]; // only for one jet finder
	if(!jets) 	return;
	ftriggList->Clear();
	fassocList->Clear();
	//Make Trigger list and assocList
	int noTrigg=0,noAssoc=0;
	for(int ii=0;ii<jets->GetEntriesFast();ii++) {
		AliJJet *jet = (AliJJet*)jets->At(ii);
		jet->SetTriggBin( fCard->GetBin(kTriggType, jet->Pt()) );
		jet->SetTrackEff( 1.0 ); 
		if( jet->GetTriggBin() > -1 ) { new ((*ftriggList)[noTrigg++]) AliJBaseTrack(*jet); }

		TObjArray *constituents = (TObjArray*)jet->GetConstituents();
		for(int jj=0;jj<constituents->GetEntriesFast();jj++) {
			AliJBaseTrack *con = (AliJBaseTrack*)constituents->At(jj);
			con->SetAssocBin( fCard->GetBin(kAssocType, con->Pt()) );
			con->SetParticleType(kJHadron);
			double pt = con->Pt();
			con->SetTrackEff( fEfficiency->GetCorrection(pt, 5, fcent) ); // only here for hybrid cut = 5
			if( con->GetAssocBin() > -1 ) { new ((*fassocList)[noAssoc++]) AliJBaseTrack(*con); }
			if(con->GetAssocBin() > -1){
				int ipta  = con->GetAssocBin();
				double effCorrection = 1.0/con->GetTrackEff();
				fHistos->fhIphiAssoc[cBin][ipta]->Fill( con->Phi(), effCorrection);
				fHistos->fhIetaAssoc[cBin][ipta]->Fill( con->Eta(), effCorrection);
			}
		}
	}
    // correlation loop
    if(fDebugMode) cout << "Start of Correlation Loop noTrigg = "<< ftriggList->GetEntriesFast()<<"\t noAssoc="<<fassocList->GetEntriesFast()<<endl;
    if(fassocList->GetEntriesFast()>0 ) fassocPool->AcceptList(fassocList, fcent, zVert, fInputList->GetEntriesFast(), fevt);
    PlayCorrelation(ftriggList, fassocList);
    fassocPool->Mix(ftriggList, kAzimuthFill, fcent, zVert, fInputList->GetEntriesFast(), fevt);
}
void AliJJetCORRAnalysis::Terminate() const{
	// comment needed
}


void AliJJetCORRAnalysis::FillHistosJets(TObjArray *Jets){
	for(int i=0; i<Jets->GetEntries();i++) {
		AliJJet *jet = (AliJJet*)Jets->At(i);
		jet->Print();
	}
}

void AliJJetCORRAnalysis::CreateHistos(){
	// Create Histograms
}


//________________________________________________________________________
void AliJJetCORRAnalysis::PlayCorrelation(TClonesArray *triggList, TClonesArray *assocList){
	int noTrigg = triggList->GetEntries();
	int noAssoc = assocList->GetEntries();

	fHistos->fhAssocMult->Fill(noAssoc);
	//------------------------------------------------------------------
	//==== Correlation Loop
	//------------------------------------------------------------------
	for(int ii=0;ii<noTrigg;ii++){ // trigger loop
		AliJBaseTrack * triggTr = (AliJBaseTrack*)triggList->At(ii);
		int iptt   = triggTr->GetTriggBin();
		if( iptt < 0 ) continue;
		double effCorr = 1.0/triggTr->GetTrackEff();
		fHistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(triggTr->Pt(), effCorr);//inclusive
		for(int jj=0;jj<noAssoc;jj++){ // assoc loop
			AliJBaseTrack  *assocTr = (AliJBaseTrack*)assocList->At(jj);
			fcorrelations->FillAzimuthHistos(kReal, cBin, zBin, triggTr, assocTr); // cBin and zBin from the members in header vis setter from TaskCode.
		}
	} // end of trigg
}
