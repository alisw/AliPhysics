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
#include "AliJHistos.h"

#include "AliJJet.h"
#include "AliJJetQAAna.h"
#include "AliJHistManager.h"

AliJJetQAAna::AliJJetQAAna():
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
   
AliJJetQAAna::AliJJetQAAna( AliJCard * card ):
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

AliJJetQAAna::AliJJetQAAna( const AliJJetQAAna & obj ):
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

AliJJetQAAna& AliJJetQAAna::operator=(const AliJJetQAAna & obj){
    if( this != &obj ){
        this->~AliJJetQAAna();
        new(this) AliJJetQAAna(obj);
    }
    return *this;
}

AliJJetQAAna::~AliJJetQAAna(){
  // destructor
	delete fHistos;
	
}

void AliJJetQAAna::UserCreateOutputObjects(){
    // comment needed
    fJetPtBins = fCard->GetVector("Jet:PtBins");
    fJetPtMinCut = fCard->Get("Jet:PtMinCut");

	fHistos = new AliJHistos(fCard);
	fHistos->CreateEventTrackHistos();
	fHistos->CreateJetHistos();
	fHistos->fHMG->Print();

	fHistos->fHMG->WriteConfig();

}

void AliJJetQAAna::ClearBeforeEvent(){

}

void AliJJetQAAna::UserExec(){
	fevt++;
	fHistos->fhZVert[cBin]->Fill(zVert);
	fHistos->fhiCentr->Fill(fcent);
	TObjArray *jets = (TObjArray*) fJetListOfList[fTargetJetIndex]; // only for one jet finder
	if(!jets) 	return;
	FillHistosJets(jets);
	
}

void AliJJetQAAna::Terminate() const{
	// comment needed
}


void AliJJetQAAna::FillHistosJets(TObjArray *Jets){
	for(int i=0; i<Jets->GetEntries();i++) {
		AliJJet *jet = (AliJJet*)Jets->At(i);
		fHistos->fhJetPt[cBin]->Fill(jet->Pt());
	}
}

void AliJJetQAAna::CreateHistos(){
	// Create Histograms
}

