/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																          *
 * Authors: Svein Lindal                                                  *
 * Version 1.0			                                                  *
 *																          *
 *																          *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is		  *
 * provided "as is" without express or implied warranty.				  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class doing conversion gamma dPhi correlations
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskdPhi.h"

#include <TH2I.h>
#include <TList.h>
#include <TChain.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisFilter.h>

#include "AliConversionTrackCuts.h"
#include "AliConversionCuts.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"
#include "AliAnaConvCorrPhoton.h"
#include "AliAnaConvCorrPion.h"
#include "AliAnaConvIsolation.h"
// Author Svein Lindal <slindal@fys.uio.no>
using namespace std;

ClassImp(AliAnalysisTaskdPhi)


//________________________________________________________________________
AliAnalysisTaskdPhi::AliAnalysisTaskdPhi(const char *name) : AliAnalysisTaskSE(name),
  fHistograms(NULL),
  fHistoGamma(NULL),
  fHistoPion(NULL),
  fV0Filter(NULL),
  fTrackCuts(NULL),
  fGammas(NULL),
  fPions(NULL),
  hMETracks(NULL), 
  hMEPhotons(NULL), 
  hMEPions(NULL),
  hMEvents(NULL),
  fPhotonCorr(NULL),
  fPionCorr(NULL), 
  fIsoAna(NULL),
  fL1(-1),
  fL2(-1),
  fDeltaAODBranchName("AliAODGammaConversion_gamma"), 
  fAxistPt(),
  fAxiscPt(),
  fAxisEta(),
  fAxisPhi(),
  fAxisCent(),
  fAxisZ(), 
  fAxisPiM()
{
  //constructor
  fAxistPt.SetNameTitle("tPtAxis", "tPt");
  fAxistPt.Set(20, 0, 100);

  fAxiscPt.SetNameTitle("cPtAxis", "cPt");
  fAxiscPt.Set(20, 0, 100);

  fAxisEta.SetNameTitle("EtaAxis", "Eta");
  fAxisEta.Set(160, -0.8, 0.8);

  fAxisPhi.SetNameTitle("PhiAxis", "Phi");
  fAxisPhi.Set(128, 0, TMath::TwoPi());

  fAxisZ.SetNameTitle("ZAxis", "Z");
  fAxisZ.Set(4, -10, 10);

  fAxisCent.SetNameTitle("CentAxis", "Cent");


  Double_t centbins[5] = {0, 10, 30, 60, 100.1};
  fAxisCent.Set(4, centbins);

  Double_t mbins[7] = {0.1, 0.11, 0.12, 0.15, 0.16, 0.18, 0.2};
  fAxisPiM.SetNameTitle("InvMassPi0", "Invariant mass");
  fAxisPiM.Set(6, mbins);

  fGammas = new TObjArray();
  fGammas->SetOwner(kFALSE);

  fPions = new TObjArray();
  fPions->SetOwner(kFALSE);

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineInput(1, TClonesArray::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}



//________________________________________________________________________
AliAnalysisTaskdPhi::~AliAnalysisTaskdPhi(){
  //destructor
  if(fPions)
	delete fPions;
  fPions = NULL;

  if(fGammas)
	delete fGammas;
  fGammas = NULL;
  
  if(fIsoAna)
  	delete fIsoAna;
  fIsoAna = NULL;

  if(fV0Filter)
	delete fV0Filter;
  fV0Filter = NULL;

  if(fHistograms)
	delete fHistograms;
  fHistograms = NULL;

if(fHistoPion)
	delete fHistoPion;
  fHistoPion = NULL;

  if(fHistoGamma)
	delete fHistoGamma;
  fHistoGamma = NULL;

}

///________________________________________________________________________
void AliAnalysisTaskdPhi::SetUpCorrObjects() {
  ///Creat corr obj
  fIsoAna = new AliAnaConvIsolation();


  fPhotonCorr = new TObjArray();
  fPionCorr = new TObjArray();
  
  TList * hPhoton = new TList();
  hPhoton->SetName("hPhotonCorr");
  hPhoton->SetOwner(kTRUE);
  fHistoGamma->Add(hPhoton);

  TList * hPion = new TList();
  hPion->SetName("hPionCorr");
  hPion->SetOwner(kTRUE);
  fHistoPion->Add(hPion);


  for(Int_t iz = 0; iz < fAxisZ.GetNbins(); iz++) {
	TObjArray * photonArray = new TObjArray();
	photonArray->SetOwner(kTRUE);
	fPhotonCorr->AddAt(photonArray, iz);

	TObjArray * pionArray = new TObjArray();
	pionArray->SetOwner(kTRUE);
	fPionCorr->AddAt(pionArray, iz);

	TList * photonList = new TList();
	photonList->SetName(Form("photon_%d", iz));
	photonList->SetOwner(kTRUE);
	hPhoton->AddAt(photonList, iz);

	TList * pionList = new TList();
	pionList->SetName(Form("pion_%d", iz));
	pionList->SetOwner(kTRUE);
	hPion->AddAt(pionList, iz);
	
	for(Int_t ic = 0; ic < fAxisCent.GetNbins(); ic++) {
	  
	  TString nameString = Form("%d_%d", iz, ic);
	  TString titleString = Form("%f < Z < %f ... %f cent %f", 
								 fAxisZ.GetBinLowEdge(iz+1), fAxisZ.GetBinUpEdge(iz+1), 
								 fAxisCent.GetBinLowEdge(ic+1), fAxisCent.GetBinUpEdge(ic+1));



	  AliAnaConvCorrPhoton * photonCorr = new AliAnaConvCorrPhoton(Form("PhotonCorr_%s", nameString.Data()), Form("photon %s", titleString.Data()));
	  photonArray->AddAt(photonCorr, ic);
	  //photonCorr->GetAxistPt().Set(fAxistPt.GetNbins(), fAxistPt.GetXbins()->GetArray());
	  //photonCorr->GetAxiscPt().Set(fAxiscPt.GetNbins(), fAxiscPt.GetXbins()->GetArray());
	  photonCorr->CreateHistograms();
	  photonList->Add(photonCorr->GetHistograms());

	  AliAnaConvCorrPion * pionCorr = new AliAnaConvCorrPion(Form("PionCorr_%s", nameString.Data()), Form("pion %s", titleString.Data()));
	  pionArray->AddAt(pionCorr, ic);
	  //pionCorr->GetAxistPt().Set(fAxistPt.GetNbins(), fAxistPt.GetXbins()->GetArray());
	  //pionCorr->GetAxiscPt().Set(fAxiscPt.GetNbins(), fAxiscPt.GetXbins()->GetArray());
	  //pionCorr->GetAxisM().Set(fAxisPiM.GetNbins(), fAxisPiM.GetXbins()->GetArray());
	  pionCorr->CreateHistograms();
	  pionList->Add(pionCorr->GetHistograms());
	}
  }
}

//________________________________________________________________________
void AliAnalysisTaskdPhi::UserCreateOutputObjects() {
  // Create histograms
  
  fHistograms = new TList();
  fHistograms->SetName("dPhi_histograms");
  fHistograms->SetOwner(kTRUE);

  fHistoGamma = new TList();
  fHistoGamma->SetName("Gamma_histo");
  fHistoGamma->SetOwner(kTRUE);

  fHistoPion = new TList();
  fHistoPion->SetName("Pion_histo");
  fHistoPion->SetOwner(kTRUE);

  
  if(fV0Filter) {
	fV0Filter->InitCutHistograms();
	fHistograms->Add(fV0Filter->GetCutHistograms());
  }
  

  AliConversionTrackCuts * tc = dynamic_cast<AliConversionTrackCuts*>(fTrackCuts);
  if(tc) fHistograms->Add(tc->CreateHistograms());

  SetUpCorrObjects();


  ///Set up ME histograms
  TList * MEHistograms = new TList();
  MEHistograms->SetName("MEHistograms");
  MEHistograms->SetOwner(kTRUE);
  fHistograms->Add(MEHistograms);

  // hMETracks = new TObjArray();
  // hMETracks->SetName("TrackArray");
  // hMETracks->SetOwner(kTRUE);
  // hMEPhotons = new TObjArray();
  // hMEPhotons->SetName("PhotonArray");
  // hMEPhotons->SetOwner(kTRUE);
  // hMEPions = new TObjArray();
  // hMEPions->SetName("PionArray");
  // hMEPions->SetOwner(kTRUE);

  // MEHistograms->Add(hMETracks);
  // MEHistograms->Add(hMEPions);
  // MEHistograms->Add(hMEPhotons);

  hMEvents = new TH2I("hMEvents", "Nevents vs centrality vertexz",
					  fAxisZ.GetNbins(), fAxisZ.GetBinLowEdge(1), fAxisZ.GetBinUpEdge(fAxisZ.GetNbins()),
 					  fAxisCent.GetNbins(), fAxisCent.GetBinLowEdge(1), fAxisCent.GetBinUpEdge(fAxisCent.GetNbins()));
  hMEvents->GetYaxis()->Set(fAxisCent.GetNbins(), fAxisCent.GetXbins()->GetArray());
  MEHistograms->Add(hMEvents);


  // TList axesList;
  // axesList.AddAt(&GetAxisEta(), 0);
  // axesList.AddAt(&GetAxisPhi(), 1);
  // axesList.AddAt(&GetAxistPt(), 2);
  // axesList.SetOwner(kFALSE);
  
  // TList piAxesList;
  // piAxesList.AddAt(&GetAxisEta(), 0);
  // piAxesList.AddAt(&GetAxisPhi(), 1);
  // piAxesList.AddAt(&GetAxistPt(), 2);
  // piAxesList.AddAt(&GetAxisPiMass(), 3);
  // piAxesList.SetOwner(kFALSE);


  TList * outAxesList = new TList();
  outAxesList->Add(&fAxisCent);
  outAxesList->Add(&fAxisZ);
  fHistograms->Add(outAxesList);

  // for(Int_t iz = 0; iz < fAxisZ.GetNbins(); iz++) {
  // 	TObjArray * trackArray = new TObjArray();
  // 	trackArray->SetName(Form("METracks_%d", iz));
  // 	trackArray->SetOwner(kTRUE);
  // 	TObjArray * photonArray = new TObjArray();
  // 	photonArray->SetName(Form("MEPhotons_%d", iz));
  // 	photonArray->SetOwner(kTRUE);
  // 	TObjArray * pionArray = new TObjArray();
  // 	pionArray->SetName(Form("MEPions_%d", iz));
  // 	pionArray->SetOwner(kTRUE);


  // 	hMEPions->AddAt(pionArray, iz);
  // 	hMETracks->AddAt(trackArray, iz);
  // 	hMEPhotons->AddAt(photonArray, iz);

  // 	for(Int_t ic = 0; ic < fAxisCent.GetNbins(); ic++) {

  // 	  TString nameString = Form("%d_%d", iz, ic);
  // 	  TString titleString = Form("%f < Z < %f ... %f cent %f", 
  // 								 fAxisZ.GetBinLowEdge(iz+1), fAxisZ.GetBinUpEdge(iz+1), 
  // 								 fAxisCent.GetBinLowEdge(ic+1), fAxisCent.GetBinUpEdge(ic+1));


  // 	  THnSparseF * trackHistogram = CreateSparse(Form("tracks_%s", nameString.Data()), 
  // 												 Form("tracks %s", titleString.Data()), &axesList );
  // 	  trackArray->AddAt(trackHistogram, ic);

  // 	  THnSparseF * photonHistogram = CreateSparse(Form("photons_%s", nameString.Data()), 
  // 												 Form("photons %s", titleString.Data()), &axesList );
  // 	  photonArray->AddAt(photonHistogram, ic);

  // 	  THnSparseF * pionHistogram = CreateSparse(Form("pions_%s", nameString.Data()), 
  // 												 Form("pions %s", titleString.Data()), &piAxesList );
  // 	  pionArray->AddAt(pionHistogram, ic);
  // 	}
  // }

  PostData(1, fHistograms);
  PostData(2, fHistoGamma);
  PostData(3, fHistoPion);

}

///________________________________________________________________________
THnSparseF * AliAnalysisTaskdPhi::CreateSparse(TString nameString, TString titleString, TList * axesList) {
  ///Create sparse
  const Int_t dim = axesList->GetSize();

  TAxis * axes[dim];
  Int_t bins[dim];
  Double_t min[dim];
  Double_t max[dim];

  for(Int_t i = 0; i<dim; i++) {
	TAxis * axis = dynamic_cast<TAxis*>(axesList->At(i));
	if(axis) {
	  axes[i] = axis;
  	} else {
	  cout << "AliAnalysisTaskdPhi::CreateSparse: Error error, all the axes are not present in axis list" << endl;
	  return NULL;
	}
  }

  for(Int_t i = 0; i<dim; i++) {
	bins[i] = axes[i]->GetNbins(); 
	min[i] = axes[i]->GetBinLowEdge(1);
	max[i] = axes[i]->GetBinUpEdge(axes[i]->GetNbins());
  }

  THnSparseF * sparse = new THnSparseF(Form("METracks_%s", nameString.Data()), 
											   Form("tracks %s", titleString.Data()), 
											   dim, bins, min, max);
  
  for(Int_t i = 0; i<dim; i++) {
	sparse->GetAxis(i)->SetNameTitle(axes[i]->GetName(), axes[i]->GetTitle() );
	if(axes[i]->GetXbins()->GetSize() > 0) {
	  sparse->SetBinEdges(i, axes[i]->GetXbins()->GetArray() );
	}
  }

  return sparse;
}

//________________________________________________________________________
void AliAnalysisTaskdPhi::UserExec(Option_t *) {
  ///User exec. 

  //if(! fV0Filter->EventIsSelected(fInputEvent)) return;


  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) {
	cout << "cout no input event handler"<<endl;
	return;
  }

  
  if ( fV0Filter && !fV0Filter->GetPIDResponse() ) {
	if ( inputHandler->GetPIDResponse() ){
	  fV0Filter->SetPIDResponse( inputHandler->GetPIDResponse() );
	} else {
	  
	  //AOD case
	  if (isAOD){
		if (!fV0Filter->GetPIDResponse()){
		  fV0Filter->InitAODpidUtil(1);
		}
	  }
	}
  }

  Double_t centrality = 0.0;
  Double_t eventPlane = 0.0;
  Double_t vertexz = fInputEvent->GetPrimaryVertex()->GetZ();
  if(isAOD) {
    AliAODHeader * header = static_cast<AliAODHeader*>(fInputEvent->GetHeader());
	centrality = header->GetCentrality();
	eventPlane = header->GetEventplane();
  } else {
	centrality = static_cast<AliESDEvent*>(fInputEvent)->GetCentrality()->GetCentralityPercentile("kV0M");
	eventPlane = fInputEvent->GetEventplane()->GetEventplane("Q");
  }


  const Int_t centBin = GetBin(fAxisCent, centrality);
  const Int_t vertexBin = GetBin(fAxisZ, vertexz);


  if(DebugLevel () > 4) {
	cout << "centrality: " << centrality <<  " " << GetBin(fAxisCent, centrality) << endl;
	cout << "vertexz: " << vertexz <<  " " << GetBin(fAxisZ, vertexz) << endl;
	cout << "eventPlane: " << eventPlane <<  " " << endl;
  }



  if(centBin < 0 || vertexBin < 0) {
	AliError("bin out of range");
	return;
  }

  fGammas->Clear();
  fPions->Clear();

  TClonesArray * aodGammas = GetConversionGammas(isAOD);
  if(!aodGammas) {
	AliError("no aod gammas found!");
	return;
  }

  if(aodGammas->GetEntriesFast() > 0) {
	if( static_cast<AliAODConversionParticle*>(aodGammas->At(0))->GetLabel(0) == fL1 && 
		static_cast<AliAODConversionParticle*>(aodGammas->At(0))->GetLabel(1) == fL2 
		) {
	  return;
	}
	fL1 = static_cast<AliAODConversionParticle*>(aodGammas->At(0))->GetLabel(0);
	fL2 = static_cast<AliAODConversionParticle*>(aodGammas->At(0))->GetLabel(1);
	//cout << aodGammas->GetEntriesFast() << " " << fInputEvent->GetNumberOfTracks() << "c" << endl;
  }

  if(DebugLevel() > 1) printf("Number of conversion gammas %d \n", aodGammas->GetEntriesFast());
  for(Int_t ig = 0; ig < aodGammas->GetEntriesFast(); ig++) {
    AliAODConversionPhoton * photon = dynamic_cast<AliAODConversionPhoton*>(aodGammas->At(ig));
    
    if(!photon) continue;
    if(!fV0Filter || fV0Filter->PhotonIsSelected(static_cast<AliConversionPhotonBase*>(photon), fInputEvent)) {
      fGammas->Add(static_cast<TObject*>(photon));
    }
  }
  
  if(DebugLevel() > 4) printf("Number of accepted gammas %d \n", fGammas->GetEntriesFast());
  hMEvents->Fill(vertexz, centrality);
  
  
  
  ///create track array
  TObjArray tracks;
  const Double_t etalim[2] = { fAxisEta.GetBinLowEdge(1), fAxisEta.GetBinUpEdge(fAxisEta.GetNbins())};
  for(Int_t iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); iTrack++) {

	AliVTrack * track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
	if(track->Pt() < fAxiscPt.GetBinLowEdge(1) ) continue;
	if(track->Eta() < etalim[0] || track->Eta() > etalim[1]) continue;

	
	if(fTrackCuts && fTrackCuts->IsSelected((track))) {
	  tracks.Add(track);
	}
  }

  Process(fGammas, &tracks, vertexBin, centBin);

  PostData(1, fHistograms);
  PostData(2, fHistoGamma);
  PostData(3, fHistoPion);
  
}


//________________________________________________________________________
void AliAnalysisTaskdPhi::Process(TObjArray * gammas, TObjArray * tracks, Int_t vertexBin, Int_t centBin) {
  ///Process stuff

  if(DebugLevel() > 4) printf("Number of accepted gammas, tracks %d  %d \n", gammas->GetEntriesFast(), tracks->GetEntriesFast());
 
  AliAnaConvCorrBase * gCorr = GetCorrObject(vertexBin, centBin, fPhotonCorr);
  AliAnaConvCorrBase * piCorr = GetCorrObject(vertexBin, centBin, fPionCorr);
  
  if(!gCorr || !piCorr) {
	AliError("corr object missing");
	return;
  }

  for(Int_t i1 = 0; i1 < gammas->GetEntriesFast(); i1++) {
   	AliAODConversionPhoton * ph1 = static_cast<AliAODConversionPhoton*>(gammas->UncheckedAt(i1));
	Int_t tIDs[4] = {ph1->GetLabel(0), ph1->GetLabel(1), -1, -1};

	Int_t leading = fIsoAna->IsLeading(static_cast<AliAODConversionParticle*>(ph1), tracks, tIDs);
	if(ph1->Pt() > fAxistPt.GetBinLowEdge(1)) {
	  gCorr->CorrelateWithTracks( static_cast<AliAODConversionParticle*>(ph1), tracks, tIDs, leading);
	}
	for(Int_t i2 = 0; i2 < i1; i2++) {
	  AliAODConversionPhoton * ph2 = static_cast<AliAODConversionPhoton*>(gammas->UncheckedAt(i2));
	  
	  if( ph2->GetTrackLabelPositive()==ph1->GetTrackLabelPositive() 
		  || ph2->GetTrackLabelNegative()==ph1->GetTrackLabelNegative()
		  || ph2->GetTrackLabelNegative()==ph1->GetTrackLabelPositive()
		  || ph2->GetTrackLabelPositive()==ph1->GetTrackLabelNegative()) {
		continue;
	  }

	  AliAODConversionMother * pion = new AliAODConversionMother(ph1, ph2);
	  pion->SetLabels(i1, i2);
	  
	  if(!fV0Filter || fV0Filter->MesonIsSelected(pion, kTRUE) ) {
		

		Int_t leadingpi = fIsoAna->IsLeading(static_cast<AliAODConversionParticle*>(pion), tracks, tIDs);
		//piCorr->FillTriggerCounters(pion, leadingpi);
		tIDs[2] = ph2->GetLabel(0);
		tIDs[3] = ph2->GetLabel(1);
		if(pion->Pt() > fAxistPt.GetBinLowEdge(1) && 
		   pion->M() > fAxisPiM.GetBinLowEdge(1) && 
		   pion->M() < fAxisPiM.GetBinUpEdge(fAxisPiM.GetNbins())) {
		  piCorr->CorrelateWithTracks(pion, tracks, tIDs, leadingpi);
		}
	  }
	}
  }
}

//________________________________________________________________________
void AliAnalysisTaskdPhi::Terminate(Option_t *) {
 
  // Draw result to the screen
  // Called once at the end of the query
}

//________________________________________________________________________
TClonesArray * AliAnalysisTaskdPhi::GetConversionGammas(Bool_t isAOD) {
  

  if(isAOD) {

	TClonesArray * gammas = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fDeltaAODBranchName.Data()));
	if(gammas) {
	  return gammas;
	}

	FindDeltaAODBranchName(fInputEvent);
	gammas = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fDeltaAODBranchName.Data()));
	return gammas;

  } else {
	TClonesArray * gammas = dynamic_cast<TClonesArray*>(GetInputData(1));
	return gammas;
  }

}

//________________________________________________________________________
void AliAnalysisTaskdPhi::FindDeltaAODBranchName(AliVEvent * event){
  ///Find aod branch
  TList *list=event->GetList();
  for(Int_t ii=0;ii<list->GetEntries();ii++){
	TString name((list->At(ii))->GetName());
	if(name.BeginsWith("GammaConv")&&name.EndsWith("gamma")){
	  fDeltaAODBranchName=name;
	  AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
	  return;
	}
  }
}
  

