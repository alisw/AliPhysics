/*************************************************************************
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
//

// Author:
// 

#include "TFile.h"

#include "AliESDInputHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliESDAD.h"
#include "AliCDMesonUtils.h"

#include "AliAnalysisTaskCEPQA.h"

//______________________________________________________________________________

AliAnalysisTaskCEPQA::AliAnalysisTaskCEPQA(const char* name):
	AliAnalysisTaskSE(name)
  , ff(0x0)
	, fESDEvent(0x0)
  , fTrigger()
	, fTree(0x0)
	, fList(0x0)
  , fHistStatus(0x0)
	, fRunNum(0x0)
	, fEvCounter(0x0)
	, fisMBOR(0x0)
	, fisMBAND(0x0)
	, fisPileup(0x0)
	, fSPDmul(0x0)
	, fVtype(0x0)
	, fVposx(0x0)
	, fVposy(0x0)
	, fVposz(0x0)
	, fisV0A(0x0)
	, fisV0C(0x0)
	, fFMDAnumHC(0x0)
	, fFMDCnumHC(0x0)
	, fADAmul(0x0)
	, fADCmul(0x0)
	, fADATCh(0x0)
	, fADCTCh(0x0)
	, fADATD(0x0)
	, fADCTD(0x0)
	, fZDCTime(0x0)
	, fZNAEne(0x0)
	, fZNCEne(0x0)
	, fZPAEne(0x0)
	, fZPCEne(0x0)
  , fisNG(0x0)
  , fisSG(0x0)
  , fisDG(0x0)
{	

	// create output slot 1
	DefineOutput(1, TTree::Class());
	DefineOutput(2, TList::Class());

}

//______________________________________________________________________________

AliAnalysisTaskCEPQA::AliAnalysisTaskCEPQA():
	AliAnalysisTaskSE()
  , ff(0x0)
	, fESDEvent(0x0)
  , fTrigger()
	, fTree(0x0)
	, fList(0x0)
  , fHistStatus(0x0)
	, fRunNum(0x0)
	, fEvCounter(0x0)
	, fisMBOR(0x0)
	, fisMBAND(0x0)
	, fisPileup(0x0)
	, fSPDmul(0x0)
	, fVtype(0x0)
	, fVposx(0x0)
	, fVposy(0x0)
	, fVposz(0x0)
	, fisV0A(0x0)
	, fisV0C(0x0)
	, fFMDAnumHC(0x0)
	, fFMDCnumHC(0x0)
	, fADAmul(0x0)
	, fADCmul(0x0)
	, fADATCh(0x0)
	, fADCTCh(0x0)
	, fADATD(0x0)
	, fADCTD(0x0)
	, fZDCTime(0x0)
	, fZNAEne(0x0)
	, fZNCEne(0x0)
	, fZPAEne(0x0)
	, fZPCEne(0x0)
  , fisNG(0x0)
  , fisSG(0x0)
  , fisDG(0x0)
{	

}

//______________________________________________________________________________

AliAnalysisTaskCEPQA::~AliAnalysisTaskCEPQA()
{

	// output tree
	if (fTree) {
		delete fTree;
		fTree = 0x0;
	}
	
	if (fList) {
		delete fList;
		fList = 0x0;
	}
	
	if (fHistStatus) delete fHistStatus;

}

//______________________________________________________________________________

void AliAnalysisTaskCEPQA::UserCreateOutputObjects()
{
	
	// the output consists of an ntuple and histograms
  // the ntuple s filled only with information from selected events
  // histograms are filled with all events
  
  // define output tree
	OpenFile(1);
	fTree = new TTree("CEP","CEPQA");
	
	// for QA
	fTree->Branch("RunNum",&fRunNum);
	fTree->Branch("EventNum",&fEvCounter);
	fTree->Branch("isMBor",&fisMBOR);
	fTree->Branch("isMBand",&fisMBAND);
	fTree->Branch("isPUp",&fisPileup);
	fTree->Branch("Vtype",&fVtype);
	fTree->Branch("Vposx",&fVposx);
	fTree->Branch("Vposy",&fVposy);
	fTree->Branch("Vposz",&fVposz);
	fTree->Branch("SPD",&fSPDmul);
	fTree->Branch("V0A",&fisV0A);
	fTree->Branch("V0C",&fisV0C);
	fTree->Branch("FMDAnumHC",&fFMDAnumHC);
	fTree->Branch("FMDCnumHC",&fFMDCnumHC);
	fTree->Branch("ADA",&fADAmul);
	fTree->Branch("ADC",&fADCmul);
	fTree->Branch("ADATCh",&fADATCh);
	fTree->Branch("ADCTCh",&fADCTCh);
	fTree->Branch("ADATD",&fADATD);
	fTree->Branch("ADCTD",&fADCTD);
	fTree->Branch("ZDCtime",&fZDCTime);
	fTree->Branch("ZNA",&fZNAEne);
	fTree->Branch("ZNC",&fZNCEne);
	fTree->Branch("ZPA",&fZPAEne);
	fTree->Branch("ZPC",&fZPCEne);
	fTree->Branch("isNG",&fisNG);
	fTree->Branch("isSG",&fisSG);
	fTree->Branch("isDG",&fisDG);
  
  // create histograms
	fList = new TList();
	fList->SetOwner();

  Int_t nlabs = 16;
  fHistStatus = new TH1D("fHistStatus","Status of event",nlabs+1,0,nlabs+1);
  const char *labels[] = {"Event","MBor","MBand","Pileup","V0A","V0C","SPD",
    "FMDA","FMDC","ADA","ADC","ZDNA","ZDNC","NoGap","SingleGap","DoubleGap"};
	for (Int_t ind=0; ind < nlabs; ind++)
    fHistStatus->GetXaxis()->SetBinLabel(ind+1,labels[ind]);
  fList->Add(fHistStatus);
  
  // creates the trigger monitoring histograms
  //ff = new TFile("axibixi.root","RECREATE");
  //TString triggerClass = "trigger_histograms_";
  //gDirectory->mkdir(triggerClass);
  //gDirectory->cd(triggerClass);
  fTrigger = new AliTriggerAnalysis();
  //fTrigger->EnableHistograms();
  
  PostOutputs();

}

//______________________________________________________________________________

void AliAnalysisTaskCEPQA::UserExec(Option_t *)
{

	// check if input event is available
	if (!CheckInput()) {
		PostOutputs();
    return;
	}
	
  // fill the trigger monitoring histograms
  // fTrigger->PrintTriggerClasses();
  // fTrigger->FillHistograms(fESDEvent);

	// get the run number
	fRunNum = fESDEvent->GetRunNumber();
	
	// increment event number
	fEvCounter++;
	
	// is there an event number?
	//printf("\n<I - TaskCEPQA> This is event %i in the file\n",
	//	fESDEvent->GetEventNumberInFile());

	// analyse the trigger
	TString ftc = fESDEvent->GetFiredTriggerClasses();
	// printf("<I - TaskCEPQA> FiredTriggerClasses:\n%s\n",ftc.Data());

	fisMBOR  = (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kMB1)) ? kTRUE : kFALSE;
  fisMBAND = (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0AND)) ? kTRUE : kFALSE;

	// check for pileup
	fisPileup = fESDEvent->IsPileupFromSPD
	(
		2,		// minContributors, default = 3
		0.8,	// minZdist, default = 0.8
		3., 	// nSigmaZdist, default = 3.
		2., 	// nSigmaDiamXY, default = 2.
		5.		// nSigmaDiamZ, default = 5.
	);
	
	// check the vertex information
	CheckVertex();

	// analyze the different detectors
	// SPD
	CheckSPDactive();
	
	// V0 (V0A and V0C)
	CheckV0active();
	
	// FMD
	CheckFMDactive();
	
	// AD (ADA and ADC)
	CheckADactive();

	// ZDC
	// ZDN (ZDNA and ZDNC)
	CheckZDNactive();
	
	// ZNP (ZPA and ZPC)
	CheckZDPactive();
  
  // is NG, SG, or DG?
  CheckGapType();

	// fill histogram Status
                    fHistStatus->Fill( 0.5);
  if (fisMBOR)      fHistStatus->Fill( 1.5);
	if (fisMBAND)     fHistStatus->Fill( 2.5);
  if (fisPileup)    fHistStatus->Fill( 3.5);
  if (fisV0A)       fHistStatus->Fill( 4.5);
  if (fisV0C)       fHistStatus->Fill( 5.5);
  if (fSPDmul>0)    fHistStatus->Fill( 6.5);
  if (fFMDAnumHC>0) fHistStatus->Fill( 7.5);
  if (fFMDCnumHC>0) fHistStatus->Fill( 8.5);
  if (fADAmul>0)    fHistStatus->Fill( 9.5);
  if (fADCmul>0)    fHistStatus->Fill(10.5);
  // ZDNA           fHistStatus->Fill(11.5);
  // ZDNC           fHistStatus->Fill(12.5);
  if (fisNG)        fHistStatus->Fill(13.5);
  if (fisSG)        fHistStatus->Fill(14.5);
  if (fisDG)        fHistStatus->Fill(15.5);
  
  // apply event filter
  if (!FilterEvent()) {
    PostOutputs();
    return;
  }
  
  // apply track filters
  //TList new goodtracks();
  //Int_t nGoodTracks = FilterTracks(goodtracks);
  
  // 2*n track events, with n=1,2,3
  
  fTree->Fill();
	PostOutputs();

  return;

}

//______________________________________________________________________________
// retrieve the vertex information
//
// there are obviously two sources of vertex which are considered
// a default and SPD
// the default is the best possible estimate and includes information from
// the ITS and TPC, if that does not exists, then the vertex estimate from
// SPD only is used
//
// fVposx: x-position of vertex
// fVposy: y-position of vertex
// fVposz: z-position of vertex
// Vtype:
//	-1: no vertex information
//	 1: vertex from SPD only
//	 2: vertex from ITS+TPC

void AliAnalysisTaskCEPQA::CheckVertex()
{

	// initialisations
	fVtype = 0;
	fVposx = -999.;
	fVposy = -999.;
	fVposz = -999.;
	
	// ITS + TPC vertex
	const AliESDVertex *vertex = fESDEvent->GetPrimaryVertexTracks();
	if (vertex) {
	
		// check also the number of tracklets/tracks used for the estimate of the
		// vertex
		if (vertex->GetNContributors() <1) {
		
			// SPD vertex
			vertex = fESDEvent->GetPrimaryVertexSPD();
			if (vertex->GetNContributors() <1) {
				fVtype = -1;
			} else fVtype = 1;
		
		} else fVtype = 2;
	
	} else fVtype  = -1;

	if (fVtype > 0) {
		fVposx = vertex->GetX();
		fVposy = vertex->GetY();
		fVposz = vertex->GetZ();
	}
	
}

//______________________________________________________________________________
// compute the hit multiplicity in the SPD detector
// fSPDmul: total number of hits in SPD

void AliAnalysisTaskCEPQA::CheckSPDactive()
{

	const AliMultiplicity *mult = fESDEvent->GetMultiplicity();

	fSPDmul = 0;
	for (Int_t iChipKey=0; iChipKey < 1200; iChipKey++) {
		if(mult->TestFastOrFiredChips(iChipKey)){
			fSPDmul++;
		}
	}
		
}

//______________________________________________________________________________
// checks for activity in either side of the V0 detector system
// fisV0A: yes/no activity in V0A
// fisV0C: yes/no activity in V0C

void AliAnalysisTaskCEPQA::CheckV0active()
{

	AliTriggerAnalysis triggerAnalysis;
	
	fisV0A =
		(triggerAnalysis.V0Trigger(fESDEvent,AliTriggerAnalysis::kASide,kTRUE) ==
		AliTriggerAnalysis::kV0BB);
	fisV0C = 
		(triggerAnalysis.V0Trigger(fESDEvent,AliTriggerAnalysis::kCSide,kTRUE) ==
		AliTriggerAnalysis::kV0BB);

}

//______________________________________________________________________________
// checks for activity in the FMD detector
// fFMDAnumHC: number of hit combinations in FMDA
// fFMDCnumHC: number of hit combinations in FMDA

void AliAnalysisTaskCEPQA::CheckFMDactive()
{

  AliCDMesonUtils::GetMultFMD(fESDEvent,fFMDAnumHC,fFMDCnumHC,NULL);
	
}

//______________________________________________________________________________
// get the multiplicity in the AD detectors
// fADAmul: multiplicity measured with ADA
// fADCmul: multiplicity measured with ADC
// fADATCh: Sum of trigger charge on ADA
// fADCTCh: Sum of trigger charge on ADC
// fADATD: ADA trigger decision
// fADCTD: ADA trigger decision
void AliAnalysisTaskCEPQA::CheckADactive()
{

	AliESDAD *fADData = fESDEvent->GetADData();
	if (!fADData) return;
	
	// multiplicities
	for (Int_t cell = 0; cell<8; cell++) {
		fADAmul += fADData->GetMultiplicityADA(cell);
		fADCmul += fADData->GetMultiplicityADC(cell);
	}
	
	// sum of trigger charge
	fADATCh = (Int_t) fADData->GetTriggerChargeA();
	fADCTCh = (Int_t) fADData->GetTriggerChargeC();

	// trigger decision
	fADATD  = fADData->GetADADecision();
	fADCTD  = fADData->GetADCDecision();

}

//______________________________________________________________________________
// read time and energy information in the ZDN detector
// fZDCTime: time information
// fZNAEne: energy measured in ZNA
// fZNCEne: energy measured in ZNC

void AliAnalysisTaskCEPQA::CheckZDNactive()
{

	AliESDZDC *fZDCData = fESDEvent->GetZDCData();
	if (!fZDCData) return;
	
	fZDCTime = fZDCData->GetZDCTimeSum();
	fZNAEne = fZDCData->GetZDCN1Energy();
	fZNCEne = fZDCData->GetZDCN2Energy();

}

//______________________________________________________________________________
// read energy information in the ZDP detector
// fZPAEne: energy measured in ZPA
// fZPCEne: energy measured in ZPC

void AliAnalysisTaskCEPQA::CheckZDPactive()
{

	AliESDZDC *fZDCData = fESDEvent->GetZDCData();
	if (!fZDCData) return;
	
	fZPAEne = fZDCData->GetZDCP1Energy();
	fZPCEne = fZDCData->GetZDCP2Energy();
	
}

//______________________________________________________________________________
// Gap type can be
// NG: No Gap
// SG: Single Gap
// DG: Double Gap

void AliAnalysisTaskCEPQA::CheckGapType()
{
  
  fisNG = kFALSE;
  fisSG = kFALSE;
  fisDG = kFALSE;
  
  fisSG = (fSPDmul>0) && (fisV0A || fisV0C) && !(fisV0A && fisV0C);
  fisDG = (fSPDmul>0) && !fisV0A && !fisV0C;
  fisNG = !fisSG && !fisDG;
  
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskCEPQA::FilterEvent()
{
  
  Bool_t isaccepted = kTRUE;
  
  // apply minimum bias trigger
  isaccepted = isaccepted && fisMBOR;
  
  // primary vertex cut
  isaccepted = isaccepted && (TMath::Abs(fVposz)<100.);
  
  // pileup
  isaccepted = isaccepted && (!fisPileup);
  return isaccepted;
  
  // number of cluster versus number of tracklet cut
  Double_t cut_slope = 4.;
  Double_t cut_offset = 65.;
  const AliMultiplicity *mult = fESDEvent->GetMultiplicity();
  Int_t ntracklets = mult->GetNumberOfTracklets();
  Int_t ncluster12 =
    fESDEvent->GetNumberOfITSClusters(0)+fESDEvent->GetNumberOfITSClusters(1);

  isaccepted = isaccepted && 
    ( ncluster12 <= (cut_offset+cut_slope*ntracklets) );
    
  // only well reconstructed events  
  // + Martin's specific selections
  // N_pure_ITSSA_track > N_trk -> reject
  // N_tracklets        > N_trk -> reject
  Int_t Ntracks = fESDEvent->GetNumberOfTracks();
  
  Int_t NpureITStracks = 0;
  for (Int_t iTrack=0; iTrack < Ntracks; iTrack++) {

    AliESDtrack* track = fESDEvent->GetTrack(iTrack);
    track->SetESDEvent(fESDEvent);

    if ((track->GetStatus() & AliESDtrack::kITSpureSA) != 0)
    {
      NpureITStracks++;
      continue;
    }
  }
  isaccepted = isaccepted && 
    (Ntracks > ntracklets) && (Ntracks > NpureITStracks);
    

  // Double gap event?
  // isaccepted = isaccepted && fisDG;
  
  return isaccepted;

}

//______________________________________________________________________________

Bool_t AliAnalysisTaskCEPQA::CheckInput()
{
	// General protection
	if (const AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler)) {
	  fESDEvent = (AliESDEvent*)esdH->GetEvent();
	}

	if (!fESDEvent) {
		printf("<E - TaskCEPQA> AliAnalysisTaskCEPQA no valid event!\n");
		return kFALSE;
	}
  
	return kTRUE;

}

//______________________________________________________________________________

void AliAnalysisTaskCEPQA::PostOutputs()
{

	// save output slot 1
	PostData(1,fTree);
	PostData(2,fList);

	return;

}

//______________________________________________________________________________

void AliAnalysisTaskCEPQA::Terminate(Option_t *)
{

  // save the trigger monitoring histograms
  if (ff) {
    TString triggerClass = "trigger_histograms_";
    gDirectory->cd(triggerClass);
    fTrigger->SaveHistograms();
    
    ff->Write();
    ff->Close();
    
  } else {
    printf("<A - AliAnalysisTaskCEPQA> Output file not open!\n");
  }

}

//______________________________________________________________________________
void AliAnalysisTaskCEPQA::Clear(Option_t *)
{


}

//______________________________________________________________________________
