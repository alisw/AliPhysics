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

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations
// author: D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <AliAnalysisManager.h>
#include <AliQnCorrectionsManager.h>
#include <AliJBaseTrack.h>
#include "AliJFlowCMBTask.h" 
static double decAcc[AliJFlowCMBTask::D_COUNT][2] = {
	{-0.8,0.8},
	{-1.5,-0.8},
	{0.8,1.5},
	{2.8,5.1},   // V0A
	{-3.7,-1.7}, // V0C
	{2.8,5.1},   // V0P+ need to do it manually
	{2.5,5.1} // Virtual dector +-
};
static const char *pdetn[4] = {"TPC","VZERO","VZEROA","VZEROC"};
static int newDetID[4] = {AliJFlowCMBTask::D_TPC,AliJFlowCMBTask::D_V0P,AliJFlowCMBTask::D_V0A,AliJFlowCMBTask::D_V0C};

//______________________________________________________________________________
AliJFlowCMBTask::AliJFlowCMBTask() :   
	AliAnalysisTaskSE("JFlowBaseTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fFlowVectorTask(NULL),
	fIsMC(kTRUE),
	fhistos(NULL),
	fCBin(-1),
	fOutput(NULL)
{
}

//______________________________________________________________________________
AliJFlowCMBTask::AliJFlowCMBTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fFlowVectorTask(NULL),
	fIsMC(kTRUE),
	fhistos(NULL),
	fCBin(-1),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJFlowCMBTask Constructor ----");
	for(uint i=0;i<D_COUNT;i++) {
		for(uint j=0;j<2;j++) fEventPlaneALICE[i][j] =-999;
	}
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJFlowCMBTask::AliJFlowCMBTask(const AliJFlowCMBTask& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fFlowVectorTask(ap.fFlowVectorTask),
	fIsMC(ap.fIsMC),
	fhistos(ap.fhistos),
	fCBin(ap.fCBin),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJFlowCMBTask COPY ----");

}

//_____________________________________________________________________________
AliJFlowCMBTask& AliJFlowCMBTask::operator = (const AliJFlowCMBTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJFlowCMBTask operator= ----");
	this->~AliJFlowCMBTask();
	new(this) AliJFlowCMBTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJFlowCMBTask::~AliJFlowCMBTask()
{
	// destructor 
	delete fOutput;
	delete fhistos;

}

//________________________________________________________________________

void AliJFlowCMBTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJFlowCMBTask::UserCreateOutPutData() \n");
	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
	if(!fIsMC) fFlowVectorTask = (AliAnalysisTaskFlowVectorCorrections*)(man->GetTask("FlowQnVectorCorrections"));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	fhistos = new AliJFlowHistos();
	fhistos->CreateEventTrackHistos();

	fhistos->fHMG->Print();


	PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJFlowCMBTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJFlowCMBTask Exec-------"<<endl;
	if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	fCBin = AliJFlowHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
	if(fCBin == -1)
		return;
	if(fIsMC) {
		TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();
		CalculateEventPlane(fInputList);
	} else {
		AliQnCorrectionsManager *fFlowVectorMgr = fFlowVectorTask->GetAliQnCorrectionsManager();
		const AliQnCorrectionsQnVector *fQnVector;
		for(UInt_t di = 0; di < sizeof(pdetn)/sizeof(pdetn[0]); ++di) {
			for(int iH=2;iH<=3;iH++) {
				fQnVector = fFlowVectorMgr->GetDetectorQnVector(pdetn[di]);
				if(fQnVector) fEventPlaneALICE[newDetID[di]][iH-2] = fQnVector->EventPlane(iH);
			}
		}
		for(UInt_t di = 0; di < sizeof(pdetn)/sizeof(pdetn[0]); ++di) {
			for(int iH=2;iH<=3;iH++) {		
				double EPref = fEventPlaneALICE[D_V0A][iH-2];
				double EP = fEventPlaneALICE[newDetID[di]][iH-2];
				fhistos->fhEPCorrInHar[fCBin][newDetID[di]][iH-2]->Fill( EP-EPref );
			}
		}
	}
}

//______________________________________________________________________________
void AliJFlowCMBTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJFlowCMBTask::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJFlowCMBTask Analysis DONE !!"<<endl; 
}

//______________________________________________________________________________
void AliJFlowCMBTask::CalculateEventPlane(TClonesArray *inList) {

	int NtracksDEC[AliJFlowCMBTask::D_COUNT]; // Num of particle in each dector
	int noTracks = inList->GetEntries();
	for(int is = 0; is < AliJFlowCMBTask::D_COUNT; is++)  { 
		NtracksDEC[is] = 0;
		for(int iH=2;iH<=3;iH++) { QvectorsEP[is][iH-2] = TComplex(0,0); }
	}
	for(int itrack=0;itrack<noTracks; itrack++){
		AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(itrack);
		double phi = trk->Phi();
		double eta = trk->Eta();
		for(int is = 0; is < AliJFlowCMBTask::D_COUNT; is++){
			if(is == AliJFlowCMBTask::D_VIRT) {
				if(decAcc[is][0]<eta && decAcc[is][1]>eta) {
					for(int iH=2;iH<=3;iH++) { QvectorsEP[is][iH-2] += TComplex(TMath::Cos(iH*phi),TMath::Sin(iH*phi));}
					NtracksDEC[is]++;
					fhistos->fh_eta[fCBin][is]->Fill(eta);
					fhistos->fh_phi[fCBin][is]->Fill(phi);
					fhistos->fh_pt[fCBin][is]->Fill(trk->Pt());
				}
			} else if(is == AliJFlowCMBTask::D_V0P) {
				if((decAcc[AliJFlowCMBTask::D_V0A][0]<eta && decAcc[AliJFlowCMBTask::D_V0A][1]>eta) || 
						(decAcc[AliJFlowCMBTask::D_V0C][0]<eta && decAcc[AliJFlowCMBTask::D_V0C][1]>eta) ) {
					for(int iH=2;iH<=3;iH++) { QvectorsEP[is][iH-2] += TComplex(TMath::Cos(iH*phi),TMath::Sin(iH*phi));}
					NtracksDEC[is]++;
					fhistos->fh_eta[fCBin][is]->Fill(eta);
					fhistos->fh_phi[fCBin][is]->Fill(phi);
					fhistos->fh_pt[fCBin][is]->Fill(trk->Pt());
				}
			} else {
				if(decAcc[is][0]<eta && decAcc[is][1]>eta) {
					for(int iH=2;iH<=3;iH++) { QvectorsEP[is][iH-2] += TComplex(TMath::Cos(iH*phi),TMath::Sin(iH*phi));}
					NtracksDEC[is]++;
					fhistos->fh_eta[fCBin][is]->Fill(eta);
					fhistos->fh_phi[fCBin][is]->Fill(phi);
					fhistos->fh_pt[fCBin][is]->Fill(trk->Pt());
				}
			}
		}
	}

	for(int is = 0; is < AliJFlowCMBTask::D_COUNT; is++){
		for(int iH=2;iH<=3;iH++) {		
			QvectorsEP[is][iH-2] /= (double)NtracksDEC[is];
			QvectorsEP[is][iH-2] = QvectorsEP[is][iH-2]/TComplex::Abs(QvectorsEP[is][iH-2]);
			fhistos->fh_EP[fCBin][is][iH-2]->Fill(QvectorsEP[is][iH-2].Theta()/double(iH));
		}
	}

	// Calculation resolution
	for(int is = 0; is < AliJFlowCMBTask::D_COUNT; is++){
		for(int iH=2;iH<=3;iH++) {		
			double EPref = QvectorsEP[D_V0A][iH-2].Theta()/double(iH);
			double EP = QvectorsEP[is][iH-2].Theta()/double(iH);
			fhistos->fhEPCorrInHar[fCBin][is][iH-2]->Fill( EP-EPref );
			//fhistos->fhEPCorr2D[fCBin][is][iH-2]->Fill(EPref,EP);
		}
	}

}
