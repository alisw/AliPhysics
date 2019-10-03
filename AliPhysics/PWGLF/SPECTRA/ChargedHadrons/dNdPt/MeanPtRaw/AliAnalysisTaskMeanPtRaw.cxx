//------------------------------------------------------------------------------
// AliAnalysisTaskMeanPtRaw 
// 
// generates one TH2F pT vs Mult raw
// 
// made by P. Luettig
// last modified 2015/05/08
//------------------------------------------------------------------------------


#include <RVersion.h>
#include "TChain.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"



#include "AliAnalysisTaskMeanPtRaw.h"

ClassImp(AliAnalysisTaskMeanPtRaw)
//________________________________________________________________________
AliAnalysisTaskMeanPtRaw::AliAnalysisTaskMeanPtRaw(const char *name) 
: AliAnalysisTaskSE(name), 
fESD(0), 
fOutputList(0),
fPtVsMultRaw(0),
fTrackCutName(0),
fESDTrackCuts(0),
fTrackCutNames(0),
fNTrackCuts(0),
fMaxVertexZ(0),
fEventCount(0)
{
  // Constructor
  fNTrackCuts = 0;
  // Define input and output slots here
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskMeanPtRaw::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  OpenFile(1, "RECREATE");
  
  fOutputList = new TList();
  fOutputList->SetOwner();
  
  Double_t binsPtDefault[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  
  Int_t binsMultPtCuts[3] =   {  3500,   81, fNTrackCuts};
  Double_t minMultPtCuts[3] = {  -0.5,    0,           0};
  Double_t maxMultPtCuts[3] = {3499.5,  200, (Double_t)fNTrackCuts};
  
  fPtVsMultRaw = new THnSparseF("fPtVsMultRaw", "raw <pT> vs mult",3, binsMultPtCuts, minMultPtCuts, maxMultPtCuts);
  fPtVsMultRaw->SetBinEdges(1, binsPtDefault);
  fPtVsMultRaw->GetAxis(0)->SetTitle("n_{acc} (raw)");
  fPtVsMultRaw->GetAxis(1)->SetTitle("p_{T} (raw) (GeV/c)");
  fPtVsMultRaw->GetAxis(2)->SetTitle("cut setting");
  
  fTrackCutName = new TH1I("fTrackCutName","fTrackCutName",10,0,10);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  fTrackCutName->SetBit(TH1::kCanRebin);
#endif
  fTrackCutName->Sumw2();
  
  for(Int_t i = 0; i < fTrackCutNames.size(); i++)
  {
	fTrackCutName->Fill(fTrackCutNames.at(i), 1);
  }
  
  // Add Histos, Profiles etc to List
  fOutputList->Add(fPtVsMultRaw);
  fOutputList->Add(fTrackCutName);
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMeanPtRaw::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  Int_t *nacc = new Int_t[fNTrackCuts];
  for(Int_t iTrackCuts = 0; iTrackCuts < fNTrackCuts; iTrackCuts++) { nacc[iTrackCuts] = 0; }
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) { Printf("ERROR: Could not get ESDInputHandler\n"); }
  else fESD = (AliESDEvent*)esdH->GetEvent();
  if (!fESD) { Printf("ERROR: fESD not available"); return; } 
  fEventCount++;
  
  if( (fESD->GetNumberOfTracks() == 0 ) )  
  {
	return;
  }
  
  if(!fESDTrackCuts.at(0)) 
  {
	Printf("ERROR: No esd track cuts available");
	return;
  }
  
  const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
  if( vtxESD->GetNContributors() < 1 ) 
  {
	return;
  }
  
  
  // Event cut on the z-Position of the vertex
  if(vtxESD->GetZ() > fMaxVertexZ || vtxESD->GetZ() < (-1.*fMaxVertexZ)) 
  {
	//     Printf("VertexZ out of range, Zv = %f",vtxESD->GetZv());
	return;
  }
  
  
  
  // loop over all reconstructed tracks
  // to get number of accepted tracks
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) 
  {
	AliESDtrack* track = fESD->GetTrack(iTracks);
	
	if (!track) 
	{
	  Printf("ERROR: Could not receive track %d", iTracks);
	  continue;
	}
	
	for(Int_t iTrackCuts = 0; iTrackCuts < fNTrackCuts; iTrackCuts++)
	{
	  if( fESDTrackCuts.at(iTrackCuts)->AcceptTrack(track) ) nacc[iTrackCuts]++;
	}
  } //track reconstruced loop
  
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) 
  {
	AliESDtrack* track = fESD->GetTrack(iTracks);
	
	if (!track) 
	{
	  Printf("ERROR: Could not receive track %d", iTracks);
	  continue;
	}
	
	Double_t pt = track->Pt();
	for(Int_t iTrackCuts = 0; iTrackCuts < fNTrackCuts; iTrackCuts++)
	{
	  Double_t dNaccPtCut[3] = {(Double_t)nacc[iTrackCuts], pt, (Double_t)iTrackCuts};
	  if( fESDTrackCuts.at(iTrackCuts)->AcceptTrack(track) ) fPtVsMultRaw->Fill(dNaccPtCut);
	}
  } //track reconstruced loop
  
  PostData(1, fOutputList);
  
  delete [] nacc;
}


void AliAnalysisTaskMeanPtRaw::AddAliESDtrackCut(AliESDtrackCuts *esdTrackCuts, const char *name) 
{ 
  fTrackCutNames.push_back(name);
  
  fESDTrackCuts.push_back(new AliESDtrackCuts(*esdTrackCuts));
  
  fNTrackCuts++;
}

//________________________________________________________________________
void AliAnalysisTaskMeanPtRaw::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

