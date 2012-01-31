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
//------------------------------------------------------------------------------
// AlidNdPtEfficiency class. 
//
// a. functionality:
// - fills generic cut histograms
// - generates cuts (selection criteria)
//
// b. data members:
// - generic cut histograms
// - control histograms
//
// Author: J.Otwinowski 18/11/2010 
//------------------------------------------------------------------------------
#include "TH1.h"
#include "TH2.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 
#include "AliTracker.h" 

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AlidNdPtBackgroundCuts.h"
#include "AlidNdPtAnalysis.h"
#include "AliPhysicsSelection.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtEfficiency.h"

using namespace std;

ClassImp(AlidNdPtEfficiency)

//_____________________________________________________________________________
  AlidNdPtEfficiency::AlidNdPtEfficiency(): AlidNdPt(),
  fAnalysisFolder(0),
  fRecMCTrackHistTPCITS(0),
  fRecMCTrackHistITSTPC(0)
{
  // default constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtEfficiency::AlidNdPtEfficiency(Char_t* name, Char_t* title): AlidNdPt(name,title),
  fAnalysisFolder(0),
  fRecMCTrackHistTPCITS(0),
  fRecMCTrackHistITSTPC(0)
{
  // constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtEfficiency::~AlidNdPtEfficiency() {
  // 
  if(fRecMCTrackHistTPCITS) delete fRecMCTrackHistTPCITS; fRecMCTrackHistTPCITS=0;
  if(fRecMCTrackHistITSTPC) delete fRecMCTrackHistITSTPC; fRecMCTrackHistITSTPC=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtEfficiency::Init(){
  //
  // Init histograms
  //
  const Int_t ptNbins = 63; 
  const Double_t ptMin = 0.; 
  const Double_t ptMax = 20.; 

  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0, 20.,25.,30.,35.,40.,50};

  // 
  // THnSparse track histograms
  //

  // TPC -> ITS matching efficiency
  // eta:phi:pt:isPrim:charge:isMatch:isTPC
  Int_t binsRecMCTrackHistTPCITS[7]=  { 30,  90,             ptNbins, 2,   3,  2,  2  };
  Double_t minRecMCTrackHistTPCITS[7]={-1.5, 0.,             ptMin,   0., -1., 0., 0. };
  Double_t maxRecMCTrackHistTPCITS[7]={ 1.5, 2.*TMath::Pi(), ptMax,   2.,  2., 2., 2. };

  fRecMCTrackHistTPCITS = new THnSparseF("fRecMCTrackHistTPCITS","eta:phi:pt:isPrim:charge:isMatch:isTPC",7,binsRecMCTrackHistTPCITS,minRecMCTrackHistTPCITS,maxRecMCTrackHistTPCITS);
  fRecMCTrackHistTPCITS->SetBinEdges(2,binsPt);
  fRecMCTrackHistTPCITS->GetAxis(0)->SetTitle("#eta");
  fRecMCTrackHistTPCITS->GetAxis(1)->SetTitle("#phi (rad)");
  fRecMCTrackHistTPCITS->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fRecMCTrackHistTPCITS->GetAxis(3)->SetTitle("isPrim");
  fRecMCTrackHistTPCITS->GetAxis(4)->SetTitle("charge");
  fRecMCTrackHistTPCITS->GetAxis(5)->SetTitle("isMatch");
  fRecMCTrackHistTPCITS->GetAxis(6)->SetTitle("isTPC");
  fRecMCTrackHistTPCITS->Sumw2();

  // ITS -> TPC matching efficiency
  // eta:phi:pt:isPrim:charge:isMatch
  Int_t binsRecMCTrackHistITSTPC[6]=  { 30,  90,             ptNbins,   2,   3,  2 };
  Double_t minRecMCTrackHistITSTPC[6]={-1.5, 0.,             ptMin,     0., -1., 0 };
  Double_t maxRecMCTrackHistITSTPC[6]={ 1.5, 2.*TMath::Pi(), ptMax,     2.,  2., 2.};

  fRecMCTrackHistITSTPC = new THnSparseF("fRecMCTrackHistITSTPC","eta:phi:pt:isPrim:charge:isMatch",6,binsRecMCTrackHistITSTPC,minRecMCTrackHistITSTPC,maxRecMCTrackHistITSTPC);
  fRecMCTrackHistITSTPC->SetBinEdges(2,binsPt);
  fRecMCTrackHistITSTPC->GetAxis(0)->SetTitle("#eta");
  fRecMCTrackHistITSTPC->GetAxis(1)->SetTitle("#phi (rad)");
  fRecMCTrackHistITSTPC->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fRecMCTrackHistITSTPC->GetAxis(3)->SetTitle("isPrim");
  fRecMCTrackHistITSTPC->GetAxis(4)->SetTitle("charge");
  fRecMCTrackHistITSTPC->GetAxis(5)->SetTitle("isMatch");
  fRecMCTrackHistITSTPC->Sumw2();

  // init output folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");
}

//_____________________________________________________________________________
void AlidNdPtEfficiency::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent)
{
  //
  // Process real and/or simulated events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts      = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts   = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

   

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
   
  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  if(IsUseMCInfo())
  {
    //
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  if(evtCuts->IsRecVertexRequired()) 
  {
     if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
        vtxESD = esdEvent->GetPrimaryVertexTPC();
    }
    else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
      vtxESD = esdEvent->GetPrimaryVertexTracks();
    }
    else {
    	return;
    }
  }
  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  TObjArray *allChargedTracks=0;
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK,isEventTriggered);

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    Bool_t isTPC = kFALSE;
    Bool_t isMatch = kFALSE;

    // TPC -> ITS prolongation efficiency
    for(Int_t iTrack=0; iTrack<entries;++iTrack) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(iTrack);
      if(!track) continue;

      isTPC = kFALSE;

      if(track->Charge()==0) continue;
      if(!track->GetTPCInnerParam()) continue;
      if(!(track->GetStatus()&AliESDtrack::kTPCrefit)) continue;

      // check loose cuts for TPC tracks
      if(!esdTrackCuts->AcceptTrack(track))  { continue; } 

      isTPC = kTRUE;
      isMatch = kFALSE;
      if( (track->GetStatus()&AliESDtrack::kITSrefit) && 
	  (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) ) 
      {
        isMatch = kTRUE;
      }

      //
      FillHistograms(track, stack, isMatch, isTPC, kFALSE);
      //if(tpcTrack) delete tpcTrack;
    } 

    //
    // ITS -> TPC prolongation efficiency
    //
    for(Int_t iTrack=0; iTrack<entries;++iTrack) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(iTrack);
      if(!track) continue;

      //
      // ITS stand alone tracks
      //
      if(!(track->GetStatus() & AliESDtrack::kITSpureSA)) continue;
      if(!(track->GetStatus() & AliESDtrack::kITSrefit)) continue;
      if(track->GetNcls(0)<4) continue;
      if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) continue;

      // Check matching with TPC only track
      for(Int_t jTrack=0; jTrack<entries;++jTrack) 
      {
        isMatch = kFALSE;

        if(iTrack==jTrack) continue;

        AliESDtrack *track2 = (AliESDtrack*)allChargedTracks->At(jTrack);
        if(!track2) continue;
        if(track2->Charge()==0) continue;
	if(!track2->GetTPCInnerParam()) continue;
        if(!(track2->GetStatus() & AliESDtrack::kTPCrefit)) continue;

        // Get TPC only tracks (must be deleted by user) 
        AliESDtrack* tpcTrack2 = AliESDtrackCuts::GetTPCOnlyTrack(esdEvent, jTrack);
        if(!tpcTrack2) continue;
        if(!tpcTrack2->RelateToVertex(vtxESD,esdEvent->GetMagneticField(),100.)) { delete tpcTrack2; continue; } 

        // check loose cuts for TPC tracks
        if(!esdTrackCuts->AcceptTrack(tpcTrack2)) { delete tpcTrack2; continue; }

        // check matching
        if (TMath::Abs(track->GetY() - tpcTrack2->GetY()) > 3) { delete tpcTrack2; continue; }
        if (TMath::Abs(track->GetSnp() - tpcTrack2->GetSnp()) > 0.2) { delete tpcTrack2; continue; }
        if (TMath::Abs(track->GetTgl() - tpcTrack2->GetTgl()) > 0.2) { delete tpcTrack2; continue; }

	isMatch = kTRUE;
	if(tpcTrack2) { 
	  delete tpcTrack2;
	}
	  break;
      }

       //
       FillHistograms(track, stack, isMatch, kFALSE, kTRUE);
    } 
  }

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;

}

//_____________________________________________________________________________
void AlidNdPtEfficiency::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Bool_t isMatch, const Bool_t isTPC,  const Bool_t isITSTPC) const
{
  //
  // Fill ESD track and MC histograms 
  //
  if(!esdTrack) return;
  Int_t charge = esdTrack->Charge();
  if(charge == 0.) return;

  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();

  //
  // Fill rec vs MC information
  //
  Bool_t isPrim = kTRUE;

  if(IsUseMCInfo()) {
    if(!stack) return;
    Int_t label = esdTrack->GetLabel(); 
    if(label < 0.) return; // fake ITS track
    TParticle* particle = stack->Particle(label);
    if(!particle) return;
    if(particle->GetPDG() && particle->GetPDG()->Charge()==0.) return;
    isPrim = stack->IsPhysicalPrimary(label);
  }

  // fill histo
  Double_t vRecMCTrackHist[6] = { eta,phi,pt,isPrim,charge,isMatch }; 
  Double_t vRecMCTrackHistTPCITS[7] = { eta,phi,pt,isPrim,charge,isMatch,isTPC }; 

  if(isITSTPC) {
    fRecMCTrackHistITSTPC->Fill(vRecMCTrackHist);
  }
  else {
    fRecMCTrackHistTPCITS->Fill(vRecMCTrackHistTPCITS);
  }
}

//_____________________________________________________________________________
Long64_t AlidNdPtEfficiency::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtEfficiency* entry = dynamic_cast<AlidNdPtEfficiency*>(obj);
    if (entry == 0) continue; 
  
    // track histo
    fRecMCTrackHistTPCITS->Add(entry->fRecMCTrackHistTPCITS);
    fRecMCTrackHistITSTPC->Add(entry->fRecMCTrackHistITSTPC);

  count++;
  }

return count;
}

//_____________________________________________________________________________
void AlidNdPtEfficiency::Analyse() 
{
  //
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;
  if(!aFolderObj) return; 

  TH1D *h1Dall = 0; 
  TH1D *h1D = 0; 
  TH1D *h1Dc = 0; 


  //
  // get cuts
  //
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts || !esdTrackCuts) {
    Error("AlidNdPtEfficiency::Analyse()", "cuts not available");
    return;
  }

  //
  // TPC->ITS efficiency
  //

  //eff vs eta
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(2,2);  
  h1Dall = (TH1D *)fRecMCTrackHistTPCITS->Projection(0);
  if(!h1Dall) return;
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistTPCITS->Projection(0);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_eta_TPCITS");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(1,2);  
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(1,2);  

  //eff vs phi
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(2,2);  
  fRecMCTrackHistTPCITS->GetAxis(0)->SetRangeUser(-0.8, 0.799);  
  h1Dall = (TH1D *)fRecMCTrackHistTPCITS->Projection(1);
  if(!h1Dall) return;
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistTPCITS->Projection(1);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_phi_TPCITS");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(1,2);  
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(1,2);  

  //eff vs pT
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(2,2);  
  fRecMCTrackHistTPCITS->GetAxis(0)->SetRangeUser(-0.8, 0.799);  
  h1Dall = (TH1D *)fRecMCTrackHistTPCITS->Projection(2);
  if(!h1Dall) return;
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistTPCITS->Projection(2);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_pT_TPCITS");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistTPCITS->GetAxis(5)->SetRange(1,2);  
  fRecMCTrackHistTPCITS->GetAxis(6)->SetRange(1,2);  


  //
  // ITS->TPC efficiency
  //

  fRecMCTrackHistITSTPC->GetAxis(0)->SetRangeUser(-1.5, 1.499);  
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(1,2);  

  //eff vs eta
  h1Dall = (TH1D *)fRecMCTrackHistITSTPC->Projection(0);
  if(!h1Dall) return;
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistITSTPC->Projection(0);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_eta_ITSTPC");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(1,2);  

  //eff vs phi
  fRecMCTrackHistITSTPC->GetAxis(0)->SetRangeUser(-0.8, 0.799);  
  h1Dall = (TH1D *)fRecMCTrackHistITSTPC->Projection(1);
  if(!h1Dall) return;
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistITSTPC->Projection(1);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_phi_ITSTPC");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(1,2);  

  //eff vs pT
  fRecMCTrackHistITSTPC->GetAxis(0)->SetRangeUser(-0.8, 0.799);  
  h1Dall = (TH1D *)fRecMCTrackHistITSTPC->Projection(2);
  if(!h1Dall) return;
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(2,2);  
  h1D = (TH1D *)fRecMCTrackHistITSTPC->Projection(2);
  if(!h1D) return;

  h1Dc = (TH1D *)h1D->Clone("eff_vs_pT_ITSTPC");
  h1Dc->Divide(h1Dall);
  aFolderObj->Add(h1Dc);
  fRecMCTrackHistITSTPC->GetAxis(5)->SetRange(1,2);  
  
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);
  if(!fAnalysisFolder) {
    if(aFolderObj) delete aFolderObj;
    return;
  }

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AlidNdPtEfficiency::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtEfficiency * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}

//_____________________________________________________________________________
TFolder* AlidNdPtEfficiency::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
