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
/* $Id: $ */

//_________________________________________________________________________
// Count events with different selections
//
// It produces a histogram with the number of events with 9 bins:
// 0: all events (that passed the physics selection if it was on)
// 1: same but cross check that event pointer did exist
// 2: passes vertex cut
// 3: passes track number cut, tracks for eta < 0.8
// 4: 3 && 2
// 5: pass VAND
// 6: 5 && 2
// 7: 5 && 3
// 8: 5 && 3 && 2
// 9: not pileup from SPD
// 10: Good vertex
// 11: 10 && 5
// 12: 10 && 3
// 13: 10 && 2
// 14: 10 && 2 && 3 && 5
// 15: 10 && 9
// 16: 9  && 5
//
// Author: Gustavo Conesa Balbastre (LPSC)
//         
//_________________________________________________________________________

#include "TH2F.h"
#include "AliAODHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskCounter.h"
ClassImp(AliAnalysisTaskCounter)

//________________________________________________________________________
AliAnalysisTaskCounter::AliAnalysisTaskCounter(const char *name) 
: AliAnalysisTaskSE(name), 
  fZVertexCut(10.),
  fTrackMultEtaCut(0.8),
  fCaloFilterPatch(kFALSE),
  fOutputContainer(0x0), 
  fESDtrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
  fTriggerAnalysis (new AliTriggerAnalysis),
  fhNEvents(0),
  fhXVertex(0),fhYVertex(0),fhZVertex(0),
  fhXGoodVertex(0),fhYGoodVertex(0),fhZGoodVertex(0)
{
  //ctor
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskCounter::AliAnalysisTaskCounter() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskCounter"),
    fZVertexCut(10.),
    fTrackMultEtaCut(0.8),
    fCaloFilterPatch(kFALSE),
    fOutputContainer(0x0), 
    fESDtrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
    fTriggerAnalysis (new AliTriggerAnalysis),
    fhNEvents(0),
    fhXVertex(0),fhYVertex(0),fhZVertex(0),
    fhXGoodVertex(0),fhYGoodVertex(0),fhZGoodVertex(0)
{
  // ctor
  DefineOutput(1, TList::Class());
}

//__________________________________________________
AliAnalysisTaskCounter::~AliAnalysisTaskCounter()
{
  //Destructor
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
  
  if(fESDtrackCuts)    delete fESDtrackCuts;
  if(fTriggerAnalysis) delete fTriggerAnalysis;
  
}


//-------------------------------------------------------------------
void AliAnalysisTaskCounter::UserCreateOutputObjects()
{
  // Init histograms
  
  fOutputContainer = new TList();
  
  fhZVertex     = new TH1F("hZVertex", " Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZVertex);

  fhZGoodVertex     = new TH1F("hZGoodVertex", " Good Z vertex distribution"   , 200 , -50 , 50  ) ;
  fhZGoodVertex->SetXTitle("v_{z} (cm)");
  fOutputContainer->Add(fhZGoodVertex);
  
  fhXVertex     = new TH1F("hXVertex", " X vertex distribution"   , 200 , -2 , 2  ) ;
  fhXVertex->SetXTitle("v_{x} (cm)");
  fOutputContainer->Add(fhXVertex);
  
  fhXGoodVertex     = new TH1F("hXGoodVertex", " Good X vertex distribution"   , 200 , -2 , 2  ) ;
  fhXGoodVertex->SetXTitle("v_{x} (cm)");
  fOutputContainer->Add(fhXGoodVertex);
  
  fhYVertex     = new TH1F("hYVertex", " Y vertex distribution"   , 200 , -2 , 2  ) ;
  fhYVertex->SetXTitle("v_{y} (cm)");
  fOutputContainer->Add(fhYVertex);
  
  fhYGoodVertex     = new TH1F("hYGoodVertex", " Good Y vertex distribution"   , 200 , -2 , 2  ) ;
  fhYGoodVertex->SetXTitle("v_{y} (cm)");
  fOutputContainer->Add(fhYGoodVertex);
  
  
  fhNEvents = new TH1I("hNEvents", "Number of analyzed events", 20, 0, 20) ;
  fhNEvents->SetXTitle("Selection");
  fhNEvents->SetYTitle("# events");
  fhNEvents->GetXaxis()->SetBinLabel(1 ,"1  = PS");
  fhNEvents->GetXaxis()->SetBinLabel(2 ,"2  = 1  & ESD");
  fhNEvents->GetXaxis()->SetBinLabel(3 ,"3  = 2  & |Z|<10");
  fhNEvents->GetXaxis()->SetBinLabel(4 ,"4  = 2  & |eta|<0.8");
  fhNEvents->GetXaxis()->SetBinLabel(5 ,"5  = 3  & 4");
  fhNEvents->GetXaxis()->SetBinLabel(6 ,"6  = 2  & V0AND");
  fhNEvents->GetXaxis()->SetBinLabel(7 ,"7  = 3  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(8 ,"8  = 4  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(9 ,"9  = 5  & 6");
  fhNEvents->GetXaxis()->SetBinLabel(10,"10 = 2  & not pileup");
  fhNEvents->GetXaxis()->SetBinLabel(11,"11 = 2  & good vertex");
  fhNEvents->GetXaxis()->SetBinLabel(12,"12 = 3  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(13,"13 = 4  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(14,"14 = 6  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(15,"15 = 9  & 11");
  fhNEvents->GetXaxis()->SetBinLabel(16,"16 = 10 & 11");
  fhNEvents->GetXaxis()->SetBinLabel(17,"17 = 6  & 10");
  fhNEvents->GetXaxis()->SetBinLabel(17,"17 = 1  & |Z|<50");  
  fhNEvents->GetXaxis()->SetBinLabel(18,"18 = Reject EMCAL");
  fhNEvents->GetXaxis()->SetBinLabel(19,"19 = 18 & 2");
  fhNEvents->GetXaxis()->SetBinLabel(20,"20 = 18 & |Z|<50");

  fOutputContainer->Add(fhNEvents);

  fOutputContainer->SetOwner(kTRUE);
  
  PostData(1,fOutputContainer);
  
}

//________________________________________________________________________
void AliAnalysisTaskCounter::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  //printf("___ Event __ %d __\n",(Int_t)Entry());
  
  fhNEvents->Fill(0.5);  
  
  AliVEvent * event = InputEvent();
  if (!event) {
    printf("AliAnalysisTaskCounter::UserExec() - ERROR: event not available \n");
    return;
  }
  AliESDEvent * esdevent = dynamic_cast<AliESDEvent*> (event);

  fhNEvents->Fill(1.5);  
    
  //Initialize bools
  Bool_t bSelectVZ    = kFALSE;
  Bool_t bV0AND       = kFALSE; 
  Bool_t bPileup      = kFALSE;
  Bool_t bGoodV       = kFALSE;
  Bool_t bSelectTrack = kFALSE;   
  Int_t  trackMult    = 0;
  
  //---------------------------------
  //Get the primary vertex, cut on Z
  //---------------------------------
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v) ;
  fhXVertex->Fill(v[0]);
  fhYVertex->Fill(v[1]);
  fhZVertex->Fill(v[2]);
  
  if(TMath::Abs(v[2]) < fZVertexCut) {
    bSelectVZ=kTRUE;
    fhNEvents->Fill(2.5);  
  }
  //else printf("Vertex out %f \n",v[2]);
  

  //--------------------------------------------------
  //Tweak for calorimeter only productions
  //--------------------------------------------------
  if(fCaloFilterPatch && !esdevent){ 
    if(event->GetNumberOfCaloClusters() > 0) {
      AliVCluster * calo = event->GetCaloCluster(0);
      if(calo->GetNLabels() == 4){
        Int_t * selection = calo->GetLabels();
        bPileup   = selection[0];
        bGoodV    = selection[1]; 
        bV0AND    = selection[2]; 
        trackMult = selection[3];
        //if(selection[0] || selection[1] || selection[2])
        //printf(" pu %d, gv %d, v0 %d, track mult %d\n ", selection[0], selection[1], selection[2], selection[3]);
        if(trackMult > 0 )  
          bSelectTrack = kFALSE;
      } else {
        //First filtered AODs, track multiplicity stored there.  
        trackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
      }
    }else{//at least one cluster
        //printf("AliAnalysisTaskCounter::UserExec() - No clusters in event\n");
        //Remove events with  vertex (0,0,0), bad vertex reconstruction
        if(TMath::Abs(v[0]) < 1.e-6 && TMath::Abs(v[1]) < 1.e-6 && TMath::Abs(v[2]) < 1.e-6) bGoodV = kFALSE;
      
        //First filtered AODs, track multiplicity stored there.  
        trackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
    }
  }
  else {
    //--------------------------------------------------
    //Count tracks, cut on number of tracks in eta < 0.8
    //--------------------------------------------------
    Int_t nTracks   = event->GetNumberOfTracks() ;
    for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
      AliVTrack * track = (AliVTrack*)event->GetTrack(itrack) ; // retrieve track from esd
      
      //Only for ESDs
      if(esdevent && !fESDtrackCuts->AcceptTrack((AliESDtrack*)track)) continue;
      
      //Do not count tracks out of acceptance cut
      if(TMath::Abs(track->Eta())< fTrackMultEtaCut) trackMult++;
    }
  }
  
  //printf("AliAnalysisTaskCounter::UserExec() - Track Mult %d \n",trackMult);
  
  //--------------------------------------------------
  // At least one track
  //--------------------------------------------------
  if (trackMult > 0) {
    bSelectTrack = kTRUE; 
                  fhNEvents->Fill(3.5);
    if(bSelectVZ) fhNEvents->Fill(4.5);
  }
  
  //---------------------------------
  // V0AND
  //---------------------------------
  if(esdevent && !fCaloFilterPatch) bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esdevent, AliTriggerAnalysis::kV0AND);
  //else if(aodevent  && !fCaloFilterPatch) bV0AND = //FIXME FOR AODs
  
  if(bV0AND)
  {
                                    fhNEvents->Fill(5.5);
    if (bSelectVZ)                  fhNEvents->Fill(6.5);
    if (bSelectTrack)               fhNEvents->Fill(7.5);
    if (bSelectVZ &&  bSelectTrack) fhNEvents->Fill(8.5);
  }
  
  //---------------------------------
  // Pileup
  //---------------------------------
  if(!fCaloFilterPatch)
    bPileup = event->IsPileupFromSPD(3, 0.8, 3., 2., 5.); //Default values, if not it does not compile
  //bPileup = event->IsPileupFromSPD(); 
  if (!bPileup){
                fhNEvents->Fill(9.5);
    if(bV0AND)  fhNEvents->Fill(16.5);
  }
  
  //---------------------------------
  // Good vertex
  //---------------------------------
  if(esdevent && !fCaloFilterPatch) bGoodV = CheckForPrimaryVertex();
  if(bGoodV) 
  {
    fhXGoodVertex->Fill(v[0]);
    fhYGoodVertex->Fill(v[1]);
    fhZGoodVertex->Fill(v[2]);
    
                     fhNEvents->Fill(10.5);
    if(bSelectVZ)    fhNEvents->Fill(11.5);
    if(bSelectTrack) fhNEvents->Fill(12.5);
    if(bV0AND)       fhNEvents->Fill(13.5);
    if(bSelectVZ && bSelectTrack && bV0AND)    
                     fhNEvents->Fill(14.5); 
    if(!bPileup)     fhNEvents->Fill(15.5); 
  }

  //printf("AliAnalysisTaskCounter::UserExec() : z vertex %d, good vertex %d, v0and %d, pile up %d, track mult %d\n ", bSelectVZ, bGoodV, bV0AND, bPileup, trackMult);
  
  // Events that could be rejected in EMCAL
  for (Int_t i = 0; i < InputEvent()->GetNumberOfCaloClusters(); i++)
  {
    AliVCluster *clus = InputEvent()->GetCaloCluster(i);
    if(clus->IsEMCAL()){    
      
      if ((clus->E() > 500 && clus->GetNCells() > 200 ) || clus->GetNCells() > 300)  {
        //printf("Counter: Reject event with cluster: E %f, ncells %d\n",clus->E(),clus->GetNCells());
                         fhNEvents->Fill(17.5); 
        if(bSelectVZ)    fhNEvents->Fill(18.5);
        if(TMath::Abs(v[2]) < 50.)  fhNEvents->Fill(19.5); 
        break;
      }
        
    }
  }
  
  PostData(1,fOutputContainer);

}

//____________________________________________________________________________
Bool_t AliAnalysisTaskCounter::CheckForPrimaryVertex(){
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //It only works for ESDs
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (InputEvent());
  if(!event) return 0;
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0) {
    return 1;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1) {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0) {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return 1;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1) {
      //      cout<<"bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return 0;
    }
  }
  return 0;
  //return fInputEvent->GetPrimaryVertex()->GetNContributors()>0;
}



//_____________________________________________________
void AliAnalysisTaskCounter::FinishTaskOutput()
{
  // Put in the output some event summary histograms
  
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputH = dynamic_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!inputH) return; 
  TH2F *histStat = dynamic_cast<TH2F*>(inputH->GetStatistics()); 
  TH2F *histBin0 = dynamic_cast<TH2F*>(inputH->GetStatistics("BIN0"));
  
  if(histStat)
    fOutputContainer->Add(histStat);
  else
    printf("AliAnalysisTaskCounter::FinishTaskOutput() - Stat histogram not available check, \n if ESDs, that AliPhysicsSelection was on, \n if AODs, if EventStat_temp.root exists \n");

  if(histBin0)
    fOutputContainer->Add(histBin0); 
  
}
