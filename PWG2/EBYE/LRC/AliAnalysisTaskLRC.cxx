/**************************************************************************
 * Author: Andrey Ivanov.                                           *
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
// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCBase objects that are processing LRC analysis
// for a given Eta window 
// This task is worcking with ESD data only

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliAnalysisTaskLRC.h"
#include <AliLRCBase.h>
#include <AliVEvent.h>
#include <AliMCEvent.h>
#include <AliESDEvent.h>     
#include <AliPID.h>

#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliStack.h>
#include <AliESDtrackCuts.h>
#include <AliPhysicsSelection.h>
#include <TParticle.h>
#include <TH1.h>
#include <TH2.h>


ClassImp(AliAnalysisTaskLRC)

//________________________________________________________________________
AliAnalysisTaskLRC::AliAnalysisTaskLRC(const char *name,Bool_t runKine) 
  : AliAnalysisTaskSE(name),fEsdTrackCuts(0),fMaxPtLimit(5.0),fMinPtLimit(0.0),fMinAceptedTracksCut(0),fMaxAceptedTracksCut(0),fCheckForkVtx(kTRUE),fCheckForVtxPosition(kFALSE),fVxMax(0),fVyMax(0),fVzMax(0),fLRCproc(0),fOutList(0),fRunKine(0),fShowEventStats(kFALSE),fShowPerTrackStats(kFALSE),fHistEventCutStats(0),fHistTrackCutStats(0),fHistVx(0),fHistVy(0),fHistVz(0),fHistEtaVsZvCoverage(0),fHistEtaVsZvCoverageAcepted(0),fHistAceptedMult(0)
{
  //Init
  fRunKine=runKine;

  // Output slot #1 writes into a TList container for common task data and QA
  DefineOutput(1, TList::Class());

  //Defining output slots for each LRC processor (required to avoid TList of TLists on merging)
  for(Int_t i=0; i < fLRCproc.GetEntries(); i++)
  {
  DefineOutput(Proc(i)->GetOutputSlotNumber(),TList::Class());
  }
}


// ---------------------------------------  Setters ------------------

  void AliAnalysisTaskLRC::SetMaxPtLimit(Double_t MaxPtLimit)
  {
   //Sets  Max Pt filter
     fMaxPtLimit=MaxPtLimit;
  }
  void AliAnalysisTaskLRC::SetMinPtLimit(Double_t MinPtLimit)
  {
     //Sets  Min Pt filter 
  fMinPtLimit=MinPtLimit;
  }
  AliLRCBase*  AliAnalysisTaskLRC::Proc(Int_t index) 
  {
  // Get Processor i
  return (dynamic_cast<AliLRCBase*> (fLRCproc.At(index)));
  }

//________________________________________________________________________
void AliAnalysisTaskLRC::UserCreateOutputObjects()
{
  //Disabling "replacing existing TH2D ...etc" warning
   
   Bool_t lTH1oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
  
   Bool_t lTH2oldStatus = TH2::AddDirectoryStatus();
   TH2::AddDirectory(kFALSE);
  
  
  //LRC processors init
  Int_t lLrcNum=fLRCproc.GetEntries();
  
  for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
    if(lrcBase)
      lrcBase->InitDataMembers();
    else continue;
  }
  

    // --------- Output list
  
  fOutList = new TList();
 
  fOutList->Add( fEsdTrackCuts );

  //Event level QA
  fHistEventCutStats = new TH1D("fHistEventCutStats","Event statistics;;N_{events}", 7,-0.5,6.5);
  TString gEventCutBinNames[7] = {"Total","No trigger","No vertex","Bad vertex position","HighMult cut","LowMult cut","Analyzed"};
  for(Int_t i = 1; i <= 7; i++)fHistEventCutStats->GetXaxis()->SetBinLabel(i,gEventCutBinNames[i-1].Data());
  fOutList->Add(fHistEventCutStats);

  //Track level QA
  fHistTrackCutStats = new TH1D("fHistTrackCutStats","Event statistics;;N_{events}", 5,-0.5,4.5);
  TString gTrackCutBinNames[5] = {"Total","AliESDtrackCuts","HighPt cut","LowPt cut","Good"};
  for(Int_t i = 1; i <= 5; i++)fHistTrackCutStats->GetXaxis()->SetBinLabel(i,gTrackCutBinNames[i-1].Data());
  fOutList->Add(fHistTrackCutStats);


  //Vertex distributions
  fHistVx = new TH1D("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fOutList->Add(fHistVx);
  fHistVy = new TH1D("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fOutList->Add(fHistVy);
  fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fOutList->Add(fHistVz);

  //Eta vs Zv
  fHistEtaVsZvCoverage = new TH2D("fHistEtaVsZvCoverage","TPC tracks ETA vs Zv;V_{z} (cm);ETA",100,-20,20,50,-2,2);
  fHistEtaVsZvCoverageAcepted = new TH2D("fHistEtaVsZvCoverage","Acepted TPC tracks ETA vs Zv;V_{z} (cm);ETA",100,-20,20,50,-2,2);
  fOutList->Add(fHistEtaVsZvCoverage);
  fOutList->Add(fHistEtaVsZvCoverageAcepted);

  fHistAceptedMult = new TH1D("fHistAceptedMult","N_{ch} - accepted tracks;N_{ch} acepted;Entries",401,-0.5,400.5);
  fOutList->Add(fHistAceptedMult);
  
  // Returning TH1,TH2 AddDirectory status 
  TH1::AddDirectory(lTH1oldStatus);
  TH2::AddDirectory(lTH2oldStatus);
    
}

//________________________________________________________________________
void AliAnalysisTaskLRC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliVEvent *event = InputEvent();
   AliESDEvent *fESD = 0x0;
   fESD=dynamic_cast<AliESDEvent*> (event) ;
   AliStack *stack = 0x0;
   AliMCEvent *eventMC = 0x0;
   if( fRunKine ) 
   {
	eventMC=MCEvent(); 
	stack=eventMC->Stack();
	Printf("Number of primaries: %d",stack->GetNprimary());
   }


  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }
  if((fESD)&&(!fEsdTrackCuts)){ AliDebug(AliLog::kError, "No ESD track cuts avalible"); return;}

   // Numbrt of processors attached
   Int_t lLrcNum=fLRCproc.GetEntries();
   

// Processing event selection 
fHistEventCutStats->Fill(0);

Bool_t lTrigger=kTRUE;
Bool_t lVertexPresent=kTRUE;
Bool_t lVertexAceptable=kTRUE;


//Trigger
   if(fESD){
	if(fShowEventStats)Printf("Trigger classes: %s:", fESD->GetFiredTriggerClasses().Data());
	
   	lTrigger=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
   	if(!lTrigger){
   		if(fShowEventStats)Printf("Rejected!");
		fHistEventCutStats->Fill(1);
 		PostData(1, fOutList);
		return;
	}
   }

// Vertex present
const AliESDVertex *vertex = fESD->GetPrimaryVertex();
 if(vertex) {
   lVertexPresent=((vertex)&&(vertex->GetNContributors() > 0)&&(vertex->GetZRes() != 0));
   if((!lVertexPresent)&&fCheckForkVtx)
     {if(fShowEventStats)Printf("No vertex");
       fHistEventCutStats->Fill(2);
       PostData(1, fOutList);
       return;
     }

   // Vertex in range 
   lVertexAceptable=(TMath::Abs(vertex->GetXv()) < fVxMax) && (TMath::Abs(vertex->GetYv()) < fVyMax); 
   if(lVertexAceptable) {
     if(fVzMax>0)   //   fVzMax < 0 -> select Zv outside selected range
       { lVertexAceptable = (TMath::Abs(vertex->GetZv()) < fVzMax);}
     else
       { lVertexAceptable = (TMath::Abs(vertex->GetZv()) > -fVzMax);}
   }
   if((!lVertexAceptable) && fCheckForVtxPosition) 
     {if(fShowEventStats)Printf("Vertex out of range");
       fHistEventCutStats->Fill(3);
       PostData(1, fOutList);
       return;
     }
   
   fHistVx->Fill(vertex->GetXv());
   fHistVy->Fill(vertex->GetYv());
   fHistVz->Fill(vertex->GetZv());
 

int lNchTrigger=0;
// Pre event loop 
if(fESD &&  !fRunKine )
	for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++)
		if(fEsdTrackCuts->AcceptTrack(fESD->GetTrack(iTracks)))lNchTrigger++;

//if(fRunKine)lNchTrigger=eventMC->GetNumberOfPrimaries();

// Nch bins cut

if( (lNchTrigger > fMaxAceptedTracksCut) && (fMaxAceptedTracksCut != 0))
	{
	fHistEventCutStats->Fill(4);
 	PostData(1, fOutList);
	return;
	}


if(lNchTrigger < fMinAceptedTracksCut)
	{	
	fHistEventCutStats->Fill(5);
	 PostData(1, fOutList);
	return;
	}
fHistAceptedMult->Fill(lNchTrigger);
fHistEventCutStats->Fill(6);


  //Track variables
  double lPt;   // Temp Pt
  double lEta;	  // Temp ETA	
  double lPhi;    // Temp Phi
  
  //Track selection counters 
  int lNaccept=0;
  int lPtOver=0;
  int lPtUnder=0;
  int lNoCharge=0;
  int lCutfail=0;
  int lDecay=0;

 for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
    if(lrcBase)
      lrcBase->StartEvent();
    else continue;
  }

    // Track loop -----------------------------------------------------------------
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
  
    AliVParticle* track = event->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    
    fHistTrackCutStats->Fill(0);
    fHistEtaVsZvCoverage->Fill(vertex->GetZv(),track->Eta());
	
    lPt=track->Pt();
    lEta=track->Eta();
    lPhi=track->Phi();
    
    // ESD track cuts
    if(fESD &&  !fRunKine )
    {
    	AliESDtrack* lESDtarck= fESD->GetTrack(iTracks);
	if(fShowPerTrackStats)printf("ESD Track N%d , Eta=%f, Pt=%f , TPC Nclusters =%d Sigma=%f ",  iTracks , lEta , lPt, lESDtarck-> GetTPCNcls(),fEsdTrackCuts->GetSigmaToVertex( lESDtarck) );
	if(!fEsdTrackCuts->AcceptTrack(lESDtarck))
		{	
		lCutfail++;
		if(fShowPerTrackStats)Printf("Rejected by cuts");
		fHistTrackCutStats->Fill(1);
		continue;
		}
	else
		{
		if(fShowPerTrackStats)Printf("OK");
		}	
    
    }
    // end of ESD tack cuts

      if(lPt>fMaxPtLimit)
	{lPtOver++; 
	fHistTrackCutStats->Fill(2);
	continue;}     // Dropping traks with hi Pt
      
	if(lPt<fMinPtLimit)
	{lPtUnder++;
	fHistTrackCutStats->Fill(3);
	continue;}   // Dropping traks with lo Pt

   
	fHistTrackCutStats->Fill(4);
 
      for(Int_t i=0; i < lLrcNum; i++)
  {
    AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
    if(lrcBase)
      lrcBase->AddTrackPtEta(lPt,lEta,lPhi);
    else continue;
  }

    lNaccept++;
    fHistEtaVsZvCoverageAcepted->Fill(vertex->GetZv(),track->Eta());
        
  } //end of track loop 
  
  
  
   for(Int_t i=0; i < lLrcNum; i++)
  {
      AliLRCBase *lrcBase = dynamic_cast<AliLRCBase*> (fLRCproc.At(i));
    if(lrcBase)
      lrcBase->FinishEvent();
    else continue;
  }

  //Debuging output of track filter
  if(fShowPerTrackStats)Printf("NofTracks= %d , accepted %d , LowPt %d, HighPt %d, LowCharge %d,  ESD cut fail %d , Decay Filer %d", event->GetNumberOfTracks(), lNaccept, lPtUnder, lPtOver ,lNoCharge , lCutfail  , lDecay);
  
   // Post output data.
  
  PostData(1, fOutList);
  
  for(Int_t i=0; i < fLRCproc.GetEntries(); i++)
  {
  PostData(Proc(i)->GetOutputSlotNumber(),Proc(i)->CreateOutput());
  }
 
 }//vertex
  
}      

//________________________________________________________________________
void AliAnalysisTaskLRC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
//fOutList = dynamic_cast<TList*> (GetOutputData(0));

//fEsdTrackCuts->DrawHistograms();

 
}


void AliAnalysisTaskLRC::AddLRCProcess(AliLRCBase *newProc)
{
// Add new AliLRCBase (Main LRC processor per ETA window) to the processing list
// Used to add new ETA window to AnalysisTask 
if(!newProc)
{Printf("ERROR:No AliLRCBase object -  NULL pointer!");
return;
}


fLRCproc.Add(newProc);
newProc->SetOutputSlotNumber(fLRCproc.GetEntries()+1);
DefineOutput(newProc->GetOutputSlotNumber(),TList::Class());
return ;
}

