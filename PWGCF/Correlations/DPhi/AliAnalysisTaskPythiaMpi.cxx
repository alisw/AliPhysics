// **************************************
// Task used for Pythia MPI studies in ZZ train 
// Output is stored in an exchance container
// *******************************************


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

#include "TChain.h"
#include "TH3.h"
#include "TProfile.h"

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliHeader.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskPythiaMpi.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVParticle.h"
#include "AliPythia.h"
#include "AliPythia8.h"
//#include "AliAnalysisHelperJetTasks.h"
#include <typeinfo>

using namespace std;

ClassImp(AliAnalysisTaskPythiaMpi)

//_________________________________________________________| Constructor

AliAnalysisTaskPythiaMpi::AliAnalysisTaskPythiaMpi(): 
  fMcEvent(0x0),  
  fMcHandler(0x0),
  fOutputList(0)
{
  //
  // Default constructor
  //
}

////_________________________________________________________| Specific Constructor

AliAnalysisTaskPythiaMpi::AliAnalysisTaskPythiaMpi(const char* name):
  AliAnalysisTaskSE(name),
  fMcEvent(0x0),  
  fMcHandler(0x0),
  fOutputList(0)
{
  //
  // Specific constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisMpi","Calling Constructor");
  // Output slot #1 writes into a TList container 
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());  
}

//_________________________________________________________| Assignment Operator

AliAnalysisTaskPythiaMpi& AliAnalysisTaskPythiaMpi::operator=(const AliAnalysisTaskPythiaMpi& c)
{
  // Assigment operator
  if(this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fMcEvent = c.fMcEvent;
    fMcHandler = c.fMcHandler;
    fOutputList = c.fOutputList;

  }
  return *this;

}
//_________________________________________________________| Copy Constructor

AliAnalysisTaskPythiaMpi::AliAnalysisTaskPythiaMpi(const AliAnalysisTaskPythiaMpi& c):
  AliAnalysisTaskSE(c),
  fMcEvent(c.fMcEvent),
  fMcHandler(c.fMcHandler), 
  fOutputList(c.fOutputList)
{
  //Copy Constructor
}

//_________________________________________________________| Destructor

AliAnalysisTaskPythiaMpi::~AliAnalysisTaskPythiaMpi()
{
  //
  // Destructor
  //
  Info("~AliAnalysisTaskPythiaMpi","Calling Destructor");
  if(fOutputList) delete fOutputList;
  
}

//_____________________________________________________________________

void AliAnalysisTaskPythiaMpi::UserCreateOutputObjects()
{
  Info("CreateOutputObjects","CreateOutputObjects of task %s",GetName());
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
 
  TH1I *fHistEvents = new TH1I("fHistNEvents","fHistNEvents",2,0,2);
  fHistEvents->GetXaxis()->SetBinLabel(1,"All events");

  TH1F *fHistPt = new TH1F("fHistPt","p_{T} distribution",39,0.15,10.0);
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c^{-1})");
  fHistPt->SetMarkerStyle(kFullCircle);

  TH1F *fHistEta = new TH1F("fHistEta","eta distribution",20,-1,1);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("dN/d#eta");
  fHistEta->SetMarkerStyle(kFullCircle);

  TH1F *fHistMpi = new TH1F("fHistMpi","MPIs distribution",100,-0.5,99.5);
  fHistMpi->GetXaxis()->SetTitle("#MPIs");
  fHistMpi->GetYaxis()->SetTitle("dN/dMPIs");
  fHistMpi->SetMarkerStyle(kFullCircle);

  TH2F *fHistMultMpi = new TH2F("fHistMultMpi","Multiplicity vs MPIs",100,-0.5,99.5,100,0,1e4);
  fHistMultMpi->GetXaxis()->SetTitle("#MPIs");
  fHistMultMpi->GetYaxis()->SetTitle("MC particles");
  fHistMultMpi->SetMarkerStyle(kFullDotSmall);

  TH2F *fHistdNdetaMpi = new TH2F("fHistdNdetaMpi","dN/d#etadMPIs",100,-0.5,99.5,100,-5,5);
  fHistdNdetaMpi->GetXaxis()->SetTitle("#MPIs");
  fHistdNdetaMpi->GetYaxis()->SetTitle("#eta");

  TH2F *fHistPtMpi = new TH2F("fHistPtMpi","fHistPtMpi",100,-0.5,99.5,39,0.15,10.0);
  fHistPtMpi->GetXaxis()->SetTitle("#MPIs");
  fHistPtMpi->GetYaxis()->SetTitle("p_{T} (GeV/c)");

  TProfile* fProfileMpiPt = new TProfile("fProfileMpiPt","fProfileMpiPt",100,-0.5,99.5,0.15,10.);
  fProfileMpiPt->GetXaxis()->SetTitle("#MPIs");
  fProfileMpiPt->GetYaxis()->SetTitle("p_{T} (GeV/c)");

  TProfile* fProfilePtMpi = new TProfile("fProfilePtMpi","fProfilePtMpi",39,0.15,10.,-0.5,99.5);
  fProfilePtMpi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fProfilePtMpi->GetYaxis()->SetTitle("#MPIs");

  TH3F *fHistTracks = new TH3F("fHistTracks","fHistTracks",20,-1.,1.,39,0.15,10.0,100,-0.5,99.5);
  fHistTracks->GetXaxis()->SetTitle("#eta");
  fHistTracks->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fHistTracks->GetZaxis()->SetTitle("#MPIs");


  fOutputList->Add(fHistEvents);
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistMpi);
  fOutputList->Add(fHistMultMpi);
  fOutputList->Add(fHistdNdetaMpi);
  fOutputList->Add(fHistPtMpi);
  fOutputList->Add(fProfileMpiPt);
  fOutputList->Add(fProfilePtMpi);
  fOutputList->Add(fHistTracks);

  PostData(1,fOutputList);

  return;
}
//_________________________________________________________________________

void AliAnalysisTaskPythiaMpi::Init()
{
   // MC handler
   if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Init() \n");
   fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); 
}

//_________________________________________________________________________

void AliAnalysisTaskPythiaMpi::LocalInit()
{
   // MC handler
   //
   // Initialization
   //

   if(fDebug > 1) printf("AnalysisTaskPythiaMpi::LocalInit() \n");
   Init();
}

//____________________________________________| User Exec

void AliAnalysisTaskPythiaMpi::UserExec(Option_t*)
{

   // handle and reset the output jet branch 
  Info("UserExec","Start of method");
   if(fDebug > 1) printf("AliAnalysisTaskPythiaMpi::UserExec \n");

   //
   // Execute analysis for current event
   //
   Init();

   TH1I* fHistEvents = (TH1I*) fOutputList->FindObject("fHistNEvents");
   fHistEvents->Fill(0.5);

  

   if(fMcHandler){
      fMcEvent = fMcHandler->MCEvent(); 

      AliGenEventHeader *fMcHeader = (fMcEvent->GenEventHeader());
      if(!fMcHeader) Printf("fMcHeader NULL!");

      AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (fMcHeader);
      if(!fMcPythiaHeader) Printf("fMcPythiaHeader NULL!");
      //Printf("PROCESS: %d",(fMcPythiaHeader->ProcessType()));

   }else{
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Handler() fMcHandler=NULL\n");
      return;
   }
   if(!fMcEvent){
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Exec()   fMcEvent=NULL \n");
      return;
   }


   Printf("MC particles: %d",fMCEvent->GetNumberOfTracks());

   //MPIs
   TH1F* fHistMpi = (TH1F*) fOutputList->FindObject("fHistMpi");
   Double_t Nmpi = (AliPythia8::Instance())->GetNMPI();
   fHistMpi->Fill(Nmpi);

   TH2F* fHistMultMpi = (TH2F*) fOutputList->FindObject("fHistMultMpi");
   fHistMultMpi->Fill(Nmpi,fMCEvent->GetNumberOfTracks());

   //loop over the tracks 

   TH1F* fHistPt = (TH1F*) fOutputList->FindObject("fHistPt");
   TH1F* fHistEta = (TH1F*) fOutputList->FindObject("fHistEta");
   TH2F* fHistdNdetaMpi = (TH2F*) fOutputList->FindObject("fHistdNdetaMpi");
   TH2F* fHistPtMpi = (TH2F*) fOutputList->FindObject("fHistPtMpi");
   TProfile* fProfileMpiPt = (TProfile*) fOutputList->FindObject("fProfileMpiPt");
   TProfile* fProfilePtMpi = (TProfile*) fOutputList->FindObject("fProfilePtMpi");
   TH3F* fHistTracks = (TH3F*) fOutputList->FindObject("fHistTracks");

   for(Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++){
     AliVParticle* track = fMCEvent->GetTrack(iTracks);
     if(!track){
       Printf("ERROR: Could not receive track %d (mc loop)",iTracks);
       continue;
     }
     
     fHistPt->Fill(track->Pt());
     fHistEta->Fill(track->Eta());

     fHistdNdetaMpi->Fill(Nmpi,track->Eta());
     fHistPtMpi->Fill(Nmpi,track->Pt());
     fProfileMpiPt->Fill(Nmpi,track->Pt());
     fProfilePtMpi->Fill(track->Pt(),Nmpi);

     fHistTracks->Fill(track->Eta(),track->Pt(),Nmpi);
   }
   
   PostData(1, fOutputList);
   return; 
   
}
//____________________________________________________________________________

void AliAnalysisTaskPythiaMpi::Terminate(Option_t*){
  //
  // Terminate analysis
  //
  Info("Terminate","Start and end of Method");
  AliAnalysisTaskSE::Terminate();

  fOutputList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: fOutputList not available\n");
    return;
  }   
    
}

