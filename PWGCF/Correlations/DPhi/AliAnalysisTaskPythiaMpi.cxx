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

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliAnalysisTaskPythiaMpi.h"
#include "AliGenEventHeader.h"
#include "AliVParticle.h"
#include "AliPythia.h"
#include "AliPythia8.h"

using namespace std;

ClassImp(AliAnalysisTaskPythiaMpi)

//_________________________________________________________| Constructor

AliAnalysisTaskPythiaMpi::AliAnalysisTaskPythiaMpi(): 
  fMcEvent(0x0),  
  fMcHandler(0x0),
  fOutputList(0),
  fHistEvents(0),
  fHistPt(0),
  fHistEta(0),
  fHistMpi(0),
  fHistMultMpi(0),
  fHistdNdetaMpi(0)
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
  fOutputList(0),
  fHistEvents(0),
  fHistPt(0),
  fHistEta(0),
  fHistMpi(0),
  fHistMultMpi(0),
  fHistdNdetaMpi(0)
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
    fHistEvents = c.fHistEvents;
    fHistPt = c.fHistPt;
    fHistEta = c.fHistEta;
    fHistMpi = c.fHistMpi;
    fHistMultMpi = c.fHistMultMpi;
    fHistdNdetaMpi = c.fHistdNdetaMpi;

  }
  return *this;

}
//_________________________________________________________| Copy Constructor

AliAnalysisTaskPythiaMpi::AliAnalysisTaskPythiaMpi(const AliAnalysisTaskPythiaMpi& c):
  AliAnalysisTaskSE(c),
  fMcEvent(c.fMcEvent),
  fMcHandler(c.fMcHandler), 
  fOutputList(c.fOutputList),
  fHistEvents(c.fHistEvents),
  fHistPt(c.fHistPt),
  fHistEta(c.fHistEta),
  fHistMpi(c.fHistMpi),
  fHistMultMpi(c.fHistMultMpi),
  fHistdNdetaMpi(c.fHistdNdetaMpi)
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
  if(fHistEvents) delete fHistEvents;
  if(fHistPt) delete fHistPt;
  if(fHistEta) delete fHistEta;
  if(fHistMpi) delete fHistMpi;
  if(fHistMultMpi) delete fHistMultMpi;
  if(fHistdNdetaMpi) delete fHistdNdetaMpi;
  
}

//_____________________________________________________________________

void AliAnalysisTaskPythiaMpi::UserCreateOutputObjects()
{
  Info("CreateOutputObjects","CreateOutputObjects of task %s",GetName());
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
 
  fHistEvents = new TH1I("fHistNEvents","fHistNEvents",2,0,2);
  fHistEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"Events with |zVtx|<10cm");

  fHistPt = new TH1F("fHistPt","p_{T} distribution",15,0.2,5.0);
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c^{-1})");
  fHistPt->SetMarkerStyle(kFullCircle);

  fHistEta = new TH1F("fHistEta","eta distribution",20,-1,1);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("dN/d#eta");
  fHistEta->SetMarkerStyle(kFullCircle);

  fHistMpi = new TH1F("fHistMpi","MPIs distribution",50,0,50);
  fHistMpi->GetXaxis()->SetTitle("#MPIs");
  fHistMpi->GetYaxis()->SetTitle("dN/dMPIs");
  fHistMpi->SetMarkerStyle(kFullCircle);

  fHistMultMpi = new TH2F("fHistMultMpi","Multiplicity vs MPIs",50,0,50,100,0,1e4);
  fHistMultMpi->GetXaxis()->SetTitle("#MPIs");
  fHistMultMpi->GetYaxis()->SetTitle("MC particles");
  fHistMultMpi->SetMarkerStyle(kFullDotSmall);

  fHistdNdetaMpi = new TH2F("fHistdNdetaMpi","dN/d#etadMPIs",50,0,50,100,-5,5);
  fHistdNdetaMpi->GetXaxis()->SetTitle("#MPIs");
  fHistdNdetaMpi->GetYaxis()->SetTitle("#eta");
  //fHistdNdetaMpi->SetMarkerStyle();


  fOutputList->Add(fHistEvents);
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistMpi);
  fOutputList->Add(fHistMultMpi);
  fOutputList->Add(fHistdNdetaMpi);

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
   fHistEvents->Fill(0.5);
   if(fMcHandler){
      fMcEvent = fMcHandler->MCEvent(); 
   }else{
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Handler() fMcHandler=NULL\n");
      return;
   }
   if(!fMcEvent){
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Exec()   fMcEvent=NULL \n");
      return;
   }

   const AliVVertex *vtxMC = fMcEvent->GetPrimaryVertex();
   Float_t zVtx = vtxMC->GetZ();
   if(TMath::Abs(zVtx)<10) fHistEvents->Fill(1.5);
   else return;

   Printf("MC particles: %d",fMCEvent->GetNumberOfTracks());


   //MPIs
   Double_t Nmpi = (AliPythia8::Instance())->GetNMPI();
   fHistMpi->Fill(Nmpi);

   fHistMultMpi->Fill(Nmpi,fMCEvent->GetNumberOfTracks());

   //loop over the tracks 

   for(Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++){
     AliVParticle* track = fMCEvent->GetTrack(iTracks);
     if(!track){
       Printf("ERROR: Could not receive track %d (mc loop)",iTracks);
       continue;
     }
     fHistPt->Fill(track->Pt());
     fHistEta->Fill(track->Eta());

     fHistdNdetaMpi->Fill(Nmpi,track->Eta());
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

