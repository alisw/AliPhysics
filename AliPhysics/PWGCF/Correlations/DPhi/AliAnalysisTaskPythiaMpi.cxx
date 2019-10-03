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
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"

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
#include "AliStack.h"
//#include "AliPythia.h"
//#include "AliPythia8.h"
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

  TProfile2D* fProfileMpiPt = new TProfile2D("fProfileMpiPt","fProfileMpiPt",100,-0.5,99.5,20,-1.,1.,0.15,10.);
  fProfileMpiPt->GetXaxis()->SetTitle("#MPIs");
  fProfileMpiPt->GetYaxis()->SetTitle("#eta");
  fProfileMpiPt->GetZaxis()->SetTitle("p_{T} (GeV/c)");

  TProfile2D* fProfilePtMpi = new TProfile2D("fProfilePtMpi","fProfilePtMpi",39,0.15,10.,20,-1.,1.,-0.5,99.5);
  fProfilePtMpi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fProfilePtMpi->GetYaxis()->SetTitle("#eta");
  fProfilePtMpi->GetZaxis()->SetTitle("#MPIs");

  TH3F *fHistTracks = new TH3F("fHistTracks","fHistTracks",20,-1.,1.,39,0.15,10.0,100,-0.5,99.5);
  fHistTracks->GetXaxis()->SetTitle("#eta");
  fHistTracks->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fHistTracks->GetZaxis()->SetTitle("#MPIs");

    //4 dim matrix
  Int_t nbins[4]   = {20, 39, 100, 4};
  Double_t minbin[4] = {-1.,           0.15,        -0.5,         1 };
  Double_t maxbin[4] = {1.,            10.0,        99.5,         4 };

  THnSparseD *fHistTracksPID = new THnSparseD("fHistTracksPID","eta:pt:mpi:pid",4,nbins,minbin,maxbin);
  fHistTracksPID->GetAxis(0)->SetTitle("#eta");
  fHistTracksPID->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  fHistTracksPID->GetAxis(2)->SetTitle("#MPI");
  fHistTracksPID->GetAxis(3)->SetTitle("PID");
  fHistTracksPID->Sumw2();

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
  fOutputList->Add(fHistTracksPID);

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

   Int_t Nmpi = 0;
   Int_t particle_id = 0;

   if(fMcHandler){
      fMcEvent = fMcHandler->MCEvent(); 

      AliGenEventHeader *fMcHeader = (fMcEvent->GenEventHeader());
      if(!fMcHeader) Printf("fMcHeader NULL!");

      AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (fMcHeader);
      if(!fMcPythiaHeader) Printf("fMcPythiaHeader NULL!");
      //Printf("PROCESS: %d",(fMcPythiaHeader->ProcessType()));

      //cout<<"TEST: "<<fMcPythiaHeader->GetXJet()<<endl;
   Nmpi = fMcPythiaHeader->GetNMPI();

   }else{
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Handler() fMcHandler=NULL\n");
      return;
   }
   if(!fMcEvent){
      if(fDebug > 1) printf("AnalysisTaskPythiaMpi::Exec()   fMcEvent=NULL \n");
      return;
   }


   Printf("MC particles: %d",fMCEvent->GetNumberOfTracks());
   Printf("MC primary particles: %d",fMCEvent->GetNumberOfPrimaries());

   //MPIs
   TH1F* fHistMpi = (TH1F*) fOutputList->FindObject("fHistMpi");
   fHistMpi->Fill(Nmpi);

   TH2F* fHistMultMpi = (TH2F*) fOutputList->FindObject("fHistMultMpi");
   fHistMultMpi->Fill(Nmpi,fMCEvent->GetNumberOfTracks());

   //loop over the tracks 

   Int_t n_phys_prim_charged = 0;

   TH1F* fHistPt = (TH1F*) fOutputList->FindObject("fHistPt");
   TH1F* fHistEta = (TH1F*) fOutputList->FindObject("fHistEta");
   TH2F* fHistdNdetaMpi = (TH2F*) fOutputList->FindObject("fHistdNdetaMpi");
   TH2F* fHistPtMpi = (TH2F*) fOutputList->FindObject("fHistPtMpi");
   TProfile2D* fProfileMpiPt = (TProfile2D*) fOutputList->FindObject("fProfileMpiPt");
   TProfile2D* fProfilePtMpi = (TProfile2D*) fOutputList->FindObject("fProfilePtMpi");
   TH3F* fHistTracks = (TH3F*) fOutputList->FindObject("fHistTracks");
   THnSparseD* fHistTracksPID = (THnSparseD*) fOutputList->FindObject("fHistTracksPID");

   for(Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++){
     AliVParticle* track = fMCEvent->GetTrack(iTracks);
     if(!track){
       Printf("ERROR: Could not receive track %d (mc loop)",iTracks);
       continue;
     }
     
     if(fMCEvent->IsPhysicalPrimary(iTracks) && track->Charge()!=0){

       fHistPt ->Fill(track->Pt());
       fHistEta->Fill(track->Eta());

       fHistdNdetaMpi->Fill(Nmpi,track->Eta());
       fHistPtMpi    ->Fill(Nmpi,track->Pt());
       fProfileMpiPt ->Fill(Nmpi,track->Eta(),track->Pt());
       fProfilePtMpi ->Fill(track->Pt(),track->Eta(),Nmpi);

       fHistTracks->Fill(track->Eta(),track->Pt(),Nmpi);
       Double_t pid = 0, weight = 1.;

       particle_id = track->PdgCode();
       if(particle_id == 111 || particle_id == 211)
	 pid = 1; //pions
       if(particle_id == 311 || particle_id == 321)
	 pid = 2; //kaons
       if(particle_id == 2212)
	 pid = 3; //protons
       if(particle_id == 3322 || particle_id == 3312)
	 pid = 4; //csis

       Double_t prop[4] = {track->Eta(),track->Pt(),(Double_t)Nmpi,pid};
       fHistTracksPID->Fill(prop, weight);
       n_phys_prim_charged++;

     }
   }
   
   Printf("MC charged physical primaries: %d",n_phys_prim_charged);

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

