/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//							   //
//	Class to analyze centrality measurements           //
//							   //
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliAnalysisTaskCentrality.h"

ClassImp(AliAnalysisTaskCentrality)


//________________________________________________________________________
AliAnalysisTaskCentrality::AliAnalysisTaskCentrality():
  AliAnalysisTaskSE(),
  fDebug(0),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fOutput(0x0),
  fhEzdc(0x0),
  fhEzem(0x0),
  fhNtracks(0x0),
  fhNtracklets(0x0),
  fhNclusters0(0x0),
  fhmultV0(0x0),
  fhmultFMD(0x0),
  fhEzemvsEzdc(0x0),
  fhNtracksvsEzdc(0x0),
  fhNtrackletsvsEzdc(0x0),
  fhNclusters0vsEzdc(0x0),
  fhmultV0vsEzdc(0x0),
  fhmultFMDvsEzdc(0x0),
  fhNtracksvsEzem(0x0),
  fhNtrackletsvsEzem(0x0),
  fhNclusters0vsEzem(0x0),
  fhmultV0vsEzem(0x0),
  fhmultFMDvsEzem(0x0),
  fhNtracksvsmultV0(0x0),
  fhNtrackletsvsmultV0(0x0),
  fhNclusters0vsmultV0(0x0),
  fhNtracksvsmultFMD(0x0),
  fhNtrackletsvsmultFMD(0x0),
  fhNclusters0vsmultFMD(0x0),
  fhmultV0vsmultFMD(0x0),
  fNev(0),	
  fBeamEnergy(0),	
  fNmyTracksgen(0),
  fxVertex(0),	
  fyVertex(0),	
  fzVertex(0),	
  fVertexer3d(0),	
  fbMC(0),	
  fNpartTargMC(0),
  fNpartProjMC(0),
  fNNColl(0),     
  fNNwColl(0),    
  fNwNColl(0),    
  fNwNwColl(0),   
  fNTracklets(0),	
  fNSingleClusters(0),
  fbZDC(0),
  fNpartZDC(0),
  fbZDCA(0),
  fNpartZDCA(0),
  fbZDCC(0),
  fNpartZDCC(0),
  fESDFlag(0),
  fZNCEnergy(0),
  fZPCEnergy(0),
  fZNAEnergy(0),
  fZPAEnergy(0),
  fZEM1Energy(0),
  fZEM2Energy(0),
  fNTracks(0),
  fNPmdTracks(0),
  fMultV0A(0),
  fMultV0C(0),
  fMultFMDA(0),
  fMultFMDC(0)
{   
   // Default constructor
}   

//________________________________________________________________________
AliAnalysisTaskCentrality::AliAnalysisTaskCentrality(const char *name):
  AliAnalysisTaskSE(name),
  fDebug(0),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fOutput(0x0),
  fhEzdc(0x0),
  fhEzem(0x0),
  fhNtracks(0x0),
  fhNtracklets(0x0),
  fhNclusters0(0x0),
  fhmultV0(0x0),
  fhmultFMD(0x0),
  fhEzemvsEzdc(0x0),
  fhNtracksvsEzdc(0x0),
  fhNtrackletsvsEzdc(0x0),
  fhNclusters0vsEzdc(0x0),
  fhmultV0vsEzdc(0x0),
  fhmultFMDvsEzdc(0x0),
  fhNtracksvsEzem(0x0),
  fhNtrackletsvsEzem(0x0),
  fhNclusters0vsEzem(0x0),
  fhmultV0vsEzem(0x0),
  fhmultFMDvsEzem(0x0),
  fhNtracksvsmultV0(0x0),
  fhNtrackletsvsmultV0(0x0),
  fhNclusters0vsmultV0(0x0),
  fhNtracksvsmultFMD(0x0),
  fhNtrackletsvsmultFMD(0x0),
  fhNclusters0vsmultFMD(0x0),
  fhmultV0vsmultFMD(0x0),
  fNev(0),	
  fBeamEnergy(0),	
  fNmyTracksgen(0),
  fxVertex(0),	
  fyVertex(0),	
  fzVertex(0),	
  fVertexer3d(0),	
  fbMC(0),	
  fNpartTargMC(0),
  fNpartProjMC(0),
  fNNColl(0),     
  fNNwColl(0),    
  fNwNColl(0),    
  fNwNwColl(0),   
  fNTracklets(0),	
  fNSingleClusters(0),
  fbZDC(0),
  fNpartZDC(0),
  fbZDCA(0),
  fNpartZDCA(0),
  fbZDCC(0),
  fNpartZDCC(0),
  fESDFlag(0),
  fZNCEnergy(0),
  fZPCEnergy(0),
  fZNAEnergy(0),
  fZPAEnergy(0),
  fZEM1Energy(0),
  fZEM2Energy(0),
  fNTracks(0),
  fNPmdTracks(0),
  fMultV0A(0),
  fMultV0C(0),
  fMultFMDA(0),
  fMultFMDC(0)
{
  // Default constructor
  
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); 

}

//________________________________________________________________________
AliAnalysisTaskCentrality& AliAnalysisTaskCentrality::operator=(const AliAnalysisTaskCentrality& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
  }
  return *this;
}

//________________________________________________________________________
AliAnalysisTaskCentrality::AliAnalysisTaskCentrality(const AliAnalysisTaskCentrality& ana):
  AliAnalysisTaskSE(ana),
  fDebug(ana.fDebug),	  
  fAnalysisInput(ana.fDebug),
  fIsMCInput(ana.fIsMCInput),
  fOutput(ana.fOutput),
  fhEzdc(ana.fhEzdc),
  fhEzem(ana.fhEzem),
  fhNtracks(ana.fhNtracks),
  fhNtracklets(ana.fhNtracklets),
  fhNclusters0(ana.fhNclusters0),
  fhmultV0(ana.fhmultV0),
  fhmultFMD(ana.fhmultFMD),
  fhEzemvsEzdc(ana.fhEzemvsEzdc),
  fhNtracksvsEzdc(ana.fhNtracksvsEzdc),
  fhNtrackletsvsEzdc(ana.fhNtrackletsvsEzdc),
  fhNclusters0vsEzdc(ana.fhNclusters0vsEzdc),
  fhmultV0vsEzdc(ana.fhmultV0vsEzdc),
  fhmultFMDvsEzdc(ana.fhmultFMDvsEzdc),
  fhNtracksvsEzem(ana.fhNtracksvsEzem),
  fhNtrackletsvsEzem(ana.fhNtrackletsvsEzem),
  fhNclusters0vsEzem(ana.fhNclusters0vsEzem),
  fhmultV0vsEzem(ana.fhmultV0vsEzem),
  fhmultFMDvsEzem(ana.fhmultFMDvsEzem),
  fhNtracksvsmultV0(ana.fhNtracksvsmultV0),
  fhNtrackletsvsmultV0(ana.fhNtrackletsvsmultV0),
  fhNclusters0vsmultV0(ana.fhNclusters0vsmultV0),
  fhNtracksvsmultFMD(ana.fhNtracksvsmultFMD),
  fhNtrackletsvsmultFMD(ana.fhNtrackletsvsmultFMD),
  fhNclusters0vsmultFMD(ana.fhNclusters0vsmultFMD),
  fhmultV0vsmultFMD(ana.fhmultV0vsmultFMD),
  fNev(ana.fNev),	
  fBeamEnergy(ana.fBeamEnergy),	
  fNmyTracksgen(ana.fNmyTracksgen),
  fxVertex(ana.fxVertex),	
  fyVertex(ana.fyVertex),	
  fzVertex(ana.fzVertex),	
  fVertexer3d(ana.fVertexer3d),	
  fbMC(ana.fbMC),	
  fNpartTargMC(ana.fNpartTargMC),
  fNpartProjMC(ana.fNpartProjMC),
  fNNColl(ana.fNNColl),     
  fNNwColl(ana.fNNwColl),    
  fNwNColl(ana.fNwNColl),    
  fNwNwColl(ana.fNwNwColl),   
  fNTracklets(ana.fNTracklets),	
  fNSingleClusters(ana.fNSingleClusters),
  fbZDC(ana.fbZDC),
  fNpartZDC(ana.fNpartZDC),
  fbZDCA(ana.fbZDCA),
  fNpartZDCA(ana.fNpartZDCA),
  fbZDCC(ana.fbZDCC),
  fNpartZDCC(ana.fNpartZDCC),
  fESDFlag(ana.fESDFlag),
  fZNCEnergy(ana.fZNCEnergy),
  fZPCEnergy(ana.fZPCEnergy),
  fZNAEnergy(ana.fZNAEnergy),
  fZPAEnergy(ana.fZPAEnergy),
  fZEM1Energy(ana.fZEM1Energy),
  fZEM2Energy(ana.fZEM2Energy),
  fNTracks(ana.fNTracks),
  fNPmdTracks(ana.fNPmdTracks),
  fMultV0A(ana.fMultV0A),
  fMultV0C(ana.fMultV0C),
  fMultFMDA(ana.fMultFMDA),
  fMultFMDC(ana.fMultFMDC)
{
  //
  // Copy Constructor	
  //
}
 
//________________________________________________________________________
 AliAnalysisTaskCentrality::~AliAnalysisTaskCentrality()
 {
   // Destructor
   if(fOutput){
     delete fOutput; fOutput=0;
   } 
 }  

//________________________________________________________________________
void AliAnalysisTaskCentrality::UserCreateOutputObjects()
{  

  // Create the output containers
  if(fDebug>1) printf("AnalysisTaskZDCpp::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fhEzdc         = new TH1F("hEzdc","hEzdc",500,0,150);
  fhEzem         = new TH1F("hEzem","hEzem",500,0,5);
  fhNtracks      = new TH1F("hNtracks","hNtracks",500,0,17000);
  fhNtracklets   = new TH1F("hNtracklets","hNtracklets",500,0,10000);
  fhNclusters0   = new TH1F("hNclusters0","hNclusters0",500,0,15000);
  fhmultV0       = new TH1F("hmultV0","hmultV0",500,0,30000);
  fhmultFMD      = new TH1F("hmultFMD","hmultFMD",500,0,24000);

  fhEzemvsEzdc         = new TProfile("hEzemvsEzdc","hEzemvsEzdc",500,0,5,"");
  fhNtracksvsEzdc      = new TProfile("hNtracksvsEzdc","hNtracksvsEzdc",500,0,17000,"");
  fhNtrackletsvsEzdc   = new TProfile("hNtrackletsvsEzdc","hNtrackletsvsEzdc",500,0,10000,"");
  fhNclusters0vsEzdc   = new TProfile("hNclusters0vsEzdc","hNclusters0vsEzdc",500,0,15000,"");
  fhmultV0vsEzdc       = new TProfile("hmultV0vsEzdc","hmultV0vsEzdc",500,0,30000,"");
  fhmultFMDvsEzdc      = new TProfile("hmultFMDvsEzdc","hmultFMDvsEzdc",500,0,24000,"");
  fhNtracksvsEzem      = new TProfile("hNtracksvsEzem","hNtracksvsEzem",500,0,17000,"");
  fhNtrackletsvsEzem   = new TProfile("hNtrackletsvsEzem","hNtrackletsvsEzem",500,0,10000,"");
  fhNclusters0vsEzem   = new TProfile("hNclusters0vsEzem","hNclusters0vsEzem",500,0,15000,"");
  fhmultV0vsEzem       = new TProfile("hmultV0vsEzem","hmultV0vsEzem",500,0,30000,"");
  fhmultFMDvsEzem      = new TProfile("hmultFMDvsEzem","hmultFMDvsEzem",500,0,24000,"");
  fhNtracksvsmultV0    = new TProfile("hNtracksvsmultV0","hNtracksvsmultV0",500,0,17000,"");      
  fhNtrackletsvsmultV0 = new TProfile("hNtrackletsvsmultV0","hNtrackletsvsmultV0",500,0,10000,"");    
  fhNclusters0vsmultV0 = new TProfile("hNclusters0vsmultV0","hNclusters0vsmultV0",500,0,15000,"");
  fhNtracksvsmultFMD   = new TProfile("hNtracksvsmultFMD","hNtracksvsmultFMD",500,0,17000,"");
  fhNtrackletsvsmultFMD= new TProfile("hNtrackletsvsmultFMD","hNtrackletsvsmultFMD",500,0,10000,"");
  fhNclusters0vsmultFMD= new TProfile("hNclusters0vsmultFMD","hNclusters0vsmultFMD",500,0,15000,"");		   
  fhmultV0vsmultFMD    = new TProfile("hmultV0vsmultFMD","hmultV0vsmultFMD",500,0,30000,"");

  fhEzdc         ->GetXaxis()->SetTitle("E_{ZDC}[TeV]");
  fhEzem         ->GetXaxis()->SetTitle("E_{ZEM}[TeV]");
  fhNtracks      ->GetXaxis()->SetTitle("N_{tracks}");
  fhNtracklets   ->GetXaxis()->SetTitle("N_{tracklets}");
  fhNclusters0   ->GetXaxis()->SetTitle("N_{clusters0}");
  fhmultV0       ->GetXaxis()->SetTitle("V0 mult");
  fhmultFMD      ->GetXaxis()->SetTitle("FMD mult");
  
  fhEzemvsEzdc         ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhNtracksvsEzdc      ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhNtrackletsvsEzdc   ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhNclusters0vsEzdc   ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhmultV0vsEzdc       ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhmultFMDvsEzdc      ->GetYaxis()->SetTitle("E_{ZDC}[TeV]");
  fhNtracksvsEzem      ->GetYaxis()->SetTitle("E_{ZEM}[TeV]");
  fhNtrackletsvsEzem   ->GetYaxis()->SetTitle("E_{ZEM}[TeV]");
  fhNclusters0vsEzem   ->GetYaxis()->SetTitle("E_{ZEM}[TeV]");
  fhmultV0vsEzem       ->GetYaxis()->SetTitle("E_{ZEM}[TeV]");
  fhmultFMDvsEzem      ->GetYaxis()->SetTitle("E_{ZEM}[TeV]");
  fhNtracksvsmultV0    ->GetYaxis()->SetTitle("V0 mult");    
  fhNtrackletsvsmultV0 ->GetYaxis()->SetTitle("V0 mult");  
  fhNclusters0vsmultV0 ->GetYaxis()->SetTitle("V0 mult");
  fhNtracksvsmultFMD   ->GetYaxis()->SetTitle("FMD mult");
  fhNtrackletsvsmultFMD->GetYaxis()->SetTitle("FMD mult");
  fhNclusters0vsmultFMD->GetYaxis()->SetTitle("FMD mult");
  fhmultV0vsmultFMD    ->GetYaxis()->SetTitle("FMD mult");
  
  fhEzemvsEzdc         ->GetXaxis()->SetTitle("E_{ZEM}[TeV]");
  fhNtracksvsEzdc      ->GetXaxis()->SetTitle("N_{tracks}");
  fhNtrackletsvsEzdc   ->GetXaxis()->SetTitle("N_{tracklets}");
  fhNclusters0vsEzdc   ->GetXaxis()->SetTitle("N_{clusters0}");
  fhmultV0vsEzdc       ->GetXaxis()->SetTitle("V0 mult");
  fhmultFMDvsEzdc      ->GetXaxis()->SetTitle("FMD mult");
  fhNtracksvsEzem      ->GetXaxis()->SetTitle("N_{tracks}");
  fhNtrackletsvsEzem   ->GetXaxis()->SetTitle("N_{tracklets}");
  fhNclusters0vsEzem   ->GetXaxis()->SetTitle("N_{clusters0}");
  fhmultV0vsEzem       ->GetXaxis()->SetTitle("V0 mult");
  fhmultFMDvsEzem      ->GetXaxis()->SetTitle("FMD mult");
  fhNtracksvsmultV0    ->GetXaxis()->SetTitle("N_{tracks}");    
  fhNtrackletsvsmultV0 ->GetXaxis()->SetTitle("N_{tracklets}");  
  fhNclusters0vsmultV0 ->GetXaxis()->SetTitle("N_{clusters0}");
  fhNtracksvsmultFMD   ->GetXaxis()->SetTitle("N_{tracks}");
  fhNtrackletsvsmultFMD->GetXaxis()->SetTitle("N_{tracklets}");
  fhNclusters0vsmultFMD->GetXaxis()->SetTitle("N_{clusters}");
  fhmultV0vsmultFMD    ->GetXaxis()->SetTitle("V0 mult");
  
  fOutput->Add(fhEzdc);
  fOutput->Add(fhEzem);
  fOutput->Add(fhNtracks);
  fOutput->Add(fhNtracklets);
  fOutput->Add(fhNclusters0);
  fOutput->Add(fhmultV0);
  fOutput->Add(fhmultFMD);

  fOutput->Add(fhEzemvsEzdc);
  fOutput->Add(fhNtracksvsEzdc);
  fOutput->Add(fhNtrackletsvsEzdc);
  fOutput->Add(fhNclusters0vsEzdc);
  fOutput->Add(fhmultV0vsEzdc);
  fOutput->Add(fhmultFMDvsEzdc);
  fOutput->Add(fhNtracksvsEzem);
  fOutput->Add(fhNtrackletsvsEzem);
  fOutput->Add(fhNclusters0vsEzem);
  fOutput->Add(fhmultV0vsEzem);
  fOutput->Add(fhmultFMDvsEzem);
  fOutput->Add(fhNtracksvsmultV0);
  fOutput->Add(fhNtrackletsvsmultV0);
  fOutput->Add(fhNclusters0vsmultV0);
  fOutput->Add(fhNtracksvsmultFMD);
  fOutput->Add(fhNtrackletsvsmultFMD);
  fOutput->Add(fhNclusters0vsmultFMD);
  fOutput->Add(fhmultV0vsmultFMD);
  
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskCentrality::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskCentrality::UserExec() \n");
  
  if(fAnalysisInput.CompareTo("ESD")==0){

    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (!esd) return;
      fNev++;

      fNTracks    = event->GetNumberOfTracks();     
      fNPmdTracks = esd->GetNumberOfPmdTracks();     

      AliESDVZERO* esdV0 = esd->GetVZEROData();
      fMultV0A=esdV0->GetMTotV0A();
      fMultV0C=esdV0->GetMTotV0C();

      if(fIsMCInput){

        AliMCEvent* mcEvent = MCEvent();
        if (!mcEvent) {
          printf("   Could not retrieve MC event!!!\n");
          return;
        }

	fNmyTracksgen = 0;
	AliStack *stack = 0x0; // needed for MC studies
	stack = MCEvent()->Stack();
	for (Int_t iTrack = 0; iTrack < MCEvent()->GetNumberOfTracks(); iTrack++) {
	  //get properties of mc particle
	  AliMCParticle* mcP = (AliMCParticle*) MCEvent()->GetTrack(iTrack);
	  // Primaries only
	  if (!(stack->IsPhysicalPrimary(mcP->Label()))) continue;
	  //charged tracks only
	  if (mcP->Particle()->GetPDG()->Charge() == 0) continue;
	  //same cuts as on ESDtracks
// 	  if(TMath::Abs(mcP->Eta())>0.9)continue;
// 	  if(mcP->Pt()<0.2)continue;
// 	  if(mcP->Pt()>200)continue;

	  fNmyTracksgen ++;
	} 

        AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
        if(!genHeader){
          printf("  Event generator header not available!!!\n");
	  return;
        }
	
	if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
          fbMC = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
          Int_t specNeutronProj = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsn();
          Int_t specProtonProj  = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsp();
          Int_t specNeutronTarg = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsn();
          Int_t specProtonTarg  = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsp();
	  fNpartTargMC = (Int_t)(208.-(specNeutronTarg+specProtonTarg));
	  fNpartProjMC = (Int_t)(208.-(specNeutronProj+specProtonProj));
	  fNNColl   = ((AliGenHijingEventHeader*) genHeader)->NN();
	  fNNwColl  = ((AliGenHijingEventHeader*) genHeader)->NNw();
	  fNwNColl  = ((AliGenHijingEventHeader*) genHeader)->NwN();
	  fNwNwColl = ((AliGenHijingEventHeader*) genHeader)->NwNw();
	}  
	
      }
      
      fBeamEnergy = esd->GetBeamEnergy();

      // ***** Trigger selection
      TString triggerClass = esd->GetFiredTriggerClasses();
      sprintf(fTrigClass,"%s",triggerClass.Data());
          
      const AliESDVertex *vertex = esd->GetPrimaryVertexSPD();
      fxVertex = vertex->GetX();
      fyVertex = vertex->GetY();
      fzVertex = vertex->GetZ();
      if(vertex->IsFromVertexer3D()) fVertexer3d = kTRUE;
      else fVertexer3d = kFALSE;
      Double_t vertex3[3];
      vertex->GetXYZ(vertex3);

      const AliMultiplicity *mult = esd->GetMultiplicity();
      fNTracklets = mult->GetNumberOfTracklets();
     
      for(Int_t ilay=0; ilay<6; ilay++){
        fNClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
      }
      fNSingleClusters = mult->GetNumberOfSingleClusters();

      for(Int_t ilay=0; ilay<2; ilay++){
        fNChips[ilay] = mult->GetNumberOfFiredChips(ilay);
      }


      AliESDFMD *fmd = esd->GetFMDData();
      Float_t totalMultA = 0;
      Float_t totalMultC = 0;
      const Float_t fFMDLowCut = 0.4;
      
      for(UShort_t det=1;det<=3;det++) {
  	Int_t nRings = (det==1 ? 1 : 2);
  	for (UShort_t ir = 0; ir < nRings; ir++) {	  
  	  Char_t   ring = (ir == 0 ? 'I' : 'O');
  	  UShort_t nsec = (ir == 0 ? 20  : 40);
  	  UShort_t nstr = (ir == 0 ? 512 : 256);
	  for(UShort_t sec =0; sec < nsec;  sec++)  {
  	    for(UShort_t strip = 0; strip < nstr; strip++) {

	      Float_t fmdmult = fmd->Multiplicity(det,ring,sec,strip);
	      if(fmdmult == 0 || fmdmult == AliESDFMD::kInvalidMult) continue;

	      Float_t nParticles=0;
		
		if(fmdmult > fFMDLowCut) {
		  nParticles = 1.;
		}
	      
	      if (det<3) totalMultA = totalMultA + nParticles;
	      else totalMultC = totalMultC + nParticles;
	      
	    }
	  }
	}
      }
      fMultFMDA = totalMultA;
      fMultFMDC = totalMultC;

      AliESDZDC *esdZDC = esd->GetESDZDC();
      fESDFlag =  esdZDC->GetESDQuality();   
      fZNCEnergy = (Float_t) (esdZDC->GetZDCN1Energy());
      fZPCEnergy = (Float_t) (esdZDC->GetZDCP1Energy());
      fZNAEnergy = (Float_t) (esdZDC->GetZDCN2Energy());
      fZPAEnergy = (Float_t) (esdZDC->GetZDCP2Energy());
      fZEM1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0));
      fZEM2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1));
            
      fbZDC = esdZDC->GetImpactParameter();
      fNpartZDC = esdZDC->GetZDCParticipants();
      fbZDCA = esdZDC->GetImpactParamSideA();
      fNpartZDCA = esdZDC->GetZDCPartSideA();
      fbZDCC = esdZDC->GetImpactParamSideC();
      fNpartZDCC = esdZDC->GetZDCPartSideC();
      
      const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
      const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
      const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
      const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();

      for(Int_t it=0; it<5; it++){
         fZNCtower[it] = (Float_t) (towZNC[it]);
         fZPCtower[it] = (Float_t) (towZPC[it]);
         fZNAtower[it] = (Float_t) (towZNA[it]); 
         fZPAtower[it] = (Float_t) (towZPA[it]);  
      }

      Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
      esdZDC->GetZNCentroidInPbPb(fBeamEnergy, xyZNC, xyZNA);
      for(Int_t it=0; it<2; it++){
         fCentrZNC[it] = xyZNC[it];
         fCentrZNA[it] = xyZNA[it];
      }

      // filling histos
      fhEzdc         ->Fill((fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhEzem         ->Fill(fZEM1Energy+fZEM2Energy);
      fhNtracks      ->Fill(fNTracks);
      fhNtracklets   ->Fill(fNTracklets);
      fhNclusters0   ->Fill(fNClusters[0]);
      fhmultV0       ->Fill(fMultV0A+fMultV0C);
      fhmultFMD      ->Fill(fMultFMDA+fMultFMDC);
      fhEzemvsEzdc         ->Fill(fZEM1Energy+fZEM2Energy, (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhNtracksvsEzdc      ->Fill(fNTracks, (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhNtrackletsvsEzdc   ->Fill(fNTracklets,  (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhNclusters0vsEzdc   ->Fill(fNClusters[0],  (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhmultV0vsEzdc       ->Fill(fMultV0A+fMultV0C,  (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhmultFMDvsEzdc      ->Fill(fMultFMDA+fMultFMDC,  (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fhNtracksvsEzem      ->Fill(fNTracks, fZEM1Energy+fZEM2Energy);
      fhNtrackletsvsEzem   ->Fill(fNTracklets, fZEM1Energy+fZEM2Energy);
      fhNclusters0vsEzem   ->Fill(fNClusters[0], fZEM1Energy+fZEM2Energy);
      fhmultV0vsEzem       ->Fill(fMultV0A+fMultV0C, fZEM1Energy+fZEM2Energy);
      fhmultFMDvsEzem      ->Fill(fMultFMDA+fMultFMDC, fZEM1Energy+fZEM2Energy);
      fhNtracksvsmultV0    ->Fill(fNTracks,fMultV0A+fMultV0C);    
      fhNtrackletsvsmultV0 ->Fill(fNTracklets,fMultV0A+fMultV0C);    
      fhNclusters0vsmultV0 ->Fill(fNClusters[0],fMultV0A+fMultV0C);    
      fhNtracksvsmultFMD   ->Fill(fNTracks,fMultFMDA+fMultFMDC);
      fhNtrackletsvsmultFMD->Fill(fNTracklets,fMultFMDA+fMultFMDC);
      fhNclusters0vsmultFMD->Fill(fNClusters[0],fMultFMDA+fMultFMDC);
      fhmultV0vsmultFMD    ->Fill(fMultV0A+fMultV0C,fMultFMDA+fMultFMDC);
  }   
  else if(fAnalysisInput.CompareTo("AOD")==0){
    //AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    // to be implemented
    printf("  AOD analysis not yet implemented!!!\n\n");
    return;
  }
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskCentrality::Terminate(Option_t */*option*/)
{
  // Terminate analysis
}
