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
//
// AliAnalysisTaskSE to make AOD centrality
// Author: Alberica Toia, CERN, Alberica.Toia@cern.ch
//
/////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODCentrality.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESD.h"
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
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliAnalysisTaskAODCentralityMaker.h"

ClassImp(AliAnalysisTaskAODCentralityMaker)


//________________________________________________________________________
AliAnalysisTaskAODCentralityMaker::AliAnalysisTaskAODCentralityMaker():
AliAnalysisTaskSE(),
fAODCentrality(0),
fDeltaAODFileName("AliAOD.Centrality.root"),
fIsMCInput        (0),
fNev              (0),
fBeamEnergy       (0),	
fNmyTracks_gen    (0),
fxVertex          (0),
fyVertex          (0),
fzVertex          (0),
fVertexer3d       (0),
fbMC 		  (0),
fNpartTargMC	  (0),
fNpartProjMC	  (0),
fNNColl     	  (0),
fNNwColl    	  (0),
fNwNColl    	  (0),
fNwNwColl   	  (0),
fNTracklets 	  (0),
fNSingleClusters  (0),
fbZDC             (0),
fNpartZDC         (0),
fbZDCA            (0),
fNpartZDCA        (0), 
fbZDCC            (0),   
fNpartZDCC        (0),     
fESDFlag 	  (0),
fZNCEnergy	  (0),
fZPCEnergy	  (0),
fZNAEnergy	  (0),
fZPAEnergy	  (0),
fZEM1Energy	  (0),
fZEM2Energy	  (0),
fNTracks    	  (0),
fNPmdTracks 	  (0),
fMultV0A    	  (0),
fMultV0C    	  (0),
fMultFMDA    	  (0),   
fMultFMDC         (0)
{
  // Default constructor
    
  for (int i=0;i<6;i++) fNClusters[i]=0;
  for (int i=0;i<2;i++) fNChips[i]=0;
  for (int i=0;i<5;i++) fZNCtower[i]=0;
  for (int i=0;i<5;i++) fZPCtower[i]=0;
  for (int i=0;i<5;i++) fZNAtower[i]=0;
  for (int i=0;i<5;i++) fZPAtower[i]=0;
  for (int i=0;i<2;i++) fCentrZNC[i]=0;
  for (int i=0;i<2;i++) fCentrZNA[i]=0;
}

//________________________________________________________________________
AliAnalysisTaskAODCentralityMaker::AliAnalysisTaskAODCentralityMaker(const char *name):
AliAnalysisTaskSE(name),
fAODCentrality(0),
fDeltaAODFileName("AliAOD.Centrality.root"),
fIsMCInput        (0),
fNev              (0),
fBeamEnergy       (0),	
fNmyTracks_gen    (0),
fxVertex          (0),
fyVertex          (0),
fzVertex          (0),
fVertexer3d       (0),
fbMC 		  (0),
fNpartTargMC	  (0),
fNpartProjMC	  (0),
fNNColl     	  (0),
fNNwColl    	  (0),
fNwNColl    	  (0),
fNwNwColl   	  (0),
fNTracklets 	  (0),
fNSingleClusters  (0),
fbZDC             (0),
fNpartZDC         (0),
fbZDCA            (0),
fNpartZDCA        (0), 
fbZDCC            (0),   
fNpartZDCC        (0),     
fESDFlag 	  (0),
fZNCEnergy	  (0),
fZPCEnergy	  (0),
fZNAEnergy	  (0),
fZPAEnergy	  (0),
fZEM1Energy	  (0),
fZEM2Energy	  (0),
fNTracks    	  (0),
fNPmdTracks 	  (0),
fMultV0A    	  (0),
fMultV0C    	  (0),
fMultFMDA    	  (0),   
fMultFMDC         (0)
{
  // Standard constructor
    
  for (int i=0;i<6;i++) fNClusters[i]=0;
  for (int i=0;i<2;i++) fNChips[i]=0;
  for (int i=0;i<5;i++) fZNCtower[i]=0;
  for (int i=0;i<5;i++) fZPCtower[i]=0;
  for (int i=0;i<5;i++) fZNAtower[i]=0;
  for (int i=0;i<5;i++) fZPAtower[i]=0;
  for (int i=0;i<2;i++) fCentrZNC[i]=0;
  for (int i=0;i<2;i++) fCentrZNA[i]=0;
}


//________________________________________________________________________
AliAnalysisTaskAODCentralityMaker::~AliAnalysisTaskAODCentralityMaker()
{
  // Destructor
}  

//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::Init()
{
  // Initialization
  if(fDebug > 1) printf("AnalysisTaskAODCentralityMaker::Init() \n");
  AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFileName.Data());

  return;
}

//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::UserCreateOutputObjects()
{
  
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskAODCentralityMaker::UserCreateOutPutData() \n");
  // Support both the case when the AOD + deltaAOD are produced in an ESD
  // analysis or if the deltaAOD is produced on an analysis on AOD's. (A.G. 27/04/09)
  if(!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }   
  TString filename = fDeltaAODFileName;
  // When running on standard AOD to produce deltas, IsStandardAOD is never set,
  // If AODEvent is NULL, new branches have to be added to the new file(s) (A.G. 15/01/10)
  if(!IsStandardAOD() && AODEvent()) filename = "";

  fAODCentrality = new AliAODCentrality();
  fAODCentrality->SetName("AODCentrality");
  AddAODBranch("AliAODCentrality", &fAODCentrality, filename);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::UserExec(Option_t */*option*/)
{
  AliVEvent* event = InputEvent();
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);

  fBeamEnergy = esd->GetBeamEnergy();
  fNTracks    = event->GetNumberOfTracks();     
  fNPmdTracks = esd->GetNumberOfPmdTracks();     
    
  // ***** V0 info
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  fMultV0A=esdV0->GetMTotV0A();
  fMultV0C=esdV0->GetMTotV0C();
    
  // ***** Trigger selection
  TString triggerClass = esd->GetFiredTriggerClasses();
  sprintf(fTrigClass,"%s",triggerClass.Data());
  
  // ***** vertex info
  const AliESDVertex *vertex = esd->GetPrimaryVertexSPD();
  fxVertex = vertex->GetX();
  fyVertex = vertex->GetY();
  fzVertex = vertex->GetZ();
  if(vertex->IsFromVertexer3D()) fVertexer3d = kTRUE;
  else fVertexer3d = kFALSE;
  Double_t vertex3[3];
  vertex->GetXYZ(vertex3);
  
  // ***** CB info (tracklets, clusters, chips)
  const AliMultiplicity *mult = esd->GetMultiplicity();
  fNTracklets = mult->GetNumberOfTracklets();
  
  for(Int_t ilay=0; ilay<6; ilay++){
    fNClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
  }
  fNSingleClusters = mult->GetNumberOfSingleClusters();
  
  for(Int_t ilay=0; ilay<2; ilay++){
    fNChips[ilay] = mult->GetNumberOfFiredChips(ilay);
  }
  
  // ***** FMD info
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
	  
	  Float_t FMDmult = fmd->Multiplicity(det,ring,sec,strip);
	  if(FMDmult == 0 || FMDmult == AliESDFMD::kInvalidMult) continue;
	  
	  Float_t nParticles=0;
	  
	  if(FMDmult > fFMDLowCut) {
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
  
  // ***** ZDC info
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
  //
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

  // ***** MC info
  if(fIsMCInput){
    
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      printf("   Could not retrieve MC event!!!\n");
      return;
    }
    
    fNmyTracks_gen = 0;
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
      
      fNmyTracks_gen ++;
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
      fNpartTargMC = Int_t (208.-(specNeutronTarg+specProtonTarg));
      fNpartProjMC = Int_t (208.-(specNeutronProj+specProtonProj));
      fNNColl   = ((AliGenHijingEventHeader*) genHeader)->NN();
      fNNwColl  = ((AliGenHijingEventHeader*) genHeader)->NNw();
      fNwNColl  = ((AliGenHijingEventHeader*) genHeader)->NwN();
      fNwNwColl = ((AliGenHijingEventHeader*) genHeader)->NwNw();
    }  
    
  }
  
  fAODCentrality->SetTrigClass         (fTrigClass[100]);
  fAODCentrality->SetxVertex	       (fxVertex       );
  fAODCentrality->SetyVertex	       (fyVertex       );
  fAODCentrality->SetzVertex	       (fzVertex       );
  fAODCentrality->SetVertexer3d        (fVertexer3d    );
  fAODCentrality->SetbMC               (fbMC);
  fAODCentrality->SetNpartTargMC       (fNpartTargMC);
  fAODCentrality->SetNpartProjMC       (fNpartProjMC);
  fAODCentrality->SetNNColl            (fNNColl);
  fAODCentrality->SetNNwColl           (fNNwColl);
  fAODCentrality->SetNwNColl           (fNwNColl);
  fAODCentrality->SetNwNwColl          (fNwNwColl);
  fAODCentrality->SetNTracklets        (fNTracklets);
  fAODCentrality->SetNSingleClusters   (fNSingleClusters);
  fAODCentrality->SetNClusters         (
					   fNClusters[0],
					   fNClusters[1],
					   fNClusters[2],
					   fNClusters[3],
					   fNClusters[4],
					   fNClusters[5]);
  fAODCentrality->SetNChips            (
					   fNChips[0],
					   fNChips[1]);
  fAODCentrality->SetbZDC              (fbZDC);
  fAODCentrality->SetNpartZDC          (fNpartZDC);
  fAODCentrality->SetbZDCA             (fbZDCA);
  fAODCentrality->SetNpartZDCA         (fNpartZDCA);
  fAODCentrality->SetbZDCC             (fbZDCC);
  fAODCentrality->SetNpartZDCC         (fNpartZDCC);
  fAODCentrality->SetESDFlag  	  (fESDFlag);
  fAODCentrality->SetZNCEnergy 	  (fZNCEnergy);
  fAODCentrality->SetZPCEnergy	  (fZPCEnergy);
  fAODCentrality->SetZNAEnergy	  (fZNAEnergy);
  fAODCentrality->SetZPAEnergy	  (fZPAEnergy);
  fAODCentrality->SetZEM1Energy	  (fZEM1Energy);
  fAODCentrality->SetZEM2Energy	  (fZEM2Energy);
  fAODCentrality->SetZNCtower          (
					   fZNCtower[0],
					   fZNCtower[1],
					   fZNCtower[2],
					   fZNCtower[3],
					   fZNCtower[4]);
  fAODCentrality->SetZPCtower          (
					 fZPCtower[0],
					 fZPCtower[1],
					 fZPCtower[2],
					 fZPCtower[3],
					 fZPCtower[4]);
 
  fAODCentrality-> SetZNAtower          (
			     fZNAtower[0],  
			     fZNAtower[1], 
			     fZNAtower[2],  
			     fZNAtower[3],  
			     fZNAtower[4]); 
  fAODCentrality-> SetZPAtower          (
					    fZPAtower[0], 
					    fZPAtower[1], 
					    fZPAtower[2], 
					    fZPAtower[3],
					    fZPAtower[4]);
  fAODCentrality-> SetCentrZNC          (
					    fCentrZNC[0],
					    fCentrZNC[1]);
  fAODCentrality-> SetCentrZNA          (
					    fCentrZNA[0],
					    fCentrZNA[1]);
  fAODCentrality-> SetNTracks           (fNTracks);
  fAODCentrality-> SetNPmdTracks        (fNPmdTracks);
  fAODCentrality-> SetMultV0A           (fMultV0A);
  fAODCentrality-> SetMultV0C           (fMultV0C);
  fAODCentrality-> SetMultFMDA          (fMultFMDA);
  fAODCentrality-> SetMultFMDC          (fMultFMDC);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskAODCentralityMaker: Terminate() \n");
}
