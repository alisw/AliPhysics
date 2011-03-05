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
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAnalysisTaskAODCentralityMaker.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskAODCentralityMaker)


//________________________________________________________________________
AliAnalysisTaskAODCentralityMaker::AliAnalysisTaskAODCentralityMaker():
AliAnalysisTaskSE(),
fAODCentrality(0),
fDeltaAODFileName("AliAOD.Centrality.root"),
fAODHeader        (0),
fIsMCInput        (0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskAODCentralityMaker::AliAnalysisTaskAODCentralityMaker(const char *name):
AliAnalysisTaskSE(name),
fAODCentrality(0),
fDeltaAODFileName("AliAOD.Centrality.root"),
fAODHeader        (0),
fIsMCInput        (0)
{
  // Standard constructor
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
  

  fAODHeader = new AliAODHeader();
  AddAODBranch("AliAODHeader", &fAODHeader, filename);
  return;
}


//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::UserExec(Option_t */*option*/)
{
// User Exec
  AliVEvent*   event = InputEvent();
  AliESDEvent* esd   = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
      AliError("No ESD Event");
      return;
  }
  
  Float_t beamEnergy = esd->GetBeamEnergy();
  Int_t   nTracks    = event->GetNumberOfTracks();     
  Int_t   nPmdTracks = esd->GetNumberOfPmdTracks();     
    
  // ***** V0 info
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  Double_t multV0A = esdV0->GetMTotV0A();
  Double_t multV0C = esdV0->GetMTotV0C();
    
  
  // ***** vertex info
  const AliESDVertex *vertex = esd->GetPrimaryVertexSPD();
  Double_t xVertex = vertex->GetX();
  Double_t yVertex = vertex->GetY();
  Double_t zVertex = vertex->GetZ();
  Bool_t vertexer3d;
  
  if(vertex->IsFromVertexer3D()) vertexer3d = kTRUE;
  else vertexer3d = kFALSE;
  Double_t vertex3[3];
  vertex->GetXYZ(vertex3);
  
  // ***** CB info (tracklets, clusters, chips)
  const AliMultiplicity *mult = esd->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  Int_t nSingleClusters;
  Int_t nClusters[6];
  
  for(Int_t ilay = 0; ilay < 6; ilay++){
    nClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
  }
  nSingleClusters = mult->GetNumberOfSingleClusters();

  Int_t nChips[2];
  for(Int_t ilay = 0; ilay < 2; ilay++){
    nChips[ilay] = mult->GetNumberOfFiredChips(ilay);
  }
  
  // ***** FMD info
  AliESDFMD *fmd = esd->GetFMDData();
  Float_t totalMultA = 0;
  Float_t totalMultC = 0;
  const Float_t fmdLowCut = 0.4;
  
  for(UShort_t det = 1;det <= 3; det++) {
      Int_t nRings = (det==1 ? 1 : 2);
      for (UShort_t ir = 0; ir < nRings; ir++) {	  
	  Char_t   ring = (ir == 0 ? 'I' : 'O');
	  UShort_t nsec = (ir == 0 ? 20  : 40);
	  UShort_t nstr = (ir == 0 ? 512 : 256);
	  for(UShort_t sec =0; sec < nsec;  sec++)  {
	      for(UShort_t strip = 0; strip < nstr; strip++) {
		  Float_t fmdMult = fmd->Multiplicity(det,ring,sec,strip);
		  if(fmdMult == 0 || fmdMult == AliESDFMD::kInvalidMult) continue;
		  Float_t nParticles=0;
		  if(fmdMult > fmdLowCut) {
		      nParticles = 1.;
		  }
	  
		  if (det<3) totalMultA = totalMultA + nParticles;
		  else totalMultC = totalMultC + nParticles;
		  
	      }
	  }
      }
  }
  Float_t multFMDA = totalMultA;
  Float_t multFMDC = totalMultC;
  
  // ***** ZDC info
  AliESDZDC *esdZDC = esd->GetESDZDC();
  UInt_t esdFlag =  esdZDC->GetESDQuality();   
  
  Float_t znCEnergy  = (Float_t) (esdZDC->GetZDCN1Energy());
  Float_t zpCEnergy  = (Float_t) (esdZDC->GetZDCP1Energy());
  Float_t znAEnergy  = (Float_t) (esdZDC->GetZDCN2Energy());
  Float_t zpAEnergy  = (Float_t) (esdZDC->GetZDCP2Energy());
  Float_t zem1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0));
  Float_t zem2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1));
  
  Double_t bZDC      = esdZDC->GetImpactParameter();
  Int_t    nPartZDC  = esdZDC->GetZDCParticipants();
  Double_t bZDCA     = esdZDC->GetImpactParamSideA();
  Int_t    nPartZDCA = esdZDC->GetZDCPartSideA();
  Double_t bZDCC     = esdZDC->GetImpactParamSideC();
  Int_t    nPartZDCC = esdZDC->GetZDCPartSideC();
  
  const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
  const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
  const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
  const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
  //
  Float_t  znCtower[5];	//  ZNC 5 tower signals
  Float_t  zpCtower[5];	//  ZPC 5 tower signals
  Float_t  znAtower[5];	//  ZNA 5 tower signals
  Float_t  zpAtower[5];	//  ZPA 5 tower signals
  Float_t  centrZNC[2];	//  centroid over ZNC
  Float_t  centrZNA[2];	//  centroid over ZNA

  for(Int_t it = 0; it < 5; it++){
    znCtower[it] = (Float_t) (towZNC[it]);
    zpCtower[it] = (Float_t) (towZPC[it]);
    znAtower[it] = (Float_t) (towZNA[it]); 
    zpAtower[it] = (Float_t) (towZPA[it]);  
  }
  
  Double_t xyZNC[2] = {-99.,-99.};
  Double_t xyZNA[2] = {-99.,-99.};

  esdZDC->GetZNCentroidInPbPb(beamEnergy, xyZNC, xyZNA);
  for(Int_t it = 0; it < 2; it++){
      centrZNC[it] = xyZNC[it];
      centrZNA[it] = xyZNA[it];
  }

  // ***** MC info
  Double_t bMC          = 0.;
  Int_t specNeutronProj = 0;
  Int_t specProtonProj  = 0;
  Int_t specNeutronTarg = 0;
  Int_t specProtonTarg  = 0;
  Int_t nPartTargMC     = 0;
  Int_t nPartProjMC     = 0;
  Int_t nnColl          = 0;
  Int_t nnwColl         = 0;
  Int_t nwNColl         = 0;
  Int_t nwNwColl        = 0;

  if(fIsMCInput){
    
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      printf("   Could not retrieve MC event!!!\n");
      return;
    }
    
    Int_t nMyTracks_gen = 0;
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
      
      nMyTracks_gen ++;
    } 
    
    AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
    if(!genHeader){
      printf("  Event generator header not available!!!\n");
      return;
    }
	

    if(genHeader->InheritsFrom(AliGenHijingEventHeader::Class())){
	bMC = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
	specNeutronProj = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsn();
	specProtonProj  = ((AliGenHijingEventHeader*) genHeader)->ProjSpectatorsp();
	specNeutronTarg = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsn();
	specProtonTarg  = ((AliGenHijingEventHeader*) genHeader)->TargSpectatorsp();
	nPartTargMC = Int_t (208.-(specNeutronTarg+specProtonTarg));
	nPartProjMC = Int_t (208.-(specNeutronProj+specProtonProj));
	nnColl   = ((AliGenHijingEventHeader*) genHeader)->NN();
	nnwColl  = ((AliGenHijingEventHeader*) genHeader)->NNw();
	nwNColl  = ((AliGenHijingEventHeader*) genHeader)->NwN();
	nwNwColl = ((AliGenHijingEventHeader*) genHeader)->NwNw();
    }  
  }
  
  fAODCentrality->SetxVertex	       (xVertex       );
  fAODCentrality->SetyVertex	       (yVertex       );
  fAODCentrality->SetzVertex	       (zVertex       );
  fAODCentrality->SetVertexer3d        (vertexer3d    );
  fAODCentrality->SetbMC               (bMC);
  fAODCentrality->SetNpartTargMC       (nPartTargMC);
  fAODCentrality->SetNpartProjMC       (nPartProjMC);
  fAODCentrality->SetNNColl            (nnColl);
  fAODCentrality->SetNNwColl           (nnwColl);
  fAODCentrality->SetNwNColl           (nwNColl);
  fAODCentrality->SetNwNwColl          (nwNwColl);
  fAODCentrality->SetNTracklets        (nTracklets);
  fAODCentrality->SetNSingleClusters   (nSingleClusters);
  fAODCentrality->SetNClusters         (
					   nClusters[0],
					   nClusters[1],
					   nClusters[2],
					   nClusters[3],
					   nClusters[4],
					   nClusters[5]);
  fAODCentrality->SetNChips            (
					   nChips[0],
					   nChips[1]);
  fAODCentrality->SetbZDC              (bZDC);
  fAODCentrality->SetNpartZDC          (nPartZDC);
  fAODCentrality->SetbZDCA             (bZDCA);
  fAODCentrality->SetNpartZDCA         (nPartZDCA);
  fAODCentrality->SetbZDCC             (bZDCC);
  fAODCentrality->SetNpartZDCC         (nPartZDCC);
  fAODCentrality->SetESDFlag  	  (esdFlag);
  fAODCentrality->SetZNCEnergy 	  (znCEnergy);
  fAODCentrality->SetZPCEnergy	  (zpCEnergy);
  fAODCentrality->SetZNAEnergy	  (znAEnergy);
  fAODCentrality->SetZPAEnergy	  (zpAEnergy);
  fAODCentrality->SetZEM1Energy	  (zem1Energy);
  fAODCentrality->SetZEM2Energy	  (zem2Energy);
  fAODCentrality->SetZNCtower          (
					   znCtower[0],
					   znCtower[1],
					   znCtower[2],
					   znCtower[3],
					   znCtower[4]);
  fAODCentrality->SetZPCtower          (
					 zpCtower[0],
					 zpCtower[1],
					 zpCtower[2],
					 zpCtower[3],
					 zpCtower[4]);
 
  fAODCentrality-> SetZNAtower          (
			     znAtower[0],  
			     znAtower[1], 
			     znAtower[2],  
			     znAtower[3],  
			     znAtower[4]); 
  fAODCentrality-> SetZPAtower          (
					    zpAtower[0], 
					    zpAtower[1], 
					    zpAtower[2], 
					    zpAtower[3],
					    zpAtower[4]);
  fAODCentrality-> SetCentrZNC          (
					    centrZNC[0],
					    centrZNC[1]);
  fAODCentrality-> SetCentrZNA          (
					    centrZNA[0],
					    centrZNA[1]);
  fAODCentrality-> SetNTracks           (nTracks);
  fAODCentrality-> SetNPmdTracks        (nPmdTracks);
  fAODCentrality-> SetMultV0A           (multV0A);
  fAODCentrality-> SetMultV0C           (multV0C);
  fAODCentrality-> SetMultFMDA          (multFMDA);
  fAODCentrality-> SetMultFMDC          (multFMDC);

//
// Header Replication
//
  AliAODHeader* hdr = AODEvent()->GetHeader();
  *fAODHeader =  *hdr;
//
  return;
}

//________________________________________________________________________
void AliAnalysisTaskAODCentralityMaker::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskAODCentralityMaker: Terminate() \n");
}
