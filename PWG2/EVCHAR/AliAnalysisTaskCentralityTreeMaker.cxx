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
#include <TH2F.h>
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
#include "AliAnalysisTaskCentralityTreeMaker.h"

ClassImp(AliAnalysisTaskCentralityTreeMaker)


//________________________________________________________________________
AliAnalysisTaskCentralityTreeMaker::AliAnalysisTaskCentralityTreeMaker():
  AliAnalysisTaskSE(),
    fDebug(0),
    fAnalysisInput("ESD"),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fEZDCvsEZEM(0x0),
    fEZDCvsNtracklets(0x0),
    fTreeFilling(kTRUE),
    fCentralityTree(0x0),
    fNev(0),
    fBeamEnergy(0), 
    fNmyTracks_gen(0),
    fxVertex(0),	 
    fyVertex(0),	 
    fzVertex(0),	 
    fVertexer3d(kFALSE),
    fbMC(0.),
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
AliAnalysisTaskCentralityTreeMaker::AliAnalysisTaskCentralityTreeMaker(const char *name):
  AliAnalysisTaskSE(name),
    fDebug(0),
    fAnalysisInput("ESD"),
    fIsMCInput(kFALSE),
    fOutput(0x0),
    fEZDCvsEZEM(0x0),
    fEZDCvsNtracklets(0x0),
    fTreeFilling(kTRUE),
    fCentralityTree(0x0),
    fNev(0),
    fBeamEnergy(0), 
    fNmyTracks_gen(0),
    fxVertex(0),	 
    fyVertex(0),	 
    fzVertex(0),	 
    fVertexer3d(kFALSE), 
    fbMC(0.),
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
  // Output slot #2 writes into a TTree container
  if(fTreeFilling) DefineOutput(2, TTree::Class()); 

}

//________________________________________________________________________
AliAnalysisTaskCentralityTreeMaker& AliAnalysisTaskCentralityTreeMaker::operator=(const AliAnalysisTaskCentralityTreeMaker& c)
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
AliAnalysisTaskCentralityTreeMaker::AliAnalysisTaskCentralityTreeMaker(const AliAnalysisTaskCentralityTreeMaker& ana):
  AliAnalysisTaskSE(ana),
  fDebug(ana.fDebug),	  
  fAnalysisInput(ana.fDebug),
  fIsMCInput(ana.fIsMCInput),
  fOutput(ana.fOutput),
  fEZDCvsEZEM(ana.fEZDCvsEZEM),
  fEZDCvsNtracklets(ana.fEZDCvsNtracklets),
  fTreeFilling(ana.fTreeFilling),
  fCentralityTree(ana.fCentralityTree),
  fNev(ana.fDebug),        
  fBeamEnergy(ana.fBeamEnergy), 
  fNmyTracks_gen(ana.fNmyTracks_gen),
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
 AliAnalysisTaskCentralityTreeMaker::~AliAnalysisTaskCentralityTreeMaker()
 {
   // Destructor
   if(fOutput){
     delete fOutput; fOutput=0;
   } 
   if(fCentralityTree){
     delete fCentralityTree; fCentralityTree=0;
   } 
 
 }  

//________________________________________________________________________
void AliAnalysisTaskCentralityTreeMaker::UserCreateOutputObjects()
{  

  // Create the output containers
  if(fDebug>1) printf("AnalysisTaskZDCpp::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fEZDCvsEZEM = new TH2F("fEZDCvsEZEM", "E_{ZDC} vs. E_{ZEM}", 100,0.,600.,100,0.,5000.);
  fEZDCvsEZEM->GetXaxis()->SetTitle("ZEM signal (GeV)");
  fEZDCvsEZEM->GetYaxis()->SetTitle("ZDC signal (TeV)");
  fOutput->Add(fEZDCvsEZEM);      
  
  fEZDCvsNtracklets = new TH2F("fEZDCvsNtracklets", "E_{ZDC} vs. N_{tracklets}", 100,0.,600.,100,0.,6000.);  
  fEZDCvsNtracklets->GetXaxis()->SetTitle("N_{tracklets}");
  fEZDCvsNtracklets->GetYaxis()->SetTitle("ZDC signal (TeV)");
  fOutput->Add(fEZDCvsNtracklets);      
  
  if(fTreeFilling){
    OpenFile(2);
    fCentralityTree = new TTree("fCentralityTree", "Centrality vs. multiplicity tree");
    //
    fCentralityTree->Branch("nev", &fNev,"nev/I");
    fCentralityTree->Branch("beamEnergy", &fBeamEnergy,"beamEnergy/F"); 
    fCentralityTree->Branch("nmyTracks_gen", &fNmyTracks_gen,"nmyTracks_gen/I"); 
    fCentralityTree->Branch("trigClass",&fTrigClass,"trigClass/C");
    fCentralityTree->Branch("xVertex", &fxVertex,"xVertex/D");
    fCentralityTree->Branch("yVertex", &fyVertex,"yVertex/D");
    fCentralityTree->Branch("zVertex", &fzVertex,"zVertex/D");
    fCentralityTree->Branch("vertexer3d", &fVertexer3d,"vertexer3d/O");
    fCentralityTree->Branch("bMC", &fbMC,"bMC/D");
    fCentralityTree->Branch("npartTargMC", &fNpartTargMC,"npartTargMC/I");
    fCentralityTree->Branch("npartProjMC", &fNpartProjMC,"npartProjMC/I");
    fCentralityTree->Branch("NNColl", &fNNColl,"NNColl/I");
    fCentralityTree->Branch("NwNColl", &fNwNColl,"NwNColl/I");
    fCentralityTree->Branch("NNwColl", &fNNwColl,"NNwColl/I");
    fCentralityTree->Branch("NwNwColl", &fNwNwColl,"NwNwColl/I");
    fCentralityTree->Branch("nTracklets", &fNTracklets,"nTracklets/I");
    fCentralityTree->Branch("nSingleClusters", &fNSingleClusters,"nSingleClusters/I");
    fCentralityTree->Branch("nClusters", fNClusters,"nClusters[6]/I");
    fCentralityTree->Branch("nChips", fNChips,"nChips[2]/I");
    fCentralityTree->Branch("bZDC", &fbZDC,"bZDC/D");
    fCentralityTree->Branch("npartZDC", &fNpartZDC,"npartZDC/I");
    fCentralityTree->Branch("bZDCA", &fbZDCA,"bZDCA/D");
    fCentralityTree->Branch("npartZDCA", &fNpartZDCA,"npartZDCA/I");
    fCentralityTree->Branch("bZDCC", &fbZDCC,"bZDCC/D");
    fCentralityTree->Branch("npartZDCC", &fNpartZDCC,"npartZDCC/I");
    fCentralityTree->Branch("esdFlag", &fESDFlag,"esdFlag/i");
    fCentralityTree->Branch("zncEnergy",  &fZNCEnergy,  "zncEnergy/F");
    fCentralityTree->Branch("zpcEnergy",  &fZPCEnergy,  "zpcEnergy/F");
    fCentralityTree->Branch("znaEnergy",  &fZNAEnergy,  "znaEnergy/F");
    fCentralityTree->Branch("zpaEnergy",  &fZPAEnergy,  "zpaEnergy/F");
    fCentralityTree->Branch("zem1Energy", &fZEM1Energy, "zem1Energy/F");
    fCentralityTree->Branch("zem2Energy", &fZEM2Energy, "zem2Energy/F");
    fCentralityTree->Branch("znctower", fZNCtower, "znctower[5]/F");
    fCentralityTree->Branch("zpctower", fZPCtower, "zpctower[5]/F");
    fCentralityTree->Branch("znatower", fZNAtower, "znatower[5]/F");
    fCentralityTree->Branch("zpatower", fZPAtower, "zpatower[5]/F");
    fCentralityTree->Branch("centrZNC", fCentrZNC, "centrZNC[2]/F");
    fCentralityTree->Branch("centrZNA", fCentrZNA, "centrZNA[2]/F");
    fCentralityTree->Branch("nTracks",    &fNTracks,"nTracks/I");
    fCentralityTree->Branch("nPmdTracks", &fNPmdTracks,"nPmdTracks/I");
    fCentralityTree->Branch("multV0A",    &fMultV0A,"multV0A/F");
    fCentralityTree->Branch("multV0C",    &fMultV0C,"multV0C/F");
    fCentralityTree->Branch("multFMDA", &fMultFMDA,"multFMDA/F");
    fCentralityTree->Branch("multFMDC", &fMultFMDC,"multFMDC/F");
  }
}

//________________________________________________________________________
void AliAnalysisTaskCentralityTreeMaker::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskCentralityTreeMaker::UserExec() \n");
  
  if(fAnalysisInput.CompareTo("ESD")==0){

    //    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    
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
	  fNpartTargMC = 208.-(specNeutronTarg+specProtonTarg);
	  fNpartProjMC = 208.-(specNeutronProj+specProtonProj);
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
      //fmd->Print();
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
	

      AliESDZDC *esdZDC = esd->GetESDZDC();
      
      fESDFlag =  esdZDC->GetESDQuality();   
      
      fZNCEnergy = (Float_t) (esdZDC->GetZDCN1Energy());
      fZPCEnergy = (Float_t) (esdZDC->GetZDCP1Energy());
      fZNAEnergy = (Float_t) (esdZDC->GetZDCN2Energy());
      fZPAEnergy = (Float_t) (esdZDC->GetZDCP2Energy());
      fZEM1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0));
      fZEM2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1));
      //
      //printf(" E_ZEM = %1.2f GeV, E_ZDC = %1.2f  TeV\n",fZEM1Energy+fZEM2Energy,  (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fEZDCvsEZEM->Fill(fZEM1Energy+fZEM2Energy, (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      fEZDCvsNtracklets->Fill(fNTracklets, (fZNCEnergy+fZPCEnergy+fZNAEnergy+fZPAEnergy)/1000.);
      PostData(1, fOutput);
      
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
      
      if(fTreeFilling){
        fCentralityTree->Fill();
        PostData(2, fCentralityTree);
      }
  }   
  else if(fAnalysisInput.CompareTo("AOD")==0){
    //AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    // to be implemented
    printf("  AOD analysis not yet implemented!!!\n\n");
    return;
  }
       
}



//________________________________________________________________________
void AliAnalysisTaskCentralityTreeMaker::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
/*  if(fDebug > 1) printf(" **** AliAnalysisTaskCentralityTreeMaker::Terminate() \n");
  
  fCentralityTree = dynamic_cast<TTree*> (GetOutputData(1));
  if(!fCentralityTree) printf("ERROR: fCentralityTree not available\n");

  //fOutput = dynamic_cast<TList*> (GetOutputData(1));
  //if(!fOutput) printf("ERROR: fOutput not available\n");
*/
}


