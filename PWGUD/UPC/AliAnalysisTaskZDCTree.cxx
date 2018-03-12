/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Based on files of the ALICE Off-line Project.                          *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 * 																		  *
 * Based on the class AliAnalysisTaskZDCTreeMaker by Chiara Oppedisano	  *																	  *
 * Created by Uliana Dmitrieva uliana.dmitrieva@cern.ch on 11/01/2018     *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//							   //
//	Class to analyze ZDC data			   //
//							   //
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>

#include "TChain.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskZDCTree.h"

ClassImp(AliAnalysisTaskZDCTree)


//________________________________________________________________________
AliAnalysisTaskZDCTree::AliAnalysisTaskZDCTree():
  AliAnalysisTaskSE(),
    fDebug(0),
    fESD(0),
    fOutput(0x0),
    fZDCTree(0x0),
    fRunNum(0),
    fIsZEM1(kFALSE),
    fIsZEM2(kFALSE),
    fIsZNC(kFALSE),
    fIsZPC(kFALSE),
    fIsZNA(kFALSE),
    fIsZPA(kFALSE),
    fIs1ZED(kFALSE),
    fIsCTRUE(kFALSE),
    fZNCEnergy(0), 
    fZPCEnergy(0),  
    fZNAEnergy(0),  
    fZPAEnergy(0),
    fZEM1Energy(0), 
    fZEM2Energy(0)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}   

//________________________________________________________________________
AliAnalysisTaskZDCTree::AliAnalysisTaskZDCTree(const char *name):
  AliAnalysisTaskSE(name),
    fDebug(0),
    fESD(0),
    fOutput(0x0),
    fZDCTree(0x0),
    fRunNum(0),
    fIsZEM1(kFALSE),
    fIsZEM2(kFALSE),
    fIsZNC(kFALSE),
    fIsZPC(kFALSE),
    fIsZNA(kFALSE),
    fIsZPA(kFALSE),
    fIs1ZED(kFALSE),
    fIsCTRUE(kFALSE), 
    fZNCEnergy(0), 
    fZPCEnergy(0),  
    fZNAEnergy(0),  
    fZPAEnergy(0),
    fZEM1Energy(0), 
    fZEM2Energy(0)


{
  // Default constructor
   
  for(Int_t itow=0; itow<5; itow++)		fZNCtower[itow]= fZPCtower[itow]= fZNAtower[itow]= fZPAtower[itow]= 0.;  
     
  for(Int_t ihit=0; ihit<4; ihit++)		fZNCTDC[ihit]= fZPCTDC[ihit]= fZNATDC[ihit]= fZNATDC[ihit]= 9999.; 
    
  // constructor
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TList::Class()); 
}
 
//________________________________________________________________________
AliAnalysisTaskZDCTree::~AliAnalysisTaskZDCTree()
{
  // Destructor
  if(fOutput){   delete fOutput;  } 

}
//________________________________________________________________________
void AliAnalysisTaskZDCTree::UserCreateOutputObjects()
{
  // Create the output containers
  if(fDebug>1) printf("AliAnalysisTaskZDCTree::UserCreateOutputObjects() \n");

	// Several histograms are more conveniently managed in a TList
	fOutput = new TList;
	fOutput->SetOwner(kTRUE);
 
	//define tree
    fZDCTree = new TTree("fZDCTree", "ZDC tree");
    
    fZDCTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
    
    fZDCTree->Branch("isZEM1",&fIsZEM1,"isZEM1/O");
    fZDCTree->Branch("isZEM2",&fIsZEM2,"isZEM2/O");
    
    fZDCTree->Branch("isZNC",&fIsZNC,"isZNC/O");
    fZDCTree->Branch("isZPC",&fIsZPC,"isZPC/O");
    fZDCTree->Branch("isZNA",&fIsZNA,"isZNA/O");
    fZDCTree->Branch("isZPA",&fIsZPA,"isZPA/O");
     
    fZDCTree->Branch("is1ZED",&fIs1ZED,"is1ZED/O");
    fZDCTree->Branch("isCTRUE",&fIsCTRUE,"isCTRUE/O");
    
    fZDCTree->Branch("zncEnergy",  &fZNCEnergy,  "zncEnergy/F");
    fZDCTree->Branch("zpcEnergy",  &fZPCEnergy,  "zpcEnergy/F");
    fZDCTree->Branch("znaEnergy",  &fZNAEnergy,  "znaEnergy/F");
    fZDCTree->Branch("zpaEnergy",  &fZPAEnergy,  "zpaEnergy/F");
    fZDCTree->Branch("zem1Energy", &fZEM1Energy, "zem1Energy/F");
    fZDCTree->Branch("zem2Energy", &fZEM2Energy, "zem2Energy/F");

    fZDCTree->Branch("znctower", fZNCtower, "znctower[5]/F");
    fZDCTree->Branch("zpctower", fZPCtower, "zpctower[5]/F");
    fZDCTree->Branch("znatower", fZNAtower, "znatower[5]/F");
    fZDCTree->Branch("zpatower", fZPAtower, "zpatower[5]/F");
 
    fZDCTree->Branch("zncTDC", fZNCTDC, "zncTDC[4]/F");
    fZDCTree->Branch("zpcTDC", fZPCTDC, "zpcTDC[4]/F");
    fZDCTree->Branch("znaTDC", fZNATDC, "znaTDC[4]/F");
    fZDCTree->Branch("zpaTDC", fZPATDC, "zpaTDC[4]/F");
  
   fOutput->Add(fZDCTree);   
    PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskZDCTree::UserExec(Option_t */*option*/)
{
	// Execute analysis for current event:
	if(fDebug>1) printf(" **** AliAnalysisTaskZDCTree::UserExec() \n");
  
      
      fESD = dynamic_cast<AliESDEvent*> (InputEvent());
      if(!fESD) return;
      
      //event info
      fRunNum = fESD->GetRunNumber();
            
      // ***** Trigger selection
      TString triggerClass = fESD->GetFiredTriggerClasses();
       
       Bool_t isZED = kFALSE;
       Bool_t isCTRUE = kFALSE;
     
      // Taking  C1ZED triggers
      if(triggerClass.Contains("C1ZED-B")) isZED = kTRUE;
      fIs1ZED = Bool_t (isZED);

      // Taking  CTRUE triggers
      if(triggerClass.Contains("CTRUE-B")) isCTRUE = kTRUE;
      fIsCTRUE = Bool_t (isCTRUE);

      if (!(isCTRUE || isZED)) return;
        
      AliESDZDC *esdZDC = fESD->GetESDZDC();
      
      fIsZEM1 = Bool_t (esdZDC->IsZEM1hit());
      fIsZEM2 = Bool_t (esdZDC->IsZEM2hit());
      fIsZNC = Bool_t (esdZDC->IsZNChit());
      fIsZPC = Bool_t (esdZDC->IsZPChit());
      fIsZNA = Bool_t (esdZDC->IsZNAhit());
      fIsZPA = Bool_t (esdZDC->IsZPAhit());
      
      fZNCEnergy = (Float_t) (esdZDC->GetZDCN1Energy());
      fZPCEnergy = (Float_t) (esdZDC->GetZDCP1Energy());
      fZNAEnergy = (Float_t) (esdZDC->GetZDCN2Energy());
      fZPAEnergy = (Float_t) (esdZDC->GetZDCP2Energy());
      fZEM1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0));
      fZEM2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1));
       
      const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
      const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
      const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
      const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
      //
      
      //
      for(Int_t it=0; it<5; it++){
         fZNCtower[it] = (Float_t) (towZNC[it]);
         fZPCtower[it] = (Float_t) (towZPC[it]);
         fZNAtower[it] = (Float_t) (towZNA[it]); 
         fZPAtower[it] = (Float_t) (towZPA[it]);  
         }
      
      //   
      for(int itdc=0; itdc<4; itdc++){
     fZNCTDC[itdc] = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNCTDCChannel(), itdc);
	 fZPCTDC[itdc] = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNCTDCChannel(), itdc);
	 fZNATDC[itdc] = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNCTDCChannel(), itdc);
	 fZPATDC[itdc] = esdZDC->GetZDCTDCCorrected(esdZDC->GetZNCTDCChannel(), itdc);
      }

      
         
  fZDCTree->Fill();

  PostData(1, fOutput);

   
}



//________________________________________________________________________
void AliAnalysisTaskZDCTree::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  printf(" **** AliAnalysisTaskZDCTree::Terminate() \n");
}

