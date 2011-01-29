/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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

/* $Id$ sjena*/

//-------------------------------------------------------------------------
//                 AliEbyEFluctuationAnalysisTask
//   This is the class to deal with the EbyE Charge Fluctuation  analysis
//               Origin: Satyajit Jena, sjena@cern.ch
//-------------------------------------------------------------------------

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVZERO.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliVertexerTracks.h"
#include "AliCentrality.h"
#include "AliEbyEFluctuationAnalysisTask.h"

ClassImp(AliEbyEFluctuationAnalysisTask)

//________________________________________________________________________
AliEbyEFluctuationAnalysisTask::AliEbyEFluctuationAnalysisTask() 
  :AliAnalysisTaskSE(),
   fTrackCut(0),
   fListQA(0),
   fListResults(0),
   fAnalysisType("ESD"),
   fAnalysisMode("TPC"),
   fCentType("M"),
   fCentBin(10),
   fCentMode("F"),
   fCentEstimator("V0M"),
   fVx(5.),
   fVy(5.),
   fVz(10.),
   fHistEventStats(0),
   fhVertex(0),
   fhVzeroMult(0),
   fCentVsMult(0),
   fhVnT(0) {
  // Dummy constructor ALWAYS needed for I/O.
  for(Int_t i = 0; i < 20; i++) fHCentrality[i] = NULL;
}


//________________________________________________________________________
AliEbyEFluctuationAnalysisTask::AliEbyEFluctuationAnalysisTask(const char *name) 
  :AliAnalysisTaskSE(name),
   fTrackCut(0),
   fListQA(0),
   fListResults(0),
   fAnalysisType("ESD"),
   fAnalysisMode("TPC"),
   fCentType("M"),
   fCentBin(5),
   fCentMode("F"),
   fCentEstimator("V0M"),
   fVx(5.),
   fVy(5.),
   fVz(10.),
   fHistEventStats(0),
   fhVertex(0),
   fhVzeroMult(0),
   fCentVsMult(0),
   fhVnT(0) {
  //Constructor
  for(Int_t i = 0; i < 20; i++) fHCentrality[i] = NULL;
  
  DefineOutput(1, TList::Class());                                         
  DefineOutput(2, TList::Class());                                         
}

//________________________________________________________________________
AliEbyEFluctuationAnalysisTask::~AliEbyEFluctuationAnalysisTask() {
  //Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if(fListQA) delete fListQA;
    if(fTrackCut) delete fTrackCut;
  }
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserCreateOutputObjects() {
  //Output objects
  fListQA = new TList();
  fListQA->SetName("fListQA");
  
  //InitTrackCuts();
  
  Int_t lbin = 100/fCentBin + 2; 
  
  TString gCutName[4] = {"Total","Offline trigger","Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats","Event statistics;;N_{events}",4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  
  fListQA->Add(fHistEventStats);
  
  fhVertex = new TH2F("hVertex","Vertexer;Bin For Different Vertexer;Value",7,0,7,1000,-25,25);
  fListQA->Add(fhVertex);
  
  fhVzeroMult = new TH2F("hVzMultilicity","",lbin,0,(Double_t)lbin,1000,0,25000);
  fListQA->Add(fhVzeroMult);
  
  fhVnT = new TH2F("hVzAndTrack","",lbin,0,(Double_t)lbin,1000,0,25000);
  fListQA->Add(fhVnT);
  
  fCentVsMult = new TH2F("hCentVsMult","Centrality vs Multiplicity[N_{+} + N_{-}];Centrality; Multiplicity[N_{+} + N_{-}]",lbin,0,(Double_t)lbin,4000,0,25000);
  
  fListQA->Add(fCentVsMult);
  
  //===========================================================//
  fListResults = new TList();
  fListResults->SetName("fListResults");
  TString histName;
  for(Int_t iBin = 0; iBin < 20 ; iBin++) {
    histName = "fHistPN"; histName += "Centrality"; histName += iBin;
    fHCentrality[iBin] = new TH2F(histName.Data(),
				  ";N_{+};N_{-}",
				  2000,0.5,2000.5,2000,0.5,2000.5);
    fListResults->Add(fHCentrality[iBin]);
  }
  
  PostData(1, fListQA); 
  PostData(2, fListResults); 
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserExec(Option_t *) {
  //Analysis part
  fHistEventStats->Fill(1); //all events 
  TString gAnalysisLevel = "ESD";
  
  //AliCentrality *centrality = esd->GetCentrality();  
  Double_t mult = 0;
  Int_t pch  = 0;
  Int_t nch  = 0;
  Int_t cent = 0;  
  
  //ESD analysis
  if(fAnalysisType.CompareTo("ESD") == 0) {
    AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); // from Task
    if (!gESD) {
      Printf("ERROR: gESD not available");
      return;
    }
       
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(!isSelected) return; 
    fHistEventStats->Fill(2); //physics selection
 
    AliESDVZERO* esdV0 = gESD->GetVZEROData();
    mult = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();

    fHistEventStats->Fill(3); //offline cut

    const AliESDVertex *vertex = GetVertex(gESD);
    if(!vertex) return;
    fHistEventStats->Fill(3); //vertexcut

    AliCentrality *centrality = gESD->GetCentrality();
    cent = centrality->GetCentralityClass5(fCentEstimator.Data());

    // cent = GetCentrality(gESD);
    if(cent < 0 || cent > 19 ) return;
    fHistEventStats->Fill(4); //Analysed Events
    
    //Track loop
    for (Int_t iTracks = 0; iTracks < gESD->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = gESD->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive track %d", iTracks);
	continue;
      }
      
      if(fTrackCut) {
	if(fAnalysisMode = "TPC") {
	  AliESDtrack *tpcOnlyTrack = fTrackCut->GetTPCOnlyTrack(gESD,iTracks);
	  if(!tpcOnlyTrack) continue;
	  if(!fTrackCut->AcceptTrack(tpcOnlyTrack)) continue;
	}
	else 
	  if(!fTrackCut->AcceptTrack(track)) continue;
      }

      //----------------
      Short_t charge = track->Charge();
      if(charge > 0) pch += 1; 
      if(charge < 0) nch += 1; 
      //---------------
    }//track loop
  }//ESD analysis

  fCentVsMult->Fill(cent,pch+nch);
  fHCentrality[cent]->Fill(pch,nch);
  
  fhVnT->Fill(cent,mult);
  fhVnT->Fill(cent,pch+nch);
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::Terminate(Option_t *) {
  //Terminate
  fListQA = dynamic_cast<TList*> (GetOutputData(1));
  if(!fListQA) { Printf("ERROR: could not retrieve TList fListQA"); return; }
}

//____________________________________________________________________//
Int_t AliEbyEFluctuationAnalysisTask::GetCentrality(AliESDEvent *esd) const {
  //Return the centrality
  Int_t cent = -1;
  
  Int_t lbin = 100/fCentBin + 1;

  AliESDVZERO* esdV0 = esd->GetVZEROData();
  Double_t mult = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
  fhVzeroMult->Fill(lbin,mult);
  
  if(fCentType.CompareTo("M") == 0) {
   
    if(fAnalysisType.CompareTo("ESD") == 0) 
      cent = FindCentralityESD(mult);
    if(fAnalysisType.CompareTo("MC") == 0) 
      cent = FindCentralityMC(mult);
    if(cent > -1)  fhVzeroMult->Fill(cent,mult);
    return cent;
  }
  

 if(fCentType.CompareTo("G") == 0) {
   //   printf("========= Inside Official Centrality Selection ======\n");
   AliCentrality *centrality = esd->GetCentrality();
   if(fCentBin == 10)
     cent = centrality->GetCentralityClass10(fCentEstimator);
   else
     cent = centrality->GetCentralityClass5(fCentEstimator);

   return cent;
 }

 return cent;
}

//____________________________________________________________________//
Int_t AliEbyEFluctuationAnalysisTask::FindCentralityESD(Double_t mult) const {
  //hardcoded centrality bins (to be removed)
  Double_t data[101] = { 250000.0,
			 17760.0,  16870.0,  16070.0,  15360.0,  14660.0, 
			 13970.0,  13320.0,  12700.0,  12140.0,  11570.0, 
			 11020.0,  10530.0,  10030.0,   9560.0,   9110.0, 
			 8680.0,   8270.0,   7870.0,   7490.0,   7120.0, 
			 6780.0,   6440.0,   6120.0,   5810.0,   5510.0, 
			 5230.0,   4950.0,   4690.0,   4440.0,   4190.0, 
			 3960.0,   3740.0,   3530.0,   3330.0,   3130.0, 
			 2940.0,   2770.0,   2600.0,   2440.0,   2290.0, 
			 2150.0,   2010.0,   1880.0,   1760.0,   1640.0, 
			 1540.0,   1440.0,   1350.0,   1260.0,   1180.0, 
			 1110.0,   1030.0,    970.0,    900.0,    840.0, 
			 790.0,    740.0,    700.0,    650.0,    610.0, 
			 580.0,    540.0,    510.0,    480.0,    450.0, 
			 430.0,    400.0,    380.0,    360.0,    340.0, 
			 320.0,    300.0,    290.0,    270.0,    250.0, 
			 240.0,    230.0,    210.0,    200.0,    190.0, 
			 180.0,    170.0,    160.0,    150.0,    150.0, 
			 140.0,    130.0,    120.0,    120.0,    110.0, 
			 100.0,    100.0,     90.0,     90.0,     80.0, 
			 80.0,     80.0,     70.0,     70.0,     70.0 };
  

 if(fCentMode.CompareTo("F") == 0) { 
   for(Int_t i = 0; i < 100; i+=fCentBin) {
      
     if( (mult < data[i]) && (mult >= data[i+fCentBin]) ) 
       return i/fCentBin;
      }
  }
  if(fCentMode.CompareTo("A") == 0) {
    for(int i = 0; i < 100; i+=fCentBin) {
      if((mult < data[0]) && (mult >data[i+fCentBin]) ) {
	return i/fCentBin; 
      }
    }
  } 
  
return -2;
  
  
}

//____________________________________________________________________//
Int_t AliEbyEFluctuationAnalysisTask::FindCentralityMC(Double_t mult) const {
  //hardcoded centrality bins (MC) ==> to be removed
  Double_t data[101] = { 250000.0, 17760.0,  16870.0,  16070.0,  15360.0,  
			 14660.0, 13970.0,  13320.0,  12700.0,  12140.0,  
			 11570.0, 11020.0,  10530.0,  10030.0,   9560.0,   9110.0, 
			 8680.0,   8270.0,   7870.0,   7490.0,   7120.0, 
			 6780.0,   6440.0,   6120.0,   5810.0,   5510.0, 
			 5230.0,   4950.0,   4690.0,   4440.0,   4190.0, 
			 3960.0,   3740.0,   3530.0,   3330.0,   3130.0, 
			 2940.0,   2770.0,   2600.0,   2440.0,   2290.0, 
			 2150.0,   2010.0,   1880.0,   1760.0,   1640.0, 
			 1540.0,   1440.0,   1350.0,   1260.0,   1180.0, 
			 1110.0,   1030.0,    970.0,    900.0,    840.0, 
			 790.0,    740.0,    700.0,    650.0,    610.0, 
			 580.0,    540.0,    510.0,    480.0,    450.0, 
			 430.0,    400.0,    380.0,    360.0,    340.0, 
			 320.0,    300.0,    290.0,    270.0,    250.0, 
			 240.0,    230.0,    210.0,    200.0,    190.0, 
			 180.0,    170.0,    160.0,    150.0,    150.0, 
			 140.0,    130.0,    120.0,    120.0,    110.0, 
			 100.0,    100.0,     90.0,     90.0,     80.0, 
			 80.0,     80.0,     70.0,     70.0,     70.0 
  };

  if(fCentMode.CompareTo("F") == 0) { 
    for(Int_t i = 0; i < 100; i+=fCentBin) {
      if( (mult < data[i]) && (mult >= data[i+fCentBin]) ) 
	return i/fCentBin;
      }
  }
  if(fCentMode.CompareTo("A") == 0) {
    for(int i = 0; i < 100; i+=fCentBin) {
      if((mult < data[0]) && (mult >= data[i+fCentBin]) ) {
	return i/fCentBin; 
      }
    }
  } 
  
  return -3;
}

//____________________________________________________________________//
const AliESDVertex* AliEbyEFluctuationAnalysisTask::GetVertex(AliESDEvent* esd) {
  //Return the primary vertex
  const AliESDVertex* vertex = 0;

  if(fAnalysisMode = "TPC") 
    vertex = esd->GetPrimaryVertexTPC();
  if(fAnalysisMode = "Hybrid") 
    vertex = esd->GetPrimaryVertexSPD();
  if(fAnalysisMode = "Global") 
    vertex = esd->GetPrimaryVertex();

  if(vertex) {
    if(vertex->GetNContributors() > 0) {
      if(vertex->GetZRes() != 0) {
	fhVertex->Fill(1,vertex->GetXv());
	fhVertex->Fill(2,vertex->GetYv());
	fhVertex->Fill(3,vertex->GetZv());
	
	if(TMath::Abs(vertex->GetXv()) < fVx) {
	  if(TMath::Abs(vertex->GetYv()) < fVy) {
	    if(TMath::Abs(vertex->GetZv()) < fVz) {
	      
	      fhVertex->Fill(4,vertex->GetXv());
   	      fhVertex->Fill(5,vertex->GetYv());
	      fhVertex->Fill(6,vertex->GetZv());
	      return vertex;	
	      
	    }
	  }
	}
      }
    }
  }
  
  return vertex;
}
