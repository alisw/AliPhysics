// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This task is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data.  

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.0
// Version: 3.0.8  may 08


#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskLRC.h"



ClassImp(AliAnalysisTaskLRC)

//________________________________________________________________________
AliAnalysisTaskLRC::AliAnalysisTaskLRC(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fHistPt(0),fHistEta(0),fOutList(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fHistNberr(0),fProfdPtB(0),fProfTestLRC(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistPtBakward(0),fHistEtaBakward(0),fHistNchBakward(0)
{
  //Init
  
  // Constructor

  SetETAWindows(-1.0,1.0,-1.0,1.0);   //Default windows (full range)
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container for all histogramms
  DefineOutput(0, TList::Class());
  
}


// ---------------------------------------  Setters ------------------

  void AliAnalysisTaskLRC::SetForwardWindow(double StartETA,double EndETA)
  {
  fStartForwardETA=StartETA;
  fEndForwardETA=EndETA;
  }
  void AliAnalysisTaskLRC::SetBackwardWindow(double StartETA,double EndETA)
  {
  fStartBakwardETA=StartETA;
  fEndBakwardETA=EndETA;
  }
  void AliAnalysisTaskLRC::SetETAWindows(double _StartForwardETA,double _EndForwardETA,double _StartBakwardETA,double _EndBakwardETA)
  {
  fStartForwardETA=_StartForwardETA;
  fEndForwardETA=_EndForwardETA;
  fStartBakwardETA=_StartBakwardETA;
  fEndBakwardETA=_EndBakwardETA;
  }

 


//________________________________________________________________________
void AliAnalysisTaskLRC::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
 
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskLRC::CreateOutputObjects()
{
  // Create histograms
  

   // Window statistics histograms
  
   // Forward
  
  fHistPtForward = new TH1D("fHistPtForward", "P_{T} distribution in Forward window", 150, 0.1, 3.1);
  fHistPtForward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtForward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtForward->SetMarkerStyle(kFullCircle);
  
  
  fHistEtaForward = new TH1D("fEtaForward", "#eta distribution in Forward window", 200, -2, 2);
  fHistEtaForward->GetXaxis()->SetTitle("ETA");
  fHistEtaForward->GetYaxis()->SetTitle("dN/ETA");
  fHistEtaForward->SetMarkerStyle(kFullCircle);
  
  
  fHistNchForward = new TH1D("fHistNchForward", "N_{ch} distribution in Forward window", 201, -0.5, 200.5);
  fHistNchForward->GetXaxis()->SetTitle("N_{ch}");
  fHistNchForward->GetYaxis()->SetTitle("dN/dN_{ch}");
  fHistNchForward->SetMarkerStyle(kFullCircle);
 
     // Bakward
  
  fHistPtBakward = new TH1D("fHistPtBakward", "P_{T} distribution in Bakward window", 150, 0.1, 3.1);
  fHistPtBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistEtaBakward = new TH1D("fEtaBakward", "#eta distribution in Bakward window", 200, -2, 2);
  fHistEtaBakward->GetXaxis()->SetTitle("ETA");
  fHistEtaBakward->GetYaxis()->SetTitle("dN/ETA");
  fHistEtaBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistNchBakward = new TH1D("fHistNchBakward", "N_{ch} distribution in Bakward window", 201, -0.5, 200.5);
  fHistNchBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistNchBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistNchBakward->SetMarkerStyle(kFullCircle);

  
    // --------- Output list
  
  fOutList = new TList();
 
   
  //Overal statistics histograms
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 150, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  
  fHistEta = new TH1F("fEta", "#eta distribution", 200, -2, 2);
  fHistEta->GetXaxis()->SetTitle("ETA");
  fHistEta->GetYaxis()->SetTitle("dN/ETA");
  fHistEta->SetMarkerStyle(kFullCircle);
    
    
    
    // -------- LRC histograms
  
  fHistNN = new TH2D("NN","NN",100,0.5,100.5,100,0.5,100.5);//1000 ������ � ���� �� �������� 
  fHistPtN = new TH2D("PtN","PtN",100,0.5,100.5,40,0,4);//4 gev - ����������� �������
  fHistPtPt = new TH2D("PtPt","PtPt",40,0,4,40,0,4);//4 gev - ����������� �������
  fHistNberr = new TH2D("nber","nber",100,0.5,100.5,1000,0,1);//��� ������� ������ ��� x ��������� � Nf
  fProfdPtB = new TProfile("dPtB","Overal multievent Pt_Backward (first bin) Pt_Backward^2 (sec. bin) ",2,0.5,2.5);  
  fProfTestLRC = new TProfile("TestLRC","Test LRC calculaion via TProfile",100,0.5,100.5);  

  
  // ---------- Adding data members to output list
  
  // Adding overal statistics
  
  fOutList->Add(fHistPt);
  fOutList->Add(fHistEta);
  
  //Adding LRC hists
  
  fOutList->Add(fHistNN);
  fOutList->Add(fHistPtN);
  fOutList->Add(fHistPtPt);
  fOutList->Add(fHistNberr);
  fOutList->Add(fProfdPtB);
  fOutList->Add(fProfTestLRC);

  
  //Adding window statistics

  fOutList->Add(fHistPtForward);
  fOutList->Add(fHistEtaForward);
  fOutList->Add(fHistNchForward);
  fOutList->Add(fHistPtBakward);
  fOutList->Add(fHistEtaBakward);
  fOutList->Add(fHistNchBakward);
  
  
  
}

//________________________________________________________________________
void AliAnalysisTaskLRC::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  

   //Event variables
  double lPtF=0;   //Forward window Pt
  double lPtB=0;   //Bakward window Pt
  
  int lNF=0;  // Forward multiplicity
  int lNB=0;  // Bakward multiplicity
  
  
  //Track variables
  double lPt;   // Temp Pt
  double lEta;	  // Temp ETA	
  
    // Track loop 
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    
    lPt=track->Pt();
    lEta=track->Eta();
    
    //  Glabal trak data
    fHistPt->Fill(lPt);
    fHistEta->Fill(lEta);
    
    //Forward window
    if( (fStartForwardETA<lEta)&&(lEta<fEndForwardETA))
    	{
    	lNF++;
    	lPtF+=lPt;
  	fHistPtForward->Fill(lPt);
  	fHistEtaForward->Fill(lEta);
	
    	}
    
    //Backward window
    if((fStartBakwardETA<lEta)&&(lEta<fEndBakwardETA))
    	{
	lNB++;
	lPtB+=lPt; 
	if(lPt<4.0)
	{
	fProfdPtB->Fill(1,lPt); 
	fProfdPtB->Fill(2,lPt*lPt);
	}
	 
  	fHistPtBakward->Fill(lPt);
  	fHistEtaBakward->Fill(lEta);
 
    	}
    
    
    
    
  } //end of track loop 
  
  
  //Filling even-total data
  fHistNN->Fill(lNF,lNB);
  
  if(lNB!=0)
  {
	lPtB=lPtB/lNB;
	fProfTestLRC->Fill(lNF,lPtB);
	fHistPtN->Fill(lNF,lPtB);
	fHistNberr->Fill(lNF,1.0/lNB);
	
	if(lNF!=0)
	{
		lPtF=lPtF/lNF;
		fHistPtPt->Fill(lPtF,lPtB);
	}
  }
	
  
 fHistNchForward->Fill(lNF);
 fHistNchBakward->Fill(lNB);
  // Post output data.
  
  PostData(0, fOutList);
}      

//________________________________________________________________________
void AliAnalysisTaskLRC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
fOutList = dynamic_cast<TList*> (GetOutputData(0));

  
 
  
}
