// This class is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.5
// Version 3.5.5


#include "Riostream.h"
#include "AliLRCProcess.h"
#include "iostream"
ClassImp(AliLRCProcess)

AliLRCProcess::AliLRCProcess(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA): fIsEventOpend(kFALSE), fIsOnline(kFALSE), fDisplayInitOnDemandWarning(kTRUE), fEventCount(0),fStartForwardETA(0), fEndForwardETA(0), fStartBakwardETA(0), fEndBakwardETA(0),fHiPt(0),fLoPt(0),fHiMult(0),fLoMult(0),fMultBins(0),fPtBins(0),fSumPtFw(0),  fSumPtBw(0),fNchFw(0), fNchBw(0),fOutList(0), fShortDef(0), fHistPt(0),fHistEta(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fProfNberr(0),fProfdPtB(0),fProfTestLRC(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistPtBakward(0),fHistEtaBakward(0),fHistNchBakward(0)
{
// Constructor with window setup makes ready-to-run processor
fEventCount=0;

SetETAWindows( _StartForwardETA, _EndForwardETA,_StartBakwardETA,_EndBakwardETA);
SetPtRange(0.0,4.0,160);
SetMultRange(0,100);     
}

Bool_t AliLRCProcess::InitDataMembers()
{
// This method is actualy creating output histogramms
// Thist method  is to be called in CreateOutputObjects method of AliAnalysisTask
cout<<" # Init for "<<fShortDef<<" this="<<this<<"\n";
if(fIsOnline)
{
Printf("Can't init data members more then one time! \n");
return kFALSE;
}
fEventCount=0;
fOutList = new TList();
fOutList->SetName(fShortDef);

   Double_t LoMult,HiMult;
   
   LoMult=fLoMult-0.5;
   HiMult=fHiMult+0.5;

   // Window statistics histograms
  
   // Forward
  
  fHistPtForward = new TH1D("fHistPtForward", "P_{T} distribution in Forward window", 100, 0.0, 5);
  fHistPtForward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtForward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtForward->SetMarkerStyle(kFullCircle);
  
  
  fHistEtaForward = new TH1D("fEtaForward", "#eta distribution in Forward window", 200, -2, 2);
  fHistEtaForward->GetXaxis()->SetTitle("ETA");
  fHistEtaForward->GetYaxis()->SetTitle("dN/ETA");
  fHistEtaForward->SetMarkerStyle(kFullCircle);
  
  
  fHistNchForward = new TH1D("fHistNchForward", "N_{ch} distribution in Forward window", fMultBins, LoMult, HiMult);
  fHistNchForward->GetXaxis()->SetTitle("N_{ch}");
  fHistNchForward->GetYaxis()->SetTitle("dN/dN_{ch}");
  fHistNchForward->SetMarkerStyle(kFullCircle);
 
     // Bakward
  
  fHistPtBakward = new TH1D("fHistPtBakward", "P_{T} distribution in Bakward window", 100, 0.0, 5);
  fHistPtBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistEtaBakward = new TH1D("fEtaBakward", "#eta distribution in Bakward window", 200, -2, 2);
  fHistEtaBakward->GetXaxis()->SetTitle("ETA");
  fHistEtaBakward->GetYaxis()->SetTitle("dN/ETA");
  fHistEtaBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistNchBakward = new TH1D("fHistNchBakward", "N_{ch} distribution in Bakward window", fMultBins, LoMult, HiMult);
  fHistNchBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistNchBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistNchBakward->SetMarkerStyle(kFullCircle);


    //Overal statistics histograms
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 100, 0.0, 5.0);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  
  fHistEta = new TH1F("fHistEta", "#eta distribution", 200, -2, 2);
  fHistEta->GetXaxis()->SetTitle("ETA");
  fHistEta->GetYaxis()->SetTitle("dN/ETA");
  fHistEta->SetMarkerStyle(kFullCircle);
    
    
    
    // -------- LRC histograms
  
  fHistNN = new TH2D("NN","NN",fMultBins, LoMult, HiMult,fMultBins, LoMult, HiMult);
  fHistPtN = new TH2D("PtN","PtN",fMultBins, LoMult, HiMult,fPtBins,fLoPt,fHiPt);
  fHistPtPt = new TH2D("PtPt","PtPt",fPtBins,fLoPt,fHiPt,fPtBins,fLoPt,fHiPt);
  fProfNberr = new TProfile("nber","nber",fMultBins, LoMult, HiMult);
  fProfdPtB = new TProfile("dPtB","Overal multievent Pt_Backward (first bin) Pt_Backward^2 (sec. bin) ",2,0.5,2.5);  
  fProfTestLRC = new TProfile("TestLRC","Test LRC calculaion via TProfile",fMultBins, LoMult, HiMult);  


   // ---------- Adding data members to output list
  
  // Adding overal statistics
  
  fOutList->Add(fHistPt);
  fOutList->Add(fHistEta);
  
  //Adding LRC hists
  
  fOutList->Add(fHistNN);
  fOutList->Add(fHistPtN);
  fOutList->Add(fHistPtPt);
  fOutList->Add(fProfNberr);
  fOutList->Add(fProfdPtB);
  fOutList->Add(fProfTestLRC);

  
  //Adding window statistics

  fOutList->Add(fHistPtForward);
  fOutList->Add(fHistEtaForward);
  fOutList->Add(fHistNchForward);
  fOutList->Add(fHistPtBakward);
  fOutList->Add(fHistEtaBakward);
  fOutList->Add(fHistNchBakward);
  
 


fIsOnline=kTRUE;
return kTRUE;
}
AliLRCProcess::~AliLRCProcess()
{
//Destructor

}

// ---------------------------------------  Setters ------------------
  void AliLRCProcess::SetShortDef()
  {
  // Creating task and output container name 
   TString lF,hF,lB,hB;
   lF+=fStartForwardETA; lF.Remove(TString::kBoth,' ');
   hF+=fEndForwardETA; hF.Remove(TString::kBoth,' ');
   lB+=fStartBakwardETA; lB.Remove(TString::kBoth,' ');
   hB+=fEndBakwardETA; hB.Remove(TString::kBoth,' ');
   
   fShortDef= "TaskLRCw"+lF+"to"+hF+"vs"+lB+"to"+hB;

  }

  void AliLRCProcess::SetForwardWindow(Double_t StartETA,Double_t EndETA)
  {
    //setter for the forward eta window
  fStartForwardETA=StartETA;
  fEndForwardETA=EndETA;
  SetShortDef();
  }
  void AliLRCProcess::SetBackwardWindow(Double_t StartETA,Double_t EndETA)
  {
    //setter for the backward eta window
  fStartBakwardETA=StartETA;
  fEndBakwardETA=EndETA;
  SetShortDef();
  }
  void AliLRCProcess::SetETAWindows(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA)
  {
  //setter for the eta windows
  fStartForwardETA=_StartForwardETA;
  fEndForwardETA=_EndForwardETA;
  fStartBakwardETA=_StartBakwardETA;
  fEndBakwardETA=_EndBakwardETA;
  SetShortDef();
  }

  void AliLRCProcess::SetPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins)
  {
  // Setter for Pt range and N bins in histos
  fLoPt=LoPt;
  fHiPt=HiPt;
  fPtBins=PtBins;  
  }
  void AliLRCProcess::SetMultRange(Int_t LoMult,Int_t HiMult,Int_t MultBins)
  {
  // Setter for multiplicity range and N bins in histos
  fLoMult=LoMult;
  fHiMult=HiMult;
  if(!MultBins)
  	{
  	fMultBins=fHiMult-fLoMult+1;
  	}else
  	{
  	fMultBins=MultBins;  
  	}
  }


//________________________________________________________________________



TList* AliLRCProcess::CreateOutput()
{
// Creates a link to output data TList
return fOutList;
}

TString AliLRCProcess::GetShortDef()
{
return fShortDef;
}

void AliLRCProcess::StartEvent()
{
// Open new Event for track by track event import
if(fIsEventOpend)
	{Printf("Event is already opened! Auto finishing ! \n");
	 cout<<fShortDef<<": event count = "<<fEventCount<<" ";
	 Printf("NchF = %i,NchB = %i \n",fNchFw,fNchBw);
	 
	 FinishEvent();
	}
	
fNchFw=0;
fSumPtFw=0;
fNchBw=0;
fSumPtBw=0;

fIsEventOpend=kTRUE;
}
void AliLRCProcess::AddTrackPtEta(Double_t Pt, Double_t Eta )
{
//Track by track event import :  Imports track to the event 

if(!fIsEventOpend)
	{Printf("Event is not opened!\n");
	return;}

    Double_t lPt=Pt;
    Double_t lEta=Eta;
    
    //  Glabal trak data
    fHistPt->Fill(lPt);
    fHistEta->Fill(lEta);
    
       
    //Forward window
    if( (fStartForwardETA<lEta)&&(lEta<fEndForwardETA))
    	{
    	fNchFw++;
    	fSumPtFw+=lPt;
  	fHistPtForward->Fill(lPt);
  	fHistEtaForward->Fill(lEta);
	
    	}
    
    //Backward window
    if((fStartBakwardETA<lEta)&&(lEta<fEndBakwardETA))
    	{
	fNchBw++;
	fSumPtBw+=lPt; 
	if(lPt<4.0)
	{
	fProfdPtB->Fill(1,lPt); 
	fProfdPtB->Fill(2,lPt*lPt);
	}
	 
  	fHistPtBakward->Fill(lPt);
  	fHistEtaBakward->Fill(lEta);
 
    	}
    
	

}
void AliLRCProcess::FinishEvent()
{
// Track by track event import : Close opened event and fill event summary histos

if(!fIsEventOpend)
	{Printf("Event is not opened!\n");
	return;}

 //Filling even-total data
  
  fHistNN->Fill(fNchFw,fNchBw);
  
  if(fNchBw!=0)
  {
	fSumPtBw=fSumPtBw/fNchBw;
	fProfTestLRC->Fill(fNchFw,fSumPtBw);
	fHistPtN->Fill(fNchFw,fSumPtBw);
	fProfNberr->Fill(fNchFw,1.0/fNchBw);
	
	if(fNchFw!=0)
	{
		fSumPtFw=fSumPtFw/fNchFw;
		fHistPtPt->Fill(fSumPtFw,fSumPtBw);
	}
  }
	
  
 fHistNchForward->Fill(fNchFw);
 fHistNchBakward->Fill(fNchBw);

fEventCount++;
fIsEventOpend=kFALSE;
//cout<<fShortDef<<": event count = "<<fEventCount<<" ";
//	 Printf("NchF = %i,NchB = %i",fNchFw,fNchBw);
}

