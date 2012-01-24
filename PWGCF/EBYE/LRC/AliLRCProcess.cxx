/**************************************************************************
 * Author: Andrey Ivanov.                                                 *
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

// This class is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

#include "Riostream.h"
#include "AliLRCProcess.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
ClassImp(AliLRCProcess)

AliLRCProcess::AliLRCProcess():fIsEventOpend(kFALSE), fIsOnline(kFALSE), fDisplayInitOnDemandWarning(kTRUE), fEventCount(0),fStartForwardETA(0), fEndForwardETA(0), fStartForwardPhi(0),fEndForwardPhi(0),fStartBakwardETA(0), fEndBakwardETA(0),fStartBackwardPhi(0),fEndBackwardPhi(0),fHiPt(0),fLoPt(0),fHiMult(0),fLoMult(0),fMultBins(0),fPtBins(0),fSumPtFw(0),  fSumPtBw(0), fSumPtBw2(0), fNchFw(0), fNchBw(0),fOutList(0), fShortDef(0),fOutputSlot(0), fHistPt(0),fHistEta(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fProfNberr(0),fProfNberrPtPt(0),fProfdPtB(0),fProfTestLRC(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistPhiForward(0),fHistPtBakward(0),fHistEtaBakward(0),fHistNchBakward(0),fHistPhiBakward(0){};

AliLRCProcess::AliLRCProcess(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA): fIsEventOpend(kFALSE), fIsOnline(kFALSE), fDisplayInitOnDemandWarning(kTRUE), fEventCount(0),fStartForwardETA(0), fEndForwardETA(0), fStartForwardPhi(0),fEndForwardPhi(0),fStartBakwardETA(0), fEndBakwardETA(0),fStartBackwardPhi(0),fEndBackwardPhi(0),fHiPt(0),fLoPt(0),fHiMult(0),fLoMult(0),fMultBins(0),fPtBins(0),fSumPtFw(0),  fSumPtBw(0), fSumPtBw2(0),fNchFw(0), fNchBw(0),fOutList(0), fShortDef(0),fOutputSlot(0), fHistPt(0),fHistEta(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fProfNberr(0),fProfNberrPtPt(0),fProfdPtB(0),fProfTestLRC(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistPhiForward(0),fHistPtBakward(0),fHistEtaBakward(0),fHistNchBakward(0),fHistPhiBakward(0)
{
// Constructor with window setup makes ready-to-run processor
fEventCount=0;

SetETAWindows( _StartForwardETA, _EndForwardETA,_StartBakwardETA,_EndBakwardETA);
SetHistPtRange(0.0,4.0,200);
SetHistMultRange(0,100);  
SetForwardWindowPhi(0,2*TMath::Pi());
SetBackwardWindowPhi(0,2*TMath::Pi());
}

Bool_t AliLRCProcess::InitDataMembers()
{
// This method is actualy creating output histogramms
// Thist method  is to be called in CreateOutputObjects method of AliAnalysisTask
//cout<<" # Init for "<<fShortDef<<" this="<<this<<"\n";
if(fIsOnline)
{
Printf("Can't init data members more then one time! \n");
return kFALSE;
}
fEventCount=0;
fOutList = new TList();
fOutList->SetName(fShortDef);

   Double_t loMult,hiMult;
   
   loMult=fLoMult-0.5;
   hiMult=fHiMult+0.5;

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
  
    
  fHistNchForward = new TH1D("fHistNchForward", "N_{ch} distribution in Forward window", fMultBins, loMult, hiMult);
  fHistNchForward->GetXaxis()->SetTitle("N_{ch}");
  fHistNchForward->GetYaxis()->SetTitle("dN/dN_{ch}");
  fHistNchForward->SetMarkerStyle(kFullCircle);
  
  
  fHistPhiForward = new TH1D("fPhiForward", "Phi distribution in Forward window", 36, 0, 2*TMath::Pi());
  fHistPhiForward->GetXaxis()->SetTitle("Phi");
  fHistPhiForward->GetYaxis()->SetTitle("dN/Phi");
  fHistPhiForward->SetMarkerStyle(kFullCircle);
  
  
 
     // Bakward
  
  fHistPtBakward = new TH1D("fHistPtBakward", "P_{T} distribution in Bakward window", 100, 0.0, 5);
  fHistPtBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistEtaBakward = new TH1D("fEtaBakward", "#eta distribution in Bakward window", 200, -2, 2);
  fHistEtaBakward->GetXaxis()->SetTitle("ETA");
  fHistEtaBakward->GetYaxis()->SetTitle("dN/ETA");
  fHistEtaBakward->SetMarkerStyle(kFullCircle);
  
  
  fHistNchBakward = new TH1D("fHistNchBakward", "N_{ch} distribution in Bakward window", fMultBins, loMult, hiMult);
  fHistNchBakward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistNchBakward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistNchBakward->SetMarkerStyle(kFullCircle);

  fHistPhiBakward = new TH1D("fPhiBakward", "#Phi distribution in Bakward window", 36, 0, 2*TMath::Pi());
  fHistPhiBakward->GetXaxis()->SetTitle("Phi");
  fHistPhiBakward->GetYaxis()->SetTitle("dN/Phi");
  fHistPhiBakward->SetMarkerStyle(kFullCircle);
  

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
  
  fHistNN = new TH2D("NN","NN",fMultBins, loMult, hiMult,fMultBins, loMult, hiMult);
  fHistPtN = new TH2D("PtN","PtN",fMultBins, loMult, hiMult,fPtBins,fLoPt,fHiPt);
  fHistPtPt = new TH2D("PtPt","PtPt",fPtBins,fLoPt,fHiPt,fPtBins,fLoPt,fHiPt);
  fProfNberr = new TProfile("nber","nber",fMultBins, loMult, hiMult);
  fProfNberrPtPt = new TProfile("nberPtPt","nberPtPt",fPtBins,fLoPt,fHiPt);
  fProfdPtB = new TProfile("dPtB","Overal multievent Pt_Backward (first bin) Pt_Backward^2 (sec. bin) ",16,0.5,16.5);  
  fProfTestLRC = new TProfile("TestLRC","Test LRC calculaion via TProfile",fMultBins, loMult, hiMult);  


   // ---------- Adding data members to output list
  
  // Adding overal statistics
  
  fOutList->Add(fHistPt);
  fOutList->Add(fHistEta);
  
  //Adding LRC hists
  
  fOutList->Add(fHistNN);
  fOutList->Add(fHistPtN);
  fOutList->Add(fHistPtPt);
  fOutList->Add(fProfNberr);
  fOutList->Add(fProfNberrPtPt);
  fOutList->Add(fProfdPtB);
  fOutList->Add(fProfTestLRC);

  
  //Adding window statistics

  fOutList->Add(fHistPtForward);
  fOutList->Add(fHistEtaForward);
  fOutList->Add(fHistNchForward);
  fOutList->Add(fHistPtBakward);
  fOutList->Add(fHistEtaBakward);
  fOutList->Add(fHistNchBakward);
  fOutList->Add(fHistPhiBakward);
  fOutList->Add(fHistPhiForward);
  
  // Adding status to dPtB
  
  fProfdPtB->Fill(3 , fStartForwardETA);
  fProfdPtB->Fill(4 , fEndForwardETA);
  fProfdPtB->Fill(5 , fStartBakwardETA);
  fProfdPtB->Fill(6 , fEndBakwardETA);
  fProfdPtB->Fill(7 , fStartForwardPhi);
  fProfdPtB->Fill(8 , fEndForwardPhi);
  fProfdPtB->Fill(9 , fStartBackwardPhi);
  fProfdPtB->Fill(10 , fEndBackwardPhi);
  



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
   char str[80];
   snprintf(str,80, "TaskLRCw%3.1fto%3.1fvs%3.1fto%3.1f",fStartForwardETA,fEndForwardETA,fStartBakwardETA,fEndBakwardETA);

   fShortDef= str;

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

  void AliLRCProcess::SetHistPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins)
  {
// Sets Pt range and number of bins for Pt axis of histos
  if(fIsOnline)
  {
    Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
    return ;
  }
  // Setter for Pt range and N bins in histos
  fLoPt=LoPt;
  fHiPt=HiPt;
  fPtBins=PtBins;  
  }
  void AliLRCProcess::SetHistMultRange(Int_t LoMult,Int_t HiMult,Int_t MultBins)
  {
  // Setter for multiplicity range and N bins in histos
  if(fIsOnline)
  {
    Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
    return ;
  }
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

  void AliLRCProcess::SetOutputSlotNumber(Int_t SlotNumber)
  {
  //Sets number of output slot for LRCProcessor
  fOutputSlot=SlotNumber;
  }

//________________________________________________________________________



TList*   AliLRCProcess::CreateOutput() const
{
// Creates a link to output data TList
return fOutList;
}

 TString  AliLRCProcess::GetShortDef() const
{
return fShortDef;
}

Int_t AliLRCProcess::GetOutputSlotNumber() const
{
// Returns number of output slot for LRCProcessor
return fOutputSlot;
}

void AliLRCProcess::StartEvent()
{
// Open new Event for track by track event import
if(fIsEventOpend)                     // Check if trying to open event more than once !
	{Printf("Event is already opened! Auto finishing ! \n");
	 cout<<fShortDef<<": event count = "<<fEventCount<<" ";
	 Printf("NchF = %i,NchB = %i \n",fNchFw,fNchBw);
	 
	 FinishEvent();
	}
if(!fIsOnline)                        // Autocreating histos if InitDataMembers was not called by hand
{ 
Printf("InitDataMembers was not called by hand ! Autocreating histos...\n");
InitDataMembers();
}

fNchFw=0;
fSumPtFw=0;
fNchBw=0;
fSumPtBw=0;
fSumPtBw2=0;

fIsEventOpend=kTRUE;
}
void AliLRCProcess::AddTrackForward(Double_t Pt, Double_t Eta ,Double_t Phi)
{
// Imports track to the event directly to Forward window
if(!fIsEventOpend)
	{Printf("Event is not opened!\n");
	return;}

fNchFw++;
fSumPtFw+=Pt;
fHistPtForward->Fill(Pt);
fHistEtaForward->Fill(Eta);
fHistPhiForward->Fill(Phi);

}
void AliLRCProcess::AddTrackBackward(Double_t Pt, Double_t Eta ,Double_t Phi)
{
// Imports track to the event directly to Backward window
if(!fIsEventOpend)
	{Printf("Event is not opened!\n");
	return;}
	
fNchBw++;
fSumPtBw+=Pt;
fSumPtBw2+=Pt*Pt;
fProfdPtB->Fill(1,Pt); 
fProfdPtB->Fill(2,Pt*Pt);
fHistPtBakward->Fill(Pt);
fHistEtaBakward->Fill(Eta);
fHistPhiBakward->Fill(Phi);
}



void AliLRCProcess::AddTrackPtEta(Double_t Pt, Double_t Eta ,Double_t Phi)
{
//Track by track event import :  Imports track to the event 

if(!fIsEventOpend)
	{Printf("Event is not opened!\n");
	return;}

    //  Glabal trak data
    fHistPt->Fill(Pt);
    fHistEta->Fill(Eta);


    //Forward window
    if( (fStartForwardETA<Eta)&&(Eta<fEndForwardETA))
    if( (fStartForwardPhi<Phi)&&(Phi<fEndForwardPhi))
	AddTrackForward(Pt,Eta,Phi);

    //Backward window
    if((fStartBakwardETA<Eta)&&(Eta<fEndBakwardETA))
    if((fStartBackwardPhi<Phi)&&(Phi<fEndBackwardPhi))
	AddTrackBackward(Pt,Eta,Phi);

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
		fProfNberrPtPt->Fill(fSumPtFw,1.0/fNchBw);
		// dPtB for PtPt
		fProfdPtB->Fill(15,fSumPtBw,fNchBw);
		fProfdPtB->Fill(16,fSumPtBw2/fNchBw,fNchBw);
	}
  }
	
  
 fHistNchForward->Fill(fNchFw);
 fHistNchBakward->Fill(fNchBw);

fEventCount++;
fIsEventOpend=kFALSE;
//cout<<fShortDef<<": event count = "<<fEventCount<<" ";
//	 Printf("NchF = %i,NchB = %i",fNchFw,fNchBw);
}

