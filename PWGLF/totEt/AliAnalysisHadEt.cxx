//Created by Christine Nattrass
//University of Tennessee at Knoxville
//
// This class is designed for the analysis of the hadronic component of 
// transverse energy.  It is used by AliAnalysisTaskHadEt.
// This gets information about the hadronic component of the transverse energy 
// from tracks reconstructed in an event
// it has daughters, AliAnalysisHadEtMonteCarlo and 
// AliAnalysisHadEtReconstructed which loop over either Monte Carlo data or 
// real data to get Et

#include "AliAnalysisHadEt.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "Rtypes.h"
#include "AliPDG.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h" 
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

using namespace std;

ClassImp(AliAnalysisHadEt);


Int_t AliAnalysisHadEt::fgnumOfEtaBins = 16;
Float_t AliAnalysisHadEt::fgEtaAxis[17]={-0.78, -0.7, -0.58, -0.46, -0.34, -0.22, -0.12, -0.06, -0.0, 0.06, 0.12, 0.22, 0.34, 0.46, 0.58, 0.7, 0.78};
Int_t AliAnalysisHadEt::fgNumOfPtBins = 111;
Float_t AliAnalysisHadEt::fgPtAxis[117]=
  {0.0,0.01,0.02,0.03,0.04, 0.05, 0.06,0.07,0.08,0.09, 0.10,0.11, .12,0.13, .14,0.15, .16,0.17, .18,0.19,
   0.2, .22, .24, .26, .28, 0.30, 0.32, .34, .36, .38, 0.40, .42, .44, .46, .48,
   0.5, .52, .54, .56, .58, 0.60, 0.62, .64, .66, .68, 0.70, .72, .74, .76, .78,
   .80, .82, .84, .86, .88, 0.90, 0.92, .94, .96, .98, 1.00,1.05, 1.1,1.15, 1.2,
  1.25, 1.3,1.35,1.40,1.45, 1.50, 1.55, 1.6,1.65, 1.7, 1.75, 1.8,1.85, 1.9,1.95,
   2.0, 2.2, 2.4, 2.6, 2.8, 3.00, 3.20, 3.4, 3.6, 3.8, 4.00, 4.2, 4.4, 4.6, 4.8,
   5.0, 5.5, 6.0, 6.5, 7.0, 7.50, 8.00, 8.5, 9.0, 9.5, 10.0,12.0,14.0,16.0,18.0,
  20.0,25.0,30.0,35.0,40.0, 45.0, 50.0};
Int_t AliAnalysisHadEt::fgNumOfPtSpectraBins = 81;
Float_t AliAnalysisHadEt::fgPtSpectraAxis[82]=
  {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200};
Float_t AliAnalysisHadEt::fgResAxis[81] = {-0.150,-0.140,-0.130,-0.120,-0.110,-0.100,-0.090,-0.080,-0.070,-0.060,
					   -0.050,-0.045,-0.040,-0.035,-0.030,-0.025,-0.024,-0.023,-0.022,-0.021,
					   -0.020,-0.019,-0.018,-0.017,-0.016,-0.015,-0.014,-0.013,-0.012,-0.011,
					   -0.010,-0.009,-0.008,-0.007,-0.006,-0.005,-0.004,-0.003,-0.002,-0.001,
					   -0.000,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
					   0.010,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,
					   0.020,0.021,0.022,0.023,0.024,0.025,0.030,0.035,0.040,0.045,
					   0.050,0.060,0.070,0.080,0.090,0.100,0.110,0.120,0.130,0.140,
					   0.150,};
// Float_t AliAnalysisHadEt::fgResAxis[31] = {-0.15,-0.14,-0.13,-0.12,-0.11,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15};
Int_t AliAnalysisHadEt::fgNumOfResBins = 80;


AliAnalysisHadEt::AliAnalysisHadEt() : AliAnalysisEtCommon()
				     ,fPIDResponse(0)
				     ,fSumEt(0)
				     ,fSumEtAcc(0)
				     ,fTotEt(0)
				     ,fTotEtAcc(0)
				     ,fTotNeutralEt(0)
				     ,fTotNeutralEtAcc(0)
				     ,fTotChargedEt(0)
				     ,fTotChargedEtAcc(0)
				     ,fMultiplicity(0)
				     ,fChargedMultiplicity(0)
				     ,fNeutralMultiplicity(0)
				     ,fhistoList(0)
				     ,fGoodEvent(0)
{//default constructor

}

AliAnalysisHadEt::~AliAnalysisHadEt()
{//destructor
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsITS;
  delete fEsdtrackCutsTPC;
  delete fPIDResponse;
}

Int_t AliAnalysisHadEt::AnalyseEvent(AliVEvent *event)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  AliAnalysisEtCommon::AnalyseEvent(event);
  ResetEventValues();
  return 0;
}


void AliAnalysisHadEt::Init()
{// clear variables, set up cuts and PDG info
  AliAnalysisEtCommon::Init();
  
}

void AliAnalysisHadEt::ResetEventValues()
{//Resets event values of et to zero
  AliAnalysisEtCommon::ResetEventValues();
  fTotEt = 0;
  fTotEtAcc = 0;
  fTotNeutralEt = 0;
  fTotNeutralEtAcc = 0;
  fTotChargedEt  = 0;
  fTotChargedEtAcc = 0;
  fMultiplicity = 0;
  fChargedMultiplicity = 0;
  fNeutralMultiplicity = 0;
  
}
void AliAnalysisHadEt::CreateEtaPtHisto2D(TString name, TString title)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();

  histoname->Append(name);
  histotitle->Append(title);

  TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),fgNumOfPtBins, fgPtAxis, fgnumOfEtaBins, fgEtaAxis);
  histo->SetYTitle("#eta");
  histo->SetXTitle("p_{T}");
  histo->SetZTitle("E_{T}");
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateResolutionPtHisto2D(TString name, TString title, TString xtitle, TString ytitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();

  histoname->Append(name);
  histotitle->Append(title);

  //TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),fgNumOfPtBins, fgPtAxis, fgNumOfResBins, fgResAxis);
  TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),fgNumOfPtBins, fgPtAxis, fgNumOfResBins, fgResAxis);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->SetZTitle("Number of entries");
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreatePtHisto1D(TString name, TString title, TString xtitle, TString ytitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();

  histoname->Append(name);
  histotitle->Append(title);

  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),fgNumOfPtBins, fgPtAxis);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreatePtSpectraHisto1D(TString name, TString title, TString xtitle, TString ytitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();

  histoname->Append(name);
  histotitle->Append(title);

  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),fgNumOfPtSpectraBins, fgPtSpectraAxis);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}

void AliAnalysisHadEt::CreateHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh)
{     //creates a 1d histogram of the given dimensions and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  histoname->Append(name);
  histotitle->Append(title);
  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
}
void AliAnalysisHadEt::CreateIntHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh)
{     //creates a 1d integer histogram and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  histoname->Append(name);
  histotitle->Append(title);
  TH1I *histo = new TH1I(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh)
{     //creates a 2d histogram and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  histoname->Append(name);
  histotitle->Append(title);
  TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh,ybins,ylow,yhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateIntHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh,Int_t ybins,Int_t ylow,Int_t yhigh)
{     //creates a 2-d integer histogram and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  histoname->Append(name);
  histotitle->Append(title);
  TH2I *histo = new TH2I(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh,ybins,ylow,yhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}

void AliAnalysisHadEt::CreateEtaHisto1D(TString name, TString title)
{     //creates 1d histogram in eta and adds it to the list of histograms to be saved
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  histoname->Append(name);
  histotitle->Append(title);
  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),fgnumOfEtaBins, fgEtaAxis);
  histo->SetYTitle("E_{T}");
  histo->SetXTitle("#eta");
  histo->Sumw2();
  fhistoList->Add(histo);
  delete histoname;
  delete histotitle;
}
void AliAnalysisHadEt::FillHisto1D(TString histname, Float_t x, Float_t weight)
{//fills a 1d histogram with the name histoname with the value x and the weight "weight"
  TH1F     *histo; 
  TString  *name   = new TString();
  name->Append(histname);       
  histo = (TH1F *)fhistoList->FindObject(name->Data()); 
  if(histo){
    histo->Fill((Double_t)x, weight);
  }
  else{cerr<<"CorrelationMaker::FillHisto1D: no histogram "<< name->Data()<<endl;}
  delete name;
}
void AliAnalysisHadEt::FillHisto2D(TString histname, Float_t x, Float_t y, Float_t weight)
{//fills a 2d histogram with the name histoname with the value x and the weight "weight"
  TH2F     *histo; 
  TString  *name   = new TString();
  name->Append(histname);       
  histo = (TH2F *)fhistoList->FindObject(name->Data()); 
  if(histo){
    histo->Fill((Double_t)x,(Double_t)y, weight);
  }
  else{cerr<<"CorrelationMaker::FillHisto2D: no histogram "<< name->Data()<<endl;}
  delete name;
}


Float_t AliAnalysisHadEt::Et(TParticle *part, float mass){//function to calculate et in the same way as it would be calculated in a calorimeter
  if(mass+1000<0.01){//if no mass given return default.  The default argument is -1000
    if(TMath::Abs(part->GetPDG(0)->PdgCode())==2212 || TMath::Abs(part->GetPDG(0)->PdgCode())==2112){
      if(part->GetPDG(0)->PdgCode()==-2212 || part->GetPDG(0)->PdgCode()==-2112){//antiproton or antineutron
	//for antinucleons we specifically want to return the kinetic energy plus twice the rest mass
	return (part->Energy()+part->GetMass())*TMath::Sin(part->Theta());
      }
      if(part->GetPDG(0)->PdgCode()==2212 || part->GetPDG(0)->PdgCode()==2112){//proton or neutron
	//for nucleons we specifically want to return the kinetic energy only
	return (part->Energy()-part->GetMass())*TMath::Sin(part->Theta());
      }
    }
    else{//otherwise go to the default
      return part->Energy()*TMath::Sin(part->Theta());
    }
  }
  else{//otherwise use the mass that was given
    return (TMath::Sqrt(TMath::Power(part->P(),2.0)+TMath::Power(mass,2.0)))*TMath::Sin(part->Theta());
  }
  return 0.0;
}
Float_t AliAnalysisHadEt::Et(Float_t p, Float_t theta, Int_t pid, Short_t charge) const {//function to calculate et in the same way as it would be calculated in a calorimeter
  if(pid==AliAnalysisEtCommon::fgPiPlusCode || pid==AliAnalysisEtCommon::fgPiMinusCode){//Nothing special for pions
    return TMath::Sqrt(p*p + fgPionMass*fgPionMass) * TMath::Sin(theta);
  }
  if(pid==AliAnalysisEtCommon::fgKPlusCode || pid==AliAnalysisEtCommon::fgKMinusCode){//Nothing special for kaons
    return TMath::Sqrt(p*p + AliAnalysisEtCommon::fgKaonMass*AliAnalysisEtCommon::fgKaonMass) * TMath::Sin(theta);
  }
  if(pid==AliAnalysisEtCommon::fgEPlusCode || pid==AliAnalysisEtCommon::fgEMinusCode){//Nothing special for electrons
    return TMath::Sqrt(p*p + AliAnalysisEtCommon::fgElectronMass*AliAnalysisEtCommon::fgElectronMass) * TMath::Sin(theta);
  }
  if(pid==AliAnalysisEtCommon::fgProtonCode || pid==AliAnalysisEtCommon::fgAntiProtonCode){//But for protons we must be careful...
    if(charge<0.0){//antiprotns: kinetic energy plus twice the rest mass
      return (TMath::Sqrt(p*p + AliAnalysisEtCommon::fgProtonMass*AliAnalysisEtCommon::fgProtonMass) + AliAnalysisEtCommon::fgProtonMass) * TMath::Sin(theta);
    }
    if(charge>0.0){//antiprotns: kinetic energy only
      return (TMath::Sqrt(p*p + AliAnalysisEtCommon::fgProtonMass*AliAnalysisEtCommon::fgProtonMass) - AliAnalysisEtCommon::fgProtonMass) * TMath::Sin(theta);
    }
  }
  //cerr<<"Uh-oh!  Et not set properly!"<<endl;
  return 0.0;
}

Float_t AliAnalysisHadEt::TrueP(float pTrec) const {
  if(pTrec>1.0) return pTrec;
  return pTrec/(1-599.334*pTrec+7285.15*pTrec*pTrec)+pTrec;
}

Int_t AliAnalysisHadEt::GetCentralityBin(Int_t numberofbins,AliCentrality *centrality){
  Int_t centralitybin = -1;
  if(numberofbins<21) centralitybin= centrality->GetCentralityClass10(fCentralityMethod);
  else{
    if(numberofbins<41) centralitybin= centrality->GetCentralityClass5(fCentralityMethod);
    else{
      Float_t centpercent = centrality->GetCentralityPercentile(fCentralityMethod);
      centralitybin= centrality->GetCentralityClass5(fCentralityMethod);
      Float_t centralitybinfloat = (centpercent/2.5);
      if(centralitybin>=0){
	centralitybin =(Int_t) (centralitybinfloat);
	//cout<<" centbin "<<centralitybin<<" centrality "<<centpercent<<" "<<centpercent/2.5 <<endl;
      }
    }
  }
  //cout<<" centrality bin "<<centralitybin<<endl;
  return centralitybin;
}
Int_t AliAnalysisHadEt::GetCentralityBin(Int_t numberofbins,AliMultSelection *centrality){
  Int_t centralitybin = -1;
  Float_t lPercentile = centrality->GetMultiplicityPercentile(fCentralityMethod,kTRUE);
  if(lPercentile<0) return centralitybin;
  if(numberofbins<21){//10% bins
    centralitybin= lPercentile/10;
  }
  else{
    if(numberofbins<41){//5% bins
    centralitybin= lPercentile/5;
    }
    else{//2.5% bins
      centralitybin= lPercentile/2.5;
    }
  }
  //cout<<" centrality bin "<<centralitybin<<" percentile "<<lPercentile<<endl;
  return centralitybin;
}
