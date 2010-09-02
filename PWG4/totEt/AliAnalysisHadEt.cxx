//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#include "AliAnalysisHadEt.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
#include "AliVEvent.h"

using namespace std;

ClassImp(AliAnalysisHadEt);


Int_t AliAnalysisHadEt::numOfEtaBins = 46;
Float_t AliAnalysisHadEt::etaAxis[47]={-0.78, -0.74, -0.7, -0.66, -0.62, -0.58, -0.54, -0.5, -0.46, -0.42, -0.38, -0.34, -0.3, -0.26, -0.22, -0.18, -0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, -0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.18, 0.22, 0.26, 0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7, 0.74, 0.78};
Int_t AliAnalysisHadEt::numOfPtBins = 111;
Float_t AliAnalysisHadEt::ptAxis[117]=
  {0.0,0.01,0.02,0.03,0.04, 0.05, 0.06,0.07,0.08,0.09, 0.10,0.11, .12,0.13, .14,0.15, .16,0.17, .18,0.19,
   0.2, .22, .24, .26, .28, 0.30, 0.32, .34, .36, .38, 0.40, .42, .44, .46, .48,
   0.5, .52, .54, .56, .58, 0.60, 0.62, .64, .66, .68, 0.70, .72, .74, .76, .78,
   .80, .82, .84, .86, .88, 0.90, 0.92, .94, .96, .98, 1.00,1.05, 1.1,1.15, 1.2,
  1.25, 1.3,1.35,1.40,1.45, 1.50, 1.55, 1.6,1.65, 1.7, 1.75, 1.8,1.85, 1.9,1.95,
   2.0, 2.2, 2.4, 2.6, 2.8, 3.00, 3.20, 3.4, 3.6, 3.8, 4.00, 4.2, 4.4, 4.6, 4.8,
   5.0, 5.5, 6.0, 6.5, 7.0, 7.50, 8.00, 8.5, 9.0, 9.5, 10.0,12.0,14.0,16.0,18.0,
  20.0,25.0,30.0,35.0,40.0, 45.0, 50.0}; 

AliAnalysisHadEt::AliAnalysisHadEt() :
        fHistogramNameSuffix("")
        ,fPdgDB(0)
        ,PiPlusCode(0)
        ,PiMinusCode(0)
        ,KPlusCode(0)
        ,KMinusCode(0)
        ,ProtonCode(0)
        ,AntiProtonCode(0)
        ,LambdaCode(0)
        ,AntiLambdaCode(0)
        ,K0SCode(0)
        ,OmegaCode(0)
        ,AntiOmegaCode(0)
        ,Xi0Code(0)
        ,AntiXi0Code(0)
        ,XiCode(0)
        ,AntiXiCode(0)
        ,SigmaCode(0)
        ,AntiSigmaCode(0)
        ,K0LCode(0)
        ,NeutronCode(0)
        ,AntiNeutronCode(0)
        ,EPlusCode(0)
        ,EMinusCode(0)
        ,PionMass(0)
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
        ,fEtaCut(EtCommonCuts::kEtaCut)
        ,fEtaCutAcc(0)
				  //,fPhiCutAccMin(0)
				  //,fPhiCutAccMax(360.)
        ,fVertexXCut(0)
        ,fVertexYCut(0)
        ,fVertexZCut(0)
				  ,fIPxyCut(0)
				  ,fIPzCut(0)
				  //,fSingleCellEnergyCut(0)
				  //,fClusterEnergyCut(EtCommonCuts::kClusterEnergyCut)
				  //,fTrackPtCut(EtCommonCuts::kTrackPtCut)
        ,esdtrackCutsITSTPC(0)
        ,esdtrackCutsTPC(0)
	,esdtrackCutsITS(0)
        ,histoList(0)
{

}

AliAnalysisHadEt::~AliAnalysisHadEt()
{

}

Int_t AliAnalysisHadEt::AnalyseEvent(AliVEvent *event)
{
  //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  cout<<"This event has "<<event->GetNumberOfTracks()<<" tracks"<<endl;
  return 0;
}

void AliAnalysisHadEt::FillOutputList()
{
}

void AliAnalysisHadEt::Init()
{
  if(!fPdgDB) fPdgDB = new TDatabasePDG();
  //the codes are defined in $ROOTSYS/etc/pdg_table.txt
  PionMass = fPdgDB->GetParticle("pi+")->Mass();
  PiPlusCode = fPdgDB->GetParticle("pi+")->PdgCode();
    PiMinusCode = fPdgDB->GetParticle("pi-")->PdgCode();
    KPlusCode = fPdgDB->GetParticle("K+")->PdgCode();
    KMinusCode = fPdgDB->GetParticle("K-")->PdgCode();
    ProtonCode = fPdgDB->GetParticle("proton")->PdgCode();
    AntiProtonCode = fPdgDB->GetParticle("antiproton")->PdgCode();
    LambdaCode = fPdgDB->GetParticle("Lambda0")->PdgCode();
    AntiLambdaCode = fPdgDB->GetParticle("Lambda0_bar")->PdgCode();
    K0SCode = fPdgDB->GetParticle("K_S0")->PdgCode();
    OmegaCode = fPdgDB->GetParticle("Omega-")->PdgCode();
    AntiOmegaCode = fPdgDB->GetParticle("Omega+")->PdgCode();
    Xi0Code = fPdgDB->GetParticle("Xi0")->PdgCode();
    AntiXi0Code = fPdgDB->GetParticle("Xi0_bar")->PdgCode();
    XiCode = fPdgDB->GetParticle("Xi-")->PdgCode();
    AntiXiCode = fPdgDB->GetParticle("Xi-_bar")->PdgCode();
    SigmaCode = fPdgDB->GetParticle("Sigma-")->PdgCode();
    AntiSigmaCode = fPdgDB->GetParticle("Sigma+")->PdgCode();
    K0LCode = fPdgDB->GetParticle("K_L0")->PdgCode();
    NeutronCode = fPdgDB->GetParticle("neutron")->PdgCode();
    AntiNeutronCode = fPdgDB->GetParticle("antineutron")->PdgCode();
    EPlusCode = fPdgDB->GetParticle("e+")->PdgCode();
    EMinusCode = fPdgDB->GetParticle("e-")->PdgCode();
}

void AliAnalysisHadEt::CreateHistograms()
{
}

void AliAnalysisHadEt::FillHistograms()
{
}

void AliAnalysisHadEt::ResetEventValues()
{
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
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();

  histoname->Append(name);
  histotitle->Append(title);
  //TH2F *h1 = new TH2F("h1", "Histogram with Gaussian random distribution", numOfPtBins, ptBinsArray, numOfEtaBins, etaBinsArray);

  TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),numOfPtBins, ptAxis, numOfEtaBins, etaAxis);
  histo->SetYTitle("#eta");
  histo->SetXTitle("p_{T}");
  histo->SetZTitle("E_{T}");
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}

void AliAnalysisHadEt::CreateHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh)
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  
  //cout<<"creating "<<name<<endl;

  histoname->Append(name);
  histotitle->Append(title);
  // printf("%s \n ",histoname->Data());
  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateIntHisto1D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh)
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  
  //cout<<"creating "<<name<<endl;

  histoname->Append(name);
  histotitle->Append(title);
  // printf("%s \n ",histoname->Data());
  TH1I *histo = new TH1I(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Float_t xlow,Float_t xhigh,Int_t ybins,Float_t ylow,Float_t yhigh)
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  
  //cout<<"creating "<<name<<endl;

  histoname->Append(name);
  histotitle->Append(title);
  // printf("%s \n ",histoname->Data());
  TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh,ybins,ylow,yhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::CreateIntHisto2D(TString name, TString title, TString xtitle, TString ytitle,Int_t xbins, Int_t xlow,Int_t xhigh,Int_t ybins,Int_t ylow,Int_t yhigh)
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  
  //cout<<"creating "<<name<<endl;

  histoname->Append(name);
  histotitle->Append(title);
  // printf("%s \n ",histoname->Data());
  TH2I *histo = new TH2I(histoname->Data(),histotitle->Data(),xbins,xlow,xhigh,ybins,ylow,yhigh);
  histo->SetYTitle(ytitle);
  histo->SetXTitle(xtitle);
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}

void AliAnalysisHadEt::CreateEtaHisto1D(TString name, TString title)
{     
  TString *histoname   = new TString();
  TString *histotitle   = new TString();
  

  histoname->Append(name);
  histotitle->Append(title);
  TH1F *histo = new TH1F(histoname->Data(),histotitle->Data(),numOfEtaBins, etaAxis);
  histo->SetYTitle("E_{T}");
  histo->SetXTitle("#eta");
  histo->Sumw2();
  histoList->Add(histo);
  delete histoname;
  delete histotitle;
    
}
void AliAnalysisHadEt::FillHisto1D(TString histname, Float_t x, Float_t weight)
{
  TH1F     *histo; 
  TString  *name   = new TString();

  name->Append(histname);       
  histo = (TH1F *)histoList->FindObject(name->Data()); 
  if(histo){
    histo->Fill((Double_t)x, weight);
  }
  else{cerr<<"CorrelationMaker::FillHisto1D: no histogram "<< name->Data()<<endl;}
  delete name;
}
void AliAnalysisHadEt::FillHisto2D(TString histname, Float_t x, Float_t y, Float_t weight)
{
  TH2F     *histo; 
  TString  *name   = new TString();
  
  name->Append(histname);       
  histo = (TH2F *)histoList->FindObject(name->Data()); 
  if(histo){
    histo->Fill((Double_t)x,(Double_t)y, weight);
  }
  else{cerr<<"CorrelationMaker::FillHisto2D: no histogram "<< name->Data()<<endl;}
  delete name;
}


Float_t AliAnalysisHadEt::Et(TParticle *part, float mass){
  if(mass == -1000){//if no mass given return default 
    if(TMath::Abs(part->GetPDG(0)->PdgCode())==2212 || TMath::Abs(part->GetPDG(0)->PdgCode())==2112){
      if(part->GetPDG(0)->PdgCode()==-2212 || part->GetPDG(0)->PdgCode()==-2112){//antiproton or antineutron
	//for antinucleons we specifically want to return the kinetic energy plus twice the rest mass
	return (part->Energy()+part->GetMass())*TMath::Sin(part->Theta());
      }
      if(part->GetPDG(0)->PdgCode()==2212 || part->GetPDG(0)->PdgCode()==2112){//antiproton or antineutron
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
