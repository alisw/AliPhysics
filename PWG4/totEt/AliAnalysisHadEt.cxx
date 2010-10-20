//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
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

using namespace std;

ClassImp(AliAnalysisHadEt);
//These are from the PDG database but by making them static the code is a bit more efficient and has no problems running with the plugin
Float_t AliAnalysisHadEt::fgPionMass = 0.13957;
Float_t AliAnalysisHadEt::fgKaonMass = 0.493677;
Float_t AliAnalysisHadEt::fgProtonMass = 0.938272;
Float_t AliAnalysisHadEt::fgElectronMass = 0.000510999;
Int_t AliAnalysisHadEt::fgPiPlusCode = 211;
Int_t AliAnalysisHadEt::fgPiMinusCode = -211;
Int_t AliAnalysisHadEt::fgKPlusCode = 321;
Int_t AliAnalysisHadEt::fgKMinusCode = -321;
Int_t AliAnalysisHadEt::fgProtonCode = 2212;
Int_t AliAnalysisHadEt::fgAntiProtonCode = -2212;
Int_t AliAnalysisHadEt::fgLambdaCode = 3122;
Int_t AliAnalysisHadEt::fgAntiLambdaCode = -3122;
Int_t AliAnalysisHadEt::fgK0SCode = 310;
Int_t AliAnalysisHadEt::fgOmegaCode = 3334;
Int_t AliAnalysisHadEt::fgAntiOmegaCode = -3334;
Int_t AliAnalysisHadEt::fgXi0Code = 3322;
Int_t AliAnalysisHadEt::fgAntiXi0Code = -3322;
Int_t AliAnalysisHadEt::fgXiCode = 3312;
Int_t AliAnalysisHadEt::fgAntiXiCode = -3312;
Int_t AliAnalysisHadEt::fgSigmaCode = 3112;
Int_t AliAnalysisHadEt::fgAntiSigmaCode = -3112;
Int_t AliAnalysisHadEt::fgK0LCode = 130;
Int_t AliAnalysisHadEt::fgNeutronCode = 2112;
Int_t AliAnalysisHadEt::fgAntiNeutronCode = -2112;
Int_t AliAnalysisHadEt::fgEPlusCode = -11;
Int_t AliAnalysisHadEt::fgEMinusCode = 11;
Int_t AliAnalysisHadEt::fgGammaCode = 22;
Int_t AliAnalysisHadEt::fgPi0Code = 111;
Int_t AliAnalysisHadEt::fgEtaCode = 221;
Int_t AliAnalysisHadEt::fgOmega0Code = 223;


// Int_t AliAnalysisHadEt::fgnumOfEtaBins = 46;
// Float_t AliAnalysisHadEt::fgEtaAxis[47]={-0.78, -0.74, -0.7, -0.66, -0.62, -0.58, -0.54, -0.5, -0.46, -0.42, -0.38, -0.34, -0.3, -0.26, -0.22, -0.18, -0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, -0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.18, 0.22, 0.26, 0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7, 0.74, 0.78};
//reduction in the number of bins
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

Float_t AliAnalysisHadEt::fgPtTPCCutOff = 0.15;
Float_t AliAnalysisHadEt::fgPtITSCutOff = 0.10;

AliAnalysisHadEt::AliAnalysisHadEt() :
        fHistogramNameSuffix("")
	,fCuts(0)
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
        ,fEsdtrackCutsITSTPC(0)
        ,fEsdtrackCutsTPC(0)
	,fEsdtrackCutsITS(0)
        ,fhistoList(0)
{//default constructor

}

AliAnalysisHadEt::~AliAnalysisHadEt()
{//destructor
  delete fCuts;
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsITS;
  delete fEsdtrackCutsTPC;
}

Int_t AliAnalysisHadEt::AnalyseEvent(AliVEvent *event)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  cout << "This event has " << event->GetNumberOfTracks() << " tracks" << endl;
  ResetEventValues();
  return 0;
}

void AliAnalysisHadEt::FillOutputList()
{//fill the output histogram list with histograms in all AliAnalysisHadEt's
}

void AliAnalysisHadEt::Init()
{// clear variables, set up cuts and PDG info

}

void AliAnalysisHadEt::CreateHistograms()
{//creates histograms included in all AliAnalysisHadEt's
}

void AliAnalysisHadEt::FillHistograms()
{//Fills histograms filled for all AliAnalysisHadEt's
}

void AliAnalysisHadEt::ResetEventValues()
{//Resets event values of et to zero
  fTotEt = 0;
  fTotEtAcc = 0;
  fTotNeutralEt = 0;
  fTotNeutralEtAcc = 0;
  fTotChargedEt  = 0;
  fTotChargedEtAcc = 0;
  fMultiplicity = 0;
  fChargedMultiplicity = 0;
  fNeutralMultiplicity = 0;
  
  if (!fCuts) { // some Init's needed
    cout << __FILE__ << ":" << __LINE__ << " : Init " << endl;
    if (!fCuts) {
      cout << " setting up Cuts " << endl;
      fCuts = new AliAnalysisEtCuts();
    }
  }
}

void AliAnalysisHadEt::SetParticleCodes()
{  //the codes are defined in $ROOTSYS/etc/pdg_table.txt
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
  if(pid==fgPiPlusCode || pid==fgPiMinusCode){//Nothing special for pions
    return TMath::Sqrt(p*p + fgPionMass*fgPionMass) * TMath::Sin(theta);
  }
  if(pid==fgKPlusCode || pid==fgKMinusCode){//Nothing special for kaons
    return TMath::Sqrt(p*p + fgKaonMass*fgKaonMass) * TMath::Sin(theta);
  }
  if(pid==fgEPlusCode || pid==fgEMinusCode){//Nothing special for electrons
    return TMath::Sqrt(p*p + fgElectronMass*fgElectronMass) * TMath::Sin(theta);
  }
  if(pid==fgProtonCode || pid==fgAntiProtonCode){//But for protons we must be careful...
    if(charge<0.0){//antiprotns: kinetic energy plus twice the rest mass
      return (TMath::Sqrt(p*p + fgProtonMass*fgProtonMass) + fgProtonMass) * TMath::Sin(theta);
    }
    if(charge>0.0){//antiprotns: kinetic energy only
      return (TMath::Sqrt(p*p + fgProtonMass*fgProtonMass) - fgProtonMass) * TMath::Sin(theta);
    }
  }
  cerr<<"Uh-oh!  Et not set properly!"<<endl;
  return 0.0;
}
