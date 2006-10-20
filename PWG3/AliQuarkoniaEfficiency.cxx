
// Root 
#include "TAxis.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TDirectory.h"

// AliRoot includes
#include "AliQuarkoniaEfficiency.h"
#include "AliLog.h"


AliQuarkoniaEfficiency::AliQuarkoniaEfficiency(Int_t quarkoniaResonance, Int_t decayChannel,
					       Int_t simParameterization):
  fEfficiencyFileName("QuarkoniaEfficiency.root"),
  fQuarkoniaResonance(quarkoniaResonance),     
  fDecayChannel(decayChannel),           
  fParameterization(simParameterization),  
  fTriggerType(kPairUnlikeApt),
  fTrigger(kFALSE),
  fEfficiency(0x0)
{

}

AliQuarkoniaEfficiency::~AliQuarkoniaEfficiency()
{
  delete fEfficiency;
}

void AliQuarkoniaEfficiency::Init()
{
  switch (fQuarkoniaResonance)  {  
  case kJpsi:
    SetTitle("Jpsi");
    break;
  case kPsiP:
    SetTitle("PsiP");
    break;
  case kUpsilon:
    SetTitle("Upsilon");
    break;
  case kUpsilonP:
    SetTitle("UpsilonP");
    break;
  case kUpsilonPP:
    SetTitle("UpsilonPP");
    break;
  case kOmega:
    SetTitle("Omega");
    break;
  case kPhi:
    SetTitle("Phi");
    break;
  }

  switch (fDecayChannel) {
  case kDimuon:
    SetName("Dimuon");
    break;
  case kDielectron:
    SetName("Dielectron");
    break;
  }

  char *param=0;
  switch (fParameterization){
  case kFlat:
    param = "Flat";
    break;
  case kCDFscaled:
    param = "CDFscaled";
    break;
  case kCDFscaledPP:
    param = "CDFscaledPP";
    break;
  }
  
  char *trig=0;
  switch (fTriggerType){
  case kSinglePlusLpt:
    trig = "SinglePlusLpt";
    break;
  case kSinglePlusHpt: 
    trig = "SinglePlusHpt";
    break;
  case kSinglePlusApt:
    trig = "SinglePlusApt";
    break;    
  case kSingleMinusLpt:
    trig = "SingleMinusLpt";
    break;    
  case kSingleMinusHpt:
    trig = "SingleMinusHpt";
    break;    
  case kSingleMinusApt:
    trig = "SingleMinusApt";
    break;    
  case kSingleUndefLpt:
    trig = "SingleUndefLpt";
    break;    
  case kSingleUndefHpt:
    trig = "SingleUndefHpt";
    break;    
  case kSingleUndefApt:
    trig = "SingleUndefApt";
    break;    
  case kPairUnlikeLpt:
    trig = "PairUnlikeLpt";
    break;    
  case kPairUnlikeHpt:
    trig = "PairUnlikeHpt";
    break;    
  case kPairUnlikeApt:
    trig = "PairUnlikeApt";
    break;    
  case kPairLikeLpt:
    trig = "PairLikeLpt";
    break;    
  case kPairLikeHpt:
    trig = "PairLikeHpt";
    break;    
  case kPairLikeApt:
    trig = "PairLikeApt";
    break;    
  }


  if(!fEfficiency) delete fEfficiency; 

  TFile efficiencyFile(fEfficiencyFileName);
  if ( efficiencyFile.IsOpen() ) {

    char quarkoniaDir[15];
    sprintf(quarkoniaDir,"%s",GetTitle());
    if (! efficiencyFile.cd(quarkoniaDir) ){
      AliError(Form("Directory %s not found in file %s \n Efficiency data for quarkonia %s and channel %s not found ",
		    quarkoniaDir,fEfficiencyFileName.Data(),GetTitle(),GetName() ));
      return;
    }
    
    char histosDir[30];
    sprintf(histosDir,"%s/%s_%s_%s",quarkoniaDir,GetTitle(),GetName(),param);
    if(! efficiencyFile.cd(histosDir) ){
      AliError(Form("Subdirectory %s/%s not found in file %s \n Efficiency data for quarkonia %s and channel %s not found ",
		    quarkoniaDir,histosDir,fEfficiencyFileName.Data(),GetTitle(),GetName() ));
      return;
    }

    char histoname[50];
    if(fTrigger) sprintf(histoname,"h%sEfficiencyPtRap_%s",GetTitle(),trig);
    else sprintf(histoname,"h%sEfficiencyPtRap",GetTitle());
    char histonameposition[99];
    sprintf(histonameposition,"%s/%s",histosDir,histoname);
    fEfficiency = (TH2F*)efficiencyFile.Get(histonameposition);

    if ( !fEfficiency ) {
      AliError(Form("Histo %s not found in file %s \n Efficiency data for quarkonia %s and channel %s not found",
		    histoname,GetTitle(), GetName() ));
    }
    else {
      fEfficiency->SetDirectory(0);
    }
    efficiencyFile.Close();

  }
  else {
    AliError(Form("File %s not found",fEfficiencyFileName.Data()));
  }

}

TH2F*  AliQuarkoniaEfficiency::GetEfficiencyHisto() const
{
  if (fEfficiency) return fEfficiency;
  else {
    AliError(Form("Efficiency data for quarkonia %s and channel %s not found",GetTitle(),GetName()));
    return 0x0;
  }
}

void  AliQuarkoniaEfficiency::GetEfficiency(Float_t rap, Float_t pT, Double_t &eff, Double_t &error)
{
  Int_t binx=0;
  Int_t biny=0;
  
  if (!fEfficiency) {
    AliError(Form("Efficiency data for quarkonia %s and channel %s not found",GetTitle(),GetName()));
  }
  else {
    if ( rap < (fEfficiency->GetXaxis())->GetXmin()  ||  
	 rap > (fEfficiency->GetXaxis())->GetXmax()  ||
	 pT  < (fEfficiency->GetYaxis())->GetXmin()  ||  
         pT  > (fEfficiency->GetYaxis())->GetXmax()   ) {
      AliInfo("Values out of range");
      eff   = 0.;
      error = 0.;
    }
    else  { 
      binx  = fEfficiency->GetXaxis()->FindBin(rap);  
      biny  = fEfficiency->GetYaxis()->FindBin(pT);
      eff   = fEfficiency->GetBinContent(binx,biny);
      error = fEfficiency->GetBinError(binx,biny);
    }
  } 
}
