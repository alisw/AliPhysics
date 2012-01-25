/* $Id$ */

// Root 
#include "TAxis.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"

// AliRoot includes
#include "AliQuarkoniaAcceptance.h"
#include "AliLog.h"



//_______________________________________________________________________
AliQuarkoniaAcceptance::AliQuarkoniaAcceptance(Int_t quarkoniaResonance, Int_t decayChannel):
  fAcceptanceFileName("$ALICE_ROOT/PWG3/QuarkoniaAcceptance.root"),
  fQuarkoniaResonance(quarkoniaResonance),
  fDecayChannel(decayChannel),
  fAcceptance(0x0)
{


}
//_______________________________________________________________________
AliQuarkoniaAcceptance::~AliQuarkoniaAcceptance()
{
  delete fAcceptance; 
}

//_______________________________________________________________________
void AliQuarkoniaAcceptance::Init()
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

  switch ( fDecayChannel) {
  case kDimuon:
    SetName("Dimuon");
    break;
  case kDielectron:
    SetName("Dielectron");
    break;
  }

  if(!fAcceptance) delete fAcceptance; 

  TFile acceptanceFile(fAcceptanceFileName);
  if ( acceptanceFile.IsOpen() ) {
    char histoname[30];
    snprintf(histoname,30,"h%s%sAccep",GetTitle(),GetName());
    fAcceptance = (TH2F*) acceptanceFile.Get(histoname);
    if ( !fAcceptance ) {
      AliError(Form("Acceptance data for quarkonia %s and channel %s not found", GetTitle(), GetName() ));
    }
    else {
      fAcceptance->SetDirectory(0);
    }
    acceptanceFile.Close();
  }
  else {
    AliError(Form("File %s not found",fAcceptanceFileName.Data()));
  }
}
//_______________________________________________________________________  
TH2F*  AliQuarkoniaAcceptance::GetAcceptanceHisto() const
{
  if (fAcceptance) return fAcceptance;
  else {
    AliError(Form("Acceptance data for quarkonia %s and channel %s not found",GetTitle(),GetName()));
    return 0x0;
  }
}
//_______________________________________________________________________  
void  AliQuarkoniaAcceptance::GetAcceptance(Float_t rap, Float_t pT, Double_t &accep, Double_t &error)
{
  Int_t binx=0;
  Int_t biny=0;

  if (!fAcceptance) {
    AliError(Form("Acceptance data for quarkonia %s and channel %s not found",GetTitle(),GetName()));
  }
  else {
    if (   rap < (fAcceptance->GetXaxis())->GetXmin()  ||  
	   rap > (fAcceptance->GetXaxis())->GetXmax()  ||
	   pT  < (fAcceptance->GetYaxis())->GetXmin()  ||  
	   pT  > (fAcceptance->GetYaxis())->GetXmax()     ) {
      AliInfo("Values out of range");
      accep = 0.;
      error = 0.;
    }
    else  { 
      binx  = fAcceptance->GetXaxis()->FindBin(rap);  
      biny  = fAcceptance->GetYaxis()->FindBin(pT);
      accep = fAcceptance->GetBinContent(binx,biny);
      error = fAcceptance->GetBinError(binx,biny);
    }
  }
} 
