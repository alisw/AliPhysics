#include "AliGFWWeights.h"
#include "TMath.h"
AliGFWWeights::AliGFWWeights():
  fDataFilled(kFALSE),
  fMCFilled(kFALSE),
  fW_data(0),
  fW_mcrec(0),
  fW_mcgen(0),
  fEffInt(0),
  fAccInt(0)
{
};
AliGFWWeights::~AliGFWWeights()
{
  delete fW_data;
  delete fW_mcrec;
  delete fW_mcgen;
  delete fEffInt;
  delete fAccInt;
  delete [] fbinsPtDefault;
};
void AliGFWWeights::Init(Bool_t AddData, Bool_t AddMC)
{
  Double_t binsPtDefault[fNbinsPtDefault+1] = {
    0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30,
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
    0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40,
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60,
    2.80, 3.00, 3.40, 3.80, 4.20, 4.60, 5.20, 5.80, 6.80,
    8.00, 12.,  16., 20};
  fbinsPtDefault = binsPtDefault;
  fW_data = new TObjArray();
  fW_data->SetName("AliGFWWeights_Data");
  fW_mcrec = new TObjArray();
  fW_mcrec->SetName("AliGFWWeights_MCRec");
  fW_mcgen = new TObjArray();
  fW_mcgen->SetName("AliGFWWeights_MCGen");
  fW_data->SetOwner(kTRUE);
  fW_mcrec->SetOwner(kTRUE);
  fW_mcgen->SetOwner(kTRUE);
  fDataFilled = kFALSE;
  fMCFilled = kFALSE;
    if(AddData) {
      //const char *tnd = GetBinName(0,0,Form("data_%s",this->GetName()));
	fW_data->Add(new TH3D("data",";#varphi;#eta;v_{z}",60,0,TMath::TwoPi(),64,-1.6,1.6,40,-10,10));
    };
    if(AddMC) {
      //const char *tnr = GetBinName(0,0,"mcrec"); //all integrated over cent. anyway
      //const char *tng = GetBinName(0,0,"mcgen"); //all integrated over cent. anyway
      fW_mcrec->Add(new TH3D("reco",";#it{p}_{T};#eta;v_{z}",fNbinsPtDefault,0,20,64,-1.6,1.6,40,-10,10));
      fW_mcgen->Add(new TH3D("gene",";#it{p}_{T};#eta;v_{z}",fNbinsPtDefault,0,20,64,-1.6,1.6,40,-10,10));
      ((TH3D*)fW_mcrec->At(fW_mcrec->GetEntries()-1))->GetXaxis()->Set(fNbinsPtDefault,fbinsPtDefault);
      ((TH3D*)fW_mcgen->At(fW_mcgen->GetEntries()-1))->GetXaxis()->Set(fNbinsPtDefault,fbinsPtDefault);
    };
};

void AliGFWWeights::Fill(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype) {
  TObjArray *tar=0;
  const char *pf="";
  if(htype==0) { tar = fW_data; pf = "data"; };
  if(htype==1) { tar = fW_mcrec; pf = "reco"; };
  if(htype==2) { tar = fW_mcgen; pf = "gene"; };
  if(!tar) return;
  TH3D *th3 = (TH3D*)tar->FindObject(pf); //pT bin 0, V0M bin 0, since all integrated
  if(!th3) {
    if(!htype) tar->Add(new TH3D(pf,";#varphi;#eta;v_{z}",60,0,TMath::TwoPi(),64,-1.6,1.6,40,-10,10)); //0,0 since all integrated
    th3 = (TH3D*)tar->At(tar->GetEntries()-1);
  };
  th3->Fill(htype?pt:phi,eta,vz);
};
Double_t AliGFWWeights::GetWeight(Double_t phi, Double_t eta, Double_t vz, Double_t pt, Double_t cent, Int_t htype) {
  TObjArray *tar=0;
  const char *pf="";
  if(htype==0) { tar = fW_data; pf = "data"; };
  if(htype==1) { tar = fW_mcrec; pf = "reco"; };
  if(htype==2) { tar = fW_mcgen; pf = "gene"; };
  if(!tar) return 1;
  TH3D *th3 = (TH3D*)tar->FindObject(pf);
  if(!th3) return 1;//-1;
  Int_t xind = th3->GetXaxis()->FindBin(htype?pt:phi);
  Int_t etaind = th3->GetYaxis()->FindBin(eta);
  Int_t vzind = th3->GetZaxis()->FindBin(vz);
  Double_t weight = th3->GetBinContent(xind, etaind, vzind);
  if(weight!=0) return 1./weight;
  return 1;
};
Double_t AliGFWWeights::FindMax(TH3D *inh, Int_t &ix, Int_t &iy, Int_t &iz) {
  Double_t maxv=inh->GetBinContent(1,1,1);
  for(Int_t i=1;i<=inh->GetNbinsX();i++)
    for(Int_t j=1;j<=inh->GetNbinsY();j++)
      for(Int_t k=1;k<=inh->GetNbinsZ();k++)
	if(inh->GetBinContent(i,j,k)>maxv) {
	  ix=i;
	  iy=j;
	  iz=k;
	  maxv=inh->GetBinContent(i,j,k);
	};
  return maxv;
};
void AliGFWWeights::MCToEfficiency() {
  if(fW_mcgen->GetEntries()<1) {
    printf("MC gen. array empty. This is probably because effs. have been calculated and the generated particle histograms have been cleared out!\n");
    return;
  };
  for(Int_t i=0;i<fW_mcrec->GetEntries();i++) {
    TH3D *hr = (TH3D*)fW_mcrec->At(i);
    TH3D *hg = (TH3D*)fW_mcgen->At(i);
    hr->Sumw2();
    hg->Sumw2();
    hr->Divide(hg);
  };
  fW_mcgen->Clear();
};
void AliGFWWeights::CreateNUA(Bool_t IntegrateOverCentAndPt) {
  if(!IntegrateOverCentAndPt) {
    printf("Method is outdated! NUA is integrated over centrality and pT. Quit now, or the behaviour will be bad\n");
    return;
  };
  TH3D *h3;
  TH1D *h1;
  if(fW_data->GetEntries()<1) return;
  if(IntegrateOverCentAndPt) {
    if(fAccInt) delete fAccInt;
    fAccInt = (TH3D*)fW_data->At(0)->Clone("IntegratedAcceptance");
    fAccInt->RebinY(2);
    fAccInt->RebinZ(5);
    fAccInt->Sumw2();
    for(Int_t etai=1;etai<=fAccInt->GetNbinsY();etai++) {
      fAccInt->GetYaxis()->SetRange(etai,etai);
      if(fAccInt->Integral()<1) continue;
      for(Int_t vzi=1;vzi<=fAccInt->GetNbinsZ();vzi++) {
	fAccInt->GetZaxis()->SetRange(vzi,vzi);
	if(fAccInt->Integral()<1) continue;
	h1 = (TH1D*)fAccInt->Project3D("x");
	Double_t maxv = h1->GetMaximum();
	for(Int_t phii=1;phii<=h1->GetNbinsX();phii++)
	  fAccInt->SetBinContent(phii,etai,vzi,fAccInt->GetBinContent(phii,etai,vzi)/maxv);
	delete h1;
      };
      fAccInt->GetZaxis()->SetRange(1,fAccInt->GetNbinsZ());
    };
    fAccInt->GetYaxis()->SetRange(1,fAccInt->GetNbinsY());
    return;
  };
};
void AliGFWWeights::CreateNUE(Bool_t IntegrateOverCentrality) {
  if(!IntegrateOverCentrality) {
    printf("Method is outdated! NUE is integrated over centrality. Quit now, or the behaviour will be bad\n");
    return;
  };
  TH3D *num=0;
  TH3D *den=0;
  if(fW_mcrec->GetEntries()<1 || fW_mcgen->GetEntries()<1) return;
  if(IntegrateOverCentrality) {
    num=(TH3D*)fW_mcrec->At(0);//->Clone(Form("temp_%s",fW_mcrec->At(0)->GetName()));
    den=(TH3D*)fW_mcgen->At(0);//->Clone(Form("temp_%s",fW_mcgen->At(0)->GetName()));
    num->Sumw2();
    den->Sumw2();
    num->RebinY(2);
    den->RebinY(2);
    num->RebinZ(5);
    den->RebinZ(5);
    fEffInt = (TH3D*)num->Clone("Efficiency_Integrated");
    fEffInt->Divide(den);
    return;
  };
};
void AliGFWWeights::ReadAndMerge(const char *filelinks) {
  FILE *flist = fopen(filelinks,"r");
  char str[150];
  Int_t nFiles=0;
  while(fscanf(flist,"%s\n",str)==1) nFiles++;
  rewind(flist);
  if(nFiles==0) {
    printf("No files to read!\n");
    return;
  };
  if(!fW_data) {
    fW_data = new TObjArray();
    fW_data->SetName("Weights_Data");
    fW_data->SetOwner(kTRUE);
  };
  if(!fW_mcrec) {
    fW_mcrec = new TObjArray();
    fW_mcrec->SetName("Weights_MCRec");
    fW_mcrec->SetOwner(kTRUE);
  };
  if(!fW_mcgen) {
    fW_mcgen = new TObjArray();
    fW_mcgen->SetName("Weights_MCGen");
    fW_mcgen->SetOwner(kTRUE);
  };
  fDataFilled = kFALSE;
  fMCFilled = kFALSE;
  TFile *tf=0;
  for(Int_t i=0;i<nFiles;i++) {
    Int_t trash = fscanf(flist,"%s\n",str);
    tf = new TFile(str,"READ");
    if(tf->IsZombie()) {
      printf("Could not open file %s!\n",str);
      tf->Close();
      continue;
    };
    TList *tl = (TList*)tf->Get("OutputList");
    AliGFWWeights *tw = (AliGFWWeights*)tl->FindObject(this->GetName());
    if(!tw) {
      printf("Could not fetch weights object from %s\n",str);
      tf->Close();
      continue;
    };
    AddArray(fW_data,tw->GetDataArray());
    AddArray(fW_mcrec,tw->GetRecArray());
    AddArray(fW_mcgen,tw->GetGenArray());
    tf->Close();
    delete tw;
  };
};
void AliGFWWeights::AddArray(TObjArray *targ, TObjArray *sour) {
  if(!sour) {
    printf("Source array does not exist!\n");
    return;
  };
  for(Int_t i=0;i<sour->GetEntries();i++) {
    TH3D *sourh = (TH3D*)sour->At(i);
    TH3D *targh = (TH3D*)targ->FindObject(sourh->GetName());
    if(!targh) {
      targh = (TH3D*)sourh->Clone(sourh->GetName());
      targh->SetDirectory(0);
      targ->Add(targh);
    } else 
      targh->Add(sourh);
  };
};
