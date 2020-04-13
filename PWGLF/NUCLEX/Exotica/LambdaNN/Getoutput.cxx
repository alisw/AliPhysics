#include <TList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TNtupleD.h>
#include <TObjString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TDatabasePDG.h>
#include <TGrid.h>
#include <TLegend.h>
#include <TStyle.h>
#include <./Getoutput.h>

const Int_t Getoutput::fgArrSize=5;
Double_t Getoutput::fgLH3Mass = 2.808921;
Double_t Getoutput::fgPionMass = 0.13957;
Double_t Getoutput::fgEleMass = 0.00051;
Double_t Getoutput::fgProtMass = 0.93827;

Getoutput::Getoutput(const char *filename,Int_t nsTri, Int_t nsPi) {
 
 TString fOutFilename; 
fOutFilename+=filename;
 fOutputList = new TList();
 fOutputList->SetOwner(kTRUE);
 fRejectLikesign = kTRUE;
 fInputList= new TList();
 fInputList->SetOwner(kTRUE);
 fIncludePidTOF=kFALSE;//FALSE;
 fUseAODCut=kFALSE;//kTRUE;
 fCaio=kFALSE;
 fIsMC=kFALSE;
 f3HPcut=3;
 f3Hcharge=2;
 TFile *grFile = TFile::Open("outAll.root"); // out2.root has 2sigma cut instead of 3 sigma cut for tritons;
 grTriPs =(TGraph*)grFile->Get(Form("tP%is",nsTri)); 
 grTriMs =(TGraph*)grFile->Get(Form("tM%is",nsTri)); 
 grPiPs =(TGraph*)grFile->Get(Form("piP%is",nsPi)); 
 grPiMs =(TGraph*)grFile->Get(Form("piM%is",nsPi)); 
}

void Getoutput::Process(){

if(f3Hcharge==2){
LoopOverV0(-1);
LoopOverV0(1);
}
else LoopOverV0(f3Hcharge);
}

void Getoutput::StoreOutputData()
{

 printf("storing...\n");
 for(Int_t elem =0; elem < fOutputList->GetEntries(); elem++){
  TObject *obj = fOutputList->At(elem);
  fOutFile->WriteObjectAny(obj,obj->ClassName(),obj->GetName());
 }

   fOutFile->Close();
 //if(fInputList) ClearInputData();
printf("end storing\n"); 
}

bool Getoutput::LoadParams(const char *paramfile)
{

 if(!gSystem->IsFileInIncludePath(paramfile))
  {
  printf("no params.txt file in %s, exiting \n",gSystem->pwd());
  return kFALSE;
 }

 FILE *infile = fopen(paramfile,"r");
 //kParMinP, kParMinPv0, kParMaxP3H, kParPiLim, kPar3hLim, 
 //kParNclusITS, kParNsigmaPID, kParNsigmaTOFmass, kParDcaTriZ, kParCosP, 
 //kParV0Dca kK0MassLow, kK0MassHigh, kLambdaMassLow, kLambdaMassHigh,kGammaMassHigh
 Int_t scanId=-1;
 TString paramNames[20]={"Pmin","Pv0Min","P3HMax","PpiMax","p3HMin",
			 "nITSclus","nSigmaPID","nSigmaTOFm2","Dca3Hz","CosThetaP",
                         "V0dca","K0mLow","K0mHigh","LmLow","LmHigh","GmHigh",
			 "3HTOFpid","IsMC","PIDyear","3HpidPlim"};
 
 scanId = fscanf(infile, "%f   %f   %f   %f  %f\n",&param[kParMinP],&param[kParMinPv0],&param[kParMaxP3H],&param[kParPiLim],&param[kPar3hLim]);
 scanId = fscanf(infile, "%f   %f   %f   %f  %f\n",&param[kParNclusITS],&param[kParNsigmaPID],&param[kParNsigmaTOFmass],&param[kParDcaTriZ],&param[kParCosP]);
 scanId = fscanf(infile, "%f   %f   %f   %f  %f   %f\n",&param[kParV0Dca],&param[kK0MassLow],&param[kK0MassHigh],&param[kLambdaMassLow],&param[kLambdaMassHigh],&param[kGammaMassHigh]);
scanId = fscanf(infile, "%f   %f   %f   %f\n",&param[kTOFpid],&param[kIsMc],&param[kPIDResponseYear],&param[k3HPlim]);
 for(Int_t i=0; i<20; i++) {
   if(i%5==0) printf("\n");
 printf("%10s=%2.2f   ",paramNames[i].Data(),param[i]);
 }
 printf("\n");
 if(param[kIsMc]>0.5) {
   SetMC();
   printf("MC is ON, ");
 }
 if(param[kTOFpid]>0.5) {
    SetTOFpid();
    printf(" TOF pid is ON, ");
 }
 if(param[k3HPlim]<0) {
    SetAODCuts();
    printf("Using AOD analysis cuts");
 } else f3HPcut = param[k3HPlim];
 printf("\n");
 return kTRUE;
}


bool Getoutput::LoadFile(const char *inputfile, const char *listName){

 if(!gGrid && !gSystem->IsFileInIncludePath(inputfile)){
  printf("no %s file in %s, exiting \n",inputfile,gSystem->pwd());
  return kFALSE;
 }


 TFile *file = TFile::Open(inputfile);
 if(!file){
  printf("no %s file",inputfile);
  return kFALSE;
 }
 TList *l = (TList*)file->Get(listName);
 if(!l){
  printf("no %s list, please check the listname \n",listName);
  return kFALSE;
 }

 //fInputList = (TList*)l->Clone("locLNNlist");
  for(Int_t ie =0; ie < l->GetEntries(); ie++){
  fInputList->AddLast(l->At(ie));
  }

 if(!fInputList){
  printf("no input list after cloning, exiting \n");
  return kFALSE;
 }
 else {
    TH1F *hIn = (TH1F *)(fInputList->FindObject("fHistEventMultiplicity"));
    TH1F *hOut = (TH1F*)(fOutputList->FindObject("hEventMultiplicity"));
    if(!hOut){
      printf("Please book histograms before loading the file\n");
      return kFALSE;
    }
    hOut->Add (hIn);
 }

 file->Delete();

    return kTRUE;
}


void Getoutput::ClearInputData(){
 printf("clear input data\n");
  if(!fInputList){
    //printf("pointer to the list %p \n",fInputList);
    fInputList->Clear();
   //delete fInputList;
  } else return;
}


bool Getoutput::EventSelectionSimple (Double_t *arr)
{
 if(arr[kV0mom] < param[kParMinPv0]) return kFALSE;
 if (arr[kDcaTriZ] > param[kParDcaTriZ]) return kFALSE;
 if (arr[kCosP] < param[kParCosP]) return kFALSE;
 if (arr[kV0Dca] > param[kParV0Dca]) return kFALSE;
}

bool Getoutput::EventSelectionAOD (Double_t *arr)
{
 if(arr[kV0mom] < param[kParMinPv0]) return kFALSE;
 if (arr[kDcaTriZ] > param[kParDcaTriZ]) return kFALSE;
 if (arr[kCosP] < param[kParCosP]) return kFALSE;
 if (arr[kV0Dca] > param[kParV0Dca]) return kFALSE;
 Int_t itsClus = (Int_t) arr[kNclusITS];
 
 if (itsClus < 100)
 {
   if (itsClus < TMath::Nint(param[kParNclusITS])) return kFALSE;
   //else printf("n its clus pion : %4i  (nITS param %3.3f)\n",itsClus,param[kParNclusITS]);
 }
 else if (itsClus % 100 < TMath::Nint(param[kParNclusITS])) return kFALSE;
 //else  printf("n its clus 3H : %4i(nITS param %3.3f \n",itsClus,param[kParNclusITS]);
 Double_t lim1 = 0.8 * (1 - arr[kV0dcaD] / 0.8);
 if (arr[kDcaTriXY] > lim1) return kFALSE;
 Double_t lim2 = 0.7 * (1 - arr[kDcaPi] / 7.);
 if (arr[kV0dcaD] > lim2) return kFALSE;
 return kTRUE;
}

Double_t Getoutput::GetInvMass (TVector3 vPos, TVector3 vNeg, Double_t mPos, Double_t mNeg)
{
 TLorentzVector v0Pos;
 v0Pos.SetXYZM (vPos.X (), vPos.Y (), vPos.Z (), mPos);
 TLorentzVector v0Neg;
 v0Neg.SetXYZM (vNeg.X (), vNeg.Y (), vNeg.Z (), mNeg);
 TLorentzVector v0Sum = v0Pos + v0Neg;
 return v0Sum.M ();
}

void Getoutput::BookOutputData(const char *name){
  printf("initialò filename %s \n",fOutFilename.Data());
 fOutFilename.ReplaceAll(".root","");
 fOutFilename+=Form("3H%i",f3Hcharge);
 if(fUseAODCut)fOutFilename+=".AODcut.";
 if(fIncludePidTOF)fOutFilename+=".3HTOFpid.";
 if(fCaio) fOutFilename+="CaioCuts.";
 fOutFilename+=name;
 if(fIsMC) fOutFilename+=Form("MC.root");
 else fOutFilename+=("root");
 printf("filename is : %s \n",fOutFilename.Data());
 fOutFile = TFile::Open(fOutFilename.Data(),"recreate"); 

 Double_t mRange[2] = { 2.94, 3.1 };
 Int_t massBins = (Int_t) ((mRange[1] - mRange[0]) / 0.0025);

 //int arrSize = 5;
 TString selTitle[fgArrSize] = { "All Candidates (PID only)",  Form ("after K0 rejection (M within [%f,%f] GeV/#it{c}^{2})", param[kK0MassLow], param[kK0MassHigh]), Form ("after #Lambda rejection (M within [%1.3f,%1.3f] GeV/#it{c}^{2})", param[kLambdaMassLow], param[kLambdaMassHigh]),
  Form ("after #gamma rejection (M within [%2.2f,%2.2f] GeV/#it{c}^{2}",0., param[kGammaMassHigh]), " after K0, #Lambda, #gamma exclusion and M>3.015" };

 TString prodVtxTitle[5] = { "All Candidates (PID only)", Form ("K0-like (M within [%f,%f] GeV/#it{c}^{2})", param[kK0MassLow],param[kK0MassHigh]), Form ("#Lambda-like (M within [%1.3f,%1.3f] GeV/#it{c}^{2})", param[kK0MassLow], param[kK0MassHigh]), Form ("#gamma-like (M within [%2.2f,%2.2f] GeV/#it{c}^{2}",0., param[kGammaMassHigh]), " after K0, #Lambda, #gamma exclusion and M>3.015" };

 for (Int_t i = 0; i < fgArrSize; i++)
 {
  hMass[i] = new TH1F (Form ("hMass%i", i), Form ("%s", selTitle[i].Data ()), massBins, mRange[0], mRange[1]);
  hMass[i]->SetXTitle ("M(#pi^{3}H) GeV/c^{2}");
  hMass[i]-> SetYTitle (Form ("entries / %2.1f MeV/c^{2}", (hMass[i]->GetXaxis ()->GetBinWidth (2) * 1000)));
  hMass[i]->GetYaxis ()->SetTitleOffset (1.5);
  hMass[i]->Sumw2 ();
  hMass[i]->SetMarkerStyle (7);
  hMass[i]->SetMarkerColor (kBlue);
  hMass[i]->SetDrawOption ("err");
  fOutputList->AddLast(hMass[i]);

  hMassTrd[i] = new TH1F (Form ("hMassTrd%i", i), Form ("%s", selTitle[i].Data ()), massBins, mRange[0], mRange[1]);
  hMassTrd[i]->SetXTitle ("M(#pi^{3}H) GeV/c^{2}");
  hMassTrd[i]-> SetYTitle (Form ("entries / %2.1f MeV/c^{2}", (hMassTrd[i]->GetXaxis ()->GetBinWidth (2) * 1000)));
  hMassTrd[i]->GetYaxis ()->SetTitleOffset (1.5);
  hMassTrd[i]->Sumw2 ();
  hMassTrd[i]->SetMarkerStyle (7);
  hMassTrd[i]->SetMarkerColor (kGreen);
  hMassTrd[i]->SetDrawOption ("err");
  fOutputList->AddLast(hMassTrd[i]);
  

  hDecayLength[i] = new TH1F (Form ("hDecayLength%i", i), Form ("%s", selTitle[i].Data ()), 300, 0, 60);
  hDecayLength[i]->SetXTitle ("d (cm)");
  hDecayLength[i]->SetYTitle (Form("entries / %2.2f cm",(hDecayLength[i]->GetXaxis ()->GetBinWidth (2))));
  fOutputList->AddLast(hDecayLength[i]);

  hProdVtx[i] = new TH1F (Form("hProdVtx%i", i), Form ("%s", prodVtxTitle[i].Data ()), 300, 0, 60);
  hProdVtx[i]->SetXTitle ("decay length (XY) (cm)");
  hProdVtx[i]->SetYTitle (Form("entries / %2.2f cm", (hProdVtx[i]->GetXaxis ()->GetBinWidth (2))));
  fOutputList->AddLast(hProdVtx[i]);

  hDiffRelP[i] = new TH2F(Form("hDiffRelP%i", i), Form("P_triton-Ppion/PLnn %s", selTitle[i].Data ()), 240,0, 6, 300, 0, 1.5);
  hDiffRelP[i]->SetXTitle("V0 p (GeV/c)");
  hDiffRelP[i]->SetYTitle("%");
  fOutputList->AddLast(hDiffRelP[i]);
 }

 hDiffRelP[fgArrSize] = new TH2F (Form("hDiffRelP%i", fgArrSize), Form ("P_triton-Ppion/PLnn  all & M>3.015"), 450, 0, 15, 500,-10, 10);
 hDiffRelP[fgArrSize]->SetXTitle ("V0 p (GeV/c)");
 hDiffRelP[fgArrSize]->SetYTitle ("%");
 fOutputList->AddLast(hDiffRelP[fgArrSize]);


 if(fIsMC){
  hMassBkg = new TH1F("hMassBkg","", massBins, mRange[0], mRange[1]);
  fOutputList->AddLast(hMassBkg);
  hMassSignal = new TH1F("hMassSignal","", massBins, mRange[0], mRange[1]);
  fOutputList->AddLast(hMassSignal);
  TString quarkCont[7]={"gamma","d","u","s","c","b","other"}; 
  for(Int_t iq=0; iq<7; iq++){
   hMassContrib[0][iq] = new TH1F(Form("hMassContrPion_%s",quarkCont[iq].Data()),Form("contribution to the invariant mass if pion comes from %s ",quarkCont[iq].Data ()), massBins, mRange[0], mRange[1]);
   hMassContrib[0][iq]->SetXTitle("M(#pi^{3}H) GeV/c^{2}");  
   hMassContrib[0][iq]->SetYTitle(Form("entries / %2.1f MeV/c^{2}",(hMassContrib[0][iq]->GetXaxis()->GetBinWidth (2) * 1000)));  
   hMassContrib[0][iq]->SetLineColor(iq+1);
   fOutputList->AddLast(hMassContrib[0][iq]);


   hMassContrib[1][iq] = new TH1F(Form("hMassContr3H_%s",quarkCont[iq].Data()),Form("contribution to the invariant mass if 3H comes from %s",quarkCont[iq].Data ()), massBins, mRange[0], mRange[1]);
   hMassContrib[1][iq]->SetXTitle("M(#pi^{3}H) GeV/c^{2}");  
   hMassContrib[1][iq]->SetYTitle(Form("entries / %2.1f MeV/c^{2}",(hMassContrib[1][iq]->GetXaxis()->GetBinWidth (2) * 1000)));  
   hMassContrib[1][iq]->SetLineColor(iq+1);
   fOutputList->AddLast(hMassContrib[1][iq]);
  }
 }

 hMassK0 = new TH1F ("hMassK0", "", (0.6 - 0.4) / 0.002, 0.4, 0.6);
 hMassK0->SetXTitle ("M(#pi^{+}#pi^{-}) GeV/#it{c}^{2}");
 hMassK0->SetYTitle (Form("entries / %i MeV/#it{c}^{2}",(Int_t) (hMassK0->GetXaxis ()->GetBinWidth (2) * 1000)));
 hMassK0->SetFillColor (kBlue);
 fOutputList->AddLast(hMassK0);
 hMassGamma = new TH1F ("hMassGamma", "", 200, 0., 1.);
 hMassGamma->SetXTitle ("M(e^{+}e^{-}) GeV/#it{c}^{2}");
 hMassGamma->SetYTitle (Form("entries / %i MeV/#it{c}^{2}",(Int_t) (hMassGamma->GetXaxis ()->GetBinWidth (2) * 1000)));
 hMassGamma->SetFillColor (kBlue);
 fOutputList->AddLast(hMassGamma);
 hMassLambda = new TH1F ("hMassLambda", "", 800, 1., 2.);
 hMassLambda->SetXTitle ("M(#pip) GeV/#it{c}^{2}");
 hMassLambda->SetYTitle (Form("entries / %i MeV/#it{c}^{2}",(Int_t) (hMassLambda->GetXaxis ()->GetBinWidth (2) * 1000)));
 hMassLambda->SetFillColor (kBlue);
 fOutputList->AddLast(hMassLambda);
 hHighM = new TH1F ("hHighM", "V0 momentum (M>3.015)", 150, 0.5, 10);
 hHighM->SetXTitle ("p GeV/#it{c}");
 hHighM->SetYTitle (Form("entries / %2.2f GeV/#it{c}", (hHighM->GetXaxis ()->GetBinWidth (2) * 1000)));
 hHighM->SetFillColor (kBlue);
 fOutputList->AddLast(hMassLambda);

 for (Int_t k = 0; k < 5; k++)
 {
  hTriP[k] = new TH1F (Form ("hTriP%i", k), Form("^{3}H momentum distribution in - v0 %s", selTitle[k].Data ()), 400, 0, 8);
  hTriP[k]->SetLineColor (kRed);
  hTriP[k]->GetXaxis()->SetTitle("p GeV/c");
  hTriP[k]->GetYaxis()->SetTitle("entries");
  fOutputList->AddLast(hTriP[k]);

  hPiP[k] = new TH1F (Form ("hPiP%i", k), Form ("#pi momentum distribution in v0 - %s", selTitle[k].Data ()), 400, 0, 8);
  hPiP[k]->GetXaxis()->SetTitle("p GeV/c");
  hPiP[k]->GetYaxis()->SetTitle("entries");
  fOutputList->AddLast(hPiP[k]);

  hV0P[k] = new TH1F (Form ("hV0P%i", k),Form("v0 momentum distribution - %s", selTitle[k].Data ()),400, 0, 8);
  hV0P[k]->SetLineColor (kGreen);
  hV0P[k]->GetXaxis()->SetTitle("p GeV/c");
  hV0P[k]->GetYaxis()->SetTitle("entries");
  fOutputList->AddLast(hV0P[k]);
 }

 triTOFmass = new TH2F ("hTOFmass","", 100, 0,10, 500, 0, 10);
 fOutputList->AddLast(triTOFmass);
 hDcaD = new TH1F ("hDcaD", "dca duaghter tracks", 240, 0, 1.2);
 fOutputList->AddLast(hDcaD);

 Int_t nArmBins[2] = {110,120};
 Double_t armVar[2] = {1.1, 0.5};
 hArmPlot = new TH2F ("hArmPlot", "", nArmBins[0], -armVar[0], armVar[0], nArmBins[1], 0, armVar[1]);
 hArmPlot->SetXTitle ("#alpha");
 hArmPlot->SetYTitle ("p#_{T}^{Arm}");
 fOutputList->AddLast(hArmPlot);


 for (Int_t j = 0; j < fgArrSize; j++)
 {
  hArmPlotSel[j] = new TH2F (Form ("hArmPlotSel%i", j), Form ("%s", selTitle[j].Data ()), nArmBins[0],-armVar[0], armVar[0], nArmBins[1], 0, armVar[1]);
  hArmPlotSel[j]->SetXTitle ("#alpha");
  hArmPlotSel[j]->SetYTitle ("p_{T}^{Arm}");
  fOutputList->AddLast(hArmPlotSel[j]);
 }

 hTPCsignalPi = new TH2F ("hTPCsignalPi", "", 240, -0.6, 0.6, 800, 0, 400);
 fOutputList->AddLast(hTPCsignalPi);
 hTPCsignalTri = new TH2F ("hTPCsignalTri", "", 1600, -8, 8, 600, 0, 600);
 fOutputList->AddLast(hTPCsignalTri);
 hTPCsignalTriAll = new TH2F ("hTPCsignalTriAll", "", 800, 0, 8, 600, 0, 600);
 fOutputList->AddLast(hTPCsignalTriAll);
 hTPCsignalPiClean = new TH2F ("hTPCsignalPiClean", "", 200, -0.5, 0.5, 800, 0, 400);
 fOutputList->AddLast(hTPCsignalPiClean);
 hTPCsignalTriClean = new TH2F ("hTPCsignalTriClean", "", 1600, -8, 8, 600, 0, 600);
 fOutputList->AddLast(hTPCsignalTriClean);

 hTPCsignalTri91Lim = new TH2F ("hTPCsignalTri91Lim", "", 1600, -8, 8, 600, 0, 600);
 fOutputList->AddLast(hTPCsignalTri91Lim);
 hTPCsignalTriTrd = new TH2F ("hTPCsignalTriTrd", "",1600, -8, 8, 600, 0, 600);
 fOutputList->AddLast(hTPCsignalTriTrd);

 const Int_t nPossible = 34;

 Int_t partId[nPossible] = { 11,13, 22, 111, 113, 130, 211, 213, 223, 310, 313, 321, 323, 411, 421, 521,511, 221, 3122, 3222, 3112,3312,443,331,-2212,-2112,1114,2214,2224,2114,431,531,-1,-999};
 TString partNames[nPossible] = { "e", "mu", "gamma", "pi0","rho0", "K_L0", "pi", "rho_ch", "omega", "K_S0", "K*0", "K","K*", "D_ch", "D0", "B_ch", "B0", "eta", "Lambda0", "Sigma+","Sigma-","Xi-bar","J/psi","eta_prime","p_bar","n_bar","Delta-","Delta0","Delta+","Delta++","D_s","B_s0","primary","none" };
 TString pType[2] = { "pion", "triton" };
 for (Int_t im = 0; im < 2; im++)
 {
  hMumCheck[im] = new TH2F (Form ("hMumCheck_%s", pType[im].Data()), Form (" %s ", pType[im].Data()),nPossible, -0.5, nPossible-0.5, nPossible/2,-0.5, nPossible/2-0.5);
  TAxis *x = hMumCheck[im]->GetXaxis (); x->SetTitle ("mother");
  TAxis *y = hMumCheck[im]->GetYaxis (); y->SetTitle ("MC truth ");
  for(Int_t misId =0; misId<nPossible/2; misId++ ) y->SetBinLabel (misId + 1, partNames[misId].Data ());
  y->SetBinLabel(nPossible/2, "none");
  for (Int_t ids = 0; ids < nPossible; ids++) x->SetBinLabel (ids + 1, partNames[ids].Data ()); 
  fOutputList->AddLast(hMumCheck[im]);
 }

 //TH1F *hEv = new TH1F ("hEventMultiplicity", "Nb of Events", 12, -0.5, 11.5);
 TH1F *hEv = new TH1F ("hEventMultiplicity", "Nb of Events", 13, -0.5, 12.5);
 hEv->GetXaxis ()->SetBinLabel (1, "All Events");
 hEv->GetXaxis ()->SetBinLabel (2, "Events w/PV");
 hEv->GetXaxis ()->SetBinLabel (3, "Events w/|Vz|<10cm");
 hEv->GetXaxis ()->SetBinLabel (4, "Central Events");
 hEv->GetXaxis ()->SetBinLabel (5, "SemiCentral Events");
 hEv->GetXaxis ()->SetBinLabel (6, "MB Events");
 hEv->GetXaxis ()->SetBinLabel (7, "Central Events  w/|Vz|<10cm");
 hEv->GetXaxis ()->SetBinLabel (8, "SemiCentral Events  w/|Vz|<10cm");
 hEv->GetXaxis ()->SetBinLabel (9, "MB Events w/|Vz|<10cm");
 hEv->GetXaxis ()->SetBinLabel (10, "Any Events");
 hEv->GetXaxis ()->SetBinLabel (11, "Any Events w/|Vz|<10cm");
 hEv->GetXaxis ()->SetBinLabel (12, "other");
 fOutputList->AddLast(hEv);

 hMonitorPlot = new TH1I("hMonitorPlot","V0 selected",12,-0.5,11.5);
 hMonitorPlot->GetXaxis()->SetBinLabel (1, "All V0");
 hMonitorPlot->GetXaxis()->SetBinLabel(2,"#pi pT cut");
 hMonitorPlot->GetXaxis()->SetBinLabel(3,"triton PT cut");
 hMonitorPlot->GetXaxis()->SetBinLabel(4,"#pi n#sigma cut");
 hMonitorPlot->GetXaxis()->SetBinLabel(5,"triton n#sigma cut");
 hMonitorPlot->GetXaxis()->SetBinLabel(6,"V0 cuts");
 hMonitorPlot->GetXaxis()->SetBinLabel(7,"TOF p removal (<4 GeV)");
 hMonitorPlot->GetXaxis()->SetBinLabel(8,"TPC pid CUT(<91) -> TOF m2");
 hMonitorPlot->GetXaxis()->SetBinLabel(9,"TOF m2 cut");
 hMonitorPlot->GetXaxis()->SetBinLabel(10,"!TPC pid CUT (>91)");
  fOutputList->AddLast(hMonitorPlot);
 ntTot = new TNtupleD ("ntTot", "ntTot", "V0Mom:piMom:triMom:v0dca:dcaTriXY:dcaV0d:decayL:decayLxy:CosThetaP:invM:bkgMassType:isHele:alphaArm:ptArm:tpcSignal:TriMom:pionMom:TRDsignal:nITSc:nSigmaElTOF");
 fOutputList->AddLast(ntTot);

 return;
}

void Getoutput::LoopOverV0(Int_t Hcharge){
 printf("start looping over V0... (charge %i) \n",Hcharge);
 //Double_t tofMassUpperLimit = 7.929 - param[kParNsigmaTOFmass]*0.3634;	//see https://aliceinfo.cern.ch/Notes/node/471 
 //Double_t tofMassLimit = 7.929 + param[kParNsigmaTOFmass]*0.3634;	// pag 12 of 44 of 2016-Nov-28-analysis_note-alice-frontpage_analysis_notes-2.pdf
 Double_t tofMassUpperLimit = 7.929 - param[kParNsigmaTOFmass]*0.6;	//see https://aliceinfo.cern.ch/Notes/node/471 
 Double_t tofMassLimit = 7.929 + param[kParNsigmaTOFmass]*0.6;	// pag 12 of 44 of 2016-Nov-28-analysis_note-alice-frontpage_analysis_notes-2.pdf
 Double_t mass[2] = { 0., 0. };
 Int_t nv0 = 0;
 Bool_t isPrint = kFALSE;

 const Int_t nPossible = 34; // for MC study only
 Int_t partId[nPossible] = { 11,13, 22, 111, 113, 130, 211, 213, 223, 310, 313, 321, 323, 411, 421, 521,511, 221, 3122, 3222, 3112,3312,443,331,-2212,-2112,1114,2214,2224,2114,431,531,-1,-999};

 TNtupleD *hn = (TNtupleD *)(fInputList->At(fInputList->GetEntries() - 1));
 if(hn) printf ("ntuple entries %i ", (int) hn->GetEntries ());
 else {
  printf("No input Ntuple found, exiting... \n");
  return;
 }
 Double_t pionMom, tritonMom, v0Mom;
 Int_t countTri=0;
 Int_t countPi=0;

 for (Int_t iv0 = 0; iv0 < hn->GetEntries (); iv0++)
 {
  if (iv0 % 500000 == 0)
  {
   printf ("analizing entry %i\n", iv0);
   isPrint = kTRUE;
  }
  nv0++;

  hn->GetEntry (iv0);
  Double_t *arr = hn->GetArgs ();
//  printf("sign %i \n",(Int_t)arr[kSign]);
 //check on the V0 daughter sign
  if (fRejectLikesign)
  {
   if (TMath::Abs (arr[kSign]) == 11) continue; // reject +1(pion charge) + 10(triton charge x 10)  = 11 (or -1 -10 = -11) <=> same charges 
   if (Hcharge == 1 && arr[kSign] == -9) continue; // if one wants only + tritons => one needs 
   else if (Hcharge == -1 && arr[kSign] == 9) continue;
  } else {
   if (TMath::Abs (arr[kSign]) == 9) continue;
  //  printf("   Entry %i -> Sign %i = 3H %i \n",iv0,(Int_t)arr[kSign],(Int_t)Hcharge);
   if (Hcharge == 1 && arr[kSign] == -11) continue;
   else if (Hcharge == -1 && arr[kSign] == 11) continue;
  }
 //printf("ok -> Entry %i -> Sign %i = 3H %i \n",iv0,(Int_t)arr[kSign],(Int_t)Hcharge);

  TVector3 pPion (arr[0], arr[1], arr[2]);
  TVector3 pTriton (arr[3], arr[4], arr[5]);
  if(!fRejectLikesign && Hcharge==1) {
  // th4e way in which they are filled, needs a switch between pion and triton at this stage
  pPion.SetXYZ(arr[3],arr[4],arr[5]);
  pTriton.SetXYZ(arr[0], arr[1], arr[2]);
 }
  pionMom = pPion.Mag();
  tritonMom = pTriton.Mag();
  TVector3 sumV0=pPion+pTriton;
 
  //if(tritonMom>arr[kV0mom]) printf("weird triton %f pion %f and sum %f, whereas stored is %f  angle %f \n",tritonMom,pionMom,sumV0.Mag(),arr[kV0mom],pPion.Angle(pTriton)*TMath::RadToDeg()); 
  Double_t pionPt = TMath::Sqrt (arr[0] * arr[0] + arr[1] * arr[1]);
  Double_t tritonPt = TMath::Sqrt (arr[3] * arr[3] + arr[4] * arr[4]);
   if(!fRejectLikesign && Hcharge==1){
   pionPt = TMath::Sqrt (arr[3] * arr[3] + arr[4] * arr[4]);
   tritonPt = TMath::Sqrt (arr[0] * arr[0] + arr[1] * arr[1]);
   }
 //printf("ok -> Sign %i = 3H %i (v0 mom %f) pion pt %f, triton pt %f (pi Lim param %f  and parMinP %f) \n",(Int_t)arr[kSign],(Int_t)Hcharge,arr[kV0mom],pionPt,tritonPt,param[kParPiLim],param[kParMinP]);
  hMonitorPlot->Fill(0);
  if (pionPt > param[kParPiLim]) continue;
//printf("ok after pion Pt\n");
   //if (pionMom < param[kParMinP]) continue;
   if (pionMom < param[kParMinP]) continue;
//printf("ok after pion min P \n");
   hMonitorPlot->Fill(1);
  if (tritonPt < param[kPar3hLim]) continue;
  //if (tritonMom < param[kPar3hLim]) continue;
   if (tritonMom > param[kParMaxP3H])  continue;
   hMonitorPlot->Fill(2);
  Double_t nsigmapi = TMath::Abs (arr[kNSPi]);
  Double_t nsigmapiTof = TMath::Abs (arr[kSigPiFromPiTof]);
  Double_t nsigmatri = TMath::Abs (arr[kNSTri]);
  //printf("BEFORE nsigmapi %f , nsigmaTri %f, limit %f \n",nsigmapi,nsigmatri,param[kParNsigmaPID]);
  if (nsigmapi > param[kParNsigmaPID]) continue;
 hMonitorPlot->Fill(3);
  if(nsigmatri > param[kParNsigmaPID]) continue;
   if(arr[kTriTPCsignal]<grTriMs->Eval(tritonMom)) continue; 
   if(arr[kTriTPCsignal]>grTriPs->Eval(tritonMom)) continue; 
   if(arr[kPiTPCsignal]<grPiMs->Eval(pionMom)) continue; 
   if(arr[kPiTPCsignal]>grPiPs->Eval(pionMom)) continue; 

 hMonitorPlot->Fill(4);
 //if(nsigmatri>2.999) printf("AFTER : nsigmapi %f , nsigmaTri %f, limit %f \n",nsigmapi,nsigmatri,param[kParNsigmaPID]);
  if (fCaio){
  if(!EventSelectionAOD(arr)) continue;
  } else if(!EventSelectionSimple(arr)) continue;
  v0Mom = arr[kV0mom];
  hArmPlot->Fill (arr[kAlphaArm], arr[kPtArm]);
  hMonitorPlot->Fill(5);
 


  hDcaD->Fill (arr[kV0dcaD]);

  Int_t signPi = 999, signTri = 999;
  if (arr[kSign] == 9)       {signTri =  1; signPi = -1;}
  else if (arr[kSign] == -9) {signTri = -1; signPi =  1;}
  else if (arr[kSign] == 11) {signTri =  1; signPi =  1;}
  else if (arr[kSign] == -11){signTri = -1; signPi = -1;}
  else {printf ("sign value not known %i\n", (Int_t)arr[kSign]);}
  if (isPrint)
  {
   printf(" (triCharge %i) : 3H sign %i | pion sign %i (stored %i) \n", f3Hcharge, signTri, signPi, (Int_t) arr[kSign]);
   isPrint = kFALSE;
  }
  //if(arr[kTriTPCsignal]<-55*tritonMom+200) printf("nsigma triton %3.3f \n",arr[kNSTri]);
  hTPCsignalPi->Fill (pionMom * signPi, arr[kPiTPCsignal]);
  hTPCsignalTri->Fill (tritonMom * signTri, arr[kTriTPCsignal]);
  hTPCsignalTriAll->Fill (tritonMom, arr[kTriTPCsignal]);
  mass[0] = fgPionMass;
  mass[1] = fgLH3Mass;

  Double_t lambdaDecay[2] = {fgPionMass,fgProtMass };

	 if(arr[kIsTrdEle]>0.5)hTPCsignalTriTrd->Fill(tritonMom * signTri,arr[kTriTPCsignal]);
    if (fIncludePidTOF)
    {
     if (tritonMom < 4 && TMath::Abs (arr[kSigPrTof]) < 3)  continue; // reject 3H identified by TOF as p
      if (tritonMom < 4 && TMath::Abs (arr[kSigPiTof]) < 3)  continue; // reject 3H identified by TOF as pion
     hMonitorPlot->Fill(6);
    }
  
    if(fUseAODCut){
   if (arr[kTriTPCsignal] < 91) // Caio's way
   //if (arr[kTriTPCsignal] < 75) // 2 sigma 's way
   {
     hMonitorPlot->Fill(7);

    Double_t mTof2 = arr[kTriTOFmass] * arr[kTriTOFmass];
    triTOFmass->Fill(tritonMom,mTof2);
    if (mTof2 < tofMassUpperLimit) continue;
    if (mTof2 > tofMassLimit) continue;
    hTPCsignalTri91Lim->Fill(tritonMom * signTri,arr[kTriTPCsignal]);
     hMonitorPlot->Fill(8);
   }
    }

    else if(tritonMom>f3HPcut){
         hMonitorPlot->Fill(7);
    Double_t mTof2 = arr[kTriTOFmass] * arr[kTriTOFmass];
    triTOFmass->Fill(tritonMom,mTof2);
    if (mTof2 < tofMassUpperLimit) continue;
    if (mTof2 > tofMassLimit) continue;
    hTPCsignalTri91Lim->Fill(tritonMom * signTri,arr[kTriTPCsignal]);
     hMonitorPlot->Fill(8);

    }  
   
    else {
      hMonitorPlot->Fill(9);
    }
   // else continue;
  // inv Mass electrons 
  hTPCsignalPiClean->Fill (pionMom * signPi, arr[kPiTPCsignal]);
  hTPCsignalTriClean->Fill (tritonMom * signTri, arr[kTriTPCsignal]);
  Double_t gamma = GetInvMass (pPion, pTriton, fgEleMass, fgEleMass);
  hMassGamma->Fill (gamma);
  // inv Mass Lambda
  Double_t lambda =  GetInvMass (pPion, pTriton, lambdaDecay[0], lambdaDecay[1]);
  hMassLambda->Fill(lambda);
  // inv Mass K0
  Double_t k0 = GetInvMass (pPion, pTriton, fgPionMass, fgPionMass);
  hMassK0->Fill (k0);
  // inv Mass LNN
  Double_t LNN = GetInvMass (pPion, pTriton, mass[0], mass[1]);

  // MC cross checks
  if ((LNN > 2.98 && LNN < 3.0 ) && fIsMC) if (arr[kMumPiPdgCode] == arr[kMumTriPdgCode] && TMath::Abs (arr[kMumPiPdgCode]) != 1010000030)  printf ("UNCORRECT :  mother %i : pi pdg %i triton pdg %i \n", (Int_t) arr[kMumPiPdgCode], (Int_t) arr[kPiPdgCode], (Int_t) arr[kTriPdgCode]);
  Int_t selection = 0;
  // flag possible bkg sources
  if (k0 > param[kK0MassLow] && k0 < param[kK0MassHigh]) selection += 1;
  if (lambda > param[kLambdaMassLow] && lambda < param[kLambdaMassHigh]) selection += 10;
  if (gamma < param[kGammaMassHigh] && gamma >= 0) selection += 100;


  TDatabasePDG *d; 
  if(fIsMC) d =  TDatabasePDG::Instance ();
  Double_t TRDsig = arr[kTRDsig]; 
  Double_t tmvaArray[20] = { v0Mom, pionPt, tritonPt, arr[kV0Dca], arr[kDcaTriXY], arr[kV0dcaD], arr[kDecayPath], arr[kDecayPathXY], arr[kCosP], LNN, (Double_t)selection, arr[kIsTrdEle],arr[kAlphaArm], arr[kPtArm],arr[kTriTPCsignal],tritonMom,pionMom,TRDsig,arr[kNclusITS],arr[kSigElTof]};
  if(!fIsMC) {
   ntTot->Fill(tmvaArray); 
  }
  else {
   // left shift in the number is a bug
   printf("Mum Pi %i    Mum Tri  %i   - Pi %i  Tri %i \n", (Int_t)arr[kMumPiPdgCode],(Int_t)arr[kMumTriPdgCode],(Int_t)arr[kPiPdgCode],(Int_t)arr[kTriPdgCode]);
   if (arr[kMumPiPdgCode] == arr[kMumTriPdgCode] && TMath::Abs (arr[kMumPiPdgCode]) == 1010000030) {
    //tmvaArray[17]=arr[kMumTriPdgCode]; 
    ntTot->Fill (tmvaArray);      // true LNN inv Mass are in that interval 
    printf("analyzing MC, Filling the ntuple (please check that TRD is correctly stored in ntTot) \n");
   // 1010000030 - LNN
   // 1000010030 - 3H
   if(TMath::Abs (arr[kMumPiPdgCode]) == 1010000030) hMassSignal->Fill(LNN);
   }
   else {
    //printf("mother pi %i, mother triton %i \n",(Int_t) arr[kMumPiPdgCode],(Int_t) arr[kMumTriPdgCode]);  
    Int_t binPi[2] = { nPossible-1,nPossible/2}; //first is mother, second is true particle
    Int_t binTri[2] = { nPossible-1,nPossible/2};
    Double_t absMumTri = TMath::Abs(arr[kMumTriPdgCode]);
    Double_t absMumPi = TMath::Abs(arr[kMumPiPdgCode]);

    for (Int_t imo = 0; imo < nPossible; imo++)
    {
     Double_t id = TMath::Abs(partId[imo]);
     if (TMath::Abs (arr[kMumPiPdgCode]) == id) binPi[0] = imo;
     if (TMath::Abs (arr[kPiPdgCode]) == id) binPi[1] = imo;
     if (TMath::Abs (arr[kMumTriPdgCode]) == id) binTri[0] = imo;
     if (TMath::Abs (arr[kTriPdgCode]) == id) binTri[1] = imo;
    }

    if( (selection!=0 && selection <12) || selection==101 || selection==110) continue; // to select invM entries excluding kaons and lambdas
    //if(LNN>2.992 && LNN<2.999){
    hMumCheck[0]->Fill (binPi[0], binPi[1]);
    hMumCheck[1]->Fill (binTri[0], binTri[1]);
    //}
    hMassBkg->Fill(LNN);
    // triton case
    Int_t selBary3H = ((Int_t)absMumTri/1000);
    if(selBary3H>6) {
      printf("weird mother %i from %i, continuing\n",(Int_t)selBary3H,(Int_t)absMumTri);
     selBary3H=6;
    }
    if(selBary3H==0){

     Int_t selMes3H = ((Int_t)absMumTri/100);
     if(selMes3H>6) {
       printf("weird mother %i from %i, continuing\n",(Int_t)selMes3H,(Int_t)absMumTri);
      selMes3H=6;
     }
     if(selMes3H==0) {
      if(absMumTri==22) {
       hMassContrib[1][0]->Fill(LNN);
       //printf("selected quark content =  gamma (%i) \n",(Int_t)absMumTri);
      }
      else {
       hMassContrib[1][6]->Fill(LNN);
       countTri++;
      }
     }
     else {
      // printf("TRITON : mother %i \n",(Int_t)absMumTri);
      hMassContrib[1][selMes3H]->Fill(LNN);
      // printf("selected quark content :%i (meson %i)\n",selMes3H,(Int_t)absMumTri);
     }
    } else {
     hMassContrib[1][selBary3H]->Fill(LNN); 
     //printf("selected quark content :%i (baryon %i)\n",selBary3H,(Int_t)absMumTri);
    }

    // pion case
    Int_t selBaryPi = ((Int_t)absMumPi/1000);
    //printf("PION : mother %i, quark content %i \n",(Int_t)absMumPi,selBaryPi);
    if(selBaryPi>6) {
      printf("weird mother %i from %i, continuing\n",(Int_t)selBaryPi,(Int_t)absMumPi);
     selBaryPi=6;
    }
    if(selBaryPi==0){
     Int_t selMesPi = ((Int_t)absMumPi/100);
     if(selMesPi>6) {
       printf("weird mother %i from %i, continuing\n",(Int_t)selMesPi,(Int_t)absMumPi);
      selMesPi=6;
     }
     if(selMesPi==0) {
      if(absMumPi==22) hMassContrib[0][0]->Fill(LNN); 
      else {
       countPi++;
       hMassContrib[0][6]->Fill(LNN); 
      }
     }
     else hMassContrib[0][selMesPi]->Fill(LNN);
    } else {
     hMassContrib[0][selBaryPi]->Fill(LNN);
    }
   } // end analysis invM bkg contribution
  }
  // end MC checks


  hMass[0]->Fill(LNN);
  if(arr[kIsTrdEle] >0.5) hMassTrd[0]->Fill(LNN); // triton not misidentified
  hArmPlotSel[0]->Fill (arr[kAlphaArm], arr[kPtArm]);

  hDecayLength[0]->Fill (arr[kDecayPathXY]);

  hTriP[0]->Fill (tritonMom);
  hPiP[0]->Fill (pionMom);
  hV0P[0]->Fill (arr[kV0mom]);
  hProdVtx[0]->Fill (arr[kDecayPathXY]);
  hDiffRelP[0]->Fill(arr[kV0mom],(tritonMom)/arr[kV0mom]); 

  if (LNN > 3.015)
  {
   hProdVtx[4]->Fill (arr[kDecayPathXY]);
   hDecayLength[4]->Fill (arr[kDecayPathXY]);
   hPiP[4]->Fill (pionMom);
   hTriP[4]->Fill (tritonMom);
   hV0P[4]->Fill (arr[kV0mom]);
   hHighM->Fill (arr[kV0mom]);
 hDiffRelP[4]->Fill(arr[kV0mom],(tritonMom)/arr[kV0mom]); 
   
  }


  // K0, L and gamma removal
  if (k0 > param[kK0MassLow] && k0 < param[kK0MassHigh])
  {
   hProdVtx[1]->Fill (arr[kDecayPathXY]);
   //printf("k0 exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  hArmPlotSel[1]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hMass[1]->Fill (LNN);
  if(arr[kIsTrdEle] >0.5) hMassTrd[1]->Fill(LNN); // triton not misidentified
  hDecayLength[1]->Fill (arr[kDecayPathXY]);
  hTriP[1]->Fill (tritonMom);
  hPiP[1]->Fill (pionMom);
  hV0P[1]->Fill (arr[kV0mom]);
 hDiffRelP[1]->Fill(arr[kV0mom],(tritonMom)/arr[kV0mom]); 

  if (lambda > param[kLambdaMassLow]&& lambda < param[kLambdaMassHigh])
  {
   hProdVtx[2]->Fill (arr[kDecayPathXY]);
   //printf("lambda exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  hArmPlotSel[2]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hMass[2]->Fill (LNN);
  if(arr[kIsTrdEle] >0.5) hMassTrd[2]->Fill(LNN); // triton not misidentified
  hDecayLength[2]->Fill (arr[kDecayPathXY]);
  hTriP[2]->Fill (tritonMom);
  hPiP[2]->Fill (pionMom);
  hV0P[2]->Fill (arr[kV0mom]);
 hDiffRelP[2]->Fill(arr[kV0mom],(tritonMom)/arr[kV0mom]); 

  if (gamma < param[kGammaMassHigh] && gamma >= 0)
  {
   hProdVtx[3]->Fill (arr[kDecayPathXY]);
   //printf("gamma exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  //  if(arr[kPtArm]>0.12) continue; //rm pt arm to reduce the bkg
  // if(TMath::Abs(arr[20])<0.5) continue; //continue; rm pt arm to reduce the bkg
  hMass[3]->Fill (LNN);
  if(arr[kIsTrdEle] >0.5) hMassTrd[3]->Fill(LNN); // triton not misidentified
  hArmPlotSel[3]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hDecayLength[3]->Fill (arr[kDecayPathXY]);
  hTriP[3]->Fill (tritonMom);
  hPiP[3]->Fill (pionMom);
  hV0P[3]->Fill (arr[kV0mom]);
 hDiffRelP[3]->Fill(arr[kV0mom],(tritonMom)/arr[kV0mom]); 

   if(iv0 == (hn->GetEntries()-1))printf("finishing the loop\n");

  }

 printf("end looping over V0...\n");

 }

 void Getoutput::DrawResults(){
  TCanvas *c = new TCanvas ("cSum", "Summary Plots", 1700, 1000);
  c->Divide (4, 2);
  for (Int_t i = 0; i < 4; i++)
  {
  if(!hMass[i]){
  printf("no hMass, continue\n");
  continue;
  }
   c->cd (1 + 2*i);  if(hMass[i]) hMass[i]->DrawCopy();
   c->cd (2 + 2*i);  if(hArmPlotSel[i])hArmPlotSel[i]->DrawCopy ("colz");
  }
 TCanvas *cDiff = new TCanvas ("cDiff", "Difference", 1700, 1000);
  cDiff->Divide (2,2);
  for(Int_t iD=1; iD<5;iD++){
    cDiff->cd(iD);
   if(!hMass[iD]) {
     printf("no hmass histo, continuing");
     continue;
     }
    hMass[iD]->DrawCopy();
    hMassTrd[iD]->SetFillColor(kGreen);
     hMassTrd[iD]->SetFillStyle(3001);
    hMassTrd[iD]->DrawCopy("samehist");
  }

  
  TCanvas *cM = new TCanvas ();
  cM->Divide (3, 1);
  cM->cd (1); if(hMassK0) hMassK0->DrawCopy();
  cM->cd (2); if(hMassLambda) hMassLambda->DrawCopy();
  cM->cd (3); if(hMassGamma) hMassGamma->DrawCopy ();

  TCanvas *cP = new TCanvas();
  cP->Divide (2, 1);
  cP->cd (1);if(hMass[3]) hMass[3]->DrawCopy ();
  cP->cd (2);if(hArmPlotSel[3]) hArmPlotSel[3]->DrawCopy ("colz");

  TCanvas *cMom = new TCanvas();
  cMom->Divide (2,2);
  for (Int_t g = 0; g < 4; g++)
  {
   cMom->cd(g + 1); if(hPiP[g]) hPiP[g]->DrawCopy();
   if(hTriP[g]) hTriP[g]->DrawCopy ("same");
   if(hV0P[g])hV0P[g]->DrawCopy("same");
   
   TLegend *l = new TLegend(0.7,0.4, 0.8,0.7);
   l->AddEntry(hPiP[g],"#pi");
   l->AddEntry(hTriP[g],"^{3}H");
   l->AddEntry(hV0P[g],"V0");
   l->Draw("same");
  }

  if(fIsMC){
   TCanvas *cCheck = new TCanvas();
   cCheck->Divide(1,2);
   cCheck->cd (1);
   if(hMumCheck[0]){
    hMumCheck[0]->DrawCopy ("colz");
    hMumCheck[0]->DrawCopy ("sametext");
   }
   cCheck->cd (2);
   if(hMumCheck[1]){
    hMumCheck[1]->DrawCopy ("colz");
    hMumCheck[1]->DrawCopy ("sametext");
   }	
  }
 }

 bool Getoutput::LoadOutputData(const char *filename){
  printf("load output data \n"); 
  gStyle->SetOptStat(10);
  TFile *f = TFile::Open(filename);
  if(!f){
   printf(" file %s not found \n",filename);
   return kFALSE;
  }
  for(Int_t i=0; i<5; i++){
   hDiffRelP[i]=(TH2F*)f->Get(Form("hDiffRelP%i", i));
   hDecayLength[i]=(TH1F*)f->Get(Form ("hDecayLength%i", i));
   hProdVtx[i]=(TH1F*)f->Get(Form("hProdVtx%i",i));
   hMass[i]=(TH1F*)f->Get(Form("hMass%i", i));
   hMassTrd[i]=(TH1F*)f->Get(Form("hMassTrd%i", i));
   hTriP[i]=(TH1F*)f->Get(Form("hTriP%i", i));
   hPiP[i]=(TH1F*)f->Get(Form("hPiP%i", i));
   hV0P[i]=(TH1F*)f->Get(Form("hV0P%i", i));
   hArmPlotSel[i]=(TH2F*)f->Get(Form("hArmPlotSel%i", i));
  }

  hMassK0 = (TH1F*)f->Get("hMassK0");
  hMassLambda = (TH1F*)f->Get("hMassLambda");
  hMassGamma = (TH1F*)f->Get("hMassGamma");
  hHighM = (TH1F*)f->Get("hHighM");
  triTOFmass = (TH2F*)f->Get("triTOFmass");
  hDcaD = (TH1F*)f->Get("hDcaD");
  hArmPlot = (TH2F*)f->Get("hArmPlot");
  hTPCsignalPi = (TH2F*)f->Get("hTPCsignalPi");
  hTPCsignalTri = (TH2F*)f->Get("hTPCsignalTri");
  hTPCsignalTri91Lim = (TH2F*)f->Get("hTPCsignalTri91Lim");
  TH1F *hEv = (TH1F*)f->Get("hEventMultiplicity"); 
  if(fIsMC){
   hMumCheck[0]=(TH2F*)f->Get("hMumCheck_pion");
   hMumCheck[1]=(TH2F*)f->Get("hMumCheck_triton");
   TString n[7]={"gamma","d","u","s","c","b","other"};
   for(Int_t j=0; j<7; j++){
    hMassContrib[0][j]=(TH1F*)f->Get(Form("hMassContrPion_%s",n[j].Data())); 
    hMassContrib[1][j]=(TH1F*)f->Get(Form("hMassContr3H_%s",n[j].Data())); 
   }
  }
  // hMassBkg = (TH1F*)f->Get("hMassBkg");
  // hMassSignal = (TH1F*)f->Get("hMassSignal");
  //TH1F *hEv;
  printf("end load output data \n"); 
  return kTRUE;
 }

