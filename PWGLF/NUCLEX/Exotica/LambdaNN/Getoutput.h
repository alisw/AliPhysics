#ifndef Getoutput_H
#define Getoutput_H

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
#include <TMath.h>

//using namespace std;

class Getoutput {
 public:
  Getoutput();
  bool LoadParams(const char *paramfile="params.txt");
  bool LoadFile(const char *inputfile="AnalysisResults.root",const char *listName="clistNtupHyper");
  void LoopOverV0(Int_t Hcharge=-1);
  void StoreOutputData(const char *filename="results.root");
  bool LoadOutputData(const char *filename="results3H1.root");
   void DrawResults();
  void ClearInputData();

  void BookOutputData();
  Double_t GetInvMass (TVector3 vPos, TVector3 vNeg, Double_t mPos, Double_t mNeg);
  bool EventSelectionAOD(Double_t *arr);

  Bool_t fIsMC;
  Bool_t fIncludePidTOF;
  Bool_t fRejectBkg; // useful in case like-sign V0 are produced
  Int_t f3Hsign;
  Double_t array[28];
  Float_t param[20];

  TList *fInputList;
  TList *fHistList;
  TList *fOutputList;
  //ClassDef(Getoutput,0);

  TH2F *hDiffRelP[6];
  TH1F *hDecayLength[5];
  TH1F *hProdVtx[5];
  TH1F *hMass[5];
  TH1F *hMassBkg;
  TH1F *hMassSignal;
  TH1F *hMassContrib[2][7]; // quark u&d, s,c,b for pion and triton
  TH1F *hMassK0;
  TH1F *hMassLambda;
  TH1F *hMassGamma;
  TH1F *hHighM;
  TH1F *hTriP[5];
  TH1F *hPiP[5];
  TH1F *hV0P[5];
  TH2F *triTOFmass;
  TH1F *hDcaD;
  TH2F *hArmPlot;
  TH2F *hArmPlotSel[5];
  TH2F *hTPCsignalPi;
  TH2F *hTPCsignalTri;
  TH2F* hTPCsignalTri91Lim;
  TH2F *hMumCheck[2];
  TH1F *hEv;

  TNtupleD *ntTot;

  static const Int_t fgArrSize;
  static Double_t fgLH3Mass;
  static Double_t fgPionMass;
  static Double_t fgEleMass;
  static Double_t fgProtMass;

};


enum { kPposx, kPposy, kPposz, kPnegx, kPnegy, kPnegz,	//0-5
 kNSPi, kNSTri, kTriTOFmass, kPiTPCsignal, kTriTPCsignal,	// 6-10
 kV0mom, kPtArm, kAlphaArm,	// 11-13
 kDcaTriXY, kDcaTriZ, kV0dcaD, kDecayPath, kDecayPathXY,	// 14-18
 kV0Dca, kCosP, kV0VtxErrSum, kSign,	// 19-22
 kDcaPi, kIsTrdEle, kSigPiFromPiTof, kSigPrTof, kSigPiTof, kNclusITS,
 kPiPdgCode, kTriPdgCode, kMumPiPdgCode, kMumTriPdgCode
};				//23-27

enum {
 kParMinP, kParMinPv0, kParMaxP3H, kParPiLim, kPar3hLim, kParNclusITS, kParNsigmaPID, kParNsigmaTOFmass, kParDcaTriZ, kParCosP, kParV0Dca,
 kK0MassLow, kK0MassHigh, kLambdaMassLow, kLambdaMassHigh,kGammaMassHigh
};
#endif

/*
   Double_t GetInvMass (TVector3 vPos, TVector3 vNeg, Double_t mPos, Double_t mNeg);
   Bool_t eventSelectionCaio (Double_t * array);

   void Getoutut(Int_t triCharge = -1, Int_t runMin = 999, Bool_t bkgRej = kTRUE)
   {

   rejectBkg = bkgRej;

   TString outFilename;
   if (rejectBkg) outFilename = Form ("plots.opposite.tritSign%i.root",triCharge);
   else outFilename = Form ("plots.same.tritSign%i.root",triCharge);
   TFile *fout = TFile::Open (outFilename.Data (), "recreate");
   TNtupleD *ntTot = new TNtupleD ("ntTot", "ntTot", "V0Mom:piMom:triMom:v0dca:dcaTriXY:dcaV0d:decayL:decayLxy:CosP:SigIIvtx:invM:selMass:pRat");
   Double_t tofMassUpperLimit = 7.929 - param[kParNsigmaTOFmass]*0.03634;	//see https://aliceinfo.cern.ch/Notes/node/471 
   Double_t tofMassLimit = 7.929 + param[kParNsigmaTOFmass]*0.03634;	// pag 12 of 44 of 2016-Nov-28-analysis_note-alice-frontpage_analysis_notes-2.pdf

//
// histogram 
//

Double_t mRange[2] = { 2.94, 3.1 };
Int_t massBins = (Int_t) ((mRange[1] - mRange[0]) / 0.0025);
Double_t k0Interval[2] = { 0.485, 0.515 };
Double_t lambdaInterval[2] = { 1.1, 1.13 };
Double_t gammaInterval[2] = { 0.0, 0.05 };

const int arrSize = 5;
TString selTitle[5] = { "All Candidates (PID only)",  Form ("after K0 rejection (M within [%f,%f] GeV/#it{c}^{2})", k0Interval[0], k0Interval[1]), Form ("after #Lambda rejection (M within [%1.3f,%1.3f] GeV/#it{c}^{2})", lambdaInterval[0], lambdaInterval[1]),
Form ("after #gamma rejection (M within [%2.2f,%2.2f] GeV/#it{c}^{2}", gammaInterval[0], gammaInterval[1]), " after K0, #Lambda, #gamma exclusion and M>3.015" };

TString prodVtxTitle[5] = { "prodVtx if M> 3.015", Form ("K0-like (M within [%f,%f] GeV/#it{c}^{2})", k0Interval[0], k0Interval[1]), Form ("#Lambda-like (M within [%1.3f,%1.3f] GeV/#it{c}^{2})", lambdaInterval[0], lambdaInterval[1]), Form ("#gamma-like (M within [%2.2f,%2.2f] GeV/#it{c}^{2}",gammaInterval[0], gammaInterval[1]), " after K0, #Lambda, #gamma exclusion and M>3.015" };

TH2F *hDiffRelP[arrSize + 1];
TH1F *hMass[arrSize], *hDecayLength[arrSize], *hProdVtx[arrSize];
for (Int_t i = 0; i < arrSize; i++)
{
hMass[i] = new TH1F (Form ("hMass%i", i), Form ("%s", selTitle[i].Data ()), massBins, mRange[0], mRange[1]);
hMass[i]->SetXTitle ("M(#pi^{3}H) GeV/c^{2}");
hMass[i]-> SetYTitle (Form ("entries / %2.1f MeV/c^{2}", (hMass[i]->GetXaxis ()->GetBinWidth (2) * 1000)));
hMass[i]->GetYaxis ()->SetTitleOffset (1.5);
hMass[i]->Sumw2 ();
hMass[i]->SetMarkerStyle (7);
hMass[i]->SetMarkerColor (kBlue);
hMass[i]->SetDrawOption ("err");

hDecayLength[i] = new TH1F (Form ("hDecayLength%i", i), Form ("%s", selTitle[i].Data ()), 300, 0, 60);
hDecayLength[i]->SetXTitle ("d (cm)");
hDecayLength[i]->SetYTitle (Form("entries / %2.2f cm",(hDecayLength[i]->GetXaxis ()->GetBinWidth (2))));

hProdVtx[i] = new TH1F (Form ("hProdVtx%i", i), Form ("%s", prodVtxTitle[i].Data ()), 300, 0, 60);
hProdVtx[i]->SetXTitle ("decay length (XY) (cm)");
hProdVtx[i]->SetYTitle (Form ("entries / %2.2f cm", (hProdVtx[i]->GetXaxis ()->GetBinWidth (2))));

hDiffRelP[i] = new TH2F (Form ("hDiffRelP%i", i), Form("P_triton-Ppion/PLnn %s", selTitle[i].Data ()), 240,0, 6, 300, 0, 1.5);
hDiffRelP[i]->SetXTitle ("V0 p (GeV/c)");
hDiffRelP[i]->SetYTitle ("%");
}


hDiffRelP[arrSize] = new TH2F (Form (Form ("hDiffRelP%i", arrSize)), Form ("P_triton-Ppion/PLnn  all & M>3.015"), 450, 0, 15, 500,-10, 10);
hDiffRelP[arrSize]->SetXTitle ("V0 p (GeV/c)");
hDiffRelP[arrSize]->SetYTitle ("%");
TH1F *hMassBkg = new TH1F("hMassBkg","", massBins, mRange[0], mRange[1]);
TH1F *hMassSignal = new TH1F("hMassSignal","", massBins, mRange[0], mRange[1]);
TH1F *hMassContrib[2][7]; // quark u&d, s,c,b for pion and triton
TString quarkCont[7]={"gamma","d","u","s","c","b","other"}; 
for(Int_t iq=0; iq<7; iq++){
hMassContrib[0][iq] = new TH1F(Form("hMassContrPion_%s",quarkCont[iq].Data()),Form("contribution to the invariant mass if pion comes from %s ",quarkCont[iq].Data ()), massBins, mRange[0], mRange[1]);
hMassContrib[0][iq]->SetXTitle("M(#pi^{3}H) GeV/c^{2}");  
hMassContrib[0][iq]->SetYTitle(Form("entries / %2.1f MeV/c^{2}",(hMassContrib[0][iq]->GetXaxis()->GetBinWidth (2) * 1000)));  
hMassContrib[0][iq]->SetLineColor(iq+1);

hMassContrib[1][iq] = new TH1F(Form("hMassContr3H_%s",quarkCont[iq].Data()),Form("contribution to the invariant mass if 3H comes from %s",quarkCont[iq].Data ()), massBins, mRange[0], mRange[1]);
hMassContrib[1][iq]->SetXTitle("M(#pi^{3}H) GeV/c^{2}");  
hMassContrib[1][iq]->SetYTitle(Form("entries / %2.1f MeV/c^{2}",(hMassContrib[1][iq]->GetXaxis()->GetBinWidth (2) * 1000)));  
hMassContrib[1][iq]->SetLineColor(iq+1);
}


TH1F *hMassK0 = new TH1F ("hMassK0", "", (0.6 - 0.4) / 0.002, 0.4, 0.6);
hMassK0->SetXTitle ("M(#pi^+#pi^-) GeV/#it{c}^{2}");
hMassK0->
SetYTitle (Form
  ("entries / %i MeV/#it{c}^{2}",
   (Int_t) (hMassK0->GetXaxis ()->GetBinWidth (2) * 1000)));
hMassK0->SetFillColor (kBlue);
TH1F *hMassGamma = new TH1F ("hMassGamma", "", 200, 0., 1.);
hMassGamma->SetXTitle ("M(e^{+}e^{-}) GeV/#it{c}^{2}");
hMassGamma->
SetYTitle (Form
  ("entries / %i MeV/#it{c}^{2}",
   (Int_t) (hMassGamma->GetXaxis ()->GetBinWidth (2) * 1000)));
hMassGamma->SetFillColor (kBlue);
TH1F *hMassLambda = new TH1F ("hMassLambda", "", 800, 1., 2.);
hMassLambda->SetXTitle ("M(#pip) GeV/#it{c}^{2}");
hMassLambda->SetYTitle (Form("entries / %i MeV/#it{c}^{2}",(Int_t) (hMassLambda->GetXaxis ()->GetBinWidth (2) * 1000)));
hMassLambda->SetFillColor (kBlue);
TH1F *hHighM = new TH1F ("hHighM", "V0 momentum (M>3.015)", 150, 0.5, 10);
hHighM->SetXTitle ("p GeV/#it{c}");
hHighM->SetYTitle (Form("entries / %2.2f GeV/#it{c}", (Int_t) (hHighM->GetXaxis ()->GetBinWidth (2) * 1000)));
hHighM->SetFillColor (kBlue);

TH1F *hTriP[5], *hPiP[5], *hV0P[5];
for (Int_t k = 0; k < 5; k++)
{
 hTriP[k] = new TH1F (Form ("hTriP%i", k), Form("^{3}H momentum distribution in - v0 %s", selTitle[k].Data ()), 400, 0, 8);
 hTriP[k]->SetLineColor (kRed);
 hPiP[k] = new TH1F (Form ("hPiP%i", k), Form ("#pi momentum distribution in v0 - %s", selTitle[k].Data ()), 400, 0, 8);
 hV0P[k] = new TH1F (Form ("hV0P%i", k),Form("v0 momentum distribution - %s", selTitle[k].Data ()),400, 0, 8);
 hV0P[k]->SetLineColor (kGreen);
 //TH1F *hPiN = new TH1F("hPiN","#pi momentum distribution in v0",400,0,8);
 //TH1F *hV0P = new TH1F("hV0P","V0 momentum distribution",400,0,8);
}
TH2F *triTOFmass = new TH2F ("hTOFmass","#Delta m (m_{TOF}-m_{3H}) with TPC ID ^3H TOF mass", 100, 0,10, 100, -5, 5);
TH1F *hDcaD = new TH1F ("hDcaD", "dca duaghter tracks", 240, 0, 1.2);

Int_t nArmBins[2] = { 110, 120 };
Double_t armVar[2] = { 1.1, 0.5 };
TH2F *hArmPlot = new TH2F ("hArmPlot", "", nArmBins[0], -armVar[0], armVar[0], nArmBins[1], 0, armVar[1]);
hArmPlot->SetXTitle ("#alpha");
hArmPlot->SetYTitle ("p#_{T}^{Arm}");
TH2F *hArmPlotSel[arrSize];
for (Int_t j = 0; j < arrSize; j++)
{
 hArmPlotSel[j] = new TH2F (Form ("hArmPlotSel%i", j), Form ("%s", selTitle[j].Data ()), nArmBins[0],-armVar[0], armVar[0], nArmBins[1], 0, armVar[1]);
 hArmPlotSel[j]->SetXTitle ("#alpha");
 hArmPlotSel[j]->SetYTitle ("p_{T}^{Arm}");
}

TH2F *hTPCsignalPi = new TH2F ("hTPCsignalPi", "", 600, -6, 6, 600, 0, 600);
TH2F *hTPCsignalTri = new TH2F ("hTPCsignalTri", "", 600, -6, 6, 600, 0, 600);
TH2F *hTPCsignalTri91Lim = new TH2F ("hTPCsignalTri91Lim", "", 600, -6, 6, 600, 0, 600);

const Int_t nPossible = 34;

Int_t partId[nPossible] = { 11,13, 22, 111, 113, 130, 211, 213, 223, 310, 313, 321, 323, 411, 421, 521,511, 221, 3122, 3222, 3112,3312,443,331,-2212,-2112,1114,2214,2224,2114,431,531,-1,-999};
TString partNames[nPossible] = { "e", "mu", "gamma", "pi0","rho0", "K_L0", "pi", "rho_ch", "omega", "K_S0", "K*0", "K","K*", "D_ch", "D0", "B_ch", "B0", "eta", "Lambda0", "Sigma+","Sigma-","Xi-bar","J/psi","eta_prime","p_bar","n_bar","Delta-","Delta0","Delta+","Delta++","D_s","B_s0","primary","none" };
TH2F *hMumCheck[2];
TString pType[2] = { "pion", "triton" };
for (Int_t im = 0; im < 2; im++)
{
 hMumCheck[im] = new TH2F (Form ("hMumCheck_%s", pType[im].Data()), Form (" %s ", pType[im].Data()),nPossible, -0.5, nPossible-0.5, nPossible/2,-0.5, nPossible/2-0.5);
 TAxis *x = hMumCheck[im]->GetXaxis (); x->SetTitle ("mother");
 TAxis *y = hMumCheck[im]->GetYaxis (); y->SetTitle ("MC truth ");
 for(Int_t misId =0; misId<nPossible/2; misId++ ) y->SetBinLabel (misId + 1, partNames[misId].Data ());
 y->SetBinLabel(nPossible/2, "none");
 for (Int_t ids = 0; ids < nPossible; ids++) x->SetBinLabel (ids + 1, partNames[ids].Data ()); 
}


Double_t LH3Mass = 2.808921;
Double_t PionMass = 0.13957;
Double_t eleMass = 0.00051;
Double_t ProtMass = 0.93827;

TH1F *hEv = new TH1F ("hEventMultiplicity", "Nb of Events", 12, -0.5, 11.5);
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


Double_t mass[2] = { 0., 0. };
Int_t nv0 = 0;


TObjArray *files = new TObjArray ();

//files->AddLast(new TObjString("LNN.Ntuple.MC.2.root"));
// files->AddLast(new TObjString("second.root"));

//files->AddLast(new TObjString("testLHC14a1b.2.root"));
//files->AddLast(new TObjString("LNN.Ntuple.MC.LHC14a6.root")); // LNN MC production

//files->AddLast (new TObjString ("LNN.Ntuple.MC.LHC14a1b.root"));

files->AddLast (new TObjString ("LNN.Ntuple.LHC12a17a_fix.0.root"));
files->AddLast (new TObjString ("LNN.Ntuple.LHC12a17a_fix.1.root"));
files->AddLast (new TObjString ("LNN.Ntuple.LHC12a17a_fix.2.root"));
files->AddLast (new TObjString ("LNN.Ntuple.LHC12a17a_fix.3.root"));
files->AddLast (new TObjString ("LNN.Ntuple.LHC12a17a_fix.4.root"));

Int_t countTri=0;
Int_t countPi=0;



Bool_t isPrint = kFALSE;
for (Int_t iFile = 0; iFile < files->GetEntries (); iFile++)
{
 TString name = ((TObjString *) (files->At (iFile)))->GetString ();
 if (!gSystem->IsFileInIncludePath (name.Data ())) continue;
 printf ("Processing file %s\n", name.Data ());
 TFile *f = TFile::Open (name.Data ());
 if (!f)
 {
  printf ("file %s not available \n", name.Data ());
  continue;
 }

 l = (TList *) f->Get ("LNNlist");
 if (!l)
 {
  printf ("-------> !!! no list in this file !!!\n");
  continue;
 }
 TH1F *h = (TH1F *) l->At (0);
 hEv->Add (h);

 TNtupleD *hn = (TNtupleD *) (l->At (l->GetEntries () - 1));
 printf ("ntuple entries %i ", (int) hn->GetEntries ());
 Double_t pionMom, tritonMom, v0Mom;

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
  TVector3 pPion (arr[0], arr[1], arr[2]);
  TVector3 pTriton (arr[3], arr[4], arr[5]);
  pionMom = pPion.Mag ();
  tritonMom = pTriton.Mag ();
  Double_t pionPt = TMath::Sqrt (arr[0] * arr[0] + arr[1] * arr[1]);
  Double_t tritonPt = TMath::Sqrt (arr[3] * arr[3] + arr[4] * arr[4]);
  if (pionPt > param[kParPiLim]) continue;
  if (tritonPt < param[kPar3hLim]) continue;

  if (!eventSelectionCaio (arr)) continue;
  v0Mom = arr[kV0mom];
  hArmPlot->Fill (arr[kAlphaArm], arr[kPtArm]);

  //check on the V0 daughter sign
  if (rejectBkg)
  {
   if (TMath::Abs (arr[kSign]) == 11) continue; 
   if (triCharge == 1 && arr[kSign] == -9) continue;
   else if (triCharge == -1 && arr[kSign] == 9) continue;
  } else {
   if (TMath::Abs (arr[kSign]) == 9) continue;
   if (triCharge == 1 && arr[kSign] == -11) continue;
   else if (triCharge == -1 && arr[kSign] == 11) continue;
  }


  Double_t nsigmapi = TMath::Abs (arr[kNSPi]);
  Double_t nsigmapiTof = TMath::Abs (arr[kSigPiFromPiTof]);
  Double_t nsigmatri = TMath::Abs (arr[kNSTri]);

  hDcaD->Fill (arr[kV0dcaD]);
  if (pionPt < param[kParMinP]) continue;
  if (tritonPt > param[kParMaxP3H])  continue;

  Int_t signPi = 999, signTri = 999;
  if (arr[kSign] == 9)
  {
   signTri = 1;
   signPi = -1;
  }
  else if (arr[kSign] == -9)
  {
   signTri = -1;
   signPi = 1;
  }
  else if (arr[kSign] == 11)
  {
   signTri = 1;
   signPi = 1;
  }
  else if (arr[kSign] == -11)
  {
   signTri = -1;
   signPi = -1;
  }
  else
  {
   printf ("sign value not known %i\n", arr[kSign]);
  }
  if (isPrint)
  {
   printf
    (" (triCharge %i) : 3H sign %i | pion sign %i (stored %i) \n", triCharge, signTri, signPi, (Int_t) arr[kSign]);
   isPrint = kFALSE;
  }

  hTPCsignalPi->Fill (pionMom * signPi, arr[kPiTPCsignal]);
  hTPCsignalTri->Fill (tritonMom * signTri, arr[kTriTPCsignal]);

  mass[0] = PionMass;
  mass[1] = LH3Mass;

  Double_t lambdaDecay[2] = {PionMass,ProtMass };


  if (nsigmapi < param[kParNSPi] && nsigmatri < param[kParNSPi])
  {
   if (arr[kTriTPCsignal] < 91) // Caio's way
   {
    if (includePidTof)
    {
     if (tritonMom < 4 && TMath::Abs (arr[kSigPrTof]) < 3)  continue; // reject 3H identified by TOF as p 
    }
    Double_t mTof2 = arr[kTriTOFmass] * arr[kTriTOFmass];
    if (mTof2 < tofMassUpperLimit) continue;
    if (mTof2 > tofMassLimit) continue;
    hTPCsignalTri91Lim->Fill(tritonMom * signTri,arr[kTriTPCsignal]);
   }

  }

  else continue;
  // inv Mass electrons 
  Double_t gamma = GetInvMass (pPion, pTriton, eleMass, eleMass);
  hMassGamma->Fill (gamma);

  // inv Mass Lambda
  Double_t lambda =
   GetInvMass (pPion, pTriton, lambdaDecay[0], lambdaDecay[1]);
  hMassLambda->Fill (lambda);

  // inv Mass K0
  Double_t k0 = GetInvMass (pPion, pTriton, PionMass, PionMass);
  hMassK0->Fill (k0);
  // inv Mass LNN
  Double_t LNN = GetInvMass (pPion, pTriton, mass[0], mass[1]);



  if (LNN > 2.98 && LNN < 3.0)
   if (arr[kMumPiPdgCode] == arr[kMumTriPdgCode]
     && TMath::Abs (arr[kMumPiPdgCode]) != 1010000030)
    printf ("UNCORRECT :  mother %i : pi pdg %i triton pdg %i \n",
      (Int_t) arr[kMumPiPdgCode], (Int_t) arr[kPiPdgCode],
      (Int_t) arr[kTriPdgCode]);
  Int_t selection = 0;
  if (k0 > k0Interval[0] && k0 < k0Interval[1]) selection += 1;
  if (lambda > lambdaInterval[0] && lambda < lambdaInterval[1]) selection += 10;
  if (gamma < gammaInterval[1] && gamma > gammaInterval[0]) selection += 100;

  TDatabasePDG *d = TDatabasePDG::Instance ();

  Double_t tmvaArray[13] = { v0Mom, pionMom, tritonMom, arr[kV0Dca], arr[kDcaTriXY],
   arr[kV0dcaD], arr[kDecayPath], arr[kDecayPathXY], arr[kCosP], arr[kV0VtxErrSum], LNN,
   selection, 0. };
  if (arr[kMumPiPdgCode] == arr[kMumTriPdgCode] && TMath::Abs (arr[kMumPiPdgCode]) == 1010000030) ntTot->Fill (tmvaArray);	// trule LNN inv Mass are in that interval 
  // 1010000030 - LNN
  // 1000010030 - 3H
  if(TMath::Abs (arr[kMumPiPdgCode]) == 1010000030) hMassSignal->Fill(LNN);
  else {
   //printf("mother pi %i, mother triton %i \n",(Int_t) arr[kMumPiPdgCode],(Int_t) arr[kMumTriPdgCode]);  
   Int_t binPi[2] = { nPossible-1,nPossible/2};	//first is mother, second is true particle
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
    printf("weird mother %i from %i, continuing\n",selBary3H,absMumTri);
    selBary3H=6;
   }
   //Int_t selMes3H = ((Int_t)absMumTri/100);
   // if(selBary3H!=0) printf("TRITON : mother %i \n",(Int_t)absMumTri);
   if(selBary3H==0){

    Int_t selMes3H = ((Int_t)absMumTri/100);
    if(selMes3H>6) {
     printf("weird mother %i from %i, continuing\n",selMes3H,absMumTri);
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
    printf("weird mother %i from %i, continuing\n",selBaryPi,absMumPi);
    selBaryPi=6;
   }

   if(selBaryPi==0){
    Int_t selMesPi = ((Int_t)absMumPi/100);
    if(selMesPi>6) {
     printf("weird mother %i from %i, continuing\n",selMesPi,absMumPi);
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



  }

  hMass[0]->Fill (LNN);
  hArmPlotSel[0]->Fill (arr[kAlphaArm], arr[kPtArm]);

  hDecayLength[0]->Fill (arr[kDecayPathXY]);

  hTriP[0]->Fill (tritonMom);
  //if(pionMom > 0.35) printf("weird!! \n"); 
  hPiP[0]->Fill (pionMom);
  hV0P[0]->Fill (arr[kV0mom]);

  if (LNN > 3.015)
  {
   hProdVtx[4]->Fill (arr[kDecayPathXY]);
   hDecayLength[4]->Fill (arr[kDecayPathXY]);
   hPiP[4]->Fill (pionMom);
   hTriP[4]->Fill (tritonMom);
   hV0P[4]->Fill (arr[kV0mom]);
  }

  // K0, L and gamma removal

  if (k0 > k0Interval[0] && k0 < k0Interval[1])
  {
   hProdVtx[1]->Fill (arr[kDecayPathXY]);
   printf("k0 exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  hArmPlotSel[1]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hMass[1]->Fill (LNN);
  hDecayLength[1]->Fill (arr[kDecayPathXY]);
  hTriP[1]->Fill (tritonMom);
  hPiP[1]->Fill (pionMom);
  hV0P[1]->Fill (arr[kV0mom]);


  if (lambda > lambdaInterval[0] && lambda < lambdaInterval[1])
  {
   hProdVtx[2]->Fill (arr[kDecayPathXY]);
   printf("lambda exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  hArmPlotSel[2]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hMass[2]->Fill (LNN);
  hDecayLength[2]->Fill (arr[kDecayPathXY]);
  hTriP[2]->Fill (tritonMom);
  hPiP[2]->Fill (pionMom);
  hV0P[2]->Fill (arr[kV0mom]);


  if (gamma < gammaInterval[1] && gamma > gammaInterval[0])
  {
   hProdVtx[3]->Fill (arr[kDecayPathXY]);
   //printf("gamma exclusion, selection = %i \n",(Int_t)selection);
   continue;
  }
  //  if(arr[kPtArm]>0.12) continue; //rm pt arm to reduce the bkg
  // if(TMath::Abs(arr[20])<0.5) continue; //continue; rm pt arm to reduce the bkg
  hMass[3]->Fill (LNN);
  hArmPlotSel[3]->Fill (arr[kAlphaArm], arr[kPtArm]);
  hDecayLength[3]->Fill (arr[kDecayPathXY]);
  hTriP[3]->Fill (tritonMom);
  hPiP[3]->Fill (pionMom);
  hV0P[3]->Fill (arr[kV0mom]);

  if (LNN > 3.015)
  {
   //if(LNN<2.98) {
   hProdVtx[4]->Fill (arr[kDecayPathXY]);
   hHighM->Fill (arr[kV0mom]);
   hTriP[4]->Fill (tritonMom);
   hPiP[4]->Fill (pionMom);
   hV0P[4]->Fill (arr[kV0mom]);
   hDecayLength[4]->Fill (arr[kDecayPathXY]);
  }

  }

  f->Close ();
 }
 //hProdVtx[0]->SetTitle("After K0, #Lambda and #gamma rejection and M > 3.02");
 TCanvas *c = new TCanvas ("cSum", "Summary Plots", 1700, 1000);
 c->Divide (4, 2);
 for (Int_t i = 0; i < 4; i++)
 {
  c->cd (1 + 2 * i);
  hMass[i]->DrawCopy ();
  c->cd (2 + 2 * i);
  hArmPlotSel[i]->DrawCopy ("colz");
 }

 TCanvas *cM = new TCanvas ();
 cM->Divide (3, 1);
 cM->cd (1);
 hMassK0->DrawCopy ();
 cM->cd (2);
 hMassLambda->DrawCopy ();
 cM->cd (3);
 hMassGamma->DrawCopy ();

 TCanvas *cP = new TCanvas ();
 cP->Divide (2, 1);
 cP->cd (1);
 hMass[3]->DrawCopy ();
 cP->cd (2);
 hArmPlotSel[3]->DrawCopy ("colz");

 TCanvas *cMom = new TCanvas ();
 cMom->Divide (2, 3);
 for (Int_t g = 0; g < 5; g++)
 {
  cMom->cd (g + 2);
  hPiP[g]->DrawCopy ();
  hTriP[g]->DrawCopy ("same");
  hV0P[g]->DrawCopy ("same");
 }

 TCanvas *cComp = new TCanvas ();
 cComp->Divide (3, 1);
 cComp->cd (1);
 hDecayLength[3]->DrawCopy ();
 hProdVtx[4]->DrawCopy ("same")->SetLineColor (kRed);
 cComp->cd (2);
 hDecayLength[3]->Add (hProdVtx[4], -1);
 hDecayLength[3]->DrawCopy ();
 cComp->cd (3);
 hHighM->DrawCopy ();

 TCanvas *cCheck = new TCanvas ();
 cCheck->Divide (1, 2);
 cCheck->cd (1);
 hMumCheck[0]->DrawCopy ("colz");
 hMumCheck[0]->DrawCopy ("sametext");
 cCheck->cd (2);
 hMumCheck[1]->DrawCopy ("colz");
 hMumCheck[1]->DrawCopy ("sametext");
 printf("missing triton %i, missing pion %i \n",countTri,countPi);
 fout->Write ();
 fout->Close ();

}

*/
