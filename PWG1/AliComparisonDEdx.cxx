//------------------------------------------------------------------------------
// Implementation of AliComparisonDEdx class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data of AliComparisonDEdx.
//  
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  //AliComparisonDEdx * compObj = (AliComparisonDEdx*)f.Get("AliComparisonDEdx");
  AliComparisonDEdx * compObj = (AliComparisonDEdx*)cOutput->FindObject("AliComparisonDEdx");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderDEdx" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_DEdx.root"."recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include <TAxis.h>
#include <TCanvas.h>
#include <TH1.h>

#include "AliComparisonDEdx.h" 
#include "AliESDEvent.h"
#include "AliESDRecInfo.h" 
#include "AliLog.h" 
#include "AliMCInfo.h" 
#include "AliMCInfoCuts.h" 
#include "AliMathBase.h"
#include "AliRecInfoCuts.h" 
#include "AliTreeDraw.h"

using namespace std;

ClassImp(AliComparisonDEdx)

//_____________________________________________________________________________
AliComparisonDEdx::AliComparisonDEdx():
//  TNamed("AliComparisonDEdx","AliComparisonDEdx"),
  AliComparisonObject("AliComparisonDEdx"),

  // dEdx 
  fDeDxHisto(0),
  
  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // default constructor	
}

//_____________________________________________________________________________
AliComparisonDEdx::AliComparisonDEdx(Char_t* name="AliComparisonDEdx", Char_t* title="AliComparisonDEdx",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliComparisonObject(name,title),

  // dEdx 
  fDeDxHisto(0),
  
  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // named constructor

  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);
  Init();
}


//_____________________________________________________________________________
AliComparisonDEdx::~AliComparisonDEdx()
{
  // destructor
  if(fDeDxHisto)  delete fDeDxHisto; fDeDxHisto=0; 

  if(fCutsRC) delete fCutsRC; fCutsRC=0;
  if(fCutsMC) delete fCutsMC; fCutsMC=0;
  
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonDEdx::Init()
{
  // Init histograms
  
  // TPC dEdx
  Int_t nPBins = 31;
    Double_t binsP[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
    Double_t pMin = 0., pMax = 10.;

    if(IsHptGenerator() == kTRUE) {
      nPBins = 100;
      pMin = 0.; pMax = 100.;
    }

   //signal:alpha:y:z:snp:tgl:ncls:pid:p
   Int_t binsQA[9]    = {600,50, 50,  50, 50, 50, 80, 5, nPBins};
   Double_t xminQA[9] = {0,  -4,-20,-250, -1, -2, 0, 0., pMin};
   Double_t xmaxQA[9] = {300, 4, 20, 250,  1,  2, 160, 5., pMax};

   fDeDxHisto = new THnSparseF("fDeDxHisto","signal:alpha:y:z:snp:tgl:ncls:pid:momentum",9,binsQA,xminQA,xmaxQA);
   if(!IsHptGenerator()) fDeDxHisto->SetBinEdges(8,binsP);

   fDeDxHisto->GetAxis(0)->SetTitle("signal");
   fDeDxHisto->GetAxis(1)->SetTitle("alpha (rad)");
   fDeDxHisto->GetAxis(2)->SetTitle("y (cm)");
   fDeDxHisto->GetAxis(3)->SetTitle("z (cm)");
   fDeDxHisto->GetAxis(4)->SetTitle("snp");
   fDeDxHisto->GetAxis(5)->SetTitle("tgl");
   fDeDxHisto->GetAxis(6)->SetTitle("ncls");
   fDeDxHisto->GetAxis(6)->SetTitle("pid");
   fDeDxHisto->GetAxis(7)->SetTitle("p (GeV/c)");
   fDeDxHisto->Sumw2();

   // Init cuts
   if(!fCutsMC) 
     AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
   if(!fCutsRC) 
     AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

   // init folder
   fAnalysisFolder = CreateFolder("folderDEdx","Analysis de/dx Folder");
}

//_____________________________________________________________________________
void AliComparisonDEdx::ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Fill dE/dx  comparison information
  
  // Check selection cuts 
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
 
  Double_t pid = -1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetEM() ) pid = 0;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetMuM() ) pid = 1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP() ) pid = 2;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP() ) pid = 3;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt() ) pid = 4;

  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  Float_t dedx = infoRC->GetESDtrack()->GetTPCsignal();
  Int_t ncls = infoRC->GetESDtrack()->GetTPCNcls();

  const AliExternalTrackParam *innerParam =  0;
  if ((innerParam = infoRC->GetESDtrack()->GetInnerParam()) == 0) return;

  Double_t pt = innerParam->Pt();
  Double_t lam = TMath::ATan2(innerParam->Pz(),innerParam->Pt());
  Double_t p = pt/TMath::Cos(lam);
  Double_t alpha = innerParam->GetAlpha();
  Double_t y = innerParam->GetY();
  Double_t z = innerParam->GetZ();
  Double_t snp = innerParam->GetSnp();
  Double_t tgl = innerParam->GetTgl();

  Double_t vDeDxHisto[9] = {dedx,alpha,y,z,snp,tgl,ncls,pid,p};
  fDeDxHisto->Fill(vDeDxHisto); 
}

//_____________________________________________________________________________
void AliComparisonDEdx::ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Fill dE/dx  comparison information
  
   AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
void AliComparisonDEdx::ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Fill dE/dx  comparison information
  
   AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
Long64_t AliComparisonDEdx::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliComparisonDEdx* entry = dynamic_cast<AliComparisonDEdx*>(obj);
    if (entry == 0) continue;

    fDeDxHisto->Add(entry->fDeDxHisto);
    count++;
  }

return count;
}

//_____________________________________________________________________________
void AliComparisonDEdx::Exec(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Process comparison information

  if(GetAnalysisMode() == 0) ProcessTPC(infoMC,infoRC);
  else if(GetAnalysisMode() == 1) ProcessTPCITS(infoMC,infoRC);
  else if(GetAnalysisMode() == 2) ProcessConstrained(infoMC,infoRC);
  else {
    printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
    return;
  }
}

//_____________________________________________________________________________
TH1F* AliComparisonDEdx::MakeResol(TH2F * his, Int_t integ, Bool_t type)
{
  // Make resolution histograms
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoI(his,&hism,integ);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliComparisonDEdx::Analyse()
{
  // Analyze comparison information and store output histograms
  // in the folder "folderDEdx"
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;

  /*
  TH1F *hiss=0;
  TGraph2D * gr=0;

  // write results in the folder 
  TCanvas * c = new TCanvas("can","TPC dedx");
  c->cd();

  hiss = comp->MakeResol(comp->fTPCSignalNormTan,4,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma_{dEdx}");
  hiss->Draw();
  hiss->SetName("TPCdEdxResolTan");

  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fTPCSignalNormTan,4,1); 
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("<dEdx>");
  hiss->Draw(); 
  hiss->SetName("TPCdEdxMeanTan");

  aFolderObj->Add(hiss);
  //
  gr = AliMathBase::MakeStat2D(comp->fTPCSignalNormTanSPt,3,1,4);
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr->GetZaxis()->SetTitle("<dEdx>");
  gr->SetName("TPCdEdxMeanTanPt_1");
  gr->GetHistogram()->Draw("colz"); 

  gr->GetHistogram()->SetName("TPCdEdxMeanTanPt_1");
  aFolderObj->Add(gr->GetHistogram());
  //
  gr = AliMathBase::MakeStat2D(comp->fTPCSignalNormTanSPt,3,1,5);
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr->GetZaxis()->SetTitle("#sigma_{dEdx}");
  gr->SetName("TPCdEdxMeanTanPt_2");
  gr->GetHistogram()->Draw("colz"); 

  gr->GetHistogram()->SetName("TPCdEdxMeanTanPt_2");
  aFolderObj->Add(gr->GetHistogram());
  */

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjrArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonDEdx::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonDEdx * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}


//_____________________________________________________________________________
TFolder* AliComparisonDEdx::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}

