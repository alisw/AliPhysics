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
  AliComparisonDEdx * compObj = (AliComparisonDEdx*)f.Get("AliComparisonDEdx");

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

#include <iostream>

#include "TFile.h"
#include "TCint.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TGraph.h"
//
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
//
#include "AliMathBase.h"
#include "AliTreeDraw.h"
//#include "TStatToolkit.h"

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonDEdx.h" 

using namespace std;

ClassImp(AliComparisonDEdx)

//_____________________________________________________________________________
AliComparisonDEdx::AliComparisonDEdx():
//  TNamed("AliComparisonDEdx","AliComparisonDEdx"),
  AliComparisonObject("AliComparisonDEdx"),

  // dEdx 
  fTPCSignalNormTan(0), 
  fTPCSignalNormSPhi(0),
  fTPCSignalNormTPhi(0), 
  //
  fTPCSignalNormTanSPhi(0),
  fTPCSignalNormTanTPhi(0),
  fTPCSignalNormTanSPt(0), 
  
  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),
  fMCPtMin(0),
  fMCAbsTanThetaMax(0),
  fMCPdgCode(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  Init();
}

//_____________________________________________________________________________
AliComparisonDEdx::~AliComparisonDEdx(){
   
  if(fTPCSignalNormTan)  delete fTPCSignalNormTan; fTPCSignalNormTan=0; 
  if(fTPCSignalNormSPhi) delete fTPCSignalNormSPhi; fTPCSignalNormSPhi=0;
  if(fTPCSignalNormTPhi) delete fTPCSignalNormTPhi; fTPCSignalNormTPhi=0;
  //
  if(fTPCSignalNormTanSPhi) delete fTPCSignalNormTanSPhi; fTPCSignalNormTanSPhi=0;
  if(fTPCSignalNormTanTPhi) delete fTPCSignalNormTanTPhi; fTPCSignalNormTanTPhi=0;
  if(fTPCSignalNormTanSPt)  delete fTPCSignalNormTanSPt; fTPCSignalNormTanSPt=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonDEdx::Init()
{
  // Init histograms
  
  // TPC dEdx
  fTPCSignalNormTan = new TH2F("CdEdxTan","CdEdxTan",50, -2,2,  40,30,70); 
  fTPCSignalNormTan->SetXTitle("tan(#theta)");
  fTPCSignalNormTan->SetYTitle("rec. dE/dx / calc. dE/dx");

  fTPCSignalNormSPhi   = new TH2F("CdEdxSPhi","CdEdxSPhi",10,0.0,1,40,30,70);
  fTPCSignalNormSPhi->SetXTitle("sin(#phi)");
  fTPCSignalNormSPhi->SetYTitle("rec. dE/dx / calc. dE/dx");

  fTPCSignalNormTPhi   = new TH2F("CdEdxTPhi","CdEdxTPhi",10,0.0,2,40,30,70); 
  fTPCSignalNormTPhi->SetXTitle("tan(#phi)");
  fTPCSignalNormTPhi->SetYTitle("rec. dE/dx / calc. dE/dx");

  fTPCSignalNormTanSPhi= new TH3F("CdEdxTanSPhi","CdEdxTanSPhi",20, -2,2, 10,0.0 ,1,  40,30,70);
  fTPCSignalNormTanSPhi->SetXTitle("tan(#theta)");
  fTPCSignalNormTanSPhi->SetYTitle("sin(#phi)");
  fTPCSignalNormTanSPhi->SetZTitle("rec. dE/dx / calc. dE/dx");

  fTPCSignalNormTanTPhi= new TH3F("CdEdxTanTPhi","CdEdxTanTPhi",20, -2,2, 10,0.0 ,1,  40,30,70);
  fTPCSignalNormTanTPhi->SetXTitle("tan(#theta)");
  fTPCSignalNormTanTPhi->SetYTitle("tan(#phi)");
  fTPCSignalNormTanTPhi->SetZTitle("rec. dE/dx / calc. dE/dx");

  fTPCSignalNormTanSPt= new TH3F("CdEdxTanSPt","CdEdxTanSPt",20, -2,2, 10,0.3 ,3, 40,30,70); 
  fTPCSignalNormTanSPt->SetXTitle("tan(#theta)");
  fTPCSignalNormTanSPt->SetYTitle("#sqrt{p_{t}}");
  fTPCSignalNormTanSPt->SetZTitle("rec. dE/dx / calc. dE/dx");

  // Init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

    // init folder
    fAnalysisFolder = CreateFolder("folderDEdx","Analysis de/dx Folder");
}

//_____________________________________________________________________________
void AliComparisonDEdx::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC){

  // Fill dE/dx  comparison information
  
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Float_t mprim = infoMC->GetPrim();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
  
  // Check selection cuts 
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;
  //if (mprim>1.4) return;
  //if (mprim<0.5) return;
  if (mprim > fCutsMC->GetMaxTPCSignal()) return;
  if (mprim < fCutsMC->GetMinTPCSignal()) return;
  if (infoRC->GetESDtrack()->GetTPCsignalN()<fCutsRC->GetMinTPCsignalN()) return;
  //
  Float_t ratio = infoRC->GetESDtrack()->GetTPCsignal()/infoMC->GetPrim();
  Float_t sphi =  infoRC->GetESDtrack()->GetInnerParam()->GetSnp();
  Float_t tphi =  sphi/TMath::Sqrt((1.-sphi)*(1.+sphi));

  if (TMath::Abs(infoMC->GetParticle().GetPdgCode()) != GetMCPdgCode()) return;
  //if (mcpt>0.5) {
  if (mcpt > GetMCPtMin()) {
    fTPCSignalNormTan->Fill(tantheta,ratio);    // only subset
  }

  //if (TMath::Abs(tantheta)<0.5){
  if (TMath::Abs(tantheta) < GetMCAbsTanThetaMax()){
    fTPCSignalNormSPhi->Fill(sphi,ratio);       // only subset
    fTPCSignalNormTPhi->Fill(tphi,ratio);       // only subset
  }
  fTPCSignalNormTanSPhi->Fill(tantheta,sphi,ratio);    
  fTPCSignalNormTanTPhi->Fill(tantheta,tphi,ratio);    
  fTPCSignalNormTanSPt->Fill(tantheta,TMath::Sqrt(mcpt),ratio);
}

//_____________________________________________________________________________
Long64_t AliComparisonDEdx::Merge(TCollection* list) 
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

    fTPCSignalNormTan->Add(entry->fTPCSignalNormTan);
    fTPCSignalNormSPhi->Add(entry->fTPCSignalNormSPhi);
    fTPCSignalNormTPhi->Add(entry->fTPCSignalNormTPhi);
    //
    fTPCSignalNormTanSPhi->Add(entry->fTPCSignalNormTanSPhi);
    fTPCSignalNormTanTPhi->Add(entry->fTPCSignalNormTanTPhi);
    fTPCSignalNormTanSPt->Add(entry->fTPCSignalNormTanSPt);

    count++;
  }

return count;
}

//_____________________________________________________________________________
void AliComparisonDEdx::Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Process comparison information
  Process(infoMC,infoRC);
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
  
  AliComparisonDEdx * comp=this;
  TObjArray *aFolderObj = new TObjArray;

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

