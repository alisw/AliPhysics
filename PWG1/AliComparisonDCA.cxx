//------------------------------------------------------------------------------
// Implementation of AliComparisonDCA class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder
// which is a data member of AliComparisonDCA.
//  
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  AliComparisonDCA * compObj = (AliComparisonDCA*)f.Get("AliComparisonDCA");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderDCA" 
  compObj->GetAnalysisFolder()->ls("*");
 
  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_DCA.root","recreate");
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
#include "AliTracker.h"   
#include "AliESDEvent.h"   
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
#include "AliESDVertex.h" 
//
#include "AliMathBase.h"
#include "AliTreeDraw.h" 

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonDCA.h" 

using namespace std;

ClassImp(AliComparisonDCA)

//_____________________________________________________________________________
AliComparisonDCA::AliComparisonDCA():
//  TNamed("AliComparisonDCA","AliComparisonDCA"),
  AliComparisonObject("AliComparisonDCA"),

  // DCA histograms
  fD0TanSPtB1(0),
  fD1TanSPtB1(0),
  fD0TanSPtL1(0),
  fD1TanSPtL1(0),
  fD0TanSPtInTPC(0),
  fD1TanSPtInTPC(0),
  fVertex(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  Init();
}

//_____________________________________________________________________________
AliComparisonDCA::~AliComparisonDCA()
{
  //
  if(fD0TanSPtB1) delete fD0TanSPtB1; fD0TanSPtB1=0;
  if(fD1TanSPtB1) delete fD1TanSPtB1; fD1TanSPtB1=0;
  if(fD0TanSPtL1) delete fD0TanSPtL1; fD0TanSPtL1=0;
  if(fD1TanSPtL1) delete fD1TanSPtL1; fD1TanSPtL1=0;
  if(fD0TanSPtInTPC) delete fD0TanSPtInTPC; fD0TanSPtInTPC=0;
  if(fD1TanSPtInTPC) delete fD1TanSPtInTPC; fD1TanSPtInTPC=0;
  if(fVertex) delete fVertex; fVertex=0;
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;

}

//_____________________________________________________________________________
void AliComparisonDCA::Init()
{
  // DCA histograms
  fD0TanSPtB1 = new TH3F("DCAyTanSPtB1","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPtB1->SetXTitle("tan(#theta)");
  fD0TanSPtB1->SetYTitle("#sqrt{p_{t}(GeV/c)}");
  fD0TanSPtB1->SetZTitle("DCA_{xy}");

  fD1TanSPtB1 = new TH3F("DCAzTanSPtB1","DCAzTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD1TanSPtB1->SetXTitle("tan(#theta)");
  fD1TanSPtB1->SetYTitle("#sqrt(p_{t}(GeV/c))");
  fD1TanSPtB1->SetZTitle("DCA_{z}");

  fD0TanSPtL1 = new TH3F("DCAyTanSPtL1","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPtL1->SetXTitle("tan(#theta)");
  fD0TanSPtL1->SetYTitle("#sqrt{p_{t}(GeV/c)}");
  fD0TanSPtL1->SetZTitle("DCA_{xy}");

  fD1TanSPtL1 = new TH3F("DCAzTanSPtL1","DCAzTanSPt",40,-2,2, 10,0.3,3, 100, -1,1);
  fD1TanSPtL1->SetXTitle("tan(#theta)");
  fD1TanSPtL1->SetYTitle("#sqrt{p_{t}(GeV/c)}");
  fD1TanSPtL1->SetZTitle("DCA_{z}");

  fD0TanSPtInTPC = new TH3F("DCAyTanSPtInTPC","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPtInTPC->SetXTitle("tan(#theta)");
  fD0TanSPtInTPC->SetYTitle("#sqrt{p_{t}(GeV/c)}");
  fD0TanSPtInTPC->SetZTitle("DCA_{xy}");

  fD1TanSPtInTPC = new TH3F("DCAzTanSPtInTPC","DCAzTanSPt",40,-2,2, 10,0.3,3, 100, -1,1);
  fD1TanSPtInTPC->SetXTitle("tan(#theta)");
  fD1TanSPtInTPC->SetYTitle("#sqrt{p_{t}(GeV/c)}");
  fD1TanSPtInTPC->SetZTitle("DCA_{z}");

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
 
  // init folder
  fAnalysisFolder = CreateFolder("folderDCA","Analysis DCA Folder");

  // vertex (0,0,0)
  fVertex = new AliESDVertex();
  fVertex->SetXv(0.0);
  fVertex->SetYv(0.0);
  fVertex->SetZv(0.0);
}

//_____________________________________________________________________________
void AliComparisonDCA::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill DCA comparison information
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Float_t spt = TMath::Sqrt(mcpt);

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

  // calculate and set prim. vertex
  fVertex->SetXv( infoMC->GetParticle().Vx() - dv[0] );
  fVertex->SetYv( infoMC->GetParticle().Vy() - dv[1] );
  fVertex->SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // calculate track parameters at vertex
  if (infoRC->GetESDtrack()->GetTPCInnerParam())
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bDCAStatus = track->PropagateToDCA(fVertex,field,kMaxD,dca,cov);

      if(bDCAStatus) {
        fD0TanSPtInTPC->Fill(tantheta,spt,dca[0]);
        fD1TanSPtInTPC->Fill(tantheta,spt,dca[1]);
	  }

	  delete track;
    }
  }

 if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)==0){
    fD0TanSPtB1->Fill(tantheta,spt,dca[0]);
    fD1TanSPtB1->Fill(tantheta,spt,dca[1]);
  }
    fD0TanSPtL1->Fill(tantheta,spt,dca[0]);
    fD1TanSPtL1->Fill(tantheta,spt,dca[1]);
}

//_____________________________________________________________________________
Long64_t AliComparisonDCA::Merge(TCollection* list) 
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
    AliComparisonDCA* entry = dynamic_cast<AliComparisonDCA*>(obj);
    if (entry == 0) continue; 
    

    fD0TanSPtB1->Add(entry->fD0TanSPtB1);
    fD1TanSPtB1->Add(entry->fD1TanSPtB1);
    fD0TanSPtL1->Add(entry->fD0TanSPtL1);
    fD1TanSPtL1->Add(entry->fD1TanSPtL1);
    fD0TanSPtInTPC->Add(entry->fD0TanSPtInTPC);
    fD1TanSPtInTPC->Add(entry->fD1TanSPtInTPC);

    count++;
  }

return count;
}

//_____________________________________________________________________________
void AliComparisonDCA::Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  // Process comparison information
  Process(infoMC,infoRC);
}

//_____________________________________________________________________________
void AliComparisonDCA::Analyse()
{
  //
  // Analyse comparison information and store output histograms
  // in the analysis folder "folderDCA" 
  //
  
  TH1::AddDirectory(kFALSE);

  TGraph2D *gr=0;
  TGraph * gr0=0;
  AliComparisonDCA * comp=this;
  TObjArray *aFolderObj = new TObjArray;

  // write results in the folder 
  // Canvas to draw analysed histograms
  TCanvas * c = new TCanvas("canDCA","DCA resolution");
  c->Divide(1,2);
  c->cd(1);
  //
  // DCA resolution
  //
  gr0 = AliMathBase::MakeStat1D(comp->fD0TanSPtB1,2,5);
  gr0->GetXaxis()->SetTitle("Tan(#theta)");
  gr0->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gr0->SetName("DCAResolTan");
  gr0->Draw("Al*");

  //if(folder) folder->Add(gr0->GetHistogram());
  aFolderObj->Add(gr0);
  //
  c->cd(2);
  gr = AliMathBase::MakeStat2D(comp->fD0TanSPtB1,4,2,5); 
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV/c)}");
  gr->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gr->SetName("DCAResolSPTTan");
  gr->GetHistogram()->Draw("colz");

  aFolderObj->Add(gr->GetHistogram());

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonDCA::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonDCA * comp=this;
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
TFolder* AliComparisonDCA::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
