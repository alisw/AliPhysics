//------------------------------------------------------------------------------
// Implementation of AliPerformanceDCA class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder
// which is a data member of AliPerformanceDCA.
//  
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  AliPerformanceDCA * compObj = (AliPerformanceDCA*)coutput->FindObject("AliPerformanceDCA");

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

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "AliPerformanceDCA.h" 
#include "AliESDEvent.h"   
#include "AliESDVertex.h" 
#include "AliLog.h" 
#include "AliMathBase.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliStack.h" 
#include "AliMCEvent.h" 
#include "AliTracker.h"   
#include "AliHeader.h"   
#include "AliGenEventHeader.h"   

using namespace std;

ClassImp(AliPerformanceDCA)

//_____________________________________________________________________________
AliPerformanceDCA::AliPerformanceDCA():
  AliPerformanceObject("AliPerformanceDCA"),

  // DCA histograms
  fDCAHisto(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  // default constructor	
  Init();
}

//_____________________________________________________________________________
AliPerformanceDCA::AliPerformanceDCA(Char_t* name="AliPerformanceDCA", Char_t* title="AliPerformanceDCA",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),

  // DCA histograms
  fDCAHisto(0),

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
AliPerformanceDCA::~AliPerformanceDCA()
{
  // destructor
  if(fDCAHisto)  delete fDCAHisto; fDCAHisto=0; 
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceDCA::Init()
{
  // DCA histograms
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 10.;

  Double_t *binsPt = 0;
  if (IsHptGenerator())  { 
    nPtBins = 100; ptMax = 100.;
    binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
  } else {
    binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
  }

  /*
  Int_t nPtBins = 31;
   Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
   Double_t ptMin = 0., ptMax = 10.;

   if(IsHptGenerator() == kTRUE) {
     nPtBins = 100;
     ptMin = 0.; ptMax = 100.;
   }
   */

   //dca_r, dca_z, eta, pt
   Int_t binsQA[4]    = {100,100,30,nPtBins};
   Double_t xminQA[4] = {-10.,-10.,-1.5,ptMin};
   Double_t xmaxQA[4] = {10.,10.,1.5,ptMax};

   fDCAHisto = new THnSparseF("fDCAHisto","dca_r:dca_z:eta:pt",4,binsQA,xminQA,xmaxQA);
   fDCAHisto->SetBinEdges(3,binsPt);

   fDCAHisto->GetAxis(0)->SetTitle("dca_r (cm)");
   fDCAHisto->GetAxis(1)->SetTitle("dca_z (cm)");
   fDCAHisto->GetAxis(2)->SetTitle("#eta");
   fDCAHisto->GetAxis(3)->SetTitle("p_{T} (GeV/c)");
   fDCAHisto->Sumw2();

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
 
  // init folder
  fAnalysisFolder = CreateFolder("folderDCA","Analysis DCA Folder");
}

//_____________________________________________________________________________
void AliPerformanceDCA::ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack)
{
  // Fill DCA comparison information
  if(!esdTrack) return;

  const AliExternalTrackParam *track = esdTrack->GetTPCInnerParam();
  if(!track) return;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

  if (esdTrack->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return; // min. nb. TPC clusters  
 
  Double_t vDCAHisto[4]={dca[0],dca[1],track->Eta(),track->Pt()};
  fDCAHisto->Fill(vDCAHisto);

  //
  // Fill rec vs MC information
  //
  if(!stack) return;

 }

//_____________________________________________________________________________
void AliPerformanceDCA::ProcessTPCITS(AliStack* const stack, AliESDtrack *const esdTrack)
{
  // Fill DCA comparison information
  if(!esdTrack) return;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParameters(dca,cov);

  if ((esdTrack->GetStatus()&AliESDtrack::kTPCrefit)==0) return; // TPC refit
  if (esdTrack->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return; // min. nb. TPC clusters  
  Int_t clusterITS[200];
  if(esdTrack->GetITSclusters(clusterITS)<fCutsRC->GetMinNClustersITS()) return;  // min. nb. ITS clusters

  Double_t vDCAHisto[4]={dca[0],dca[1],esdTrack->Eta(),esdTrack->Pt()};
  fDCAHisto->Fill(vDCAHisto);

  //
  // Fill rec vs MC information
  //
  if(!stack) return;

}

void AliPerformanceDCA::ProcessConstrained(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill DCA comparison information
  
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
Long64_t AliPerformanceDCA::Merge(TCollection* const list) 
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
    AliPerformanceDCA* entry = dynamic_cast<AliPerformanceDCA*>(obj);
    if (entry == 0) continue; 

    fDCAHisto->Add(entry->fDCAHisto);
    count++;
  }

return count;
}

//_____________________________________________________________________________
void AliPerformanceDCA::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, const Bool_t bUseMC)
{
  // Process comparison information
  if(!esdEvent) 
  {
      AliDebug(AliLog::kError, "esdEvent not available");
      return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  if(bUseMC)
  {
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

  } // end bUseMC

  //  Process events
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if(!track) continue;

    if(GetAnalysisMode() == 0) ProcessTPC(stack,track);
    else if(GetAnalysisMode() == 1) ProcessTPCITS(stack,track);
    else if(GetAnalysisMode() == 2) ProcessConstrained(stack,track);
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }
}

//_____________________________________________________________________________
void AliPerformanceDCA::Analyse()
{
  //
  // Analyse comparison information and store output histograms
  // in the analysis folder "folderDCA" 
  //
  
  TH1::AddDirectory(kFALSE);
  TH1 *h1D=0;
  TH2 *h2D=0;
  //TH3 *h3D=0;
  TObjArray *aFolderObj = new TObjArray;
  char title[256];

  // set pt measurable range 
  fDCAHisto->GetAxis(3)->SetRangeUser(0.10,10.);

  //
  h2D = fDCAHisto->Projection(0,1); // inverse projection convention
  h2D->SetName("dca_r_vs_dca_z");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(1)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(0)->GetTitle());
  sprintf(title,"%s vs %s",fDCAHisto->GetAxis(0)->GetTitle(),fDCAHisto->GetAxis(1)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  //
  h2D = fDCAHisto->Projection(0,2);
  h2D->SetName("dca_r_vs_eta");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(0)->GetTitle());
  sprintf(title,"%s vs %s",fDCAHisto->GetAxis(2)->GetTitle(),fDCAHisto->GetAxis(0)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  h1D = MakeStat1D(h2D,0,0);
  h1D->SetName("mean_dca_r_vs_eta");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h1D->GetYaxis()->SetTitle("mean_dca_r (cm)");
  sprintf(title," mean_dca_r (cm) vs %s",fDCAHisto->GetAxis(2)->GetTitle());
  h1D->SetTitle(title);
  aFolderObj->Add(h1D);

  h1D = MakeStat1D(h2D,0,1);
  h1D->SetName("rms_dca_r_vs_eta");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h1D->GetYaxis()->SetTitle("rms_dca_r (cm)");
  sprintf(title," rms_dca_r (cm) vs %s",fDCAHisto->GetAxis(2)->GetTitle());
  h1D->SetTitle(title);
  aFolderObj->Add(h1D);

  //
  h2D = fDCAHisto->Projection(0,3);
  h2D->SetName("dca_r_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(0)->GetTitle());
  sprintf(title,"%s vs %s",fDCAHisto->GetAxis(0)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  h2D->SetBit(TH1::kLogX);
  aFolderObj->Add(h2D);

  h1D = MakeStat1D(h2D,0,0);
  h1D->SetName("mean_dca_r_vs_pt");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h1D->GetYaxis()->SetTitle("mean_dca_r (cm)");
  sprintf(title,"mean_dca_r (cm) vs %s",fDCAHisto->GetAxis(3)->GetTitle());
  h1D->SetTitle(title);
  h1D->SetBit(TH1::kLogX);
  aFolderObj->Add(h1D);

  h1D = MakeStat1D(h2D,0,1);
  h1D->SetName("rms_dca_r_vs_pt");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h1D->GetYaxis()->SetTitle("rms_dca_r (cm)");
  sprintf(title,"rms_dca_r (cm) vs %s",fDCAHisto->GetAxis(3)->GetTitle());
  h1D->SetTitle(title);
  h1D->SetBit(TH1::kLogX);
  aFolderObj->Add(h1D);
   
  // 
  h2D = fDCAHisto->Projection(1,2);
  h2D->SetName("dca_z_vs_eta");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(1)->GetTitle());
  sprintf(title,"%s vs %s",fDCAHisto->GetAxis(1)->GetTitle(),fDCAHisto->GetAxis(2)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  h1D = MakeStat1D(h2D,0,0);
  h1D->SetName("mean_dca_z_vs_eta");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h1D->GetYaxis()->SetTitle("mean_dca_z (cm)");
  sprintf(title,"mean_dca_z (cm) vs %s",fDCAHisto->GetAxis(2)->GetTitle());
  h1D->SetTitle(title);
  aFolderObj->Add(h1D);

  h1D = MakeStat1D(h2D,0,1);
  h1D->SetName("rms_dca_z_vs_eta");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h1D->GetYaxis()->SetTitle("rms_dca_z (cm)");
  sprintf(title,"rms_dca_z (cm) vs %s",fDCAHisto->GetAxis(2)->GetTitle());
  h1D->SetTitle(title);
  aFolderObj->Add(h1D);

  //
  h2D = fDCAHisto->Projection(1,3);
  h2D->SetName("dca_z_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(1)->GetTitle());
  sprintf(title,"%s vs %s",fDCAHisto->GetAxis(1)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  h2D->SetBit(TH1::kLogX);
  aFolderObj->Add(h2D);

  h1D = MakeStat1D(h2D,0,0);
  h1D->SetName("mean_dca_z_vs_pt");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h1D->GetYaxis()->SetTitle("mean_dca_z (cm)");
  sprintf(title,"mean_dca_z (cm) vs %s",fDCAHisto->GetAxis(3)->GetTitle());
  h1D->SetTitle(title);
  h1D->SetBit(TH1::kLogX);
  aFolderObj->Add(h1D);

  h1D = MakeStat1D(h2D,0,1);
  h1D->SetName("rms_dca_z_vs_pt");
  h1D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h1D->GetYaxis()->SetTitle("rms_dca_z (cm)");
  sprintf(title,"rms_dca_z (cm) vs %s",fDCAHisto->GetAxis(3)->GetTitle());
  h1D->SetTitle(title);
  h1D->SetBit(TH1::kLogX);
  aFolderObj->Add(h1D);
   
  /*
  h3D = fDCAHisto->Projection(2,3,0); // normal 3D projection convention
  h3D->SetName("dca_r_vs_eta_vs_pt");

  h2D = MakeStat2D(h3D,0,0,0);
  h2D->SetName("mean_dca_r_vs_eta_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetZaxis()->SetTitle("mean_dca_r (cm)");
  sprintf(title,"mean_dca_r (cm) vs %s vs %s",fDCAHisto->GetAxis(2)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  h2D = MakeStat2D(h3D,0,0,1);
  h2D->SetName("rms_dca_r_vs_eta_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetZaxis()->SetTitle("rms_dca_r (cm)");
  sprintf(title,"rms_dca_r (cm) vs %s vs %s",fDCAHisto->GetAxis(2)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  //
  h3D = fDCAHisto->Projection(2,3,1);
  h3D->SetName("dca_z_vs_eta_vs_pt");

  h2D = MakeStat2D(h3D,0,0,0);
  h2D->SetName("mean_dca_z_vs_eta_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetZaxis()->SetTitle("mean_dca_z (cm)");
  sprintf(title,"mean_dca_z (cm) vs %s vs %s",fDCAHisto->GetAxis(2)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);

  h2D = MakeStat2D(h3D,0,0,1);
  h2D->SetName("rms_dca_z_vs_eta_vs_pt");
  h2D->GetXaxis()->SetTitle(fDCAHisto->GetAxis(2)->GetTitle());
  h2D->GetYaxis()->SetTitle(fDCAHisto->GetAxis(3)->GetTitle());
  h2D->GetZaxis()->SetTitle("rms_dca_z (cm)");
  sprintf(title,"rms_dca_z (cm) vs %s vs %s",fDCAHisto->GetAxis(2)->GetTitle(),fDCAHisto->GetAxis(3)->GetTitle());
  h2D->SetTitle(title);
  aFolderObj->Add(h2D);
  */







  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TH1F* AliPerformanceDCA::MakeStat1D(TH2 *hist, Int_t delta0, Int_t type) 
{
  // Return TH1F histogram 
  // delta - number of bins to integrate
  // with mean (type == 0) or RMS (type==1) 

  char hname[256];
  const char* suffix = "_stat1d";
  sprintf(hname,"%s%s",hist->GetName(),suffix);
  TAxis* xaxis = hist->GetXaxis();
  Int_t  nbinx = xaxis->GetNbins();

  TH1F *hnew = (TH1F*)hist->ProjectionX()->Clone();
  hnew->SetName(hname);

  char name[256];
  for (Int_t ix=0; ix<=nbinx;ix++) {
    sprintf(name,"%s_%d",hist->GetName(),ix);
    TH1 *projection = hist->ProjectionY(name,ix-delta0,ix+delta0);

    Float_t stat= 0., stat_err =0.;
    if (type==0) { stat = projection->GetMean(); stat_err = projection->GetMeanError(); } 
    if (type==1) { stat = projection->GetRMS(); stat_err = projection->GetRMSError(); }
 
    hnew->SetBinContent(ix, stat);
    hnew->SetBinError(ix, stat_err);
  }
  
return hnew;
}

//_____________________________________________________________________________
TH2F* AliPerformanceDCA::MakeStat2D(TH3 *hist, Int_t delta0, Int_t delta1, Int_t type) 
{
  // Return TH1F histogram 
  // delta0 - number of bins to integrate in x
  // delta1 - number of bins to integrate in y
  // with mean (type==0) or RMS (type==1) 

  char hname[256];
  const char* suffix = "_stat2d";
  sprintf(hname,"%s%s",hist->GetName(),suffix);

  TAxis* xaxis = hist->GetXaxis();
  Int_t  nbinx = xaxis->GetNbins(); 

  TH2F *hnew = (TH2F*)hist->Project3D("yx")->Clone();
  hnew->SetName(hname);

  TAxis* yaxis = hist->GetYaxis();
  Int_t  nbiny = yaxis->GetNbins(); 

  char name[256];
  for (Int_t ix=0; ix<=nbinx;ix++) {
    for (Int_t iy=0; iy<=nbiny;iy++) {
      sprintf(name,"%s_%d_%d",hist->GetName(),ix,iy);
      TH1 *projection = hist->ProjectionZ(name,ix-delta0,ix+delta0,iy-delta1,iy+delta1);

      Float_t stat= 0., stat_err =0.;
      if (type==0) { stat = projection->GetMean(); stat_err = projection->GetMeanError(); } 
      if (type==1) { stat = projection->GetRMS(); stat_err = projection->GetRMSError(); }
     
      hnew->SetBinContent(ix,iy,stat);
      hnew->SetBinError(ix,iy,stat_err);
    }
  }
  
return hnew;
}

//_____________________________________________________________________________
TFolder* AliPerformanceDCA::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceDCA * comp=this;
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
TFolder* AliPerformanceDCA::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
