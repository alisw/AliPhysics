//------------------------------------------------------------------------------
// Implementation of AliPerformanceTPC class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformanceTPC.
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceTPC * compObj = (AliPerformanceTPC*)coutput->FindObject("AliPerformanceTPC");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderTPC" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_TPC.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TPostScript.h"

#include "AliPerformanceTPC.h" 
#include "AliESDEvent.h" 
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliStack.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTreeDraw.h" 

using namespace std;

ClassImp(AliPerformanceTPC)

//_____________________________________________________________________________
AliPerformanceTPC::AliPerformanceTPC():
  AliPerformanceObject("AliPerformanceTPC"),
  fTPCHisto(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  Init();
}

//_____________________________________________________________________________
AliPerformanceTPC::AliPerformanceTPC(Char_t* name="AliPerformanceTPC", Char_t* title="AliPerformanceTPC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),
  fTPCHisto(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  // named constructor	
  // 
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);

  Init();
}

//_____________________________________________________________________________
AliPerformanceTPC::~AliPerformanceTPC()
{
  // destructor
   
  if(fTPCHisto) delete fTPCHisto; fTPCHisto=0;     
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceTPC::Init(){
  //
  // histogram bining
  //

  // set pt bins
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

  // nTPCClust:chi2PerTPCClust:nTPCClustFindRatio:eta:phi:pt
  Int_t binsTPCHisto[6]={160,100,100,30,144,nPtBins};
  Double_t minTPCHisto[6]={0., 0., 0., -1.5, 0., ptMin};
  Double_t maxTPCHisto[6]={160.,10.,1.2, 1.5, 2.*TMath::Pi(), ptMax};

  fTPCHisto = new THnSparseF("fTPCHisto","nClust:chi2PerClust:nClust/nFindableClust:eta:phi:pt",6,binsTPCHisto,minTPCHisto,maxTPCHisto);
  fTPCHisto->SetBinEdges(5,binsPt);

  fTPCHisto->GetAxis(0)->SetTitle("nClust");
  fTPCHisto->GetAxis(1)->SetTitle("chi2PerClust");
  fTPCHisto->GetAxis(2)->SetTitle("nClust/nFindableClust");
  fTPCHisto->GetAxis(3)->SetTitle("#eta");
  fTPCHisto->GetAxis(4)->SetTitle("#phi (rad)");
  fTPCHisto->GetAxis(5)->SetTitle("p_{T} (GeV/c)");
  fTPCHisto->Sumw2();

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderTPC","Analysis Resolution Folder");
}

//_____________________________________________________________________________
void AliPerformanceTPC::ProcessTPC(AliStack* const stack, AliESDtrack *const esdTrack)
{
  if(!esdTrack) return;

  // Fill TPC only resolution comparison information 
  const AliExternalTrackParam *track = esdTrack->GetTPCInnerParam();
  if(!track) return;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

  //Float_t q = esdTrack->Charge();
  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();
  Int_t nClust = esdTrack->GetTPCclusters(0);
  Int_t nFindableClust = esdTrack->GetTPCNclsF();

  Float_t chi2PerCluster = 0.;
  if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2()/Float_t(nClust);

  Float_t clustPerFindClust = 0.;
  if(nFindableClust>0.) clustPerFindClust = Float_t(nClust)/nFindableClust;
  
  if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ()) {
    Double_t vTPCHisto[6] = {nClust,chi2PerCluster,clustPerFindClust,eta,phi,pt};
    fTPCHisto->Fill(vTPCHisto); 
  }
 
  //
  // Fill rec vs MC information
  //
  if(!stack) return;

}

//_____________________________________________________________________________
void AliPerformanceTPC::ProcessTPCITS(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill comparison information (TPC+ITS) 
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}
 
//_____________________________________________________________________________
void AliPerformanceTPC::ProcessConstrained(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill comparison information (constarained parameters) 
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}
 
//_____________________________________________________________________________
void AliPerformanceTPC::Exec(AliMCEvent* const mcEvent, AliESDEvent* const esdEvent, const Bool_t bUseMC)
{
  // Process comparison information 
  //
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
    AliESDtrack*track = esdEvent->GetTrack(iTrack);
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
void AliPerformanceTPC::Analyse() {
  //
  // Analyse comparison information and store output histograms
  // in the folder "folderTPC"
  //
  TH1::AddDirectory(kFALSE);
  //TH1F *h=0;
  TH2F *h2D=0;
  TObjArray *aFolderObj = new TObjArray;

  char name[256];
  char title[256];
  for(Int_t i=0; i<5; i++) 
  {
    for(Int_t j=i+1; j<6; j++) 
    {
      if(j==5) fTPCHisto->GetAxis(5)->SetRangeUser(0.1,10.);
      h2D = (TH2F*)fTPCHisto->Projection(i,j);
      sprintf(name,"h_tpc_%d_vs_%d",i,j);
      h2D->SetName(name);
      h2D->GetXaxis()->SetTitle(fTPCHisto->GetAxis(j)->GetTitle());
      h2D->GetYaxis()->SetTitle(fTPCHisto->GetAxis(i)->GetTitle());
      sprintf(title,"%s vs %s",fTPCHisto->GetAxis(j)->GetTitle(),fTPCHisto->GetAxis(i)->GetTitle());
      h2D->SetTitle(title);

      if(j==5) h2D->SetBit(TH1::kLogX);
      aFolderObj->Add(h2D);
    }  
  }

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliPerformanceTPC::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceTPC * comp=this;
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
Long64_t AliPerformanceTPC::Merge(TCollection* const list) 
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
  AliPerformanceTPC* entry = dynamic_cast<AliPerformanceTPC*>(obj);
  if (entry == 0) continue; 

  fTPCHisto->Add(entry->fTPCHisto);

  count++;
  }

return count;
}

//_____________________________________________________________________________
TFolder* AliPerformanceTPC::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
