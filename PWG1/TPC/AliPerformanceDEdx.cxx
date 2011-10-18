//------------------------------------------------------------------------------
// Implementation of AliPerformanceDEdx class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data of AliPerformanceDEdx.
//  
// Author: J.Otwinowski 04/02/2008 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  //AliPerformanceDEdx * compObj = (AliPerformanceDEdx*)f.Get("AliPerformanceDEdx");
  AliPerformanceDEdx * compObj = (AliPerformanceDEdx*)coutput->FindObject("AliPerformanceDEdx");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderDEdx" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_DEdx.root","recreate");
  //compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TDirectory.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSystem.h"
#include "TChain.h"

#include "AliPerformanceDEdx.h"
#include "AliPerformanceTPC.h"
#include "AliTPCPerformanceSummary.h"
#include "AliESDEvent.h"
#include "AliTracker.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliLog.h" 
#include "AliMCInfoCuts.h" 
#include "AliMathBase.h"
#include "AliRecInfoCuts.h" 
#include "AliTreeDraw.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

using namespace std;

ClassImp(AliPerformanceDEdx)

Bool_t AliPerformanceDEdx::fgMergeTHnSparse = kFALSE;
Bool_t AliPerformanceDEdx::fgUseMergeTHnSparse = kFALSE;

//_____________________________________________________________________________
AliPerformanceDEdx::AliPerformanceDEdx():
  AliPerformanceObject("AliPerformanceDEdx"),

  // dEdx 
  fDeDxHisto(0),
  fFolderObj(0),
  
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
AliPerformanceDEdx::AliPerformanceDEdx(Char_t* name="AliPerformanceDEdx", Char_t* title="AliPerformanceDEdx",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),

  // dEdx 
  fDeDxHisto(0),
  fFolderObj(0),
  
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
AliPerformanceDEdx::~AliPerformanceDEdx()
{
  // destructor
  if(fDeDxHisto)  delete fDeDxHisto; fDeDxHisto=0; 
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::Init()
{
  // Init histograms
  
  // TPC dEdx
  // set p bins
  Int_t nPBins = 50;
  Double_t pMin = 1.e-2, pMax = 20.;

  Double_t *binsP = 0;

  if (IsHptGenerator())  { 
        pMax = 100.;
  } 
   binsP = CreateLogAxis(nPBins,pMin,pMax);


  /*
  Int_t nPBins = 31;
    Double_t binsP[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
    Double_t pMin = 0., pMax = 10.;

    if(IsHptGenerator() == kTRUE) {
      nPBins = 100;
      pMin = 0.; pMax = 100.;
    }
   */

   //dedx:alpha:y:z:snp:tgl:ncls:p
   //dedx:phi:y:z:snp:tgl:ncls:p
   //Int_t binsQA[8]    = {300, 50, 50,  50, 50, 50, 80, nPBins};
   //Double_t xminQA[8] = {0, -4,-20,-250, -1, -2, 0, pMin};
   //Double_t xmaxQA[8] = {300, 4, 20, 250,  1,  2, 160, pMax};
  Int_t binsQA[10]    = {300, 144, 50,  50, 50, 50, 80, nPBins, 160,100};
  Double_t xminQA[10] = {0, -TMath::Pi(),-20,-250, -1, -2, 0, pMin, 0., 0.};
  Double_t xmaxQA[10] = {300, TMath::Pi(), 20, 250,  1,  2, 160, pMax ,160., 1.};

   //fDeDxHisto = new THnSparseF("fDeDxHisto","dedx:alpha:y:z:snp:tgl:ncls:momentum",8,binsQA,xminQA,xmaxQA);
   fDeDxHisto = new THnSparseF("fDeDxHisto","dedx:phi:y:z:snp:tgl:ncls:momentum:TPCSignalN:clsF",10,binsQA,xminQA,xmaxQA);
   fDeDxHisto->SetBinEdges(7,binsP);

   fDeDxHisto->GetAxis(0)->SetTitle("dedx (a.u.)");
   fDeDxHisto->GetAxis(1)->SetTitle("#phi (rad)");
   fDeDxHisto->GetAxis(2)->SetTitle("y (cm)");
   fDeDxHisto->GetAxis(3)->SetTitle("z (cm)");
   fDeDxHisto->GetAxis(4)->SetTitle("sin#phi");
   fDeDxHisto->GetAxis(5)->SetTitle("tan#lambda");
   fDeDxHisto->GetAxis(6)->SetTitle("ncls");
   fDeDxHisto->GetAxis(7)->SetTitle("p (GeV/c)");
   fDeDxHisto->GetAxis(8)->SetTitle("number of cls used for dEdx");
   fDeDxHisto->GetAxis(9)->SetTitle("number of cls found over findable");
   //fDeDxHisto->Sumw2();

   // Init cuts
   if(!fCutsMC) {
     AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
   }
   if(!fCutsRC) {
     AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
   }

   // init folder
   fAnalysisFolder = CreateFolder("folderDEdx","Analysis de/dx Folder");
   
   // save merge status in object
   fMergeTHnSparseObj = fgMergeTHnSparse;

}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessTPC(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill dE/dx  comparison information
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessInnerTPC(AliStack* const stack, AliESDtrack *const esdTrack, AliESDEvent* const esdEvent)
{
 //
 // Fill TPC track information at inner TPC wall
 // 
  if(!esdEvent) return;
  if(!esdTrack) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliESDVertex *vtxESD = esdEvent->GetPrimaryVertexTracks();
    Double_t x[3]; esdTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = kFALSE;
    if(fabs(b[2])>0.000001)
      isOK = esdTrack->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
    //    Bool_t isOK = esdTrack->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }

  // get external param. at inner TPC wall
  const AliExternalTrackParam *innerParam =  esdTrack->GetInnerParam();
  if(!innerParam) return;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

  if((esdTrack->GetStatus()&AliESDtrack::kTPCrefit)==0) return; // TPC refit

  //
  // select primaries
  //
  Double_t dcaToVertex = -1;
  if( fCutsRC->GetDCAToVertex2D() ) 
  {
      dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC->GetMaxDCAToVertexXY()/fCutsRC->GetMaxDCAToVertexXY()                    + dca[1]*dca[1]/fCutsRC->GetMaxDCAToVertexZ()/fCutsRC->GetMaxDCAToVertexZ()); 
  }
  if(fCutsRC->GetDCAToVertex2D() && dcaToVertex > 1) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC->GetMaxDCAToVertexXY()) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC->GetMaxDCAToVertexZ()) return;

  Float_t dedx = esdTrack->GetTPCsignal();
  Int_t ncls = esdTrack->GetTPCNcls();
  Int_t TPCSignalN = esdTrack->GetTPCsignalN();
  Float_t nClsF = esdTrack->GetTPCClusterInfo(2,0);


  Double_t pt = innerParam->Pt();
  Double_t lam = TMath::ATan2(innerParam->Pz(),innerParam->Pt());
  Double_t p = pt/TMath::Cos(lam);
  //Double_t alpha = innerParam->GetAlpha();
  Double_t phi = TMath::ATan2(innerParam->Py(),innerParam->Px());
  //if(phi<0.) phi += 2.*TMath::Phi();
  Double_t y = innerParam->GetY();
  Double_t z = innerParam->GetZ();
  Double_t snp = innerParam->GetSnp();
  Double_t tgl = innerParam->GetTgl();

  //Double_t vDeDxHisto[8] = {dedx,alpha,y,z,snp,tgl,ncls,p};
  Double_t vDeDxHisto[10] = {dedx,phi,y,z,snp,tgl,ncls,p,TPCSignalN,nClsF};
  fDeDxHisto->Fill(vDeDxHisto); 

  if(!stack) return;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessTPCITS(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill dE/dx  comparison information
  
   AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessConstrained(AliStack* const /*stack*/, AliESDtrack *const /*esdTrack*/)
{
  // Fill dE/dx  comparison information
  
   AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
Long64_t AliPerformanceDEdx::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;
  
  Bool_t merge = ((fgUseMergeTHnSparse && fgMergeTHnSparse) || (!fgUseMergeTHnSparse && fMergeTHnSparseObj));

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  TObjArray* objArrayList = 0;
  objArrayList = new TObjArray();

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliPerformanceDEdx* entry = dynamic_cast<AliPerformanceDEdx*>(obj);
    if (entry == 0) continue; 
    if (merge) {
        if ((fDeDxHisto) && (entry->fDeDxHisto)) { fDeDxHisto->Add(entry->fDeDxHisto); }        
    }
    // the analysisfolder is only merged if present
    if (entry->fFolderObj) { objArrayList->Add(entry->fFolderObj); }

    count++;
  }
  if (fFolderObj) { fFolderObj->Merge(objArrayList); } 
  // to signal that track histos were not merged: reset
  if (!merge) { fDeDxHisto->Reset(); }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;  

return count;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend)
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

  // use ESD friends
  if(bUseESDfriend) {
    if(!esdFriend) {
      AliDebug(AliLog::kError, "esdFriend not available");
      return;
    }
  }

  // trigger
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = esdEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) return; 
  }

  // get event vertex
  const AliESDVertex *vtxESD = NULL;
  if( IsUseTrackVertex() ) 
  { 
    // track vertex
    vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    // TPC track vertex
    vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  if(vtxESD && (vtxESD->GetStatus()<=0)) return;

  //  Process events
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if(!track) continue;

    if(GetAnalysisMode() == 0) ProcessTPC(stack,track);
    else if(GetAnalysisMode() == 1) ProcessTPCITS(stack,track);
    else if(GetAnalysisMode() == 2) ProcessConstrained(stack,track);
    else if(GetAnalysisMode() == 3) ProcessInnerTPC(stack,track,esdEvent);
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }
}

//_____________________________________________________________________________
void AliPerformanceDEdx::Analyse()
{
  // Analyze comparison information and store output histograms
  // in the folder "folderDEdx"
  //
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kFALSE);
  TH1F *h1D=0;
  TH2F *h2D=0;
  TObjArray *aFolderObj = new TObjArray;
  TString selString;

  char name[256];
  char title[256];

  for(Int_t i=1; i<10; i++) { 
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, i);
  }

  // resolution histograms for mips
  //dedx:phi:y:z:snp:tgl:ncls:p
  fDeDxHisto->GetAxis(2)->SetRangeUser(-15.,14.999);
  fDeDxHisto->GetAxis(3)->SetRangeUser(-120.,119.999);
  fDeDxHisto->GetAxis(4)->SetRangeUser(-0.4, 0.399);
  fDeDxHisto->GetAxis(5)->SetRangeUser(-0.9,0.89);
  fDeDxHisto->GetAxis(6)->SetRangeUser(60.,160.);
  fDeDxHisto->GetAxis(7)->SetRangeUser(0.4,0.499);
  fDeDxHisto->GetAxis(8)->SetRangeUser(60.,160.);
  
  
  selString = "mipsres";
  AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, &selString);

  //
  TObjArray *arr[9] = {0};
  TF1 *f1[9] = {0};
  
  for(Int_t i=0; i<9; i++) 
  { 
    arr[i] = new TObjArray;
    f1[i] = new TF1("gaus","gaus");
    //printf("i %d \n",i);

    h2D = (TH2F*)fDeDxHisto->Projection(0,i+1);

    f1[i]->SetRange(40,60); // should be pion peak
    h2D->FitSlicesY(f1[i],0,-1,10,"QNR",arr[i]); // gaus fit of pion peak

    h1D = (TH1F*)arr[i]->At(1);
    snprintf(name,256,"mean_dedx_mips_vs_%d",i);
    h1D->SetName(name);
    snprintf(title,256,"%s vs %s","mean_dedx_mips (a.u.)",fDeDxHisto->GetAxis(i)->GetTitle());
    h1D->SetTitle(title);
    h1D->GetXaxis()->SetTitle(fDeDxHisto->GetAxis(i)->GetTitle());
    h1D->GetYaxis()->SetTitle("mean_dedx_mips (a.u.)");
    //h1D->SetMinimum(40);
    //h1D->SetMaximum(60);

    aFolderObj->Add(h1D);

    h1D = (TH1F*)arr[i]->At(2);
    snprintf(name,256,"res_dedx_mips_vs_%d",i);
    h1D->SetName(name);
    snprintf(title,256,"%s vs %s","res_dedx_mips (a.u)",fDeDxHisto->GetAxis(i)->GetTitle());
    h1D->SetTitle(title);
    h1D->GetXaxis()->SetTitle(fDeDxHisto->GetAxis(i)->GetTitle());
    h1D->GetYaxis()->SetTitle("res_dedx_mips (a.u.)");
    //h1D->SetMinimum(0);
    //h1D->SetMaximum(6);

    aFolderObj->Add(h1D);
  }

    // select MIPs (version from AliTPCPerfomanceSummary)
    fDeDxHisto->GetAxis(0)->SetRangeUser(35,60);
    fDeDxHisto->GetAxis(2)->SetRangeUser(-20,20);
    fDeDxHisto->GetAxis(3)->SetRangeUser(-250,250);
    fDeDxHisto->GetAxis(4)->SetRangeUser(-1, 1);
    fDeDxHisto->GetAxis(5)->SetRangeUser(-1,1);
    fDeDxHisto->GetAxis(6)->SetRangeUser(80,160);
    fDeDxHisto->GetAxis(7)->SetRangeUser(0.4,0.55);
    fDeDxHisto->GetAxis(8)->SetRangeUser(80,160);
    fDeDxHisto->GetAxis(9)->SetRangeUser(0.5,1.);

    selString = "mips";
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, &selString);
    
    selString = "mips_C";
    fDeDxHisto->GetAxis(5)->SetRangeUser(-3,0);
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 5, &selString);
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 1, &selString);
    
    selString = "mips_A";
    fDeDxHisto->GetAxis(5)->SetRangeUser(0,3);
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 5, &selString);    
    AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 1, &selString);
    
    //restore cuts
    fDeDxHisto->GetAxis(0)->SetRangeUser(0,300);
    fDeDxHisto->GetAxis(2)->SetRangeUser(-20,20);
    fDeDxHisto->GetAxis(3)->SetRangeUser(-250,250);
    fDeDxHisto->GetAxis(4)->SetRangeUser(-1, 1);
    fDeDxHisto->GetAxis(5)->SetRangeUser(-3,3);
    fDeDxHisto->GetAxis(6)->SetRangeUser(0,160);
    fDeDxHisto->GetAxis(7)->SetRangeUser(1e-2,10);    
    fDeDxHisto->GetAxis(8)->SetRangeUser(0,160);
    fDeDxHisto->GetAxis(9)->SetRangeUser(0.,1.);




    printf("exportToFolder\n");
    // export objects to analysis folder
    fAnalysisFolder = ExportToFolder(aFolderObj);
    if (fFolderObj) delete fFolderObj;
    fFolderObj = aFolderObj;
    aFolderObj=0;


  for(Int_t i=0;i<9;i++) { 
    if(f1[i]) delete f1[i]; f1[i]=0;
  }

}

//_____________________________________________________________________________
TFolder* AliPerformanceDEdx::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceDEdx * comp=this;
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
TFolder* AliPerformanceDEdx::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}

//_____________________________________________________________________________
TTree* AliPerformanceDEdx::CreateSummary()
{
    // implementaion removed, switched back to use AliPerformanceSummary (now called in AliPerformanceTask)
    return 0;
}

