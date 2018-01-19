//------------------------------------------------------------------------------
// Implementation of AliPerformanceDEdx class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data of AliPerformanceDEdx.
//  
// Author: J.Otwinowski   04/02/2008 
// Changes by M.Knichel   15/10/2010
// Changes by J.Salzwedel 15/10/2014
//------------------------------------------------------------------------------

/*AliPerformanceDEdx.cxx
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/Macros/LoadMyLibs.C");
  LoadMyLibs();cd /hera/alice/atarant/train/trunk/atarant_spectra/qa/

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
#include "AliVEvent.h"
#include "AliTracker.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
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
AliPerformanceDEdx::AliPerformanceDEdx(TRootIOCtor* b):
 AliPerformanceObject(b),
  // dEdx 
  fDeDxHisto(0),
  // histogram folder 
  fAnalysisFolder(0),
  fFolderObj(0),
  h_tpc_dedx_mips_0(NULL),
  h_tpc_dedx_mipsele_0(0),
  h_tpc_dedx_mips_c_0_5(0),
  h_tpc_dedx_mips_a_0_5(0),
  h_tpc_dedx_mips_c_0_1(0),
  h_tpc_dedx_mips_a_0_1(0)

{
  // io constructor
}

//_____________________________________________________________________________
AliPerformanceDEdx::AliPerformanceDEdx(const Char_t* name, const Char_t* title, Int_t analysisMode, Bool_t hptGenerator, Bool_t useSparse):
 AliPerformanceObject(name,title),

  // dEdx 
  fDeDxHisto(0),
  // histogram folder 
  fAnalysisFolder(0),
  fFolderObj(0),
  h_tpc_dedx_mipsele_0(0),
  h_tpc_dedx_mips_c_0_5(0),
  h_tpc_dedx_mips_a_0_5(0),
  h_tpc_dedx_mips_c_0_1(0),
  h_tpc_dedx_mips_a_0_1(0)
{
  // named constructor

  fUseSparse = useSparse;
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);
  Init();
}


//_____________________________________________________________________________
AliPerformanceDEdx::~AliPerformanceDEdx()
{
  // destructor
  delete fDeDxHisto;

  if (fFolderObj && fAnalysisFolder && !fAnalysisFolder->IsOwner()) {
    fFolderObj->Delete();
  } //delete the registered non-sparse histograms

  delete fFolderObj;
  delete fAnalysisFolder;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::Init()
{

    if(!fUseSparse) fFolderObj = new TObjArray;

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

   //Int_t binsQA[8]    = {300, 50, 50,  50, 50, 50, 80, nPBins};
   //Double_t xminQA[8] = {0, -4,-20,-250, -1, -2, 0, pMin};
   //Double_t xmaxQA[8] = {300, 4, 20, 250,  1,  2, 160, pMax};
 // signal:phi:y:z:snp:tgl:ncls:p:nclsDEdx:nclsF
  Int_t binsQA[10]    = {300, 144, 50,  50, 50, 50, 80, nPBins, 160, 80};
  Double_t xminQA[10] = {0, -TMath::Pi(),-20,-250, -1, -2, 0, pMin, 0., 0.};
  Double_t xmaxQA[10] = {300, TMath::Pi(), 20, 250,  1,  2, 160, pMax ,160., 1.};

    if(fUseSparse){
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
    }
    else{
    
        h_tpc_dedx_mips_0 = new TH1D("h_tpc_dedx_mips_0","",binsQA[0],xminQA[0],xmaxQA[0]);
        h_tpc_dedx_mips_0->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mipsele_0 = new TH1D("h_tpc_dedx_mipsele_0","",binsQA[0],xminQA[0],xmaxQA[0]);
        h_tpc_dedx_mipsele_0->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mipsele_0->SetYTitle("tan#lambda");
        h_tpc_dedx_mips_c_0_5 = new TH2D("h_tpc_dedx_mips_c_0_5","",binsQA[0],xminQA[0],xmaxQA[0],binsQA[5],xminQA[5],xmaxQA[5]);
        h_tpc_dedx_mips_c_0_5->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mips_c_0_5->SetYTitle("tan#lambda");
        h_tpc_dedx_mips_a_0_5 = new TH2D("h_tpc_dedx_mips_a_0_5","",binsQA[0],xminQA[0],xmaxQA[0],binsQA[5],xminQA[5],xmaxQA[5]);
        h_tpc_dedx_mips_a_0_5->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mips_a_0_5->SetYTitle("tan#lambda");
        h_tpc_dedx_mips_c_0_1 = new TH2D("h_tpc_dedx_mips_c_0_1","",binsQA[0],xminQA[0],xmaxQA[0],binsQA[1],xminQA[1],xmaxQA[1]);
        h_tpc_dedx_mips_c_0_1->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mips_c_0_1->SetYTitle("#phi (rad)");
        h_tpc_dedx_mips_a_0_1 = new TH2D("h_tpc_dedx_mips_a_0_1","",binsQA[0],xminQA[0],xmaxQA[0],binsQA[1],xminQA[1],xmaxQA[1]);
        h_tpc_dedx_mips_a_0_1->SetXTitle("dedx (a.u.)");
        h_tpc_dedx_mips_a_0_1->SetYTitle("#phi (rad)");
        
        fFolderObj->Add(h_tpc_dedx_mips_0);
        fFolderObj->Add(h_tpc_dedx_mipsele_0);
        fFolderObj->Add(h_tpc_dedx_mips_c_0_5);
        fFolderObj->Add(h_tpc_dedx_mips_a_0_5);
        fFolderObj->Add(h_tpc_dedx_mips_c_0_1);
        fFolderObj->Add(h_tpc_dedx_mips_a_0_1);
        
    }
    
   //fDeDxHisto->Sumw2();

   // init folder
   fAnalysisFolder = CreateFolder("folderDEdx","Analysis de/dx Folder");

   // save merge status in object
   fMergeTHnSparseObj = fgMergeTHnSparse;

}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessTPC(AliMCEvent* const /*mcev*/, AliVTrack *const /*vTrack*/)
{
  // Fill dE/dx  comparison information
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessInnerTPC(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent)
{
    //
    // Fill TPC track information at inner TPC wall
    // Only ESD events store TPC track information at inner TPC wall
    //
    if(!vEvent || !vTrack) return;
    
    
    AliExternalTrackParam trackParams;
    vTrack->GetTrackParam(trackParams);
    AliExternalTrackParam *etpTrack = &trackParams;

    AliExternalTrackParam innerTPCtrackParams;
    vTrack->GetTrackParamTPCInner(innerTPCtrackParams);
    AliExternalTrackParam* innerTPCparam = &innerTPCtrackParams;

    Double_t dca[2] = {0.,0.};
    Double_t cov[3] = {0.,0.,0.};
    if( IsUseTrackVertex() ) {
        // Relate TPC inner params to prim. vertex
        AliESDVertex vertex;
        vEvent->GetPrimaryVertexTracks(vertex);
        const AliVVertex *vVertex = &vertex;
        Double_t x[3]; etpTrack->GetXYZ(x);
        Double_t b[3]; AliTracker::GetBxByBz(x,b);
        Bool_t isOK=kFALSE;
        if(fabs(b[2])>0.000001)
            isOK = innerTPCparam->RelateToVVertexBxByBzDCA(vVertex, b, kVeryBig, NULL, dca, cov);
        if(!isOK) return;
    }
    
    AliExternalTrackParam innerTrackParams;
    vTrack->GetTrackParamIp(innerTrackParams);
    AliExternalTrackParam* innerParam = &innerTrackParams;

    if((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
    
    //
    // select primaries
    //
    Double_t dcaToVertex = -1;
    if( fCutsRC.GetDCAToVertex2D() )
    {
        dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC.GetMaxDCAToVertexXY()/fCutsRC.GetMaxDCAToVertexXY()+dca[1]*dca[1]/fCutsRC.GetMaxDCAToVertexZ()/fCutsRC.GetMaxDCAToVertexZ());
    }
    if(fCutsRC.GetDCAToVertex2D() && dcaToVertex > 1) return;
    if(!fCutsRC.GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC.GetMaxDCAToVertexXY()) return;
    if(!fCutsRC.GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC.GetMaxDCAToVertexZ()) return;
    
    Float_t dedx = vTrack->GetTPCsignal();
    Int_t ncls = vTrack->GetTPCNcls();
    Int_t TPCSignalN = vTrack->GetTPCsignalN();
    //Float_t nCrossedRows = vTrack->GetTPCClusterInfo(2,1);
    Float_t nClsF = vTrack->GetTPCClusterInfo(2,0);
    
    
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
    
    //fill thnspars here coud add oroc mdedium long..............Atti
    //you should select which pad leng here
    // http://svnweb.cern.ch/world/wsvn/AliRoot/trunk/STEER/STEERBase/AliTPCdEdxInfo.h
    // fTPCsignalRegion[4];
    
    //Double_t vDeDxHisto[10] = {dedx,phi,y,z,snp,tgl,ncls,p,TPCSignalN,nCrossedRows};
    Double_t vDeDxHisto[10] = {dedx,phi,y,z,snp,tgl,Double_t(ncls),p,Double_t(TPCSignalN),nClsF};
    if(fUseSparse) fDeDxHisto->Fill(vDeDxHisto);
    else  FilldEdxHisotgram(vDeDxHisto);
    
    if(!mcev) return;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessTPCITS(AliMCEvent* const /*mcev*/, AliVTrack *const /*vTrack*/)
{
  // Fill dE/dx  comparison information
  
   AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
void AliPerformanceDEdx::ProcessConstrained(AliMCEvent* const /*mcev*/, AliVTrack *const /*vTrack*/)
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
  if (!merge) {
      if(fDeDxHisto) fDeDxHisto->Reset();
  }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;  

return count;
}

//_____________________________________________________________________________
void AliPerformanceDEdx::Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend)
{
  // Process comparison information 
  //
  if(!vEvent)
  {
      AliDebug(AliLog::kError, "esdEvent not available");
      return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
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

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

  } // end bUseMC

  // use V friends
  if(bUseVfriend) {
    if(!vFriendEvent) {
      AliDebug(AliLog::kError, "vFriend not available");
      return;
    }
  }

  // trigger
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = vEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) return; 
  }

  // get event vertex
  const AliVVertex *vVertex = NULL;
  if( IsUseTrackVertex() ) 
  { 
    // track vertex
    vVertex = vEvent->GetPrimaryVertexTracks();
  } else {
    // TPC track vertex
    vVertex = vEvent->GetPrimaryVertexTPC();
  }
  if(vVertex && (vVertex->GetStatus()<=0)) return;
  if(!vVertex) {
    printf("ERROR: Could not determine primary vertex");
    return;
  }
  
  //  Process events
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  {
    AliVParticle *particle = vEvent->GetTrack(iTrack);
    if(!particle) continue;
    AliVTrack *track = dynamic_cast<AliVTrack*>(particle);
    if(!track) continue;

    if(GetAnalysisMode() == 0) ProcessTPC(mcEvent,track);
    else if(GetAnalysisMode() == 1) ProcessTPCITS(mcEvent,track);
    else if(GetAnalysisMode() == 2) ProcessConstrained(mcEvent,track);
    else if(GetAnalysisMode() == 3) ProcessInnerTPC(mcEvent,track,vEvent);
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
  //Atti h_tpc_dedx_mips_0 
  //fai fit con range p(.32,.38) and dEdx(65- 120 or 100) e ripeti cosa fatta per pion e fai trending della media e res, poio la loro differenza
  //fai dedx vs lamda ma for e e pion separati
  //
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kFALSE);
    if(fUseSparse){
      TH1F *h1D=0;
      TH2F *h2D=0;
      TObjArray *aFolderObj = new TObjArray;
      TString selString;

      char name[256];
      char title[256];

      for(Int_t i=1; i<10; i++) { 
        AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, i);
      }

        AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 6, 7);
        AddProjection(aFolderObj, "dedx", fDeDxHisto, 7, 8, 9);
        AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, 8, 9);
        AddProjection(aFolderObj, "dedx", fDeDxHisto, 6, 8, 9);

      // resolution histograms for mips
      //-> signal:phi:y:z:snp:tgl:ncls:p:nclsDEdx:nclsF
      fDeDxHisto->GetAxis(2)->SetRangeUser(-15.,14.999);
      fDeDxHisto->GetAxis(3)->SetRangeUser(-120.,119.999);
      fDeDxHisto->GetAxis(4)->SetRangeUser(-0.4, 0.399);
      fDeDxHisto->GetAxis(5)->SetRangeUser(-0.9,0.89);
      fDeDxHisto->GetAxis(6)->SetRangeUser(60.,160.);
      fDeDxHisto->GetAxis(7)->SetRangeUser(0.4,0.499); //p
      fDeDxHisto->GetAxis(8)->SetRangeUser(60.,160.);
      
     
      selString = "mipsres";
      AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, &selString);

      //
      TObjArray *arr[10] = {0};
      TF1 *f1[10] = {0};
      
      for(Int_t i=1; i<10; i++) 
      { 
        arr[i] = new TObjArray;
        f1[i] = new TF1("gaus","gaus");
        //printf("i %d \n",i);

        h2D = (TH2F*)fDeDxHisto->Projection(0,i);

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
        fDeDxHisto->GetAxis(2)->SetRangeUser(-20,19.999);
        fDeDxHisto->GetAxis(3)->SetRangeUser(-250,249.999);
        fDeDxHisto->GetAxis(4)->SetRangeUser(-1, 0.99);
        fDeDxHisto->GetAxis(5)->SetRangeUser(-1,0.99);
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
        
        //////////////////////////////////////// atti new start
        // 
        // select (version from AliTPCPerfomanceSummary) electrons                                                                                                      
        fDeDxHisto->GetAxis(0)->SetRangeUser(70,100); //dedx for electrons
        fDeDxHisto->GetAxis(2)->SetRangeUser(-20,19.999);
        fDeDxHisto->GetAxis(3)->SetRangeUser(-250,249.999);
        fDeDxHisto->GetAxis(4)->SetRangeUser(-1, 0.99);
        fDeDxHisto->GetAxis(5)->SetRangeUser(-1,0.99);
        fDeDxHisto->GetAxis(6)->SetRangeUser(80,160);
        fDeDxHisto->GetAxis(7)->SetRangeUser(0.32,0.38); //momenta for electrons
        fDeDxHisto->GetAxis(8)->SetRangeUser(80,160);
        fDeDxHisto->GetAxis(9)->SetRangeUser(0.5,1.);

        selString = "mipsele";
        AddProjection(aFolderObj, "dedx", fDeDxHisto, 0, &selString);
        //////////////////////////////////////// atti new stop

        //restore cuts
        for (Int_t i=0; i<fDeDxHisto->GetNdimensions(); i++) {
          fDeDxHisto->GetAxis(i)->SetRange(1,fDeDxHisto->GetAxis(i)->GetNbins());
        }

        printf("exportToFolder\n");
        // export objects to analysis folder
        fAnalysisFolder = ExportToFolder(aFolderObj);
        if (fFolderObj) delete fFolderObj;
        fFolderObj = aFolderObj;
        aFolderObj=0;


      for(Int_t i=0;i<10;i++) { 
        if(f1[i]) delete f1[i]; f1[i]=0;
      }
    }
    else {
        printf("exportToFolder\n");
        fAnalysisFolder = ExportToFolder(fFolderObj);
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

void AliPerformanceDEdx::FilldEdxHisotgram(double *vDeDxHisto){
    
    if(!vDeDxHisto) return;
    if(vDeDxHisto[0] >= 35 && vDeDxHisto[0] < 60)
    if(vDeDxHisto[2] >= -20 && vDeDxHisto[2] < 20)
    if(vDeDxHisto[3] >= -250 && vDeDxHisto[3] < 250)
    if(vDeDxHisto[4] >= -1 && vDeDxHisto[4] < 1)
    if(vDeDxHisto[6] >= 80 && vDeDxHisto[6] < 160)
    if(vDeDxHisto[7] >= 0.4 && vDeDxHisto[7] < 0.55)
    if(vDeDxHisto[8] >= 80 && vDeDxHisto[8] < 160)
    if(vDeDxHisto[9] >= 0.5 && vDeDxHisto[9] < 1.){
        if(vDeDxHisto[5] > -1 && vDeDxHisto[5] < 1) h_tpc_dedx_mips_0->Fill(vDeDxHisto[0]);
        if(vDeDxHisto[5] > -3 && vDeDxHisto[5] < 0){
           if(h_tpc_dedx_mips_c_0_5) h_tpc_dedx_mips_c_0_5->Fill(vDeDxHisto[0],vDeDxHisto[5]);
           if(h_tpc_dedx_mips_c_0_1) h_tpc_dedx_mips_c_0_1->Fill(vDeDxHisto[0],vDeDxHisto[1]);
        }
        else if(vDeDxHisto[5] > 0 && vDeDxHisto[5] < 3){
            if(h_tpc_dedx_mips_a_0_5) h_tpc_dedx_mips_a_0_5->Fill(vDeDxHisto[0],vDeDxHisto[5]);
            if(h_tpc_dedx_mips_a_0_1) h_tpc_dedx_mips_a_0_1->Fill(vDeDxHisto[0],vDeDxHisto[1]);

        }
    }
    if(vDeDxHisto[0] >= 70 && vDeDxHisto[0] < 100) /*dedx for electrons*/
    if(vDeDxHisto[2] >= -20 && vDeDxHisto[2] < 19.999)
    if(vDeDxHisto[3] >= -250 && vDeDxHisto[3] < 249.999)
    if(vDeDxHisto[4] >= -1 && vDeDxHisto[4] < 0.99)
    if(vDeDxHisto[5] >= -1 && vDeDxHisto[5] < 0.99)
    if(vDeDxHisto[6] >= 80 && vDeDxHisto[6] < 160)
    if(vDeDxHisto[7] >= 0.32 && vDeDxHisto[7] < 0.38)  /*momenta for electrons*/
    if(vDeDxHisto[8] >= 80 && vDeDxHisto[8] < 160)
    if(vDeDxHisto[9] >= 0.5 && vDeDxHisto[9] < 1.)
    if(h_tpc_dedx_mipsele_0) h_tpc_dedx_mipsele_0->Fill(vDeDxHisto[0]);

}

void AliPerformanceDEdx::ResetOutputData(){

    if(fUseSparse){
        if(fDeDxHisto) fDeDxHisto->Reset("ICE");
    }
    else{
        if(h_tpc_dedx_mips_0) h_tpc_dedx_mips_0->Reset("ICE");
        if(h_tpc_dedx_mipsele_0) h_tpc_dedx_mipsele_0->Reset("ICE");
        if(h_tpc_dedx_mips_c_0_5) h_tpc_dedx_mips_c_0_5->Reset("ICE");
        if(h_tpc_dedx_mips_a_0_5) h_tpc_dedx_mips_a_0_5->Reset("ICE");
        if(h_tpc_dedx_mips_c_0_1) h_tpc_dedx_mips_c_0_1->Reset("ICE");
        if(h_tpc_dedx_mips_a_0_1) h_tpc_dedx_mips_a_0_1->Reset("ICE");
    }
    
}

//_____________________________________________________________________________
TCollection* AliPerformanceDEdx::GetListOfDrawableObjects()
{
  TObjArray* tmp = fFolderObj;
  fFolderObj = NULL;
  if (fAnalysisFolder) { fAnalysisFolder->SetOwner(kFALSE); }
  return tmp;
}

