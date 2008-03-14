//------------------------------------------------------------------------------
// Implementation of AliComparisonDCA class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms) are stored in the output picture_dca.root file.
//  
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
  // after running analysis, read the file, and get component
  gSystem->Load("libPWG1.so");
  TFile f("Output.root");
  AliComparisonDCA * comp = (AliComparisonDCA*)f.Get("AliComparisonDCA");

  // analyse comparison data (output stored in pictures_dca.root)
  comp->Analyse();

  // TPC track length parameterisation
  TF1 fl("fl","((min(250./(abs(x+0.000001)),250)-90))",0,2);  // length function
  TF1 fl2("fl2","[0]/((min(250./(abs(x+0.000001)),250)-90))^[1]",0,2);
  fl2.SetParameter(1,1);
  fl2.SetParameter(0,1);
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
  TNamed("AliComparisonDCA","AliComparisonDCA"),

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
  fCutsMC(0)  
{
  InitHisto();
  InitCuts();

  // vertex (0,0,0)
  fVertex = new AliESDVertex();
  fVertex->SetXv(0.0);
  fVertex->SetYv(0.0);
  fVertex->SetZv(0.0);
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
}

//_____________________________________________________________________________
void AliComparisonDCA::InitHisto()
{
  // DCA histograms
  fD0TanSPtB1 = new TH3F("DCAyTanSPtB1","DCAyTanSPt",20,1,2, 10,0.3,2, 100,-4,4);
  fD0TanSPtB1->SetXTitle("tan(#theta)");
  fD0TanSPtB1->SetYTitle("sqrt(p_{t})");
  fD0TanSPtB1->SetZTitle("DCA_{xy}");

  fD1TanSPtB1 = new TH3F("DCAzTanSPtB1","DCAzTanSPt",20,1,2, 10,0.3,2, 100,-4,4);
  fD1TanSPtB1->SetXTitle("tan(#theta)");
  fD1TanSPtB1->SetYTitle("sqrt(p_{t})");
  fD1TanSPtB1->SetZTitle("DCA_{z}");

  fD0TanSPtL1 = new TH3F("DCAyTanSPtL1","DCAyTanSPt",20,0,1, 10,0.3,2, 100,-0.1,0.1);
  fD0TanSPtL1->SetXTitle("tan(#theta)");
  fD0TanSPtL1->SetYTitle("sqrt(p_{t})");
  fD0TanSPtL1->SetZTitle("DCA_{xy}");

  fD1TanSPtL1 = new TH3F("DCAzTanSPtL1","DCAzTanSPt",20,0,1, 10,0.3,2, 100, -0.1,0.1);
  fD1TanSPtL1->SetXTitle("tan(#theta)");
  fD1TanSPtL1->SetYTitle("sqrt(p_{t})");
  fD1TanSPtL1->SetZTitle("DCA_{z}");

  fD0TanSPtInTPC = new TH3F("DCAyTanSPtInTPC","DCAyTanSPt",100,0,100, 10,0.3,2, 100,-0.1,0.1);
  fD0TanSPtInTPC->SetXTitle("tan(#theta)");
  fD0TanSPtInTPC->SetYTitle("sqrt(p_{t})");
  fD0TanSPtInTPC->SetZTitle("DCA_{xy}");

  fD1TanSPtInTPC = new TH3F("DCAzTanSPtInTPC","DCAzTanSPt",100,0,100, 10,0.3,2, 100, -0.1,0.1);
  fD1TanSPtInTPC->SetXTitle("tan(#theta)");
  fD1TanSPtInTPC->SetYTitle("sqrt(p_{t})");
  fD1TanSPtInTPC->SetZTitle("DCA_{z}");
}

//_____________________________________________________________________________
void AliComparisonDCA::InitCuts()
{
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
}

//_____________________________________________________________________________
void AliComparisonDCA::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill DCA comparison information
   
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = infoMC->GetParticle().R()<fCutsMC->GetMaxR() && TMath::Abs(infoMC->GetParticle().Vz())<fCutsMC->GetMaxVz() ;
  Float_t spt = TMath::Sqrt(mcpt);

  AliExternalTrackParam *track = 0;
  Double_t kRadius    = 3.0;      // beam pipe radius
  Double_t kMaxStep   = 5.0;      // max step
  Double_t field      = 0.4;      // mag. field
  Double_t kMaxD      = 123456.0; // max distance

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  const Double_t* dv;

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;

  // calculate and set prim. vertex
  dv = infoMC->GetVDist(); // distance to Prim. vertex

  fVertex->SetXv( infoMC->GetParticle().Vx() - dv[0] );
  fVertex->SetYv( infoMC->GetParticle().Vy() - dv[1] );
  fVertex->SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // calculate track parameters at vertex

  if (infoRC->GetESDtrack()->GetTPCInnerParam())
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      AliTracker::PropagateTrackTo(track,kRadius,infoMC->GetMass(),kMaxStep,kTRUE);
      track->PropagateToDCA(fVertex,field,kMaxD,dca,cov);

      fD0TanSPtInTPC->Fill(tantheta,spt,dca[0]);
      fD1TanSPtInTPC->Fill(tantheta,spt,dca[1]);

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
    if (entry == 0) { 
      Error("Add","Attempt to add object of class: %s to a %s",
	  obj->ClassName(),this->ClassName());
      return -1;
    }

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
  // Analyse output histograms
  
  AliComparisonDCA * comp=this;
  TGraph2D * gr=0;
  TGraph * gr0=0;

  TFile *fp = new TFile("pictures_dca.root","recreate");
  fp->cd();

  TCanvas * c = new TCanvas("DCA","DCA resloution");
  c->cd();

  // DCA resolution
  gr0 = AliMathBase::MakeStat1D(comp->fD0TanSPtB1,2,5);
  gr0->GetXaxis()->SetTitle("Tan(#theta)");
  gr0->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gr0->Write("DCAResolTan");
  //
  gr = AliMathBase::MakeStat2D(comp->fD0TanSPtB1,4,2,5); 
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gr->GetHistogram()->Write("DCAResolSPTTan");

  gPad->Clear();
  gr0->Draw("al*");

  fp->Close();
}
