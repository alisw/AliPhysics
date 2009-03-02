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
  //AliComparisonDCA * compObj = (AliComparisonDCA*)f.Get("AliComparisonDCA");
  AliComparisonDCA * compObj = (AliComparisonDCA*)cOutput->FindObject("AliComparisonDCA");

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

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TGraph.h"

#include "AliTracker.h"   
#include "AliESDEvent.h"   
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
#include "AliESDVertex.h" 
#include "AliMathBase.h"

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonDCA.h" 

using namespace std;

ClassImp(AliComparisonDCA)

//_____________________________________________________________________________
AliComparisonDCA::AliComparisonDCA():
  AliComparisonObject("AliComparisonDCA"),

  // DCA histograms
  fDCAHisto(0),
  /*
  fD0TanSPtTPCITS(0),
  fD1TanSPtTPCITS(0),
  fD0TanSPt(0),
  fD1TanSPt(0),
  fD0TanSPtTPC(0),
  fD1TanSPtTPC(0),
  */

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  // default constructor	
}

//_____________________________________________________________________________
AliComparisonDCA::AliComparisonDCA(Char_t* name="AliComparisonDCA", Char_t* title="AliComparisonDCA",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliComparisonObject(name,title),

  // DCA histograms
  fDCAHisto(0),
  /*
  fD0TanSPtTPCITS(0),
  fD1TanSPtTPCITS(0),
  fD0TanSPt(0),
  fD1TanSPt(0),
  fD0TanSPtTPC(0),
  fD1TanSPtTPC(0),
  */

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
AliComparisonDCA::~AliComparisonDCA()
{
  // destructor
  if(fDCAHisto)  delete fDCAHisto; fDCAHisto=0; 
  /*
  if(fD0TanSPtTPCITS) delete fD0TanSPtTPCITS; fD0TanSPtTPCITS=0;
  if(fD1TanSPtTPCITS) delete fD1TanSPtTPCITS; fD1TanSPtTPCITS=0;
  if(fD0TanSPt) delete fD0TanSPt; fD0TanSPt=0;
  if(fD1TanSPt) delete fD1TanSPt; fD1TanSPt=0;
  if(fD0TanSPtTPC) delete fD0TanSPtTPC; fD0TanSPtTPC=0;
  if(fD1TanSPtTPC) delete fD1TanSPtTPC; fD1TanSPtTPC=0;
  */
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;

}

//_____________________________________________________________________________
void AliComparisonDCA::Init()
{
  // DCA histograms

 Int_t nPBins = 31;
    Double_t binsP[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
    Double_t pMin = 0., pMax = 10.;

    if(IsHptGenerator() == kTRUE) {
      nPBins = 100;
      pMin = 0.; pMax = 100.;
    }

   //dca_r, dca_z, eta, pt
   Int_t binsQA[4]    = {100,100,20,nPBins};
   Double_t xminQA[4] = {-10.,-10.,-1., pMin};
   Double_t xmaxQA[4] = {10.,10.,1., pMax};

   fDCAHisto = new THnSparseF("fDCAHisto","dca_r:dca_z:eta:pt",4,binsQA,xminQA,xmaxQA);
   if(!IsHptGenerator()) fDCAHisto->SetBinEdges(3,binsP);

   fDCAHisto->GetAxis(0)->SetTitle("dca_r (cm)");
   fDCAHisto->GetAxis(1)->SetTitle("dca_z (cm)");
   fDCAHisto->GetAxis(2)->SetTitle("eta");
   fDCAHisto->GetAxis(3)->SetTitle("pt (GeV/c)");
   fDCAHisto->Sumw2();
	
  /*	
  fD0TanSPtTPCITS = new TH3F("DCAyTanSPtTPCITS","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPtTPCITS->SetXTitle("tan(#theta)");
  fD0TanSPtTPCITS->SetYTitle("#sqrt{p_{t}(GeV)}");
  fD0TanSPtTPCITS->SetZTitle("DCA_{xy}");

  fD1TanSPtTPCITS = new TH3F("DCAzTanSPtTPCITS","DCAzTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD1TanSPtTPCITS->SetXTitle("tan(#theta)");
  fD1TanSPtTPCITS->SetYTitle("#sqrt(p_{t}(GeV))");
  fD1TanSPtTPCITS->SetZTitle("DCA_{z}");

  fD0TanSPt = new TH3F("DCAyTanSPt","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPt->SetXTitle("tan(#theta)");
  fD0TanSPt->SetYTitle("#sqrt{p_{t}(GeV)}");
  fD0TanSPt->SetZTitle("DCA_{xy}");

  fD1TanSPt = new TH3F("DCAzTanSPt","DCAzTanSPt",40,-2,2, 10,0.3,3, 100, -1,1);
  fD1TanSPt->SetXTitle("tan(#theta)");
  fD1TanSPt->SetYTitle("#sqrt{p_{t}(GeV)}");
  fD1TanSPt->SetZTitle("DCA_{z}");

  fD0TanSPtTPC = new TH3F("DCAyTanSPtTPC","DCAyTanSPt",40,-2,2, 10,0.3,3, 100,-1,1);
  fD0TanSPtTPC->SetXTitle("tan(#theta)");
  fD0TanSPtTPC->SetYTitle("#sqrt{p_{t}(GeV)}");
  fD0TanSPtTPC->SetZTitle("DCA_{xy}");

  fD1TanSPtTPC = new TH3F("DCAzTanSPtTPC","DCAzTanSPt",40,-2,2, 10,0.3,3, 100, -1,1);
  fD1TanSPtTPC->SetXTitle("tan(#theta)");
  fD1TanSPtTPC->SetYTitle("#sqrt{p_{t}(GeV)}");
  fD1TanSPtTPC->SetZTitle("DCA_{z}");
  */

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
 
  // init folder
  fAnalysisFolder = CreateFolder("folderDCA","Analysis DCA Folder");
}

//_____________________________________________________________________________
void AliComparisonDCA::ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo * const infoRC)
{
  // Fill DCA comparison information
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance

  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t mceta = infoMC->GetParticle().Eta();
  //Float_t spt = TMath::Sqrt(mcpt);

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;
  //if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;

  // calculate and set prim. vertex
  AliESDVertex vertexMC;
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // calculate track parameters at vertex
  if (infoRC->GetESDtrack()->GetTPCInnerParam())
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);

      if(bDCAStatus) {
	 Double_t vDCAHisto[4]={dca[0],dca[1],mceta,mcpt};
	 fDCAHisto->Fill(vDCAHisto);
      }
    delete track;
    }
  }
}

//_____________________________________________________________________________
void AliComparisonDCA::ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo * const infoRC)
{
  // Fill DCA comparison information
  Int_t clusterITS[200];
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t mceta = infoMC->GetParticle().Eta();
  //Float_t spt = TMath::Sqrt(mcpt);

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;
  //if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;

  infoRC->GetESDtrack()->GetImpactParameters(dca,cov);

  // ITS + TPC
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS())
  {
    Double_t vDCAHisto[4]={dca[0],dca[1],mceta,mcpt};
    fDCAHisto->Fill(vDCAHisto);
  }
}

void AliComparisonDCA::ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo * const infoRC)
{
  // Fill DCA comparison information
  
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}

//_____________________________________________________________________________
Long64_t AliComparisonDCA::Merge(TCollection* const list) 
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

    fDCAHisto->Add(entry->fDCAHisto);
    /*
    fD0TanSPtTPCITS->Add(entry->fD0TanSPtTPCITS);
    fD1TanSPtTPCITS->Add(entry->fD1TanSPtTPCITS);
    fD0TanSPt->Add(entry->fD0TanSPt);
    fD1TanSPt->Add(entry->fD1TanSPt);
    fD0TanSPtTPC->Add(entry->fD0TanSPtTPC);
    fD1TanSPtTPC->Add(entry->fD1TanSPtTPC);
    */

    count++;
  }

return count;
}

//_____________________________________________________________________________
void AliComparisonDCA::Exec(AliMCInfo* const infoMC, AliESDRecInfo * const infoRC)
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
void AliComparisonDCA::Analyse()
{
  //
  // Analyse comparison information and store output histograms
  // in the analysis folder "folderDCA" 
  //
  
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;

  /*
  TGraph * gr[8]= { 0,0,0,0,0,0,0,0 };
  TGraph2D *gr2[8]= { 0,0,0,0,0,0,0,0};
  AliComparisonDCA * comp=this;

  // write results in the folder 
  // Canvas to draw analysed histograms
  TCanvas * c = new TCanvas("canDCA","DCA resolution");
  c->Divide(2,4);
  //
  // DCA resolution
  //
  c->cd(1);
  gr[0] = AliMathBase::MakeStat1D(comp->fD0TanSPtTPC,2,5);
  gr[0]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[0]->GetYaxis()->SetTitle("#sigmaDCA_xy (cm)");
  gr[0]->SetName("DCAXYResolTanTPC");
  gr[0]->SetTitle("resol. DCA_xy (TPC only)");
  gr[0]->Draw("Al*");

  aFolderObj->Add(gr[0]);

  c->cd(2);
  gr[1] = AliMathBase::MakeStat1D(comp->fD1TanSPtTPC,2,5);
  gr[1]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[1]->GetYaxis()->SetTitle("#sigmaDCA_z (cm)");
  gr[1]->SetName("DCAZResolTanTPC");
  gr[1]->SetTitle("resol. DCA_z (TPC only)");
  gr[1]->Draw("Al*");

  aFolderObj->Add(gr[1]);

  c->cd(3);
  gr[2] = AliMathBase::MakeStat1D(comp->fD0TanSPtTPCITS,2,5);
  gr[2]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[2]->GetYaxis()->SetTitle("#sigmaDCA_xy (cm)");
  gr[2]->SetName("DCAXYResolTanTPCITS");
  gr[2]->SetTitle("resol. DCA_xy (TPC+ITS)");
  gr[2]->Draw("Al*");

  aFolderObj->Add(gr[2]);

  c->cd(4);
  gr[3] = AliMathBase::MakeStat1D(comp->fD1TanSPtTPCITS,2,5);
  gr[3]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[3]->GetYaxis()->SetTitle("#sigmaDCA_z (cm)");
  gr[3]->SetName("DCAZResolTanTPCITS");
  gr[3]->SetTitle("resol. DCA_z (TPC+ITS)");
  gr[3]->Draw("Al*");

  aFolderObj->Add(gr[3]);

  //
  // DCA mean value
  //
  c->cd(5);
  gr[4] = AliMathBase::MakeStat1D(comp->fD0TanSPtTPC,2,4);
  gr[4]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[4]->GetYaxis()->SetTitle("mean DCA_xy (cm)");
  gr[4]->SetName("DCAXYMeanTanTPC");
  gr[4]->SetTitle("mean DCA_xy (TPC only)");
  gr[4]->Draw("Al*");

  aFolderObj->Add(gr[4]);

  c->cd(6);
  gr[5] = AliMathBase::MakeStat1D(comp->fD1TanSPtTPC,2,4);
  gr[5]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[5]->GetYaxis()->SetTitle("mean DCA_z (cm)");
  gr[5]->SetName("DCAZMeanTanTPC");
  gr[5]->SetTitle("mean DCA_z (TPC only)");
  gr[5]->Draw("Al*");

  aFolderObj->Add(gr[5]);

  c->cd(7);
  gr[6] = AliMathBase::MakeStat1D(comp->fD0TanSPtTPCITS,2,4);
  gr[6]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[6]->GetYaxis()->SetTitle("mean DCA_xy (cm)");
  gr[6]->SetName("DCAXYMeanTanTPCITS");
  gr[6]->SetTitle("mean DCA_xy (TPC+ITS)");
  gr[6]->Draw("Al*");

  aFolderObj->Add(gr[6]);

  c->cd(8);
  gr[7] = AliMathBase::MakeStat1D(comp->fD1TanSPtTPCITS,2,4);
  gr[7]->GetXaxis()->SetTitle("Tan(#theta)");
  gr[7]->GetYaxis()->SetTitle("mean DCA_z (cm)");
  gr[7]->SetName("DCAZMeanTanTPCITS");
  gr[7]->SetTitle("mean DCA_z (TPC+ITS)");
  gr[7]->Draw("Al*");

  aFolderObj->Add(gr[7]);

  // 2D DCA resolution 
  TCanvas * c1 = new TCanvas("canDCA1","2D DCA resolution");
  c1->Divide(2,4);

  // TPC only
  c1->cd(1);
  gr2[0] = AliMathBase::MakeStat2D(comp->fD0TanSPtTPC,4,2,5); 
  gr2[0]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[0]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[0]->GetZaxis()->SetTitle("#sigmaDCA_xy (cm)");
  gr2[0]->SetName("DCAXYResolSPTTanTPC");
  gr2[0]->SetTitle("#sigma DCA_xy (TPC only)");
  gr2[0]->GetHistogram()->Draw("colz");

  gr2[0]->GetHistogram()->SetName("DCAXYResolSPTTanTPC");
  aFolderObj->Add(gr2[0]->GetHistogram());

  c1->cd(2);
  gr2[1] = AliMathBase::MakeStat2D(comp->fD1TanSPtTPC,4,2,5); 
  gr2[1]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[1]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[1]->GetZaxis()->SetTitle("#sigmaDCA_z (cm)");
  gr2[1]->SetName("DCAZResolSPTTanTPC");
  gr2[1]->SetTitle("#sigma DCA_z (TPC only)");
  gr2[1]->GetHistogram()->Draw("colz");

  gr2[1]->GetHistogram()->SetName("DCAZResolSPTTanTPC");
  aFolderObj->Add(gr2[1]->GetHistogram());

  // TPC+ITS
  c1->cd(3);
  gr2[2] = AliMathBase::MakeStat2D(comp->fD0TanSPtTPCITS,4,2,5); 
  gr2[2]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[2]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[2]->GetZaxis()->SetTitle("#sigmaDCA_xy (cm)");
  gr2[2]->SetName("DCAXYResolSPTTanTPCITS");
  gr2[2]->SetTitle("#sigma DCA_xy (TPC+ITS)");
  gr2[2]->GetHistogram()->Draw("colz");

  gr2[2]->GetHistogram()->SetName("DCAXYResolSPTTanTPCITS");
  aFolderObj->Add(gr2[2]->GetHistogram());

  c1->cd(4);
  gr2[3] = AliMathBase::MakeStat2D(comp->fD1TanSPtTPCITS,4,2,5); 
  gr2[3]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[3]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[3]->GetZaxis()->SetTitle("#sigmaDCA_z (cm)");
  gr2[3]->SetName("DCAZResolSPTTanTPCITS");
  gr2[3]->SetTitle("#sigma DCA_z (TPC+ITS)");
  gr2[3]->GetHistogram()->Draw("colz");

  gr2[3]->GetHistogram()->SetName("DCAZResolSPTTanTPCITS");
  aFolderObj->Add(gr2[3]->GetHistogram());

  // 2D DCA mean value  
  c1->cd(5);
  gr2[4] = AliMathBase::MakeStat2D(comp->fD0TanSPtTPC,4,2,4); 
  gr2[4]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[4]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[4]->GetZaxis()->SetTitle("mean DCA_xy (cm)");
  gr2[4]->SetName("DCAXYMeanSPTTanTPC");
  gr2[4]->SetTitle("mean DCA_xy (TPC only)");
  gr2[4]->GetHistogram()->Draw("colz");

  gr2[4]->GetHistogram()->SetName("DCAXYMeanSPTTanTPC");
  aFolderObj->Add(gr2[4]->GetHistogram());

  c1->cd(6);
  gr2[5] = AliMathBase::MakeStat2D(comp->fD1TanSPtTPC,4,2,4); 
  gr2[5]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[5]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[5]->GetZaxis()->SetTitle("mean DCA_z (cm)");
  gr2[5]->SetName("DCAZMeanSPTTanTPC");
  gr2[5]->SetTitle("mean DCA_z (TPC only)");
  gr2[5]->GetHistogram()->Draw("colz");

  gr2[5]->GetHistogram()->SetName("DCAZMeanSPTTanTPC");
  aFolderObj->Add(gr2[5]->GetHistogram());

  c1->cd(7);
  gr2[6] = AliMathBase::MakeStat2D(comp->fD0TanSPtTPCITS,4,2,4); 
  gr2[6]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[6]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[6]->GetZaxis()->SetTitle("mean DCA_xy (cm)");
  gr2[6]->SetName("DCAXYMeanSPTTanTPCITS");
  gr2[6]->SetTitle("mean DCA_xy (TPC+ITS)");
  gr2[6]->GetHistogram()->Draw("colz");

  gr2[6]->GetHistogram()->SetName("DCAXYMeanSPTTanTPCITS");
  aFolderObj->Add(gr2[6]->GetHistogram());

  c1->cd(8);
  gr2[7] = AliMathBase::MakeStat2D(comp->fD1TanSPtTPCITS,4,2,4); 
  gr2[7]->GetXaxis()->SetTitle("Tan(#theta)");
  gr2[7]->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr2[7]->GetZaxis()->SetTitle("mean DCA_z (cm)");
  gr2[7]->SetName("DCAZMeanSPTTanTPCITS");
  gr2[7]->SetTitle("mean DCA_z (TPC+ITS)");
  gr2[7]->GetHistogram()->Draw("colz");

  gr2[7]->GetHistogram()->SetName("DCAZMeanSPTTanTPCITS");
  aFolderObj->Add(gr2[7]->GetHistogram());

  */
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
