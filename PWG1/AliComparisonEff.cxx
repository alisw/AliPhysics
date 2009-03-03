//------------------------------------------------------------------------------
// Implementation of AliComparisonEff class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data member of AliComparisonEff.
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  //AliComparisonEff * compObj = (AliComparisonEff*)f.Get("AliComparisonEff");
  AliComparisonEff * compObj = (AliComparisonEff*)cOutput->FindObject("AliComparisonEff");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderEff" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Eff.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include <TAxis.h>
#include <TH1D.h>

// 
#include "AliESDtrack.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
#include "AliESDVertex.h" 
#include "AliExternalTrackParam.h" 
#include "AliTracker.h" 
#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonEff.h" 

using namespace std;

ClassImp(AliComparisonEff)

//_____________________________________________________________________________
AliComparisonEff::AliComparisonEff():
  AliComparisonObject("AliComparisonEff"),

  // histograms
 
  fEffHisto(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // default consttructor	
}

//_____________________________________________________________________________
AliComparisonEff::AliComparisonEff(Char_t* name="AliComparisonEff",Char_t*title="AliComparisonEff",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliComparisonObject(name,title),

  // histograms
  fEffHisto(0),

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
AliComparisonEff::~AliComparisonEff()
{
// destructor

  if(fEffHisto)  delete  fEffHisto; fEffHisto=0;

  if(fCutsRC) delete fCutsRC; fCutsRC=0;
  if(fCutsMC) delete fCutsMC; fCutsMC=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonEff::Init()
{
  // Init histograms
  //
  Int_t nPtBins = 31;
  Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
  Double_t ptMin = 0., ptMax = 10.;

  if(IsHptGenerator() == kTRUE) {
    nPtBins = 100;
    ptMin = 0.; ptMax = 100.;
  }

  //mceta:mcphi:mcpt:pid:isPrim:recStatus:findable
  Int_t binsEffHisto[7]={30,144,nPtBins,5,2,2,2};
  Double_t minEffHisto[7]={-1.5,0.,ptMin,0.,0.,0.,0.};
  Double_t maxEffHisto[7]={ 1.5,2.*TMath::Pi(), ptMax,5.,2.,2.,2.};

  fEffHisto = new THnSparseF("fEffHisto","mceta:mcphi:mcpt:pid:isPrim:recStatus:findable",7,binsEffHisto,minEffHisto,maxEffHisto);
  if(!IsHptGenerator()) fEffHisto->SetBinEdges(2,binsPt);

  fEffHisto->GetAxis(0)->SetTitle("eta");
  fEffHisto->GetAxis(1)->SetTitle("phi (rad)");
  fEffHisto->GetAxis(2)->SetTitle("pt (GeV/c)");
  fEffHisto->GetAxis(3)->SetTitle("pid");
  fEffHisto->GetAxis(4)->SetTitle("isPrim");
  fEffHisto->GetAxis(5)->SetTitle("recStatus");
  fEffHisto->GetAxis(6)->SetTitle("findable");
  fEffHisto->Sumw2();

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderEff","Analysis Efficiency Folder");
}

//_____________________________________________________________________________
void AliComparisonEff::ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC)
{
  // Fill efficiency comparison information

  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  AliESDVertex vertexMC;  // MC primary vertex

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t mceta = infoMC->GetParticle().Eta();
  Float_t mcphi = infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();

  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
  Bool_t findable = (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());
  Bool_t recStatus = kFALSE;
 
  // calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );
  
  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 

  // transform Pdg to Pid
  // Pdg convension is different for hadrons and leptons 
  // (e.g. K+/K- = 321/-321; e+/e- = -11/11 ) 
  Double_t pid = -1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetEM() ) pid = 0; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetMuM() ) pid = 1; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP() ) pid = 2; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP() ) pid = 3; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt() ) pid = 4; 

  if (infoRC->GetESDtrack() && infoRC->GetESDtrack()->GetTPCInnerParam()) 
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0)
    {
      Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);
      if(bDCAStatus) {
        if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
        {
          recStatus = kTRUE;
        }
      }
    delete track;
    }
  }

  // Fill histograms
  Double_t vEffHisto[7] = { mceta, mcphi, mcpt, pid, isPrim, recStatus, findable }; 
  fEffHisto->Fill(vEffHisto);
}

//_____________________________________________________________________________
void AliComparisonEff::ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC)
{
  // Fill efficiency comparison information
 
  Int_t clusterITS[200];
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t mceta = infoMC->GetParticle().Eta();
  Float_t mcphi = infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
  Bool_t findable = (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());
  Bool_t recStatus =kFALSE;
 
  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 

  // transform Pdg to Pid
  // Pdg convension is different for hadrons and leptons 
  // (e.g. K+/K- = 321/-321; e+/e- = -11/11 ) 
  Double_t pid = -1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetEM() ) pid = 0; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetMuM() ) pid = 1; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP() ) pid = 2; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP() ) pid = 3; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt() ) pid = 4; 

  if(!infoRC->GetESDtrack()) return;
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS())
  {
    infoRC->GetESDtrack()->GetImpactParameters(dca,cov);
    if(TMath::Abs(dca[0]) < fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1]) < fCutsRC->GetMaxDCAToVertexZ()) 
    {
       recStatus =(infoRC->GetStatus(1)==3);
    } 
  }

  // fill histograms
  Double_t vEffHisto[7] = { mceta, mcphi, mcpt, pid, isPrim, recStatus, findable }; 
  fEffHisto->Fill(vEffHisto);
}

//_____________________________________________________________________________
void AliComparisonEff::ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC)
{
  // Fill efficiency comparison information
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  AliESDVertex vertexMC;  // MC primary vertex

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t mceta = infoMC->GetParticle().Eta();
  Float_t mcphi = infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();

  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
  Bool_t findable = (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());
  Bool_t recStatus = kFALSE;
 
  // calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );
  
  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 

  // transform Pdg to Pid
  // Pdg convension is different for hadrons and leptons 
  // (e.g. K+/K- = 321/-321; e+/e- = -11/11 ) 
  Double_t pid = -1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetEM() ) pid = 0; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetMuM() ) pid = 1; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP() ) pid = 2; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP() ) pid = 3; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt() ) pid = 4; 

  // constrained parameters resolution
  if (!infoRC->GetESDtrack()) return;
  const AliExternalTrackParam * cparam = infoRC->GetESDtrack()->GetConstrainedParam();
  if(!cparam) return;

  if ((track = new AliExternalTrackParam(*cparam)) != 0)
  {
    Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);
    if(bDCAStatus) {
      if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
      {
        recStatus =  (infoRC->GetStatus(1)!=3);
      }
    }
  delete track;
  }

  // Fill histograms
  Double_t vEffHisto[7] = { mceta, mcphi, mcpt, pid, isPrim, recStatus, findable }; 
  fEffHisto->Fill(vEffHisto);
}

//_____________________________________________________________________________
void AliComparisonEff::Exec(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC)
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
Long64_t AliComparisonEff::Merge(TCollection* const list) 
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
    AliComparisonEff* entry = dynamic_cast<AliComparisonEff*>(obj);
    if (entry == 0) continue; 
  
     fEffHisto->Add(entry->fEffHisto);
  count++;
  }

return count;
}
 
//_____________________________________________________________________________
void AliComparisonEff::Analyse() 
{
  // Analyse comparison information and store output histograms
  // in the folder "folderEff" 
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;

  // efficiency vs pt
  fEffHisto->GetAxis(0)->SetRangeUser(-0.9,0.9); // eta range
  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);    // primary tracks

  // rec efficiency vs pt
  TH1D *ptAll = fEffHisto->Projection(2);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *ptRec = fEffHisto->Projection(2);
  ptRec->Divide(ptAll);
  ptRec->SetName("ptRecEff");

  aFolderObj->Add(ptRec);


  // rec efficiency vs pid vs pt
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions
  TH1D *ptAllPi = fEffHisto->Projection(2);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecPi = fEffHisto->Projection(2);
  ptRecPi->Divide(ptAllPi);
  ptRecPi->SetName("ptRecEffPi");

  aFolderObj->Add(ptRecPi);

  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *ptAllK = fEffHisto->Projection(2);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecK = fEffHisto->Projection(2);
  ptRecK->Divide(ptAllK);
  ptRecK->SetName("ptRecEffK");

  aFolderObj->Add(ptRecK);

  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *ptAllP = fEffHisto->Projection(2);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecP = fEffHisto->Projection(2);
  ptRecP->Divide(ptAllP);
  ptRecP->SetName("ptRecEffP");

  aFolderObj->Add(ptRecP);
  
  // findable efficiency 
  fEffHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(6)->SetRangeUser(1.,1.); // findable
  TH1D *ptAllF = fEffHisto->Projection(2);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);
  fEffHisto->GetAxis(6)->SetRangeUser(1.,1.);
  TH1D *ptRecF = fEffHisto->Projection(2); // rec findable
  ptRecF->Divide(ptAllF);
  ptRecF->SetName("ptRecFindableEff");

  aFolderObj->Add(ptRecF);

  //
  // efficiency vs eta
  //

  fEffHisto->GetAxis(0)->SetRangeUser(-1.5,1.5); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.2,10.); // pt range
  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);   // primary tracks
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.);   // all
  fEffHisto->GetAxis(6)->SetRangeUser(0.,1.);   // all

  // rec efficiency vs eta
  TH1D *etaAll = fEffHisto->Projection(0);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *etaRec = fEffHisto->Projection(0);
  etaRec->Divide(etaAll);
  etaRec->SetName("etaRecEff");

  aFolderObj->Add(etaRec);

  // rec efficiency vs pid vs eta
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions
  TH1D *etaAllPi = fEffHisto->Projection(0);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecPi = fEffHisto->Projection(0);
  etaRecPi->Divide(etaAllPi);
  etaRecPi->SetName("etaRecEffPi");

  aFolderObj->Add(etaRecPi);

  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *etaAllK = fEffHisto->Projection(0);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecK = fEffHisto->Projection(0);
  etaRecK->Divide(etaAllK);
  etaRecK->SetName("etaRecEffK");

  aFolderObj->Add(etaRecK);

  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *etaAllP = fEffHisto->Projection(0);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecP = fEffHisto->Projection(0);
  etaRecP->Divide(etaAllP);
  etaRecP->SetName("etaRecEffP");

  aFolderObj->Add(etaRecP);
  
  //
  // findable efficiency 
  //
  fEffHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(6)->SetRangeUser(1.,1.); // findable
  TH1D *etaAllF = fEffHisto->Projection(0);

  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);
  fEffHisto->GetAxis(6)->SetRangeUser(1.,1.);
  TH1D *etaRecF = fEffHisto->Projection(0);   // rec findable
  etaRecF->Divide(etaAllF);
  etaRecF->SetName("etaRecFindableEff");

  aFolderObj->Add(etaRecF);

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonEff::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonEff * comp=this;
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
TFolder* AliComparisonEff::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
