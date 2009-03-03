//------------------------------------------------------------------------------
// Implementation of AliComparisonRes class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliComparisonRes.
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  //AliComparisonRes * compObj = (AliComparisonRes*)f.Get("AliComparisonRes");
  AliComparisonRes * compObj = (AliComparisonRes*)cOutput->FindObject("AliComparisonRes");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderRes" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Res.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TAxis.h"

#include "AliComparisonRes.h" 
#include "AliESDRecInfo.h" 
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h" 
#include "AliMCInfo.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTreeDraw.h" 

using namespace std;

ClassImp(AliComparisonRes)

//_____________________________________________________________________________
AliComparisonRes::AliComparisonRes():
  AliComparisonObject("AliComparisonRes"),
  fResolHisto(0),
  fPullHisto(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  //Init();
}

//_____________________________________________________________________________
AliComparisonRes::AliComparisonRes(Char_t* name="AliComparisonRes", Char_t* title="AliComparisonRes",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE):
//AliComparisonRes::AliComparisonRes(Char_t* name, Char_t* title,Int_t analysisMode,Bool_t hptGenerator):
  AliComparisonObject(name,title),
  fResolHisto(0),
  fPullHisto(0),

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
AliComparisonRes::~AliComparisonRes()
{
  // destructor
   
  if(fResolHisto) delete fResolHisto; fResolHisto=0;     
  if(fPullHisto)  delete fPullHisto;  fPullHisto=0;     

  if(fCutsRC) delete fCutsRC; fCutsRC=0;
  if(fCutsMC) delete fCutsMC; fCutsMC=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonRes::Init(){

  //
  // histogram bining
  //

  // set pt bins
   Int_t nPtBins = 31;
   Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
  Double_t ptMin = 0., ptMax = 10.; 


  if(IsHptGenerator() == kTRUE) {
    nPtBins = 100;
    ptMin = 0.; ptMax = 100.; 
  }

  // res_y:res_z:res_phi,res_lambda:res_1pt:y:z:eta:phi:pt
  Int_t binsResolHisto[10]={100,100,100,100,100,50,100,30,144,nPtBins};
  Double_t minResolHisto[10]={-1.,-1.,-0.03,-0.03,-0.2,-0.03,-20.,-1.5,0.,ptMin};
  Double_t maxResolHisto[10]={ 1., 1., 0.03, 0.03, 0.2, 0.03, 20., 1.5,2.*TMath::Pi(), ptMax};


  fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_1pt:y:z:eta:phi:pt",10,binsResolHisto,minResolHisto,maxResolHisto);
  if(!IsHptGenerator()) fResolHisto->SetBinEdges(9,binsPt);

  fResolHisto->GetAxis(0)->SetTitle("res_y (cm)");
  fResolHisto->GetAxis(1)->SetTitle("res_z (cm)");
  fResolHisto->GetAxis(2)->SetTitle("res_phi (rad)");
  fResolHisto->GetAxis(3)->SetTitle("res_lambda (rad)");
  fResolHisto->GetAxis(4)->SetTitle("res_1pt (GeV/c)^{-1}");
  fResolHisto->GetAxis(5)->SetTitle("y (cm)");
  fResolHisto->GetAxis(6)->SetTitle("z (cm)");
  fResolHisto->GetAxis(7)->SetTitle("eta");
  fResolHisto->GetAxis(8)->SetTitle("phi (rad)");
  fResolHisto->GetAxis(9)->SetTitle("pt (GeV/c)");
  fResolHisto->Sumw2();

  //pull_y:pull_z:pull_lambda:pull_phi:pull_1pt:y:z:eta:phi:pt
  Int_t binsPullHisto[10]={100,100,100,100,100,50,100,30,144,nPtBins};
  Double_t minPullHisto[10]={-5.,-5.,-5.,-5.,-5.,-0.03, -20.,-1.5, 0., ptMin};
  Double_t maxPullHisto[10]={ 5., 5., 5., 5., 5., 0.03, 20., 1.5, 2.*TMath::Pi(),ptMax};

  fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_phi:pull_lambda:pull_1pt:y:z:eta:phi:pt",10,binsPullHisto,minPullHisto,maxPullHisto);
  if(!IsHptGenerator()) fPullHisto->SetBinEdges(9,binsPt);
  
  fPullHisto->GetAxis(0)->SetTitle("res_y (cm)");
  fPullHisto->GetAxis(1)->SetTitle("res_z (cm)");
  fPullHisto->GetAxis(2)->SetTitle("res_phi (rad)");
  fPullHisto->GetAxis(3)->SetTitle("res_lambda (rad)");
  fPullHisto->GetAxis(4)->SetTitle("res_1pt (GeV/c)^{-1}");
  fPullHisto->GetAxis(5)->SetTitle("y (rad)");
  fPullHisto->GetAxis(6)->SetTitle("z (rad)");
  fPullHisto->GetAxis(7)->SetTitle("eta");
  fPullHisto->GetAxis(8)->SetTitle("phi (rad)");
  fPullHisto->GetAxis(9)->SetTitle("pt (GeV/c)");
  fPullHisto->Sumw2();

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderRes","Analysis Resolution Folder");
}

//_____________________________________________________________________________
void AliComparisonRes::ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Fill resolution comparison information 
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  AliESDVertex vertexMC;  // MC primary vertex

  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t delta1PtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 
  //Float_t delta1Pt2TPC, deltaY1PtTPC, deltaZ1PtTPC, deltaPhi1PtTPC, deltaTheta1PtTPC; 

  Float_t mceta =  infoMC->GetParticle().Eta();
  Float_t mcphi =  infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = infoMC->GetParticle().Pt();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  //if (infoRC->GetStatus(1)!=3) return; // TPC refit
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  // calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // calculate track parameters at vertex
  const AliExternalTrackParam *innerTPC =  0;
  if ((innerTPC = infoRC->GetESDtrack()->GetTPCInnerParam()) != 0)
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);

      // Fill parametrisation histograms (only TPC track)
      if(bDCAStatus) 
      {
          // select primaries
          if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ()) 
	  { 

            deltaYTPC= track->GetY()-infoMC->GetParticle().Vy();
            deltaZTPC = track->GetZ()-infoMC->GetParticle().Vz();
            deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
            deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
	    delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;

            pullYTPC= (track->GetY()-infoMC->GetParticle().Vy()) / TMath::Sqrt(track->GetSigmaY2());
            pullZTPC = (track->GetZ()-infoMC->GetParticle().Vz()) / TMath::Sqrt(track->GetSigmaZ2());
            pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
            pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 
	    pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());

            Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,delta1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};
	    fResolHisto->Fill(vResolHisto);

            Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};
	    fPullHisto->Fill(vPullHisto);
	  }

          /*
            delta1Pt2TPC = (1/mcpt-innerTPC->OneOverPt())/TMath::Power(1+1/mcpt,2);       
            deltaY1PtTPC= (infoMC->GetParticle().Vy()-innerTPC->GetY()) / (0.2+1/mcpt);
            deltaZ1PtTPC = (infoMC->GetParticle().Vz()-innerTPC->GetZ()) / (0.2+1/mcpt);
            deltaPhi1PtTPC = deltaPhiTPC   / (0.1+1/mcpt);
            deltaTheta1PtTPC = deltaLambdaTPC / (0.1+1/mcpt);

	    fPhiEtaPtTPC->Fill(TMath::ATan2(innerTPC->Py(),innerTPC->Px()), innerTPC->Eta(), innerTPC->Pt()); 

            f1Pt2ResolS1PtTPC->Fill(s1mcpt,delta1Pt2TPC);
            fYResolS1PtTPC->Fill(s1mcpt,deltaY1PtTPC);
            fZResolS1PtTPC->Fill(s1mcpt,deltaZ1PtTPC);
            fPhiResolS1PtTPC->Fill(s1mcpt,deltaPhi1PtTPC);
            fThetaResolS1PtTPC->Fill(s1mcpt,deltaTheta1PtTPC);

            fPtResolPtTPC->Fill(mcpt,deltaPtTPC);
            fPtPullPtTPC->Fill(mcpt,pullPtTPC);
            fPhiResolEtaTPC->Fill(eta,deltaPhiTPC);
            fPhiPullEtaTPC->Fill(eta,pullPhiTPC);
            fLambdaResolEtaTPC->Fill(eta,deltaLambdaTPC);
            fLambdaPullEtaTPC->Fill(eta,pullLambdaTPC);
          */
	}
    delete track;
    }
  }
}

//_____________________________________________________________________________
void AliComparisonRes::ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  Int_t clusterITS[200];
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t delta1PtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  Float_t mceta =  infoMC->GetParticle().Eta();
  Float_t mcphi =  infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = infoMC->GetParticle().Pt();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return; // TPC refit
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  infoRC->GetESDtrack()->GetImpactParameters(dca,cov);

  // select primaries
  if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ()) 
  { 
    deltaYTPC= infoRC->GetESDtrack()->GetY()-infoMC->GetParticle().Vy();
    deltaZTPC = infoRC->GetESDtrack()->GetZ()-infoMC->GetParticle().Vz();
    deltaLambdaTPC = TMath::ATan2(infoRC->GetESDtrack()->Pz(),infoRC->GetESDtrack()->Pt())-TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
    deltaPhiTPC = TMath::ATan2(infoRC->GetESDtrack()->Py(),infoRC->GetESDtrack()->Px())-TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
    delta1PtTPC = (infoRC->GetESDtrack()->OneOverPt()-1./mcpt)*mcpt;

    pullYTPC= (infoRC->GetESDtrack()->GetY()-infoMC->GetParticle().Vy()) / TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaY2());
    pullZTPC = (infoRC->GetESDtrack()->GetZ()-infoMC->GetParticle().Vz()) / TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaZ2());
    pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaTgl2());
    pullPhiTPC = deltaPhiTPC / TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaSnp2()); 
    pull1PtTPC = (infoRC->GetESDtrack()->OneOverPt()-1./mcpt) / TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,delta1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};
    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};

    // TPC and ITS clusters in the system
    if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS()) 
    {
      fResolHisto->Fill(vResolHisto);
      fPullHisto->Fill(vPullHisto);
    }
  }
}

//_____________________________________________________________________________
void AliComparisonRes::ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC)
{
  // Fill resolution comparison information (constarained parameters) 
  //
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  AliESDVertex vertexMC;  // MC primary vertex

  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t delta1PtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  Float_t mceta =  infoMC->GetParticle().Eta();
  Float_t mcphi =  infoMC->GetParticle().Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = infoMC->GetParticle().Pt();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return; // TPC refit
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  // constrained parameters resolution
  const AliExternalTrackParam * cparam = infoRC->GetESDtrack()->GetConstrainedParam();
  if(!cparam) return;

  if ((track = new AliExternalTrackParam(*cparam)) != 0)
  {
    Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);
    if(bDCAStatus) {
      if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
      {
        deltaYTPC= track->GetY()-infoMC->GetParticle().Vy();
        deltaZTPC = track->GetZ()-infoMC->GetParticle().Vz();
        deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
        deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
        delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;

        pullYTPC= (track->GetY()-infoMC->GetParticle().Vy()) / TMath::Sqrt(track->GetSigmaY2());
        pullZTPC = (track->GetZ()-infoMC->GetParticle().Vz()) / TMath::Sqrt(track->GetSigmaZ2());
        pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
        pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 
        pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());

        Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,delta1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};
        fResolHisto->Fill(vResolHisto);

        Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,infoMC->GetParticle().Vy(),infoMC->GetParticle().Vz(),mceta,mcphi,mcpt};
        fPullHisto->Fill(vPullHisto);
      }
    }
  delete track;
  }
}
 
//_____________________________________________________________________________
void AliComparisonRes::Exec(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC){
  
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
TH1F* AliComparisonRes::MakeResol(TH2F * his, Int_t integ, Bool_t type){
  // Create resolution histograms
 
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoI(his,&hism,integ);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliComparisonRes::Analyse(){
  // Analyse comparison information and store output histograms
  // in the folder "folderRes"
  //
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH2F *h2D=0;
  TObjArray *aFolderObj = new TObjArray;

  // write results in the folder 
  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  c->cd();

  char name[256];
  char res_title[256] = {"res_y:res_z:res_lambda:res_phi:res_1pt:y:z:eta:phi:pt"} ;
  char pull_title[256] = {"pull_y:pull_z:pull_lambda:pull_phi:pull_1pt:y:z:eta:phi:pt"};

  for(Int_t i=0; i<5; i++) 
  {
    for(Int_t j=5; j<10; j++) 
    {
      if(j==7) fResolHisto->GetAxis(9)->SetRangeUser(0.2,10.); // pt threshold
      if(j==9) fResolHisto->GetAxis(7)->SetRangeUser(-0.9,0.9); // eta window

      h2D = (TH2F*)fResolHisto->Projection(i,j);
      h = AliComparisonRes::MakeResol(h2D,1,0);
      sprintf(name,"h_res_%d_vs_%d",i,j);
      h->SetName(name);
      h->SetTitle(res_title);

      aFolderObj->Add(h);

      h = AliComparisonRes::MakeResol(h2D,1,1);
      sprintf(name,"h_mean_res_%d_vs_%d",i,j);
      h->SetName(name);
      h->SetTitle(res_title);

      aFolderObj->Add(h);

      if(j==7) fResolHisto->GetAxis(9)->SetRangeUser(0.0,10.);
      if(j==9) fResolHisto->GetAxis(7)->SetRangeUser(-1.5,1.5);

      //
      if(j==7) fPullHisto->GetAxis(9)->SetRangeUser(0.2,10.);
      if(j==9) fPullHisto->GetAxis(7)->SetRangeUser(-0.9,0.9);

      h2D = (TH2F*)fPullHisto->Projection(i,j);
      h = AliComparisonRes::MakeResol(h2D,1,0);
      sprintf(name,"h_pull_%d_vs_%d",i,j);
      h->SetName(name);
      h->SetTitle(pull_title);

      aFolderObj->Add(h);

      h = AliComparisonRes::MakeResol(h2D,1,1);
      sprintf(name,"h_mean_pull_%d_vs_%d",i,j);
      h->SetName(name);
      h->SetTitle(pull_title);

      aFolderObj->Add(h);

      if(j==7) fPullHisto->GetAxis(9)->SetRangeUser(0.0,10.);
      if(j==9) fPullHisto->GetAxis(7)->SetRangeUser(-1.5,1.5);
    }
  }

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonRes::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonRes * comp=this;
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
Long64_t AliComparisonRes::Merge(TCollection* const list) 
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
  AliComparisonRes* entry = dynamic_cast<AliComparisonRes*>(obj);
  if (entry == 0) continue; 

  fResolHisto->Add(entry->fResolHisto);
  fPullHisto->Add(entry->fPullHisto);

  count++;
  }

return count;
}

//_____________________________________________________________________________
TFolder* AliComparisonRes::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
