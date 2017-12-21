//------------------------------------------------------------------------------
// Implementation of AliPerformanceRes class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformanceRes.
//
// Author: J.Otwinowski 04/02/2008
// Changes by J.Salzwedel 23/10/2014
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceRes * compObj = (AliPerformanceRes*)coutput->FindObject("AliPerformanceRes");
 
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
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"

#include "AliPerformanceRes.h" 
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliVfriendTrack.h"
#include "AliVfriendEvent.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTreeDraw.h" 

using namespace std;

ClassImp(AliPerformanceRes)
Double_t AliPerformanceRes::fgkMergeEntriesCut=5000000.; //5*10**6 tracks (small default to keep default memory foorprint low)

//_____________________________________________________________________________
AliPerformanceRes::AliPerformanceRes(TRootIOCtor* b):
  AliPerformanceObject(b),
  fResolHisto(0),
  fPullHisto(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // io constructor	
}

//_____________________________________________________________________________
AliPerformanceRes::AliPerformanceRes(const Char_t* name, const Char_t* title, Int_t analysisMode, Bool_t hptGenerator):
  AliPerformanceObject(name,title),
  fResolHisto(0),
  fPullHisto(0),

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
AliPerformanceRes::~AliPerformanceRes()
{
  // destructor
   
  if(fResolHisto) delete fResolHisto; fResolHisto=0;     
  if(fPullHisto)  delete fPullHisto;  fPullHisto=0;     
  
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceRes::Init(){

  //
  // histogram bining
  //

  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 15.e-3, ptMax = 20.;

  Double_t *binsPt = 0;

  if (IsHptGenerator())  { 
        ptMax = 100.;
  } 
  binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);

  Double_t yMin = -0.02, yMax = 0.02;
  Double_t zMin = -12.0, zMax = 12.0;
  if(GetAnalysisMode() == 3) { // TrackRef coordinate system
    yMin = -100.; yMax = 100.; 
    zMin = -100.; zMax = 100.; 
  }

  // res_y:res_z:res_phi,res_lambda:res_pt:y:z:eta:phi:pt
  Int_t binsResolHisto[10]={100,100,100,100,100,25,50,144,30,nPtBins};
  Double_t minResolHisto[10]={-1.,-1.,-0.03,-0.03,-0.2, yMin, zMin, 0., -1.5, ptMin};
  Double_t maxResolHisto[10]={ 1., 1., 0.03, 0.03, 0.2, yMax, zMax, 2.*TMath::Pi(), 1.5, ptMax};

  fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_pt:y:z:phi:eta:pt",10,binsResolHisto,minResolHisto,maxResolHisto);

  //fResolHisto->SetBinEdges(9,binsPt);
  fResolHisto->GetAxis(9)->Set(nPtBins,binsPt);

  fResolHisto->GetAxis(0)->SetTitle("y-y_{mc} (cm)");
  fResolHisto->GetAxis(1)->SetTitle("z-z_{mc} (cm)");
  fResolHisto->GetAxis(2)->SetTitle("#phi-#phi_{mc} (rad)");
  fResolHisto->GetAxis(3)->SetTitle("#lambda-#lambda_{mc} (rad)");
  fResolHisto->GetAxis(4)->SetTitle("(p_{T}/p_{Tmc}-1)");
  fResolHisto->GetAxis(5)->SetTitle("y_{mc} (cm)");
  fResolHisto->GetAxis(6)->SetTitle("z_{mc} (cm)");
  fResolHisto->GetAxis(7)->SetTitle("#phi_{mc} (rad)");
  fResolHisto->GetAxis(8)->SetTitle("#eta_{mc}");
  fResolHisto->GetAxis(9)->SetTitle("p_{Tmc} (GeV/c)");
  fResolHisto->Sumw2();

  ////pull_y:pull_z:pull_phi:pull_lambda:pull_1pt:y:z:eta:phi:pt
  //Int_t binsPullHisto[10]={100,100,100,100,100,50,50,30,144,nPtBins};
  //Double_t minPullHisto[10]={-5.,-5.,-5.,-5.,-5.,yMin, zMin,-1.5, 0., ptMin};
  //Double_t maxPullHisto[10]={ 5., 5., 5., 5., 5., yMax, zMax, 1.5, 2.*TMath::Pi(),ptMax};
  //fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_phi:pull_lambda:pull_1pt:y:z:eta:phi:pt",10,binsPullHisto,minPullHisto,maxPullHisto);

  //pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt
  //Int_t binsPullHisto[10]={100,100,100,100,100,50,50,50,50,nPtBins};
  //Double_t minPullHisto[10]={-5.,-5.,-5.,-5.,-5.,yMin, zMin,-1., -2.0, ptMin};
  //Double_t maxPullHisto[10]={ 5., 5., 5., 5., 5., yMax, zMax, 1., 2.0, ptMax};
  Int_t binsPullHisto[10]={100,100,100,100,100,50,50,50,50,20};
  Double_t minPullHisto[10]={-5.,-5.,-5.,-5.,-5.,yMin, zMin,-1., -2.0, 0.};
  Double_t maxPullHisto[10]={ 5., 5., 5., 5., 5., yMax, zMax, 1., 2.0, 10.};
  fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt",10,binsPullHisto,minPullHisto,maxPullHisto);

  /*
  if(!IsHptGenerator()) fPullHisto->SetBinEdges(9,bins1Pt);
  fPullHisto->GetAxis(0)->SetTitle("(y-y_{mc})/#sigma");
  fPullHisto->GetAxis(1)->SetTitle("(z-z_{mc})/#sigma");
  fPullHisto->GetAxis(2)->SetTitle("(#phi-#phi_{mc})/#sigma");
  fPullHisto->GetAxis(3)->SetTitle("(#lambda-#lambda_{mc})/#sigma");
  fPullHisto->GetAxis(4)->SetTitle("(p_{Tmc}/p_{T}-1)/#sigma");
  fPullHisto->GetAxis(5)->SetTitle("y_{mc} (cm)");
  fPullHisto->GetAxis(6)->SetTitle("z_{mc} (cm)");
  fPullHisto->GetAxis(7)->SetTitle("#eta_{mc}");
  fPullHisto->GetAxis(8)->SetTitle("#phi_{mc} (rad)");
  fPullHisto->GetAxis(9)->SetTitle("p_{Tmc} (GeV/c)");
  fPullHisto->Sumw2();
  */

  fPullHisto->GetAxis(9)->Set(nPtBins,binsPt);

  fPullHisto->GetAxis(0)->SetTitle("(y-y_{mc})/#sigma");
  fPullHisto->GetAxis(1)->SetTitle("(z-z_{mc})/#sigma");
  fPullHisto->GetAxis(2)->SetTitle("(sin#phi-sin#phi_{mc})/#sigma");
  fPullHisto->GetAxis(3)->SetTitle("(tan#lambda-tan#lambda_{mc})/#sigma");
  fPullHisto->GetAxis(4)->SetTitle("(p_{Tmc}/p_{T}-1)/#sigma");
  fPullHisto->GetAxis(5)->SetTitle("y_{mc} (cm)");
  fPullHisto->GetAxis(6)->SetTitle("z_{mc} (cm)");
  fPullHisto->GetAxis(7)->SetTitle("sin#phi_{mc}");
  fPullHisto->GetAxis(8)->SetTitle("tan#lambda_{mc}");
  fPullHisto->GetAxis(9)->SetTitle("1/p_{Tmc} (GeV/c)^{-1}");
  fPullHisto->Sumw2();

  // init folder
  fAnalysisFolder = CreateFolder("folderRes","Analysis Resolution Folder");
}

//_____________________________________________________________________________
void AliPerformanceRes::ProcessTPC(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent)
{
  if(!vEvent || !vTrack) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexTPCBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }

  // Fill TPC only resolution comparison information 
  const AliExternalTrackParam* tmpTrack = vTrack->GetTPCInnerParam();
  if(!tmpTrack) return;

  AliExternalTrackParam track = *tmpTrack;

  Float_t dca[2] = {0.,0.};// dca_xy, dca_z, 
  Float_t cov[3] = {0.,0.,0.}; //sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParametersTPC(dca,cov);
 
  //
  // Fill rec vs MC information
  //
  if(!mcev) return;
  Int_t label = vTrack->GetTPCLabel(); //Use TPC-only label for TPC-only resolution analysis
  if (label <= 0) return;
  TParticle* particle = ((AliMCParticle*)mcev->GetTrack(label))->Particle();
  if(!particle) return;
  if(!particle->GetPDG()) return;
  if(particle->GetPDG()->Charge()==0) return;
  //printf("charge %d \n",particle->GetPDG()->Charge());
  
  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) return;

  // exclude electrons
  if (fCutsMC.GetEM()==TMath::Abs(particle->GetPdgCode())) return;

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  Float_t mceta =  particle->Eta();
  Float_t mcphi =  particle->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = particle->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(particle->Py(),particle->Px())); 
  Float_t mctgl = TMath::Tan(TMath::ATan2(particle->Pz(),particle->Pt()));

  // nb. TPC clusters cut
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return;

  // select primaries
  if(TMath::Abs(dca[0])<fCutsRC.GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC.GetMaxDCAToVertexZ()) 
  { 
    if(mcpt == 0) return;
    double Bz = vEvent->GetMagneticField();

    Double_t mclocal[4]; //Rotated x,y,px,py mc-coordinates - the MC data should be rotated since the track is propagated best along x
    Double_t c = TMath::Cos(track.GetAlpha());
    Double_t s = TMath::Sin(track.GetAlpha());
    Double_t x = particle->Vx();
    Double_t y = particle->Vy();
    mclocal[0] = x*c + y*s;
    mclocal[1] =-x*s + y*c;
    Double_t px = particle->Px();
    Double_t py = particle->Py();
    mclocal[2] = px*c + py*s;
    mclocal[3] =-px*s + py*c;
    Float_t mcsnplocal = TMath::Sin(TMath::ATan2(mclocal[3],mclocal[2])); 


    track.AliExternalTrackParam::PropagateTo(mclocal[0],Bz);

    deltaYTPC= track.GetY()-mclocal[1];
    deltaZTPC = track.GetZ()-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(track.Pz(),track.Pt())-TMath::ATan2(particle->Pz(),particle->Pt());
	//See comments in ProcessInnerTPC for remarks on local and global momentum coordinates for deltaPhi / pullSnp calculation
    deltaPhiTPC = TMath::ATan2(track.Py(),track.Px())-TMath::ATan2(particle->Py(),particle->Px());
    //delta1PtTPC = (track.OneOverPt()-1./mcpt)*mcpt;
    deltaPtTPC = (track.Pt()-mcpt) / mcpt;

    pullYTPC= deltaYTPC / TMath::Sqrt(track.GetSigmaY2());
    pullZTPC = deltaZTPC / TMath::Sqrt(track.GetSigmaZ2());
 
    //Double_t sigma_lambda = 1./(1.+track.GetTgl()*track.GetTgl()) * TMath::Sqrt(track.GetSigmaTgl2()); 
    //Double_t sigma_phi = 1./TMath::Sqrt(1-track.GetSnp()*track.GetSnp()) * TMath::Sqrt(track.GetSigmaSnp2());
    pullPhiTPC = (track.GetSnp() - mcsnplocal) / TMath::Sqrt(track.GetSigmaSnp2());
    pullLambdaTPC = (track.GetTgl() - mctgl) / TMath::Sqrt(track.GetSigmaTgl2());

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track.GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track.GetSigmaSnp2()); 
    if (mcpt) pull1PtTPC = (track.OneOverPt()-1./mcpt) / TMath::Sqrt(track.GetSigma1Pt2());
    else pull1PtTPC = 0.; 

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,particle->Vy(),particle->Vz(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);
  }
}

//_____________________________________________________________________________
void AliPerformanceRes::ProcessTPCITS(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent)
{
  // Fill resolution comparison information (TPC+ITS)
  if(!vTrack) return;
  if(!vEvent) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }
  
  //
  //Get DCA information
  //
  Float_t dca[2] = {0.,0.}; // dca_xy, dca_z
  Float_t cov[3] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParameters(dca,cov);
  
  //
  // Fill rec vs MC information
  //
  if(!mcev) return;

  Int_t label = TMath::Abs(vTrack->GetLabel()); //Use global label for combined resolution analysis
  TParticle* particle = ((AliMCParticle*)mcev->GetTrack(label))->Particle();
  if(!particle) return;
  if(!particle->GetPDG()) return;
  if(particle->GetPDG()->Charge()==0) return;
  //printf("charge %d \n",particle->GetPDG()->Charge());


  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) return;

  // exclude electrons
  if (fCutsMC.GetEM()==TMath::Abs(particle->GetPdgCode())) return;

  Float_t mceta =  particle->Eta();
  Float_t mcphi =  particle->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = particle->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(particle->Py(),particle->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(particle->Pz(),particle->Pt()));

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return; // min. nb. TPC clusters
  if(vTrack->GetITSclusters(0)<fCutsRC.GetMinNClustersITS()) return;  // min. nb. ITS clusters
  
  Float_t deltaPtTPC=0, deltaYTPC=0, deltaZTPC=0, deltaPhiTPC=0, deltaLambdaTPC=0; 
  Float_t pull1PtTPC=0, pullYTPC=0, pullZTPC=0, pullPhiTPC=0, pullLambdaTPC=0; 

  // select primaries
  if(TMath::Abs(dca[0])<fCutsRC.GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC.GetMaxDCAToVertexZ()) 
  { 
    if(mcpt == 0) return;
    Double_t trackYvert = vTrack->Yv(); // ****  should initialize to zero
    Double_t trackZvert = vTrack->Zv();
    Double_t trackP[3] = {0.,0.,0.};
    vTrack->GetPxPyPz(trackP);
    Double_t trackPt = sqrt(trackP[0]*trackP[0] + trackP[1]*trackP[1]);
    deltaYTPC= trackYvert-particle->Vy();
    deltaZTPC = trackZvert-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(trackP[2],trackPt)-TMath::ATan2(particle->Pz(),particle->Pt());
    deltaPhiTPC = TMath::ATan2(trackP[1],trackP[0])-TMath::ATan2(particle->Py(),particle->Px());
    deltaPtTPC = (trackPt-mcpt) / mcpt;

    AliExternalTrackParam extTrackParam;
    extTrackParam.CopyFromVTrack(vTrack);
    pullYTPC= (trackYvert-particle->Vy()) / TMath::Sqrt(extTrackParam.GetSigmaY2());
    pullZTPC = (trackZvert-particle->Vz()) / TMath::Sqrt(extTrackParam.GetSigmaZ2());
    pullPhiTPC = (extTrackParam.GetSnp() - mcsnp) / TMath::Sqrt(extTrackParam.GetSigmaSnp2());
    pullLambdaTPC = (extTrackParam.GetTgl() - mctgl) / TMath::Sqrt(extTrackParam.GetSigmaTgl2());
    if (mcpt) pull1PtTPC = (extTrackParam.OneOverPt()-1./mcpt) / TMath::Sqrt(extTrackParam.GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,particle->Vy(),particle->Vz(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);

   
    /*
    deltaYTPC= esdTrack->GetY()-particle->Vy();
    deltaZTPC = esdTrack->GetZ()-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(esdTrack->Pz(),esdTrack->Pt())-TMath::ATan2(particle->Pz(),particle->Pt());
    deltaPhiTPC = TMath::ATan2(esdTrack->Py(),esdTrack->Px())-TMath::ATan2(particle->Py(),particle->Px());
    delta1PtTPC = (esdTrack->OneOverPt()-1./mcpt)*mcpt;

    pullYTPC= (esdTrack->GetY()-particle->Vy()) / TMath::Sqrt(esdTrack->GetSigmaY2());
    pullZTPC = (esdTrack->GetZ()-particle->Vz()) / TMath::Sqrt(esdTrack->GetSigmaZ2());

    Double_t sigma_lambda = 1./(1.+esdTrack->GetTgl()*esdTrack->GetTgl()) * TMath::Sqrt(esdTrack->GetSigmaTgl2()); 
    Double_t sigma_phi = 1./TMath::Sqrt(1-esdTrack->GetSnp()*esdTrack->GetSnp()) * TMath::Sqrt(esdTrack->GetSigmaSnp2());
    pullLambdaTPC = deltaLambdaTPC / sigma_lambda;
    pullPhiTPC = deltaPhiTPC / sigma_phi;

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(esdTrack->GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(esdTrack->GetSigmaSnp2()); 
    if (mcpt) pull1PtTPC = (esdTrack->OneOverPt()-1./mcpt) / TMath::Sqrt(esdTrack->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,delta1PtTPC,particle->Vy(),particle->Vz(),mceta,mcphi,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mceta,mcphi,mcpt};
    fPullHisto->Fill(vPullHisto);
    */
  }
}

//_____________________________________________________________________________
void AliPerformanceRes::ProcessConstrained(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent)
{
  // Fill resolution comparison information (constrained parameters) 
  //
  if(!vEvent) return;
  if(!vTrack) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }


  const AliExternalTrackParam * track = vTrack->GetConstrainedParam();
  if(!track) return;

  Float_t dca[2] = {0.,0.}; // dca_xy, dca_z
  Float_t cov[3] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParameters(dca,cov);
 
  //
  // Fill rec vs MC information
  //
  if(!mcev) return;

  Int_t label = TMath::Abs(vTrack->GetLabel()); 
  TParticle* particle = ((AliMCParticle*)mcev->GetTrack(label))->Particle();
  if(!particle) return;
  if(!particle->GetPDG()) return;
  if(particle->GetPDG()->Charge()==0) return;
  //printf("charge %d \n",particle->GetPDG()->Charge());

  // Only 5 charged particle species (e,mu,pi,K,p)
  if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) return;

  // exclude electrons
  if (fCutsMC.GetEM()==TMath::Abs(particle->GetPdgCode())) return;

  Float_t mceta =  particle->Eta();
  Float_t mcphi =  particle->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = particle->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(particle->Py(),particle->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(particle->Pz(),particle->Pt()));

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return; // min. nb. TPC clusters

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  // select primaries
  if(TMath::Abs(dca[0])<fCutsRC.GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC.GetMaxDCAToVertexZ()) 
  { 

    if(mcpt == 0) return;
    
    deltaYTPC= track->GetY()-particle->Vy();
    deltaZTPC = track->GetZ()-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(particle->Pz(),particle->Pt());
    deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(particle->Py(),particle->Px());
    //delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;
    deltaPtTPC = (track->Pt()-mcpt) / mcpt;

    pullYTPC= (track->GetY()-particle->Vy()) / TMath::Sqrt(track->GetSigmaY2());
    pullZTPC = (track->GetZ()-particle->Vz()) / TMath::Sqrt(track->GetSigmaZ2());
 
    //Double_t sigma_lambda = 1./(1.+track->GetTgl()*track->GetTgl()) * TMath::Sqrt(track->GetSigmaTgl2()); 
    //Double_t sigma_phi = 1./TMath::Sqrt(1-track->GetSnp()*track->GetSnp()) * TMath::Sqrt(track->GetSigmaSnp2());
    pullPhiTPC = (track->GetSnp() - mcsnp) / TMath::Sqrt(track->GetSigmaSnp2());
    pullLambdaTPC = (track->GetTgl() - mctgl) / TMath::Sqrt(track->GetSigmaTgl2());

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 
    if (mcpt) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,particle->Vy(),particle->Vz(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);

    /*

    deltaYTPC= track->GetY()-particle->Vy();
    deltaZTPC = track->GetZ()-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(particle->Pz(),particle->Pt());
    deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(particle->Py(),particle->Px());
    delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;

    pullYTPC= (track->GetY()-particle->Vy()) / TMath::Sqrt(track->GetSigmaY2());
    pullZTPC = (track->GetZ()-particle->Vz()) / TMath::Sqrt(track->GetSigmaZ2());

    Double_t sigma_lambda = 1./(1.+track->GetTgl()*track->GetTgl()) * TMath::Sqrt(track->GetSigmaTgl2()); 
    Double_t sigma_phi = 1./TMath::Sqrt(1-track->GetSnp()*track->GetSnp()) * TMath::Sqrt(track->GetSigmaSnp2());
    pullLambdaTPC = deltaLambdaTPC / sigma_lambda;
    pullPhiTPC = deltaPhiTPC / sigma_phi;

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 

    if (mcpt) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,delta1PtTPC,particle->Vy(),particle->Vz(),mceta,mcphi,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mceta,mcphi,mcpt};
    fPullHisto->Fill(vPullHisto);

    */
  }
}
 
//_____________________________________________________________________________
void AliPerformanceRes::ProcessInnerTPC(AliMCEvent *const mcEvent, AliVTrack *const vTrack, AliVEvent* const vEvent)
{
  //
  // Fill resolution comparison information (inner params at TPC reference point) 
  //
  if(!vEvent) return;
  if(!vTrack) return;
  
  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexTPCBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }

  const AliExternalTrackParam * innerParam = vTrack->GetInnerParam();
  if(!innerParam) return;

  // create new AliExternalTrackParam
  AliExternalTrackParam *track = new AliExternalTrackParam(*innerParam);
  if(!track) return;

  Float_t dca[2] = {0.,0.}; // dca_xy, dca_z
  Float_t cov[3] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParametersTPC(dca,cov);
 
  //
  // Fill rec vs MC information
  //
  if(!mcEvent) return;

  Int_t label = vTrack->GetTPCLabel(); //Use TPC-only label for TPC-only resolution analysis
  if (label <= 0) return;
  AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(label);
  if(!mcParticle) return;

  // get the first TPC track reference
  AliTrackReference *ref0 = GetFirstTPCTrackRef(mcParticle);
  if(!ref0) return;

  // Only 5 charged particle species (e,mu,pi,K,p)
  TParticle *particle = mcParticle->Particle();
  if(!particle) return;

  if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) return;

  // exclude electrons
  if (fCutsMC.GetEM()==TMath::Abs(particle->GetPdgCode())) return;

  Double_t mclocal[4]; //Rotated x,y,Px,Py mc-coordinates - the MC data should be rotated since the track is propagated best along x
  Double_t c = TMath::Cos(track->GetAlpha());
  Double_t s = TMath::Sin(track->GetAlpha());
  Double_t x = ref0->X();
  Double_t y = ref0->Y();
  mclocal[0] = x*c + y*s;
  mclocal[1] =-x*s + y*c;
  Double_t px = ref0->Px();
  Double_t py = ref0->Py();
  mclocal[2] = px*c + py*s;
  mclocal[3] =-px*s + py*c;

  //Double_t xyz[3] = {ref0->X(),ref0->Y(),ref0->Z()};
  // propagate track to the radius of the first track reference within TPC
  //Double_t trRadius = TMath::Sqrt(xyz[1] * xyz[1] + xyz[0] * xyz[0]);
  Double_t field[3]; track->GetBxByBz(field);
  if (TGeoGlobalMagField::Instance()->GetField() == NULL) {
    Error("ProcessInnerTPC", "Magnetic Field not set");
  }
  Bool_t isOK = track->PropagateToBxByBz(mclocal[0],field);
  if(!isOK) {return;}
  //track->AliExternalTrackParam::PropagateTo(mclocal[0],esdEvent->GetMagneticField());  //Use different propagation since above methods returns zero B field and does not work
  
  Float_t mceta =  -TMath::Log(TMath::Tan(0.5 * ref0->Theta()));
  Float_t mcphi =  ref0->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = ref0->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(ref0->Py(),ref0->Px()));
  Float_t mcsnplocal = TMath::Sin(TMath::ATan2(mclocal[3],mclocal[2]));
  Float_t mctgl = TMath::Tan(TMath::ATan2(ref0->Pz(),ref0->Pt()));

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return; // min. nb. TPC clusters

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  // select primaries
  Bool_t isPrimary;
  if ( IsUseTrackVertex() )
  {
    isPrimary = TMath::Abs(dca[0])<fCutsRC.GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC.GetMaxDCAToVertexZ();
  }
  else
  {
	  //If Track vertex is not used the above check does not work, hence we use the MC reference track
    isPrimary = mcEvent->IsPhysicalPrimary(label);
  }
  if(isPrimary) 
  { 
    if(mcpt == 0) return;
    
    deltaYTPC= track->GetY()-mclocal[1];
    deltaZTPC = track->GetZ()-ref0->Z();
    deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(ref0->Pz(),ref0->Pt());
	//track->Px() etc returns momentum in global coordinates, hence mc momentum must not be rotated.
    deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(ref0->Py(),ref0->Px());
    //delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;
    deltaPtTPC = (track->Pt()-mcpt) / mcpt;

    pullYTPC= deltaYTPC / TMath::Sqrt(track->GetSigmaY2());
    pullZTPC = deltaZTPC / TMath::Sqrt(track->GetSigmaZ2());
 
    //Double_t sigma_lambda = 1./(1.+track->GetTgl()*track->GetTgl()) * TMath::Sqrt(track->GetSigmaTgl2()); 
    //Double_t sigma_phi = 1./TMath::Sqrt(1-track->GetSnp()*track->GetSnp()) * TMath::Sqrt(track->GetSigmaSnp2());
	//track->GetSnp returns value in local coordinates, hence here, in contrast to deltaPhiTPC, the mc momentum must be rotated
    pullPhiTPC = (track->GetSnp() - mcsnplocal) / TMath::Sqrt(track->GetSigmaSnp2());
    pullLambdaTPC = (track->GetTgl() - mctgl) / TMath::Sqrt(track->GetSigmaTgl2());

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 
    if (mcpt) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,ref0->Y(),ref0->Z(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,ref0->Y(),ref0->Z(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);
  }

  if(track) delete track;
}

//_____________________________________________________________________________
void AliPerformanceRes::ProcessOuterTPC(AliMCEvent *const mcEvent, AliVTrack *const vTrack, const AliVfriendTrack *const friendTrack, AliVEvent* const vEvent)
{
  //
  // Fill resolution comparison information (outer params at TPC reference point) 
  //
  if(!friendTrack) return;
  if(!vEvent) return;
  if(!vTrack) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexTPCBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }

  const AliExternalTrackParam * outerParam = friendTrack->GetTPCOut();
  if(!outerParam) return;

  // create new AliExternalTrackParam
  AliExternalTrackParam *track = new AliExternalTrackParam(*outerParam);
  if(!track) return;

  Float_t dca[2] = {0.,0.}; // dca_xy, dca_z
  Float_t cov[3] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParametersTPC(dca,cov);
 
  //
  // Fill rec vs MC information
  //
  if(!mcEvent) return;

  Int_t label = vTrack->GetTPCLabel(); //Use TPC-only label for TPC-only resolution analysis
  if (label <= 0) return;
  AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(label);
  if(!mcParticle) return;

  // get the last TPC track reference
  AliTrackReference *ref0 = GetLastTPCTrackRef(mcParticle);
  if(!ref0) return;

  // Only 5 charged particle species (e,mu,pi,K,p)
  TParticle *particle = mcParticle->Particle();
  if(!particle) return;

  if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) return;

  // exclude electrons
  if (fCutsMC.GetEM()==TMath::Abs(particle->GetPdgCode())) return;

  // calculate alpha angle
  Double_t xyz[3] = {ref0->X(),ref0->Y(),ref0->Z()};
  Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
  //if (alpha<0) alpha += TMath::TwoPi();

  // rotate outer track to local coordinate system
  // and propagate to the radius of the last track reference of TPC
  Double_t trRadius = TMath::Sqrt(xyz[1] * xyz[1] + xyz[0] * xyz[0]);
  //Bool_t isOK = track->Propagate(alpha,trRadius,AliTracker::GetBz());
  Double_t field[3]; track->GetBxByBz(field);
  Bool_t isOK = track->PropagateBxByBz(alpha,trRadius,field);
  if(!isOK) return;

  Float_t mceta =  -TMath::Log(TMath::Tan(0.5 * ref0->Theta()));
  Float_t mcphi =  ref0->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = ref0->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(ref0->Py(),ref0->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(ref0->Pz(),ref0->Pt()));

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return; // min. nb. TPC clusters

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  // select primaries
  if(TMath::Abs(dca[0])<fCutsRC.GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC.GetMaxDCAToVertexZ()) 
  { 
    if(mcpt == 0) return;
    
    deltaYTPC= track->GetY(); // already rotated
    deltaZTPC = track->GetZ()-ref0->Z();
    deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(ref0->Pz(),ref0->Pt());
    deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(ref0->Py(),ref0->Px());
    //delta1PtTPC = (track->OneOverPt()-1./mcpt)*mcpt;
    deltaPtTPC = (track->Pt()-mcpt) / mcpt;

    pullYTPC= track->GetY() / TMath::Sqrt(track->GetSigmaY2());
    pullZTPC = (track->GetZ()-ref0->Z()) / TMath::Sqrt(track->GetSigmaZ2());
 
    //Double_t sigma_lambda = 1./(1.+track->GetTgl()*track->GetTgl()) * TMath::Sqrt(track->GetSigmaTgl2()); 
    //Double_t sigma_phi = 1./TMath::Sqrt(1-track->GetSnp()*track->GetSnp()) * TMath::Sqrt(track->GetSigmaSnp2());
    pullPhiTPC = (track->GetSnp() - mcsnp) / TMath::Sqrt(track->GetSigmaSnp2());
    pullLambdaTPC = (track->GetTgl() - mctgl) / TMath::Sqrt(track->GetSigmaTgl2());

    //pullLambdaTPC = deltaLambdaTPC / TMath::Sqrt(track->GetSigmaTgl2());
    //pullPhiTPC = deltaPhiTPC / TMath::Sqrt(track->GetSigmaSnp2()); 
    if (mcpt) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,ref0->Y(),ref0->Z(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,ref0->Y(),ref0->Z(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);
  }

  if(track) delete track;
}

//_____________________________________________________________________________
AliTrackReference * AliPerformanceRes::GetFirstTPCTrackRef(AliMCParticle *mcParticle) 
{
  // return first TPC track reference 
  if(!mcParticle) return 0;

  // find first track reference 
  // check direction to select proper reference point for looping tracks
  Int_t nTrackRef = mcParticle->GetNumberOfTrackReferences();
  AliTrackReference *ref = 0;
  AliTrackReference *refIn = 0;
  for (Int_t iref = 0; iref < nTrackRef; iref++) { 
    ref = mcParticle->GetTrackReference(iref);
    if(ref && (ref->DetectorId()==AliTrackReference::kTPC))
    {
      Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
      if(dir < 0.) break;

      refIn = ref;
      break;
    }
  }

return refIn;
}

//_____________________________________________________________________________
AliTrackReference * AliPerformanceRes::GetLastTPCTrackRef(AliMCParticle *mcParticle) 
{
  // return last TPC track reference 
  if(!mcParticle) return 0;

  // find last track reference 
  // check direction to select proper reference point for looping tracks
  Int_t nTrackRef = mcParticle->GetNumberOfTrackReferences();
  AliTrackReference *ref = 0;
  AliTrackReference *refOut = 0;
  for (Int_t iref = 0; iref < nTrackRef; iref++) { 
    ref = mcParticle->GetTrackReference(iref);
    if(ref && (ref->DetectorId()==AliTrackReference::kTPC)) { 
      Float_t dir=ref->X()*ref->Px()+ref->Y()*ref->Py();
      if(dir< 0.0) break;
      refOut = ref;
    }
  }

return refOut;
}

//_____________________________________________________________________________
void AliPerformanceRes::Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vfriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend)
{
  // Process comparison information
  if(!vEvent) 
  {
    Error("Exec","vEvent not available");
    return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  TArrayF vtxMC(3);
  
  if(bUseMC)
  {
    if(!mcEvent) {
      Error("Exec","mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      Error("Exec","Header not available");
      return;
    }
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      Error("Exec","Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);
  } 
  else {
    Error("Exec","MC information required!");
    return;
  }
  
  // use V friends
  if(bUseVfriend) {
    if(!vfriendEvent) {
      Error("Exec","vfriendEvent not available");
      return;
    }
  }

  // get event vertex
  const AliVVertex *vVertex = NULL;
  if( IsUseTrackVertex() ) 
  { 
    // track vertex
    vVertex = vEvent->GetPrimaryVertexTracks();
    if(vVertex && (vVertex->GetStatus()<=0)) return;
  }
  // Coverity - removed else branch as vtxESD is not further used in method
  //  else {  
  //    // TPC track vertex
  //    vtxESD = esdEvent->GetPrimaryVertexTPC();
  //  }
 



  //  Process events
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  { 
    AliVTrack *vTrack = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!vTrack) continue;
    
    const AliVfriendTrack *friendTrack= NULL;
    
    Int_t label = TMath::Abs(vTrack->GetLabel());
    /* // RS this check is not needed
    if ( label > mcEvent->GetNumberOfTracks() ) 
    {
      ULong_t status = track->GetStatus();
      printf ("Error : ESD MCLabel %d - StackSize %d - Status %lu \n",
	       track->GetLabel(), stack->GetNtrack(), status );
      printf(" NCluster %d \n", track->GetTPCclusters(0) );
//      if ((status&AliESDtrack::kTPCrefit)== 0 ) printf("   kTPCrefit \n");
//      if ((status&AliESDtrack::kTPCin)== 0 )    printf("   kTPCin \n");
//      if ((status&AliESDtrack::kTPCout)== 0 )   printf("   kTPCout \n");
//      if ((status&AliESDtrack::kTRDrefit)== 0 ) printf("   kTRDrefit \n");
//      if ((status&AliESDtrack::kTRDin)== 0 )    printf("   kTRDin \n");
//      if ((status&AliESDtrack::kTRDout)== 0 )   printf("   kTRDout \n");
//      if ((status&AliESDtrack::kITSrefit)== 0 ) printf("   kITSrefit \n");
//      if ((status&AliESDtrack::kITSin)== 0 )    printf("   kITSin \n");
//      if ((status&AliESDtrack::kITSout)== 0 )   printf("   kITSout \n");

      continue;
    }
*/
    if (label == 0) continue;		//Cannot distinguish between track or fake track
    if (vTrack->GetLabel() < 0) continue; //Do not consider fake tracks

    if(GetAnalysisMode() == 0) ProcessTPC(mcEvent,vTrack,vEvent);
    else if(GetAnalysisMode() == 1) ProcessTPCITS(mcEvent,vTrack,vEvent);
    else if(GetAnalysisMode() == 2) ProcessConstrained(mcEvent,vTrack,vEvent);
    else if(GetAnalysisMode() == 3) ProcessInnerTPC(mcEvent,vTrack,vEvent);
    else if(GetAnalysisMode() == 4) {
      if(bUseVfriend && vfriendEvent && vfriendEvent->TestSkipBit()==kFALSE && iTrack<vfriendEvent->GetNumberOfTracks()) {
	friendTrack=vfriendEvent->GetTrack(iTrack);
	if(!friendTrack) continue;
      }
      ProcessOuterTPC(mcEvent,vTrack,friendTrack,vEvent);
    }
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }
}

//_____________________________________________________________________________
TH1F* AliPerformanceRes::MakeResol(TH2F * his, Int_t integ, Bool_t type, Int_t cut){
  // Create resolution histograms

  //Debugging: Turning off creation/display of histograms because they
  //seem to get stuck in a loop
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  bool shouldDrawBinFits = kFALSE; //default = kTRUE
  hisr = AliTreeDraw::CreateResHistoII(his,&hism,integ,shouldDrawBinFits,cut);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliPerformanceRes::Analyse() {
  // Analyse comparison information and store output histograms
  // in the folder "folderRes"
  //
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH2F *h2D=0;
  TObjArray *aFolderObj = new TObjArray;
  if(!aFolderObj) return;

  // write results in the folder 
  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  c->cd();

  char name[256];
  char title[256];

  for(Int_t i=0; i<5; i++) 
  {
    for(Int_t j=5; j<10; j++) 
    {
      //Ranges
      if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
      else fResolHisto->GetAxis(8)->SetRangeUser(-1.5,1.49);
      if (j!=9) fResolHisto->GetAxis(9)->SetRangeUser(0.1,19.99);            // pt threshold
      else fResolHisto->GetAxis(9)->SetRangeUser(0.015,19.99);
      if(GetAnalysisMode() == 3) fResolHisto->GetAxis(5)->SetRangeUser(-80.,79.99); // y range

      if(j!=8) fPullHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
      else  fPullHisto->GetAxis(8)->SetRangeUser(-1.5,1.49);      // eta window
      if (j!=9) fPullHisto->GetAxis(9)->SetRangeUser(0.1,19.99);            // pt threshold
      else fPullHisto->GetAxis(9)->SetRangeUser(0.015,19.99);
      if(GetAnalysisMode() == 3) fPullHisto->GetAxis(5)->SetRangeUser(-80.,79.99); // y range
      
      //Resolutions

      h2D = (TH2F*)fResolHisto->Projection(i,j);

      h = AliPerformanceRes::MakeResol(h2D,1,0,100);
      snprintf(name,256,"h_res_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fResolHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fResolHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fResolHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      h = AliPerformanceRes::MakeResol(h2D,1,1,100);
      //h = (TH1F*)arr->At(1);
      snprintf(name,256,"h_mean_res_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fResolHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fResolHisto->GetAxis(i)->GetTitle(),"(mean)");
      h->GetYaxis()->SetTitle(title);

      snprintf(title,256,"%s vs %s",title,fResolHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      fResolHisto->GetAxis(8)->SetRangeUser(-1.5,1.5);
      //fResolHisto->GetAxis(9)->SetRangeUser(0.1,100.);

      //Pulls

      h2D = (TH2F*)fPullHisto->Projection(i,j);

      h = AliPerformanceRes::MakeResol(h2D,1,0,100);
      snprintf(name,256,"h_pull_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      //if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      h = AliPerformanceRes::MakeResol(h2D,1,1,100);
      snprintf(name,256,"h_mean_pull_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(mean)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      //if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);
    }
  }

  for (Int_t i = 0;i < fResolHisto->GetNdimensions();i++)
  {
    fResolHisto->GetAxis(i)->UnZoom();				//Reset Range
  }
  for (Int_t i = 0;i < fPullHisto->GetNdimensions();i++)
  {
    fPullHisto->GetAxis(i)->UnZoom();				//Reset Range
  }

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliPerformanceRes::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceRes * comp=this;
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
Long64_t AliPerformanceRes::Merge(TCollection* const list) 
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
  AliPerformanceRes* entry = dynamic_cast<AliPerformanceRes*>(obj);
  if (entry == 0) continue; 
  if (fResolHisto->GetEntries()<fgkMergeEntriesCut){
    fResolHisto->Add(entry->fResolHisto);  
    fPullHisto->Add(entry->fPullHisto);
  }

  count++;
  }

return count;
}

//_____________________________________________________________________________
TFolder* AliPerformanceRes::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
