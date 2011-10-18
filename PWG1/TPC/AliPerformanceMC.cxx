//------------------------------------------------------------------------------
// Implementation of AliPerformanceMC class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformanceMC.
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceMC * compObj = (AliPerformanceMC*)coutput->FindObject("AliPerformanceMC");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderRes" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_MC.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"

#include "AliPerformanceMC.h" 
#include "AliESDEvent.h" 
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliESDfriend.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliStack.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTreeDraw.h" 
#include "AliTPCseed.h" 

using namespace std;

ClassImp(AliPerformanceMC)

//_____________________________________________________________________________
AliPerformanceMC::AliPerformanceMC():
  AliPerformanceObject("AliPerformanceMC"),
  fResolHisto(0),
  fPullHisto(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  Init();
}

//_____________________________________________________________________________
AliPerformanceMC::AliPerformanceMC(Char_t* name="AliPerformanceMC", Char_t* title="AliPerformanceMC",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),
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
AliPerformanceMC::~AliPerformanceMC()
{
  // destructor
   
  if(fResolHisto) delete fResolHisto; fResolHisto=0;     
  if(fPullHisto)  delete fPullHisto;  fPullHisto=0;     
  
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceMC::Init(){

  //
  // histogram bining
  //

  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 20.;

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
  Int_t binsResolHisto[10]={100,100,100,100,100,25,50,90,30,nPtBins};
  Double_t minResolHisto[10]={-0.2,-0.2,-0.002,-0.002,-0.002, yMin, zMin, 0., -1.5, ptMin};
  Double_t maxResolHisto[10]={ 0.2, 0.2, 0.002, 0.002, 0.002, yMax, zMax, 2.*TMath::Pi(), 1.5, ptMax};

  fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_pt:y:z:phi:eta:pt",10,binsResolHisto,minResolHisto,maxResolHisto);
  fResolHisto->SetBinEdges(9,binsPt);

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

  //pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt
  Int_t binsPullHisto[10]={100,100,100,100,100,50,50,50,50,nPtBins};
  Double_t minPullHisto[10]={-5.,-5.,-5.,-5.,-5., yMin, zMin,-1., -2.0, ptMin};
  Double_t maxPullHisto[10]={ 5., 5., 5., 5., 5., yMax, zMax, 1., 2.0, ptMax};
  fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt",10,binsPullHisto,minPullHisto,maxPullHisto);

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

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderMC","Analysis MC Folder");
}

//_____________________________________________________________________________
void AliPerformanceMC::Exec(AliMCEvent* const mcEvent, AliESDEvent *const /*esdEvent*/, AliESDfriend *const /*esdFriend*/, const Bool_t bUseMC, const Bool_t /*bUseESDfriend*/)
{
  // Process pure MC information 
  //
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
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
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      Error("Exec","Stack not available");
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
  
  Int_t nPart = mcEvent->GetNumberOfTracks();
  if (nPart==0) return;

  //TParticle * part = new TParticle;
  //TClonesArray * trefs = new TClonesArray("AliTrackReference");
  TParticle *part=0;
  TClonesArray *trefs=0;

  //  Process events
  for (Int_t iPart = 0; iPart < nPart; iPart++) 
  { 
    Int_t status = mcEvent->GetParticleAndTR(iPart, part, trefs);
    if (status<0) continue;
    if(!part) continue;
    if(!trefs) continue;

    AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(iPart);
    if(!mcParticle) continue;

    TParticle *particle = mcParticle->Particle();
    if(!particle) continue;

    Bool_t prim = stack->IsPhysicalPrimary(iPart);

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC->IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) break;

    // exclude electrons
    if (fCutsMC->GetEP()==TMath::Abs(particle->GetPdgCode())) return;

    AliTrackReference *refTPCIn = 0;
    AliTrackReference *refTPCOut = 0;
    AliTrackReference *refITSIn = 0;
    AliTrackReference *refITSOut = 0;
    AliTrackReference *refTRDIn = 0;
    AliTrackReference *refTRDOut = 0;

    Float_t rmax = 0.;
    Int_t nTrackRef = mcParticle->GetNumberOfTrackReferences();
    if(nTrackRef < 5) continue;

    for (Int_t iref = 0; iref < nTrackRef; iref++) 
    {
      AliTrackReference *ref = mcParticle->GetTrackReference(iref);
      if(!ref) continue;

      Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
      if(dir < 0.) break; // looping tracks (return back)
      if (ref->R()<rmax) break;

      // pt threshold
      if(ref->Pt() < fCutsRC->GetPtMin()) break;

      // TPC
      if(ref->DetectorId()==AliTrackReference::kTPC)
      {
	if(!refTPCIn) {
	  refTPCIn = ref;
	}
	else {
          if(ref->R() > refTPCIn->R())
	  refTPCOut = ref;
	}
      }
      // ITS
      if(ref->DetectorId()==AliTrackReference::kITS)
      {
	if(!refITSIn) {
	  refITSIn = ref;
	}
	else {
          if(ref->R() > refITSIn->R())
	  refITSOut = ref;
	}
      }
      // TRD
      if(ref->DetectorId()==AliTrackReference::kTRD)
      {
	if(!refTRDIn) {
	  refTRDIn = ref;
	}
	else {
          if(ref->R() > refTRDIn->R())
	  refTRDOut = ref;
	}
      }
      if (ref->R()>rmax) rmax=ref->R();
    }

    if(GetAnalysisMode() == 0) {
      // Propagate refTPCIn to DCA to prim. vertex and make comparison
      // only primary particles
      if(refTPCIn) { 
        if(prim) ProcessTPC(refTPCIn,particle); 
      }
    }

    if(GetAnalysisMode() == 3) {
      // Propagate refTPCOut to refTPCIn using AliTrack::PropagateTrackToBxByBz()
      //
      if(refTPCIn && refTPCOut) 
        ProcessInnerTPC(refTPCIn,refTPCOut,particle); 
    }

    if(GetAnalysisMode() == 4) { 
      // propagate refTPCIn to refTPCOut using AliExternalTrackParam::PropagateToBxByBz()
      ProcessOuterTPCExt(part,trefs); 
    }
  }

  //if(part) delete part;
  //if(trefs) delete trefs;
}

//_____________________________________________________________________________
void AliPerformanceMC::ProcessTPC(AliTrackReference* const refIn, TParticle *const particle) {
  //
  // Test propagation from refIn to DCA
  //
  if(!refIn) return;
  if(!particle) return;

  AliExternalTrackParam *track = 0;
  Double_t radius = 2.8; // cm
  Double_t mass = particle->GetMass();
  Double_t step=1;
  //
  track=MakeTrack(refIn,particle);
  if (!track) return;
  //
  //AliTracker::PropagateTrackTo(track, radius+step, mass, step, kTRUE,0.99);
  //AliTracker::PropagateTrackTo(track, radius+0.5, mass, step*0.1, kTRUE,0.99);
  AliTracker::PropagateTrackToBxByBz(track, radius+step, mass, step, kTRUE,0.99);
  AliTracker::PropagateTrackToBxByBz(track, radius+0.5, mass, step*0.1, kTRUE,0.99);
  Double_t xyz[3] = {particle->Vx(),particle->Vy(),particle->Vz()};
  Double_t sxyz[3]={0.0,0.0,0.0};
  AliESDVertex vertex(xyz,sxyz);
  Double_t dca[2], cov[3];
  //Bool_t isOK = track->PropagateToDCA(&vertex,AliTracker::GetBz(),10,dca,cov);
  Double_t field[3];  track->GetBxByBz(field); 
  Bool_t isOK = track->PropagateToDCABxByBz(&vertex,field,10,dca,cov);
  if(!isOK) return;

  // Fill histogram
  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC, pullYTPC, pullZTPC, pullPhiTPC, pullLambdaTPC; 

  Float_t mceta =  particle->Eta();
  Float_t mcphi =  particle->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = particle->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(particle->Py(),particle->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(particle->Pz(),particle->Pt()));

  //if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ()) 
  //{ 
    if(mcpt == 0) return;
    
    //deltaYTPC= track->GetY()-particle->Vy();
    deltaYTPC= track->GetY(); // it is already rotated
    deltaZTPC = track->GetZ()-particle->Vz();
    deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(particle->Pz(),particle->Pt());
    deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(particle->Py(),particle->Px());
    deltaPtTPC = (track->Pt()-mcpt) / mcpt;

    /*
    printf("------------------------------ \n");
    printf("trY %f, trZ %f, trTheta %f ,trPhi %f, trPt %f \n", track->GetY(), track->GetZ(), TMath::ATan2(track->Pz(),track->Pt()), TMath::ATan2(track->Py(),track->Px()), track->Pt() ); 
    printf("partY %f, partZ %f, partTheta %f ,partPhi %f, partPt %f \n", particle->Vy(), particle->Vz(), TMath::ATan2(particle->Pz(),particle->Pt()), TMath::ATan2(particle->Py(),particle->Px()), mcpt ); 
    */


    pullYTPC= (track->GetY()-particle->Vy()) / TMath::Sqrt(track->GetSigmaY2());
    pullZTPC = (track->GetZ()-particle->Vz()) / TMath::Sqrt(track->GetSigmaZ2());
 
    pullPhiTPC = (track->GetSnp() - mcsnp) / TMath::Sqrt(track->GetSigmaSnp2());
    pullLambdaTPC = (track->GetTgl() - mctgl) / TMath::Sqrt(track->GetSigmaTgl2());

    if (mcpt) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());
    else pull1PtTPC = 0.;

    Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,particle->Vy(),particle->Vz(),mcphi,mceta,mcpt};
    fResolHisto->Fill(vResolHisto);

    Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,particle->Vy(),particle->Vz(),mcsnp,mctgl,1./mcpt};
    fPullHisto->Fill(vPullHisto);
  //}
  if(track) delete track;
}

//_____________________________________________________________________________
void AliPerformanceMC::ProcessInnerTPC(AliTrackReference* const refIn,  AliTrackReference *const  refOut, TParticle *const particle) {
  //
  // Test propagation from Out to In
  //
  if(!refIn) return;
  if(!refOut) return;
  if(!particle) return;

  AliExternalTrackParam *track=MakeTrack(refOut,particle);
  if (!track) return;

  Double_t mass = particle->GetMass();
  Double_t step=1;

  Float_t alphaIn= TMath::ATan2(refIn->Y(),refIn->X());
  //Float_t alphaOut= TMath::ATan2(refOut->Y(),refOut->X());
  //track->Rotate(alphaIn-alphaOut);
  track->Rotate(alphaIn);
  //printf("alphaIn %f, alphaOut %f \n",alphaIn,alphaOut);

  //
  //Bool_t isOK = AliTracker::PropagateTrackTo(track, refIn->R(), mass, step, kTRUE,0.99);
  //Bool_t isOK = AliTracker::PropagateTrackTo(track, refIn->R(), mass, step, kFALSE,0.99);
  Bool_t isOK = AliTracker::PropagateTrackToBxByBz(track, refIn->R(), mass, step, kFALSE,0.99);
  if(!isOK) return;

  // calculate alpha angle
  //Double_t xyz[3] = {refIn->X(),refIn->Y(),refIn0->Z()};
  //Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
  //if (alpha<0) alpha += TMath::TwoPi();

  // rotate outer track to inner 
  // and propagate to the radius of the first track referenco of TPC
  //Double_t trRadius = TMath::Sqrt(xyz[1] * xyz[1] + xyz[0] * xyz[0]);
  //Bool_t isOK = track->Propagate(alpha,trRadius,AliTracker::GetBz());
  //if(!isOK) return;
  
  Float_t mceta =  -TMath::Log(TMath::Tan(0.5 * refIn->Theta()));
  Float_t mcphi =  refIn->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = refIn->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(refIn->Py(),refIn->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(refIn->Pz(),refIn->Pt()));

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC=0., pullYTPC=0., pullZTPC=0., pullPhiTPC=0., pullLambdaTPC=0.; 

  if(mcpt == 0) return;

  //deltaYTPC = track->GetY()-refIn->Y();
  deltaYTPC = track->GetY();
  deltaZTPC = track->GetZ()-refIn->Z();
  deltaLambdaTPC = TMath::ATan2(track->Pz(),track->Pt())-TMath::ATan2(refIn->Pz(),refIn->Pt());
  deltaPhiTPC = TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(refIn->Py(),refIn->Px());
  deltaPtTPC = (track->Pt()-mcpt) / mcpt;

  /*
  printf("------------------------------ \n");
  printf("trX %f, trY %f, trZ %f, trTheta %f ,trPhi %f, trPt %f \n", track->GetX(), track->GetY(), track->GetZ(), TMath::ATan2(track->Pz(),track->Pt()), TMath::ATan2(track->Py(),track->Px()), track->Pt() ); 
  printf("partX %f, partY %f, partZ %f, partTheta %f ,partPhi %f, partPt %f \n",refIn->X(),  refIn->Y(), refIn->Z(), TMath::ATan2(refIn->Pz(),refIn->Pt()), TMath::ATan2(refIn->Py(),refIn->Px()), mcpt ); 
  */

  if(TMath::Sqrt(track->GetSigmaY2())) pullYTPC = track->GetY() / TMath::Sqrt(track->GetSigmaY2());
  if(TMath::Sqrt(track->GetSigmaZ2())) pullZTPC = (track->GetZ()-refIn->Z()) / TMath::Sqrt(track->GetSigmaZ2());
  if(TMath::Sqrt(track->GetSigmaSnp2())) pullPhiTPC = (track->GetSnp() - mcsnp) / TMath::Sqrt(track->GetSigmaSnp2());
  if(TMath::Sqrt(track->GetSigmaTgl2())) pullLambdaTPC = (track->GetTgl() - mctgl) / TMath::Sqrt(track->GetSigmaTgl2());
  if(TMath::Sqrt(track->GetSigma1Pt2())) pull1PtTPC = (track->OneOverPt()-1./mcpt) / TMath::Sqrt(track->GetSigma1Pt2());

  //printf("pullYTPC %f,pullZTPC %f,pullPhiTPC %f,pullLambdaTPC %f,pull1PtTPC %f,refIn->Y() %f,refIn->Z() %f,mcsnp %f,mctgl %f,1./mcpt %f \n",pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,refIn->Y(),refIn->Z(),mcsnp,mctgl,1./mcpt);

  Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,refIn->Y(),refIn->Z(),mcphi,mceta,mcpt};
  fResolHisto->Fill(vResolHisto);

  Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,refIn->Y(),refIn->Z(),mcsnp,mctgl,1./mcpt};
  fPullHisto->Fill(vPullHisto);

  if(track) delete track;
}

//_____________________________________________________________________________
void AliPerformanceMC::ProcessOuterTPCExt(TParticle *const part, TClonesArray * const trefs)
{
  const Int_t kMinRefs=5;
  Int_t nrefs = trefs->GetEntries();
  if (nrefs<kMinRefs) return; // we should have enough references
  Int_t iref0 =-1;
  Int_t iref1 =-1;
  
  for (Int_t iref=0; iref<nrefs; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    if (!ref) continue;    
    Float_t dir = ref->X()*ref->Px()+ref->Y()*ref->Py();
    if (dir<0) break;
    if (ref->DetectorId()!=AliTrackReference::kTPC) continue;
    if (iref0<0) iref0 = iref;
    iref1 = iref;    
  }
  if (iref1-iref0<kMinRefs) return;
  Double_t covar[15];
  for (Int_t icov=0; icov<15; icov++) covar[icov]=0;
  covar[0]=1; 
  covar[2]=1; 
  covar[5]=1;
  covar[9]=1;
  covar[14]=1;

  AliTrackReference * refIn = (AliTrackReference*)trefs->At(iref0);
  AliTrackReference * refOut = (AliTrackReference*)trefs->At(iref1);
  AliExternalTrackParam *paramPropagate= MakeTrack(refIn,part);
  AliExternalTrackParam *paramUpdate   = MakeTrack(refIn,part);
  paramUpdate->AddCovariance(covar);

  Bool_t isOKP=kTRUE;
  Bool_t isOKU=kTRUE;
  AliMagF * field = (AliMagF*) TGeoGlobalMagField::Instance()->GetField();
  if(!field) return;
  for (Int_t iref = iref0; iref<=iref1; iref++){
    AliTrackReference * ref = (AliTrackReference*)trefs->At(iref);
    if(!ref) continue;
    Float_t alphaC= TMath::ATan2(ref->Y(),ref->X());
    Double_t pos[3] = {ref->X(), ref->Y(), ref->Z()};
    Double_t mag[3];
    field->Field(pos,mag);
    isOKP&=paramPropagate->Rotate(alphaC);
    isOKU&=paramUpdate->Rotate(alphaC);
    for (Float_t xref= paramPropagate->GetX(); xref<ref->R(); xref++){
      //isOKP&=paramPropagate->PropagateTo(xref, mag[2]);
      //isOKU&=paramUpdate->PropagateTo(xref, mag[2]);
      isOKP&=paramPropagate->PropagateToBxByBz(xref, mag);
      isOKU&=paramUpdate->PropagateToBxByBz(xref, mag);
    }
    //isOKP&=paramPropagate->PropagateTo(ref->R(), mag[2]);
    //isOKU&=paramUpdate->PropagateTo(ref->R(), mag[2]);
    isOKP&=paramPropagate->PropagateToBxByBz(ref->R(), mag);
    isOKU&=paramUpdate->PropagateToBxByBz(ref->R(), mag);
    Double_t clpos[2] = {0, ref->Z()};
    Double_t clcov[3] = { 0.005,0,0.005};
    isOKU&= paramUpdate->Update(clpos, clcov);  
  }

  Float_t mceta =  -TMath::Log(TMath::Tan(0.5 * refOut->Theta()));
  Float_t mcphi =  refOut->Phi();
  if(mcphi<0) mcphi += 2.*TMath::Pi();
  Float_t mcpt = refOut->Pt();
  Float_t mcsnp = TMath::Sin(TMath::ATan2(refOut->Py(),refOut->Px()));
  Float_t mctgl = TMath::Tan(TMath::ATan2(refOut->Pz(),refOut->Pt()));

  Float_t deltaPtTPC, deltaYTPC, deltaZTPC, deltaPhiTPC, deltaLambdaTPC; 
  Float_t pull1PtTPC=0., pullYTPC=0., pullZTPC=0., pullPhiTPC=0., pullLambdaTPC=0.; 

  if(mcpt == 0) return;

  //deltaYTPC = track->GetY()-refIn->Y();
  deltaYTPC = paramPropagate->GetY();
  deltaZTPC = paramPropagate->GetZ()-refOut->Z();
  deltaLambdaTPC = TMath::ATan2(paramPropagate->Pz(),paramPropagate->Pt())-TMath::ATan2(refOut->Pz(),refOut->Pt());
  deltaPhiTPC = TMath::ATan2(paramPropagate->Py(),paramPropagate->Px())-TMath::ATan2(refOut->Py(),refOut->Px());
  deltaPtTPC = (paramPropagate->Pt()-mcpt) / mcpt;

  if(TMath::Sqrt(paramPropagate->GetSigmaY2())) pullYTPC = paramPropagate->GetY() / TMath::Sqrt(paramPropagate->GetSigmaY2());
  if(TMath::Sqrt(paramPropagate->GetSigmaZ2())) pullZTPC = (paramPropagate->GetZ()-refOut->Z()) / TMath::Sqrt(paramPropagate->GetSigmaZ2());
  if(TMath::Sqrt(paramPropagate->GetSigmaSnp2())) pullPhiTPC = (paramPropagate->GetSnp() - mcsnp) / TMath::Sqrt(paramPropagate->GetSigmaSnp2());
  if(TMath::Sqrt(paramPropagate->GetSigmaTgl2())) pullLambdaTPC = (paramPropagate->GetTgl() - mctgl) / TMath::Sqrt(paramPropagate->GetSigmaTgl2());
  if(TMath::Sqrt(paramPropagate->GetSigma1Pt2())) pull1PtTPC = (paramPropagate->OneOverPt()-1./mcpt) / TMath::Sqrt(paramPropagate->GetSigma1Pt2());

  //printf("pullYTPC %f,pullZTPC %f,pullPhiTPC %f,pullLambdaTPC %f,pull1PtTPC %f,refIn->Y() %f,refIn->Z() %f,mcsnp %f,mctgl %f,1./mcpt %f \n",pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,refIn->Y(),refIn->Z(),mcsnp,mctgl,1./mcpt);

  Double_t vResolHisto[10] = {deltaYTPC,deltaZTPC,deltaPhiTPC,deltaLambdaTPC,deltaPtTPC,refIn->Y(),refIn->Z(),mcphi,mceta,mcpt};
  fResolHisto->Fill(vResolHisto);

  Double_t vPullHisto[10] = {pullYTPC,pullZTPC,pullPhiTPC,pullLambdaTPC,pull1PtTPC,refIn->Y(),refIn->Z(),mcsnp,mctgl,1./mcpt};
  fPullHisto->Fill(vPullHisto);
}
 
//_____________________________________________________________________________
AliExternalTrackParam * AliPerformanceMC::MakeTrack(const AliTrackReference* const ref, TParticle* const part)
{
  //
  // Make track out of the track ref
  // part - TParticle used to determine chargr
  // the covariance matrix - equal 0 - starting from ideal MC position
  Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
  Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
  Double_t cv[21];
  for (Int_t i=0; i<21;i++) cv[i]=0;
  if (!part->GetPDG()) return 0;
  AliExternalTrackParam * param = new AliExternalTrackParam(xyz,pxyz,cv,TMath::Nint(part->GetPDG()->Charge()/3.));
  return param;
}

//_____________________________________________________________________________
TH1F* AliPerformanceMC::MakeResol(TH2F * his, Int_t integ, Bool_t type, Int_t cut){
  // Create resolution histograms
 
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoII(his,&hism,integ,kTRUE,cut);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliPerformanceMC::Analyse() {
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
      //if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
      fResolHisto->GetAxis(9)->SetRangeUser(0.16,10.); // pt threshold
      if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(0.,0.9); // eta window
      //fResolHisto->GetAxis(9)->SetRangeUser(0.16,3.); // pt window
      if(GetAnalysisMode() == 3) fResolHisto->GetAxis(5)->SetRangeUser(-80.,80.); // y range

      h2D = (TH2F*)fResolHisto->Projection(i,j);

      h = AliPerformanceMC::MakeResol(h2D,1,0,20);
      snprintf(name,256,"h_res_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fResolHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fResolHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fResolHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      h = AliPerformanceMC::MakeResol(h2D,1,1,20);
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
      fResolHisto->GetAxis(9)->SetRangeUser(0.0,10.);

      //
      if(j!=8) fPullHisto->GetAxis(8)->SetRangeUser(-1.0,0.99); // theta window
      else  fPullHisto->GetAxis(8)->SetRangeUser(-1.5,1.5); // theta window
      fPullHisto->GetAxis(9)->SetRangeUser(0.0,10.);  // pt threshold
      if(GetAnalysisMode() == 3) fPullHisto->GetAxis(5)->SetRangeUser(-80.,80.); // y range

      h2D = (TH2F*)fPullHisto->Projection(i,j);

      h = AliPerformanceMC::MakeResol(h2D,1,0,20);
      snprintf(name,256,"h_pull_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      //if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      h = AliPerformanceMC::MakeResol(h2D,1,1,20);
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

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliPerformanceMC::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceMC * comp=this;
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
Long64_t AliPerformanceMC::Merge(TCollection* const list) 
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
  AliPerformanceMC* entry = dynamic_cast<AliPerformanceMC*>(obj);
  if (entry == 0) continue; 

  fResolHisto->Add(entry->fResolHisto);
  fPullHisto->Add(entry->fPullHisto);

  count++;
  }

return count;
}

//_____________________________________________________________________________
TFolder* AliPerformanceMC::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
