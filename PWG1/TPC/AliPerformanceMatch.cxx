//------------------------------------------------------------------------------
// Implementation of AliPerformanceMatch class. It checks the track matching 
// quality (residuals, pulls) bewteen tracking detectors TPC-ITS and TPC-TRD 
// with TPC as the reference detector. In addition, the ITS and TRD track 
// reconstruction efficiency with respect to TPC is calculated.
//
// The matchinig quality parameters are stored in the THnSparse histograms. 
// The analysis of these histograms to extract reference information is done in 
// the Analyse() function (post-processing).  
//  
// The histograms with reference information can be exported to the ROOT folder.
//
// Author: J.Otwinowski 17/10/2009 
// Changes by M.Knichel 22/10/2010 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceMatch * compObj = (AliPerformanceMatch*)coutput->FindObject("AliPerformanceMatchTPCTRD");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderMatch" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Match.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"

#include "AliPerformanceMatch.h" 
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
#include "AliTRDtrackV1.h" 
#include "AliTreeDraw.h" 

using namespace std;

ClassImp(AliPerformanceMatch)

Bool_t AliPerformanceMatch::fgMergeTHnSparse = kFALSE;

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch():
  AliPerformanceObject("AliPerformanceMatch"),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fFolderObj(0),
  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  
  fUseHLT(0)
{
  Init();
}

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch(Char_t* name="AliPerformanceMatch", Char_t* title="AliPerformanceMatch",Int_t analysisMode=0,Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fFolderObj(0),
  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  fUseHLT(0)
{
  // named constructor	
  // 
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);

  Init();
}

//_____________________________________________________________________________
AliPerformanceMatch::~AliPerformanceMatch()
{
  // destructor
   
  if(fResolHisto) delete fResolHisto; fResolHisto=0;     
  if(fPullHisto)  delete fPullHisto;  fPullHisto=0;
  if(fTrackingEffHisto) delete fTrackingEffHisto; fTrackingEffHisto = 0x0;
  
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
  if(fFolderObj) delete fFolderObj; fFolderObj=0;
}

//_____________________________________________________________________________
void AliPerformanceMatch::Init(){

  //
  // Make performance histogrms
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

  Double_t yMin = -0.02, yMax = 0.02;
  Double_t zMin = -12.0, zMax = 12.0;
  if(GetAnalysisMode() == 0 || GetAnalysisMode() == 1) { 
    yMin = -100.; yMax = 100.; 
    zMin = -100.; zMax = 100.; 
  }

  // res_y:res_z:res_phi,res_lambda:res_pt:y:z:eta:phi:pt:isRec
  Int_t binsResolHisto[11]={100,100,100,100,100,25,50,90,30,nPtBins,2};
  Double_t minResolHisto[11]={-1.,-1.,-0.03,-0.03,-0.2, yMin, zMin, 0., -1.5, ptMin,0};
  Double_t maxResolHisto[11]={ 1., 1., 0.03, 0.03, 0.2, yMax, zMax, 2.*TMath::Pi(), 1.5, ptMax,2};

  fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_pt:y:z:phi:eta:pt:isRec",11,binsResolHisto,minResolHisto,maxResolHisto);
  fResolHisto->SetBinEdges(9,binsPt);

  fResolHisto->GetAxis(0)->SetTitle("y-y_{ref} (cm)");
  fResolHisto->GetAxis(1)->SetTitle("z-z_{ref} (cm)");
  fResolHisto->GetAxis(2)->SetTitle("#phi-#phi_{ref} (rad)");
  fResolHisto->GetAxis(3)->SetTitle("#lambda-#lambda_{ref} (rad)");
  fResolHisto->GetAxis(4)->SetTitle("(p_{T}/p_{Tref}-1)");
  fResolHisto->GetAxis(5)->SetTitle("y_{ref} (cm)");
  fResolHisto->GetAxis(6)->SetTitle("z_{ref} (cm)");
  fResolHisto->GetAxis(7)->SetTitle("#phi_{ref} (rad)");
  fResolHisto->GetAxis(8)->SetTitle("#eta_{ref}");
  fResolHisto->GetAxis(9)->SetTitle("p_{Tref} (GeV/c)");
  fResolHisto->GetAxis(10)->SetTitle("isReconstructed");
  fResolHisto->Sumw2();

  //pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt:isRec
  Int_t binsPullHisto[11]={100,100,100,100,100,50,50,50,50,nPtBins,2};
  Double_t minPullHisto[11]={-5.,-5.,-5.,-5.,-5., yMin, zMin,-1.,-2.0, ptMin,0};
  Double_t maxPullHisto[11]={ 5., 5., 5., 5., 5., yMax, zMax, 1., 2.0, ptMax,2};
  fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt:isRec",11,binsPullHisto,minPullHisto,maxPullHisto);

  fPullHisto->GetAxis(0)->SetTitle("(y-y_{ref})/#sigma");
  fPullHisto->GetAxis(1)->SetTitle("(z-z_{ref})/#sigma");
  fPullHisto->GetAxis(2)->SetTitle("(sin#phi-sin#phi_{ref})/#sigma");
  fPullHisto->GetAxis(3)->SetTitle("(tan#lambda-tan#lambda_{ref})/#sigma");
  fPullHisto->GetAxis(4)->SetTitle("(p_{Tref}/p_{T}-1)/#sigma");
  fPullHisto->GetAxis(5)->SetTitle("y_{ref} (cm)");
  fPullHisto->GetAxis(6)->SetTitle("z_{ref} (cm)");
  fPullHisto->GetAxis(7)->SetTitle("sin#phi_{ref}");
  fPullHisto->GetAxis(8)->SetTitle("tan#lambda_{ref}");
  fPullHisto->GetAxis(9)->SetTitle("1/p_{Tref} (GeV/c)^{-1}");
  fPullHisto->GetAxis(10)->SetTitle("isReconstructed");
  fPullHisto->Sumw2();

  // -> has match:y:z:snp:tgl:phi:pt:ITSclusters
  Int_t binsTrackingEffHisto[8]    = { 2,    50, 100, 50, 50, 90,           100,  7   };
  Double_t minTrackingEffHisto[8]  = {-0.5, -25, -50, -1, -2, 0.,            0,   -0.5 };
  Double_t maxTrackingEffHisto[8]  = { 1.5,  25,  50,  1,  2, 2*TMath::Pi(), 20,   6.5 };
  
  fTrackingEffHisto = new THnSparseF("fTrackingEffHisto","has match:y:z:snp:tgl:phi:pt:ITSclusters",8,binsTrackingEffHisto,minTrackingEffHisto,maxTrackingEffHisto);
  fTrackingEffHisto->GetAxis(0)->SetTitle("IsMatching");
  fTrackingEffHisto->GetAxis(1)->SetTitle("local y (cm)");
  fTrackingEffHisto->GetAxis(2)->SetTitle("z (cm)");
  fTrackingEffHisto->GetAxis(3)->SetTitle("sin(#phi)");
  fTrackingEffHisto->GetAxis(4)->SetTitle("tan(#lambda)");
  fTrackingEffHisto->GetAxis(5)->SetTitle("phi (rad)");
  fTrackingEffHisto->GetAxis(6)->SetTitle("p_{T}");
  fTrackingEffHisto->GetAxis(7)->SetTitle("number of ITS clusters");
  fTrackingEffHisto->Sumw2();

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderMatch","Analysis Matching Folder");
}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessITSTPC(Int_t iTrack, AliESDEvent *const esdEvent, AliStack* /*const stack*/, AliESDtrack *const esdTrack, AliESDfriendTrack *const esdFriendTrack)
{
  //
  // addition to standard analysis - check if ITS stand-alone tracks have a match in the TPC
  // Origin: A. Kalwait
  // Modified: J. Otwinowski
  if(!esdEvent) return;
  if(!esdTrack) return;
  if(!esdFriendTrack) return;

  // ITS stand alone tracks with SPD layers 
  if(!(esdTrack->GetStatus() & AliESDtrack::kITSpureSA)) return;
  if(!(esdTrack->GetStatus() & AliESDtrack::kITSrefit)) return;
  if(esdTrack->GetNcls(0)<fCutsRC->GetMinNClustersITS()) return;
  if(!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1)) return;

  Bool_t hasMatch = kFALSE;
    for (Int_t jTrack = 0; jTrack < esdEvent->GetNumberOfTracks(); jTrack++) {
      // loop for all those tracks and check if the corresponding TPC track is found
      if (jTrack==iTrack) continue;
      AliESDtrack *trackTPC = esdEvent->GetTrack(jTrack);
      if (!trackTPC) continue;
      if (!trackTPC->GetTPCInnerParam()) continue;

      // TPC nClust/track after first tracking pass
      if(trackTPC->GetTPCNclsIter1()<fCutsRC->GetMinNClustersTPC()) continue; 

      AliExternalTrackParam *innerTPC = new AliExternalTrackParam(*(trackTPC->GetTPCInnerParam()));
      if(!innerTPC) continue;

      Double_t x[3]; trackTPC->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);
      Double_t dca[2],cov[3];
      Bool_t isTPCOK = innerTPC->PropagateToDCABxByBz(esdEvent->GetPrimaryVertexSPD(),b,kVeryBig,dca,cov);
      if(!isTPCOK) { 
        if(innerTPC) delete innerTPC;  
        continue;
      }

      //
      // select primaries
      //
      Double_t dcaToVertex = -1;
      if( fCutsRC->GetDCAToVertex2D() ) 
      {
        dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC->GetMaxDCAToVertexXY()/fCutsRC->GetMaxDCAToVertexXY()+dca[1]*dca[1]/fCutsRC->GetMaxDCAToVertexZ()/fCutsRC->GetMaxDCAToVertexZ()); 
      }
      if(fCutsRC->GetDCAToVertex2D() && dcaToVertex > 1) { 
        delete innerTPC; continue; }
      if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC->GetMaxDCAToVertexXY()) { 
        delete innerTPC; continue; }
      if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC->GetMaxDCAToVertexZ()) { 
        delete innerTPC; continue; }

      // rough checking if they match
      /*
      printf("+++++++++++++++++++++++++++++++++++++++++++++++ ");
      printf("esdTrack->GetY() %e, esdTrack->GetSnp() %e, esdTrack->GetTgl() %e \n", 
              esdTrack->GetY(), esdTrack->GetSnp(), esdTrack->GetTgl());
      printf("innerTPC->GetY() %e, innerTPC->GetSnp() %e, innerTPC->GetTgl() %e \n", 
              innerTPC->GetY() , innerTPC->GetSnp() , innerTPC->GetTgl());
      */
      if (TMath::Abs(esdTrack->GetY() - innerTPC->GetY()) > 3) { 
        delete innerTPC; continue; 
      }
      if (TMath::Abs(esdTrack->GetSnp() - innerTPC->GetSnp()) > 0.2) { 
        delete innerTPC; continue; 
      }
      if (TMath::Abs(esdTrack->GetTgl() - innerTPC->GetTgl()) > 0.2) { 
        delete innerTPC; continue; 
      }

      hasMatch = kTRUE;
      if(innerTPC) delete innerTPC;
    }
    //has match:y:z:snp:tgl:phi:pt:ITSclusters
    Double_t vecTrackingEff[8] = { hasMatch,esdTrack->GetY(),esdTrack->GetZ(),esdTrack->GetSnp(),esdTrack->GetTgl(),esdTrack->Phi(), esdTrack->Pt(),esdTrack->GetITSclusters(0) };
    fTrackingEffHisto->Fill(vecTrackingEff);
    
}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessTPCITS(AliStack* /*const stack*/, AliESDtrack *const esdTrack, AliESDfriendTrack *const esdFriendTrack)
{
  //
  // Match TPC and ITS min-bias tracks
  // at radius between detectors
  //
  if(!esdTrack) return;
  if(!esdFriendTrack) return;
   
  //
  // Propagate tracks to the radius between TPC-ITS
  // using B-field and material budget
  //
  Double_t radius = fCutsRC->GetTPCITSMatchingRadius();
  Double_t mass = esdTrack->GetMass();
  Double_t step=1.0; // cm

  //
  // Propagate TPCinner (reference detector)
  //
  Bool_t isTPCOK=kFALSE;
  AliExternalTrackParam *innerTPC=NULL;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

  //
  // select primaries
  //
  Double_t dcaToVertex = -1;
  if( fCutsRC->GetDCAToVertex2D() ) 
  {
      dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC->GetMaxDCAToVertexXY()/fCutsRC->GetMaxDCAToVertexXY() + dca[1]*dca[1]/fCutsRC->GetMaxDCAToVertexZ()/fCutsRC->GetMaxDCAToVertexZ()); 
  }
  if(fCutsRC->GetDCAToVertex2D() && dcaToVertex > 1) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC->GetMaxDCAToVertexXY()) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC->GetMaxDCAToVertexZ()) return;

  if( (esdTrack->GetTPCNclsIter1()>fCutsRC->GetMinNClustersTPC()) && 
      (esdTrack->GetTPCInnerParam()) &&
      (innerTPC=new AliExternalTrackParam(*(esdTrack->GetTPCInnerParam())))) 
  {
     isTPCOK = AliTracker::PropagateTrackToBxByBz(innerTPC,radius,mass,step,kTRUE);
  }
  if(!isTPCOK) { 
   if(innerTPC) delete innerTPC;  
   return;
  }

  //
  // Propagate ITSouter
  //
  Bool_t isITSOK=kFALSE;
  AliExternalTrackParam *outerITS=NULL;

  if( (esdTrack->GetNcls(0)>fCutsRC->GetMinNClustersITS()) &&
      (esdFriendTrack->GetITSOut()) &&
      (outerITS=new AliExternalTrackParam(*(esdFriendTrack->GetITSOut())))) 
  {
    isITSOK = AliTracker::PropagateTrackToBxByBz(outerITS,radius,mass,step,kTRUE);
  }

  //
  // Fill histograms (TPC reference detector)
  //
  if(isTPCOK)
    FillHistograms(innerTPC,outerITS,isITSOK);

  if(outerITS) delete outerITS;  
  if(innerTPC) delete innerTPC;  
}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessTPCTRD(AliStack* /*const stack*/, AliESDtrack *const esdTrack, AliESDfriendTrack *const esdFriendTrack)
{
  //
  // Match TPC and TRD min-bias tracks
  // at radius between detectors. TPC is the reference detector.
  //
  if(!esdTrack) return;
  if(!esdFriendTrack) return;
  
  //
  // Propagate tracks to the radius between TPC-TRD
  // using B-field and material budget
  //
  Double_t radius = fCutsRC->GetTPCTRDMatchingRadius();
  Double_t mass = esdTrack->GetMass();
  Double_t step=1.0; // cm

  //
  // Propagate TPCouter (reference detector)
  //
  Bool_t isTPCOK=kFALSE;
  AliExternalTrackParam *outerTPC=NULL;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

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

  if( (esdTrack->GetTPCNclsIter1()>fCutsRC->GetMinNClustersTPC()) && 
      (esdFriendTrack->GetTPCOut()) &&
      (outerTPC=new AliExternalTrackParam(*(esdFriendTrack->GetTPCOut())))) 
  {
     isTPCOK = AliTracker::PropagateTrackToBxByBz(outerTPC,radius,mass,step,kTRUE);
  }
  if(!isTPCOK)  { 
    if(outerTPC) delete outerTPC;  
    return;
  }

  //
  // Propagate TRDinner
  //
  Bool_t isTRDOK = kFALSE;
  AliExternalTrackParam *innerTRD=NULL;

  // get TRD track
  AliTRDtrackV1 *trdTrack=NULL; //esdFriendTrack = fESDfriend->GetTrack(itrk);
  TObject *calObject=NULL;
  Int_t icalib = 0;
  while((calObject = esdFriendTrack->GetCalibObject(icalib++))) {
    if(strcmp(calObject->IsA()->GetName(),"AliTRDtrackV1") != 0) continue; // Look for the TRDtrack
    if(!(trdTrack = dynamic_cast<AliTRDtrackV1*>(calObject))) break;
  }

  if( (trdTrack) &&
      (trdTrack->GetNumberOfTracklets()>fCutsRC->GetMinNTrackletsTRD()) &&
      (trdTrack->GetTracklet(0)) &&
      (esdFriendTrack->GetTRDIn()) &&
      (innerTRD = new AliExternalTrackParam(*(esdFriendTrack->GetTRDIn())))) 
  {
    isTRDOK = AliTracker::PropagateTrackToBxByBz(innerTRD,radius,mass,step,kTRUE);
  }

  //
  // Fill histograms (TPC reference detector)
  //
  if(isTPCOK)
    FillHistograms(outerTPC,innerTRD,isTRDOK);

  if(outerTPC) delete outerTPC;  
  if(innerTRD) delete innerTRD;  
}

//_____________________________________________________________________________
void AliPerformanceMatch::FillHistograms(AliExternalTrackParam *const refParam, AliExternalTrackParam *const param, Bool_t isRec) 
{
  //
  // fill performance histograms 
  // (TPC always as reference)
  //
  if(!refParam) return;

  //
  // Deltas (dy,dz,dphi,dtheta,dpt)
  //
  Float_t delta[5] = {0};
  if(param && isRec) {
    delta[0] = param->GetY()-refParam->GetY();
    delta[1] = param->GetZ()-refParam->GetZ();
    delta[2] = TMath::ATan2(param->Py(),param->Px())-TMath::ATan2(refParam->Py(),refParam->Px());
    delta[3] = TMath::ATan2(param->Pz(),param->Pt())-TMath::ATan2(refParam->Pz(),refParam->Pt());
    if(refParam->Pt()) delta[4] = (param->Pt()-refParam->Pt())/refParam->Pt();
  }
  // 
  // Pulls (y,z,snp,tanl,1/pt)
  //
  Float_t sigma[5] = {0};
  Float_t pull[5] = {0};
  if(param && isRec) {
    sigma[0] = TMath::Sqrt(param->GetSigmaY2()+refParam->GetSigmaY2());   
    sigma[1] = TMath::Sqrt(param->GetSigmaZ2()+refParam->GetSigmaZ2());
    sigma[2] = TMath::Sqrt(param->GetSigmaSnp2()+refParam->GetSigmaSnp2());
    sigma[3] = TMath::Sqrt(param->GetSigmaTgl2()+refParam->GetSigmaTgl2());
    sigma[4] = TMath::Sqrt(param->GetSigma1Pt2()+refParam->GetSigma1Pt2());
    if(sigma[0]) pull[0] = delta[0] / sigma[0]; 
    if(sigma[1]) pull[1] = delta[1] / sigma[1]; 
    if(sigma[2]) pull[2] = (param->GetSnp()-refParam->GetSnp()) / sigma[2]; 
    if(sigma[3]) pull[3] = (param->GetTgl()-refParam->GetTgl()) / sigma[3]; 
    if(sigma[4]) pull[4] = (param->OneOverPt()-refParam->OneOverPt()) / sigma[4]; 
  }

  // Fill histograms
  Double_t vResolHisto[11] = {delta[0],delta[1],delta[2],delta[3],delta[4],refParam->GetY(),refParam->GetZ(),refParam->Phi(),refParam->Eta(),refParam->Pt(),isRec};
  fResolHisto->Fill(vResolHisto);

  Double_t vPullHisto[11] = {pull[0],pull[1],pull[2],pull[3],pull[4],refParam->GetY(),refParam->GetZ(),refParam->GetSnp(),refParam->GetTgl(),refParam->OneOverPt(),isRec};
  fPullHisto->Fill(vPullHisto);
}

//_____________________________________________________________________________
void AliPerformanceMatch::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend)
{
  // Process comparison information 
  //
  if(!esdEvent) 
  {
    Error("Exec","esdEvent not available");
    return;
  }
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
  
  // use ESD friends
  if(bUseESDfriend) {
    if(!esdFriend) {
      Error("Exec","esdFriend not available");
      return;
    }
  }

  // trigger
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = esdEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) return; 
  }

  // get TPC event vertex
  const AliESDVertex *vtxESD = esdEvent->GetPrimaryVertexTPC();
  if(vtxESD && (vtxESD->GetStatus()<=0)) return;

  //  Process events
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if(!track) continue;

    AliESDfriendTrack *friendTrack=0;
    if(bUseESDfriend) {
      friendTrack=esdFriend->GetTrack(iTrack);
      if(!friendTrack) continue;
    }

    if(GetAnalysisMode() == 0) ProcessTPCITS(stack,track,friendTrack);
    else if(GetAnalysisMode() == 1) ProcessTPCTRD(stack,track,friendTrack);
    else if(GetAnalysisMode() == 2) ProcessITSTPC(iTrack,esdEvent,stack,track,friendTrack);
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }
}

//_____________________________________________________________________________
TH1F* AliPerformanceMatch::MakeResol(TH2F * his, Int_t integ, Bool_t type, Int_t cut){
  // Create resolution histograms
 
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoII(his,&hism,integ,kTRUE,cut);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliPerformanceMatch::Analyse() {
  // Analyse comparison information and store output histograms
  // in the folder "folderMatch"
  //
  TString selString;
  /*
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH1F *h2=0;
  TH2F *h2D=0;
  */
  TObjArray *aFolderObj = new TObjArray;

  // write results in the folder 
  // TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  // c->cd();

  // char name[256];
  // char title[256];

  if(GetAnalysisMode()==0 || GetAnalysisMode()==1) { 

  fResolHisto->GetAxis(10)->SetRangeUser(1.0,2.0); // only reconstructed
  fPullHisto->GetAxis(10)->SetRangeUser(1.0,2.0);  // only reconstructed
  for(Int_t i=0; i<5; i++) 
  {
    for(Int_t j=5; j<10; j++) 
    {
      //if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
      if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(0.0,0.89); // eta window
      else fResolHisto->GetAxis(8)->SetRangeUser(-1.5,1.49);
      fResolHisto->GetAxis(9)->SetRangeUser(0.1,100.); // pt threshold
      
      selString = "resol";
      AddProjection(aFolderObj, "match", fResolHisto, i, j, &selString);

      /*
      h2D = (TH2F*)fResolHisto->Projection(i,j);
      h = AliPerformanceMatch::MakeResol(h2D,1,0,100);
      snprintf(name,256,"h_res_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fResolHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fResolHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fResolHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      //if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);
      */
      /*
      h = AliPerformanceMatch::MakeResol(h2D,1,1,100);
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
      */
      //
      if(j!=8) fPullHisto->GetAxis(8)->SetRangeUser(0.0,0.89); // eta window
      else  fPullHisto->GetAxis(8)->SetRangeUser(-1.5,1.49); // eta window
      fPullHisto->GetAxis(9)->SetRangeUser(0.1,100.);  // pt threshold
      selString = "pull";
      AddProjection(aFolderObj, "match", fPullHisto, i, j, &selString);

      /*
      h2D = (TH2F*)fPullHisto->Projection(i,j);
      
      h = AliPerformanceMatch::MakeResol(h2D,1,0,100);
      snprintf(name,256,"h_pull_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(resolution)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);

      h = AliPerformanceMatch::MakeResol(h2D,1,1,100);
      snprintf(name,256,"h_mean_pull_%d_vs_%d",i,j);
      h->SetName(name);

      h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
      snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(mean)");
      h->GetYaxis()->SetTitle(title);
      snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
      h->SetTitle(title);

      //if(j==9) h->SetBit(TH1::kLogX);    
      aFolderObj->Add(h);
      */
    }
  }

  //
  // Efficiency plots
  //
  for(Int_t i=5; i<10; i++) 
  {
    if(i!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
    else fResolHisto->GetAxis(8)->SetRangeUser(-1.5,1.49);
    fResolHisto->GetAxis(9)->SetRangeUser(0.1,100.); // pt threshold

    fResolHisto->GetAxis(10)->SetRange(1,fResolHisto->GetAxis(10)->GetNbins()); // all 
    selString = "eff_all";
    AddProjection(aFolderObj, "match", fResolHisto, i, &selString);
    // h = (TH1F*)fResolHisto->Projection(i);

    fResolHisto->GetAxis(10)->SetRange(2,2); // only reconstructed
    selString = "eff_rec";
    AddProjection(aFolderObj, "match", fResolHisto, i, &selString);
    //h2 = (TH1F*)fResolHisto->Projection(i);

    /*
    TH1F* h2c = (TH1F*)h2->Clone();
    h2c->Divide(h2,h,1,1,"B");
 
    snprintf(name,256,"h_eff_%d",i);
    h2c->SetName(name);

    h2c->GetXaxis()->SetTitle(h2c->GetXaxis()->GetTitle());
    h2c->GetYaxis()->SetTitle("efficiency");
    h2c->SetTitle("matching effciency");

    aFolderObj->Add(h2c);
    */
  }

  }
  
  // 
  // TPC efficiency wrt ITS
  //
  if(GetAnalysisMode()==2) { 
    selString = "trackingeff";
    AddProjection(aFolderObj, "match", fTrackingEffHisto, 0, &selString);

    // h = (TH1F*)fTrackingEffHisto->Projection(0);
    // aFolderObj->Add(h);

    for(Int_t i=1; i<7; i++) 
    {
      //
      // 
      // calculate efficiency 
      //

      // all ITS standalone tracks
      fTrackingEffHisto->GetAxis(0)->SetRange(1,fTrackingEffHisto->GetAxis(0)->GetNbins());
      //h = (TH1F*)fTrackingEffHisto->Projection(i);
      selString = "trackingeff_all";
      AddProjection(aFolderObj, "match", fTrackingEffHisto, i, &selString);

      // TPC tracks which has matching with TPC
      fTrackingEffHisto->GetAxis(0)->SetRange(2,2);
      //h2 = (TH1F*)fTrackingEffHisto->Projection(i);
      selString = "trackingeff_tpc";
      AddProjection(aFolderObj, "match", fTrackingEffHisto, i, &selString);

      /*
      TH1F* h2c = (TH1F*)h2->Clone();
      h2c->Divide(h2,h,1,1,"B");
 
      snprintf(name,256,"h_TPC_eff_%d",i);
      h2c->SetName(name);

      h2c->GetXaxis()->SetTitle(h2c->GetXaxis()->GetTitle());
      h2c->GetYaxis()->SetTitle("efficiency");
      h2c->SetTitle("TPC effciency wrt ITS");
      
      aFolderObj->Add(h2c);
      */
    }

  }
  printf("exportToFolder\n");
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(fFolderObj) delete fFolderObj;
  fFolderObj = aFolderObj;  
    aFolderObj=0;  
  
}

//_____________________________________________________________________________
void AliPerformanceMatch::AnalyseFinal() {
  if (!fFolderObj) {   
      printf("AliPerformanceMatch: no projections available to analyse\n");
      return;
  }
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH1F *h2=0;
  TH2F *h2D=0;
  TObjArray *aFolderObj = fFolderObj;

  // write results in the folder 
  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  c->cd();

  char name[256];
  char title[256];

  if(GetAnalysisMode()==0 || GetAnalysisMode()==1) { 

  for(Int_t i=0; i<5; i++) 
  {
    for(Int_t j=5; j<10; j++) 
    {
      snprintf(name,256,"h_tpc_match_resol_%d_%d",i,j);
      h2D = dynamic_cast<TH2F*>(aFolderObj->FindObject(name));
      if (h2D) {
          h = AliPerformanceMatch::MakeResol(h2D,1,0,100);
          snprintf(name,256,"h_res_%d_vs_%d",i,j);
          h->SetName(name);          
          h->GetXaxis()->SetTitle(fResolHisto->GetAxis(j)->GetTitle());
          snprintf(title,256,"%s %s",fResolHisto->GetAxis(i)->GetTitle(),"(resolution)");
          h->GetYaxis()->SetTitle(title);
          snprintf(title,256,"%s vs %s",title,fResolHisto->GetAxis(j)->GetTitle());
          h->SetTitle(title);

          //if(j==9) h->SetBit(TH1::kLogX);    
          aFolderObj->Add(h);

          h = AliPerformanceMatch::MakeResol(h2D,1,1,100);
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
          h=0;
      } 
      //
      snprintf(name,256,"h_tpc_match_pull_%d_%d",i,j);
      h2D = dynamic_cast<TH2F*>(aFolderObj->FindObject(name));      
      if (h2D) {
          h = AliPerformanceMatch::MakeResol(h2D,1,0,100);
          snprintf(name,256,"h_pull_%d_vs_%d",i,j);
          h->SetName(name);

          h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
          snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(resolution)");
          h->GetYaxis()->SetTitle(title);
          snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
          h->SetTitle(title);

          if(j==9) h->SetBit(TH1::kLogX);    
          aFolderObj->Add(h);

          h = AliPerformanceMatch::MakeResol(h2D,1,1,100);
          snprintf(name,256,"h_mean_pull_%d_vs_%d",i,j);
          h->SetName(name);

          h->GetXaxis()->SetTitle(fPullHisto->GetAxis(j)->GetTitle());
          snprintf(title,256,"%s %s",fPullHisto->GetAxis(i)->GetTitle(),"(mean)");
          h->GetYaxis()->SetTitle(title);
          snprintf(title,256,"%s vs %s",title,fPullHisto->GetAxis(j)->GetTitle());
          h->SetTitle(title);

           //if(j==9) h->SetBit(TH1::kLogX);    
           aFolderObj->Add(h);
           h2D=0;
      }
    }
  }

  //
  // Efficiency plots
  //
  for(Int_t i=5; i<10; i++) 
  {
    snprintf(name,256,"h_tpc_match_eff_all_%d",i);
    h = dynamic_cast<TH1F*>(aFolderObj->FindObject(name));    

    
    snprintf(name,256,"h_tpc_match_eff_rec_%d",i);
    h2 = dynamic_cast<TH1F*>(aFolderObj->FindObject(name));    

    if (h && h2) {
        TH1F* h2c = (TH1F*)h2->Clone();
        h2c->Divide(h2,h,1,1,"B");
 
        snprintf(name,256,"h_eff_%d",i);
        h2c->SetName(name);

        h2c->GetXaxis()->SetTitle(h2c->GetXaxis()->GetTitle());
        h2c->GetYaxis()->SetTitle("efficiency");
        h2c->SetTitle("matching effciency");

        aFolderObj->Add(h2c);
        h=0;
        h2=0;
    }
  }

  }
  
  // 
  // TPC efficiency wrt ITS
  //
  if(GetAnalysisMode()==2) { 

    for(Int_t i=1; i<7; i++) 
    {
      //
      // 
      // calculate efficiency 
      //

      // all ITS standalone tracks                  
      snprintf(name,256,"h_tpc_match_trackingeff_all_%d",i);
      h = dynamic_cast<TH1F*>(aFolderObj->FindObject(name));
  

      // TPC tracks which has matching with TPC
      snprintf(name,256,"h_tpc_match_trackingeff_tpc_%d",i);
      h2 = dynamic_cast<TH1F*>(aFolderObj->FindObject(name));

      if (h && h2) {
          TH1F* h2c = (TH1F*)h2->Clone();
          h2c->Divide(h2,h,1,1,"B");
 
          snprintf(name,256,"h_TPC_eff_%d",i);
          h2c->SetName(name);

          h2c->GetXaxis()->SetTitle(h2c->GetXaxis()->GetTitle());
          h2c->GetYaxis()->SetTitle("efficiency");
          h2c->SetTitle("TPC effciency wrt ITS");

          aFolderObj->Add(h2c);
          h=0;
          h2=0;
          h2c=0;
    }

  }
  printf("exportToFolder\n");
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(fFolderObj) delete fFolderObj;
  fFolderObj = aFolderObj;  
    aFolderObj=0;  
  }
  
}


//_____________________________________________________________________________
TFolder* AliPerformanceMatch::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceMatch * comp=this;
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
TFolder* AliPerformanceMatch::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}

//_____________________________________________________________________________
Long64_t AliPerformanceMatch::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  TObjArray* objArrayList = 0;
  objArrayList = new TObjArray();

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliPerformanceMatch* entry = dynamic_cast<AliPerformanceMatch*>(obj);
    if (entry == 0) continue; 
    if (fgMergeTHnSparse) {
        if ((fResolHisto) && (entry->fResolHisto)) { fResolHisto->Add(entry->fResolHisto); }
        if ((fPullHisto) && (entry->fPullHisto)) { fPullHisto->Add(entry->fPullHisto); }
        if ((fTrackingEffHisto) && (entry->fTrackingEffHisto)) { fTrackingEffHisto->Add(entry->fTrackingEffHisto); }
    }
    // the analysisfolder is only merged if present
    if (entry->fFolderObj) { objArrayList->Add(entry->fFolderObj); }

    count++;
  }
  if (fFolderObj) { fFolderObj->Merge(objArrayList); } 
  // to signal that track histos were not merged: reset
  if (!fgMergeTHnSparse) { fResolHisto->Reset(); fPullHisto->Reset(); fTrackingEffHisto->Reset(); }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;
return count;
}
