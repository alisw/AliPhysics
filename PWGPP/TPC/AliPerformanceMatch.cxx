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
// Author: J.Otwinowski   17/10/2009 
// Changes by M.Knichel   22/10/2010
// Changes by J.Salzwedel 14/10/2014
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/Macros/LoadMyLibs.C");
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
#include "TH3.h"
#include "TAxis.h"
#include "TF1.h"

#include "AliPerformanceMatch.h" 
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliVfriendEvent.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTRDtrackV1.h" 
#include "AliTreeDraw.h" 
#include "AliFlatESDTrack.h"

using namespace std;

ClassImp(AliPerformanceMatch)

Bool_t AliPerformanceMatch::fgMergeTHnSparse = kFALSE;
Bool_t AliPerformanceMatch::fgUseMergeTHnSparse = kFALSE;

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch(TRootIOCtor*):
  AliPerformanceObject(),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fTPCConstrain(0),
  fFolderObj(0),
  fAnalysisFolder(0),
  fUseHLT(0),
  h_tpc_match_trackingeff_all_2_3(NULL),
  h_tpc_match_trackingeff_tpc_2_3(NULL),
  h_tpc_match_pull_2_7(NULL),
  h_tpc_match_pull_4_7(NULL),
  h_tpc_match_pull_0_7(NULL),
  h_tpc_match_pull_1_7(NULL),
  h_tpc_match_pull_3_7(NULL),
  h_tpc_constrain_tpc_0_2_3(NULL)

{
  // io constructor
}

//_____________________________________________________________________________
AliPerformanceMatch::AliPerformanceMatch(const Char_t* name, const Char_t* title, Int_t analysisMode, Bool_t hptGenerator, Bool_t useSparse):
  AliPerformanceObject(name,title),
  fResolHisto(0),
  fPullHisto(0),
  fTrackingEffHisto(0),
  fTPCConstrain(0),
  fFolderObj(0),

  // histogram folder 
  fAnalysisFolder(0),
  fUseHLT(0),
  h_tpc_match_trackingeff_all_2_3(NULL),
  h_tpc_match_trackingeff_tpc_2_3(NULL),
  h_tpc_match_pull_2_7(NULL),
  h_tpc_match_pull_4_7(NULL),
  h_tpc_match_pull_0_7(NULL),
  h_tpc_match_pull_1_7(NULL),
  h_tpc_match_pull_3_7(NULL),
  h_tpc_constrain_tpc_0_2_3(NULL)
{
  // named constructor	
  // 
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);
  fUseSparse = useSparse;
  Init();
}

//_____________________________________________________________________________
AliPerformanceMatch::~AliPerformanceMatch()
{
  // destructor

  delete fResolHisto;
  delete fPullHisto;
  delete fTrackingEffHisto;
  delete fTPCConstrain;

  if (fFolderObj && fAnalysisFolder && !fAnalysisFolder->IsOwner()) {
    fFolderObj->Delete();
  } //delete the registered non-sparse histograms

  delete fFolderObj;
  delete fAnalysisFolder;
}

//_____________________________________________________________________________
void AliPerformanceMatch::Init(){

  //
  // Make performance histogrms
  //

  // set pt bins
  
    if(!fUseSparse) fFolderObj = new TObjArray;
    
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 20.;
 
  Double_t *binsPt = 0;

  if (IsHptGenerator())  { 
        ptMax = 100.;
  } 
   binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);


  // Double_t yMin = -0.02, yMax = 0.02;
  // Double_t zMin = -12.0, zMax = 12.0;
  // if(GetAnalysisMode() == 0 || GetAnalysisMode() == 1) { 
  //   yMin = -100.; yMax = 100.; 
  //   zMin = -100.; zMax = 100.; 
  // }

  //
  //init ITS TPC Mactching
  //
  if(GetAnalysisMode()==1){
    // res_y:res_z:res_phi,res_lambda:res_pt:y:z:eta:phi:pt:isRec
    Int_t binsResolHisto[9]={100,100,100,100,100,90,30,nPtBins,2};
    Double_t minResolHisto[9]={-1.,-1.,-0.03,-0.03,-0.2, 0., -1.5, ptMin,0};
    Double_t maxResolHisto[9]={ 1., 1., 0.03, 0.03, 0.2, 2.*TMath::Pi(), 1.5, ptMax,2};
    
    fResolHisto = new THnSparseF("fResolHisto","res_y:res_z:res_phi:res_lambda:res_pt:phi:eta:pt:isRec",9,binsResolHisto,minResolHisto,maxResolHisto);
    fResolHisto->SetBinEdges(7,binsPt);
    
    fResolHisto->GetAxis(0)->SetTitle("y-y_{ref} (cm)");
    fResolHisto->GetAxis(1)->SetTitle("z-z_{ref} (cm)");
    fResolHisto->GetAxis(2)->SetTitle("#phi-#phi_{ref} (rad)");
    fResolHisto->GetAxis(3)->SetTitle("#lambda-#lambda_{ref} (rad)");
    fResolHisto->GetAxis(4)->SetTitle("(p_{T}/p_{Tref}-1)");
    fResolHisto->GetAxis(5)->SetTitle("#phi_{ref} (rad)");
    fResolHisto->GetAxis(6)->SetTitle("#eta_{ref}");
    fResolHisto->GetAxis(7)->SetTitle("p_{Tref} (GeV/c)");
    fResolHisto->GetAxis(8)->SetTitle("isReconstructed");
    fResolHisto->Sumw2();
    
    //pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:y:z:snp:tgl:1pt:isRec
    Int_t binsPullHisto[9]={100,100,100,100,100,90,30,nPtBins,2};
    Double_t minPullHisto[9]={-5.,-5.,-5.,-5.,-5., 0.,-1.5, ptMin,0};
    Double_t maxPullHisto[9]={ 5., 5., 5., 5., 5., 2.*TMath::Pi(), 1.5, ptMax,2};
      if(fUseSparse){
        fPullHisto = new THnSparseF("fPullHisto","pull_y:pull_z:pull_snp:pull_tgl:pull_1pt:phi:eta:1pt:isRec",9,binsPullHisto,minPullHisto,maxPullHisto);
        fPullHisto->GetAxis(0)->SetTitle("(y-y_{ref})/#sigma");
        fPullHisto->GetAxis(1)->SetTitle("(z-z_{ref})/#sigma");
        fPullHisto->GetAxis(2)->SetTitle("(sin#phi-sin#phi_{ref})/#sigma");
        fPullHisto->GetAxis(3)->SetTitle("(tan#lambda-tan#lambda_{ref})/#sigma");
        fPullHisto->GetAxis(4)->SetTitle("(p_{Tref}/p_{T}-1)/#sigma");
        fPullHisto->GetAxis(5)->SetTitle("#phi_{ref} (rad)");
        fPullHisto->GetAxis(6)->SetTitle("#eta_{ref}");
        fPullHisto->GetAxis(7)->SetTitle("1/p_{Tref} (GeV/c)^{-1}");
        fPullHisto->GetAxis(8)->SetTitle("isReconstructed");
        fPullHisto->Sumw2();
      }
      else{
          h_tpc_match_pull_2_7 = new TH2D("h_tpc_match_pull_2_7","",binsPullHisto[2],minPullHisto[2],maxPullHisto[2],binsPullHisto[7],minPullHisto[7],maxPullHisto[7]);
          h_tpc_match_pull_2_7->SetXTitle("(sin#phi-sin#phi_{ref})/#sigma");
          h_tpc_match_pull_2_7->SetYTitle("1/p_{Tref} (GeV/c)^{-1}");
          
          h_tpc_match_pull_4_7 = new TH2D("h_tpc_match_pull_2_7","",binsPullHisto[4],minPullHisto[4],maxPullHisto[4],binsPullHisto[7],minPullHisto[7],maxPullHisto[7]);
          h_tpc_match_pull_4_7->SetXTitle("(p_{Tref}/p_{T}-1)/#sigma");
          h_tpc_match_pull_4_7->SetYTitle("1/p_{Tref} (GeV/c)^{-1}");

          h_tpc_match_pull_0_7 = new TH2D("h_tpc_match_pull_2_7","",binsPullHisto[0],minPullHisto[0],maxPullHisto[0],binsPullHisto[7],minPullHisto[7],maxPullHisto[7]);
          h_tpc_match_pull_0_7->SetXTitle("(y-y_{ref})/#sigma");
          h_tpc_match_pull_0_7->SetYTitle("1/p_{Tref} (GeV/c)^{-1}");

          h_tpc_match_pull_1_7 = new TH2D("h_tpc_match_pull_2_7","",binsPullHisto[1],minPullHisto[1],maxPullHisto[1],binsPullHisto[7],minPullHisto[7],maxPullHisto[7]);
          h_tpc_match_pull_1_7->SetXTitle("1/p_{Tref} (GeV/c)^{-1}");
          h_tpc_match_pull_1_7->SetYTitle("1/p_{Tref} (GeV/c)^{-1}");

          h_tpc_match_pull_3_7 = new TH2D("h_tpc_match_pull_2_7","",binsPullHisto[3],minPullHisto[3],maxPullHisto[3],binsPullHisto[7],minPullHisto[7],maxPullHisto[7]);
          h_tpc_match_pull_3_7->SetXTitle("(tan#lambda-tan#lambda_{ref})/#sigma");
          h_tpc_match_pull_3_7->SetYTitle("1/p_{Tref} (GeV/c)^{-1}");
      
          fFolderObj->Add(h_tpc_match_pull_2_7);
          fFolderObj->Add(h_tpc_match_pull_4_7);
          fFolderObj->Add(h_tpc_match_pull_0_7);
          fFolderObj->Add(h_tpc_match_pull_1_7);
          fFolderObj->Add(h_tpc_match_pull_3_7);
      }
  }
  
  //
  //TPC  ITS  Mactching
  // 
  if(GetAnalysisMode()==0){
    // -> has match:y:z:snp:tgl:phi:pt:ITSclusters
    Int_t binsTrackingEffHisto[5]    = { 2,   90,       nPtBins, 30,    7   };
    Double_t minTrackingEffHisto[5]  = {-0.5, 0.,          ptMin, -1.5, -0.5 };
    Double_t maxTrackingEffHisto[5]  = { 1.5, 2*TMath::Pi(), ptMax, 1.5,  6.5 };
      if(fUseSparse){
        fTrackingEffHisto = new THnSparseF("fTrackingEffHisto","has match:phi:pt:eta:ITSclusters",5,binsTrackingEffHisto,minTrackingEffHisto,maxTrackingEffHisto);
        fTrackingEffHisto->SetBinEdges(2,binsPt);
        fTrackingEffHisto->GetAxis(0)->SetTitle("IsMatching");
        fTrackingEffHisto->GetAxis(1)->SetTitle("phi (rad)");
        fTrackingEffHisto->GetAxis(2)->SetTitle("p_{T}");
        fTrackingEffHisto->GetAxis(3)->SetTitle("eta");
        fTrackingEffHisto->GetAxis(4)->SetTitle("number of ITS clusters");
              fTrackingEffHisto->Sumw2();
      }
      else{
          h_tpc_match_trackingeff_all_2_3 = new TH2D("h_tpc_match_trackingeff_all_2_3","",binsTrackingEffHisto[2],minTrackingEffHisto[2],maxTrackingEffHisto[2],binsTrackingEffHisto[3],minTrackingEffHisto[3],maxTrackingEffHisto[3]);
          h_tpc_match_trackingeff_all_2_3->SetXTitle("p_{T}");
          h_tpc_match_trackingeff_all_2_3->SetYTitle("eta");
          
          h_tpc_match_trackingeff_tpc_2_3 = new TH2D("h_tpc_match_trackingeff_tpc_2_3","",binsTrackingEffHisto[2],minTrackingEffHisto[2],maxTrackingEffHisto[2],binsTrackingEffHisto[3],minTrackingEffHisto[3],maxTrackingEffHisto[3]);
          h_tpc_match_trackingeff_tpc_2_3->SetXTitle("p_{T}");
          h_tpc_match_trackingeff_tpc_2_3->SetYTitle("eta");
          
          fFolderObj->Add(h_tpc_match_trackingeff_all_2_3);
          fFolderObj->Add(h_tpc_match_trackingeff_tpc_2_3);
      }
  }

  //
  //TPC constraining to vertex
  //
  if(GetAnalysisMode()==2){
    //initilization of fTPCConstrain
    Int_t  binsTPCConstrain[4] = {100, 90,            nPtBins, 30};
    Double_t minTPCConstrain[4] = {-5, 0,             ptMin,   -1.5};
    Double_t maxTPCConstrain[4] = {5,  2*TMath::Pi(), ptMax,  1.5};
    
    if(fUseSparse){
        fTPCConstrain = new THnSparseF("fTPCConstrain","pull_phi:phi:pt:eta",4, binsTPCConstrain,minTPCConstrain,maxTPCConstrain);
        fTPCConstrain->SetBinEdges(2,binsPt);
        fTPCConstrain->GetAxis(0)->SetTitle("(#phi-#phi_{ref})/#sigma");
        fTPCConstrain->GetAxis(1)->SetTitle("phi (rad)");
        fTPCConstrain->GetAxis(2)->SetTitle("p_{T}");
        fTPCConstrain->GetAxis(3)->SetTitle("eta");
        fTPCConstrain->Sumw2();
      }
    else{
        h_tpc_constrain_tpc_0_2_3 = new TH3D("h_tpc_constrain_tpc_0_2_3","",binsTPCConstrain[0],minTPCConstrain[0],maxTPCConstrain[0],binsTPCConstrain[2],minTPCConstrain[2],maxTPCConstrain[2],binsTPCConstrain[3],minTPCConstrain[3],maxTPCConstrain[3]);
        h_tpc_constrain_tpc_0_2_3->SetXTitle("(#phi-#phi_{ref})/#sigma");
        h_tpc_constrain_tpc_0_2_3->SetYTitle("p_{T}");
        h_tpc_constrain_tpc_0_2_3->SetZTitle("eta");
        
        fFolderObj->Add(h_tpc_constrain_tpc_0_2_3);
        
    }
  }

  // init folder
  fAnalysisFolder = CreateFolder("folderMatch","Analysis Matching Folder");
  
   // save merge status in object
  fMergeTHnSparseObj = fgMergeTHnSparse;

  delete [] binsPt;

}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessITSTPC(Int_t iTrack, AliVEvent *const vEvent, AliMCEvent* /*const mcev*/, AliVTrack *const vTrack)
{
  //
  // addition to standard analysis - check if ITS stand-alone tracks have a match in the TPC
  // Origin: A. Kalwait
  // Modified: J. Otwinowski
  if(!vEvent && !vTrack) return;

  // ITS stand alone tracks with SPD layers 
  if(!(vTrack->GetStatus() & AliVTrack::kITSpureSA)) return;
  if(!(vTrack->GetStatus() & AliVTrack::kITSrefit)) return;
  if(vTrack->GetNcls(0)<4) return;
  if(!vTrack->HasPointOnITSLayer(0) && !vTrack->HasPointOnITSLayer(1)) return;

  const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
  AliESDtrack* tpcTrack2 = NULL;
  Bool_t hasMatch = kFALSE;
  for (Int_t jTrack = 0; jTrack < vEvent->GetNumberOfTracks(); jTrack++) {
    // loop for all those tracks and check if the corresponding TPC track is found
    if (jTrack==iTrack) continue;
    AliVTrack *vTrackTPC = (AliVTrack*)vEvent->GetTrack(jTrack);
    if (!vTrackTPC) continue;
    if (!vTrackTPC->GetTPCInnerParam()) continue;
    if(!(vTrackTPC->GetStatus() & AliVTrack::kTPCrefit)) continue;
    
    // TPC nClust/track after first tracking pass
    // if(vTrackTPC->GetTPCNclsIter1()<fCutsRC.GetMinNClustersTPC()) continue;
    tpcTrack2 = AliESDtrackCuts::GetTPCOnlyTrackFromVEvent(vEvent, jTrack);
    if(!tpcTrack2) continue;
    //TODO: this will only work offline with ESD tracks! needs some work to make it work online as well.
    if(!tpcTrack2->RelateToVVertex(vVertex,vEvent->GetMagneticField(),100.)) { delete tpcTrack2; tpcTrack2=0; continue; } 
    
    if(!fCutsRC.AcceptVTrack(tpcTrack2)) { delete tpcTrack2; tpcTrack2=0; continue; }
    // check matching
    if (TMath::Abs(vTrack->GetY() - tpcTrack2->GetY()) > 3) { delete tpcTrack2; tpcTrack2=0; continue; }
    if (TMath::Abs(vTrack->GetSnp() - tpcTrack2->GetSnp()) > 0.2) { delete tpcTrack2; tpcTrack2=0; continue; }
    if (TMath::Abs(vTrack->GetTgl() - tpcTrack2->GetTgl()) > 0.2) { delete tpcTrack2; tpcTrack2=0; continue; }
    
    hasMatch=kTRUE;
    break;
  }
  
  FillHistograms(tpcTrack2,vTrack,hasMatch);     
  delete tpcTrack2;
}

//_____________________________________________________________________________
void AliPerformanceMatch::ProcessTPCITS(AliMCEvent* /*const mcev*/, AliVEvent *const vEvent,AliVTrack *const vTrack)
{
  
  //
  // Match TPC and ITS min-bias tracks
  // at radius between detectors
  //
  if(!vTrack) return;
  if(!vEvent) return;
  AliFlatESDTrack *flatTrack = dynamic_cast<AliFlatESDTrack*>(vTrack);

  Bool_t isTPC = kFALSE;
  Bool_t isMatch = kFALSE;

  if(!flatTrack){
      if(vTrack->Charge()==0) return;
      if(!vTrack->GetTPCInnerParam()) return;
      if(!fCutsRC.AcceptVTrack(vTrack)) return;
  }
  else{
      if(!fCutsRC.AcceptFTrack(vTrack,vEvent)) return;
  }
  if(!(vTrack->GetStatus()&AliVTrack::kTPCrefit)) return;

  isTPC = kTRUE;
  if( (vTrack->GetStatus()&AliVTrack::kITSrefit))
    isMatch = kTRUE;
    
    AliExternalTrackParam trackParams;
    vTrack->GetTrackParam(trackParams);
    AliExternalTrackParam *etpTrack = &trackParams;
    
  if(isTPC){
    Double_t vecTrackingEff[5] = { static_cast<Double_t>(isMatch),etpTrack->Phi(), etpTrack->Pt(),etpTrack->Eta(),static_cast<Double_t>(vTrack->GetITSclusters(0)) };
    if(fUseSparse) fTrackingEffHisto->Fill(vecTrackingEff);
    else{
        if(vecTrackingEff[0] > -0.5) h_tpc_match_trackingeff_all_2_3->Fill(vecTrackingEff[2],vecTrackingEff[3]);
        if(vecTrackingEff[0] > 0.5) h_tpc_match_trackingeff_tpc_2_3->Fill(vecTrackingEff[2],vecTrackingEff[3]);
    }
      
  }
}

//_____________________________________________________________________________
/*void AliPerformanceMatch::ProcessTPCTRD(AliMCEvent* , AliVTrack *const vTrack, AliVfriendTrack *const vFriendTrack)
{
  return;
}*/

void AliPerformanceMatch::ProcessTPCConstrain(AliMCEvent* /*const mcev*/, AliVEvent *const vEvent, AliVTrack *const vTrack){
  //
  // Contrain TPC inner track to the vertex
  // then compare to the global tracks
  //
    if(!vEvent && !vTrack) return;
    if(!vTrack->GetTPCInnerParam()) return;
  
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    AliFlatESDTrack *flatTrack = dynamic_cast<AliFlatESDTrack*>(vTrack);

    if(!flatTrack){
        if(vTrack->Charge()==0) return;
        if(!fCutsRC.AcceptVTrack(vTrack)) return;
    }
    else{
        if(!fCutsRC.AcceptFTrack(vTrack,vEvent)) return;
    }
    if(!(vTrack->GetStatus()&AliVTrack::kITSrefit)) return;
    if(!(vTrack->GetStatus()&AliVTrack::kTPCrefit)) return;

    AliExternalTrackParam trackParams, trackParamsTPCInner;
    vTrack->GetTrackParam(trackParams);
    AliExternalTrackParam *etpTrack = &trackParams;
    Double_t x[3]; etpTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = kFALSE;

    if (vTrack->GetTrackParamTPCInner(trackParamsTPCInner)<0) { return; }
    AliExternalTrackParam * TPCinner = &trackParamsTPCInner;
  
    isOK = TPCinner->ConstrainToVertex(vVertex, b);

    // transform to the track reference frame 
    isOK = TPCinner->Rotate(etpTrack->GetAlpha());
    if(!isOK) return;
    isOK = TPCinner->PropagateTo(etpTrack->GetX(),vEvent->GetMagneticField());
    if(!isOK) return;

    Double_t sigmaPhi=0,deltaPhi=0,pullPhi=0;
    deltaPhi = TPCinner->GetSnp() - etpTrack->GetSnp();
    sigmaPhi = TMath::Sqrt(vTrack->GetSigmaSnp2()+TPCinner->GetSigmaSnp2());
    if(sigmaPhi!=0)
    pullPhi = deltaPhi/sigmaPhi;

    Double_t vTPCConstrain[4] = {pullPhi,etpTrack->Phi(),etpTrack->Pt(),etpTrack->Eta()};
    if(fUseSparse) fTPCConstrain->Fill(vTPCConstrain);
    else {
        h_tpc_constrain_tpc_0_2_3->Fill(vTPCConstrain[0],vTPCConstrain[2],vTPCConstrain[3]);
    }

    return;
}
//_____________________________________________________________________________
void AliPerformanceMatch::FillHistograms(AliVTrack *const refParamVTrack, AliVTrack *const paramVTrack, Bool_t isRec) 
{
  //
  // fill performance histograms 
  // (TPC always as reference)
  //

  
  if(!refParamVTrack || !paramVTrack || !isRec) return;

  //this should be done with a copy ctor?
  AliExternalTrackParam refParam_;
  AliExternalTrackParam param_;
  refParam_.CopyFromVTrack(refParamVTrack);
  param_.CopyFromVTrack(paramVTrack);

  AliExternalTrackParam *refParam = &refParam_;
  AliExternalTrackParam *param = &param_;
  
  if(!refParam) return;
  if(!param) return;
  if(!isRec) return;
  
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
  Double_t vResolHisto[9] = {delta[0],delta[1],delta[2],delta[3],delta[4],refParam->Phi(),refParam->Eta(),refParam->Pt(),static_cast<Double_t>(isRec)};
  Double_t vPullHisto[9] = {pull[0],pull[1],pull[2],pull[3],pull[4],refParam->Phi(),refParam->Eta(),refParam->OneOverPt(),static_cast<Double_t>(isRec)};
    if(fabs(vPullHisto[4])<5){
        if(fUseSparse){
            fResolHisto->Fill(vResolHisto);
            fPullHisto->Fill(vPullHisto);
        }
        else {
            if(vPullHisto[6] > 0. && vPullHisto[6] < 1.49)
            if(vPullHisto[7] > 0.01 && vPullHisto[7] < 10)
            if(vPullHisto[8] > 1.0 && vPullHisto[8] < 2.0){
                if (h_tpc_match_pull_2_7) h_tpc_match_pull_2_7->Fill(vPullHisto[2],vPullHisto[7]);
                if (h_tpc_match_pull_4_7) h_tpc_match_pull_4_7->Fill(vPullHisto[4],vPullHisto[7]);
                if (h_tpc_match_pull_0_7) h_tpc_match_pull_0_7->Fill(vPullHisto[0],vPullHisto[7]);
                if (h_tpc_match_pull_1_7) h_tpc_match_pull_1_7->Fill(vPullHisto[1],vPullHisto[7]);
                if (h_tpc_match_pull_3_7) h_tpc_match_pull_3_7->Fill(vPullHisto[3],vPullHisto[7]);
            }
            
        }
    }
}
//_____________________________________________________________________________
void AliPerformanceMatch::Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent */*const vfriendEvent*/, const Bool_t bUseMC, const Bool_t /*bUseVfriend*/)
{
  // Process comparison information 
  //
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

  // trigger
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = vEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) return; 
  }

  // get TPC event vertex
    AliESDVertex vertex;
    vEvent->GetPrimaryVertex(vertex);
    if(!(vertex.GetStatus())) return;

  //  Process events
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  {
    AliVTrack *track = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!track) continue;

    /*AliESDfriendTrack *friendTrack=0;
    if(bUseESDfriend) {
      friendTrack=esdFriend->GetTrack(iTrack);
      if(!friendTrack) continue;
      }*/

    if(GetAnalysisMode() == 0){
      if(!IsUseTOFBunchCrossing()){
	ProcessTPCITS(mcEvent,vEvent,track);
      }
      else
	if( track->GetTOFBunchCrossing(vEvent->GetMagneticField())==0) {
	  ProcessTPCITS(mcEvent,vEvent,track);
	}
    }
    /* else if(GetAnalysisMode() == 2) ProcessTPCTRD(mcev,track,friendTrack);*/
    else if(GetAnalysisMode() == 1) {ProcessITSTPC(iTrack,vEvent,mcEvent,track);}
    else if(GetAnalysisMode() == 2){
      if(!IsUseTOFBunchCrossing()){
	ProcessTPCConstrain(mcEvent,vEvent,track);
      }
      else
	if( track->GetTOFBunchCrossing(vEvent->GetMagneticField())==0) {
	  ProcessTPCConstrain(mcEvent,vEvent,track);
	}
    }
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }
  }

}

//_____________________________________________________________________________
TH1F* AliPerformanceMatch::MakeResol(TH2F * his, Int_t integ, Bool_t type, Int_t cut){
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
void AliPerformanceMatch::Analyse() {
  // Analyse comparison information and store output histograms
  // in the folder "folderMatch"
  //
  /*
  TH1::AddDirectory(kFALSE);
  TH1F *h=0;
  TH1F *h2=0;
  TH2F *h2D=0;
  */
    if(fUseSparse){
        TString selString;
        TObjArray *aFolderObj = new TObjArray;

      // write results in the folder 
      // TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
      // c->cd();

      // char name[256];
      // char title[256];

      if(GetAnalysisMode()==1) { 

      fResolHisto->GetAxis(8)->SetRangeUser(1.0,2.0); // only reconstructed
      fPullHisto->GetAxis(8)->SetRangeUser(1.0,2.0);  // only reconstructed
      for(Int_t i=0; i<5; i++) 
        {
          for(Int_t j=5; j<8; j++) 
        {
          //if(j!=8) fResolHisto->GetAxis(8)->SetRangeUser(-0.9,0.89); // eta window
          if(j!=6) fResolHisto->GetAxis(6)->SetRangeUser(0.0,1.49); // eta window
          else fResolHisto->GetAxis(6)->SetRangeUser(-1.5,1.49);
          fResolHisto->GetAxis(7)->SetRangeUser(0.01,10.); // pt threshold
          
          selString = "resol";
          AddProjection(aFolderObj, "match", fResolHisto, i, j, &selString);
          
          
          if(j!=6) fPullHisto->GetAxis(6)->SetRangeUser(0.0,1.49); // eta window
          else  fPullHisto->GetAxis(6)->SetRangeUser(-1.5,1.49); // eta window
          fPullHisto->GetAxis(7)->SetRangeUser(0.01,10.);  // pt threshold
          selString = "pull";
          AddProjection(aFolderObj, "match", fPullHisto, i, j, &selString);
          
        }
        }
      }
      // 
      // TPC efficiency wrt ITS
      //
      if(GetAnalysisMode()==0) { 
        selString = "trackingeff";
        AddProjection(aFolderObj, "match", fTrackingEffHisto, 0, &selString);
        
        for(Int_t i=1; i<5; i++) 
          {
        // all ITS standalone tracks
        fTrackingEffHisto->GetAxis(0)->SetRange(1,fTrackingEffHisto->GetAxis(0)->GetNbins());
        selString = "trackingeff_all";
        AddProjection(aFolderObj, "match", fTrackingEffHisto, i, 3,&selString);
        
          // TPC tracks which has matching with TPC
        fTrackingEffHisto->GetAxis(0)->SetRange(2,2);
        selString = "trackingeff_tpc";
        AddProjection(aFolderObj, "match", fTrackingEffHisto, i, 3,&selString);
          }
      }

      //
      //TPC constrained to vertex
      //
      if(GetAnalysisMode()==2) { 
        selString = "tpc";
        //    for(Int_t i=1; i<4; i++)
        AddProjection(aFolderObj, "constrain", fTPCConstrain, 0, 1, 2, &selString);
        AddProjection(aFolderObj, "constrain", fTPCConstrain, 0, 1, 3, &selString);
        AddProjection(aFolderObj, "constrain", fTPCConstrain, 0, 2, 3, &selString);
      }
      
      printf("exportToFolder\n");
      fAnalysisFolder = ExportToFolder(aFolderObj);
      
      // delete only TObjArray
      if(fFolderObj) delete fFolderObj;
      fFolderObj = aFolderObj;  
      aFolderObj=0;  
    }
    else {
        printf("exportToFolder\n");
        fAnalysisFolder = ExportToFolder(fFolderObj);
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
  
  Bool_t merge = ((fgUseMergeTHnSparse && fgMergeTHnSparse) || (!fgUseMergeTHnSparse && fMergeTHnSparseObj));

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
    if (merge) {
        if ((fResolHisto) && (entry->fResolHisto)) { fResolHisto->Add(entry->fResolHisto); }
        if ((fPullHisto) && (entry->fPullHisto)) { fPullHisto->Add(entry->fPullHisto); }
        if ((fTrackingEffHisto) && (entry->fTrackingEffHisto)) { fTrackingEffHisto->Add(entry->fTrackingEffHisto); }

        if ((fTPCConstrain) && (entry->fTPCConstrain)) { fTPCConstrain->Add(entry->fTPCConstrain); }
    }
    // the analysisfolder is only merged if present
    if (entry->fFolderObj) { objArrayList->Add(entry->fFolderObj); }

    count++;
  }
  if (fFolderObj) { fFolderObj->Merge(objArrayList); } 
  // to signal that track histos were not merged: reset
  if (!merge) {
      if(fResolHisto) fResolHisto->Reset();
      if(fPullHisto) fPullHisto->Reset();
      if(fTrackingEffHisto) fTrackingEffHisto->Reset();
      if(fTPCConstrain) fTPCConstrain->Reset();
  }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;
return count;
}

void AliPerformanceMatch::ResetOutputData(){

    if(fUseSparse){
        if(fResolHisto) fResolHisto->Reset("ICE");
        if(fPullHisto) fPullHisto->Reset("ICE");
        if(fTrackingEffHisto) fTrackingEffHisto->Reset("ICE");
        if(fTPCConstrain) fTPCConstrain->Reset("ICE");
    }
    else{
        if(h_tpc_match_trackingeff_all_2_3) h_tpc_match_trackingeff_all_2_3->Reset("ICE");
        if(h_tpc_match_trackingeff_tpc_2_3) h_tpc_match_trackingeff_tpc_2_3->Reset("ICE");
        if(h_tpc_match_pull_2_7) h_tpc_match_pull_2_7->Reset("ICE");
        if(h_tpc_match_pull_4_7) h_tpc_match_pull_4_7->Reset("ICE");
        if(h_tpc_match_pull_0_7) h_tpc_match_pull_0_7->Reset("ICE");
        if(h_tpc_match_pull_1_7) h_tpc_match_pull_1_7->Reset("ICE");
        if(h_tpc_match_pull_3_7) h_tpc_match_pull_3_7->Reset("ICE");
        if(h_tpc_constrain_tpc_0_2_3) h_tpc_constrain_tpc_0_2_3->Reset("ICE");
    }

}

//_____________________________________________________________________________
TCollection* AliPerformanceMatch::GetListOfDrawableObjects() 
{
  TObjArray* tmp = fFolderObj;
  fFolderObj = NULL;
  if (fAnalysisFolder) { fAnalysisFolder->SetOwner(kFALSE); }
  return tmp;
}

