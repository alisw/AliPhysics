// charge dependent mixed harmonics
// author: Jaap Onderwaater, j.onderwaater@gsi.de
//         2014/Oct

#include <iostream>
#include "AliHistogramManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsAxes.h"
#include "AliQnCorrectionsSteps.h"
///#include "AliCMEVarManager.h"
//#include "AliCMEReducedVarManager.h"
#include "AliChargeOnePwrtRP.h"
#include "AliChargeTwoPwrtRP.h"
#include "AliCMEAnalysisCuts.h"
#include <THashList.h>
#include <TClonesArray.h>
#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedAnalysisTaskTwoPwrtRP.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedEventInfo.h"
#include "AliReducedFMDInfo.h"
#include "TCutG.h"
#include "TRandom3.h"

#define FILL AliReducedVarManager

using std::cout;
using std::endl;


ClassImp(AliReducedAnalysisTaskTwoPwrtRP)



  //_________________________________________________________________________________
  AliReducedAnalysisTaskTwoPwrtRP::AliReducedAnalysisTaskTwoPwrtRP() :
    AliReducedAnalysisTaskSE(),
    fInitialized(kFALSE),
    fFriendPath(""),
    fListHistos(),
    fListHistosQn(),
    fOffset(0),
    fNfiles(1),
    fQvecTPC(0x0),
    fQvecVZEROA(0x0),
    fQvecVZEROC(0x0),
    fQvecFMDA(0x0),
    fQvecFMDC(0x0),
    fNOnePwrtRP(0),
    fNTwoPwrtRP(0),
    fNCtbins(-1),
    fNPtbins(-1),
    fNEtabins(-1),
    fNDPtbins(-1),
    fNDEtabins(-1),
    fNMPtbins(-1),
    fNMEtabins(-1),
    fMinHarmonic(1000),
    fMaxHarmonic(0),
    fTracksTPCstandalone(kTRUE),
    fRemoveOutliers(kFALSE),
    fNevents(0)
{
  //
  // Default constructor
  //
}



//_________________________________________________________________________________
AliReducedAnalysisTaskTwoPwrtRP::AliReducedAnalysisTaskTwoPwrtRP(const char* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fInitialized(kFALSE),
  fFriendPath(""),
  fListHistos(),
  fListHistosQn(),
  fOffset(0),
  fNfiles(1),
  fQvecTPC(0x0),
  fQvecVZEROA(0x0),
  fQvecVZEROC(0x0),
  fQvecFMDA(0x0),
  fQvecFMDC(0x0),
  fNOnePwrtRP(0),
  fNTwoPwrtRP(0),
  fNCtbins(-1),
  fNPtbins(-1),
  fNEtabins(-1),
  fNDPtbins(-1),
  fNDEtabins(-1),
  fNMPtbins(-1),
  fNMEtabins(-1),
  fMinHarmonic(1000),
  fMaxHarmonic(0),
  fTracksTPCstandalone(kTRUE),
  fRemoveOutliers(kFALSE),
  fNevents(0)
{
  //
  // Constructor
  //

  for(Int_t i=0; i<10; i++) fTrackingQuality[i] = 0x0;
  for(Int_t i=0; i<AliChargeOnePwrtRP::Nqvectors; i++) fEPcor[i] = 0x0;
  //fTask = new AliReducedAnalysisTaskTwoPwrtRP();

  //DefineInput(0,TChain::Class());
  //localHelper::GetInstance()->DefineInput(0,AliReducedEvent::Class());
  //localHelper::GetInstance()->DefineInput(1,TList::Class());
  //DefineInput(0,TTree::Class());
  //localHelper::GetInstance()->DefineOutput(1, TList::Class());   // Calibration histograms
  //localHelper::GetInstance()->DefineOutput(2, TList::Class());   // QA histograms

  //fListHistos.SetName("CalibrationHistograms");
  //fListHistos.SetOwner(kFALSE);
}




//_______________________________________________________________________________
AliReducedAnalysisTaskTwoPwrtRP::~AliReducedAnalysisTaskTwoPwrtRP()
{
  //
  // De-Constructor
  //
  for(Int_t g=0; g<7; g++) delete fAxes[g];
  for(Int_t ic=0; ic<fNOnePwrtRP; ic++) if(fFlow[ic]) delete fFlow[ic];
  for(Int_t icom=0; icom<fNTwoPwrtRP; icom++)  if(fChargeCorrelation[icom]) delete fChargeCorrelation[icom];
  for(Int_t ih=0; ih<AliChargeOnePwrtRP::Nqvectors; ih++) if(fEPcor[ih]) delete fEPcor[ih];




}





//_________________________________________________________________________________
void AliReducedAnalysisTaskTwoPwrtRP::Init()
{
  //
  // Add all histogram manager histogram lists to the output TList (called in AddTask)
  //

  if (!fListHistos.IsEmpty()) return; //already initialised

  fListHistos.SetName("ChargeCorrelations");
  fListHistos.SetOwner();
  fListHistosQn.SetName("QnCorrelations");
  fListHistosQn.SetOwner();

  TAxis trackaxis[4];
  Double_t  tCtbins[][2]  = {{ 0.0, 3}, {10.0, 2}, {80.0,7}};
  Double_t  tPtbins[][2]  = {{ 0.0, 4}, {2.0, 8}, {3.0, 2}, {5.0, 1}};
  Double_t tEtabins[][2]  = {{-0.8, 2}, {0.8, 8}};
  Double_t tPhibins[][2]  = {{ 0.0, 2}, {2.*TMath::Pi(), 50}};

  Double_t tdEtabins[][2]  = {{0.0, 2}, {1.6, 16}};
  Double_t tdPhibins[][2]  = {{-TMath::Pi(), 2}, {TMath::Pi(), 100}};

  trackaxis[0] =  AliQnCorrectionsAxes::MakeAxis(tCtbins);
  trackaxis[1] =  AliQnCorrectionsAxes::MakeAxis(tPtbins);
  trackaxis[2] =  AliQnCorrectionsAxes::MakeAxis(tEtabins);
  trackaxis[3] =  AliQnCorrectionsAxes::MakeAxis(tPhibins);
  fTrackDistribution=AliHistogramManager::CreateHistogram("tracks","tracks",4,trackaxis);

  TAxis thndphideta[2];
  thndphideta[0] =  AliQnCorrectionsAxes::MakeAxis(tdPhibins);
  thndphideta[1] =  AliQnCorrectionsAxes::MakeAxis(tdEtabins);

  TAxis axis[7];
  for(Int_t i=0; i<= fAxes[0]->GetNbins(); i++) {SetBinEdgeCt(  i, fAxes[0]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[1]->GetNbins(); i++) {SetBinEdgePt(  i, fAxes[1]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[2]->GetNbins(); i++) {SetBinEdgeMPt( i, fAxes[2]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[3]->GetNbins(); i++) {SetBinEdgeDPt( i, fAxes[3]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[4]->GetNbins(); i++) {SetBinEdgeEta( i, fAxes[4]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[5]->GetNbins(); i++) {SetBinEdgeMEta(i, fAxes[5]->GetBinUpEdge(i)); };
  for(Int_t i=0; i<= fAxes[6]->GetNbins(); i++) {SetBinEdgeDEta(i, fAxes[6]->GetBinUpEdge(i)); };

  TString chargecombisname[4] = {"pp","mp","pm","mm"};
  TString harmonicsname[4] = {"11","22","33","44"};
  TString centname[9] = {"cent00x05","cent05x10","cent10x20","cent20x30","cent30x40","cent40x50","cent50x60","cent60x70","cent70x80"};
  for(Int_t i=0; i<4; i++){
    for(Int_t j=0; j<4; j++){
      fProfileCent[j][i][0] = new TProfile(Form("c2_%s_%s_cent",harmonicsname[j].Data(),chargecombisname[i].Data()),Form("c2_%s_%s_cent",harmonicsname[j].Data(),chargecombisname[i].Data()), fAxes[0]->GetNbins(), fCtbinning);
      fProfileCent[j][i][1] = new TProfile(Form("s2_%s_%s_cent",harmonicsname[j].Data(),chargecombisname[i].Data()),Form("s2_%s_%s_cent",harmonicsname[j].Data(),chargecombisname[i].Data()), fAxes[0]->GetNbins(), fCtbinning);
      for(Int_t k=0; k<9;k++){
        fProfileMPt[k][j][i][0] = new TProfile(Form("c2_%s_%s_mPt_%s",harmonicsname[j].Data(),chargecombisname[i].Data(),centname[k].Data()),Form("c2_%s_%s_mPt_%s",harmonicsname[j].Data(),chargecombisname[i].Data(),centname[k].Data()), fAxes[1]->GetNbins(), fPtbinning);
        fProfileMPt[k][j][i][1] = new TProfile(Form("s2_%s_%s_mPt_%s",harmonicsname[j].Data(),chargecombisname[i].Data(),centname[k].Data()),Form("s2_%s_%s_mPt_%s",harmonicsname[j].Data(),chargecombisname[i].Data(),centname[k].Data()), fAxes[1]->GetNbins(), fPtbinning);

      }}}
  for(Int_t k=0; k<9;k++){
    fNDPhiDEta[k] = AliHistogramManager::CreateHistogram(Form("c2_aa_dPhidEta_%s_entries",centname[k].Data()),Form("c2_aa_dPhidEta_%s",centname[k].Data()), 2, thndphideta);

  }
  //MakeBinMapping();

  TAxis* axes[3];
  for(Int_t ic=0; ic<fNOnePwrtRP; ic++){
    axes[0] = fAxes[fFlow[ic]->GetVarX()];
    axes[1] = fAxes[fFlow[ic]->GetVarY()];
    axes[2] = fAxes[fFlow[ic]->GetVarZ()];
    fFlow[ic]->Init(axes[0],axes[1],axes[2]);
  }
  for(Int_t ic=0; ic<fNTwoPwrtRP; ic++){
    axes[0] = fAxes[fChargeCorrelation[ic]->GetVarX()];
    axes[1] = fAxes[fChargeCorrelation[ic]->GetVarY()];
    axes[2] = fAxes[fChargeCorrelation[ic]->GetVarZ()];
    fChargeCorrelation[ic]->Init(axes[0],axes[1],axes[2]);
  }

  fQvecTPC    = (AliQnCorrectionsQnVector*)  fInputObjects[0];
  fQvecVZEROA = (AliQnCorrectionsQnVector*)  fInputObjects[1];
  fQvecVZEROC = (AliQnCorrectionsQnVector*)  fInputObjects[2];
  fQvecFMDA   = (AliQnCorrectionsQnVector*)  fInputObjects[3];
  fQvecFMDC   = (AliQnCorrectionsQnVector*)  fInputObjects[4];

  fEPcor[0] = new EventPlaneCorrelation(fAxes[0],"VZEROA","VZEROC","TPC");
  fEPcor[0]->SetEventPlanes(fQvecVZEROA,fQvecVZEROC,fQvecTPC);
  fEPcor[1] = new EventPlaneCorrelation(fAxes[0],"FMDA","FMDC","TPC");
  fEPcor[1]->SetEventPlanes(fQvecFMDA,fQvecFMDC,fQvecTPC);

  fQvectors[0] = fQvecVZEROA;
  fQvectors[1] = fQvecVZEROC;
  fQvectors[2] = fQvecFMDA;
  fQvectors[3] = fQvecFMDC;



  fListHistos.Add(fTrackDistribution);

  for(Int_t i=0; i<4; i++){
    fListHistos.Add(fNDPhiDEta[i]);
  }

  for(Int_t i=0; i<4; i++){
    for(Int_t j=0; j<4; j++){
      fListHistos.Add(fProfileCent[j][i][0]);
      fListHistos.Add(fProfileCent[j][i][1]);
      for(Int_t k=0; k<9;k++){
        fListHistos.Add(fProfileMPt[k][j][i][0]);
        fListHistos.Add(fProfileMPt[k][j][i][1]);
      }}}

  //THashList* histos = new THashList();
  //TList histos=fListHistos;
  for(Int_t ic=0; ic<fNOnePwrtRP; ic++){
    for(Int_t iep=0; iep<AliChargeOnePwrtRP::Nqvectors; iep++){
      for(Int_t ich=0; ich<3; ich++){
        for(Int_t icomp=0; icomp<2; icomp++){
          for(Int_t ih=fFlow[ic]->GetMinHarmonic(); ih<=fFlow[ic]->GetMaxHarmonic(); ih++){
            if(fFlow[ic]){
              if(fFlow[ic]->GetTHn(iep, ich, ih, icomp)){
                fListHistos.Add(fFlow[ic]->GetTHn(iep, ich, ih, icomp));
                fListHistos.Add(fFlow[ic]->GetTHn2(iep, ich, ih, icomp));
                fListHistos.Add(fFlow[ic]->GetSumQ(iep, ich, ih, icomp));
                fListHistos.Add(fFlow[ic]->GetSumQ2(iep, ich, ih, icomp));
              }
            }
          }
        }
        if(fFlow[ic]) if(fFlow[ic]->GetTHnMult(iep,ich)) fListHistos.Add(fFlow[ic]->GetTHnMult(iep,ich));
        if(fFlow[ic]) if(fFlow[ic]->GetTHnMult(iep,ich)) fListHistos.Add(fFlow[ic]->GetSumQMult(iep,ich));

      }
    }
  }

  //fListHistos.SetName("ChargeCorrelations");



  for(Int_t icom=0; icom<fNTwoPwrtRP; icom++){
    for(Int_t ic=0; ic<4; ic++){
      for(Int_t ih=0; ih<AliChargeTwoPwrtRP::N2p; ih++){
      //for(Int_t ih=0; ih<1; ih++){
        fListHistos.Add(fChargeCorrelation[icom]->Get2pCorrelationTHn(ic,ih,0));
        fListHistos.Add(fChargeCorrelation[icom]->Get2pCorrelationTHn(ic,ih,1));
      }
      for(Int_t ih=0; ih<AliChargeTwoPwrtRP::N3p; ih++){
        for(Int_t ip=0; ip<AliChargeOnePwrtRP::Nqvectors; ip++) {
          fListHistos.Add(fChargeCorrelation[icom]->Get3pCorrelationTHn(ic,ih,ip,0));
          fListHistos.Add(fChargeCorrelation[icom]->Get3pCorrelationTHn(ic,ih,ip,1));
        }
      }
      fListHistos.Add(fChargeCorrelation[icom]->GetCorrelationTHnMult(ic));
  }}



  fHistosManager->AddToOutputList(&fListHistos);

  //THashList* epcor = new THashList();
  //epcor->SetName("QnCorrelations");
  //epcor->SetOwner(kTRUE);
  //TList* epcor=fListHistos;
  for(Int_t iep=0; iep<AliChargeOnePwrtRP::Nqvectors; iep++){
    if(!fEPcor[iep]) continue;
    for(Int_t icor=0; icor<3; icor++){
      for(Int_t ih=fMinHarmonic; ih<=fMaxHarmonic; ih++){
        for(Int_t icomp=0; icomp<4; icomp++){
          if(!fEPcor[iep]->CorrelationProfile(icor,ih,icomp)) continue;
          fListHistosQn.Add(fEPcor[iep]->CorrelationProfile(icor,ih,icomp));
        }}}
  }


  fHistosManager->AddToOutputList(&fListHistosQn);
}


//________________________________________________________________________________________________________
void AliReducedAnalysisTaskTwoPwrtRP::Process(){
  //
  // Main loop. Called for every event
  //
  Int_t nPosC1=0;
  Int_t nPosC2=0;
  Int_t nPosC3=0;

  //if(fNevents>1000) return;
  if(fNevents>200) return;

  AliReducedEventInfo* event = (AliReducedEventInfo*)fEvent;
  //if(event->CentralityVZERO()<79.0) return;

  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i) {fValues[i]=-9999.;}

  FILL::FillEventInfo(event, fValues);
  fValues[AliReducedVarManager::kNtracksSelected]=0;
  fValues[AliReducedVarManager::kNtracksTotal]=0;

  AliReducedTrackInfo* trackM = 0x0;
  TClonesArray* trackListM = event->GetTracks();
  TIter nextTrackM(trackListM);
  Double_t mTPC=0.,mGlobal=0.;
  while((trackM=static_cast<AliReducedTrackInfo*>(nextTrackM()))) {
    if(!trackM) continue;
    FILL::FillTrackInfo(trackM,fValues);
    if(fTPCtrackcuts->IsSelected(fValues))    mTPC=fValues[AliReducedVarManager::kNtracksTotal]++;
    if(fGlobaltrackcuts->IsSelected(fValues)) mGlobal=fValues[AliReducedVarManager::kNtracksSelected]++;
  }

  fHistosManager->FillHistClass("Event_NoCuts", fValues);

  for(UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
    fHistosManager->FillHistClass("OnlineTriggers_NoCuts", fValues);
  }

  //cout<<"----------"<<endl;
  //cout<<fValues[AliReducedVarManager::kCentVZERO]<<endl;
  //cout<<fValues[AliReducedVarManager::kCentVZEROmTPC]<<endl;
  //cout<<fValues[AliReducedVarManager::kVtxZ]<<endl;
  //cout<<fValues[AliReducedVarManager::kCentQuality]<<endl;

  if(!fEventCuts->IsSelected(fValues)) return;

  //if(!IsEventSelected(event)) return;

  if(fRemoveOutliers){
    // Reduced tree implementation of Mikolaj's cut
    if(mTPC<(1.095*mGlobal-29.32)||mTPC>(1.282*mGlobal+73.60)) return;
    if(mGlobal<(0.780*mTPC-57.43)||mGlobal>(0.913*mTPC+26.77)) return;

    // Remove outliers from centrality VZERO-TPC correlation
    TCutG *cutg = new TCutG("CUTG",15);
    cutg->SetVarX("Centrality(VZERO)-Centrality(TPC)");
    cutg->SetVarY("");
    cutg->SetTitle("Graph");
    cutg->SetFillColor(1);
    cutg->SetPoint(0,0.0,-4.17698);
    cutg->SetPoint(1,13.8925,-6.73886);
    cutg->SetPoint(2,47.6522,-9.91337);
    cutg->SetPoint(3,64.3297,-12.9765);
    cutg->SetPoint(4,75.5019,-17.3762);
    cutg->SetPoint(5,82.3025,-12.0297);
    cutg->SetPoint(6,80.4404,14.0903);
    cutg->SetPoint(7,60.6056,8.8552);
    cutg->SetPoint(8,24.4171,6.18193);
    cutg->SetPoint(9,3.20596,5.01238);
    cutg->SetPoint(10,-0.0323842,1.28094);
    cutg->SetPoint(11,0.0,-4.17698);
    if(!cutg->IsInside(fValues[AliReducedVarManager::kCentVZERO],fValues[AliReducedVarManager::kCentVZEROmTPC])) return;
  }

  fNevents++;

  fHistosManager->FillHistClass("Event_AfterCuts", fValues);

  for(UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
    fHistosManager->FillHistClass("OnlineTriggers_AfterCuts", fValues);
  }


  for(Int_t ip=0; ip<AliChargeOnePwrtRP::Nqvectors; ip++){
      if((fQvectors[ip]->CheckEventPlaneStatus(2,AliQnCorrectionsConstants::kUndefined)||!fQvectors[ip]->CheckEventPlaneStatus(2,AliQnCorrectionsConstants::kRecentering)||fQvectors[ip]->SumOfWeights()<1E-3)){
      fHistosManager->FillHistClass("Event_NoEventPlane", fValues);
      return;
    }
  }


  Float_t Qn[AliChargeOnePwrtRP::Nharmonics][AliChargeOnePwrtRP::Nqvectors][2];
  for(Int_t ip=0; ip<AliChargeOnePwrtRP::Nqvectors; ip++) {
    for(Int_t ih=fMinHarmonic; ih<=fMaxHarmonic; ih++){
      Qn[ih-1][ip][0] = fQvectors[ip]->QxNorm(ih);
      if(ip<2&&ih>3) Qn[ih-1][ip][0]=0.0;  // Set VZERO 4th harmonic x-component to 0
      Qn[ih-1][ip][1] = fQvectors[ip]->QyNorm(ih);
  }}

  Float_t trackMapX[20000][AliChargeOnePwrtRP::Nharmonics] = {{-1.}};
  Float_t trackMapY[20000][AliChargeOnePwrtRP::Nharmonics] = {{-1.}};
  Float_t trackMapPhi[20000] = {-1.};
  Float_t trackMapPtEtaCharge[20000][3] = {{-1.}};
  ULong_t  trackMapFlag[20000] = {0};


  for(Int_t ic=0; ic<fNOnePwrtRP; ic++) fFlow[ic]->Clear();
  for(Int_t ic=0; ic<fNTwoPwrtRP; ic++) fChargeCorrelation[ic]->Clear();




    //// Track efficiency //////////////

  Double_t PbeffDCA[48];
  PbeffDCA[0] = 0.89596; PbeffDCA[1] = 0.92504; PbeffDCA[2] = 0.944326; PbeffDCA[3] = 0.948261; PbeffDCA[4] = 0.948522;
  PbeffDCA[5] = 0.946941; PbeffDCA[6] = 0.948537; PbeffDCA[7] = 0.95194; PbeffDCA[8] = 0.95228; PbeffDCA[9] = 0.954775;
  PbeffDCA[10] = 0.956055; PbeffDCA[11] = 0.955758; PbeffDCA[12] = 0.951268; PbeffDCA[13] = 0.952081; PbeffDCA[14] = 0.947992;
  PbeffDCA[15] = 0.943755; PbeffDCA[16] = 0.938106; PbeffDCA[17] = 0.932583; PbeffDCA[18] = 0.927206; PbeffDCA[19] = 0.921112;
  PbeffDCA[20] = 0.913127; PbeffDCA[21] = 0.910508; PbeffDCA[22] = 0.905838; PbeffDCA[23] = 0.90243; PbeffDCA[24] = 0.897893;
  PbeffDCA[25] = 0.895459; PbeffDCA[26] = 0.89308; PbeffDCA[27] = 0.895013; PbeffDCA[28] = 0.890872; PbeffDCA[29] = 0.892046;
  PbeffDCA[30] = 0.890133; PbeffDCA[31] = 0.894786; PbeffDCA[32] = 0.891886; PbeffDCA[33] = 0.893649; PbeffDCA[34] = 0.891418;
  PbeffDCA[35] = 0.897408; PbeffDCA[36] = 0.890273; PbeffDCA[37] = 0.896526; PbeffDCA[38] = 0.892236; PbeffDCA[39] = 0.895075;
  PbeffDCA[40] = 0.89516; PbeffDCA[41] = 0.889562; PbeffDCA[42] = 0.890777; PbeffDCA[43] = 0.892925; PbeffDCA[44] = 0.890358;
  PbeffDCA[45] = 0.897014; PbeffDCA[46] = 0.886445; PbeffDCA[47] = 0.890928;

  Double_t  fe = PbeffDCA[46];
  for(int i=0;i<48;i++) PbeffDCA[i]=(PbeffDCA[i]-fe)/PbeffDCA[i];

  TAxis* pbins = new TAxis(48,0.2,5.0);

  Double_t w[16][10] = {{0.970977,0.99772,0.99419,0.978897,0.971534,0.97598,0.989627,0.998296,1,0.933612},
{0.954112,0.983011,0.982031,0.968427,0.965237,0.967845,0.975416,0.974368,1,0.935514},
{0.957119,0.987402,0.987601,0.976029,0.969505,0.970452,0.987786,0.997328,1,0.920526},
{0.939055,0.969236,0.97094,0.960891,0.956266,0.961443,0.962705,0.987657,1,0.928272},
{0.940339,0.972623,0.973738,0.961084,0.961272,0.967403,0.969807,1,0.989227,0.920163},
{0.948144,0.98088,0.982423,0.97495,0.970065,0.982265,0.991422,1,0.986479,0.983271},
{0.92695,0.960166,0.961328,0.954786,0.951977,0.958845,0.968716,0.975792,1,0.914216},
{0.93233,0.968171,0.970805,0.962505,0.966148,0.969439,0.982002,0.980514,1,0.892224},
{0.932492,0.9683,0.972112,0.967793,0.955948,0.960787,0.990176,0.992178,1,0.927829},
{0.940215,0.975212,0.979227,0.979416,0.976191,0.995886,0.989413,0.987496,1,0.964744},
{0.93745,0.976152,0.976806,0.97805,0.975628,0.983384,1,0.998986,0.999524,0.972516},
{0.920143,0.958541,0.960134,0.953553,0.966066,0.982041,1,0.976215,0.997212,0.850735},
{0.897359,0.93182,0.93689,0.9368,0.941417,0.947915,0.967267,1,0.997503,0.948772},
{0.924567,0.957018,0.963272,0.953267,0.96477,0.962868,0.974634,0.997557,0.983571,1},
{0.884835,0.91816,0.92338,0.932512,0.923256,0.92942,0.941799,0.88556,1,0.777778},
{0.82124,0.859637,0.85144,0.864,0.836804,0.896039,0.842593,0.942029,1,0.525253}
};

  TAxis* ax_cent = new TAxis( 20,0,100);
  TAxis* ax_pt   = new TAxis( 10,0,5);

  Int_t t_centbin=0;
  Int_t t_ptbin=0;

  t_centbin=ax_cent->FindBin(fValues[AliReducedVarManager::kCentVZERO])-1;


  TRandom3* rndm = new TRandom3();

  const Int_t minOneFlag=20;
  Int_t minTwoFlag=minOneFlag+fNOnePwrtRP;
  Int_t diffTwoFlag=minOneFlag+fNOnePwrtRP +fNTwoPwrtRP;
  Float_t dphi;
  Float_t pi = TMath::Pi();
  Double_t dphideta[2];
  Float_t ceta2;
  AliReducedTrackInfo* track = 0x0;
  TClonesArray* trackList = event->GetTracks();
  trackList->Randomize(2);
  TIter nextTrack(trackList);
  Bool_t keepTrack=kFALSE;
  Int_t itrack=-1,itrack2=-1;
  AliChargeOnePwrtRP* correlation=0x0;
  Double_t trackDist[4];
  AliReducedTrackInfo* track2 = 0x0;
  Int_t ccentbin = fProfileCent[0][0][0]->GetXaxis()->FindBin(event->CentralityVZERO())-1;
  Double_t cphi,cphi2;
  //cout<<"peek 4 "<<fChargeCorrelation[0]->Get2pCorrelationTHn(0,1,0)->GetName()<<endl;
  for(Int_t citrack=0; citrack<trackList->GetEntriesFast(); citrack++){
    track = (AliReducedTrackInfo*) trackList->At(citrack);

  //while((track=static_cast<AliReducedTrackInfo*>(nextTrack()))) {
    //if(!track) continue;
    itrack++;
    itrack2++;
    //if(itrack>500) continue;
    keepTrack=kFALSE;

    FILL::FillTrackInfo(track,fValues);
    Float_t cphi=(fTracksTPCstandalone ? track->PhiTPC() : track->Phi() );
    Float_t ceta=(fTracksTPCstandalone ? track->EtaTPC() : track->Eta() );
    if(fTrackCutsOnePwrtRP[0]->IsSelected(fValues)) {
    fHistosManager->FillHistClass("TrackQA_Profile_TPC_TOF_ITS",fValues);
    //TIter nextTrack2(trackList);
    //while((track2=static_cast<AliReducedTrackInfo*>(nextTrack2()))) {
      //if(!track2) continue;
    for(Int_t citrack2=citrack+1; citrack2<trackList->GetEntriesFast(); citrack2++){
      track2 = (AliReducedTrackInfo*) trackList->At(citrack2);
      if(track->TrackId()==track2->TrackId()) continue;
      FILL::FillTrackInfo(track2,fValues);
      if(!fTrackCutsOnePwrtRP[0]->IsSelected(fValues)) continue;
      cphi2=(fTracksTPCstandalone ? track2->PhiTPC() : track2->Phi() );
      ceta2=(fTracksTPCstandalone ? track2->EtaTPC() : track2->Eta() );
      dphi=cphi-cphi2;
      if(dphi<-pi) dphi+=2*pi;
      if(dphi>pi)  dphi-=2*pi;

      dphideta[0] = dphi;
      dphideta[1] = TMath::Abs(ceta-ceta2);

      fNDPhiDEta[ccentbin]->Fill(dphideta);

      Float_t cmpt = (track->Pt()+track2->Pt())/2.;
        if(track2->Charge()==1&&track->Charge()==1){
          //nPosC1++;
          for(Int_t ih=1; ih<5; ih++){
            //cout<<TMath::Abs(dphi)<<"  "<<pi/((Double_t) ih)<<endl;
            if(TMath::Abs(dphi)<pi/((Double_t) ih)) continue;
            fProfileCent[ih-1][0][0]->Fill(event->CentralityVZERO(), TMath::Cos(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][0][0]->Fill(cmpt, TMath::Cos(ih*(cphi2-cphi)));
            fProfileCent[ih-1][0][1]->Fill(event->CentralityVZERO(), TMath::Sin(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][0][1]->Fill(cmpt, TMath::Sin(ih*(cphi2-cphi)));
          }
        }
        else if(track2->Charge()==1&&track->Charge()==-1){
          for(Int_t ih=1; ih<5; ih++){
            if(TMath::Abs(dphi)<pi/((Double_t) ih)) continue;
            fProfileCent[ih-1][1][0]->Fill(event->CentralityVZERO(), TMath::Cos(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][1][0]->Fill(cmpt, TMath::Cos(ih*(cphi2-cphi)));
            fProfileCent[ih-1][1][1]->Fill(event->CentralityVZERO(), TMath::Sin(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][1][1]->Fill(cmpt, TMath::Sin(ih*(cphi2-cphi)));
          }
        }
        else if(track2->Charge()==-1&&track->Charge()==1){
          for(Int_t ih=1; ih<5; ih++){
            if(TMath::Abs(dphi)<pi/((Double_t) ih)) continue;
            fProfileCent[ih-1][2][0]->Fill(event->CentralityVZERO(), TMath::Cos(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][2][0]->Fill(cmpt, TMath::Cos(ih*(cphi2-cphi)));
            fProfileCent[ih-1][2][1]->Fill(event->CentralityVZERO(), TMath::Sin(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][2][1]->Fill(cmpt, TMath::Sin(ih*(cphi2-cphi)));
          }
        }
        else if(track2->Charge()==-1&&track->Charge()==-1){
          for(Int_t ih=1; ih<5; ih++){
            if(TMath::Abs(dphi)<pi/((Double_t) ih)) continue;
            fProfileCent[ih-1][3][0]->Fill(event->CentralityVZERO(), TMath::Cos(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][3][0]->Fill(cmpt, TMath::Cos(ih*(cphi2-cphi)));
            fProfileCent[ih-1][3][1]->Fill(event->CentralityVZERO(), TMath::Sin(ih*(cphi2-cphi)));
            fProfileMPt[ccentbin][ih-1][3][1]->Fill(cmpt, TMath::Sin(ih*(cphi2-cphi)));
          }
        }



    }};


    FILL::FillTrackInfo(track,fValues);



    t_ptbin=ax_pt->FindBin(fValues[AliReducedVarManager::kPtTPC])-1;

    Double_t weight = w[t_centbin][t_ptbin];

    trackDist[0] = fValues[AliReducedVarManager::kCentVZERO];
    trackDist[1] = track->Pt();
    trackDist[2] = track->Eta();
    trackDist[3] = track->Phi();

    if(fTrackCutsOnePwrtRP[0]->IsSelected(fValues)) fTrackDistribution->Fill(trackDist);

    for(Int_t flag=minOneFlag; flag<(diffTwoFlag+fNTwoPwrtRP); flag++)  track->UnsetFlag(flag);  // these bits will be used as track+pid selection flags

    //if(rndm->Rndm()>weight) {itrack2--; continue;}

    Int_t pb = pbins->FindBin(track->PtTPC())-1;
    if(rndm->Rndm()>PbeffDCA[pb]) fValues[AliReducedVarManager::kEMCALmatchedEnergy]=1.0;
    else fValues[AliReducedVarManager::kEMCALmatchedEnergy]=0.0;


    for(Int_t ic=0; ic<fNOnePwrtRP; ic++){
      if(fTrackCutsOnePwrtRP[ic]->IsSelected(fValues)) {
        track->SetFlag(minOneFlag+ic);
        fHistosManager->FillHistClass(Form("TrackQA_%s_TPC_TOF_ITS",fTrackCutsOnePwrtRP[ic]->GetName()),fValues);
        keepTrack=kTRUE;
        //nPosC3++;
      }
    }


    for(Int_t ic=0; ic<fNTwoPwrtRP; ic++){
      if(fTrackCutsTwoPwrtRP[ic][0]->IsSelected(fValues)) {
        track->SetFlag(minTwoFlag+ic);
        keepTrack=kTRUE;
        if(track->Charge()==1) nPosC3++;
      }
      if(fTrackCutsTwoPwrtRP[ic][1]->IsSelected(fValues)) {
        track->SetFlag(diffTwoFlag+ic);
        keepTrack=kTRUE;
      }
    }




    if(keepTrack){


      Float_t phit=(fTracksTPCstandalone ? track->PhiTPC() : track->Phi() );

      // store some values, will save time later
      trackMapPtEtaCharge[itrack2][0]=(fTracksTPCstandalone ? track->PtTPC() : track->Pt() );
      trackMapPtEtaCharge[itrack2][1]=(fTracksTPCstandalone ? track->EtaTPC(): track->Eta() );
      trackMapPtEtaCharge[itrack2][2]=track->Charge();
      trackMapPhi[itrack2]=phit;
      trackMapFlag[itrack2]=track->GetQualityFlag();

      for(Int_t ih=1; ih<=fMaxHarmonic; ih++){
        trackMapX[itrack2][ih-1]=TMath::Cos(ih*phit);
        trackMapY[itrack2][ih-1]=TMath::Sin(ih*phit);
      }
    }
    else {itrack2--;}
  }


  //cout<<"single "<<nPosC1<<endl;
  //return;

  Int_t nTracks=itrack2+1;


  if(nTracks<2) return;

  Int_t CtBin=0;
  Float_t cent = fValues[AliReducedVarManager::kCentVZERO];
  while(cent>fCtbinning[CtBin]) {CtBin++;}
  CtBin--;

  for(Int_t iep=0; iep<AliChargeOnePwrtRP::Nqvectors; iep++) if(fEPcor[iep]) fEPcor[iep]->FillCorrelations(cent);

  Int_t n,m;
  Int_t ncor;
  Int_t dims[7]={0};
  Int_t dimsMax[7] = {1, fNPtbins,fNMPtbins,fNDPtbins,fNEtabins,fNMEtabins,fNDEtabins}; // number of bins for event variables should be set to 1



  //Float_t somevalues[6]; // SetVars(i);   // 0: Ct, 1: 1Pt, 2: mPt, 3: dPt, 4: 1Eta 5: mEta,  6: dEta
  Float_t pt1, eta1, phi1, charge1;
  Float_t pt2, eta2, phi2, charge2;
  Float_t mpt, dpt, meta, deta, charge;
  ULong_t f1, f2;
  Int_t Pt1Bin=0, Eta1Bin=0;
  Int_t Pt2Bin=0, Eta2Bin=0;
  Int_t mPtBin=0, dPtBin=0, mEtaBin=0, dEtaBin=0;
  //Float_t addCors[AliChargeTwoPwrtRP::N2p+AliChargeTwoPwrtRP::N3p][2];
  //Float_t* addCors[] = AliChargeTwoPwrtRP::fTwoPCorrelationArray;
  Float_t cor[AliChargeTwoPwrtRP::N2p][4];
  Float_t qncor[AliChargeOnePwrtRP::Nqvectors][4];
  Int_t it,it2;
  for(it=0; it<nTracks; it++){
    pt1    =trackMapPtEtaCharge[it][0];
    eta1   =trackMapPtEtaCharge[it][1];
    charge1=trackMapPtEtaCharge[it][2];
    f1     =trackMapFlag[it];


    Pt1Bin=0;
    while(pt1>fPtbinning[Pt1Bin]) Pt1Bin++;
    Pt1Bin--;
    Eta1Bin=0;
    while(eta1>fEtabinning[Eta1Bin]) Eta1Bin++;
    Eta1Bin--;

    for(Int_t ic=0; ic<fNOnePwrtRP; ic++){
      fFlow[ic]->SetBin(0, Pt1Bin,0,0,Eta1Bin,0,0);
      //cout<<CtBin<<"  "<<Pt1Bin<<"  "<<Eta1Bin<<"  "<<fFlow[ic]->Bin()<<"  "<<trackMapX[it][0]<<"  "<<trackMapX[it][1]<<endl;
      fFlow[ic]->AddTrack( fFlow[ic]->Charge(charge1), trackMapX[it], trackMapY[it]);
    }


    for(it2=0; it2<nTracks; it2++){
      if(it==it2) continue;
      pt2    =trackMapPtEtaCharge[it2][0];
      eta2   =trackMapPtEtaCharge[it2][1];
      charge2=trackMapPtEtaCharge[it2][2];
      f2     =trackMapFlag[it2];
      mpt=(pt1+pt2)/2.;
      dpt=pt1-pt2+1E-6;
      if(dpt<0.) dpt*=-1.;
      meta=(eta1+eta2)/2.;
      deta=eta1-eta2+1E-6;
      if(deta<0.) deta*=-1.;

      dphi=trackMapPhi[it2]-trackMapPhi[it];
      if(dphi<-pi) dphi+=2*pi;
      if(dphi>pi)  dphi-=2*pi;
      if(TMath::Abs(dphi)<pi/3.) continue;

      mPtBin=0;
      while(mpt>fMPtbinning[mPtBin]) mPtBin++;
      mPtBin--;
      dPtBin=0;
      while(dpt>fDPtbinning[dPtBin]) dPtBin++;
      dPtBin--;
      mEtaBin=0;
      while(meta>fMEtabinning[mEtaBin]) mEtaBin++;
      mEtaBin--;
      dEtaBin=0;
      while(deta>fDEtabinning[dEtaBin]) dEtaBin++;
      dEtaBin--;


      cor[0][0] = trackMapX[it][0]*trackMapX[it2][0]; // cos(phi1)*cos(phi2)
      cor[0][1] = trackMapX[it][0]*trackMapY[it2][0]; // cos(phi1)*sin(phi2)
      cor[0][2] = trackMapY[it][0]*trackMapX[it2][0]; // sin(phi1)*cos(phi2)
      cor[0][3] = trackMapY[it][0]*trackMapY[it2][0]; // sin(phi1)*sin(phi2)

      cor[1][0] = trackMapX[it][1]*trackMapX[it2][1]; // cos(2*phi1)*cos(2*phi2)
      cor[1][1] = trackMapX[it][1]*trackMapY[it2][1]; // cos(2*phi1)*sin(2*phi2)
      cor[1][2] = trackMapY[it][1]*trackMapX[it2][1]; // sin(2*phi1)*cos(2*phi2)
      cor[1][3] = trackMapY[it][1]*trackMapY[it2][1]; // sin(2*phi1)*sin(2*phi2)

      cor[2][0] = trackMapX[it][2]*trackMapX[it2][2]; // cos(3*phi1)*cos(3*phi2)
      cor[2][1] = trackMapX[it][2]*trackMapY[it2][2]; // cos(3*phi1)*sin(3*phi2)
      cor[2][2] = trackMapY[it][2]*trackMapX[it2][2]; // sin(3*phi1)*cos(3*phi2)
      cor[2][3] = trackMapY[it][2]*trackMapY[it2][2]; // sin(3*phi1)*sin(3*phi2)

      cor[3][0] = trackMapX[it][3]*trackMapX[it2][3]; // cos(4*phi1)*cos(4*phi2)
      cor[3][1] = trackMapX[it][3]*trackMapY[it2][3]; // cos(4*phi1)*sin(4*phi2)
      cor[3][2] = trackMapY[it][3]*trackMapX[it2][3]; // sin(4*phi1)*cos(4*phi2)
      cor[3][3] = trackMapY[it][3]*trackMapY[it2][3]; // sin(4*phi1)*sin(4*phi2)




    //for(Int_t i=0; i<AliChargeTwoPwrtRP::N2p; i++){
    for(Int_t i=0; i<4; i++){
      AliChargeTwoPwrtRP::fPair[i][0]= cor[i][0]+cor[i][3];  // cos(n(phi1-phi2))
      AliChargeTwoPwrtRP::fPair[i][1]=(cor[i][0]+cor[i][3])*(cor[i][0]+cor[i][3]);
      AliChargeTwoPwrtRP::fPair[i+4][0]= cor[i][2]-cor[i][1];  // sin(n(phi1-phi2))
      AliChargeTwoPwrtRP::fPair[i+4][1]=(cor[i][2]-cor[i][1])*(cor[i][2]-cor[i][1]);
    }
    for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++){

        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+j][0] =  ( cor[0][0]-cor[0][3])*Qn[1][j][0] +  ( cor[0][1]+cor[0][2])*Qn[1][j][1];  // cos(phi1+phi2-2Psi2A);
        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+j][1] = (( cor[0][0]-cor[0][3])*Qn[1][j][0] +  ( cor[0][1]+cor[0][2])*Qn[1][j][1])
                                                                * (( cor[0][0]-cor[0][3])*Qn[1][j][0] +  ( cor[0][1]+cor[0][2])*Qn[1][j][1]);

        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+AliChargeOnePwrtRP::Nqvectors+j][0] =  ( cor[1][0]-cor[1][3])*Qn[3][j][0] +  ( cor[1][1]+cor[1][2])*Qn[3][j][1]; // cos(2(phi1+phi2-2Psi4A))
        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+AliChargeOnePwrtRP::Nqvectors+j][1] = (( cor[1][0]-cor[1][3])*Qn[3][j][0] +  ( cor[1][1]+cor[1][2])*Qn[3][j][1])
                                                                                              * (( cor[1][0]-cor[1][3])*Qn[3][j][0] +  ( cor[1][1]+cor[1][2])*Qn[3][j][1]);

        qncor[j][0] = trackMapX[it2][1]*Qn[1][j][0]; // cos(phi2)*cos(2*EP2)
        qncor[j][1] = trackMapX[it2][1]*Qn[1][j][1]; // cos(phi2)*sin(2*EP2)
        qncor[j][2] = trackMapY[it2][1]*Qn[1][j][0]; // sin(phi2)*cos(2*EP2)
        qncor[j][3] = trackMapY[it2][1]*Qn[1][j][1]; // sin(phi2)*sin(2*EP2)

        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+2*AliChargeOnePwrtRP::Nqvectors+j][0] =     ( cor[0][0]+cor[0][3])*(qncor[j][0]+qncor[j][3]) +  ( cor[0][1]-cor[0][2])*(qncor[j][1]-qncor[j][2]); // cos(phi1-3phi2+2Psi2A)
        AliChargeTwoPwrtRP::fPair[AliChargeTwoPwrtRP::N2p+2*AliChargeOnePwrtRP::Nqvectors+j][1] =    (( cor[0][0]+cor[0][3])*(qncor[j][0]+qncor[j][3]) +  ( cor[0][1]-cor[0][2])*(qncor[j][1]-qncor[j][2]))
                                                                                                *    (( cor[0][0]+cor[0][3])*(qncor[j][0]+qncor[j][3]) +  ( cor[0][1]-cor[0][2])*(qncor[j][1]-qncor[j][2])) ;
    }


      if(charge1==1&&charge2==1) nPosC2++;



      AliChargeTwoPwrtRP::SetCharge(charge1,charge2);
      AliChargeTwoPwrtRP::SetBin(0,Pt1Bin,mPtBin,dPtBin,Eta1Bin,mEtaBin,dEtaBin);

      //if(fChargeCorrelation[0]->Bin()>104856) cout<<Pt1Bin<<"  "<<mPtBin<<"  "<<dPtBin<<"  "<<Eta1Bin<<"  "<<mEtaBin<<"  "<<dEtaBin<<"  "<<pt1<<"  "<<pt2<<"  "<<dpt<<"  "<<eta1<<"  "<<eta2<<"  "<<it<<"  "<<it2<<"  "<<nTracks<<endl;

      for(Int_t ic=0; ic<fNTwoPwrtRP; ic++){
        if((trackMapFlag[it]&(ULong_t(1)<<((UShort_t)(minTwoFlag+ ic))))
          &&(trackMapFlag[it2]&(ULong_t(1)<<((UShort_t)(diffTwoFlag+ ic))))){
            fChargeCorrelation[ic]->GetTwoParticleCorrelation();
        }
      }




    }
  };


  dims[0] = CtBin;
  dimsMax[0] = fNCtbins;


  AliChargeOnePwrtRP* oneP;
  for(Int_t ic=0; ic<fNOnePwrtRP; ic++){
    oneP=fFlow[ic];

    for(Int_t ibinx=0; ibinx<oneP->GetNbinsX(); ibinx++){
      for(Int_t ibiny=0; ibiny<oneP->GetNbinsY(); ibiny++){
        for(Int_t ibinz=0; ibinz<oneP->GetNbinsZ(); ibinz++){

          Int_t xbin=(oneP->GetVarX()==-1 ? 0 :dims[oneP->GetVarX()]+1);
          Int_t ybin=(oneP->GetVarY()==-1 ? 0 :dims[oneP->GetVarY()]+1);
          Int_t zbin=(oneP->GetVarZ()==-1 ? 0 :dims[oneP->GetVarZ()]+1);

          if(oneP->GetTrackVarMap(0)==0) xbin=ibinx+1;
          if(oneP->GetTrackVarMap(0)==1) ybin=ibinx+1;
          if(oneP->GetTrackVarMap(0)==2) zbin=ibinx+1;
          if(oneP->GetTrackVarMap(1)==1) ybin=ibiny+1;
          if(oneP->GetTrackVarMap(1)==2) zbin=ibiny+1;
          if(oneP->GetTrackVarMap(2)==2) zbin=ibinz+1;

          Int_t nbinx = (oneP->GetVarX()==-1 ? 1 : dimsMax[oneP->GetVarX()]+2);
          Int_t nbiny = (oneP->GetVarY()==-1 ? 1 : dimsMax[oneP->GetVarY()]+2);
          Int_t nbinz = (oneP->GetVarZ()==-1 ? 1 : dimsMax[oneP->GetVarZ()]+2);


          for(Int_t ip=0; ip<AliChargeOnePwrtRP::Nqvectors; ip++)  oneP->FillCorrelation( ibinx*oneP->GetNbinsY()*oneP->GetNbinsZ()+ibiny*oneP->GetNbinsZ()+ibinz, xbin*nbiny*nbinz+ybin*nbinz+zbin, fQvectors[ip], ip);
        }
      }
    }
  }

  AliChargeTwoPwrtRP* twoP;
  for(Int_t ic=0; ic<fNTwoPwrtRP; ic++){
    twoP=fChargeCorrelation[ic];

    for(Int_t ibinx=0; ibinx<twoP->GetNbinsX(); ibinx++){
      for(Int_t ibiny=0; ibiny<twoP->GetNbinsY(); ibiny++){
        for(Int_t ibinz=0; ibinz<twoP->GetNbinsZ(); ibinz++){

          Int_t xbin=(twoP->GetVarX()==-1 ? 0 :dims[twoP->GetVarX()]+1);
          Int_t ybin=(twoP->GetVarY()==-1 ? 0 :dims[twoP->GetVarY()]+1);
          Int_t zbin=(twoP->GetVarZ()==-1 ? 0 :dims[twoP->GetVarZ()]+1);

          if(twoP->GetTrackVarMap(0)==0) xbin=ibinx+1;
          if(twoP->GetTrackVarMap(0)==1) ybin=ibinx+1;
          if(twoP->GetTrackVarMap(0)==2) zbin=ibinx+1;
          if(twoP->GetTrackVarMap(1)==1) ybin=ibiny+1;
          if(twoP->GetTrackVarMap(1)==2) zbin=ibiny+1;
          if(twoP->GetTrackVarMap(2)==2) zbin=ibinz+1;

          Int_t nbinx = (twoP->GetVarX()==-1 ? 1 : dimsMax[twoP->GetVarX()]+2);
          Int_t nbiny = (twoP->GetVarY()==-1 ? 1 : dimsMax[twoP->GetVarY()]+2);
          Int_t nbinz = (twoP->GetVarZ()==-1 ? 1 : dimsMax[twoP->GetVarZ()]+2);

          //cout<<"bin "<<ibinx<<"  "<<ibiny<<"  "<<ibinz<<"  "<<xbin<<"  "<<ybin<<"  "<<zbin<<"  "<<nbinx<<"  "<<nbiny<<"  "<<nbinz<<"  "<<ibinx*fChargeCorrelation[ic]->GetNbinsY()*fChargeCorrelation[ic]->GetNbinsZ()+ibiny*fChargeCorrelation[ic]->GetNbinsZ()+ibinz<<endl;
          twoP->FillCorrelation(xbin*nbiny*nbinz+ybin*nbinz+zbin, ibinx*fChargeCorrelation[ic]->GetNbinsY()*fChargeCorrelation[ic]->GetNbinsZ()+ibiny*fChargeCorrelation[ic]->GetNbinsZ()+ibinz);

        }
      }

    }

  }

  return;







  }



  //__________________________________________________________________
  void AliReducedAnalysisTaskTwoPwrtRP::Finish()
  {
    //
    // Finish Task
    //

  }


//__________________________________________________________________
void AliReducedAnalysisTaskTwoPwrtRP::MakeBinMapping()
{
  // using bin maps turned out to be much slower than using while loops

  for(Int_t iax=1; iax<7; iax++){
    Float_t max = fAxes[iax]->GetBinUpEdge(fAxes[iax]->GetNbins());
    Float_t min = fAxes[iax]->GetBinUpEdge(0);
    Int_t nbin= (Int_t) ((max-min)*100+10E-6);
    min+=0.005;
    for(Int_t i=0; i< nbin; i++) {
      Int_t bin=fAxes[iax]->FindBin(min)-1;
      if(iax==1) fPtBinningMap[min*100+10E-6]=bin;
      if(iax==2) fMPtBinningMap[min*100+10E-6]=bin;
      if(iax==3) fDPtBinningMap[min*100+10E-6]=bin;
      if(iax==4) fEtaBinningMap[min*100+10E-6]=bin;
      if(iax==5) fMEtaBinningMap[min*100+10E-6]=bin;
      if(iax==6) fDEtaBinningMap[min*100+10E-6]=bin;
      min+=0.01;
    }
  }
}
