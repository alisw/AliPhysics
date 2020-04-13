/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//new
//-----------------------------------------------------------------
//                 AliAnalysisTaskHelium3Pi class
//-----------------------------------------------------------------

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TChain.h>
#include <Riostream.h>
#include <TTree.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TString.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliMultSelection.h"
//#include <AliVTrack.h"
#include "AliAnalysisTaskHelium3Pi.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHelium3Pi)


//________________________________________________________________________

AliAnalysisTaskHelium3Pi::AliAnalysisTaskHelium3Pi() 
: AliAnalysisTaskSE(),
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fDataType("PbPb"),
  fYear(2011),
  fVzmax(10),
  fApplyFlatten(kTRUE),
  fFill3Hetree(kTRUE),
  fDoFlow(kTRUE),
  fESDevent(0),                         //! 
  fevent(0),   
  fListHist(0), 
  fHistEventMultiplicity(0),
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCent(0),
  fHistTrackMultiplicitySemiCent(0),
  fHistTrackMultiplicityMB(0),
  fHistTrackMultiplicityINT7(0),
  fHistTrackMultiplicityPVCent(0),
  fHistTrackMultiplicityPVSemiCent(0),
  fHistTrackMultiplicityPVMB(0),
  fHistTrackMultiplicityPVINT7(0),
  fhBB(0),
  fhTOF(0),
  fhMassTOF(0),
  fhBBPions(0),
  fhBBHe(0),
  hQVzAQVzCvsCentrality(0),
  hqEPCvsCentrality(0), 
  hqEPAvsCentrality(0),
  hqEPvsCentrality(0),
  fNtuple1(0),
  teventtype(0),
  tTrackNumber(0),
  tpercentile(0),
  txPrimaryVertex(0),
  tyPrimaryVertex(0),
  tzPrimaryVertex(0),
  txSecondaryVertex(0),
  tySecondaryVertex(0),
  tzSecondaryVertex(0),
  tdcaTracks(0),
  tCosPointingAngle(0),
  tDCAV0toPrimaryVertex(0),
  tHeSign(0),
  tHepInTPC(0),
  tHeTPCsignal(0),
  tDcaHeToPrimVertex(0),
  tHeEta(0),
  tmomHex(0),
  tmomHey(0),
  tmomHez(0),
  tmomHeAtSVx(0),
  tmomHeAtSVy(0),
  tmomHeAtSVz(0),
  tHeTPCNcls(0),
  tHeimpactXY(0),
  tHeimpactZ(0),
  tHeITSClusterMap(0),
  tIsHeITSRefit(0),
  tPionSign(0),
  tPionpInTPC(0),
  tPionTPCsignal(0),
  tDcaPionToPrimVertex(0),
  tPionEta(0),
  tmomPionx(0),
  tmomPiony(0),
  tmomPionz(0),
  tmomNegPionAtSVx(0),
  tmomNegPionAtSVy(0),
  tmomNegPionAtSVz(0),
  tPionTPCNcls(0),
  tPionimpactXY(0),
  tPionimpactZ(0),
  tPionITSClusterMap(0),
  tIsPiITSRefit(0),
  txn(0),
  txp(0),
  tuqV0A(0),
  tuqV0C(0),
  fNtuple4(0),
  tHeleventtype(0),
  tHelpercentile(0),
  tHelSign(0),
  tHelpinTPC(0),
  tHelGetTPCsignal(0),
  tHelPx(0),
  tHelPy(0),
  tHelPz(0),
  tHelEta(0),
  tHelisTOF(0),
  tHelTOFpull(0),
  tHeMass(0),
  tHelimpactXY(0),
  tHelimpactZ(0),
  tHelmapITS(0),
  tHelBetaTOF(0),
  tHelIsITSrefit(0),
  fESDtrackCuts(0),
  fPIDResponse(0)
  
{
  printf("Dummy Constructor");
  
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireITSStandAlone(kFALSE);
  fESDtrackCuts->SetRequireITSPureStandAlone(kFALSE);
      
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(60);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);

}

//________________________________________________________________________
AliAnalysisTaskHelium3Pi::AliAnalysisTaskHelium3Pi(TString name) 
  : AliAnalysisTaskSE(name),
    fAnalysisType("ESD"), 
    fCollidingSystems(0), 
    fDataType("PbPb"),
    fYear(2011),
    fVzmax(10),
    fApplyFlatten(kTRUE),
    fFill3Hetree(kTRUE),
    fDoFlow(kTRUE),
    fESDevent(0),                         //! 
    fevent(0),    
    fListHist(0), 
    fHistEventMultiplicity(0),
    fHistTrackMultiplicity(0),
    fHistTrackMultiplicityCent(0),
    fHistTrackMultiplicitySemiCent(0),
    fHistTrackMultiplicityMB(0),
    fHistTrackMultiplicityINT7(0),
    fHistTrackMultiplicityPVCent(0),
    fHistTrackMultiplicityPVSemiCent(0),
    fHistTrackMultiplicityPVMB(0),
    fHistTrackMultiplicityPVINT7(0),
    fhBB(0),
    fhTOF(0),
    fhMassTOF(0),
    fhBBPions(0),
    fhBBHe(0),
    hQVzAQVzCvsCentrality(0),
    hqEPCvsCentrality(0), 
    hqEPAvsCentrality(0),
    hqEPvsCentrality(0),
    fNtuple1(0),
    teventtype(0),
    tTrackNumber(0),
    tpercentile(0),
    txPrimaryVertex(0),
    tyPrimaryVertex(0),
    tzPrimaryVertex(0),
    txSecondaryVertex(0),
    tySecondaryVertex(0),
    tzSecondaryVertex(0),
    tdcaTracks(0),
    tCosPointingAngle(0),
    tDCAV0toPrimaryVertex(0),
    tHeSign(0),
    tHepInTPC(0),
    tHeTPCsignal(0),
    tDcaHeToPrimVertex(0),
    tHeEta(0),
    tmomHex(0),
    tmomHey(0),
    tmomHez(0),
    tmomHeAtSVx(0),
    tmomHeAtSVy(0),
    tmomHeAtSVz(0),
    tHeTPCNcls(0),
    tHeimpactXY(0),
    tHeimpactZ(0),
    tHeITSClusterMap(0),
    tIsHeITSRefit(0),
    tPionSign(0),
    tPionpInTPC(0),
    tPionTPCsignal(0),
    tDcaPionToPrimVertex(0),
    tPionEta(0),
    tmomPionx(0),
    tmomPiony(0),
    tmomPionz(0),
    tmomNegPionAtSVx(0),
    tmomNegPionAtSVy(0),
    tmomNegPionAtSVz(0),
    tPionTPCNcls(0),
    tPionimpactXY(0),
    tPionimpactZ(0),
    tPionITSClusterMap(0),
    tIsPiITSRefit(0),
    txn(0),
    txp(0),
    tuqV0A(0),
    tuqV0C(0),
    fNtuple4(0),
    tHeleventtype(0),
    tHelpercentile(0),
    tHelSign(0),
    tHelpinTPC(0),
    tHelGetTPCsignal(0),
    tHelPx(0),
    tHelPy(0),
    tHelPz(0),
    tHelEta(0),
    tHelisTOF(0),
    tHelTOFpull(0),
    tHeMass(0),
    tHelimpactXY(0),
    tHelimpactZ(0),
    tHelmapITS(0),
    tHelBetaTOF(0),
    tHelIsITSrefit(0),
    fESDtrackCuts(0),
    fPIDResponse(0)
{					  
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container ()
 
  printf("Real Constructor");
 
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireITSStandAlone(kFALSE);
  fESDtrackCuts->SetRequireITSPureStandAlone(kFALSE);
      
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(60);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
   
  //DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());

 

}
//_______________________________________________________
AliAnalysisTaskHelium3Pi::~AliAnalysisTaskHelium3Pi() 
{ 
  // Destructor
  if (fListHist) {
    delete fListHist;
    fListHist = 0;
  }
  
  if (fESDtrackCuts) delete fESDtrackCuts;
  if(fNtuple1) delete fNtuple1;
  if(fNtuple4) delete fNtuple4;
}
//=================DEFINITION BETHE BLOCH==============================

Double_t AliAnalysisTaskHelium3Pi::BetheBloch(Double_t betaGamma,Double_t charge,Bool_t isPbPb) {

  Double_t kp1, kp2, kp3, kp4, kp5;
  
  if(isPbPb){

    //    pass1 2011
    kp1 = 4.7*charge*charge;
    kp2 = 8.98482806165147636e+00;
    kp3 = 1.54000000000000005e-05;
    kp4 = 2.30445734159456084e+00;
    kp5 = 2.25624744086878559e+00;

  }
  
  else{

    // to be defined ...
    //pass1 2011
    kp1 = 4.7*charge*charge;
    kp2 = 8.98482806165147636e+00;
    kp3 = 1.54000000000000005e-05;
    kp4 = 2.30445734159456084e+00;
    kp5 = 2.25624744086878559e+00;

  }

  Double_t beta = betaGamma / TMath::Sqrt(1.0 + betaGamma * betaGamma);
  
  Double_t aa = TMath::Power(beta, kp4);
  Double_t bb = TMath::Power(1.0 / betaGamma, kp5);
  
  bb = TMath::Log(kp3 + bb);
  
  Double_t out = (kp2 - aa - bb) * kp1 / aa;

  return out;
 
}



//____________________________________________________________________

// Flattening of the centrality distribution
// true =  means skip event

Bool_t AliAnalysisTaskHelium3Pi::Flatten(Float_t cent) {
  float prob[13] = {
    0.855566,0.846964,0.829618,0.829259,0.830984,
    0.85094,0.844346,0.851818,0.874758,1,
    0.374767,0.650491,0.946963
  };
  if (cent >= 13.f) return kFALSE;
  else return gRandom->Rndm() > prob[int(cent)];
}


//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskHelium3Pi::UserCreateOutputObjects()
{

  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!

  if(! fHistEventMultiplicity ){
    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 14 , -0.5, 13.5);
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1 ,"All Events");          //0
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2 ,"Events w/PV");         //1
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3 ,"Events w/|Vz|<10cm");  //2
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4 ,"Central Events");      //3
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5 ,"SemiCentral Events");  //4
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6 ,"MB Events");           //5
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7 ,"INT7 Events");         //6
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8 ,"Central Events  w/|Vz|<10cm"); //7
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9 ,"SemiCentral Events  w/|Vz|<10cm"); //8
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10 ,"MB Events w/|Vz|<10cm");  //9
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"INT7 Events w/|Vz|<10cm"); //10
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"Any Events");              //11
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(13,"Any Events w/|Vz|<10cm");  //12

    fListHist->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity   = new TH2F( "fHistTrackMultiplicity" , "Nb of Tracks", 2500,0, 25000,210,-1,104);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicity);
  } 

  if(! fHistTrackMultiplicityCent ){
    fHistTrackMultiplicityCent   = new TH2F( "fHistTrackMultiplicityCent", "Nb of Tracks Central Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityCent->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityCent);
  } 

  if(! fHistTrackMultiplicitySemiCent ){
    fHistTrackMultiplicitySemiCent   = new TH2F( "fHistTrackMultiplicitySemiCent" , "Nb of Tracks SemiCentral Events", 2500,0, 25000 ,210,-1,104);
    fHistTrackMultiplicitySemiCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicitySemiCent->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicitySemiCent);
  } 
 
  if(! fHistTrackMultiplicityMB ){
    fHistTrackMultiplicityMB   = new TH2F( "fHistTrackMultiplicityMB" , "Nb of Tracks MBral Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityMB->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityMB);
  } 

  if(! fHistTrackMultiplicityINT7 ){
    fHistTrackMultiplicityINT7   = new TH2F( "fHistTrackMultiplicityINT7" , "Nb of Tracks INT7ral Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityINT7->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityINT7->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityINT7);
  } 
  
  if(! fHistTrackMultiplicityPVCent ){
    fHistTrackMultiplicityPVCent   = new TH2F( "fHistTrackMultiplicityPVCent" , "Nb of Tracks Central Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityPVCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVCent->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityPVCent);
  } 

  if(! fHistTrackMultiplicityPVSemiCent ){
    fHistTrackMultiplicityPVSemiCent   = new TH2F( "fHistTrackMultiplicityPVSemiCent" , "Nb of Tracks SemiCentral Events", 2500,0, 25000 ,210,-1,104);
    fHistTrackMultiplicityPVSemiCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVSemiCent->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityPVSemiCent);
  } 
 
  if(! fHistTrackMultiplicityPVMB ){
    fHistTrackMultiplicityPVMB   = new TH2F( "fHistTrackMultiplicityPVMB" , "Nb of Tracks MBral Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityPVMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVMB->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityPVMB);
  } 

  if(! fHistTrackMultiplicityPVINT7 ){
    fHistTrackMultiplicityPVINT7   = new TH2F( "fHistTrackMultiplicityPVINT7" , "Nb of Tracks INT7ral Events", 2500,0, 25000,210,-1,104 );
    fHistTrackMultiplicityPVINT7->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVINT7->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicityPVINT7);
  } 
  
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 120,-6,6,150,0,1500);
    fhBB->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBB->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fhBB);
  }

  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 120,-6,6,100,0,1.2);
    fhTOF->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhTOF->GetYaxis()->SetTitle("#beta");
    fListHist->Add(fhTOF);
  }

  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 600,-5,5);
    fhMassTOF->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    fListHist->Add(fhMassTOF);
  }

  if(! fhBBPions ){
    fhBBPions = new TH2F( "fhBBPions" , "Bethe-Bloch TPC Pions" , 120,-6,6,150,0,1500);
    fhBBPions->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBPions->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fhBBPions);
  }
  
  if(! fhBBHe ){
    fhBBHe = new TH2F( "fhBBHe" , "Bethe-Bloch TPC He" , 240,-6,6,300,0,1500);
    fhBBHe->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBHe->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fhBBHe);
  }

  hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",1000,-5,5,105,0,105);
  fListHist->Add(hQVzAQVzCvsCentrality);

  hqEPCvsCentrality   = new TH2F("hqEPCvsCentrality","hqEPCvsCentrality",100,-5,5,105,0,105);
  hqEPAvsCentrality   = new TH2F("hqEPAvsCentrality","hqEPAvsCentrality",100,-5,5,105,0,105);
  hqEPvsCentrality    = new TH2F("hqEPvsCentrality" ,"hqEPvsCentrality" ,100,-5,5,105,0,105);
  
  fListHist->Add(hqEPCvsCentrality);  
  fListHist->Add(hqEPAvsCentrality);  
  fListHist->Add(hqEPvsCentrality );  

  if(! fNtuple1 ) {
    
    fNtuple1 = new TTree("fNtuple1","fNtuple1");
    
    fNtuple1->Branch("teventtype"           ,&teventtype           ,"teventtype/F");
    fNtuple1->Branch("tTrackNumber"         ,&tTrackNumber         ,"tTrackNumber/F");
    fNtuple1->Branch("tpercentile"          ,&tpercentile          ,"tpercentile/F");
    fNtuple1->Branch("txPrimaryVertex"      ,&txPrimaryVertex      ,"txPrimaryVertex/F");
    fNtuple1->Branch("tyPrimaryVertex"      ,&tyPrimaryVertex      ,"tyPrimaryVertex/F");
    fNtuple1->Branch("tzPrimaryVertex"      ,&tzPrimaryVertex      ,"tzPrimaryVertex/F");
    fNtuple1->Branch("txSecondaryVertex"    ,&txSecondaryVertex    ,"txSecondaryVertex/F");
    fNtuple1->Branch("tySecondaryVertex"    ,&tySecondaryVertex    ,"tySecondaryVertex/F");
    fNtuple1->Branch("tzSecondaryVertex"    ,&tzSecondaryVertex    ,"tzSecondaryVertex/F");
    fNtuple1->Branch("tdcaTracks"           ,&tdcaTracks           ,"tdcaTracks/F");
    fNtuple1->Branch("tCosPointingAngle"    ,&tCosPointingAngle    ,"tCosPointingAngle/F");
    fNtuple1->Branch("tDCAV0toPrimaryVertex",&tDCAV0toPrimaryVertex,"tDCAV0toPrimaryVertex/F");
    fNtuple1->Branch("tHeSign"              ,&tHeSign              ,"tHeSign/F");
    fNtuple1->Branch("tHepInTPC"            ,&tHepInTPC            ,"tHepInTPC/F");
    fNtuple1->Branch("tHeTPCsignal"         ,&tHeTPCsignal         ,"tHeTPCsignal/F");
    fNtuple1->Branch("tDcaHeToPrimVertex"   ,&tDcaHeToPrimVertex   ,"tDcaHeToPrimVertex/F");
    fNtuple1->Branch("tHeEta"               ,&tHeEta               ,"tHeEta/F");
    fNtuple1->Branch("tmomHex"              ,&tmomHex              ,"tmomHex/F");
    fNtuple1->Branch("tmomHey"              ,&tmomHey              ,"tmomHey/F");
    fNtuple1->Branch("tmomHez"              ,&tmomHez              ,"tmomHez/F");
    fNtuple1->Branch("tmomHeAtSVx"          ,&tmomHeAtSVx          ,"tmomHeAtSVx/F");
    fNtuple1->Branch("tmomHeAtSVy"          ,&tmomHeAtSVy          ,"tmomHeAtSVy/F");
    fNtuple1->Branch("tmomHeAtSVz"          ,&tmomHeAtSVz          ,"tmomHeAtSVz/F");
    fNtuple1->Branch("tHeTPCNcls"           ,&tHeTPCNcls           ,"tHeTPCNcls/F");
    fNtuple1->Branch("tHeimpactXY"          ,&tHeimpactXY          ,"tHeimpactXY/F");
    fNtuple1->Branch("tHeimpactZ"           ,&tHeimpactZ           ,"tHeimpactZ/F");
    fNtuple1->Branch("tHeITSClusterMap"     ,&tHeITSClusterMap     ,"tHeITSClusterMap/F");
    fNtuple1->Branch("tIsHeITSRefit"        ,&tIsHeITSRefit        ,"tIsHeITSRefit/F");
    fNtuple1->Branch("tPionSign"            ,&tPionSign            ,"tPionSign/F");
    fNtuple1->Branch("tPionpInTPC"          ,&tPionpInTPC          ,"tPionpInTPC/F");
    fNtuple1->Branch("tPionTPCsignal"       ,&tPionTPCsignal       ,"tPionTPCsignal/F");
    fNtuple1->Branch("tDcaPionToPrimVertex" ,&tDcaPionToPrimVertex ,"tDcaPionToPrimVertex/F");
    fNtuple1->Branch("tPionEta"             ,&tPionEta             ,"tPionEta/F");
    fNtuple1->Branch("tmomPionx"            ,&tmomPionx            ,"tmomPionx/F");
    fNtuple1->Branch("tmomPiony"            ,&tmomPiony            ,"tmomPiony/F");
    fNtuple1->Branch("tmomPionz"            ,&tmomPionz            ,"tmomPionz/F");
    fNtuple1->Branch("tmomNegPionAtSVx"     ,&tmomNegPionAtSVx     ,"tmomNegPionAtSVx/F");
    fNtuple1->Branch("tmomNegPionAtSVy"     ,&tmomNegPionAtSVy     ,"tmomNegPionAtSVy/F");
    fNtuple1->Branch("tmomNegPionAtSVz"     ,&tmomNegPionAtSVz     ,"tmomNegPionAtSVz/F");
    fNtuple1->Branch("tPionTPCNcls"         ,&tPionTPCNcls         ,"tPionTPCNcls/F");
    fNtuple1->Branch("tPionimpactXY"        ,&tPionimpactXY        ,"tPionimpactXY/F");
    fNtuple1->Branch("tPionimpactZ"         ,&tPionimpactZ         ,"tPionimpactZ/F");
    fNtuple1->Branch("tPionITSClusterMap"   ,&tPionITSClusterMap   ,"tPionITSClusterMap/F");
    fNtuple1->Branch("tIsPiITSRefit"        ,&tIsPiITSRefit        ,"tIsPiITSRefit/F");
    fNtuple1->Branch("txn"                  ,&txn                  ,"txn/F");
    fNtuple1->Branch("txp"                  ,&txp                  ,"txp/F");
    fNtuple1->Branch("tuqV0A"               ,&tuqV0A               ,"tuqV0A/F");
    fNtuple1->Branch("tuqV0C"               ,&tuqV0C               ,"tuqV0C/F");
      
  }
  
  if(! fNtuple4 ) {
    
    fNtuple4 = new TTree("fNtuple4","fNtuple4");
      
    fNtuple4->Branch("tHeleventtype"        ,&tHeleventtype        ,"tHeleventtype/F");
    fNtuple4->Branch("tHelpercentile"       ,&tHelpercentile       ,"tHelpercentile/F");
    fNtuple4->Branch("tHelSign"             ,&tHelSign             ,"tHelSign/F");
    fNtuple4->Branch("tHelpinTPC"           ,&tHelpinTPC           ,"tHelpinTPC/F");
    fNtuple4->Branch("tHelGetTPCsignal"     ,&tHelGetTPCsignal     ,"tHelGetTPCsignal/F");
    fNtuple4->Branch("tHelPx"               ,&tHelPx               ,"tHelPx/F");
    fNtuple4->Branch("tHelPy"               ,&tHelPy               ,"tHelPy/F");
    fNtuple4->Branch("tHelPz"               ,&tHelPz               ,"tHelPz/F");
    fNtuple4->Branch("tHelEta"              ,&tHelEta              ,"tHelEta/F");
    fNtuple4->Branch("tHelisTOF"            ,&tHelisTOF            ,"tHelisTOF/F");
    fNtuple4->Branch("tHelTOFpull"          ,&tHelTOFpull          ,"tHelTOFpull/F");
    fNtuple4->Branch("tHeMass"              ,&tHeMass              ,"tHeMass/F");
    fNtuple4->Branch("tHelimpactXY"         ,&tHelimpactXY         ,"tHelimpactXY/F");
    fNtuple4->Branch("tHelimpactZ"          ,&tHelimpactZ          ,"tHelimpactZ/F");
    fNtuple4->Branch("tHelmapITS"           ,&tHelmapITS           ,"tHelmapITS/F");
    fNtuple4->Branch("tHelBetaTOF"          ,&tHelBetaTOF          ,"tHelBetaTOF/F");
    fNtuple4->Branch("tHelIsITSrefit"       ,&tHelIsITSrefit       ,"tHelIsITSrefit/F");
   
  } 

  PostData(1,  fListHist);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple4);
}// end UserCreateOutputObjects


//====================== USER EXEC ========================

void AliAnalysisTaskHelium3Pi::UserExec(Option_t *) 
{
  cout<<"Enter the user exec"<<endl;
  
  // cout<< fAnalysisType<<endl;       
  // cout<< fCollidingSystems<<endl;   
  // cout<< fDataType<<endl;           
  // cout<< fYear<<endl;               
  // cout<< fVzmax<<endl;              
  // cout<< fApplyFlatten<<endl;       
  // cout<< fFill3Hetree<<endl;        
  // cout<< fDoFlow<<endl;                 
  
  //_______________________________________________________________________
  
  //!*********************!//
  //!  Define variables   !//
  //!*********************!//

  Double_t pinTPC=0.,/*poutTPC=0.,*/TPCSignal=0.;
  // Double_t xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;
  Double_t massTOF=0.,timeTOF=0.,/*trackLenghtTOF=0.,*/betaTOF=0.,gamma=0.;

  ULong_t  status=0;
  ULong_t  statusPi=0;

  Bool_t   isTPC=kFALSE,hasTOF=kFALSE,IsHeITSRefit=kFALSE,IsPiITSRefit=kFALSE ;

  Float_t nSigmaNegPion=0.;
  Float_t nSigma3He=0.;

  Double_t cutNSigma = 3;
  // Double_t bbtheoM=0.,bbtheo=0.;
  // Double_t zNathashaNeg=0;
  // Double_t zNathashaPos=0;
  Double_t fPos[3]={0.,0.,0.};
  Double_t tPhi = -999.;

  //Masses
  Double_t        Helium3Mass = 2.80839; 
  //  Double_t        Helium3Mass = 2.80894; //tri-mass
  Double_t        PionMass    = 0.13957; 


  //Flow variables
  
  Double_t qxEPa = -999., qyEPa = -999.;
  Double_t qxEPc = -999., qyEPc = -999.;
  Double_t qxEP =  -999., qyEP = -999.;
  
  Double_t evPlAngV0A = -999.;
  Double_t evPlAngV0C = -999.;
  Double_t evPlAngV0  = -999.;

  Float_t  uqV0A = -999.;
  Float_t  uqV0C = -999.; 
  
  // TLORENTZ vectors
  
  TLorentzVector  vPion,vHelium,vSum;

  //!----------------------------------------------------------------

  //! A set of very loose parameters for cuts 
  
  Double_t fgChi2max=33.;     //! max chi2
  Double_t fgDNmin=0.05;      //! min imp parameter for the 1st daughter = 500um
  Double_t fgDCAmax=1.0;      //! max DCA between the daughter tracks in cm
  Double_t fgCPAmin=0.99;     //! min cosine of V0's pointing angle  
  //  Double_t fgRmin=0.2;    //! min radius of the fiducial volume //original
  Double_t fgRmin=0.1;        //! min radius of the fiducial volume = 1 mm 
  Double_t fgRmax=200.;       //! max radius of the fiducial volume = 2 m

  //------------------------------------------

  // Main loop
  // Called for EACH event
  Info("AliAnalysisTaskHelium3Pi","Starting UserExec");  

  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }

  if(fAnalysisType == "ESD"){
    fESDevent = dynamic_cast<AliESDEvent*>(event);
    if (!fESDevent) {
      AliError("Cannot get the ESD event");
      return;
    }  
    fevent = fESDevent;
  }
  else{
    AliError("Cannot get any event");
    return;
  }
  
  fHistEventMultiplicity->Fill(0);
  
  Double_t lMagneticField=fevent->GetMagneticField();
  Int_t TrackNumber = -1;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  //----------------------------------------
  // Centrality  
 
  Float_t percentile = -999;
  AliMultSelection *fMultSelection = 0x0; 
  
  if(fYear < 2015){
    AliCentrality *centrality = fevent->GetCentrality();
    percentile = centrality->GetCentralityPercentile("V0M");
    if(fApplyFlatten == kTRUE){
      if(fYear == 2011){
	if(Flatten(percentile))return;
      }
    }
  }
  else{
    fMultSelection = (AliMultSelection * ) fevent->FindListObject("MultSelection");
    if(!fMultSelection){
      AliWarning("AliMultSelection object not found!");
      PostData(1,fListHist);
      PostData(2,fNtuple1);
      PostData(3,fNtuple4);
      return;
    }
    else{
      percentile = fMultSelection->GetMultiplicityPercentile("V0M");
    } 
    
  }
  


  // AliCentrality *centrality = fevent->GetCentrality();
  // Float_t percentile = centrality->GetCentralityPercentile("V0M");
  // if(fApplyFlatten == kTRUE){
  //   if(fYear == 2011){
  //     if(Flatten(percentile))return;
  //   }
  // }
  TrackNumber = fevent->GetNumberOfTracks();
  if (TrackNumber<2) return;  

  //---------------------------------------------
  // PID
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse(); 
  
  //===========================================
  
  Int_t eventtype=-99;
  
  Bool_t isSelectedCentral     = (inputHandler->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (inputHandler->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (inputHandler->IsEventSelected() & AliVEvent::kMB);
  Bool_t isSelectedINT7        = (inputHandler->IsEventSelected() & AliVEvent::kINT7);
  Bool_t isSelectedAny         = (inputHandler->IsEventSelected() & AliVEvent::kAny);
 
  if(isSelectedCentral){
    fHistEventMultiplicity->Fill(3);
    fHistTrackMultiplicityCent->Fill(TrackNumber,percentile); 
    eventtype=1;
  }

  if(isSelectedSemiCentral){
    fHistEventMultiplicity->Fill(4);
    fHistTrackMultiplicitySemiCent->Fill(TrackNumber,percentile); 
    eventtype=2;
  }

  if(isSelectedMB){
    fHistEventMultiplicity->Fill(5);
    fHistTrackMultiplicityMB->Fill(TrackNumber,percentile); 
    eventtype=3;
  }

  if(isSelectedINT7){
    fHistEventMultiplicity->Fill(6);
    fHistTrackMultiplicityINT7->Fill(TrackNumber,percentile); 
    eventtype=4;
  }
 
  if(!isSelectedCentral && !isSelectedSemiCentral && !isSelectedMB && !isSelectedINT7 && isSelectedAny){
    fHistEventMultiplicity->Fill(11);
    fHistTrackMultiplicity->Fill(TrackNumber,percentile); //tracce per evento
    eventtype=5;
  }

  //if(isSelectedCentral || isSelectedSemiCentral || isSelectedMB || isSelectedAny){
  if(eventtype ==1  || eventtype ==2  || eventtype==3  || eventtype==4 || eventtype==5){
    
    // ANALISYS
    // Primary vertex cut

    const AliVVertex* vertexmain = fevent->GetPrimaryVertex();
    if (!vertexmain){
      AliWarning("No prim. vertex in ESD... return!");

      PostData(1,fListHist);
      PostData(2,fNtuple1);
      PostData(3,fNtuple4);
      return;
    }
    
    vertexmain->GetXYZ( lBestPrimaryVtxPos );
    
    fHistEventMultiplicity->Fill(1); // analyzed events with PV
 
    if((TMath::Abs(lBestPrimaryVtxPos[2])) > fVzmax) return;
    
    //------------- Pile-up
    Bool_t isPileUpSpd=kFALSE;
    isPileUpSpd=fESDevent->IsPileupFromSPD();
       
    if(isPileUpSpd){  
      PostData(1,fListHist);
      PostData(2,fNtuple1);
      PostData(3,fNtuple4);
      return;
    }
    //--------------------------


    if(eventtype==1){
      fHistTrackMultiplicityPVCent->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(7); 
    }
    
    if(eventtype==2){
      fHistTrackMultiplicityPVSemiCent->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(8); 
    }
    
    if(eventtype==3){
      fHistTrackMultiplicityPVMB->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(9); 
    }

    if(eventtype==4){
      fHistEventMultiplicity->Fill(10); 
    }
    if(eventtype==5){
      fHistEventMultiplicity->Fill(12); 
    }
    
    
    fHistEventMultiplicity->Fill(2);
    
    //Fill Flow variable
    if(fDoFlow == kTRUE){
      AliEventplane *pl = fevent->GetEventplane();
      if(!pl ){
	AliError("AliAnalysisTaskSENucleiv2SP::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
	fHistEventMultiplicity->Fill(12);
      }
      
      qxEPa = -999., qyEPa = -999.;
      qxEPc = -999., qyEPc = -999.;
      qxEP =  -999. , qyEP = -999.;
      
      evPlAngV0A = pl->CalculateVZEROEventPlane(fevent, 8, 2, qxEPa, qyEPa);
      evPlAngV0C = pl->CalculateVZEROEventPlane(fevent, 9, 2, qxEPc, qyEPc);
      evPlAngV0  = pl->CalculateVZEROEventPlane(fevent,10, 2, qxEP,  qyEP);

      Double_t  QV0AQV0C = qxEPa * qxEPc + qyEPa*qyEPc;
      hQVzAQVzCvsCentrality->Fill(QV0AQV0C,percentile);

      hqEPCvsCentrality  ->Fill(TMath::Sqrt(qxEPa*qxEPa+qyEPa*qyEPa) , percentile); 
      hqEPAvsCentrality  ->Fill(TMath::Sqrt(qxEPc*qxEPc+qyEPc*qyEPc) , percentile); 
      hqEPvsCentrality   ->Fill(TMath::Sqrt(qxEP *qxEP +qyEP *qyEP ) , percentile); 
      
    }
    
    //Find Pair candidates
    
    TArrayI PionsTPC(TrackNumber);        //Neg pions
    Int_t nPionsTPC=0;
    
    TArrayI HeTPC(TrackNumber);        //helium3
    Int_t nHeTPC=0;
    
    // const Double_t speedOfLight =  TMath::C()*1E2*1E-12; // cm/ps
    
    Float_t impactXY=-999, impactZ=-999;
    Float_t impactXYpi=-999, impactZpi=-999;

    // Double_t ptcExp  = -999;
    // Double_t expbeta = -999;
    Double_t pullTOF = -999; 

    // AliESDtrack *track;

    //*************************************************************

    for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
     
     
      AliESDtrack *track=fESDevent->GetTrack(j);       
      if(!track)continue;
      
      if (track->Eta() < -0.9 || track->Eta() > 0.9) continue;
      
      Bool_t trkFlag = 0;
      
      trkFlag = fESDtrackCuts->AcceptTrack(track);
	
     
      if(!trkFlag)
       	continue;
      
      // cout<<track<<endl;
      status  = (ULong_t)track->GetStatus();
      isTPC   = (((status) & AliESDtrack::kTPCin)  != 0);
     
      Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
      hasTOF     = kFALSE;
      if (hasTOFout) hasTOF = kTRUE;
      Float_t trackLenghtTOF = track->GetIntegratedLength(); 
      if (trackLenghtTOF < 350.) hasTOF = kFALSE;

      UInt_t mapITS=track->GetITSClusterMap();
      
      //----------------------------------------------
      //****** Cuts from  AliV0Vertex.cxx *************
     
      Double_t d=-999;
     
      d=track->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lMagneticField);
      
      //Cut on d moved into the PID loops
      
      //---- (Usefull) Stuff
      
      //      TPCSignal=atrack->GetTPCsignal(); 
      TPCSignal=track->GetTPCsignal(); 
      
      if (TPCSignal<10)continue;
      if (TPCSignal>1000)continue;
      /*
      //      if(!isTPC)continue;
      if(fAnalysisType == "ESD"){
	if(!esdtrack->GetTPCInnerParam())continue;
      
	AliExternalTrackParam trackIn(*esdtrack->GetInnerParam()); 
	pinTPC= trackIn.GetP(); 
      }

      else  if(fAnalysisType == "AODD"){
	pinTPC= aodtrack->GetTPCMomentum();
      }
      poutTPC=pinTPC;
      */
      pinTPC = track->GetTPCmomentum();
      //      cout<<TPCSignal<<endl;
      
      //   if((status) & (AliVTrack::kITSrefit!=0)){
      fhBB->Fill(pinTPC*track->Charge(),TPCSignal);
	// }
      
      Double_t p    = track->P();
      timeTOF = track->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);      // ps
      
      if(hasTOF){
	
	betaTOF= trackLenghtTOF/(timeTOF * 2.99792457999999984e-02);
	//	cout<<trackLenghtTOF<<" "<<betaTOF<<" "<<pinTPC<<endl;
	fhTOF->Fill(pinTPC*track->GetSign(),betaTOF);
	gamma =  1/TMath::Sqrt(1 - betaTOF*betaTOF);
	massTOF = pinTPC/TMath::Sqrt(gamma*gamma - 1);
      }
      
      
      nSigmaNegPion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) 2));
      
      //2 is pion
      
      if ( (nSigmaNegPion < cutNSigma)){ 
	
	if (TMath::Abs(d)>fgRmax) continue;
	
	fhBBPions->Fill(pinTPC*track->GetSign(),TPCSignal);
	
	if(pinTPC<3.){
	  PionsTPC[nPionsTPC++]=j;
	}
      }
    
      Bool_t isHeITSrefit=((status) & (AliESDtrack::kITSrefit));
      
      nSigma3He  = TMath::Abs((fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType) 7)));
      
      if(nSigma3He  < 3.) {
	
	if(pinTPC>0.6){ 
	  //cout<<hasTOF<<" massTOF: "<<massTOF<<endl;
	  if(hasTOF){
	    // if(TMath::Abs(massTOF) > 5.0)continue;
	    // if(TMath::Abs(massTOF) < 1.8 )continue;
	    fhMassTOF->Fill(massTOF*track->Charge());
	    // cout<<massTOF<<endl;
	  }
	
	//------------------------------
	
	fhBBHe->Fill(pinTPC*track->Charge(),TPCSignal);
	// if(pinTPC<1.2)
	//   HeTPC[nHeTPC++]=j;
	// else
	//   if(hasTOF && TMath::Abs(pullTOF)<=3)
	//     HeTPC[nHeTPC++]=j;
	
	//atrack->GetImpactParameters(impactXY, impactZ);
	
	track->GetImpactParameters(impactXY, impactZ);
	  
	tHeleventtype	 =(Float_t)eventtype;
	tHelpercentile   =(Float_t)percentile;
	tHelSign	 =(Float_t)track->GetSign();
	tHelpinTPC	 =(Float_t)pinTPC;
	tHelGetTPCsignal =(Float_t)track->GetTPCsignal();
	tHelPx	         =(Float_t)track->Px();
	tHelPy	         =(Float_t)track->Py();
	tHelPz	         =(Float_t)track->Pz();
	tHelEta	         =(Float_t)track->Eta();
	tHelisTOF	 =(Float_t)hasTOF;
	tHelTOFpull	 =(Float_t)pullTOF;
	tHeMass          =(Float_t)massTOF;
	tHelimpactXY	 =(Float_t)impactXY;
	tHelimpactZ	 =(Float_t)impactZ;
	tHelmapITS	 =(Float_t)mapITS;
	tHelBetaTOF      =(Float_t)betaTOF;
	tHelIsITSrefit   =(Float_t)isHeITSrefit;

	if(fFill3Hetree == kTRUE)
	  fNtuple4->Fill();
	
	//this cut was outside the loop
	if (TMath::Abs(d)>fgRmax) continue;
	HeTPC[nHeTPC++]=j;
	
	}
      }
    }  //! track
    
    PionsTPC.Set(nPionsTPC);
    HeTPC.Set(nHeTPC);
    
    Double_t        DcaHeToPrimVertex=0;
    Double_t        DcaPionToPrimVertex=0;
    
    impactXY=-999, impactZ=-999;
    impactXYpi=-999, impactZpi=-999;
    
    // Tracks
    
    // Vettors for il PxPyPz
    
    Double_t momPionVett[3];
    for(Int_t i=0;i<3;i++)momPionVett[i]=0;
    
    Double_t momHeVett[3];
    for(Int_t i=0;i<3;i++)momHeVett[i]=0;
    
    //At SV
    
    Double_t momPionVettAt[3];
    for(Int_t i=0;i<3;i++)momPionVettAt[i]=0;
    
    Double_t momHeVettAt[3];
    for(Int_t i=0;i<3;i++)momHeVettAt[i]=0;
    
    AliESDtrack  *PionTrack = 0x0;
    AliESDtrack  *HeTrack = 0x0;
    //---------------   LOOP PAIRS   ----------------
    
    for (Int_t k=0; k < nPionsTPC; k++) {                           //! Pions Loop
      
      DcaPionToPrimVertex=0.;
      DcaHeToPrimVertex=0;
      
      Int_t PionIdx=PionsTPC[k];
      
      PionTrack =fESDevent->GetTrack(PionIdx);

      //      AliVTrack* PionTrack = (AliVTrack*) fevent->GetTrack(PionIdx);
      //AliESDtrack *PionTrack = static_cast<AliESDtrack *>(aPionTrack);
      //  PionTrack = static_cast<AliESDtrack *>(aPionTrack);
      //PionTrack = new AliESDtrack((AliVTrack*)aPionTrack);--
          
      statusPi = (ULong_t)PionTrack->GetStatus();

      IsPiITSRefit = ((statusPi) & (AliESDtrack::kITSrefit)); 
     
 
      if (PionTrack) {
	DcaPionToPrimVertex = TMath::Abs(PionTrack->GetD(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1],lMagneticField)); //OK
	}
      
      if(DcaPionToPrimVertex<0.2)continue; 
      
      
      AliExternalTrackParam trackInPion(*PionTrack);  

      // AliExternalTrackParam trackInPion;
      // trackInPion.CopyFromVTrack(aPionTrack);

      for (Int_t i=0; i<nHeTPC; i++){                               //! Helium Loop
	
	Int_t HeIdx=HeTPC[i];
	
	HeTrack	= fESDevent->GetTrack(HeIdx);

	IsHeITSRefit = (status & AliESDtrack::kITSrefit); 
	
	if (HeTrack){
	  if(fAnalysisType == "ESD"){  
	    DcaHeToPrimVertex = TMath::Abs(HeTrack->GetD(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1],lMagneticField)); //OK
	  }
	}
	AliExternalTrackParam trackInHe(*HeTrack); 

	// AliExternalTrackParam trackInHe;
	// trackInHe.CopyFromVTrack(aHeTrack);

	if ( DcaPionToPrimVertex < fgDNmin)                //OK
	  if ( DcaHeToPrimVertex < fgDNmin) continue;    //OK
	
	//if(DcaHeToPrimVertex<0.1)continue; 
	
	Double_t xn, xp;
	Double_t dca=0.;
	
	dca= PionTrack->GetDCA(HeTrack,lMagneticField,xn,xp); //!dca (Neg to Pos)
	
	if (dca > fgDCAmax) continue;
	if ((xn+xp) > 2*fgRmax) continue;
	if ((xn+xp) < 2*fgRmin) continue;
	
	//CORRECTION from AliV0Vertex
	
	Bool_t corrected=kFALSE;
	if ((trackInPion.GetX() > 3.) && (xn < 3.)) {
	  //correct for the beam pipe material
	  corrected=kTRUE;
	}
	if ((trackInHe.GetX() > 3.) && (xp < 3.)) {
	  //correct for the beam pipe material
	  corrected=kTRUE;
	}
	if (corrected) {
	  dca=trackInPion.GetDCA(&trackInHe,lMagneticField,xn,xp);
	  if (dca > fgDCAmax) continue;
	  if ((xn+xp) > 2*fgRmax) continue;
	  if ((xn+xp) < 2*fgRmin) continue;
	}
	
	//=============================================//
	// Make "V0" with found tracks                 //
	//=============================================//
	
	trackInPion.PropagateTo(xn,lMagneticField); 
	trackInHe.PropagateTo(xp,lMagneticField);
	
	AliESDv0 vertex(trackInPion,PionIdx,trackInHe,HeIdx);
	if (vertex.GetChi2V0() > fgChi2max) continue;
	
	Float_t CosPointingAngle=vertex.GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]); //PointingAngle
	if (CosPointingAngle < fgCPAmin) continue;
	
	vertex.SetDcaV0Daughters(dca);
	vertex.SetV0CosineOfPointingAngle(CosPointingAngle);
	
	fPos[0]=vertex.Xv();
	fPos[1]=vertex.Yv(); 
	fPos[2]=vertex.Zv(); 
	
	HeTrack->PxPyPz(momHeVett);
	PionTrack->PxPyPz(momPionVett); 
	
	Double_t raggio=TMath::Sqrt(fPos[0]*fPos[0]+fPos[1]*fPos[1]+fPos[2]*fPos[2]);
	HeTrack->GetPxPyPzAt(raggio,lMagneticField,momHeVettAt);
	PionTrack->GetPxPyPzAt(raggio,lMagneticField,momPionVettAt); 
	
	//------------------------------------------------------------------------//
	HeTrack->GetImpactParameters(impactXY, impactZ);
	PionTrack->GetImpactParameters(impactXYpi, impactZpi);
	
	if(vertex.GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2])>3) continue;	
	//save only pairs up to 3.2 GeV/c2

	vHelium.SetXYZM(2*momHeVettAt[0],2*momHeVettAt[1],2*momHeVettAt[2],Helium3Mass); 
	vPion.SetXYZM(momPionVettAt[0],momPionVettAt[1],momPionVettAt[2],PionMass);       
	vSum=vHelium+vPion;
	
	if(vSum.M()>3.2)
	  continue;

	if(fDoFlow){
	  tPhi = vSum.Phi();
	  uqV0A = -999;
	  uqV0C = -999; 
	  
	  uqV0A = TMath::Cos(2*tPhi)*qxEPa+TMath::Sin(2*tPhi)*qyEPa;
	  uqV0C = TMath::Cos(2*tPhi)*qxEPc+TMath::Sin(2*tPhi)*qyEPc;
	}

	// Int_t  fIdxInt[200]; //dummy array
	// Int_t nClustersTPCHe = HeTrack->GetTPCclusters(fIdxInt);
	// Int_t nClustersTPCPi = PionTrack->GetTPCclusters(fIdxInt);
	
	//----------------------------------------------------------------------//

	teventtype		=(Float_t)eventtype;
	tTrackNumber		=(Float_t)TrackNumber;
	tpercentile		=(Float_t)percentile;
	txPrimaryVertex	   	=(Float_t)lBestPrimaryVtxPos[0]; //PRIMARY
	tyPrimaryVertex	   	=(Float_t)lBestPrimaryVtxPos[1];
	tzPrimaryVertex	   	=(Float_t)lBestPrimaryVtxPos[2];
	txSecondaryVertex	=(Float_t)fPos[0]; //SECONDARY
	tySecondaryVertex	=(Float_t)fPos[1];
	tzSecondaryVertex	=(Float_t)fPos[2];
	tdcaTracks		=(Float_t)dca;           //between 2 tracks
	tCosPointingAngle	=(Float_t)CosPointingAngle;          //cosPointingAngle da V0
	tDCAV0toPrimaryVertex	=(Float_t)vertex.GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
	tHeSign		   	=(Float_t)HeTrack->GetSign(); //He
	tHepInTPC		=(Float_t)trackInHe.GetP();
	tHeTPCsignal		=(Float_t)HeTrack->GetTPCsignal();
	tDcaHeToPrimVertex	=(Float_t)DcaHeToPrimVertex;
	tHeEta		   	=(Float_t)HeTrack->Eta();
	tmomHex		   	=(Float_t)momHeVett[0];
	tmomHey		   	=(Float_t)momHeVett[1];
	tmomHez		   	=(Float_t)momHeVett[2];
	tmomHeAtSVx		=(Float_t)momHeVettAt[0];
	tmomHeAtSVy		=(Float_t)momHeVettAt[1];
	tmomHeAtSVz		=(Float_t)momHeVettAt[2];
	tHeTPCNcls		=(Float_t)HeTrack->GetTPCNcls();
	tHeimpactXY		=(Float_t)impactXY;
	tHeimpactZ		=(Float_t)impactZ;
	tHeITSClusterMap	=(Float_t)HeTrack->GetITSClusterMap();
	tIsHeITSRefit		=(Float_t)IsHeITSRefit;
	tPionSign		=(Float_t)PionTrack->GetSign(); //Pion
	tPionpInTPC		=(Float_t)trackInPion.GetP();
	tPionTPCsignal	   	=(Float_t)PionTrack->GetTPCsignal();
	tDcaPionToPrimVertex	=(Float_t)DcaPionToPrimVertex;
	tPionEta		=(Float_t)PionTrack->Eta();
	tmomPionx		=(Float_t)momPionVett[0];
	tmomPiony		=(Float_t)momPionVett[1];
	tmomPionz		=(Float_t)momPionVett[2];
	tmomNegPionAtSVx	=(Float_t)momPionVettAt[0];
	tmomNegPionAtSVy	=(Float_t)momPionVettAt[1];
	tmomNegPionAtSVz	=(Float_t)momPionVettAt[2];
	tPionTPCNcls		=(Float_t)PionTrack->GetTPCNcls();
	tPionimpactXY		=(Float_t)impactXYpi;
	tPionimpactZ		=(Float_t)impactZpi;
	tPionITSClusterMap	=(Float_t)PionTrack->GetITSClusterMap();
	tIsPiITSRefit		=(Float_t)IsPiITSRefit;
	txn			=(Float_t)xn;
	txp			=(Float_t)xp;
	tuqV0A		   	=(Float_t)uqV0A;
	tuqV0C                  =(Float_t)uqV0C;
	
	fNtuple1->Fill();  
	vertex.Delete();
      }// positive TPC
      
    } //negative tpc
    
  }
  
  PostData(1,fListHist);
  PostData(2,fNtuple1);
  PostData(3,fNtuple4);
  
} //end userexec


//________________________________________________________________________

void AliAnalysisTaskHelium3Pi::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}


