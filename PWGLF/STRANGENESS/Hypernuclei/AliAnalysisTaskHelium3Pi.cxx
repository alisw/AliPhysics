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
//
//-----------------------------------------------------------------
//                 AliAnalysisTaskHelium3Pi class
//-----------------------------------------------------------------


class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0; 
class AliCascadeVertexer;

#include <iostream>
#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TCutG.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAnalysisTaskHelium3Pi.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "TString.h"
#include <TDatime.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
//#include <AliVTrack.h>

const Int_t AliAnalysisTaskHelium3Pi::fgNrot = 15;

ClassImp(AliAnalysisTaskHelium3Pi)

//________________________________________________________________________
AliAnalysisTaskHelium3Pi::AliAnalysisTaskHelium3Pi() 
: AliAnalysisTaskSE(),
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fESDtrackCuts(0),
  fDataType("REAL"),
  fListHist(0), 
  fHistEventMultiplicity(0),         
  fHistTrackMultiplicity(0),      
  fHistTrackMultiplicityCent(0),      
  fHistTrackMultiplicitySemiCent(0),  
  fHistTrackMultiplicityMB(0),        
  fHistTrackMultiplicityPVCent(0),      
  fHistTrackMultiplicityPVSemiCent(0),  
  fHistTrackMultiplicityPVMB(0),        
  fhBB(0),    
  fhTOF(0),   
  fhMassTOF(0),
  fhBBPions(0),
  fhBBHe(0),   
  fhNaPos(0),  
  fhNaNeg(0),  
  fBetavsTPCsignalPos(0),  
  fBetavsTPCsignalNeg(0),  
  fNtuple1(0),
  trunNumber(0),
  tbunchcross(0),
  torbit(0),
  tperiod(0),
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
  tchi2He(0),
  tchi2Pi(0),
  fNtuple4(0),
  tHelrunNumber(0),
  tHelBCNumber(0),
  tHelOrbitNumber(0),
  tHelPeriodNumber(0),
  tHeleventtype(0),
  tHelisHeITSrefit(0),
  tHelpercentile(0),
  tHelSign(0),
  tHelpinTPC(0),
  tHelGetTPCsignal(0),
  tHelPx(0),
  tHelPy(0),
  tHelPz(0),
  tHelEta(0),
  tHelisTOF(0),
  tHelpoutTPC(0),
  tHeltimeTOF(0),
  tHeltrackLenghtTOF(0),
  tHelimpactXY(0),
  tHelimpactZ(0),
  tHelmapITS(0),
  tHelTPCNcls(0),
  tHelTRDsignal(0),
  tHelxPrimaryVertex(0),
  tHelyPrimaryVertex(0),
  tHelzPrimaryVertex(0),
  tHelchi2PerClusterTPC(0),
  fPIDResponse(0)
  
{
  // Dummy Constructor
  
  //  printf("Dummy Constructor");
  
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
AliAnalysisTaskHelium3Pi::AliAnalysisTaskHelium3Pi(const char *name) 
: AliAnalysisTaskSE(name), 
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fESDtrackCuts(0),
  fDataType("REAL"),
  fListHist(0), 
  fHistEventMultiplicity(0),         
  fHistTrackMultiplicity(0),      
  fHistTrackMultiplicityCent(0),      
  fHistTrackMultiplicitySemiCent(0),  
  fHistTrackMultiplicityMB(0),        
  fHistTrackMultiplicityPVCent(0),      
  fHistTrackMultiplicityPVSemiCent(0),  
  fHistTrackMultiplicityPVMB(0),        
  fhBB(0),    
  fhTOF(0),   
  fhMassTOF(0),
  fhBBPions(0),
  fhBBHe(0),   
  fhNaPos(0),  
  fhNaNeg(0),  
  fBetavsTPCsignalPos(0),  
  fBetavsTPCsignalNeg(0),  
  fNtuple1(0),
  trunNumber(0),
  tbunchcross(0),
  torbit(0),
  tperiod(0),
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
  tchi2He(0),
  tchi2Pi(0),
  fNtuple4(0),
  tHelrunNumber(0),
  tHelBCNumber(0),
  tHelOrbitNumber(0),
  tHelPeriodNumber(0),
  tHeleventtype(0),
  tHelisHeITSrefit(0),
  tHelpercentile(0),
  tHelSign(0),
  tHelpinTPC(0),
  tHelGetTPCsignal(0),
  tHelPx(0),
  tHelPy(0),
  tHelPz(0),
  tHelEta(0),
  tHelisTOF(0),
  tHelpoutTPC(0),
  tHeltimeTOF(0),
  tHeltrackLenghtTOF(0),
  tHelimpactXY(0),
  tHelimpactZ(0),
  tHelmapITS(0),
  tHelTPCNcls(0),
  tHelTRDsignal(0),
  tHelxPrimaryVertex(0),
  tHelyPrimaryVertex(0),
  tHelzPrimaryVertex(0),
  tHelchi2PerClusterTPC(0),
  fPIDResponse(0)
{					  

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container ()

  DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());

  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireITSStandAlone(kFALSE);
  fESDtrackCuts->SetRequireITSPureStandAlone(kFALSE);
      
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(60);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);

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

//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskHelium3Pi::UserCreateOutputObjects()
{

  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!

  if(! fHistEventMultiplicity ){
    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 10 , 0, 10);
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"SemiCentral Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"MB Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"Central Events  w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"SemiCentral Events  w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"MB Events   w/|Vz|<10cm");

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
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 300,0 ,5);
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
    fhBBHe = new TH2F( "fhBBHe" , "Bethe-Bloch TPC He" , 120,-6,6,150,0,1500);
    fhBBHe->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBHe->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fhBBHe);
  }
  
  if(! fhNaPos ){
    fhNaPos = new TH2F( "fhNaPos" , "Distribution Pos" , 500,0,5,40,-10,10);
    fhNaPos->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhNaPos->GetYaxis()->SetTitle("(TPCSignal-bbtheo)/bbtheo (He)");
    fListHist->Add(fhNaPos);
  }
  
  if(! fhNaNeg ){
    fhNaNeg = new TH2F( "fhNaNeg" , "Distribution Neg" , 500,0,5,40,-10,10);
    fhNaNeg->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhNaNeg->GetYaxis()->SetTitle("(TPCSignal-bbtheo)/bbtheo (He)");
    fListHist->Add(fhNaNeg);
  }

  if(! fBetavsTPCsignalPos ){
    fBetavsTPCsignalPos = new TH2F("fBetavsTPCsignalPos","fBetavsTPCsignalPos",100,0,1.2,150,0,1500);
    fBetavsTPCsignalPos->GetXaxis()->SetTitle("#beta");
    fBetavsTPCsignalPos->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fBetavsTPCsignalPos);
  }
  
  if(! fBetavsTPCsignalNeg ){
    fBetavsTPCsignalNeg = new TH2F("fBetavsTPCsignalNeg","fBetavsTPCsignalNeg",100,0,1.2,150,0,1500);
    fBetavsTPCsignalNeg->GetXaxis()->SetTitle("#beta");
    fBetavsTPCsignalNeg->GetYaxis()->SetTitle("TPC Signal");
    fListHist->Add(fBetavsTPCsignalNeg);
  }
  

  
  if(! fNtuple1 ) {
    
    fNtuple1 = new TTree("fNtuple1","fNtuple1");
    
    fNtuple1->Branch("trunNumber"           ,&trunNumber           ,"trunNumber/F");
    fNtuple1->Branch("tbunchcross"          ,&tbunchcross          ,"tbunchcross/F");
    fNtuple1->Branch("torbit"               ,&torbit               ,"torbit/F");
    fNtuple1->Branch("tperiod"              ,&tperiod              ,"tperiod/F");
    fNtuple1->Branch("teventtype"           ,&teventtype           ,"teventtype/F");
    fNtuple1->Branch("tTrackNumber"         ,&tTrackNumber         ,"tTrackNumber/F");
    fNtuple1->Branch("tpercentile"          ,&tpercentile          ,"tpercentile/F") ;
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
    fNtuple1->Branch("tchi2He"              ,&tchi2He              ,"tchi2He/F");
    fNtuple1->Branch("tchi2Pi"              ,&tchi2Pi              ,"tchi2Pi/F");
    
  }
  
  if(! fNtuple4 ) {
 
    fNtuple4 = new TTree("fNtuple4","fNtuple4");
    
    fNtuple4->Branch("tHelrunNumber"        ,&tHelrunNumber        ,"tHelrunNumber/F");
    fNtuple4->Branch("tHelBCNumber"         ,&tHelBCNumber         ,"tHelBCNumber/F");
    fNtuple4->Branch("tHelOrbitNumber"      ,&tHelOrbitNumber      ,"tHelOrbitNumber/F");
    fNtuple4->Branch("tHelPeriodNumber"     ,&tHelPeriodNumber     ,"tHelPeriodNumber/F");
    fNtuple4->Branch("tHeleventtype"        ,&tHeleventtype        ,"tHeleventtype/F");
    fNtuple4->Branch("tHelisHeITSrefit"     ,&tHelisHeITSrefit     ,"tHelisHeITSrefit/F");
    fNtuple4->Branch("tHelpercentile"       ,&tHelpercentile       ,"tHelpercentile/F");
    fNtuple4->Branch("tHelSign"             ,&tHelSign             ,"tHelSign/F");
    fNtuple4->Branch("tHelpinTPC"           ,&tHelpinTPC           ,"tHelpinTPC/F");
    fNtuple4->Branch("tHelGetTPCsignal"     ,&tHelGetTPCsignal     ,"tHelGetTPCsignal/F");
    fNtuple4->Branch("tHelPx"               ,&tHelPx               ,"tHelPx/F");
    fNtuple4->Branch("tHelPy"               ,&tHelPy               ,"tHelPy/F");
    fNtuple4->Branch("tHelPz"               ,&tHelPz               ,"tHelPz/F");
    fNtuple4->Branch("tHelEta"              ,&tHelEta              ,"tHelEta/F");
    fNtuple4->Branch("tHelisTOF"            ,&tHelisTOF            ,"tHelisTOF/F");
    fNtuple4->Branch("tHelpoutTPC"          ,&tHelpoutTPC          ,"tHelpoutTPC/F");
    fNtuple4->Branch("tHeltimeTOF"          ,&tHeltimeTOF          ,"tHeltimeTOF/F");
    fNtuple4->Branch("tHeltrackLenghtTOF"   ,&tHeltrackLenghtTOF   ,"tHeltrackLenghtTOF/F");
    fNtuple4->Branch("tHelimpactXY"         ,&tHelimpactXY         ,"tHelimpactXY/F");
    fNtuple4->Branch("tHelimpactZ"          ,&tHelimpactZ          ,"tHelimpactZ/F");
    fNtuple4->Branch("tHelmapITS"           ,&tHelmapITS           ,"tHelmapITS/F");
    fNtuple4->Branch("tHelTPCNcls"          ,&tHelTPCNcls          ,"tHelTPCNcls/F");
    fNtuple4->Branch("tHelTRDsignal"        ,&tHelTRDsignal        ,"tHelTRDsignal/F");
    fNtuple4->Branch("tHelxPrimaryVertex"   ,&tHelxPrimaryVertex   ,"tHelxPrimaryVertex/F");
    fNtuple4->Branch("tHelyPrimaryVertex"   ,&tHelyPrimaryVertex   ,"tHelyPrimaryVertex/F");
    fNtuple4->Branch("tHelzPrimaryVertex"   ,&tHelzPrimaryVertex   ,"tHelzPrimaryVertex/F");
    fNtuple4->Branch("tHelchi2PerClusterTPC",&tHelchi2PerClusterTPC,"tHelchi2PerClusterTPC/F");
    

  } 

  PostData(1,  fListHist);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple4);
}// end UserCreateOutputObjects


//====================== USER EXEC ========================

void AliAnalysisTaskHelium3Pi::UserExec(Option_t *) 
{
  //_______________________________________________________________________
  
  //!*********************!//
  //!  Define variables   !//
  //!*********************!//

  Double_t pinTPC=0.,poutTPC=0.,TPCSignal=0.;
  Double_t xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;
  Double_t massTOF=0.,timeTOF=0.,trackLenghtTOF=0.,betaTOF=0.;

  ULong_t  status=0;
  //  ULong_t  statusT=0;
  ULong_t  statusPi=0;

  Bool_t   isTPC=kFALSE,isTOF=kFALSE,IsHeITSRefit=kFALSE,IsPiITSRefit=kFALSE ;

  Float_t nSigmaNegPion=0.;

  Double_t cutNSigma = 3;
  Double_t bbtheoM=0.,bbtheo=0.;
  Double_t zNathashaNeg=0;
  Double_t zNathashaPos=0;
  Double_t fPos[3]={0.,0.,0.};
  Double_t runNumber=0.;
  //  Double_t evNumber=0.;

  Double_t BCNumber=0.;
  Double_t OrbitNumber=0.;
  Double_t PeriodNumber=0.;

  Double_t        Helium3Mass = 2.80839; 
  Double_t        PionMass    = 0.13957; 
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

  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
    
  Info("AliAnalysisTaskHelium3Pi","Starting UserExec");  

  SetDataType("REAL");
  
  // create pointer to event
  AliESDEvent* lESDevent = dynamic_cast<AliESDEvent*>(event);
  if (!lESDevent) {
    AliError("Cannot get the ESD event");
    return;
  }  

  fHistEventMultiplicity->Fill(0);

  Double_t lMagneticField=lESDevent->GetMagneticField();
  Int_t TrackNumber = -1;

   
  //*****************//  
  //*   Centrality  *//
  //*****************//
 
  AliCentrality *centrality = lESDevent->GetCentrality();
  Float_t percentile=centrality->GetCentralityPercentile("V0M");

  TrackNumber = lESDevent->GetNumberOfTracks();
  if (TrackNumber<2) return;  

  fHistTrackMultiplicity->Fill(TrackNumber,percentile); //tracce per evento

  //****************************************
  
  // PID
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse(); // data member di tipo "const AliPIDResponse *fPIDResponse;"
  //  cout<<"fPIDResponse "<<fPIDResponse<<endl;
  //===========================================

  Int_t eventtype=-99;
  
  Bool_t isSelectedCentral     = (inputHandler->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (inputHandler->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (inputHandler->IsEventSelected() & AliVEvent::kMB);
 
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
  
  if(isSelectedCentral || isSelectedSemiCentral || isSelectedMB){
    
    // ANALISYS
    
    // Primary vertex cut
    
    const AliESDVertex *vtx = lESDevent->GetPrimaryVertexTracks();
    
    if(vtx->GetNContributors()<1) {
      
      // SPD vertex cut
      vtx = lESDevent->GetPrimaryVertexSPD();
      
      if(vtx->GetNContributors()<1) {
	Info("AliAnalysisTaskHelium3Pi","No good vertex, skip event");
	return; // NO GOOD VERTEX, SKIP EVENT 
      }
    }
    
    fHistEventMultiplicity->Fill(1); // analyzed events with PV
 
    xPrimaryVertex=vtx->GetXv();
    yPrimaryVertex=vtx->GetYv();
    zPrimaryVertex=vtx->GetZv();  
    
    if(TMath::Abs(zPrimaryVertex)>10) return;
    
    if(eventtype==1){
      fHistTrackMultiplicityPVCent->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(6); 
    }
    
    if(eventtype==2){
      fHistTrackMultiplicityPVSemiCent->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(7); 
    }
    
    if(eventtype==3){
      fHistTrackMultiplicityPVMB->Fill(TrackNumber,percentile); 
      fHistEventMultiplicity->Fill(8); 
    }
    
    
    fHistEventMultiplicity->Fill(2);
    
    //Find Pair candidates
    
    TArrayI PionsTPC(TrackNumber);        //Neg pions
    Int_t nPionsTPC=0;
    
    TArrayI HeTPC(TrackNumber);        //helium3
    Int_t nHeTPC=0;
    
    const Double_t speedOfLight =  TMath::C()*1E2*1E-12; // cm/ps
    
    Float_t impactXY=-999, impactZ=-999;
    Float_t impactXYpi=-999, impactZpi=-999;

    
    //*************************************************************
    
    runNumber = lESDevent->GetRunNumber();
    BCNumber    = lESDevent->GetBunchCrossNumber();
    OrbitNumber = lESDevent->GetOrbitNumber();
    PeriodNumber= lESDevent->GetPeriodNumber();
    
    //*************************************************************

    for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
      
      AliESDtrack *esdtrack=lESDevent->GetTrack(j);
      //      AliVTrack*  esdtrack= (AliVTrack *) fEvent->GetTrack(iT);

     
      if(!esdtrack) { 
	AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
	continue; 
      }

      // ************** Track cuts ****************
      
      if (!fESDtrackCuts->AcceptTrack(esdtrack)) continue;

      
      status  = (ULong_t)esdtrack->GetStatus();
      isTPC   = (((status) & AliESDtrack::kTPCin)  != 0);
      isTOF   = ((((status) & AliESDtrack::kTOFout) != 0) && (((status) & AliESDtrack::kTIME) != 0));

      
      UInt_t mapITS=esdtrack->GetITSClusterMap();
            
      //----------------------------------------------
      
      //****** Cuts from  AliV0Vertex.cxx *************
      
      Double_t d=esdtrack->GetD(xPrimaryVertex,yPrimaryVertex,lMagneticField);
      //    if (TMath::Abs(d)<fgDPmin) continue;
      if (TMath::Abs(d)>fgRmax) continue;
      
      //---- (Usefull) Stuff
      
      TPCSignal=esdtrack->GetTPCsignal(); 
      
      if (TPCSignal<10)continue;
      if (TPCSignal>1000)continue;
      
      if(!isTPC)continue;
      if(!esdtrack->GetTPCInnerParam())continue;
      
      AliExternalTrackParam trackIn(*esdtrack->GetInnerParam()); 
      pinTPC= trackIn.GetP(); 
      
      //pinTPC= esdtrack->GetTPCMomentum();

      poutTPC=pinTPC;
      
      
      if((status) & (AliESDtrack::kITSrefit!=0)){
	fhBB->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
      }
      
      timeTOF=esdtrack->GetTOFsignal();                 // ps
      trackLenghtTOF= esdtrack->GetIntegratedLength();  // cm
      
      if(isTOF){
	
	if(!esdtrack->GetOuterParam())continue;    
	
	AliExternalTrackParam trackOut(*esdtrack->GetOuterParam()); 
	
	poutTPC = trackOut.GetP(); 
	
	betaTOF= (trackLenghtTOF/timeTOF)/2.99792458e-2;
	
	fhTOF->Fill(poutTPC*esdtrack->GetSign(),betaTOF);
	
	Double_t mass2=(poutTPC*poutTPC)*((((speedOfLight*speedOfLight)*(timeTOF*timeTOF))-(trackLenghtTOF*trackLenghtTOF))/(trackLenghtTOF*trackLenghtTOF));
	if(mass2>0) massTOF=TMath::Sqrt(mass2);
	fhMassTOF->Fill(massTOF);
	
	if(esdtrack->GetSign() < 0.)fBetavsTPCsignalNeg->Fill(betaTOF,TPCSignal);
	if(esdtrack->GetSign() > 0.)fBetavsTPCsignalPos->Fill(betaTOF,TPCSignal);
	
      }
      
      //pass2
      
      // bbtheo =BetheBloch((2*pinTPC)/3.,2,kTRUE);    //! OK
      // bbtheoM=(1 - 0.08*5)*bbtheo;                  //! OK 
      // bbtheoP=(1 + 0.08*5)*bbtheo;                  //! OK
      

      bbtheo = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)0,esdtrack,(AliPID::EParticleType) 7);
      
      if(esdtrack->GetSign()<0){
	zNathashaNeg=bbtheo;//(TPCSignal-bbtheo)/bbtheo;
	//	cout<<"BBtheo 1 :"<<zNathashaNeg<<endl;
	fhNaNeg->Fill(pinTPC,zNathashaNeg); 
      }
      
      if(esdtrack->GetSign() > 0.){
       	zNathashaPos=bbtheo;//(TPCSignal-bbtheo)/bbtheo;
	fhNaPos->Fill(pinTPC,zNathashaPos); 
      }
      
      nSigmaNegPion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType) 2));
      //2 is pion
      
      if ( (nSigmaNegPion < cutNSigma)){ 
	
	//	cout<<"Nsigma pi: "<<nSigmaNegPion<<endl;
	
	fhBBPions->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
	
	if(pinTPC<3.){
	  PionsTPC[nPionsTPC++]=j;
	}
      }
    
      //      nSigmaNegPion=(fPIDResponse->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType) 2));
      
      bbtheoM = TMath::Abs((fPIDResponse->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType) 7)));
      
      //      if( TPCSignal > bbtheoM ) {
      //      if( bbtheoM > -3.) {
      if( bbtheoM < 3.) {
	
	if(pinTPC>0.6){
	  
	  fhBBHe->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
	  HeTPC[nHeTPC++]=j;
	  
	  Bool_t isHeITSrefit=((status) & (AliESDtrack::kITSrefit));
	  
	  esdtrack->GetImpactParameters(impactXY, impactZ);
	  
	  Int_t  fIdxInt[200]; //dummy array
	  Int_t nClustersTPC = esdtrack->GetTPCclusters(fIdxInt);
	  
	  Float_t chi2PerClusterTPC = esdtrack->GetTPCchi2()/(Float_t)(nClustersTPC);
	  
	  tHelrunNumber	        =(Float_t)runNumber;
	  tHelBCNumber	        =(Float_t)BCNumber;
	  tHelOrbitNumber       =(Float_t)OrbitNumber;
	  tHelPeriodNumber      =(Float_t)PeriodNumber;
	  tHeleventtype	        =(Float_t)eventtype;
	  tHelisHeITSrefit      =(Float_t)isHeITSrefit;
	  tHelpercentile        =(Float_t)percentile;
	  tHelSign	        =(Float_t)esdtrack->GetSign();
	  tHelpinTPC	        =(Float_t)pinTPC;
	  tHelGetTPCsignal      =(Float_t)esdtrack->GetTPCsignal();
	  tHelPx	        =(Float_t)esdtrack->Px();
	  tHelPy	        =(Float_t)esdtrack->Py();
	  tHelPz	        =(Float_t)esdtrack->Pz();
	  tHelEta	        =(Float_t)esdtrack->Eta();
	  tHelisTOF	        =(Float_t)isTOF;
	  tHelpoutTPC	        =(Float_t)poutTPC;
	  tHeltimeTOF	        =(Float_t)timeTOF;
	  tHeltrackLenghtTOF    =(Float_t)trackLenghtTOF;
	  tHelimpactXY	        =(Float_t)impactXY;
	  tHelimpactZ	        =(Float_t)impactZ;
	  tHelmapITS	        =(Float_t)mapITS;
	  tHelTPCNcls	        =(Float_t)esdtrack->GetTPCNcls();
	  tHelTRDsignal	        =(Float_t)esdtrack->GetTRDsignal();
	  tHelxPrimaryVertex    =(Float_t)xPrimaryVertex;
	  tHelyPrimaryVertex    =(Float_t)yPrimaryVertex;
	  tHelzPrimaryVertex    =(Float_t)zPrimaryVertex;
	  tHelchi2PerClusterTPC =(Float_t)chi2PerClusterTPC;            
	  	  
	  fNtuple4->Fill();
	}
      }
    }  //! track
	  
    PionsTPC.Set(nPionsTPC);
    HeTPC.Set(nHeTPC);
    
    Double_t        DcaHeToPrimVertex=0;
    Double_t        DcaPionToPrimVertex=0;
    
    impactXY=-999, impactZ=-999;
    impactXYpi=-999, impactZpi=-999;
    
    // Track 
    
    AliESDtrack  *PionTrack = 0x0;
    AliESDtrack  *HeTrack = 0x0;
    
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
    
    //---------------   LOOP PAIRS   ----------------
    
    for (Int_t k=0; k < nPionsTPC; k++) {                           //! Pions Loop
      
      DcaPionToPrimVertex=0.;
      DcaHeToPrimVertex=0;
      
      Int_t PionIdx=PionsTPC[k];
      
      PionTrack=lESDevent->GetTrack(PionIdx);
      
      statusPi = (ULong_t)PionTrack->GetStatus();
      //      isTOFPi  = ((((statusPi) & (AliESDtrack::kTOFout)) != 0) && (((statusPi) & (AliESDtrack::kTIME)) != 0));
      IsPiITSRefit = ((statusPi) & (AliESDtrack::kITSrefit)); 
      
      if (PionTrack) 
	DcaPionToPrimVertex = TMath::Abs(PionTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField)); //OK
      
      if(DcaPionToPrimVertex<0.2)continue; 
      
      AliExternalTrackParam trackInPion(*PionTrack);  
      
      for (Int_t i=0; i<nHeTPC; i++){                               //! Helium Loop
	
	Int_t HeIdx=HeTPC[i];
	
	HeTrack=lESDevent->GetTrack(HeIdx);
	
	//	statusT= (ULong_t)HeTrack->GetStatus();
	//	isTOFHe   = (((statusT & AliESDtrack::kTOFout) != 0) && ((statusT & AliESDtrack::kTIME) != 0));
	IsHeITSRefit = (status & AliESDtrack::kITSrefit); 
	
	if (HeTrack) 
	  DcaHeToPrimVertex = TMath::Abs(HeTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField)); //OK
	
	AliExternalTrackParam trackInHe(*HeTrack); 
    
	if ( DcaPionToPrimVertex < fgDNmin)                //OK
	  if ( DcaHeToPrimVertex < fgDNmin) continue;    //OK
	
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
	
	Float_t CosPointingAngle=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex); //PointingAngle
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
	
	if(vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex)>3) continue;
	
	//salvo solo fino a 3.1 GeV/c2
	
	vHelium.SetXYZM(2*momHeVettAt[0],2*momHeVettAt[1],2*momHeVettAt[2],Helium3Mass); 
	vPion.SetXYZM(momPionVettAt[0],momPionVettAt[1],momPionVettAt[2],PionMass);       
	vSum=vHelium+vPion;
	
	if(vSum.M()>3.2)
	  continue;

	Int_t  fIdxInt[200]; //dummy array
	Int_t nClustersTPCHe = HeTrack->GetTPCclusters(fIdxInt);
	Int_t nClustersTPCPi = PionTrack->GetTPCclusters(fIdxInt);
	
	//----------------------------------------------------------------------//

	trunNumber		=(Float_t)runNumber;
	tbunchcross		=(Float_t)BCNumber;
	torbit		   	=(Float_t)OrbitNumber;
	tperiod		   	=(Float_t)PeriodNumber;
	teventtype		=(Float_t)eventtype;
	tTrackNumber		=(Float_t)TrackNumber;
	tpercentile		=(Float_t)percentile;
	txPrimaryVertex	   	=(Float_t)xPrimaryVertex; //PRIMARY
	tyPrimaryVertex	   	=(Float_t)yPrimaryVertex;
	tzPrimaryVertex	   	=(Float_t)zPrimaryVertex;
	txSecondaryVertex	=(Float_t)fPos[0]; //SECONDARY
	tySecondaryVertex	=(Float_t)fPos[1];
	tzSecondaryVertex	=(Float_t)fPos[2];
	tdcaTracks		=(Float_t)dca;           //between 2 tracks
	tCosPointingAngle	=(Float_t)CosPointingAngle;          //cosPointingAngle da V0
	tDCAV0toPrimaryVertex	=(Float_t)vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
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
	tchi2He		   	=(Float_t)HeTrack->GetTPCchi2()/(Float_t)(nClustersTPCHe);
	tchi2Pi                 =(Float_t)PionTrack->GetTPCchi2()/(Float_t)(nClustersTPCPi);
		

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


