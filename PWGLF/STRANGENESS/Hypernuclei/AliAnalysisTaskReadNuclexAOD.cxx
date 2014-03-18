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
//                 AliAnalysisTaskNUclexAOD class AOD
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
#include "AliAnalysisTaskReadNuclexAOD.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "TString.h"
#include <TDatime.h>
#include <TRandom3.h>
#include <AliAODPid.h>
#include <AliAODHandler.h>
#include <AliAODRecoDecay.h>
#include <AliAODRecoDecayLF.h>
#include <AliAODRecoDecayLF2Prong.h>

ClassImp(AliAnalysisTaskReadNuclexAOD)

//________________________________________________________________________
AliAnalysisTaskReadNuclexAOD::AliAnalysisTaskReadNuclexAOD() 
  : AliAnalysisTaskSE(),
  fAnalysisType("AOD"), 
  fCollidingSystems(0), 
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
  fhBBPions(0),
  fhBBHe(0),
  aodTree(0x0),
  Nuclei(0x0),
  SecondaryVertices(0x0),
  DaughterTracks(0x0),
  nevent(0x0),
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
  t3LHlrunNumber(0),
  t3LHlBCNumber(0),
  t3LHlOrbitNumber(0),
  t3LHlPeriodNumber(0),
  t3LHleventtype(0),
  t3LHlpercentile(0),
  t3LHlPx(0),
  t3LHlPy(0),
  t3LHlPz(0),
  t3LHlEta(0),
  t3LHlY(0),
  t3LHlM(0),
  t3LHlxPrimaryVertex(0),
  t3LHlyPrimaryVertex(0),
  t3LHlzPrimaryVertex(0),
  t3LHlp0px(0),    
  t3LHlp0py(0),
  t3LHlp0pz(0),
  t3LHlp0ch(0),
  t3LHlp1px(0),
  t3LHlp1py(0),
  t3LHlp1pz(0),
  t3LHlp1ch(0)

{
  // Dummy Constructor 
}

//________________________________________________________________________
AliAnalysisTaskReadNuclexAOD::AliAnalysisTaskReadNuclexAOD(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAnalysisType("AOD"), 
    fCollidingSystems(0), 
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
    fhBBPions(0),
    fhBBHe(0),
    aodTree(0x0),
    Nuclei(0x0),
    SecondaryVertices(0x0),
    DaughterTracks(0x0),
    nevent(0x0),
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
    t3LHlrunNumber(0),
    t3LHlBCNumber(0),
    t3LHlOrbitNumber(0),
    t3LHlPeriodNumber(0),
    t3LHleventtype(0),
    t3LHlpercentile(0),
    t3LHlPx(0),
    t3LHlPy(0),
    t3LHlPz(0),
    t3LHlEta(0),
    t3LHlY(0),
    t3LHlM(0),
    t3LHlxPrimaryVertex(0),
    t3LHlyPrimaryVertex(0),
    t3LHlzPrimaryVertex(0),
    t3LHlp0px(0),    
    t3LHlp0py(0),
    t3LHlp0pz(0),
    t3LHlp0ch(0),
    t3LHlp1px(0),
    t3LHlp1py(0),
    t3LHlp1pz(0),
    t3LHlp1ch(0)

{
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
  
  DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class()); 
  DefineOutput(3, TTree::Class()); 

}
//_______________________________________________________
AliAnalysisTaskReadNuclexAOD::~AliAnalysisTaskReadNuclexAOD() 
{ 
  // Destructor
  if (fListHist) {
    delete fListHist;
    fListHist = 0;
  }
  
  if(fNtuple1) delete fNtuple1;
  if(fNtuple4) delete fNtuple4;
}
//=================DEFINITION BETHE BLOCH==============================

Double_t AliAnalysisTaskReadNuclexAOD::BetheBloch(Double_t betaGamma,Double_t charge,Bool_t isPbPb) {

  Double_t kp1, kp2, kp3, kp4, kp5;
  
  if(isPbPb){


 
    //pass2 2011
    kp1 = 4.7*charge*charge;
    kp2 = 8.98482806165147636e+00;
    kp3 = 1.54000000000000005e-05;
    kp4 = 2.30445734159456084e+00;
    kp5 = 2.25624744086878559e+00;

  }
  
  else{


   //pass2 2011
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

void AliAnalysisTaskReadNuclexAOD::Init() 
{
  //  cout<<"----------------------------------> Vediamo"<<endl;
  
  // AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  // cout<<aodHandler<<endl;
  // TTree *aodTree = aodHandler->GetTree();
  // //  aodTree = aodHandler->GetTree();
  // cout<<"aodTree: "<<aodTree<<endl;
  // aodTree->SetBranchAddress("Nuclei",           &Nuclei);
  // aodTree->SetBranchAddress("SecondaryVertices",&SecondaryVertices);
  // aodTree->SetBranchAddress("DaughterTracks",   &DaughterTracks);


}

//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskReadNuclexAOD::UserCreateOutputObjects()
{
  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!

  
  if(! fHistEventMultiplicity ){
    
    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 13 , 0, 13);
     
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Central Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"SemiCentral Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"MB Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"Central Events  w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"SemiCentral Events  w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"MB Events   w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"Sec. Vert. Tree");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"Sec. Vert. Rec");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"Event Not taken");
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

    fNtuple4->Branch("t3LHlrunNumber"        ,&t3LHlrunNumber        ,"t3LHlrunNumber/F");
    fNtuple4->Branch("t3LHlBCNumber"         ,&t3LHlBCNumber         ,"t3LHlBCNumber/F");
    fNtuple4->Branch("t3LHlOrbitNumber"      ,&t3LHlOrbitNumber      ,"t3LHlOrbitNumber/F");
    fNtuple4->Branch("t3LHlPeriodNumber"     ,&t3LHlPeriodNumber     ,"t3LHlPeriodNumber/F");
    fNtuple4->Branch("t3LHleventtype"        ,&t3LHleventtype        ,"t3LHleventtype/F");
    fNtuple4->Branch("t3LHlpercentile"       ,&t3LHlpercentile       ,"t3LHlpercentile/F");
    fNtuple4->Branch("t3LHlPx"               ,&t3LHlPx               ,"t3LHlPx/F");
    fNtuple4->Branch("t3LHlPy"               ,&t3LHlPy               ,"t3LHlPy/F");
    fNtuple4->Branch("t3LHlPz"               ,&t3LHlPz               ,"t3LHlPz/F");
    fNtuple4->Branch("t3LHlEta"              ,&t3LHlEta              ,"t3LHlEta/F");
    fNtuple4->Branch("t3LHlY"                ,&t3LHlY                ,"t3LHlY/F");
    fNtuple4->Branch("t3LHlM"                ,&t3LHlM                ,"t3LHlM/F");
    fNtuple4->Branch("t3LHlxPrimaryVertex"   ,&t3LHlxPrimaryVertex   ,"t3LHlxPrimaryVertex/F");
    fNtuple4->Branch("t3LHlyPrimaryVertex"   ,&t3LHlyPrimaryVertex   ,"t3LHlyPrimaryVertex/F");
    fNtuple4->Branch("t3LHlzPrimaryVertex"   ,&t3LHlzPrimaryVertex   ,"t3LHlzPrimaryVertex/F");
    fNtuple4->Branch("t3LHlp0px"             ,&t3LHlp0px             ,"t3LHlp0px/F");
    fNtuple4->Branch("t3LHlp0py"             ,&t3LHlp0py             ,"t3LHlp0py/F");
    fNtuple4->Branch("t3LHlp0pz"             ,&t3LHlp0pz             ,"t3LHlp0pz/F");
    fNtuple4->Branch("t3LHlp0ch"             ,&t3LHlp0ch             ,"t3LHlp0ch/F");
    fNtuple4->Branch("t3LHlp1px"             ,&t3LHlp1px             ,"t3LHlp1px/F");
    fNtuple4->Branch("t3LHlp1py"             ,&t3LHlp1py             ,"t3LHlp1py/F");
    fNtuple4->Branch("t3LHlp1pz"             ,&t3LHlp1pz             ,"t3LHlp1pz/F");
    fNtuple4->Branch("t3LHlp1ch"             ,&t3LHlp1ch             ,"t3LHlp1ch/F");


  } 

  
  PostData(1,  fListHist);
  PostData(2,  fNtuple1);
  PostData(3,  fNtuple4);
}// end UserCreateOutputObjects



//====================== USER EXEC ========================

void AliAnalysisTaskReadNuclexAOD::UserExec(Option_t *) 
{
  //_______________________________________________________________________

  //!*********************!//
  //!  Define variables   !//
  //!*********************!//

  Double_t  xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;

  ULong_t  statusT;
  ULong_t  statusPi;
  
  Bool_t   IsHeITSRefit=kFALSE,IsPiITSRefit=kFALSE ;

  Double_t fPos[3]={0.,0.,0.};

  Double_t runNumber=0.;
  Double_t BCNumber=0.;
  Double_t OrbitNumber=0.;
  Double_t PeriodNumber=0.;

  Double_t        Helium3Mass = 2.80839; 
  Double_t        PionMass    = 0.13957; 

  // TLORENTZ vectors
  
  TLorentzVector  vPion,vHelium,vSum;
  TLorentzVector  vHyp;
  TLorentzVector  vp0,vp1,vS;
  //!----------------------------------------------------------------
  //! A set of very loose parameters for cuts 
  
  //  Double_t fgChi2max=33.;     //! max chi2
  Double_t fgDNmin=0.05;      //! min imp parameter for the 1st daughter = 500um
  //  Double_t fgDPmin=0.05;      //! min imp parameter for the 2nd daughter = 500um
  Double_t fgDCAmax=1.;       //! max DCA between the daughter tracks in cm
  //  Double_t fgCPAmin=0.99;      //! min cosine of V0's pointing angle
  Double_t fgRmin=0.1;        //! min radius of the fiducial volume = 1 mm 
  Double_t fgRmax=200.;       //! max radius of the fiducial volume = 2 m
  
  //!-----------------------------------------------------------------

  // Main loop
  // Called for EACH event
  
  nevent++;  //--> This is important

  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  Info("AliAnalysisTaskReadNuclexAOD","Starting UserExec");  

  SetDataType("REAL");
  
  AliAODEvent *lAODevent=(AliAODEvent *)InputEvent();	
  if (!lAODevent) {
    Printf("ERROR: aod not available");
    //    fHistLog->Fill(98);
    return;
  }
 
  fHistEventMultiplicity->Fill(0.5);

  //*************************************************************
  
  runNumber    = lAODevent->GetRunNumber();
  BCNumber     = lAODevent->GetBunchCrossNumber();
  OrbitNumber  = lAODevent->GetOrbitNumber();
  PeriodNumber = lAODevent->GetPeriodNumber();
   
  AliCentrality *centrality = lAODevent->GetCentrality();
  centrality = lAODevent->GetCentrality();
    
  Float_t percentile=centrality->GetCentralityPercentile("V0M");
  Float_t refMult = lAODevent->GetHeader()->GetRefMultiplicity();
  
  fHistTrackMultiplicity->Fill(refMult,percentile); //tracce per evento

  //  cout<<"Centrality: "<<percentile<<" number of track: "<<lAODevent->GetHeader()->GetRefMultiplicity()<<endl;
  //*************************************************************
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  
  Int_t eventtype=-99;
  
  Bool_t isSelectedCentral     = (inputHandler->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (inputHandler->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (inputHandler->IsEventSelected() & AliVEvent::kMB);
 
  if(isSelectedCentral){
    fHistEventMultiplicity->Fill(3.5);
    fHistTrackMultiplicityCent->Fill(refMult,percentile); 
    eventtype=1;
  }

  if(isSelectedSemiCentral){
    fHistEventMultiplicity->Fill(4.5);
    fHistTrackMultiplicitySemiCent->Fill(refMult,percentile); 
    eventtype=2;
  }

  if(isSelectedMB){
    fHistEventMultiplicity->Fill(5.5);
    fHistTrackMultiplicityMB->Fill(refMult,percentile); 
    eventtype=3;
  }

  if(eventtype==-99)return;
  // cout<<"eventtype: "<<eventtype<<endl;
  
  // if(!isSelectedCentral)return;
  // if(!isSelectedSemiCentral)return;
  // if(!isSelectedMB)return;
  
  cout<<"eventtype: "<<eventtype<<endl;
    
  //-------------------------------
  
  Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
  Double_t lMagneticField                 = -10.; 
  
  AliAODVertex* primvtx = lAODevent->GetPrimaryVertex();
  if(!primvtx){
    cout<<"No PV"<<endl;
    return;
  }
  fHistEventMultiplicity->Fill(1.5);
  primvtx->GetXYZ( lBestPrimaryVtxPos );

  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
    return; 
  }
  fHistEventMultiplicity->Fill(2.5);


  if(eventtype==1){
    fHistTrackMultiplicityPVCent->Fill(refMult,percentile); 
    fHistEventMultiplicity->Fill(6.5); 
  }
  
  if(eventtype==2){
    fHistTrackMultiplicityPVSemiCent->Fill(refMult,percentile); 
    fHistEventMultiplicity->Fill(7.5); 
  }
  
  if(eventtype==3){
    fHistTrackMultiplicityPVMB->Fill(refMult,percentile); 
    fHistEventMultiplicity->Fill(8.5); 
  }
  
  lMagneticField = lAODevent->GetMagneticField();

  xPrimaryVertex = lBestPrimaryVtxPos[0];
  yPrimaryVertex = lBestPrimaryVtxPos[1];
  zPrimaryVertex = lBestPrimaryVtxPos[2];

  //--------------------------------------------------------

  AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

  aodTree = aodHandler->GetTree();
  aodTree->SetBranchAddress("Nuclei",&Nuclei);
  aodTree->SetBranchAddress("SecondaryVertices",&SecondaryVertices);
  aodTree->SetBranchAddress("DaughterTracks",&DaughterTracks);

  //---------------------------------------------------------

  aodTree->GetEntry(nevent);

  if(Nuclei-> GetEntriesFast()>0)
    fHistEventMultiplicity->Fill(2);

  if(SecondaryVertices-> GetEntriesFast()>0)
    fHistEventMultiplicity->Fill(3);
    
  for(Int_t i=0; i < SecondaryVertices-> GetEntriesFast() ; i++){

    AliAODRecoDecayLF2Prong *d = (AliAODRecoDecayLF2Prong*)SecondaryVertices->UncheckedAt(i);
    //cout<<"\n\n\nPT: "<<d->Pt()<<endl<<endl<<endl<<endl;
    if(!d) continue;
   
    fHistEventMultiplicity->Fill(10.5);
    
    //    Double_t tDecayVertexV0[3]; 
    // if(d->GetXYZ(tDecayVertexV0))
    
    Double_t tV0mom[3];
    d->PxPyPz(tV0mom); 
    
    // TLORENTZ vectors

    //    cout<<"px prong 1: "<<d->PxProng(0)<<endl;

    vHyp.SetXYZM(tV0mom[0],tV0mom[1],tV0mom[2],2.991);
    

    vp0.SetXYZM(d->PxProng(0),d->PyProng(0),d->PzProng(0),PionMass);
    vp1.SetXYZM(d->PxProng(1),d->PyProng(1),d->PzProng(1),Helium3Mass);
    vS = vp0 + vp1;

    //    cout<<"Momemnto pt: "<<lPt<<endl;
    
    t3LHlrunNumber	   =(Float_t)runNumber;
    t3LHlBCNumber	   =(Float_t)BCNumber;
    t3LHlOrbitNumber       =(Float_t)OrbitNumber;
    t3LHlPeriodNumber      =(Float_t)PeriodNumber;
    t3LHleventtype	   =(Float_t)eventtype;
    t3LHlpercentile        =(Float_t)percentile;
    t3LHlPx	           =(Float_t)tV0mom[0];
    t3LHlPy	           =(Float_t)tV0mom[1];
    t3LHlPz	           =(Float_t)tV0mom[2];
    t3LHlEta	           =(Float_t)d->Eta();
    t3LHlY	           =(Float_t)vHyp.Rapidity();
    t3LHlM	           =(Float_t)vS.M();
    t3LHlxPrimaryVertex    =(Float_t)xPrimaryVertex;
    t3LHlyPrimaryVertex    =(Float_t)yPrimaryVertex;
    t3LHlzPrimaryVertex    =(Float_t)zPrimaryVertex;

    t3LHlp0px    =(Float_t)d->PxProng(0);
    t3LHlp0py    =(Float_t)d->PyProng(0);
    t3LHlp0pz    =(Float_t)d->PzProng(0);
    t3LHlp0ch    =(Float_t)d->ChargeProng(0);
    t3LHlp1px    =(Float_t)d->PxProng(1);
    t3LHlp1py    =(Float_t)d->PyProng(1);
    t3LHlp1pz    =(Float_t)d->PzProng(1);
    t3LHlp1ch    =(Float_t)d->ChargeProng(1);

    fNtuple4->Fill();
    
  }
  
  
  //----------------------------------------
  // Loop on tracks like-AnalysisTaskHelium3Pi.cxx
  
  Int_t TrackNumber= DaughterTracks-> GetEntriesFast();

  TArrayI Track0(TrackNumber);        //Pions                                                                          
  Int_t nTrack0=0;
 
  TArrayI Track1(TrackNumber);        //Helium3
  Int_t nTrack1=0;
 

  Float_t impactXY=-999, impactZ=-999;
  Float_t impactXYpi=-999, impactZpi=-999;
  
  for(Int_t j=0; j<TrackNumber; j++){
    //  cout<<"eventtype inside: "<<eventtype<<endl;
    //    cout<<"Inside loop tracks"<<endl;
    AliVTrack *track = (AliVTrack*)DaughterTracks->UncheckedAt(j);
    AliAODTrack *aodtrack =(AliAODTrack*)track;
    
    fhBB->Fill(aodtrack->GetDetPid()->GetTPCmomentum()*aodtrack->Charge(),aodtrack->GetTPCsignal());
    
    //-----------------------------------------------------------
    //Track cuts
    // if(aodtrack->GetTPCNcls() < 80 )continue;
    // if(aodtrack->Chi2perNDF() > 5 )continue;

    // if (!aodtrack->IsOn(AliAODTrack::kTPCrefit)) continue;
    // if (!aodtrack->IsOn(AliAODTrack::kTPCin)) continue;
    // if (aodtrack->IsOn(AliAODTrack::kITSpureSA)) continue;

    //---------------------------------------------------------------
     
    if (aodtrack->GetTPCsignal()<100){
      Track0[nTrack0++]=j;
      fhBBPions->Fill(aodtrack->GetDetPid()->GetTPCmomentum()*aodtrack->Charge(),aodtrack->GetTPCsignal());
    }

    else{
      Track1[nTrack1++]=j;
      fhBBHe->Fill(aodtrack->GetDetPid()->GetTPCmomentum()*aodtrack->Charge(),aodtrack->GetTPCsignal());
    }
 
  }

  Track0.Set(nTrack0);
  Track1.Set(nTrack1);
 

  //cout<<"nTrack0: "<<nTrack0<<endl;
  //cout<<"nTrack1: "<<nTrack1<<endl;
  
  
  //cout<<"Track loops..."<<endl;
   
  Double_t        DcaHeToPrimVertex=0;
  Double_t        DcaPionToPrimVertex=0;
  
  impactXY=-999, impactZ=-999;
  impactXYpi=-999, impactZpi=-999;
  
 
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
  

  //Loop on the tracks

  for (Int_t k=0; k <nTrack0 ; k++) {                           //! Pion Tracks Loop

    DcaPionToPrimVertex=0.;
    DcaHeToPrimVertex=0;
    
    Int_t Track0idx=Track0[k];
    
    AliVTrack *trackpi = (AliVTrack*)DaughterTracks->UncheckedAt(Track0idx);
    PionTrack = new AliESDtrack(trackpi); 
    
    statusPi = (ULong_t)trackpi->GetStatus();
    IsPiITSRefit = ((statusPi) & (AliESDtrack::kITSrefit)); 
    
    if (PionTrack) 
      DcaPionToPrimVertex = TMath::Abs(PionTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField));

    //cout<<"DcaPionToPrimVertex: "<<DcaPionToPrimVertex<<endl;

    if(DcaPionToPrimVertex<0.2)continue; 
    AliExternalTrackParam trackInPion(*PionTrack);  

    for (Int_t i=0; i<nTrack1; i++){                            //! He Tracks Loop
      
      Int_t Track1idx=Track1[i];
      
      AliVTrack *trackhe = (AliVTrack*)DaughterTracks->UncheckedAt(Track1idx);
      HeTrack = new AliESDtrack(trackhe);
      
      statusT= (ULong_t)HeTrack->GetStatus();
      IsHeITSRefit = (statusT & AliESDtrack::kITSrefit); 

      if (HeTrack) 
	DcaHeToPrimVertex = TMath::Abs(HeTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField)); //OK
      
      AliExternalTrackParam trackInHe(*HeTrack); 
    
      if ( DcaPionToPrimVertex < fgDNmin)                //OK
	if ( DcaHeToPrimVertex < fgDNmin) continue;    //OK
      
      Double_t xn, xp;
      Double_t dca=0.;
      
      dca= PionTrack->GetDCA(HeTrack,lMagneticField,xn,xp); //!dca (Neg to Pos)
      
      //cout<<"dca tracks: "<<dca<<endl;

      if (dca > fgDCAmax) continue;
      if ((xn+xp) > 2*fgRmax) continue;
      if ((xn+xp) < 2*fgRmin) continue;
      
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
	
	AliESDv0 vertex(trackInPion,Track0idx,trackInHe,Track1idx);
	//	if (vertex.GetChi2V0() > fgChi2max) continue;
	
	//cout<<"Made v0"<<endl;

	Float_t CosPointingAngle=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex); //PointingAngle
	//cout<<"cpa: "<<CosPointingAngle<<endl;
	//	if (CosPointingAngle < fgCPAmin) continue;
	
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
	
	//	if(vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex)>3) continue;
	
	//salvo solo fino a 3.1 GeV/c2
	
	vHelium.SetXYZM(2*momHeVettAt[0],2*momHeVettAt[1],2*momHeVettAt[2],Helium3Mass); 
	vPion.SetXYZM(momPionVettAt[0],momPionVettAt[1],momPionVettAt[2],PionMass);       
	vSum=vHelium+vPion;
	
	//	if(vSum.M()>3.2)
	//  continue;
	
	Int_t  fIdxInt[200]; //dummy array
	Int_t nClustersTPCHe = HeTrack->GetTPCclusters(fIdxInt);
	Int_t nClustersTPCPi = PionTrack->GetTPCclusters(fIdxInt);
	
	//cout<<"Cos PA: "<<CosPointingAngle<<endl;
	fHistEventMultiplicity->Fill(11.5);
	//-------------------------------------------------------------
	
	trunNumber		=(Float_t)runNumber;
	tbunchcross		=(Float_t)BCNumber;
	torbit		   	=(Float_t)OrbitNumber;
	tperiod		   	=(Float_t)PeriodNumber;
	teventtype		=(Float_t)eventtype;
	tTrackNumber		=(Float_t)TrackNumber;
	tpercentile		=(Float_t)percentile;
	txPrimaryVertex	   	=(Float_t)xPrimaryVertex;            //PRIMARY
	tyPrimaryVertex	   	=(Float_t)yPrimaryVertex;
	tzPrimaryVertex	   	=(Float_t)zPrimaryVertex;
	txSecondaryVertex	=(Float_t)fPos[0];                   //SECONDARY
	tySecondaryVertex	=(Float_t)fPos[1];
	tzSecondaryVertex	=(Float_t)fPos[2];
	tdcaTracks		=(Float_t)dca;                       //between 2 tracks
	tCosPointingAngle	=(Float_t)CosPointingAngle;          //cosPointingAngle da V0
	tDCAV0toPrimaryVertex	=(Float_t)vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
	tHeSign		   	=(Float_t)HeTrack->GetSign();        //He
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
		
    }    
    
  }
  
  PostData(1,fListHist);
  PostData(2,fNtuple1);
  PostData(3,fNtuple4);
    
} //end userexec


//________________________________________________________________________

void AliAnalysisTaskReadNuclexAOD::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}


