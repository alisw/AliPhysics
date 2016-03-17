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
#include "AliAnalysisTaskHelium3PiMC.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "TString.h"
#include <TDatime.h>
#include <TRandom3.h>
using std::endl;
using std::cout;

const Int_t AliAnalysisTaskHelium3PiMC::fgNrot = 15;


ClassImp(AliAnalysisTaskHelium3PiMC)

//________________________________________________________________________
AliAnalysisTaskHelium3PiMC::AliAnalysisTaskHelium3PiMC() 
: AliAnalysisTaskSE(),
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fDataType("SIM"),
  fListHistCascade(0), 
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fHistMCMultiplicityTracks(0),
  fHistMCEta(0), 
  fHistMCPt(0), 
  fHistMCTheta(0), 
  fHistMCDecayPosition(0),
  fHistMCDecayRho(0), 
  fhRigidityHevsMomPiMC(0),
  fhRigidityHevsMomPiRec(0),
  fhInvMassMC(0),
  fhInvMassMum(0),
  fhInvMassRec(0),
  fhInvMassRec1(0),
  fhInvMassRec2(0), 
  fhInvMassRec3(0),
  fhInvMassRec4(0),
  fhInvMassRec5(0),
  fhInvMassRec6(0),
  fhInvMassRec7(0),
  fhHeMCRigidity(0),
  fhPioneMC(0),
  hBBTPCnoCuts(0),
  fhBBTPC(0),
  fhBBTPCNegativePions(0),
  fhBBTPCPositivePions(0),
  fhBBTPCHe3(0),
  fHistProvaDCA(0),
  fHistPercentileVsTrackNumber(0),
  fhHeDCAXY(0),
  fhHeDCAZ(0),
  fhPiDCAXY(0),
  fhPiDCAZ(0),
  hITSClusterMap(0),
  fNtuple1(0),
  fNtuple2(0)

{
  // Dummy Constructor(0); 
}

//________________________________________________________________________
AliAnalysisTaskHelium3PiMC::AliAnalysisTaskHelium3PiMC(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAnalysisType("ESD"), 
    fCollidingSystems(0), 
    fDataType("SIM"),
    fListHistCascade(0), 
    fHistEventMultiplicity(0), 
    fHistTrackMultiplicity(0),
    fHistMCMultiplicityTracks(0),
    fHistMCEta(0), 
    fHistMCPt(0), 
    fHistMCTheta(0), 
    fHistMCDecayPosition(0),
    fHistMCDecayRho(0), 
    fhRigidityHevsMomPiMC(0),
    fhRigidityHevsMomPiRec(0),
    fhInvMassMC(0),
    fhInvMassMum(0),
    fhInvMassRec(0),
    fhInvMassRec1(0),
    fhInvMassRec2(0), 
    fhInvMassRec3(0),
    fhInvMassRec4(0),
    fhInvMassRec5(0),
    fhInvMassRec6(0),
    fhInvMassRec7(0),
    fhHeMCRigidity(0),
    fhPioneMC(0),
    hBBTPCnoCuts(0),
    fhBBTPC(0),
    fhBBTPCNegativePions(0),
    fhBBTPCPositivePions(0),
    fhBBTPCHe3(0),
    fHistProvaDCA(0),
    fHistPercentileVsTrackNumber(0),
    fhHeDCAXY(0),
    fhHeDCAZ(0),
    fhPiDCAXY(0),
    fhPiDCAZ(0),
    hITSClusterMap(0),
    fNtuple1(0),
    fNtuple2(0)
  

{
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
}
//_______________________________________________________
AliAnalysisTaskHelium3PiMC::~AliAnalysisTaskHelium3PiMC() 
{ 
  // Destructor
  if (fListHistCascade) {
    delete fListHistCascade;
    fListHistCascade = 0;
  }

}
//=================DEFINITION BETHE BLOCH==============================

Double_t AliAnalysisTaskHelium3PiMC::BetheBloch(Double_t betaGamma,Double_t charge,Bool_t isPbPb) {

  Double_t kp1, kp2, kp3, kp4, kp5;
  
  if(isPbPb){

    //pass2 //to be checked
    kp1 = 5.2*charge*charge;
    kp2 = 8.98482806165147636e+00;
    kp3 = 1.54000000000000005e-05;
    kp4 = 2.30445734159456084e+00;
    kp5 = 2.25624744086878559e+00;


  }
  
  else{
    
    //pass2 // to be defined
    kp1 = 5.2*charge*charge;
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

void AliAnalysisTaskHelium3PiMC::UserCreateOutputObjects()
{
  fListHistCascade = new TList();
  fListHistCascade->SetOwner();  // IMPORTANT!

  if(! fHistEventMultiplicity ){
    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 6 , -1, 5 );
    fHistEventMultiplicity->GetXaxis()->SetTitle("Event Type");
    fListHistCascade->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){

    fHistTrackMultiplicity   = new TH1F( "fHistTrackMultiplicity" , "Nb of Tracks" , 25000,0, 25000 );
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fListHistCascade->Add(fHistTrackMultiplicity);
  } 
 
  if(! fHistMCMultiplicityTracks){ 
    fHistMCMultiplicityTracks =new TH1F("fHistMCMultiplicityTracks","MC Multiplicity Tracks",1000,0,1000); 
    fHistMCMultiplicityTracks->GetXaxis()->SetTitle("MC Number of tracks");
    fListHistCascade->Add(fHistMCMultiplicityTracks); 
  }
  if(!fHistMCEta ){ 
    fHistMCEta=new TH1F("fHistMCEta","MC eta",1000,-3,3);                                     
    fHistMCEta->GetXaxis()->SetTitle("Injected Eta");
    fListHistCascade->Add(fHistMCEta);
  } 
  if(! fHistMCPt){ 
    fHistMCPt =new TH1F("fHistMCPt","MC pt",1000,0,20);     
    fHistMCPt->GetXaxis()->SetTitle("Injected Pt");
    fListHistCascade->Add(fHistMCPt); 
  }  
  if(!fHistMCTheta ){ 
    fHistMCTheta=new TH1F("fHistMCTheta","MC theta",1000,-6,6);                                   
    fHistMCTheta->GetXaxis()->SetTitle("Injected Theta");
    fListHistCascade->Add(fHistMCTheta); 
  }
  if(!fHistMCDecayPosition){ 
    fHistMCDecayPosition =new TH1F("fHistMCDecayPosition","MC Decay Position",10000,0,1000);  
    fHistMCDecayPosition->GetXaxis()->SetTitle("Decay Position");
    fListHistCascade->Add(fHistMCDecayPosition); 
  }    
  if(!fHistMCDecayRho ){ 
    fHistMCDecayRho  =new TH1F("fHistMCDecayRho","MC decay position 3d",10000,0,1000);  
    fHistMCDecayRho->GetXaxis()->SetTitle("Decay rho");
    fListHistCascade->Add(fHistMCDecayRho); 
  } 

  if(!fhRigidityHevsMomPiMC ){ 
    fhRigidityHevsMomPiMC=new TH2F("fhRigidityHevsMomPiMC","Rigidity He vs Mom Pi MC",20,0,10,300,0,30);
    fhRigidityHevsMomPiMC->GetXaxis()->SetTitle("He3 Rigidity");
    fhRigidityHevsMomPiMC->GetYaxis()->SetTitle("Pi momentum");
    fListHistCascade->Add(fhRigidityHevsMomPiMC); 
  }

  if(! fhRigidityHevsMomPiRec){ 
    fhRigidityHevsMomPiRec=new TH2F("fhRigidityHevsMomPiRec","Rigidity He vs Mom Pi Rec",20,0,10,300,0,30);
    fhRigidityHevsMomPiRec->GetXaxis()->SetTitle("He3 Rigidity");
    fhRigidityHevsMomPiRec->GetYaxis()->SetTitle("Pi momentum");
    fListHistCascade->Add(fhRigidityHevsMomPiRec); 
  }
  
  if(!fhInvMassMC){
    fhInvMassMC=new TH1F("fhInvMassMC","fhInvMassMC",800,2.,6.);
    fhInvMassMC->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassMC); 
  }
  
  if(!fhInvMassMum){
    fhInvMassMum=new TH1F("fhInvMassMum","fhInvMassMum",800,2.,6.);
    fhInvMassMum->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassMum); 
  }
  
  if(!fhInvMassRec){
    fhInvMassRec=new TH1F("fhInvMassRec","fhInvMassRec",800,2.,6.);
    fhInvMassRec->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec);
  }

  if(!fhInvMassRec1){
    fhInvMassRec1=new TH1F("fhInvMassRec1","No Altri tagli",800,2.,6.);
    fhInvMassRec1->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec1);
  }
  if(!fhInvMassRec2){
    fhInvMassRec2=new TH1F("fhInvMassRec2","DCA pi > 0.1",800,2.,6.);
    fhInvMassRec2->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec2);
  }
  if(!fhInvMassRec3){
    fhInvMassRec3=new TH1F("fhInvMassRec3","DCA He > 0.05",800,2.,6.);
    fhInvMassRec3->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec3);
  }
  if(!fhInvMassRec4){
    fhInvMassRec4=new TH1F("fhInvMassRec4","DCA tracks < 1 cm",800,2.,6.);
    fhInvMassRec4->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec4);
  }
  if(!fhInvMassRec5){
    fhInvMassRec5=new TH1F("fhInvMassRec5","Condizione xn+xp",800,2.,6.);
    fhInvMassRec5->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec5);
  }
  if(!fhInvMassRec6){
    fhInvMassRec6=new TH1F("fhInvMassRec6","Ho fatto V0 ",800,2.,6.);
    fhInvMassRec6->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec6);
  }
  if(!fhInvMassRec7){
    fhInvMassRec7=new TH1F("fhInvMassRec7","V0+Taglio CPA",800,2.,6.);
    fhInvMassRec7->GetXaxis()->SetTitle("(He3,#pi) InvMass");
    fListHistCascade->Add(fhInvMassRec7);
  }

  if(!fhHeMCRigidity ){ 
    fhHeMCRigidity=new TH1F("fhHeMCRigidity","He3 rigidity distribution",200,0,20);
    fhHeMCRigidity->GetXaxis()->SetTitle("He3 rigidity");
    fListHistCascade->Add(fhHeMCRigidity); 
  }
  if(!fhPioneMC ){ 
    fhPioneMC=new TH1F("hPioneMC","Pion mom distribution",200,0,50);
    fhPioneMC->GetXaxis()->SetTitle("Pion momentum");
    fListHistCascade->Add(fhPioneMC); 
  }
  
  if(!hBBTPCnoCuts ){ 
    hBBTPCnoCuts=new TH2F("hBBTPCnoCuts","scatterPlot TPC no cuts",2000,-10,10,1000,0,3000);
    hBBTPCnoCuts->GetXaxis()->SetTitle("p/Z (GeV/#it{c})");
    hBBTPCnoCuts->GetYaxis()->SetTitle("TPC Signal (a.u)");
    fListHistCascade->Add(hBBTPCnoCuts); 
  }
  if(!fhBBTPC ){ 
    fhBBTPC=new TH2F("fhBBTPC","scatterPlot TPC",2000,-10,10,1000,0,3000);
    fhBBTPC->GetXaxis()->SetTitle("p/Z (GeV/#it{c})");
    fhBBTPC->GetYaxis()->SetTitle("TPC Signal (a.u)");
    fListHistCascade->Add(fhBBTPC); 
  }
  if(!fhBBTPCNegativePions ){ 
    fhBBTPCNegativePions=new TH2F("fhBBTPCNegativePions","scatterPlot Neg Pions",2000,-10,10,1000,0,3000);
    fhBBTPCNegativePions->GetXaxis()->SetTitle("p/Z (GeV/#it{c})");
    fhBBTPCNegativePions->GetYaxis()->SetTitle("TPC Signal (a.u)");
    fListHistCascade->Add(fhBBTPCNegativePions); 
  }
  if(!fhBBTPCPositivePions ){ 
    fhBBTPCPositivePions=new TH2F("fhBBTPCPositivePions","scatterPlot Pos Pions",2000,-10,10,1000,0,3000);
    fhBBTPCPositivePions->GetXaxis()->SetTitle("p/Z (GeV/#it{c})");
    fhBBTPCPositivePions->GetYaxis()->SetTitle("TPC Signal (a.u)");
    fListHistCascade->Add(fhBBTPCPositivePions); 
  }
  if(!fhBBTPCHe3 ){ 
    fhBBTPCHe3=new TH2F("fhBBTPCHe3","scatterPlot TPC -  He3",2000,-10,10,1000,0,3000);
    fhBBTPCHe3->GetXaxis()->SetTitle("p/Z (GeV/#it{c})");
    fhBBTPCHe3->GetYaxis()->SetTitle("TPC Signal (a.u)");
    fListHistCascade->Add(fhBBTPCHe3); 
  }
  if(!fHistProvaDCA ){ 
    fHistProvaDCA=new TH2F("fHistProvaDCA","fHistProvaDCA",1000,-50,50,1000,0,100);
    fHistProvaDCA->GetXaxis()->SetTitle("xn+xp");
    fHistProvaDCA->GetYaxis()->SetTitle("dca tracks");
    fListHistCascade->Add(fHistProvaDCA); 
  }
  
  if(!hITSClusterMap ){ 
    hITSClusterMap=new TH1F("hITSClusterMap","hITSClusterMap",65,-1,64);
    fListHistCascade->Add(hITSClusterMap); 
  }

  if(!fHistPercentileVsTrackNumber){
    fHistPercentileVsTrackNumber=new TH2F("fHistPercentileVsTrackNumber","fHistPercentileVsTrackNumber",120,-3,117,2500,0,25000);
    fHistPercentileVsTrackNumber->GetXaxis()->SetTitle("Percentile");
    fHistPercentileVsTrackNumber->GetYaxis()->SetTitle("Tracks Number");
    fListHistCascade->Add(fHistPercentileVsTrackNumber); 
  }

  if(!fhHeDCAXY){
    fhHeDCAXY=new TH1F("fhHeDCAXY","fhHeDCAXY",800,-4,4);
    fListHistCascade->Add(fhHeDCAXY); 
  }
  if(!fhHeDCAZ){
    fhHeDCAZ=new TH1F("fhHeDCAZ","fhHeDCAZ",800,-30,30);
    fListHistCascade->Add(fhHeDCAZ); 
  }
  if(!fhPiDCAXY){
    fhPiDCAXY=new TH1F("fhPiDCAXY","fhPiDCAXY",800,-4,4);
    fListHistCascade->Add(fhPiDCAXY); 
  }
  if(!fhPiDCAZ){
    fhPiDCAZ=new TH1F("fhPiDCAZ","fhPiDCAZ",800,-30,30);
    fListHistCascade->Add(fhPiDCAZ); 
  }

  if(! fNtuple1 ) {
    fNtuple1 = new TNtuple("fNtuple1","Ntuple1","runNumber:evNumber:TrackNumber:percentile:xPrimaryVertex:yPrimaryVertex:zPrimaryVertex:xSecondaryVertex:ySecondaryVertex:zSecondaryVertex:dcaTracks:CosPointingAngle:DCAV0toPrimaryVertex:HeSign:HepInTPC:HeTPCsignal:DcaHeToPrimVertex:HeEta:momHex:momHey:momHez:momHeAtSVx:momHeAtSVy:momHeAtSVz:HeTPCNcls:HeimpactXY:HeimpactZ:isTOFHe:HeBeta:HeITSClusterMap:IsHeITSRefit:PionSign:PionpInTPC:PionTPCsignal:DcaPionToPrimVertex:PionEta:momPionx:momPiony:momPionz:momNegPionAtSVx:momNegPionAtSVy:momNegPionAtSVz:PionTPCNcls:PionimpactXY:PionimpactZ:isTOFPion:PionBeta:PionITSClusterMap:IsPiITSRefit:PDGCodeNeg:PDCCodePos:motherPDGNeg:motherPDGPos:labelPi:labelHe:mumidNeg:mumidPos");
  
    fListHistCascade->Add(fNtuple1);
  }
  
  if(! fNtuple2 ) {
    
    fNtuple2 = new TNtuple("fNtuple2","Ntuple2","run:event:iMC:Centrality:PVx:PVy:PVz:PDGcodeMum:MotherIndex:SVxD0:SVyD0:SVzD0:SVxD1:SVyD1:SVzD1:SV3d:EtaMum:YMum:ThetaMum:PhiMum:PxMum:PyMum:PzMum:PdgDaughter0:PdgDaughter1:PxD0:PyD0:PzD0:PxD1:PyD1:PzD1");
    
    fListHistCascade->Add(fNtuple2);
  } 
    
  PostData(1,fListHistCascade);

}// end UserCreateOutputObjects



//====================== USER EXEC ========================

void AliAnalysisTaskHelium3PiMC::UserExec(Option_t *) 
{
  //_______________________________________________________________________
  
  //!*********************!//
  //!  Define variables   !//
  //!*********************!//
  Float_t vett1[60];
  for(Int_t i=0;i<60;i++) vett1[i]=0;

  Float_t vett2[40];
  for(Int_t i=0;i<40;i++) vett2[i]=0;


  Double_t  pinTPC=0.,TPCSignal=0.;
  Double_t xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;

  ULong_t  status=0;
  ULong_t  statusT=0;
  ULong_t  statusPi=0;

  Bool_t   isTPC=kFALSE,isTOFHe3=kFALSE,isTOFPi=kFALSE;

  Double_t fPos[3]={0.,0.,0.};
  Double_t runNumber=0.;
  Double_t evNumber=0.;
 
  Int_t id0           = 0, id1          = 0;
  Double_t mcDecayPosXD0 = 0, mcDecayPosYD0 = 0, mcDecayPosR = 0, mcDecayPosZD0 =0, mcDecayPosRho=0.;
  Double_t mcDecayPosXD1 = 0, mcDecayPosYD1 = 0, mcDecayPosZD1 =0;

  Double_t lEtaCurrentPart =0., lPtCurrentPart = 0.,lThetaCurrentPart = 0., lPhiCurrentPart = 0.;
  Int_t iCurrentMother = 0;
  Double_t mcPosX = 0.,	mcPosY = 0.,mcPosZ = 0.;

  Double_t lPdgCurrentDaughter0 = 0, lPdgCurrentDaughter1= 0., /*lPdgCurrentMother=0.,*/lPdgCurrentDaughter =0;

  Double_t PxD0 = 0, PyD0 = 0,PzD0 = 0;
  Double_t PxD1 = 0, PyD1 = 0,PzD1 = 0;

  Int_t lNbMCPart           = 0;
  Int_t lPdgcodeCurrentPart = 0;
  //!----------------------------------------------------------------

  //! A set of very loose parameters for cuts 
  
  Double_t fgChi2max=33.;     //! max chi2
  Double_t fgDNmin=0.05;      //! min imp parameter for the 1st daughter = 500um
  Double_t fgDCAmax=1.;       //! max DCA between the daughter tracks in cm
  Double_t fgCPAmin=0.9;      //! min cosine of V0's pointing angle
  Double_t fgRmin=0.1;        //! min radius of the fiducial volume = 1 mm 
  Double_t fgRmax=200.;       //! max radius of the fiducial volume = 2 m

  //------------------------------------------
  // create pointer to event

  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }



//   AliVEvent *lESDevent = InputEvent();
//   if (!lESDevent) { 
//     Printf("ERROR: Could not retrieve event"); 
//     return; 
//  }
    
  Info("AliAnalysisTaskHelium3PiMC","Starting UserExec");  

  SetDataType("SIM");
  AliStack *stack = 0;
  if(fDataType == "SIM") {
    
    // Main loop
    // Called for EACH event
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent) { 
      Printf("ERROR: Could not retrieve MC event"); 
      return; 
    }
    
    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
    
    // set up a stack for use in check for primary/stable particles
    stack = mcEvent->Stack();
    if( !stack ) { Printf( "Stack not available"); return; }
  }
  
  else{
    Printf( "This Task Works Only on Simulation");
    return;
  }
  AliESDEvent *lESDevent = 0;

  //********************************** Connect to the InputEvent ******//
  
  //Int_t TrackNumber = 0;
  if(fAnalysisType == "ESD"){
    lESDevent = dynamic_cast<AliESDEvent*>(event);
    if (!lESDevent) {
      Printf("ERROR: lESDevent not available \n");
      return;
    }
  }

  else{
    Printf("This Analysis Works Only for ESD\n");
    return;
  }

  //*****************//  
  //*   Centrality  *//
  //*****************//

  AliCentrality *centrality = (AliCentrality*)lESDevent->GetCentrality();
 
  Float_t percentile=centrality->GetCentralityPercentile("V0M");
  
  //------------------------------

  runNumber = lESDevent->GetRunNumber();
  evNumber =lESDevent->GetEventNumberInFile();

  //---------------------

  //  Int_t primary = stack->GetNprimary();
  Int_t label =-1;

  lNbMCPart    = stack->GetNtrack();
      
  fHistMCMultiplicityTracks->Fill(lNbMCPart);     //histo

  TArrayD MomPionsMC(lNbMCPart);        //Neg pions
  Int_t nPionsMC=0;
  TArrayD MomHeMC(lNbMCPart);        //helium3
  Int_t nHeMC=0;

  //------ Trimomento pion
  TArrayD PxPionsMC(lNbMCPart);       
  Int_t nPxPionsMC=0;
  TArrayD PyPionsMC(lNbMCPart);       
  Int_t nPyPionsMC=0;
  TArrayD PzPionsMC(lNbMCPart);       
  Int_t nPzPionsMC=0;
  //------ Trimomento He
  TArrayD PxHeMC(lNbMCPart);       
  Int_t nPxHeMC=0;
  TArrayD PyHeMC(lNbMCPart);       
  Int_t nPyHeMC=0;
  TArrayD PzHeMC(lNbMCPart);       
  Int_t nPzHeMC=0;

  //Mass Definition

  Double_t        Helium3Mass = 2.80839; 
  Double_t        PionMass = 0.13957;    
  
  TLorentzVector  vPionMC,vHeliumMC,vSumMC;
  TLorentzVector  vPionMum,vHeliumMum,vSumMum;
  TLorentzVector  vPionRec,vHeliumRec,vSumRec;
  Bool_t isTwoBody=kFALSE;

  for (Int_t iMC=0; iMC<stack->GetNtrack(); iMC++)
    {
      TParticle *p0 = stack->Particle(iMC);
      
      if (!p0) {
	//Printf("ERROR: particle with label %d not found in stack (mc loop)", iMc);
	continue;
      }
      
      lPdgcodeCurrentPart  = p0->GetPdgCode();	
     
      if(lPdgcodeCurrentPart == 1000020030 || lPdgcodeCurrentPart == -1000020030 ){
	
	MomHeMC[nHeMC++]=p0->P();
	
        PxHeMC[nPxHeMC++]=p0->Px();
	PyHeMC[nPyHeMC++]=p0->Py();
	PzHeMC[nPzHeMC++]=p0->Pz();
	
	fhHeMCRigidity->Fill(p0->P()/2);
      }

      if(lPdgcodeCurrentPart == 211 || lPdgcodeCurrentPart == -211 ){

	MomPionsMC[nPionsMC++]=p0->P();
	
	PxPionsMC[nPxPionsMC++]=p0->Px();
	PyPionsMC[nPyPionsMC++]=p0->Py();
	PzPionsMC[nPzPionsMC++]=p0->Pz();
	
	fhPioneMC->Fill(p0->P());
      }
      
      if ( lPdgcodeCurrentPart == 1010010030 || lPdgcodeCurrentPart == -1010010030 ){

	lEtaCurrentPart   = p0->Eta();
	lPtCurrentPart    = p0->Pt();
	lThetaCurrentPart = p0->Theta();
	lPhiCurrentPart   = p0->Phi();
	iCurrentMother    = p0->GetFirstMother();
	
	fHistMCEta->Fill(lEtaCurrentPart);    
	fHistMCPt->Fill(lPtCurrentPart);      
	fHistMCTheta->Fill(lThetaCurrentPart);
	
	//	if (iCurrentMother == -1){lPdgCurrentMother=0; } else {lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();}
	 
	mcPosX = p0->Vx();
	mcPosY = p0->Vy();
	mcPosZ = p0->Vz();

	isTwoBody=kFALSE;
	
	for(Int_t i=p0->GetFirstDaughter(); i<= p0->GetLastDaughter(); i++){
	  TParticle *pDaughter = stack->Particle(i);
	  lPdgCurrentDaughter= pDaughter->GetPdgCode();
	  cout<<lPdgCurrentDaughter<<endl;
	  if(lPdgCurrentDaughter == 1000020030 || lPdgCurrentDaughter ==-1000020030 ){
	    isTwoBody=kTRUE;
	    
	  }
	}
	
	if(isTwoBody){
	  
	  for(Int_t i=p0->GetFirstDaughter(); i<= p0->GetLastDaughter(); i++){
	  
	    TParticle *pDaughter = stack->Particle(i);
	  
	    lPdgCurrentDaughter= pDaughter->GetPdgCode();
	  
	    if(lPdgCurrentDaughter == 211 || lPdgCurrentDaughter == -211 ){
	      id0 = i;
	    }
	    
	    if(lPdgCurrentDaughter == 1000020030 || lPdgCurrentDaughter == -1000020030 ){
	      id1 = i;
	    }
	  }
	  
	  TParticle *pDaughter0 = stack->Particle(id0);
	  TParticle *pDaughter1 = stack->Particle(id1);
	  lPdgCurrentDaughter0 = pDaughter0->GetPdgCode();
	  lPdgCurrentDaughter1 = pDaughter1->GetPdgCode();
	  
	  // Decay Radius and Production Radius
	  
	  if ( id0 <= lNbMCPart && id0 > 0 && id1 <= lNbMCPart && id1 > 0) {
	    
	    lPdgCurrentDaughter0 = pDaughter0->GetPdgCode();
	    lPdgCurrentDaughter1 = pDaughter1->GetPdgCode();
	    
	    PxD0 = pDaughter0->Px();
	    PyD0 = pDaughter0->Py();
	    PzD0 = pDaughter0->Pz();
	    
	    PxD1 = pDaughter1->Px();
	    PyD1 = pDaughter1->Py();
	    PzD1 = pDaughter1->Pz();
	    
	    mcDecayPosXD0 = pDaughter0->Vx();
	    mcDecayPosYD0 = pDaughter0->Vy();
	    mcDecayPosZD0 = pDaughter0->Vz();
	    
	    mcDecayPosXD1 = pDaughter0->Vx();
	    mcDecayPosYD1 = pDaughter0->Vy();
	    mcDecayPosZD1 = pDaughter0->Vz();
	    
	    mcDecayPosR = TMath::Sqrt(mcDecayPosXD0*mcDecayPosXD0+mcDecayPosYD0*mcDecayPosYD0);
	    fHistMCDecayPosition->Fill(mcDecayPosR);  
	    
	    mcDecayPosRho = TMath::Sqrt(mcDecayPosXD0*mcDecayPosXD0+mcDecayPosYD0*mcDecayPosYD0+mcDecayPosZD0*mcDecayPosZD0);
	    fHistMCDecayRho->Fill(mcDecayPosRho);  
	    
	    //---- Initial mass Test
	    
	    vHeliumMum.SetXYZM(PxD1,PyD1,PzD1,Helium3Mass); 
	    vPionMum.SetXYZM(PxD0,PyD0,PzD0,PionMass);       
	    vSumMum=vHeliumMum+vPionMum;
	    
	    fhInvMassMum->Fill(vSumMum.M());
	    
	    //Ntupla hyper triton
	    
	    vett2[0]=(Float_t)lESDevent->GetRunNumber();
	    vett2[1]=(Float_t)lESDevent->GetEventNumberInFile();
	    vett2[2]=(Float_t)iMC;
	    vett2[3]=(Float_t)percentile;
	    vett2[4]=(Float_t)mcPosX;
	    vett2[5]=(Float_t)mcPosY;
	    vett2[6]=(Float_t)mcPosZ;
	    vett2[7]=(Float_t)lPdgcodeCurrentPart;
	    vett2[8]=(Float_t)iCurrentMother;
	    vett2[9]=(Float_t)mcDecayPosXD0;
	    vett2[10]=(Float_t)mcDecayPosYD0;
	    vett2[11]=(Float_t)mcDecayPosZD0;
	    vett2[12]=(Float_t)mcDecayPosXD1;
	    vett2[13]=(Float_t)mcDecayPosYD1;
	    vett2[14]=(Float_t)mcDecayPosZD1;
	    vett2[15]=(Float_t)mcDecayPosRho;
	    vett2[16]=(Float_t)lEtaCurrentPart;
	    vett2[17]=(Float_t)p0->Y();
	    vett2[18]=(Float_t)lThetaCurrentPart;
	    vett2[19]=(Float_t)lPhiCurrentPart;
	    vett2[20]=(Float_t)p0->Px();
	    vett2[21]=(Float_t)p0->Py();
	    vett2[22]=(Float_t)p0->Pz();
	    vett2[23]=(Float_t)lPdgCurrentDaughter0;
	    vett2[24]=(Float_t)lPdgCurrentDaughter1;
	    vett2[25]=(Float_t)PxD0; //pion
	    vett2[26]=(Float_t)PyD0;
	    vett2[27]=(Float_t)PzD0;
	    vett2[28]=(Float_t)PxD1; //He3
	    vett2[29]=(Float_t)PyD1;
	    vett2[30]=(Float_t)PzD1;
	    
	    fNtuple2->Fill(vett2);
	    
	  }//if check daughters index
	}//is twobody
      } // if mother
    } // Kinetic Track loop ends here 
    
  // Loop phase - space
  
  Double_t HeMomMC =0;
  Double_t PionMomMC=0;
  Double_t PxHeMc=0, PyHeMc=0,PzHeMc=0;
  Double_t PxPionMc=0, PyPionMc=0,PzPionMc=0;

  for(Int_t l=0; l < nHeMC; l++){
      
    HeMomMC=MomHeMC[l];

    PxHeMc=PxHeMC[l];
    PyHeMc=PyHeMC[l];
    PzHeMc=PzHeMC[l];

    for(Int_t k=0; k < nPionsMC; k++){
   
      PionMomMC=MomPionsMC[k];
      
      PxPionMc=PxPionsMC[k];
      PyPionMc=PyPionsMC[k];
      PzPionMc=PzPionsMC[k];
      
      fhRigidityHevsMomPiMC->Fill(HeMomMC/2,PionMomMC);

      vHeliumMC.SetXYZM(PxHeMc,PyHeMc,PzHeMc,Helium3Mass); 
      vPionMC.SetXYZM(PxPionMc,PyPionMc,PzPionMc,PionMass);       
      vSumMC=vHeliumMC+vPionMC;
      
      fhInvMassMC->Fill(vSumMC.M());

    }
    
  } // end loop phase space

  //-------------- RECONSTRUCTION -------------------

  fHistEventMultiplicity->Fill(0);
  
  Double_t lMagneticField=lESDevent->GetMagneticField();

  Int_t TrackNumber = -1;

  // ANALISYS reconstructed tracks
  
  // Primary vertex cut
  
  const AliESDVertex *vtx = lESDevent->GetPrimaryVertexTracks();
  
  if(vtx->GetNContributors()<1) {
    
    // SPD vertex cut
    vtx = lESDevent->GetPrimaryVertexSPD();
    
    if(vtx->GetNContributors()<1) {
      Info("AliAnalysisTaskHelium3PiMC","No good vertex, skip event");
      return; // NO GOOD VERTEX, SKIP EVENT 
    }
  }

  fHistEventMultiplicity->Fill(1); // analyzed events with PV
 
  xPrimaryVertex=vtx->GetX();
  yPrimaryVertex=vtx->GetY();
  zPrimaryVertex=vtx->GetZ();  

  TrackNumber = lESDevent->GetNumberOfTracks();
  fHistTrackMultiplicity->Fill(TrackNumber); //tracce per evento

  fHistPercentileVsTrackNumber->Fill(percentile,TrackNumber);

  if (TrackNumber<2) return;  

  fHistEventMultiplicity->Fill(2);
    
  //Find Pair candidates
    
  TArrayI PionsTPC(TrackNumber);        //Neg pions
  Int_t nPionsTPC=0;
  
  TArrayI HeTPC(TrackNumber);        //helium3
  Int_t nHeTPC=0;
 
  // find pairs candidates phase daughter tracks rec

  TArrayD MomPionsRec(TrackNumber);        //Neg pions
  Int_t nPionsRec=0;
  
  TArrayD MomHeRec(TrackNumber);        //helium3
  Int_t nHeRec=0;

  //------ Trimomento pion
  TArrayD PxPionsRec(TrackNumber);       
  Int_t nPxPionsRec=0;
  TArrayD PyPionsRec(TrackNumber);       
  Int_t nPyPionsRec=0;
  TArrayD PzPionsRec(TrackNumber);       
  Int_t nPzPionsRec=0;

  //------ Trimomento He
  TArrayD PxHeRec(TrackNumber);       
  Int_t nPxHeRec=0;
  TArrayD PyHeRec(TrackNumber);       
  Int_t nPyHeRec=0;
  TArrayD PzHeRec(TrackNumber);       
  Int_t nPzHeRec=0;

  Float_t impactXY=-999, impactZ=-999;
  Float_t impactXYpi=-999, impactZpi=-999;
  
  Int_t PDGCodePos=0;
  Int_t PDGCodeNeg=0;
  Int_t motherPDGNeg=0;
  Int_t motherPDGPos=0;

  //!   SELECTIONS:
  //! - No ITSpureSA
  //! - ITSrefit
  //! - TPCrefit
  //! - ITSmap !=0

  // ******************* Track Cuts Definitions ********************//
  
  //! ITS

  AliESDtrackCuts* esdtrackCutsITS = new AliESDtrackCuts("esdtrackCutsITS");
  esdtrackCutsITS->SetRequireITSStandAlone(kFALSE);
  esdtrackCutsITS->SetRequireITSPureStandAlone(kFALSE);

  //! TPC
  
  Int_t    minclsTPC=60;
  Double_t maxchi2perTPCcl=5.;
  
  AliESDtrackCuts* esdtrackCutsTPC = new AliESDtrackCuts("esdtrackCutsTPC");
  esdtrackCutsTPC->SetRequireTPCRefit(kTRUE);
  esdtrackCutsTPC->SetAcceptKinkDaughters(kFALSE);
  esdtrackCutsTPC->SetMinNClustersTPC(minclsTPC);
  esdtrackCutsTPC->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
        
  //*************************************************************
  
  for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
    
    AliESDtrack *esdtrack=lESDevent->GetTrack(j);
   
    if(!esdtrack) { 
      AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
      continue; 
    }

    hBBTPCnoCuts->Fill(esdtrack->GetSign()*esdtrack->P(),esdtrack->GetTPCsignal());

    // ************** Track cuts ****************
    
    status  = (ULong_t)esdtrack->GetStatus();
    
    isTPC   = (((status) & (AliESDtrack::kTPCin))  != 0);
    
    Bool_t IsTrackAcceptedTPC =  esdtrackCutsTPC->AcceptTrack(esdtrack);
    Bool_t IsTrackAcceptedITS =  esdtrackCutsITS->AcceptTrack(esdtrack);
    
    if (!(IsTrackAcceptedTPC && IsTrackAcceptedITS)) continue;

    //----------------------------------------------
    
    //****** Cuts from  AliV0Vertex.cxx *************
     
    Double_t d=esdtrack->GetD(xPrimaryVertex,yPrimaryVertex,lMagneticField);
    //    if (TMath::Abs(d)<fgDPmin) continue;
    if (TMath::Abs(d)>fgRmax) continue;
    
    //---- (Usefull) Stuff
    
    TPCSignal=esdtrack->GetTPCsignal(); 
    
    if (TPCSignal<10)continue;
    
    if(!isTPC)continue;

    if(!esdtrack->GetTPCInnerParam())continue;
    
    AliExternalTrackParam trackIn(*esdtrack->GetInnerParam()); 
    pinTPC = trackIn.GetP(); 
        
    fhBBTPC->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
    
    d=esdtrack->GetD(xPrimaryVertex,yPrimaryVertex,lMagneticField);
    //    if (TMath::Abs(d)<fgDPmin) continue;
    if (TMath::Abs(d)>fgRmax) continue;

    label = TMath::Abs(esdtrack->GetLabel());
    
    if (label>=10000000) {
      // Underlying event. 10000000 is the
      // value of fkMASKSTEP in AliRunDigitizer
      cout <<"Strange, there should be no underlying event"<<endl;
    }
    
    else {
      if (label>=lNbMCPart) {
	cout <<"Strange, label outside the range"<< endl;
	continue;
      }
    }
    
    TParticle * part = stack->Particle(label);
    	
    Int_t PDGCode=part->GetPdgCode();
    
    if(PDGCode==-211)
      fhBBTPCNegativePions->Fill(esdtrack->GetSign()*esdtrack->P(),esdtrack->GetTPCsignal());
    
    if(PDGCode==+211)
      fhBBTPCPositivePions->Fill(esdtrack->GetSign()*esdtrack->P(),esdtrack->GetTPCsignal());


    //    if(PDGCode == 211){
 	  
    if(PDGCode==-211 || PDGCode==+211){
      
      PionsTPC[nPionsTPC++]=j;
    
      esdtrack->GetImpactParameters(impactXY, impactZ);
      fhPiDCAXY->Fill(impactXY);
      fhPiDCAZ->Fill(impactZ);
      
      MomPionsRec[nPionsRec++]=esdtrack->P();

      PxPionsRec[nPxPionsRec++]=esdtrack->Px();
      PyPionsRec[nPyPionsRec++]=esdtrack->Py();
      PzPionsRec[nPzPionsRec++]=esdtrack->Pz();
	
    }
     	
    if(PDGCode==1000020030 ||PDGCode==-1000020030 ){


      HeTPC[nHeTPC++]=j;

      fhBBTPCHe3->Fill(esdtrack->GetSign()*esdtrack->P(),esdtrack->GetTPCsignal());

      esdtrack->GetImpactParameters(impactXY, impactZ);
      fhHeDCAXY->Fill(impactXY);
      fhHeDCAZ->Fill(impactZ);
      
      MomHeRec[nHeRec++]=esdtrack->P();

      PxHeRec[nPxHeRec++]=esdtrack->Px();
      PyHeRec[nPyHeRec++]=esdtrack->Py();
      PzHeRec[nPzHeRec++]=esdtrack->Pz();
	
    }  
      
  }  //     ! track

  //-------------- LOOP pairs 1 -------------
  // Fill phase space and inva mass before cuts
  
  Double_t HeMomRec =0;
  Double_t PionMomRec=0;
  Double_t PxHeReco=0, PyHeReco=0,PzHeReco=0;
  Double_t PxPionReco=0, PyPionReco=0,PzPionReco=0;

  for(Int_t l=0; l < nHeRec; l++){
      
    HeMomRec=MomHeRec[l];

    PxHeReco=PxHeRec[l];
    PyHeReco=PyHeRec[l];
    PzHeReco=PzHeRec[l];

    for(Int_t k=0; k < nPionsRec; k++){
   
      PionMomRec=MomPionsRec[k];
      
      PxPionReco=PxPionsRec[k];
      PyPionReco=PyPionsRec[k];
      PzPionReco=PzPionsRec[k];
      
      fhRigidityHevsMomPiRec->Fill(HeMomRec,PionMomRec);

      vHeliumRec.SetXYZM(2*PxHeReco,2*PyHeReco,2*PzHeReco,Helium3Mass); 
      vPionRec.SetXYZM(PxPionReco,PyPionReco,PzPionReco,PionMass);       
      vSumRec=vHeliumRec+vPionRec;
      
      fhInvMassRec->Fill(vSumRec.M());
    }
    
  } // fine loop phase space

  //--------------- LOOP PAIRS ----------------------//
  
  Double_t        DcaHeToPrimVertex=0;
  Double_t        DcaPionToPrimVertex=0;
  
  impactXY=-999, impactZ=-999;
  impactXYpi=-999, impactZpi=-999;
  
  
  // Vettori per il PxPyPz
  
  Double_t momPionVett[3];
  for(Int_t i=0;i<3;i++)momPionVett[i]=0;
    
  Double_t momHeVett[3];
  for(Int_t i=0;i<3;i++)momHeVett[i]=0;
  
  //At SV
  
  Double_t momPionVettAt[3];
  for(Int_t i=0;i<3;i++)momPionVettAt[i]=0;
  
  Double_t momHeVettAt[3];
  for(Int_t i=0;i<3;i++)momHeVettAt[i]=0;
    
  Bool_t   IsHeITSRefit,IsPiITSRefit ;
  
  //----------- My 2nd Vertex Finder
    
  for (Int_t k=0; k < nPionsTPC; k++) {                           //! Pions Loop
    
    DcaPionToPrimVertex=0.;
    DcaHeToPrimVertex=0;
    
    Int_t PionIdx=PionsTPC[k];
    
    AliESDtrack  *PionTrack=lESDevent->GetTrack(PionIdx);
    
    statusPi = (ULong_t)PionTrack->GetStatus();
    IsPiITSRefit = ((statusPi) & (AliESDtrack::kITSrefit)); 
    
    Int_t labelPi = TMath::Abs(PionTrack->GetLabel());
    TParticle * partNeg = stack->Particle(labelPi);
    PDGCodeNeg=partNeg->GetPdgCode();
    
    Int_t mumidNeg = partNeg->GetFirstMother();
    if(mumidNeg>-1){
      TParticle *motherNeg=(TParticle*)stack->Particle(mumidNeg);
      motherPDGNeg = motherNeg->GetPdgCode();
    }
    
    if (PionTrack) 
      DcaPionToPrimVertex = TMath::Abs(PionTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField)); //OK
  
    if(DcaPionToPrimVertex<0.1)continue;     

    AliExternalTrackParam trackInPion(*PionTrack);  
    
    for (Int_t i=0; i<nHeTPC; i++){                               //! Helium Loop
      
      Int_t HeIdx=HeTPC[i];
      
      AliESDtrack  *HeTrack=lESDevent->GetTrack(HeIdx);
      
      statusT= (ULong_t)HeTrack->GetStatus();
      IsHeITSRefit = ((statusT) & (AliESDtrack::kITSrefit)); 
      
      Int_t labelHe = TMath::Abs(HeTrack->GetLabel());
      TParticle * partPos = stack->Particle(labelHe);
      PDGCodePos=partPos->GetPdgCode();
      
      Int_t mumidPos = partPos->GetFirstMother();
      if(mumidPos>-1){
	TParticle *motherPos=(TParticle*)stack->Particle(mumidPos);
	motherPDGPos = motherPos->GetPdgCode();
      }
      
      if (HeTrack) 
	DcaHeToPrimVertex = TMath::Abs(HeTrack->GetD(xPrimaryVertex, yPrimaryVertex,lMagneticField)); //OK
      
      AliExternalTrackParam trackInHe(*HeTrack); 
     
      HeTrack->PxPyPz(momHeVett);
      PionTrack->PxPyPz(momPionVett);   
     
      vHeliumRec.SetXYZM(2*momHeVett[0],2*momHeVett[1],2*momHeVett[2],Helium3Mass); 
      vPionRec.SetXYZM(momPionVett[0],momPionVett[1],momPionVett[2],PionMass);       
      vSumRec=vHeliumRec+vPionRec;
      
      fhInvMassRec1->Fill(vSumRec.M());

      fhInvMassRec2->Fill(vSumRec.M());
      
      if ( DcaPionToPrimVertex < fgDNmin)                //OK
	if ( DcaHeToPrimVertex < fgDNmin) continue;    //OK

      fhInvMassRec3->Fill(vSumRec.M());

      Double_t xn, xp;
      Double_t dca=0.;
      
      dca= PionTrack->GetDCA(HeTrack,lMagneticField,xn,xp); //!distance between two tracks (Neg to Pos)
      fHistProvaDCA->Fill(xn-xp,dca);
      if (dca > fgDCAmax) continue;

      fhInvMassRec4->Fill(vSumRec.M());
       
      if ((xn+xp) > 2*fgRmax) continue;
      if ((xn+xp) < 2*fgRmin) continue;
      fhInvMassRec5->Fill(vSumRec.M());
      
      //CORREZIONE da AliV0Vertex
      
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
      // Make  a  "V0" with Tracks                   //
      //=============================================//
      
      trackInPion.PropagateTo(xn,lMagneticField); 
      trackInHe.PropagateTo(xp,lMagneticField);
            
      AliESDv0 vertex(trackInPion,PionIdx,trackInHe,HeIdx);
      if (vertex.GetChi2V0() > fgChi2max) continue;
      fhInvMassRec6->Fill(vSumRec.M());

      Float_t CosPointingAngle=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex); //PointingAngle
      if (CosPointingAngle < fgCPAmin) continue;
    
      fhInvMassRec7->Fill(vSumRec.M());

      vertex.SetDcaV0Daughters(dca);
      vertex.SetV0CosineOfPointingAngle(CosPointingAngle);

      fPos[0]=vertex.Xv();
      fPos[1]=vertex.Yv(); 
      fPos[2]=vertex.Zv(); 
      

      
      Double_t raggio=TMath::Sqrt(fPos[0]*fPos[0]+fPos[1]*fPos[1]+fPos[2]*fPos[2]);
      HeTrack->GetPxPyPzAt(raggio,lMagneticField,momHeVettAt);
      PionTrack->GetPxPyPzAt(raggio,lMagneticField,momPionVettAt); 
      
      //------------------------------------------------------------------------//

      HeTrack->GetImpactParameters(impactXY, impactZ);
      
      PionTrack->GetImpactParameters(impactXYpi, impactZpi);
      
      Float_t timeTOFHe= HeTrack->GetTOFsignal();                 // ps
      Float_t trackLenghtTOFHe= HeTrack->GetIntegratedLength();  // cm

      Float_t timeTOFPi= PionTrack->GetTOFsignal();                 // ps
      Float_t trackLenghtTOFPi= PionTrack->GetIntegratedLength();  // cm

      //----------------------------------------------------------------------//

      vett1[0]=(Float_t)runNumber;
      vett1[1]=(Float_t)evNumber;
      vett1[2]=(Float_t)lNbMCPart;
      vett1[3]=(Float_t)percentile;
      vett1[4]=(Float_t)xPrimaryVertex; //PRIMARY
      vett1[5]=(Float_t)yPrimaryVertex;
      vett1[6]=(Float_t)zPrimaryVertex;
      vett1[7]=(Float_t)fPos[0]; //SECONDARY
      vett1[8]=(Float_t)fPos[1];
      vett1[9]=(Float_t)fPos[2];
      vett1[10]=(Float_t)dca;           //between 2 tracks
      vett1[11]=(Float_t)CosPointingAngle;          //cosPointingAngle da V0
      vett1[12]=(Float_t)vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
      vett1[13]=(Float_t)HeTrack->GetSign(); //He
      vett1[14]=(Float_t)trackInHe.GetP();
      vett1[15]=(Float_t)HeTrack->GetTPCsignal();
      vett1[16]=(Float_t)DcaHeToPrimVertex;
      vett1[17]=(Float_t)HeTrack->Eta();
      vett1[18]=(Float_t)momHeVett[0];
      vett1[19]=(Float_t)momHeVett[1];
      vett1[20]=(Float_t)momHeVett[2];
      vett1[21]=(Float_t)momHeVettAt[0];
      vett1[22]=(Float_t)momHeVettAt[1];
      vett1[23]=(Float_t)momHeVettAt[2];
      vett1[24]=(Float_t)HeTrack->GetTPCNcls();
      vett1[25]=(Float_t)impactXY;
      vett1[26]=(Float_t)impactZ;
      vett1[27]=(Float_t)isTOFHe3;
      vett1[28]=(Float_t)(trackLenghtTOFHe/timeTOFHe)/2.99792458e-2;
      vett1[29]=(Float_t)HeTrack->GetITSClusterMap();
      vett1[30]=(Float_t)IsHeITSRefit;
      vett1[31]=(Float_t)PionTrack->GetSign(); //Pion
      vett1[32]=(Float_t)trackInPion.GetP();
      vett1[33]=(Float_t)PionTrack->GetTPCsignal();
      vett1[34]=(Float_t)DcaPionToPrimVertex;
      vett1[35]=(Float_t)PionTrack->Eta();
      vett1[36]=(Float_t)momPionVett[0];
      vett1[37]=(Float_t)momPionVett[1];
      vett1[38]=(Float_t)momPionVett[2];
      vett1[39]=(Float_t)momPionVettAt[0];
      vett1[40]=(Float_t)momPionVettAt[1];
      vett1[41]=(Float_t)momPionVettAt[2];
      vett1[42]=(Float_t)PionTrack->GetTPCNcls();
      vett1[43]=(Float_t)impactXYpi;
      vett1[44]=(Float_t)impactZpi;
      vett1[45]=(Float_t)isTOFPi;
      vett1[46]=(Float_t)(trackLenghtTOFPi/timeTOFPi)/2.99792458e-2;
      vett1[47]=(Float_t)PionTrack->GetITSClusterMap();
      vett1[48]=(Float_t)IsPiITSRefit;
      vett1[49]=(Float_t)PDGCodeNeg;
      vett1[50]=(Float_t)PDGCodePos;
      vett1[51]=(Float_t)motherPDGNeg;
      vett1[52]=(Float_t)motherPDGPos;
      vett1[53]=(Float_t)labelPi;
      vett1[54]=(Float_t)labelHe;
      vett1[55]=(Float_t)mumidNeg;
      vett1[56]=(Float_t)mumidPos;

      fNtuple1->Fill(vett1);  

    }// positive TPC
    
  } //negative tpc

  PostData(1,fListHistCascade);
  
}// end userexec


//________________________________________________________________________
//
void AliAnalysisTaskHelium3PiMC::Terminate(Option_t *) 
{
  //  Draw result to the screen
  //   Called once at the end of the query
}

