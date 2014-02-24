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
  fListHistCascade(0), 
  fHistEventMultiplicity(0),         
  fHistTrackMultiplicity(0),      
  fHistTrackMultiplicityCent(0),      
  fHistTrackMultiplicitySemiCent(0),  
  fHistTrackMultiplicityMB(0),        
  fHistTrackMultiplicityPVCent(0),      
  fHistTrackMultiplicityPVSemiCent(0),  
  fHistTrackMultiplicityPVMB(0),        
  fHistMult(0),
  fhBB(0),    
  fhTOF(0),   
  fhMassTOF(0),
  fhBBPions(0),
  fhBBHe(0),   
  fhNaPos(0),  
  fhNaNeg(0),  
  fBetavsTPCsignalPos(0),  
  fBetavsTPCsignalNeg(0),  
  fHelium3TOF(0),   
  fNtuple1(0),
  fNtuple4(0),
  fPIDResponse(0)

{
  // Dummy Constructor 
}

//________________________________________________________________________
AliAnalysisTaskHelium3Pi::AliAnalysisTaskHelium3Pi(const char *name) 
: AliAnalysisTaskSE(name), 
  fAnalysisType("ESD"), 
  fCollidingSystems(0), 
  fESDtrackCuts(0),
  fDataType("REAL"),
  fListHistCascade(0), 
  fHistEventMultiplicity(0),    
  fHistTrackMultiplicity(0),            
  fHistTrackMultiplicityCent(0),      
  fHistTrackMultiplicitySemiCent(0),  
  fHistTrackMultiplicityMB(0),        
  fHistTrackMultiplicityPVCent(0),      
  fHistTrackMultiplicityPVSemiCent(0),  
  fHistTrackMultiplicityPVMB(0),        
  fHistMult(0),
  fhBB(0),    
  fhTOF(0),   
  fhMassTOF(0),
  fhBBPions(0),
  fhBBHe(0),   
  fhNaPos(0),  
  fhNaNeg(0),  
  fBetavsTPCsignalPos(0),  
  fBetavsTPCsignalNeg(0),  
  fHelium3TOF(0),                       
  fNtuple1(0),
  fNtuple4(0),
  fPIDResponse(0)  
  
{
  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container (Cascade)

  DefineOutput(1, TList::Class());

  // Int_t    minclsTPC=60;
  // Double_t maxchi2perTPCcl=5.;
  
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireITSStandAlone(kFALSE);
  fESDtrackCuts->SetRequireITSPureStandAlone(kFALSE);
      
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fESDtrackCuts->SetMinNClustersTPC(minclsTPC);
  // fESDtrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
  fESDtrackCuts->SetMinNClustersTPC(60);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  
}
//_______________________________________________________
AliAnalysisTaskHelium3Pi::~AliAnalysisTaskHelium3Pi() 
{ 
  // Destructor
  if (fListHistCascade) {
    delete fListHistCascade;
    fListHistCascade = 0;
  }
  
  if (fESDtrackCuts) delete fESDtrackCuts;
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
//___________________________________________

Bool_t AliAnalysisTaskHelium3Pi::IsTrackAccepted(AliVTrack *track){

  Bool_t testTrackCuts = kFALSE;

  Int_t    minclsTPC=60;
  Double_t maxchi2perTPCcl=5.;
  
  AliESDtrackCuts *fEsdTrackCuts = new AliESDtrackCuts("esdtrackCuts");
  fEsdTrackCuts->SetRequireITSStandAlone(kFALSE);
  fEsdTrackCuts->SetRequireITSPureStandAlone(kFALSE);
    
  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
  fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fEsdTrackCuts->SetMinNClustersTPC(minclsTPC);
  fEsdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);

  AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
 
  //if (fESDtrackCuts->AcceptTrack(esdtrack)) testTrackCuts = kTRUE;
  if (fEsdTrackCuts->AcceptTrack(esdtrack)) testTrackCuts = kTRUE;
  
  // Bool_t IsTrackAccepted =  fEsdTrackCuts->AcceptTrack(fEsdTrackCuts);
  // if (IsTrackAccepted) 
  //   return kTRUE;
  return testTrackCuts;

  delete fEsdTrackCuts;
  
}
//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskHelium3Pi::UserCreateOutputObjects()
{

  fListHistCascade = new TList();
  fListHistCascade->SetOwner();  // IMPORTANT!

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

    fListHistCascade->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity   = new TH2F( "fHistTrackMultiplicity" , "Nb of Tracks", 25000,0, 25000,210,-1,104);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicity);
  } 

  if(! fHistTrackMultiplicityCent ){
    fHistTrackMultiplicityCent   = new TH2F( "fHistTrackMultiplicityCent", "Nb of Tracks Central Events", 25000,0, 25000,210,-1,104 );
    fHistTrackMultiplicityCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityCent->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityCent);
  } 

  if(! fHistTrackMultiplicitySemiCent ){
    fHistTrackMultiplicitySemiCent   = new TH2F( "fHistTrackMultiplicitySemiCent" , "Nb of Tracks SemiCentral Events", 25000,0, 25000 ,210,-1,104);
    fHistTrackMultiplicitySemiCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicitySemiCent->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicitySemiCent);
  } 
 
  if(! fHistTrackMultiplicityMB ){
    fHistTrackMultiplicityMB   = new TH2F( "fHistTrackMultiplicityMB" , "Nb of Tracks MBral Events", 25000,0, 25000,210,-1,104 );
    fHistTrackMultiplicityMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityMB->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityMB);
  } 

  if(! fHistTrackMultiplicityPVCent ){
    fHistTrackMultiplicityPVCent   = new TH2F( "fHistTrackMultiplicityPVCent" , "Nb of Tracks Central Events", 25000,0, 25000,210,-1,104 );
    fHistTrackMultiplicityPVCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVCent->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityPVCent);
  } 

  if(! fHistTrackMultiplicityPVSemiCent ){
    fHistTrackMultiplicityPVSemiCent   = new TH2F( "fHistTrackMultiplicityPVSemiCent" , "Nb of Tracks SemiCentral Events", 25000,0, 25000 ,210,-1,104);
    fHistTrackMultiplicityPVSemiCent->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVSemiCent->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityPVSemiCent);
  } 
 
  if(! fHistTrackMultiplicityPVMB ){
    fHistTrackMultiplicityPVMB   = new TH2F( "fHistTrackMultiplicityPVMB" , "Nb of Tracks MBral Events", 25000,0, 25000,210,-1,104 );
    fHistTrackMultiplicityPVMB->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicityPVMB->GetYaxis()->SetTitle("Percentile");
    fListHistCascade->Add(fHistTrackMultiplicityPVMB);
  } 

  if(! fHistMult){
    fHistMult=new TH1F ("fHistMult","Number neg-pos", 10, -1,9);
    fHistMult->GetXaxis()->SetTitle("Type of tracks");
    fListHistCascade->Add(fHistMult);
  }
 
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 1000,-6,6,1500,0,5000);
    fhBB->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBB->GetYaxis()->SetTitle("TPC Signal");
    fListHistCascade->Add(fhBB);
  }

  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 1000,-6,6,1000,0,1.2);
    fhTOF->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhTOF->GetYaxis()->SetTitle("#beta");
    fListHistCascade->Add(fhTOF);
  }

  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 300,0 ,5);
    fhMassTOF->GetXaxis()->SetTitle("Mass (GeV/#it{c}^{2})");
    fListHistCascade->Add(fhMassTOF);
  }

  if(! fhBBPions ){
    fhBBPions = new TH2F( "fhBBPions" , "Bethe-Bloch TPC Pions" , 1000,-6,6,1500,0,5000);
    fhBBPions->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBPions->GetYaxis()->SetTitle("TPC Signal");
    fListHistCascade->Add(fhBBPions);
  }
  
  if(! fhBBHe ){
    fhBBHe = new TH2F( "fhBBHe" , "Bethe-Bloch TPC He" , 1000,-6,6,1500,0,5000);
    fhBBHe->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhBBHe->GetYaxis()->SetTitle("TPC Signal");
    fListHistCascade->Add(fhBBHe);
  }
  
  if(! fhNaPos ){
    fhNaPos = new TH2F( "fhNaPos" , "Distribution Pos" , 500,0,5,500,-10,10);
    fhNaPos->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhNaPos->GetYaxis()->SetTitle("(TPCSignal-bbtheo)/bbtheo (He)");
    fListHistCascade->Add(fhNaPos);
  }
  
  if(! fhNaNeg ){
    fhNaNeg = new TH2F( "fhNaNeg" , "Distribution Neg" , 500,0,5,500,-10,10);
    fhNaNeg->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fhNaNeg->GetYaxis()->SetTitle("(TPCSignal-bbtheo)/bbtheo (He)");
    fListHistCascade->Add(fhNaNeg);
  }

  if(! fBetavsTPCsignalPos ){
    fBetavsTPCsignalPos = new TH2F("fBetavsTPCsignalPos","fBetavsTPCsignalPos",1000,0,1.2,1500,0,5000);
    fBetavsTPCsignalPos->GetXaxis()->SetTitle("#beta");
    fBetavsTPCsignalPos->GetYaxis()->SetTitle("TPC Signal");
    fListHistCascade->Add(fBetavsTPCsignalPos);
  }
  
  if(! fBetavsTPCsignalNeg ){
    fBetavsTPCsignalNeg = new TH2F("fBetavsTPCsignalNeg","fBetavsTPCsignalNeg",1000,0,1.2,1500,0,5000);
    fBetavsTPCsignalNeg->GetXaxis()->SetTitle("#beta");
    fBetavsTPCsignalNeg->GetYaxis()->SetTitle("TPC Signal");
    fListHistCascade->Add(fBetavsTPCsignalNeg);
  }
  
  if(! fHelium3TOF){
    fHelium3TOF = new TH2F("fHelium3massTOF","Helium3 beta vs p/z",1000,0,6,1000,0,1.2);
    fHelium3TOF->GetXaxis()->SetTitle("p/z (GeV/#it{c})");
    fHelium3TOF->GetYaxis()->SetTitle("#beta");
    fListHistCascade->Add(fHelium3TOF);
  }
  
  if(! fNtuple1 ) {
    fNtuple1 = new TNtuple("fNtuple1","Ntuple1","runNumber:bunchcross:orbit:period:eventtype:TrackNumber:percentile:xPrimaryVertex:yPrimaryVertex:zPrimaryVertex:xSecondaryVertex:ySecondaryVertex:zSecondaryVertex:dcaTracks:CosPointingAngle:DCAV0toPrimaryVertex:HeSign:HepInTPC:HeTPCsignal:DcaHeToPrimVertex:HeEta:momHex:momHey:momHez:momHeAtSVx:momHeAtSVy:momHeAtSVz:HeTPCNcls:HeimpactXY:HeimpactZ:HeITSClusterMap:IsHeITSRefit:PionSign:PionpInTPC:PionTPCsignal:DcaPionToPrimVertex:PionEta:momPionx:momPiony:momPionz:momNegPionAtSVx:momNegPionAtSVy:momNegPionAtSVz:PionTPCNcls:PionimpactXY:PionimpactZ:PionITSClusterMap:IsPiITSRefit:xn:xp:chi2He:chi2Pi");
    
    fListHistCascade->Add(fNtuple1);
  }
  
  if(! fNtuple4 ) {
    fNtuple4 = new TNtuple("fNtuple4","Ntuple4","runNumber:BCNumber:OrbitNumber:PeriodNumber:eventtype:isHeITSrefit:percentile:Sign:pinTPC:GetTPCsignal:Px:Py:Pz:Eta:isTOF:poutTPC:timeTOF:trackLenghtTOF:impactXY:impactZ:mapITS:TPCNcls:TRDsignal:xPrimaryVertex:yPrimaryVertex:zPrimaryVertex:chi2PerClusterTPC");
    fListHistCascade->Add(fNtuple4);
  } 

  PostData(1,  fListHistCascade);

}// end UserCreateOutputObjects


//====================== USER EXEC ========================

void AliAnalysisTaskHelium3Pi::UserExec(Option_t *) 
{
  //_______________________________________________________________________
  
  //!*********************!//
  //!  Define variables   !//
  //!*********************!//

  Float_t vett1[52];
  for(Int_t i=0;i<52;i++) vett1[i]=0;

  Float_t vett4[40];
  for(Int_t i=0;i<40;i++) vett4[i]=0;
  
  Double_t pinTPC=0.,poutTPC=0.,TPCSignal=0.;
  Double_t xPrimaryVertex=0.,yPrimaryVertex=0.,zPrimaryVertex=0.;
  Double_t massTOF=0.,timeTOF=0.,trackLenghtTOF=0.,betaTOF=0.;

  ULong_t  status=0;
  //  ULong_t  statusT=0;
  ULong_t  statusPi=0;

  Bool_t   isTPC=kFALSE,isTOF=kFALSE,IsHeITSRefit=kFALSE,IsPiITSRefit=kFALSE ;

  Float_t nSigmaNegPion=0.;
  // Float_t nSigmaNegPion1=0.;
  // Float_t nSigmaNegPion2=0.;

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

  //===========================================

  Int_t eventtype=-99;
  
  Bool_t isSelectedCentral     = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
  Bool_t isSelectedSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral);
  Bool_t isSelectedMB          = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
 
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

    Bool_t testTrackCuts = kFALSE;

    //!   SELECTIONS:
    //! - No ITSpureSA
    //! - ITSrefit
    //! - TPCrefit
    //! - ITSmap !=0
    
    // // ******************* Track Cuts Definitions ********************//
  
    // //! ITS

    // AliESDtrackCuts* esdtrackCutsITS = new AliESDtrackCuts("esdtrackCutsITS");
    // esdtrackCutsITS->SetRequireITSStandAlone(kFALSE);
    // esdtrackCutsITS->SetRequireITSPureStandAlone(kFALSE);
    
    // //! TPC
    
    // Int_t    minclsTPC=60;
    // Double_t maxchi2perTPCcl=5.;
    
    // AliESDtrackCuts* esdtrackCutsTPC = new AliESDtrackCuts("esdtrackCutsTPC");
    // esdtrackCutsTPC->SetRequireTPCRefit(kTRUE);
    // esdtrackCutsTPC->SetAcceptKinkDaughters(kFALSE);
    // esdtrackCutsTPC->SetMinNClustersTPC(minclsTPC);
    // esdtrackCutsTPC->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
    
    // //*********************************************+
    
    runNumber = lESDevent->GetRunNumber();
    //    evNumber  = lESDevent->GetEventNumberInFile();
    
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
      
      testTrackCuts = IsTrackAccepted(esdtrack);
      if(testTrackCuts == kFALSE)continue;

      //      cout<<"Is ttrack accepted: "<<(testTrackCuts)<<endl;
      //      if(IsTrackAccepted(esdtrack)==kFALSE)continue;

      // ************** Track cuts ****************
      
      status  = (ULong_t)esdtrack->GetStatus();
      
      isTPC   = (((status) & AliESDtrack::kTPCin)  != 0);
      isTOF   = ((((status) & AliESDtrack::kTOFout) != 0) && (((status) & AliESDtrack::kTIME) != 0));

      // if (!IsTrackAccepted(esdtrack))continue;
      // Bool_t IsTrackAcceptedTPC =  esdtrackCutsTPC->AcceptTrack(esdtrack);
      // Bool_t IsTrackAcceptedITS =  esdtrackCutsITS->AcceptTrack(esdtrack);
      // if (!(IsTrackAcceptedTPC && IsTrackAcceptedITS)) continue;
      
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
      
      fHistMult->Fill(0);
      
      if((status) & (AliESDtrack::kITSrefit!=0)){
	fHistMult->Fill(1);
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
      
      fHistMult->Fill(2);
      
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
      
      bbtheoM = (fPIDResponse->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType) 7));
      
      //      if( TPCSignal > bbtheoM ) {
      if( bbtheoM > -3.) {
	
	if(pinTPC>0.6){
	  
	  fhBBHe->Fill(pinTPC*esdtrack->GetSign(),TPCSignal);
	  HeTPC[nHeTPC++]=j;
	  
	  Bool_t isHeITSrefit=((status) & (AliESDtrack::kITSrefit));
	  
	  esdtrack->GetImpactParameters(impactXY, impactZ);
	  
	  Int_t  fIdxInt[200]; //dummy array
	  Int_t nClustersTPC = esdtrack->GetTPCclusters(fIdxInt);
	  
	  Float_t chi2PerClusterTPC = esdtrack->GetTPCchi2()/(Float_t)(nClustersTPC);
	  
	  vett4[0] =(Float_t)runNumber;
	  vett4[1] =(Float_t)BCNumber;
	  vett4[2] =(Float_t)OrbitNumber;
	  vett4[3] =(Float_t)PeriodNumber;
	  vett4[4] =(Float_t)eventtype;
	  vett4[5] =(Float_t)isHeITSrefit;
	  vett4[6] =(Float_t)percentile;
	  vett4[7] =(Float_t)esdtrack->GetSign();
	  vett4[8] =(Float_t)pinTPC;
	  vett4[9] =(Float_t)esdtrack->GetTPCsignal();
	  vett4[10]=(Float_t)esdtrack->Px();
	  vett4[11]=(Float_t)esdtrack->Py();
	  vett4[12]=(Float_t)esdtrack->Pz();
	  vett4[13]=(Float_t)esdtrack->Eta();
	  vett4[14]=(Float_t)isTOF;
	  vett4[15]=(Float_t)poutTPC;
	  vett4[16]=(Float_t)timeTOF;
	  vett4[17]=(Float_t)trackLenghtTOF;
	  vett4[18]=(Float_t)impactXY;
	  vett4[19]=(Float_t)impactZ;
	  vett4[20]=(Float_t)mapITS;
	  vett4[21]=(Float_t)esdtrack->GetTPCNcls();
	  vett4[22]=(Float_t)esdtrack->GetTRDsignal();
	  vett4[23]=(Float_t)xPrimaryVertex;
	  vett4[24]=(Float_t)yPrimaryVertex;
	  vett4[25]=(Float_t)zPrimaryVertex;
	  vett4[26]=(Float_t)chi2PerClusterTPC;
	  
	  fNtuple4->Fill(vett4);
	}
      }
    }       //! track
    
    
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
      
      if(DcaPionToPrimVertex<0.2)continue; //qui
      
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
	
	vett1[0] =(Float_t)runNumber;
	vett1[1] =(Float_t)BCNumber;
	vett1[2] =(Float_t)OrbitNumber;
	vett1[3] =(Float_t)PeriodNumber;
	vett1[4] =(Float_t)eventtype;
	vett1[5] =(Float_t)TrackNumber;
	vett1[6] =(Float_t)percentile;
	vett1[7] =(Float_t)xPrimaryVertex; //PRIMARY
	vett1[8] =(Float_t)yPrimaryVertex;
	vett1[9] =(Float_t)zPrimaryVertex;
	vett1[10]=(Float_t)fPos[0]; //SECONDARY
	vett1[11]=(Float_t)fPos[1];
	vett1[12]=(Float_t)fPos[2];
	vett1[13]=(Float_t)dca;           //between 2 tracks
	vett1[14]=(Float_t)CosPointingAngle;          //cosPointingAngle da V0
	vett1[15]=(Float_t)vertex.GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
	vett1[16]=(Float_t)HeTrack->GetSign(); //He
	vett1[17]=(Float_t)trackInHe.GetP();
	vett1[18]=(Float_t)HeTrack->GetTPCsignal();
	vett1[19]=(Float_t)DcaHeToPrimVertex;
	vett1[20]=(Float_t)HeTrack->Eta();
	vett1[21]=(Float_t)momHeVett[0];
	vett1[22]=(Float_t)momHeVett[1];
	vett1[23]=(Float_t)momHeVett[2];
	vett1[24]=(Float_t)momHeVettAt[0];
	vett1[25]=(Float_t)momHeVettAt[1];
	vett1[26]=(Float_t)momHeVettAt[2];
	vett1[27]=(Float_t)HeTrack->GetTPCNcls();
	vett1[28]=(Float_t)impactXY;
	vett1[29]=(Float_t)impactZ;
	vett1[30]=(Float_t)HeTrack->GetITSClusterMap();
	vett1[31]=(Float_t)IsHeITSRefit;
	vett1[32]=(Float_t)PionTrack->GetSign(); //Pion
	vett1[33]=(Float_t)trackInPion.GetP();
	vett1[34]=(Float_t)PionTrack->GetTPCsignal();
	vett1[35]=(Float_t)DcaPionToPrimVertex;
	vett1[36]=(Float_t)PionTrack->Eta();
	vett1[37]=(Float_t)momPionVett[0];
	vett1[38]=(Float_t)momPionVett[1];
	vett1[39]=(Float_t)momPionVett[2];
	vett1[40]=(Float_t)momPionVettAt[0];
	vett1[41]=(Float_t)momPionVettAt[1];
	vett1[42]=(Float_t)momPionVettAt[2];
	vett1[43]=(Float_t)PionTrack->GetTPCNcls();
	vett1[44]=(Float_t)impactXYpi;
	vett1[45]=(Float_t)impactZpi;
	vett1[46]=(Float_t)PionTrack->GetITSClusterMap();
	vett1[47]=(Float_t)IsPiITSRefit;
	vett1[48]=(Float_t)xn;
	vett1[49]=(Float_t)xp;
	vett1[50]=(Float_t)HeTrack->GetTPCchi2()/(Float_t)(nClustersTPCHe);
	vett1[51]=(Float_t)PionTrack->GetTPCchi2()/(Float_t)(nClustersTPCPi);

	fNtuple1->Fill(vett1);  
	vertex.Delete();
      }// positive TPC
      
    } //negative tpc
    
  }
  
  PostData(1,fListHistCascade);
  
} //end userexec


//________________________________________________________________________

void AliAnalysisTaskHelium3Pi::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}


