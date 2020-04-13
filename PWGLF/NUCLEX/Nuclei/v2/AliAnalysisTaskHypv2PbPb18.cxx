/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright noticxse appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//                 AliAnalysisTaskHypv2PbPb18 class
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

#include <iostream>
#include <TGrid.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include <TRandom3.h>
#include "TFile.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliOADBContainer.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskHypv2PbPb18.h"

ClassImp(AliAnalysisTaskHypv2PbPb18)

using std::cout;
using std::endl;

//_____________________________________________________________________________
AliAnalysisTaskHypv2PbPb18::AliAnalysisTaskHypv2PbPb18():
  AliAnalysisTaskSE(),                       //! 
  fESDevent(0),                         //! 
  fevent(0),   
  fRun(-1),
  fMultV0(0),
  fQxnmV0A(0),
  fQynmV0A(0),
  fQxnsV0A(0),
  fQynsV0A(0),
  fQxnmV0C(0),
  fQynmV0C(0),
  fQxnsV0C(0),
  fQynsV0C(0),
  fRecPass(0),
  fCenCalV0(0), //da qui --> centrality selection
  fFilterBit(4),
  fptc(1),     
  fVzmax(10),
  fPeriod(1),
  fNHarm(2),
  fListHist(0), 
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fhBB(0),
  fhBBHyp(0),
  fhBBAHyp(0),
  EPVzAvsCentrality(0), 
  EPVzCvsCentrality(0), 
  hQVzAQVzCvsCentrality(0),
  hQVzAQTPCvsCentrality(0),
  hQVzCQTPCvsCentrality(0),
  hQxVzAvsCentrality(0),
  hQyVzAvsCentrality(0),
  hQxVzCvsCentrality(0),
  hQyVzCvsCentrality(0),
  hCos2DeltaTPCVzAvsCentrality(0),
  hCos2DeltaTPCVzCvsCentrality(0),
  hCos2DeltaVzAVzCvsCentrality(0),
  hCos2DeltaVzATPCvsCentrality(0),
  hCos2DeltaVzCTPCvsCentrality(0),
  hCos2DeltaVzCVzAvsCentrality(0),
  eventtype(-999),
  ftree(0),           
  iEvent(0),
  zVertex(0),
  centrality(0),
  px_Daughter1(0),
  py_Daughter1(0),
  pz_Daughter1(0),
  q_Daughter1(0),
  dcaxy_Daughter1(0),
  nTPC_Clusters_Daughter1(0),
  nTPC_Clusters_dEdx_Daughter1(0),
  chi2_TPC_Daughter1(0),
  nSigmaTPC_He3_Daughter1(0),
  nSigmaTPC_Pion_Daughter1(0),
  px_Daughter2(0),
  py_Daughter2(0),
  pz_Daughter2(0),
  q_Daughter2(0),
  dcaxy_Daughter2(0),
  nTPC_Clusters_Daughter2(0),
  nTPC_Clusters_dEdx_Daughter2(0),
  chi2_TPC_Daughter2(0),
  nSigmaTPC_He3_Daughter2(0),
  nSigmaTPC_Pion_Daughter2(0),
  isOnTheFlyV0(0),
  cosPointingAngle(0),
  dcaV0Daughters(0),
  radius(0),
  chi2V0(0),
  decayLength(0),
  deltaphiV0A(0),
  deltaphiV0C(0),
  uqV0A(0),
  uqV0C(0),
  fPIDResponse(0),
  fESDtrackCuts_Pos(0),
  fESDtrackCuts_Neg(0),
  fESDtrackCutsEP(0),
  fEventCuts(0)
{
  cout<<"Dummy constructor"<<endl;
  fESDtrackCuts_Pos = new AliESDtrackCuts("fESDtrackCuts_Pos","fESDtrackCuts_Pos");
  fESDtrackCuts_Neg = new AliESDtrackCuts("fESDtrackCuts_Neg","fESDtrackCuts_Neg");
  fESDtrackCutsEP   = new AliESDtrackCuts("AliESDtrackCutsEP","AliESDtrackCutsEP");
  Initialize();
}

//______________________________________________________________________________
AliAnalysisTaskHypv2PbPb18::AliAnalysisTaskHypv2PbPb18(const char *name):
    AliAnalysisTaskSE(name),                   //! 
    fESDevent(0),                         //! 
    fevent(0),   
    fRun(-1),
    fMultV0(0),
    fQxnmV0A(0),
    fQynmV0A(0),
    fQxnsV0A(0),
    fQynsV0A(0),
    fQxnmV0C(0),
    fQynmV0C(0),
    fQxnsV0C(0),
    fQynsV0C(0),
    fRecPass(0),
    fCenCalV0(0), //da qui
    fFilterBit(4),
    fptc(1),     
    fVzmax(10),
    fPeriod(1),
    fNHarm(2),
    fListHist(0), 
    fHistEventMultiplicity(0), 
    fHistTrackMultiplicity(0),
    fhBB(0),
    fhBBHyp(0),
    fhBBAHyp(0),
    EPVzAvsCentrality(0), 
    EPVzCvsCentrality(0), 
    hQVzAQVzCvsCentrality(0),
    hQVzAQTPCvsCentrality(0),
    hQVzCQTPCvsCentrality(0),
    hQxVzAvsCentrality(0),
    hQyVzAvsCentrality(0),
    hQxVzCvsCentrality(0),
    hQyVzCvsCentrality(0),
    hCos2DeltaTPCVzAvsCentrality(0),
    hCos2DeltaTPCVzCvsCentrality(0),
    hCos2DeltaVzAVzCvsCentrality(0),
    hCos2DeltaVzATPCvsCentrality(0),
    hCos2DeltaVzCTPCvsCentrality(0),
    hCos2DeltaVzCVzAvsCentrality(0),			     
    eventtype(-999),
    ftree(0),           
    iEvent(0),
    zVertex(0),
    centrality(0),
    px_Daughter1(0),
    py_Daughter1(0),
    pz_Daughter1(0),
    q_Daughter1(0),
    dcaxy_Daughter1(0),
    nTPC_Clusters_Daughter1(0),
    nTPC_Clusters_dEdx_Daughter1(0),
    chi2_TPC_Daughter1(0),
    nSigmaTPC_He3_Daughter1(0),
    nSigmaTPC_Pion_Daughter1(0),
    px_Daughter2(0),
    py_Daughter2(0),
    pz_Daughter2(0),
    q_Daughter2(0),
    dcaxy_Daughter2(0),
    nTPC_Clusters_Daughter2(0),
    nTPC_Clusters_dEdx_Daughter2(0),
    chi2_TPC_Daughter2(0),
    nSigmaTPC_He3_Daughter2(0),
    nSigmaTPC_Pion_Daughter2(0),
    isOnTheFlyV0(0),
    cosPointingAngle(0),
    dcaV0Daughters(0),
    radius(0),
    chi2V0(0),
    decayLength(0),
    deltaphiV0A(0),
    deltaphiV0C(0),
    uqV0A(0),
    uqV0C(0),
    fPIDResponse(0),
    fESDtrackCuts_Pos(0),
    fESDtrackCuts_Neg(0),
    fESDtrackCutsEP(0),
    fEventCuts(0)
{
  
  //
  cout<<"Real constructor"<<endl;
  fESDtrackCuts_Pos = new AliESDtrackCuts("fESDtrackCuts_Pos","fESDtrackCuts_Pos");
  fESDtrackCuts_Neg = new AliESDtrackCuts("fESDtrackCuts_Neg","fESDtrackCuts_Neg");
  fESDtrackCutsEP   = new AliESDtrackCuts("AliESDtrackCutsEP","AliESDtrackCutsEP");
  Initialize();
  //  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class()); 
}
//------------------------------------------
void AliAnalysisTaskHypv2PbPb18::Initialize()
{
  //
  // updating parameters in case of changes
  //
  // fESDtrackCuts_Pos = new AliESDtrackCuts("fESDtrackCuts_Pos","fESDtrackCuts_Pos");
  // fESDtrackCuts_Neg = new AliESDtrackCuts("fESDtrackCuts_Neg","fESDtrackCuts_Neg");
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(fisPrimCut,kTRUE);
  // fESDtrackCuts->SetMaxDCAToVertexXY(3);
  // fESDtrackCuts->SetMaxDCAToVertexZ(2);
  // fESDtrackCuts->SetEtaRange(-0.8,0.8);
  
  fESDtrackCutsEP = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); 
  
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::PassedBasicTrackQualityCuts_Pos (AliESDtrack *track)  {
    
    fESDtrackCuts_Pos -> SetAcceptKinkDaughters(false);
    fESDtrackCuts_Pos -> SetMinNClustersTPC(50);
    fESDtrackCuts_Pos -> SetRequireTPCRefit(true);
    fESDtrackCuts_Pos -> SetMaxChi2PerClusterTPC(10.0);
    fESDtrackCuts_Pos -> SetEtaRange (-1.0,1.0);
    if ( !fESDtrackCuts_Pos->AcceptTrack (track) ) return false;
    return true;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::PassedBasicTrackQualityCuts_Neg (AliESDtrack *track)  {
    
    fESDtrackCuts_Neg -> SetAcceptKinkDaughters(false);
    fESDtrackCuts_Neg -> SetMinNClustersTPC(50);
    fESDtrackCuts_Neg -> SetRequireTPCRefit(true);
    fESDtrackCuts_Neg -> SetMaxChi2PerClusterTPC(10.0);
    fESDtrackCuts_Neg -> SetEtaRange (-1.0,1.0);
    if ( !fESDtrackCuts_Neg->AcceptTrack (track) ) return false;
    return true;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::PassedMinimalQualityCutsV0 (AliESDv0 *V0)  {
    
    //Basic Cuts
    if (V0->GetDcaV0Daughters()>2.0) return false;
    if (V0->GetRr()<3.0) return false;
    if (V0->GetV0CosineOfPointingAngle()<0.9) return false;
    
    return true;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::IsHyperTritonCandidate (AliESDv0 *V0)  {
    
    //Get V0 Daughters
    AliESDtrack *trackPos = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *trackNeg = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Pair Requirements
    if ( (!IsPionCandidate (trackPos)) && (!IsPionCandidate (trackNeg))) return false;
    if ( (!Is3HeCandidate  (trackPos)) && (!Is3HeCandidate  (trackNeg))) return false;
    if ( IsPionCandidate   (trackPos)  && (!Is3HeCandidate  (trackNeg))) return false;
    if ( Is3HeCandidate    (trackPos)  && (!IsPionCandidate (trackNeg))) return false;
    

    //Momentum Components of V0 Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
    
    //Hypertriton
    if (Is3HeCandidate (trackPos) && IsPionCandidate (trackNeg)) {
        
        TVector3 P1 (2.0*posMomentum[0],2.0*posMomentum[1],2.0*posMomentum[2]);
        TVector3 P2 (negMomentum[0],negMomentum[1],negMomentum[2]);
        Double_t m = InvariantMassHypertriton (P1,P2);
	

        if (m>3.1){
	  fhBBHyp->Fill(trackPos->GetTPCmomentum()*trackPos->Charge(),trackPos->GetTPCsignal());
	  fhBBHyp->Fill(trackNeg->GetTPCmomentum()*trackNeg->Charge(),trackNeg->GetTPCsignal());
	  return false;
	}
    }
    //Anti-Hypertriton
    if (IsPionCandidate (trackPos) && Is3HeCandidate (trackNeg)) {
        
        TVector3 P1 (2.0*negMomentum[0],2.0*negMomentum[1],2.0*negMomentum[2]);
        TVector3 P2 (posMomentum[0],posMomentum[1],posMomentum[2]);
        Double_t m = InvariantMassHypertriton (P1,P2);

        if (m>3.1) {
	  fhBBAHyp->Fill(trackPos->GetTPCmomentum()*trackPos->Charge(),trackPos->GetTPCsignal());
	  fhBBAHyp->Fill(trackNeg->GetTPCmomentum()*trackNeg->Charge(),trackNeg->GetTPCsignal());
	  return false;
	}
    }
    
    return true;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::IsPionCandidate (AliESDtrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
    if (TMath::Abs(nsigmaTPC) > 4.0) return false;
    if (track->Pt()>1.3) return false;
    
    Double_t dca_xy = GetTransverseDCA (track);
    if (TMath::Abs (dca_xy) < 0.1) return false;
    
    return true;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskHypv2PbPb18::Is3HeCandidate (AliESDtrack *track)  {
  Double_t fParamHe3[5]   = { 1.74962, 27.4992, 4.00313e-15, 2.48485, 8.31768};
  //Variables
  Double_t p = track->GetInnerParam()->GetP();
  Double_t mass = AliPID::ParticleMass (AliPID::kHe3);
  Double_t dEdx_au = track->GetTPCsignal();

  //Expected dE/dx for 3He
  Float_t hel3Exp = 4.0*AliExternalTrackParam::BetheBlochAleph(2.0*p/mass,fParamHe3[0],fParamHe3[1],fParamHe3[2],fParamHe3[3],fParamHe3[4]);
  Double_t sigma = 0.07;//dE/dx Resolution for 3He (7%)
  Double_t nSigmaHe3  = (dEdx_au - hel3Exp)/(sigma*hel3Exp);

  fhBB->Fill(track->GetTPCmomentum()*track->Charge(),track->GetTPCsignal());
  
  if (TMath::Abs(nSigmaHe3)>4.0) return false;
  /*
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTPC) > 4.0) return false;
 */
  return true;
}
//_______________________________________________________________________________________
Double_t AliAnalysisTaskHypv2PbPb18::InvariantMassHypertriton (TVector3 P1, TVector3 P2)  {
    
    //Mass Daughters
    Double_t m3He  = AliPID::ParticleMass (AliPID::kHe3);
    Double_t mpion = AliPID::ParticleMass (AliPID::kPion);
    
    //Invariant Mass Calculation
    TVector3 P = P1 + P2;
    
    Double_t E1 = TMath::Sqrt(m3He*m3He + P1.Mag2());
    Double_t E2 = TMath::Sqrt(mpion*mpion + P2.Mag2());
    Double_t m  = TMath::Sqrt( (E1+E2)*(E1+E2) - P.Mag2() );
    
    return m;
}
//_______________________________________________________________________________________
Double_t AliAnalysisTaskHypv2PbPb18::GetDecayLengthV0 (AliESDv0 *V0)  {
    
    
    //Initialization
    Double_t decayLengthV0 = 0;
    
    //Secondary Vertex Position
    Double_t secVertex[3] = { 0.0, 0.0, 0.0 };
    V0->GetXYZ(secVertex[0],secVertex[1],secVertex[2]);
    
    //Primary Vertex Position
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t primVertex[3] = { 0.0, 0.0, 0.0 };
    vertex->GetXYZ(primVertex);
    
    //Decay Length
    Double_t Dx = primVertex[0]-secVertex[0];
    Double_t Dy = primVertex[1]-secVertex[1];
    Double_t Dz = primVertex[2]-secVertex[2];
    decayLengthV0 = TMath::Sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
    
    return decayLengthV0;
}
//_______________________________________________________________________________________
Double_t AliAnalysisTaskHypv2PbPb18::GetTransverseDCA (AliESDtrack *track)  {
    
    /*
     Double_t impactParameter[2];
     track -> GetImpactParameters(impactParameter[0],impactParameter[1]);
     */
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_______________________________________________________________________________________
Double_t AliAnalysisTaskHypv2PbPb18::GetDCA (AliESDtrack *track)  {
    
    /*
     Double_t impactParameter[2];
     track -> GetImpactParameters(impactParameter[0],impactParameter[1]);
     */
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCA = TMath::Sqrt(impactParameter[0]*impactParameter[0]+impactParameter[1]*impactParameter[1]);
    
    return DCA;
}
//________________________________________________________________________
Float_t AliAnalysisTaskHypv2PbPb18::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
   return result;
}


//_____________________________________________________________________________
AliAnalysisTaskHypv2PbPb18::~AliAnalysisTaskHypv2PbPb18()
{

  
}

//______________________________________________________________________________
void AliAnalysisTaskHypv2PbPb18::UserCreateOutputObjects()
{ 
  //-------------------------------------------------------
  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!
  
  if(! fHistEventMultiplicity ){

    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 12 , 0.5,12.5);
    
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/good PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Events wo pileup");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"kINT7");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"HM V0");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"HM SPD");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"kCentral");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"kSemiCentral");

    fListHist->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity  = new TH2F( "fHistTrackMultiplicity", "Nb of Tracks MB Events |Vz| < 10", 250,0, 5000,105,0,105);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicity);
  } 
 
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBB);
  }
  
  if(! fhBBHyp ){
    fhBBHyp = new TH2F( "fhBBHyp" , "BetheBlochTPC - Hypertriton" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBBHyp);
  }

  if(! fhBBAHyp ){
    fhBBAHyp = new TH2F( "fhBBAHyp" , "BetheBlochTPC - AHypertriton" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBBAHyp);
  }
 
  
  EPVzAvsCentrality  = new TH2D("EPVzAvsCentrality" , "EPVzAvsCentrality" , 80,-TMath::Pi(),TMath::Pi(), 105,0,105);
  EPVzCvsCentrality  = new TH2D("EPVzCvsCentrality" , "EPVzCvsCentrality" , 80,-TMath::Pi(),TMath::Pi(), 105,0,105);

  fListHist->Add(EPVzAvsCentrality);
  fListHist->Add(EPVzCvsCentrality);
  

  if(fNHarm < 3){
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",1000,-100,100,105,0,105);
    hQVzAQTPCvsCentrality = new TH2F("hQVzAQTPCvsCentrality","hQVzAQTPCvsCentrality",1000,-100,100,105,0,105);
    hQVzCQTPCvsCentrality = new TH2F("hQVzCQTPCvsCentrality","hQVzCQTPCvsCentrality",1000,-100,100,105,0,105);
  }
  else{
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",5000,-5000,5000,105,0,105);
    hQVzAQTPCvsCentrality = new TH2F("hQVzAQTPCvsCentrality","hQVzAQTPCvsCentrality",5000,-5000,5000,105,0,105);
    hQVzCQTPCvsCentrality = new TH2F("hQVzCQTPCvsCentrality","hQVzCQTPCvsCentrality",5000,-5000,5000,105,0,105);
  }
  fListHist->Add(hQVzAQVzCvsCentrality);
  fListHist->Add(hQVzAQTPCvsCentrality);
  fListHist->Add(hQVzCQTPCvsCentrality);
  
  if(fNHarm < 3){
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",100,-20,20,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",100,-20,20,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",100,-20,20,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",100,-20,20,105,0,105);
  }
  
  else{
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",2000,-500,500,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",2000,-500,500,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",2000,-500,500,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",2000,-500,500,105,0,105);
  }

  fListHist->Add(hQxVzAvsCentrality);
  fListHist->Add(hQyVzAvsCentrality);
  fListHist->Add(hQxVzCvsCentrality);
  fListHist->Add(hQyVzCvsCentrality);
 
   hCos2DeltaTPCVzAvsCentrality   = new TH2F("hCos2DeltaTPCVzAvsCentrality"  ,"hCos2DeltaTPCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaTPCVzCvsCentrality   = new TH2F("hCos2DeltaTPCVzCvsCentrality"  ,"hCos2DeltaTPCVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzAVzCvsCentrality   = new TH2F("hCos2DeltaVzAVzCvsCentrality"  ,"hCos2DeltaVzAVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzATPCvsCentrality   = new TH2F("hCos2DeltaVzATPCvsCentrality"  ,"hCos2DeltaVzATPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCTPCvsCentrality   = new TH2F("hCos2DeltaVzCTPCvsCentrality"  ,"hCos2DeltaVzCTPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCVzAvsCentrality   = new TH2F("hCos2DeltaVzCVzAvsCentrality"  ,"hCos2DeltaVzCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);

  fListHist->Add(hCos2DeltaTPCVzAvsCentrality);
  fListHist->Add(hCos2DeltaTPCVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzAVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzATPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCTPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCVzAvsCentrality);


  if(!ftree){

    //Reduced Tree HyperTriton
    ftree = new TTree("ftree","ftree");
    ftree -> Branch("iEvent",&iEvent,"iEvent/I");
    ftree -> Branch("zVertex",&zVertex,"zVertex/D");
    ftree -> Branch("centrality",&centrality,"centrality/D");
    ftree -> Branch("px_Daughter1",&px_Daughter1,"px_Daughter1/D");
    ftree -> Branch("py_Daughter1",&py_Daughter1,"py_Daughter1/D");
    ftree -> Branch("pz_Daughter1",&pz_Daughter1,"pz_Daughter1/D");
    ftree -> Branch("q_Daughter1",&q_Daughter1,"q_Daughter1/I");
    ftree -> Branch("dcaxy_Daughter1",&dcaxy_Daughter1,"dcaxy_Daughter1/D");
    ftree -> Branch("nTPC_Clusters_Daughter1",&nTPC_Clusters_Daughter1,"nTPC_Clusters_Daughter1/I");
    ftree -> Branch("nTPC_Clusters_dEdx_Daughter1",&nTPC_Clusters_dEdx_Daughter1,"nTPC_Clusters_dEdx_Daughter1/I");
    ftree -> Branch("chi2_TPC_Daughter1",&chi2_TPC_Daughter1,"chi2_TPC_Daughter1/D");
    ftree -> Branch("nSigmaTPC_He3_Daughter1",&nSigmaTPC_He3_Daughter1,"nSigmaTPC_He3_Daughter1/D");
    ftree -> Branch("nSigmaTPC_Pion_Daughter1",&nSigmaTPC_Pion_Daughter1,"nSigmaTPC_Pion_Daughter1/D");
    ftree -> Branch("px_Daughter2",&px_Daughter2,"px_Daughter2/D");
    ftree -> Branch("py_Daughter2",&py_Daughter2,"py_Daughter2/D");
    ftree -> Branch("pz_Daughter2",&pz_Daughter2,"pz_Daughter2/D");
    ftree -> Branch("q_Daughter2",&q_Daughter2,"q_Daughter2/I");
    ftree -> Branch("dcaxy_Daughter2",&dcaxy_Daughter2,"dcaxy_Daughter2/D");
    ftree -> Branch("nTPC_Clusters_Daughter2",&nTPC_Clusters_Daughter2,"nTPC_Clusters_Daughter2/I");
    ftree -> Branch("nTPC_Clusters_dEdx_Daughter2",&nTPC_Clusters_dEdx_Daughter2,"nTPC_Clusters_dEdx_Daughter2/I");
    ftree -> Branch("chi2_TPC_Daughter2",&chi2_TPC_Daughter2,"chi2_TPC_Daughter2/D");
    ftree -> Branch("nSigmaTPC_He3_Daughter2",&nSigmaTPC_He3_Daughter2,"nSigmaTPC_He3_Daughter2/D");
    ftree -> Branch("nSigmaTPC_Pion_Daughter2",&nSigmaTPC_Pion_Daughter2,"nSigmaTPC_Pion_Daughter2/D");
    ftree -> Branch("isOnTheFlyV0",&isOnTheFlyV0,"isOnTheFlyV0/I");
    ftree -> Branch("cosPointingAngle",&cosPointingAngle,"cosPointingAngle/D");
    ftree -> Branch("dcaV0Daughters",&dcaV0Daughters,"dcaV0Daughters/D");
    ftree -> Branch("radius",&radius,"radius/D");
    ftree -> Branch("chi2V0",&chi2V0,"chi2V0/D");
    ftree -> Branch("decayLength",&decayLength,"decayLength/D");
    ftree -> Branch("deltaphiV0A"           ,&deltaphiV0A           ,"deltaphiV0A/D"         );
    ftree -> Branch("deltaphiV0C"           ,&deltaphiV0C           ,"deltaphiV0C/D"         );
    ftree -> Branch("uqV0A"  ,&uqV0A  ,"uqV0A/D");
    ftree -> Branch("uqV0C"  ,&uqV0C  ,"uqV0C/D");

    // fOutputList -> Add(ftree);
  }

  fEventCuts.AddQAplotsToList(fListHist);

  PostData(1,  fListHist);
  PostData(2,  ftree);  

}

//______________________________________________________________________________
void AliAnalysisTaskHypv2PbPb18::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
 
  Info("AliAnalysisTaskHypv2PbPb18","Starting UserExec");  
  fHistEventMultiplicity->Fill(1);
  AliVEvent *event = InputEvent();

  fESDevent = dynamic_cast<AliESDEvent*>(event);

  if (!fESDevent) {
    AliError("Cannot get the AOD event");
      return;
  }  
  fevent = fESDevent;

  if(!fevent || !fevent->GetHeader()){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }
      
  Int_t run = fevent->GetRunNumber();

  //cout<<"-------------------->RUN number "<<run<<endl;
  
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);
    fRun = run;
  }
  
  //  cout<<"Done info calib"<<endl;

  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fevent)) {
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }
  
  fHistEventMultiplicity->Fill(2);

  // Primary vertex cut
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};

  const AliVVertex* vertexmain = fevent->GetPrimaryVertex();
  if (!vertexmain){
    AliWarning("No prim. vertex in ESD... return!");
    
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }

  //  cout<<"Open vertex"<<endl;

  
  vertexmain->GetXYZ( lBestPrimaryVtxPos );
 
  if((TMath::Abs(lBestPrimaryVtxPos[2])) > fVzmax) return;
  fHistEventMultiplicity->Fill(3);

  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fESDevent->IsPileupFromSPD();
  
  if(isPileUpSpd){  
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }
  
  //  cout<<"is pileup"<<endl;

  fHistEventMultiplicity->Fill(4); // analyzed events with PV w/o pile up
  
  /*
  //event cuts
  // 1. primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fevent->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1) return;
  fHistEventMultiplicity->Fill(2);
 
  // 2. SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fevent->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return;
  fHistEventMultiplicity->Fill(3);
  
  // 3. pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(fevent);

  if(isPileupFromMV)return;
  fHistEventMultiplicity->Fill(4);
  
  // 4. cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) >  fVzmax)
    return;
  fHistEventMultiplicity->Fill(5);
  */
  
  // 5. Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  

  Bool_t isTriggerSelected = kFALSE;
 
  Bool_t isSelectedINT7        = fSelectMask& AliVEvent::kINT7;
  Bool_t isSelectedCentral     = fSelectMask& AliVEvent::kCentral;
  Bool_t isSelectedSemiCentral = fSelectMask& AliVEvent::kSemiCentral;
  Bool_t isSelectedHMV0        = fSelectMask& AliVEvent::kHighMultV0;
  Bool_t isSelectedHMSPD       = fSelectMask& AliVEvent::kHighMultSPD;
  
  //Int_t eventtype = -999;
  if(isSelectedINT7){
    eventtype = 1;
    fHistEventMultiplicity->Fill(6);
  }
  if(isSelectedHMV0){
    eventtype = 2;
    fHistEventMultiplicity->Fill(7);
  }
  if(isSelectedHMSPD){
    eventtype = 3;
    fHistEventMultiplicity->Fill(8);
  }
  if(isSelectedCentral){
    eventtype = 4;
    fHistEventMultiplicity->Fill(9);
  }
  if(isSelectedSemiCentral){
    eventtype = 5;
    fHistEventMultiplicity->Fill(10);
  }

  //  if(eventtype!=1 && eventtype!=2 && eventtype!=3)return;
  if(eventtype!=1 && eventtype!=4 && eventtype!=5)return;
  
  //  cout<<"event selected"<<endl;

  // get the PID response
  fPIDResponse=inputHandler->GetPIDResponse(); 
 

  //Analysis
  //cout<<"Start event"<<endl;
  Analyze(fevent,lBestPrimaryVtxPos[2],eventtype);
  //  Analyze(fevent);
     
}

//________________________________________________________________________
void AliAnalysisTaskHypv2PbPb18::Analyze(AliVEvent* esd, Double_t vz, Int_t evttype)
//void AliAnalysisTaskHypv2PbPb18::Analyze(AliVEvent* esd)
{  
  //  cout<<"OPen analysize"<<endl;

  //new vertex selection
  const AliESDVertex* vtTrc = (AliESDVertex*)esd->GetPrimaryVertex();
  const AliESDVertex* vtSPD = (AliESDVertex*)esd->GetPrimaryVertexSPD();

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  
  double dz = vtTrc->GetZ() - vtSPD->GetZ();
  
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;
    
  if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
    return; // bad vertexing
  
  //Centrality
  Float_t v0Centr    = -100.;
  Float_t cl1Centr   = -100.;
  Float_t cl0Centr   = -100.;
    
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*) esd->FindListObject("MultSelection");
  if( !MultSelection) {
    AliWarning("AliMultSelection object not found!");
    return;
  } else {
    v0Centr  = MultSelection->GetMultiplicityPercentile("V0M");
    cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
    cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
  }
  
  //  cout<<"v0Centr: "<<v0Centr<<endl;
  
  if (v0Centr >= 90. || v0Centr < 0)
    return; 
  
  AliESDVZERO* esdV0 = (AliESDVZERO*)esd->GetVZEROData();
  Float_t multV0a = esdV0->GetMTotV0A();
  Float_t multV0c = esdV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  //centrality estimator : 0 = V0M, 1 = CL1, 2 = CL2
  if (fCenCalV0 == 1)
    v0Centr = cl0Centr;
  else if (fCenCalV0 == 2)
    v0Centr = cl1Centr;

  /*
    Short_t centrCode = -10;
    if ((v0Centr >= 0) && (v0Centr < 5.))
    centrCode = 0;
    else if ((v0Centr >= 5.) && (v0Centr < 10.))
    centrCode = 1;
    else if ((v0Centr >= 10.) && (v0Centr < 20.))
    centrCode = 2;
    else if ((v0Centr >= 20.) && (v0Centr < 30.))
    centrCode = 3;
    else if ((v0Centr >= 30.) && (v0Centr < 40.))
    centrCode = 4;
    else if ((v0Centr >= 40.) && (v0Centr < 50.))
    centrCode = 5;
    else if ((v0Centr >= 50.) && (v0Centr < 60.))
    centrCode = 6;
    else if ((v0Centr >= 60.) && (v0Centr < 70.))
    centrCode = 7;
    else if ((v0Centr >= 70.) && (v0Centr < 80.))
    centrCode = 8;
    else if ((v0Centr >= 80.) && (v0Centr < 90.))
    centrCode = 9;
    
    if (centrCode < 0)
    return;
  */

  Int_t iCen = Int_t(v0Centr);
  
  // Int_t iCentSPD = Int_t(cl1Centr);
  // if (iCentSPD >= 90)
  //   return;
    
  // Int_t iCen = Int_t(v0Centr);
  // if (iCen >= 90)
  //   return;
  
  fHistEventMultiplicity->Fill(11);
  
  //V0 info
  Double_t Qxan = 0, Qyan = 0;
  Double_t Qxcn = 0, Qycn = 0;
  Double_t sumMa = 0, sumMc = 0;
  
  //  cout<<"qui"<<endl;

  for (Int_t iV0 = 0; iV0 < 64; iV0++) {
    
    Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
    Float_t multv0 = esdV0->GetMultiplicity(iV0);
    
    if (iV0 < 32){
      
      Double_t multCorC = -10;
      
      if (iV0 < 8)
	multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
      else if (iV0 >= 8 && iV0 < 16)
	multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
      else if (iV0 >= 16 && iV0 < 24)
	multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
      else if (iV0 >= 24 && iV0 < 32)
	multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);
      
      if (multCorC < 0){
	cout<<"Problem with multiplicity in V0C"<<endl;
	continue;
      }
      
      Qxcn += TMath::Cos(fNHarm*phiV0) * multCorC;
      Qycn += TMath::Sin(fNHarm*phiV0) * multCorC;
      
      
      // if (fIsQ2Ana){
      //     Qxc2ese += TMath::Cos(2.*phiV0) * multCorC;
      //     Qyc2ese += TMath::Sin(2.*phiV0) * multCorC;
      // } else {
      //     Qxc3ese += TMath::Cos(3.*phiV0) * multCorC;
      //     Qyc3ese += TMath::Sin(3.*phiV0) * multCorC;
      // }
      
      
      sumMc = sumMc + multCorC;
      
    } else {
      
      Double_t multCorA = -10;
      
      if (iV0 >= 32 && iV0 < 40)
	multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
      else if (iV0 >= 40 && iV0 < 48)
	multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
      else if (iV0 >= 48 && iV0 < 56)
	multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
      else if (iV0 >= 56 && iV0 < 64)
	multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);
      
      if (multCorA < 0){
	cout<<"Problem with multiplicity in V0A"<<endl;
	continue;
      }
      
      Qxan += TMath::Cos(fNHarm*phiV0) * multCorA;
      Qyan += TMath::Sin(fNHarm*phiV0) * multCorA;
      
      sumMa = sumMa + multCorA;
      
    }
    
  }

  //cout<<"qua"<<endl;

  if (sumMa < 0 || sumMc < 0)
    return;

  Double_t QxanCor = Qxan;
  Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen+1))/fQynsV0A->GetBinContent(iCen+1);
  Double_t QxcnCor = Qxcn;
  Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen+1))/fQynsV0C->GetBinContent(iCen+1);
  
  if (fNHarm != 4.){
    QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen+1))/fQxnsV0A->GetBinContent(iCen+1);
    QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen+1))/fQxnsV0C->GetBinContent(iCen+1);
  }

   
  // cout<<"quo"<<endl;
  
  // cout<<"iCen+1 "<<iCen+1 <<endl;
  
  
  // Double_t QxanCor = Qxan;
  // Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen+1))/fQynsV0A->GetBinContent(iCen+1);
  // Double_t QxcnCor = Qxcn;
  // Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen+1))/fQynsV0C->GetBinContent(iCen+1);
  
  // if (fNHarm != 4.){
  //   QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen+1))/fQxnsV0A->GetBinContent(iCen+1);
  //   QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen+1))/fQxnsV0C->GetBinContent(iCen+1);
  // }
  
  Double_t evPlAngV0A = TMath::ATan2(QyanCor, QxanCor)/fNHarm;
  Double_t evPlAngV0C = TMath::ATan2(QycnCor, QxcnCor)/fNHarm;
 
  EPVzAvsCentrality  ->Fill(evPlAngV0A  , iCen); 
  EPVzCvsCentrality  ->Fill(evPlAngV0C  , iCen); 

  //  cout<<"evPlAngV0C"<<  evPlAngV0C<<endl;
  //Qtpc --RAMONA
    
  
  const Int_t nTracks = esd->GetNumberOfTracks();
  // cout<<"TPC ev plane "<<nTracks<<endl;
  Double_t Qxtn = 0, Qytn = 0;
  
  for (Int_t it1 = 0; it1 < nTracks; it1++) {
    AliESDtrack* esdTrk1 = (AliESDtrack*)esd->GetTrack(it1);
    
    if (!esdTrk1){
      delete esdTrk1;
      continue;
    }
    
    if(!fESDtrackCutsEP->AcceptTrack((AliESDtrack*)esdTrk1))
      continue; //tpc only track

    if(!PassedBasicTrackQualityCuts_Pos (esdTrk1))continue;

    if (TMath::Abs(esdTrk1->Eta()) < 0.8 && esdTrk1->GetTPCNcls() >= 70 && esdTrk1->Pt() >= 0.2 && esdTrk1->Pt() < 3.){
      Qxtn += TMath::Cos(fNHarm*esdTrk1->Phi());
      Qytn += TMath::Sin(fNHarm*esdTrk1->Phi());
    }
  }
  
  //  cout<<"Qxtn: "<<Qxtn<<endl;
  // TBC
  Double_t evPlAngTPC = TMath::ATan2(Qytn, Qxtn)/fNHarm;
  
  hCos2DeltaTPCVzAvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngTPC - evPlAngV0A)) , iCen);
  hCos2DeltaTPCVzCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngTPC - evPlAngV0C)) , iCen);
  hCos2DeltaVzAVzCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0A - evPlAngV0C)) , iCen);
  hCos2DeltaVzATPCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0A - evPlAngTPC)) , iCen);
  hCos2DeltaVzCTPCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0C - evPlAngTPC)) , iCen);
  hCos2DeltaVzCVzAvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0C - evPlAngV0A)) , iCen);

  //Scalar Product -- Resolutions
  
  Double_t corV0AV0Cvn = QxanCor*QxcnCor + QyanCor*QycnCor;
  Double_t corV0ATPCvn = QxanCor*Qxtn + QyanCor*Qytn;
  Double_t corV0CTPCvn = QxcnCor*Qxtn + QycnCor*Qytn;
  
  // cout<<"corV0AV0Cvn "<<corV0AV0Cvn <<endl;
  // cout<<"corV0ATPCvn "<<corV0ATPCvn <<endl;
  // cout<<"corV0CTPCvn "<<corV0CTPCvn <<endl;
  
  hQVzAQVzCvsCentrality->Fill(corV0AV0Cvn,iCen);
  hQVzAQTPCvsCentrality->Fill(corV0ATPCvn,iCen);
  hQVzCQTPCvsCentrality->Fill(corV0CTPCvn,iCen);
  
  // cout<<"corV0AV0Cvn "<<corV0AV0Cvn<<endl;

  //NUA correction
 
  hQxVzAvsCentrality->Fill(QxanCor,iCen);
  hQyVzAvsCentrality->Fill(QyanCor,iCen);
  hQxVzCvsCentrality->Fill(QxcnCor,iCen);
  hQyVzCvsCentrality->Fill(QycnCor,iCen);

  //----------------------------------------------------
  // from here my analysis starts
  // Hypertriton loop

  Int_t TrackNumber = esd->GetNumberOfTracks();
  fHistTrackMultiplicity->Fill(TrackNumber,iCen); //tracce per evento

  //Loop Over Reconstructed V0s
  for ( Int_t iV0=0 ; iV0<fESDevent->GetNumberOfV0s() ; iV0++ ) {
        
    //Get V0 Candidate
    AliESDv0 *V0 = (AliESDv0*)fESDevent->GetV0(iV0);
    if (!V0) continue;
    if ( V0->GetOnFlyStatus()) isOnTheFlyV0=1;//V0-Online
    if (!V0->GetOnFlyStatus()) isOnTheFlyV0=0;//V0-Offline
    if (!PassedMinimalQualityCutsV0(V0)) continue;
        
    //Get V0 Daughters
    AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    if (!posTrack) continue;
    if (!negTrack) continue;
    if (posTrack->Charge() == negTrack->Charge()) continue;
    if (posTrack->GetID()  == negTrack->GetID() ) continue;
    
    //Quality Requirements
    if (!PassedBasicTrackQualityCuts_Pos (posTrack)) continue;
    if (!PassedBasicTrackQualityCuts_Neg (negTrack)) continue;
    
    //Hypertriton Candidate Selection
    if (!IsHyperTritonCandidate(V0)) continue;
    
    
    //Momentum Components of V0 Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);

    iEvent     = evttype;        
    zVertex    = vz;
    centrality = iCen;

    //Daughter1 (Positive Charge)
    px_Daughter1                    = posMomentum[0];
    py_Daughter1                    = posMomentum[1];
    pz_Daughter1                    = posMomentum[2];
    q_Daughter1                     = (Int_t) posTrack -> Charge();
    dcaxy_Daughter1                 = GetTransverseDCA (posTrack);
    nTPC_Clusters_Daughter1         = posTrack -> GetTPCNcls();
    nTPC_Clusters_dEdx_Daughter1    = posTrack -> GetTPCsignalN();
    chi2_TPC_Daughter1              = posTrack -> GetTPCchi2();
    nSigmaTPC_He3_Daughter1         = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kHe3);
    nSigmaTPC_Pion_Daughter1        = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kPion);
    
    //Daughter2  (Negative Charge)
    px_Daughter2                    = negMomentum[0];
    py_Daughter2                    = negMomentum[1];
    pz_Daughter2                    = negMomentum[2];
    q_Daughter2                     = (Int_t) negTrack -> Charge();
    dcaxy_Daughter2                 = GetTransverseDCA (negTrack);
    nTPC_Clusters_Daughter2         = negTrack -> GetTPCNcls();
    nTPC_Clusters_dEdx_Daughter2    = negTrack -> GetTPCsignalN();
    chi2_TPC_Daughter2              = negTrack -> GetTPCchi2();
    nSigmaTPC_He3_Daughter2         = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kHe3);
    nSigmaTPC_Pion_Daughter2        = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kPion);
    
    
    //Pair Variables
    cosPointingAngle = V0->GetV0CosineOfPointingAngle();
    dcaV0Daughters   = V0->GetDcaV0Daughters();
    radius           = V0->GetRr();
    chi2V0           = V0->GetChi2V0();
    decayLength      = GetDecayLengthV0 (V0);

    //flow variables
    deltaphiV0A=GetPhi0Pi(V0->Phi()-evPlAngV0A);
    deltaphiV0C=GetPhi0Pi(V0->Phi()-evPlAngV0C);
    
    // Scalar Product
    uqV0A = TMath::Cos(fNHarm*V0->Phi())*QxanCor+TMath::Sin(fNHarm*V0->Phi())*QyanCor;
    uqV0C = TMath::Cos(fNHarm*V0->Phi())*QxcnCor+TMath::Sin(fNHarm*V0->Phi())*QycnCor;

    
    //Fill Reduced Tree
    ftree -> Fill();
  }
  
  PostData(1, fListHist);
  PostData(2, ftree);
  
}

//_____________________________________________________________________________
void AliAnalysisTaskHypv2PbPb18::OpenInfoCalbration(Int_t run )
{

  //foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibpPb2016/calibV0NoSDD.root");
  
  TFile* foadb = 0;
  if (!gGrid) TGrid::Connect("alien");
  if (fPeriod == 0)
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibPbPb2018/calibV0Run2Vtx10P118q.root");
  //TFile::Open("calibV0Run2Vtx10P118q.root");
  else
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibPbPb2018/calibV0Run2Vtx10P118r.root");
  //foadb = TFile::Open("calibV0Run2Vtx10P118r.root");
  
  
  if(!foadb){
    printf("OADB V0 calibration file cannot be opened\n");
    return;
  }
    
  
  AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return;
  }
  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
    return;
  }
  fMultV0 = ((TH1D*) cont->GetObject(run));
    

    
  AliOADBContainer* contQxnam = 0;
  if (fNHarm == 2.)
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
  else
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");
    
  if(!contQxnam){
    printf("OADB object fqxanm is not available in the file\n");
    return;
  }
  if(!(contQxnam->GetObject(run))){
    printf("OADB object fqxanm is not available for run %i\n", run);
    return;
  }
  fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));
    
    
    
  AliOADBContainer* contQynam = 0;
  if (fNHarm == 2.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
  else if (fNHarm == 3.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
  else if (fNHarm == 4.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya4m");
    
  if(!contQynam){
    printf("OADB object fqyanm is not available in the file\n");
    return;
  }
  if(!(contQynam->GetObject(run))){
    printf("OADB object fqyanm is not available for run %i\n", run);
    return;
  }
  fQynmV0A = ((TH1D*) contQynam->GetObject(run));
    
    
    
  AliOADBContainer* contQxnas = 0;
  if (fNHarm == 2.)
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
  else
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");
    
  if(!contQxnas){
    printf("OADB object fqxans is not available in the file\n");
    return;
  }
  if(!(contQxnas->GetObject(run))){
    printf("OADB object fqxans is not available for run %i\n", run);
    return;
  }
  fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));
    
    
    
  AliOADBContainer* contQynas = 0;
  if (fNHarm == 2.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
  else if (fNHarm == 3.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
  else if (fNHarm == 4.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya4s");
    
  if(!contQynas){
    printf("OADB object fqyans is not available in the file\n");
    return;
  }
  if(!(contQynas->GetObject(run))){
    printf("OADB object fqyans is not available for run %i\n", run);
    return;
  }
  fQynsV0A = ((TH1D*) contQynas->GetObject(run));
    
    
    
  AliOADBContainer* contQxncm = 0;
  if (fNHarm == 2.)
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
  else
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");
    
  if(!contQxncm){
    printf("OADB object fqxcnm is not available in the file\n");
    return;
  }
  if(!(contQxncm->GetObject(run))){
    printf("OADB object fqxcnm is not available for run %i\n", run);
    return;
  }
  fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));
    
    
    
  AliOADBContainer* contQyncm = 0;
  if (fNHarm == 2.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
  else if (fNHarm == 3.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
  else if (fNHarm == 4.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");
    
  if(!contQyncm){
    printf("OADB object fqyc2m is not available in the file\n");
    return;
  }
  if(!(contQyncm->GetObject(run))){
    printf("OADB object fqyc2m is not available for run %i\n", run);
    return;
  }
  fQynmV0C = ((TH1D*) contQyncm->GetObject(run));
    


  AliOADBContainer* contQxncs = 0;
  if (fNHarm == 2.)
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
  else
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");
    
  if(!contQxncs){
    printf("OADB object fqxc2s is not available in the file\n");
    return;
  }
  if(!(contQxncs->GetObject(run))){
    printf("OADB object fqxc2s is not available for run %i\n", run);
    return;
  }
  fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));
    
    
    
  AliOADBContainer* contQyncs = 0;
  if (fNHarm == 2.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
  else if (fNHarm == 3.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
  else if (fNHarm == 4.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");
    
  if(!contQyncs){
    printf("OADB object fqycnm is not available in the file\n");
    return;
  }
  if(!(contQyncs->GetObject(run))){
    printf("OADB object fqycns is not available for run %i\n", run);
    return;
  }
  fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
        
        
        
  /*    
	AliOADBContainer* contQx2cmESE = (AliOADBContainer*) foadb->Get("fqxc2m");
	if(!contQx2cmESE){
        printf("OADB object fqxc2m is not available in the file\n");
        return;
	}
	if(!(contQx2cmESE->GetObject(run))){
        printf("OADB object fqxc2m is not available for run %i\n", run);
        return;
	}
	fQx2mV0CESE = ((TH1D*) contQx2cmESE->GetObject(run));
        
        
	AliOADBContainer* contQy2cmESE = (AliOADBContainer*) foadb->Get("fqyc2m");
	if(!contQy2cmESE){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
	}
	if(!(contQy2cmESE->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
	}
	fQy2mV0CESE = ((TH1D*) contQy2cmESE->GetObject(run));
        
        
        
        
	AliOADBContainer* contQx3cmESE = (AliOADBContainer*) foadb->Get("fqxc3m");
	if(!contQx3cmESE){
        printf("OADB object fqxc3m is not available in the file\n");
        return;
	}
	if(!(contQx3cmESE->GetObject(run))){
        printf("OADB object fqxc3m is not available for run %i\n", run);
        return;
	}
	fQx3mV0CESE = ((TH1D*) contQx3cmESE->GetObject(run));
        
        
	AliOADBContainer* contQy3cmESE = (AliOADBContainer*) foadb->Get("fqyc3m");
	if(!contQy3cmESE){
        printf("OADB object fqyc3m is not available in the file\n");
        return;
	}
	if(!(contQy3cmESE->GetObject(run))){
        printf("OADB object fqyc3m is not available for run %i\n", run);
        return;
	}
	fQy3mV0CESE = ((TH1D*) contQy3cmESE->GetObject(run));
  */

}
//_____________________________________________________________________________
void AliAnalysisTaskHypv2PbPb18::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
