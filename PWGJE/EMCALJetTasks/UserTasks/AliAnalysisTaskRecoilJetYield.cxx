//
// My analysis task
//
// Author: Harry Andrews

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
//#include "AliPythiaInfo.h"
#include "TRandom3.h"
#include "AliPicoTrack.h"
#include "AliEmcalJetFinder.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskRecoilJetYield.h"
#include "AliMultSelection.h"

//Globals



Double_t Eta_Upper=1.00;
Double_t Eta_Lower=-1.00;
Int_t Eta_Bins_1=100;
Double_t Phi_Upper=2*(TMath::Pi());
Double_t Phi_Lower=(-1*(TMath::Pi()));
Int_t Phi_Bins_1=100;
Float_t Percentile_1 = 300; 




using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskRecoilJetYield)

//________________________________________________________________________
AliAnalysisTaskRecoilJetYield::AliAnalysisTaskRecoilJetYield() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRecoilJetYield", kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fPtMinTriggerHadron(20.),
  fPtMaxTriggerHadron(50.),
  fRecoilAngularWindow(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
  fBeta_SD(0),
  fZCut(0.1),
  fRho(1e-6),
  fRhoParam(0),
  fNsubMeasure(kFALSE),
  fDoSoftDrop(kFALSE),
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetCounter(0x0),
  fhNumberOfJetTracks(0x0),
  fTreeJetInfo(0),
  fhJetArea(0x0),
  fhTrackPt(0x0),
  fhGroomedPtvJetPt(0x0),
  fhPtTriggerHadron(0x0), 
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0)

{
  for(Int_t i=0;i<nBranch;i++){
    fJetInfoVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskRecoilJetYield::AliAnalysisTaskRecoilJetYield(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fPtMinTriggerHadron(20.),
  fPtMaxTriggerHadron(50.),
  fRecoilAngularWindow(0.6),
  fSemigoodCorrect(0),
  fHolePos(0),
  fHoleWidth(0), 
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fSubJetAlgorithm(0),
  fSubJetRadius(0.1),
  fSubJetMinPt(1),
  fJetRadius(0.4),
  fRMatched(0.2),
  fSharedFractionPtMin(0.5),
  fDerivSubtrOrder(0),
  fFullTree(kFALSE),
  fBeta_SD(0),
  fZCut(0.1),
  fRho(1e-6),
  fRhoParam(0),
  fNsubMeasure(kFALSE),
  fDoSoftDrop(kFALSE),
  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),
  fhJetMass(0x0),
  fhJetRadius(0x0),
  fhJetCounter(0x0),
  fhNumberOfJetTracks(0x0),
  fTreeJetInfo(0),
  fhJetArea(0x0),
  fhTrackPt(0x0),
  fhGroomedPtvJetPt(0x0),
  fhPtTriggerHadron(0x0), 
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0)
  
{
  // Standard constructor.
  for(Int_t i=0;i<nBranch;i++){
    fJetInfoVar[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskRecoilJetYield::~AliAnalysisTaskRecoilJetYield()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskRecoilJetYield::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJetInfo = new TTree(nameoutput, nameoutput);
  
  if (!fFullTree){
    const Int_t nVarMin = 26;
    TString *fJetInfoVarNames = new TString [nVarMin];
  
    fJetInfoVarNames[0] = "Pt";
    fJetInfoVarNames[1] = "Pt_Truth";
    fJetInfoVarNames[2] = "Phi";
    fJetInfoVarNames[3] = "Phi_Truth";
    fJetInfoVarNames[4] = "Eta";
    fJetInfoVarNames[5] = "Eta_Truth";
    fJetInfoVarNames[6] = "SymParam";
    fJetInfoVarNames[7] = "SymParam_Truth";
    fJetInfoVarNames[8] = "Tau1";
    fJetInfoVarNames[9] = "Tau1_Truth";
    fJetInfoVarNames[10] = "Tau2";
    fJetInfoVarNames[11] = "Tau2_Truth";
    fJetInfoVarNames[12] = "Mass";
    fJetInfoVarNames[13] = "Mass_Truth";
    fJetInfoVarNames[14] = "NTracks";
    fJetInfoVarNames[15] = "NTracks_Truth";
    fJetInfoVarNames[16] = "JetPtNoCut";
    fJetInfoVarNames[17] = "JetPtNoCut_Truth";
    fJetInfoVarNames[18] = "PTD";
    fJetInfoVarNames[19] = "PTD_Truth";
    fJetInfoVarNames[20] = "Angularity";
    fJetInfoVarNames[21] = "Angularity_Truth";
    fJetInfoVarNames[22] = "Centrality";
    fJetInfoVarNames[23] = "Centrality_Truth";
    fJetInfoVarNames[24] = "DelR";
    fJetInfoVarNames[25] = "DelR_Truth";

  
    
    for(Int_t ivar=0; ivar < nVarMin; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeJetInfo->Branch(fJetInfoVarNames[ivar].Data(), &fJetInfoVar[ivar], Form("%s/D", fJetInfoVarNames[ivar].Data()));
    }
  }

  if (fJetShapeType==AliAnalysisTaskRecoilJetYield::kData){
    
    fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );   
    fOutput->Add(fhJetPt);
    fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
    fOutput->Add(fhJetPhi);
    fhJetEta= new TH1F("fhJetEta", "Jet Eta", Eta_Bins_1, Eta_Lower, Eta_Upper);
    fOutput->Add(fhJetEta);
    fhJetMass= new TH1F("fhJetMass", "Jet Mass", 4000,-0.5, 39.5);
    fOutput->Add(fhJetMass);
    fhJetRadius= new TH1F("fhJetRadius", "Jet Radius", 100, -0.05,0.995);
    fOutput->Add(fhJetRadius);
    fhNumberOfJetTracks= new TH1F("fhNumberOfJetTracks", "Number of Tracks within a Jet", 300, -0.5,299.5);
    fOutput->Add(fhNumberOfJetTracks);
    fhJetCounter= new TH1F("fhJetCounter", "Jet Counter", 150, -0.5, 149.5);
    fOutput->Add(fhJetCounter);
    fhJetArea= new TH1F("fhJetArea", "Jet Area", 400,-0.5, 1.95);
    fOutput->Add(fhJetArea);
    fhTrackPt= new TH1F("fhTrackPt", "Track Pt",600,-0.5,59.5);   
    fOutput->Add(fhTrackPt);
    fhGroomedPtvJetPt= new TH2F("fhGroomedPtvJetPt","Groomed Jet p_{T} v Original Jet p_{T}",150,0,150,150,0,150);
    fOutput->Add(fhGroomedPtvJetPt);

    if(fJetSelection == kRecoil){
      fhPtTriggerHadron= new TH1F("fhPtTriggerHadron", "fhPtTriggerHadron",1500,-0.5,149.5);  
      fOutput->Add(fhPtTriggerHadron);
      fh2PtTriggerHadronJet= new TH2F("fh2PtTriggerHadronJet", "fh2PtTriggerHadronJet",1500,-0.5,149.5,1500,-0.5,149.5);  
      fOutput->Add(fh2PtTriggerHadronJet);
      fhPhiTriggerHadronJet= new TH1F("fhPhiTriggerHadronJet", "fhPhiTriggerHadronJet",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
      fOutput->Add(fhPhiTriggerHadronJet);
      fhPhiTriggerHadronEventPlane= new TH1F("fhPhiTriggerHadronEventPlane", "fhPhiTriggerHadronEventPlane",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
      fOutput->Add(fhPhiTriggerHadronEventPlane);
      fhPhiTriggerHadronEventPlaneTPC= new TH1F("fhPhiTriggerHadronEventPlaneTPC", "fhPhiTriggerHadronEventPlaneTPC",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));  
      fOutput->Add(fhPhiTriggerHadronEventPlaneTPC);
      
  
    }
  }
  PostData(1,fOutput);
  PostData(2,fTreeJetInfo);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecoilJetYield::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecoilJetYield::FillHistograms()
{
  AliMultSelection *MultSelection = 0x0;
  AliAODEvent *fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Error("UserExec","AOD not available");
    return 0;
  }
  MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  if(!MultSelection) {
    AliWarning("AliMultSelection object not found!");
  }
  else{
    Percentile_1 = MultSelection->GetMultiplicityPercentile("V0M");
  }
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin)) return 0;
  }

  ////////////////////Find Trigger Hadron////////////////////
  AliAODTrack *TriggerHadron = 0x0;
  if (fJetSelection == kRecoil) {
    //you have to set a flag and the limits of the pT interval for your trigger
    Int_t TriggerHadronLabel = SelectTriggerHadron(fPtMinTriggerHadron, fPtMaxTriggerHadron);
    if (TriggerHadronLabel==-99999) return 0;  //Trigger Hadron Not Found
    AliTrackContainer *PartCont =NULL;
    if (fJetShapeSub==kConstSub) PartCont = GetTrackContainer(1);
    else PartCont = GetTrackContainer(0);
    TClonesArray *TrackArray = PartCont->GetArray();
    TriggerHadron = static_cast<AliAODTrack*>(TrackArray->At(TriggerHadronLabel));
    if (!TriggerHadron) return 0;//No trigger hadron with label found   
    if(fSemigoodCorrect){
      Double_t HoleDistance=RelativePhi(TriggerHadron->Phi(),fHolePos);
      if(TMath::Abs(HoleDistance)+fHoleWidth+fJetRadius>TMath::Pi()-fRecoilAngularWindow) return 0;
    }
    fhPtTriggerHadron->Fill(TriggerHadron->Pt()); //Needed for per trigger Normalisation
    fhPhiTriggerHadronEventPlane->Fill(TMath::Abs(RelativePhiEventPlane(fEPV0,TriggerHadron->Phi()))); //fEPV0 is the event plane from AliAnalysisTaskEmcal
    fhPhiTriggerHadronEventPlaneTPC->Fill(TMath::Abs(RelativePhiEventPlane(((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane(),TriggerHadron->Phi()))); //TPC event plane 
  }


  ////////////////////kData////////////////////
  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kData){
    AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                         
    AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event 
    Int_t JetCounter=0; //Counts number of jets in event  
    Double_t JetPhi=0;
    Bool_t EventCounter=kFALSE;
    Double_t JetPt_ForThreshold=0;
    AliEmcalJet *Jet2= NULL;
    if(JetCont) {
      JetCont->ResetCurrentID();
      while((Jet1=JetCont->GetNextAcceptJet())) {
	if(!Jet1) continue;
	fJetInfoVar[16]=Jet1->Pt();
	fJetInfoVar[17]=0;
	if (fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) JetPt_ForThreshold = Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	else JetPt_ForThreshold = Jet1->Pt();
	if(JetPt_ForThreshold<fPtThreshold) {
	  continue;
	}
	else {
	  Float_t RecoilDeltaPhi = 0.;
	  if (fJetSelection == kRecoil){
	    RecoilDeltaPhi = RelativePhi(TriggerHadron->Phi(), Jet1->Phi()); 
	    if (TMath::Abs(RecoilDeltaPhi) < (TMath::Pi() - fRecoilAngularWindow)) continue;  //accept the jet only if it overlaps with the recoil phi area of the trigger
	    fh2PtTriggerHadronJet->Fill(TriggerHadron->Pt(), Jet1->Pt());
	    fhPhiTriggerHadronJet->Fill(RelativePhi(TriggerHadron->Phi(), Jet1->Phi()));
	  }
          if (!EventCounter){
            EventCounter=kTRUE;
          }
	  JetCounter++;
	  fhJetPt->Fill(Jet1->Pt());
	  fhJetArea->Fill(Jet1->Area());
	  JetPhi=Jet1->Phi();
	  if(JetPhi < -1*TMath::Pi()) JetPhi += (2*TMath::Pi());
	  else if (JetPhi > TMath::Pi()) JetPhi -= (2*TMath::Pi());
	  fhJetPhi->Fill(JetPhi);
	  fhJetEta->Fill(Jet1->Eta());
	  fhJetMass->Fill(Jet1->M());
	  fhJetRadius->Fill(TMath::Sqrt((Jet1->Area()/TMath::Pi()))); //Radius of Jets per event
          fhNumberOfJetTracks->Fill(Jet1->GetNumberOfTracks());
	  if(fJetShapeSub==kNoSub || fJetShapeSub==kDerivSub) fJetInfoVar[0]= Jet1->Pt()-(GetRhoVal(0)*Jet1->Area());
	  else fJetInfoVar[0]=Jet1->Pt();
	  fJetInfoVar[1]=0;
	  fJetInfoVar[2]=JetPhi;
	  fJetInfoVar[3]=0;
	  fJetInfoVar[4]=Jet1->Eta();
	  if(fDoSoftDrop) {
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kFALSE);
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kTRUE);
	  }
	  else{
	    fJetInfoVar[6]=0;
	    fJetInfoVar[7]=0;
	    fJetInfoVar[8]=0;
	    fJetInfoVar[9]=0;
	    fJetInfoVar[10]=0;
	    fJetInfoVar[11]=0;
	    fJetInfoVar[24]=0;
	    fJetInfoVar[25]=0;
	  }		    
	  fJetInfoVar[12]=Jet1->M();
	  fJetInfoVar[13]=0;
	  fJetInfoVar[14]=Jet1->GetNumberOfConstituents();
	  fJetInfoVar[15]=0;
	  fJetInfoVar[18]=PTD(Jet1,0);
	  fJetInfoVar[19]=0;
	  fJetInfoVar[20]=Angularity(Jet1,0);
	  fJetInfoVar[21]=0;
	  fJetInfoVar[22]=fCent;
	  fJetInfoVar[23]=0;
	  fTreeJetInfo->Fill();
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
    }
  }

  ////////////////////kEmbedded////////////////////
  if(fJetShapeType == kDetEmbPart){
    AliEmcalJet *JetHybridS = NULL; //Subtracted hybrid Jet  
    AliEmcalJet *JetHybridUS = NULL; //Unsubtracted Hybrid Jet                                                                                                                     
    AliEmcalJet *JetPythDet = NULL; //Detector Level Pythia Jet
    AliEmcalJet *JetPythTrue = NULL; //Particle Level Pyhtia Jet
    AliJetContainer *JetContHybridS= GetJetContainer(0); //Jet Container for Subtracted Hybrid Jets 
    AliJetContainer *JetContHybridUS= GetJetContainer(1); //Jet Container for Unsubtracted Hybrid Jets                                                                                  
    AliJetContainer *JetContPythDet= GetJetContainer(2); //Jet Container for Detector Level Pyhtia Jets 
    AliJetContainer *JetContPythTrue= GetJetContainer(3); //Jet Container for Particle Level Pythia Jets

    Bool_t JetsMatched = kFALSE;
    Double_t JetPtThreshold;

    if(fJetShapeSub==kConstSub){
      while((JetHybridS = JetContHybridS->GetNextAcceptJet())){
	fJetInfoVar[16]=JetHybridS->Pt();
	fJetInfoVar[17]=0;
	if (fJetShapeSub==kConstSub) JetPtThreshold=JetHybridS->Pt();
	else JetPtThreshold=JetHybridS->Pt()-(GetRhoVal(0)*JetHybridS->Area());
        if ( (!JetHybridS) || (JetPtThreshold<fPtThreshold)) continue;
	Int_t JetNumber=-1;
	for(Int_t i = 0; i<JetContHybridUS->GetNJets(); i++) {
	  JetHybridUS = JetContHybridUS->GetJet(i);
	  if(JetHybridUS->GetLabel()==JetHybridS->GetLabel()) {
	    JetNumber=i;
	  }
	}
	if(JetNumber==-1) continue;
	JetHybridUS=JetContHybridUS->GetJet(JetNumber);
	if (JetContHybridUS->AliJetContainer::GetFractionSharedPt(JetHybridUS)<fSharedFractionPtMin) continue;
	JetPythDet=JetHybridUS->ClosestJet();

	if(!(fJetShapeSub==kConstSub)) JetHybridUS = JetHybridS->ClosestJet();
        if (!JetHybridUS) {
          Printf("Unsubtracted embedded jet does not exist, returning");
          continue;
        }
	if (!JetPythDet) continue;
	JetPythTrue=JetPythDet->ClosestJet();
	if(!JetPythTrue) continue;
	JetsMatched=kTRUE;

	fJetInfoVar[0]=JetHybridS->Pt();
	fJetInfoVar[1]=JetPythTrue->Pt();
	fJetInfoVar[2]=JetHybridS->Phi();
	fJetInfoVar[3]=JetPythTrue->Phi();
	fJetInfoVar[4]=JetHybridS->Eta();
	fJetInfoVar[5]=JetPythTrue->Eta();
	if(fDoSoftDrop) {
	  SoftDrop(JetHybridS,JetContHybridS,fZCut,fBeta_SD,kFALSE);
	  SoftDrop(JetPythTrue,JetContPythTrue,fZCut,fBeta_SD,kTRUE);
	}
	else{
	  fJetInfoVar[6]=0;
	  fJetInfoVar[7]=0;
	  fJetInfoVar[8]=0;
	  fJetInfoVar[9]=0;
	  fJetInfoVar[10]=0;
	  fJetInfoVar[11]=0;
	  fJetInfoVar[24]=0;
	  fJetInfoVar[25]=0;
	}		    
	fJetInfoVar[12]=JetHybridS->M();
	fJetInfoVar[13]=JetPythTrue->M();
	fJetInfoVar[14]=JetHybridS->GetNumberOfConstituents();
	fJetInfoVar[15]=JetPythTrue->GetNumberOfConstituents();
	fJetInfoVar[18]=PTD(JetHybridS,0);
	fJetInfoVar[19]=PTD(JetPythTrue,0);
	fJetInfoVar[20]=Angularity(JetHybridS,0);
	fJetInfoVar[21]=Angularity(JetPythTrue,0);
	fJetInfoVar[22]=fCent;
	fJetInfoVar[23]=fCent;
	fTreeJetInfo->Fill();
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRecoilJetYield::RelativePhiEventPlane(Double_t EventPlane, Double_t Phi){

  if(Phi < -1*TMath::Pi()) Phi += (2*TMath::Pi());
  else if (Phi > TMath::Pi()) Phi -= (2*TMath::Pi());
  Double_t DeltaPhi=Phi-EventPlane;
  if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
  else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  return DeltaPhi;
}
//________________________________________________________________________
Double_t AliAnalysisTaskRecoilJetYield::RelativePhi(Double_t Phi1, Double_t Phi2){

  if(Phi1 < -1*TMath::Pi()) Phi1 += (2*TMath::Pi());                                                           
  else if (Phi1 > TMath::Pi()) Phi1 -= (2*TMath::Pi());
  if(Phi2 < -1*TMath::Pi()) Phi2 += (2*TMath::Pi());
  else if (Phi2 > TMath::Pi()) Phi2 -= (2*TMath::Pi());
  Double_t DeltaPhi=Phi2-Phi1;
  if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
  else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  return DeltaPhi;
}




//--------------------------------------------------------------------------
Double_t AliAnalysisTaskRecoilJetYield::Angularity(AliEmcalJet *Jet, Int_t JetContNb){
  
  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  Double_t Angularity_Numerator=0;                                                                          
  Double_t Angularity_Denominator=0;
  AliVParticle *Particle=0x0;
  Double_t DeltaPhi=0;


  for (Int_t i=0; i< Jet->GetNumberOfTracks(); i++){  //loops through all tracks                                                                              
    Particle = static_cast<AliVParticle*>(Jet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));                                                                 
    if(!Particle) continue;
    DeltaPhi=RelativePhi(Jet->Phi(),Particle->Phi());
    Angularity_Numerator=Angularity_Numerator+(Particle->Pt()*TMath::Sqrt(((Particle->Eta()-Jet->Eta())*(Particle->Eta()-Jet->Eta()))+(DeltaPhi*DeltaPhi)));
    Angularity_Denominator= Angularity_Denominator+Particle->Pt();
  }
  if(Angularity_Denominator!=0) return Angularity_Numerator/Angularity_Denominator;
  else return -1;

}



//--------------------------------------------------------------------------                                                                                                         
Double_t AliAnalysisTaskRecoilJetYield::PTD(AliEmcalJet *Jet, Int_t JetContNb){

  AliJetContainer *JetCont = GetJetContainer(JetContNb);
  Double_t PTD_Numerator=0;  //Reset these values                                                                                                                            
  Double_t PTD_Denominator=0;
  AliVParticle *Particle=0x0;
  Double_t DeltaPhi=0;
  for (Int_t i=0; i< Jet->GetNumberOfTracks(); i++){  //loops through all tracks                                                                             
    Particle = static_cast<AliVParticle*>(Jet->TrackAt(i, JetCont->GetParticleContainer()->GetArray()));
    if(!Particle) continue;
    PTD_Numerator=PTD_Numerator+(Particle->Pt()*Particle->Pt());
    PTD_Denominator=PTD_Denominator+Particle->Pt();
  }
  if(PTD_Denominator!=0) return TMath::Sqrt(PTD_Numerator)/PTD_Denominator;
  else return -1;

}

//_________________________________________________________________________
 void AliAnalysisTaskRecoilJetYield::SoftDrop(AliEmcalJet *fJet,AliJetContainer *fJetCont, double zcut, double beta, Bool_t fTruthJet){
  std::vector<fastjet::PseudoJet>        fInputVectors;
  Double_t JetInvMass=0, PseudJetInvMass=0, TrackMom = 0, TrackEnergy = 0;
  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  Double_t JetEta=fJet->Eta(),JetPhi=fJet->Phi();
  Double_t FJTrackEta[9999],FJTrackPhi[9999],FJTrackPt[9999],EmcalJetTrackEta[9999],EmcalJetTrackPhi[9999],EmcalJetTrackPt[9999];
  UShort_t FJNTracks=0,EmcalJetNTracks=0;
  if (fTrackCont) for (Int_t i=0; i<fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      JetInvMass += fTrk->M();
      if (!fTrk) continue;
      fastjet::PseudoJet PseudoTracks(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      TrackMom += TMath::Sqrt(TMath::Power(fTrk->Px(),2)+TMath::Power(fTrk->Py(),2)+TMath::Power(fTrk->Pz(),2));
      TrackEnergy += fTrk->E();
      PseudoTracks.set_user_index(fJet->TrackAt(i)+100);
      PseudJetInvMass += PseudoTracks.m();
      fInputVectors.push_back(PseudoTracks);
      EmcalJetTrackEta[i]=fTrk->Eta();
      EmcalJetTrackPhi[i]=fTrk->Phi();
      EmcalJetTrackPt[i]=fTrk->Pt();
      EmcalJetNTracks++;
    }
  fastjet::JetDefinition                *fJetDef;         
  fastjet::ClusterSequence              *fClustSeqSA;
  


  fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fJetRadius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 

  try {
    fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    //return -1;
  }

  std::vector<fastjet::PseudoJet>       fOutputJets;
  fOutputJets.clear();
  fOutputJets=fClustSeqSA->inclusive_jets(0);
 

  std::vector<fastjet::PseudoJet> jet_constituents = fOutputJets[0].constituents();
  Float_t NSubjettinessResult1, NSubjettinessResult2, NSubBeta = 1, R0=0.4;
  std::vector<fastjet::PseudoJet> Subjet_Axes;
  fastjet::PseudoJet SubJet1_Axis,SubJet2_Axis;
  if(jet_constituents.size()>=2){
    fastjet::contrib::Nsubjettiness nSub1(1,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedMeasure(NSubBeta,R0));
    fastjet::contrib::Nsubjettiness nSub2(2,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedMeasure(NSubBeta,R0));

    NSubjettinessResult1 = nSub1.result(fOutputJets[0]);
    NSubjettinessResult2 = nSub2.result(fOutputJets[0]);
    Subjet_Axes = nSub2.currentAxes();
    SubJet1_Axis = Subjet_Axes[0];
    SubJet2_Axis = Subjet_Axes[1];

    Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
    Double_t SubJet2_Eta=SubJet2_Axis.pseudorapidity();
    Double_t SubJet1_Phi=SubJet1_Axis.phi();
    if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
    else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
    Double_t SubJet2_Phi=SubJet2_Axis.phi();
    if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
    else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());

    Double_t DeltaPhiSubJets,DeltaEtaSubJets;
    DeltaPhiSubJets=SubJet1_Phi-SubJet2_Phi;
    DeltaEtaSubJets=SubJet1_Eta-SubJet2_Eta;
    if(DeltaPhiSubJets < -1*TMath::Pi()) DeltaPhiSubJets += (2*TMath::Pi());
    else if (DeltaPhiSubJets > TMath::Pi()) DeltaPhiSubJets -= (2*TMath::Pi());
    Double_t DelR=0;
    DelR = TMath::Power(TMath::Power(DeltaPhiSubJets,2)+TMath::Power(DeltaEtaSubJets,2),0.5);

    if(!fTruthJet) fJetInfoVar[24]=DelR;
    else fJetInfoVar[25]=DelR;
  }
  
  fastjet::contrib::SoftDrop softdrop(beta, zcut);
  fastjet::PseudoJet finaljet = softdrop(fOutputJets[0]);
  if(!fTruthJet){
    fJetInfoVar[8]=NSubjettinessResult1;
    fJetInfoVar[10]=NSubjettinessResult2;
  }
  else {
    fJetInfoVar[9]=NSubjettinessResult1;
    fJetInfoVar[11]=NSubjettinessResult2;
  }
  AliEmcalJet* jet = new AliEmcalJet(finaljet.perp(), finaljet.eta(), finaljet.phi(), finaljet.m());
  std::vector<fastjet::PseudoJet> fSDTracks=finaljet.constituents();
  Double_t FastjetTrackDelR,EmcalTrackDelR;
  for(Int_t i=0;i<fJet->GetNumberOfConstituents();i++){
    if(i<=finaljet.constituents().size()){
      FastjetTrackDelR = TMath::Sqrt(TMath::Power(fSDTracks[i].eta()-JetEta,2)+TMath::Power(fSDTracks[i].phi()-JetPhi,2));
      FJTrackEta[i]=fSDTracks[i].eta();
      FJTrackPhi[i]=fSDTracks[i].phi();
      FJTrackPt[i]=fSDTracks[i].perp();
      FJNTracks++;
    }
    AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
    EmcalTrackDelR = TMath::Sqrt(TMath::Power(fTrk->Eta()-JetEta,2)+TMath::Power(fTrk->Phi()-JetPhi,2));       
  }
  Int_t NDroppedTracks = fJet->GetNumberOfTracks()-finaljet.constituents().size();
  Int_t nConstituents(fClustSeqSA->constituents(finaljet).size());
  jet->SetNumberOfTracks(nConstituents);
  Double_t SymParam, Mu, DeltaR;
  SymParam=(finaljet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
  Mu=(finaljet.structure_of<fastjet::contrib::SoftDrop>().mu());
  DeltaR=(finaljet.structure_of<fastjet::contrib::SoftDrop>().delta_R());
  fhGroomedPtvJetPt->Fill(finaljet.perp(),fJet->Pt());
  if(!fTruthJet) fJetInfoVar[6]=SymParam;
  else fJetInfoVar[7]=SymParam;
  return;

  
}

//________________________________________________________________________
Int_t AliAnalysisTaskRecoilJetYield::SelectTriggerHadron(Float_t PtMin, Float_t PtMax){
  AliTrackContainer *PartCont = NULL;
  if (fJetShapeSub==kConstSub) PartCont = GetTrackContainer(1);
  else PartCont = GetTrackContainer(0);
  TClonesArray *TracksArray = PartCont->GetArray();
  if(!PartCont || !TracksArray) return -99999;
  AliAODTrack *Track = 0x0;
  Int_t Trigger_Index[10000];
  for (Int_t i=0; i<100; i++) Trigger_Index[i] = 0;
  Int_t Trigger_Counter = 0;
  for(Int_t i=0; i < TracksArray->GetEntriesFast(); i++){  
    if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
      if (!Track) continue;
      fhTrackPt->Fill(Track->Pt());
      if(TMath::Abs(Track->Eta())>0.9) continue;
      if (Track->Pt()<0.15) continue;
      if ((Track->Pt() >= PtMin) && (Track->Pt()< PtMax)) {
	Trigger_Index[Trigger_Counter] = i;
	Trigger_Counter++;
      }
    }
  }
  if (Trigger_Counter == 0) return -99999;
  Int_t RandomNumber = 0, Index = 0 ; 
  TRandom3* Random = new TRandom3(0); 
  RandomNumber = Random->Integer(Trigger_Counter);
  Index = Trigger_Index[RandomNumber];
  return Index; 
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecoilJetYield::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskRecoilJetYield::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
  
  
}
