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
  fhDroppedBranches(0x0),
  fhPtTriggerHadron(0x0), 
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0),
  fhDetJetPt_Incl(0x0),
  fhDetJetPt_Matched(0x0),
  fReclusterAlgo(0),
  fSubMatching(kFALSE)

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
  fhDroppedBranches(0x0),
  fhPtTriggerHadron(0x0), 
  fh2PtTriggerHadronJet(0x0),
  fhPhiTriggerHadronJet(0x0),
  fhPhiTriggerHadronEventPlane(0x0),
  fhPhiTriggerHadronEventPlaneTPC(0x0),
  fhDetJetPt_Incl(0x0),
  fhDetJetPt_Matched(0x0),
  fReclusterAlgo(0),
  fSubMatching(kFALSE)
  
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
    const Int_t nVarMin = 18;
    TString *fJetInfoVarNames = new TString [nVarMin];
  
    fJetInfoVarNames[0] = "Pt";
    fJetInfoVarNames[1] = "Pt_Truth";
    fJetInfoVarNames[2] = "SymParam";
    fJetInfoVarNames[3] = "SymParam_Truth";
    fJetInfoVarNames[4] = "Tau1";
    fJetInfoVarNames[5] = "Tau1_Truth";
    fJetInfoVarNames[6] = "Tau2";
    fJetInfoVarNames[7] = "Tau2_Truth";
    fJetInfoVarNames[8] = "PTD";
    fJetInfoVarNames[9] = "PTD_Truth";
    fJetInfoVarNames[10] = "Angularity";
    fJetInfoVarNames[11] = "Angularity_Truth";
    fJetInfoVarNames[12] = "DelR";
    fJetInfoVarNames[13] = "DelR_Truth";
    fJetInfoVarNames[14] = "N_Groomed_Branches";
    fJetInfoVarNames[15] = "N_Groomed_Branches_Truth";
    fJetInfoVarNames[16] = "Groomed_Jet_Pt";
    fJetInfoVarNames[17] = "Groomed_Jet_Pt_Truth";

  
    
    for(Int_t ivar=0; ivar < nVarMin; ivar++){
      cout<<"looping over variables"<<endl;
      fTreeJetInfo->Branch(fJetInfoVarNames[ivar].Data(), &fJetInfoVar[ivar], Form("%s/D", fJetInfoVarNames[ivar].Data()));
    }
  }

  if (fJetShapeType==AliAnalysisTaskRecoilJetYield::kData || fJetShapeType==AliAnalysisTaskRecoilJetYield::kDetEmbPart){
    
    fhJetPt= new TH1F("fhJetPt", "Jet Pt",150,-0.5,149.5 );   
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
    fhDroppedBranches= new TH1F("fhDroppedBranches","Number of Softdropped branches",50,0,50);
    fOutput->Add(fhDroppedBranches);
    fhDetJetPt_Incl= new TH1F("fhDetJetPt_Incl", "Jet Pt",200,-0.5,199.5 );   
    fOutput->Add(fhDetJetPt_Incl);
    fhDetJetPt_Matched= new TH1F("fhDetJetPt_Matched", "Jet Pt",200,-0.5,199.5 );   
    fOutput->Add(fhDetJetPt_Matched);

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
  //cout<<"End of UserCreateOutputObjects"<<endl;
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
  //AliMultSelection *MultSelection = 0x0;
  AliAODEvent *fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Error("UserExec","AOD not available");
    return 0;
  }
  /*MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  if(!MultSelection) {
    AliWarning("AliMultSelection object not found!");
  }
  else{
    Percentile_1 = MultSelection->GetMultiplicityPercentile("V0M");
  }
  */
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
    AliParticleContainer *PartContMC = NULL;
    if (fJetShapeSub==kConstSub){
      if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) PartContMC = GetParticleContainer(1);
      else PartCont = GetTrackContainer(1);
    }
    else {
      if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) PartContMC = GetParticleContainer(0);
      else PartCont = GetTrackContainer(0);
    }
    TClonesArray *TrackArray = NULL;
    TClonesArray *TrackArrayMC = NULL;
    if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) TrackArrayMC = PartContMC->GetArray();
    else TrackArray = PartCont->GetArray();
    if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) TriggerHadron = static_cast<AliAODTrack*>(TrackArrayMC->At(TriggerHadronLabel));
    else TriggerHadron = static_cast<AliAODTrack*>(TrackArray->At(TriggerHadronLabel));
    if (!TriggerHadron) return 0;//No trigger hadron with label found   
    if(fSemigoodCorrect){
      Double_t HoleDistance=RelativePhi(TriggerHadron->Phi(),fHolePos);
      if(TMath::Abs(HoleDistance)+fHoleWidth+fJetRadius>TMath::Pi()-fRecoilAngularWindow) return 0;
    }
    fhPtTriggerHadron->Fill(TriggerHadron->Pt()); //Needed for per trigger Normalisation
    if (fJetShapeType != AliAnalysisTaskRecoilJetYield::kGenOnTheFly) fhPhiTriggerHadronEventPlane->Fill(TMath::Abs(RelativePhiEventPlane(fEPV0,TriggerHadron->Phi()))); //fEPV0 is the event plane from AliAnalysisTaskEmcal
    if (fJetShapeType != AliAnalysisTaskRecoilJetYield::kGenOnTheFly) fhPhiTriggerHadronEventPlaneTPC->Fill(TMath::Abs(RelativePhiEventPlane(((AliVAODHeader*)(InputEvent()->GetHeader()))->GetEventplane(),TriggerHadron->Phi()))); //TPC event plane 
  }

  ////////////////////kData////////////////////
  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kData){
    // cout<<"entering kData"<<endl;

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
	  if(fDoSoftDrop) {
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kFALSE);
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kTRUE);
	  }
	  else{
	    fJetInfoVar[2]=0;
	    fJetInfoVar[3]=0;
	    fJetInfoVar[4]=0;
	    fJetInfoVar[5]=0;
	    fJetInfoVar[6]=0;
	    fJetInfoVar[7]=0;
	    fJetInfoVar[12]=0;
	    fJetInfoVar[13]=0;
	    fJetInfoVar[14]=0;
	    fJetInfoVar[15]=0;
	    fJetInfoVar[16]=0;
	    fJetInfoVar[17]=0;
	  }		    
	  fJetInfoVar[8]=PTD(Jet1,0);
	  fJetInfoVar[9]=0;
	  fJetInfoVar[10]=Angularity(Jet1,0);
	  fJetInfoVar[11]=0;
	  fTreeJetInfo->Fill();
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
    }
  }

  ////////////////////kEmbedded////////////////////
  if(fJetShapeType == kDetEmbPart){
    AliEmcalJet *JetHybridS = NULL; //Subtracted hybrid Jet  
    AliEmcalJet *JetHybridUS = NULL; //Unsubtracted Hybrid Jet     //For matching SubtractedHybrid->DetPythia this jet container is also Subtracted Hybrid                                                                                                                
    AliEmcalJet *JetPythDet = NULL; //Detector Level Pythia Jet
    AliEmcalJet *JetPythTrue = NULL; //Particle Level Pyhtia Jet
    AliJetContainer *JetContHybridS= GetJetContainer(0); //Jet Container for Subtracted Hybrid Jets 
    AliJetContainer *JetContHybridUS= GetJetContainer(1); //Jet Container for Unsubtracted Hybrid Jets                                                                                  
    AliJetContainer *JetContPythDet= GetJetContainer(2); //Jet Container for Detector Level Pyhtia Jets 
    AliJetContainer *JetContPythTrue= GetJetContainer(3); //Jet Container for Particle Level Pythia Jets

  

    Bool_t JetsMatched = kFALSE;
    Double_t JetPtThreshold;
    JetContHybridS->ResetCurrentID();
    JetContHybridUS->ResetCurrentID();
    JetContPythDet->ResetCurrentID();
    JetContPythTrue->ResetCurrentID();

    while((JetPythDet = JetContPythDet->GetNextAcceptJet())){
      fhDetJetPt_Incl->Fill(JetPythDet->Pt()); //Fill histogram with all Detector level jets
    }
      
      
    while((JetHybridS = JetContHybridS->GetNextAcceptJet())){
      if (fJetShapeSub==kConstSub) JetPtThreshold=JetHybridS->Pt();
      else JetPtThreshold=JetHybridS->Pt()-(GetRhoVal(0)*JetHybridS->Area());
      if ( (!JetHybridS) || (JetPtThreshold<fPtThreshold)) continue;
      if (fJetShapeSub==kConstSub){
	Int_t JetNumber=-1;
	for(Int_t i = 0; i<JetContHybridUS->GetNJets(); i++) {
	  JetHybridUS = JetContHybridUS->GetJet(i);
	  if (!JetHybridUS) continue;
	    
	  if(JetHybridUS->GetLabel()==JetHybridS->GetLabel()) {
	    JetNumber=i;
	  }
	}
	if(JetNumber==-1) continue;
	JetHybridUS=JetContHybridUS->GetJet(JetNumber);
	//if(JetHybridUS) cout<<"Matched to jet i = "<< JetNumber<<endl;
	if (JetContHybridUS->AliJetContainer::GetFractionSharedPt(JetHybridUS)<fSharedFractionPtMin && !fSubMatching) {
	  //cout<<"Fraction shared pt below cut = "<<JetContHybridUS->AliJetContainer::GetFractionSharedPt(JetHybridUS)<<endl;
	  continue;
	}
	else if(GetFractionSharedPt_SubMatching(JetHybridUS,JetContHybridUS->GetParticleContainer())<fSharedFractionPtMin && fSubMatching){
	  continue;
	}
	JetPythDet=JetHybridUS->ClosestJet();
      }
      if(!(fJetShapeSub==kConstSub)){
	if (JetContHybridS->AliJetContainer::GetFractionSharedPt(JetHybridS)<fSharedFractionPtMin) continue;
	    JetPythDet = JetHybridS->ClosestJet();
      }
      if (!JetHybridUS) {
	Printf("Unsubtracted embedded jet does not exist, returning");
	continue;
      }
      if (!JetPythDet) continue;
      fhDetJetPt_Matched->Fill(JetPythDet->Pt()); //Fill only matched detector level jets for tagging efficiency comparison
      JetPythTrue=JetPythDet->ClosestJet();
      if(!JetPythTrue) continue;
      JetsMatched=kTRUE;

      if(fJetShapeSub==kConstSub) fJetInfoVar[0]=JetHybridS->Pt();
      else fJetInfoVar[0]=JetHybridS->Pt()-(GetRhoVal(0)*JetHybridS->Area());
      fJetInfoVar[1]=JetPythTrue->Pt();
      if(fDoSoftDrop) {
	SoftDrop(JetHybridS,JetContHybridS,fZCut,fBeta_SD,kFALSE);
	SoftDrop(JetPythTrue,JetContPythTrue,fZCut,fBeta_SD,kTRUE);
      }
      else{
	fJetInfoVar[2]=0;
	fJetInfoVar[3]=0;
	fJetInfoVar[4]=0;
	fJetInfoVar[5]=0;
	fJetInfoVar[6]=0;
	fJetInfoVar[7]=0;
	fJetInfoVar[12]=0;
	fJetInfoVar[13]=0;
	fJetInfoVar[14]=0;
	fJetInfoVar[15]=0;
	fJetInfoVar[16]=0;
	fJetInfoVar[17]=0;
      }		    
      fJetInfoVar[8]=PTD(JetHybridS,0);
      fJetInfoVar[9]=PTD(JetPythTrue,0);
      fJetInfoVar[10]=Angularity(JetHybridS,0);
      fJetInfoVar[11]=Angularity(JetPythTrue,0);
      fTreeJetInfo->Fill();
    }
    
  }

  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly){
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
	  if(fDoSoftDrop) {
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kFALSE);
	    SoftDrop(Jet1,JetCont,fZCut,fBeta_SD,kTRUE);
	  }
	  else{
	    fJetInfoVar[2]=0;
	    fJetInfoVar[3]=0;
	    fJetInfoVar[4]=0;
	    fJetInfoVar[5]=0;
	    fJetInfoVar[6]=0;
	    fJetInfoVar[7]=0;
	    fJetInfoVar[12]=0;
	    fJetInfoVar[13]=0;
	    fJetInfoVar[14]=0;
	    fJetInfoVar[15]=0;
	    fJetInfoVar[16]=0;
	    fJetInfoVar[17]=0;
	  }		    
	  fJetInfoVar[8]=PTD(Jet1,0);
	  fJetInfoVar[9]=0;
	  fJetInfoVar[10]=Angularity(Jet1,0);
	  fJetInfoVar[11]=0;
	  fTreeJetInfo->Fill();
	}
      }
      fhJetCounter->Fill(JetCounter); //Number of Jets in Each Event
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
  Float_t NSubjettinessResult[3], NSubBeta = 1, R0=0.4;
  std::vector<fastjet::PseudoJet> Subjet_Axes;
  fastjet::PseudoJet SubJet1_Axis,SubJet2_Axis;
  Double_t DelR=-5;
	
  
  for(Int_t j=1; j<3; j++){
    if(jet_constituents.size() < j){
      if(!fTruthJet){
	if(j==1) fJetInfoVar[4]=-5;
	else if(j==2) fJetInfoVar[6]=-5;
      }
      else {
	if(j==1) fJetInfoVar[5]=-5;
	else if(j==2) fJetInfoVar[7]=-5;
      }
      continue;
    }
    else{
      fastjet::contrib::Nsubjettiness nSub(j,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedMeasure(NSubBeta,R0));
      NSubjettinessResult[j] = nSub.result(fOutputJets[0]);
      if(j==2){ //find subjet axis to calculate delR
	Subjet_Axes = nSub.currentAxes();
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
	
	DelR = TMath::Power(TMath::Power(DeltaPhiSubJets,2)+TMath::Power(DeltaEtaSubJets,2),0.5);
      }
    }
  }
    if(!fTruthJet){
      fJetInfoVar[4]=NSubjettinessResult[1];
      fJetInfoVar[6]=NSubjettinessResult[2];
      // fJetInfoVar[12]=DelR;
    }
    else {
      fJetInfoVar[5]=NSubjettinessResult[1];
      fJetInfoVar[7]=NSubjettinessResult[2];
      // fJetInfoVar[13]=DelR;
    }

    
  
  fastjet::contrib::SoftDrop softdrop(beta, zcut);
  //fastjet::contrib::SoftDrop softdrop_antikt(beta,zcut);
  softdrop.set_verbose_structure(kTRUE);
  //fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, 0.4);
  // fastjet::contrib::Recluster *antiKT_Recluster(jet_def_akt);
  fastjet::contrib::Recluster *recluster;
  if(fReclusterAlgo == 2) recluster = new fastjet::contrib::Recluster(fastjet::kt_algorithm,1,true);
  if(fReclusterAlgo == 1) recluster = new fastjet::contrib::Recluster(fastjet::antikt_algorithm,1,true);
  if(fReclusterAlgo == 0) recluster = new fastjet::contrib::Recluster(fastjet::cambridge_algorithm,1,true);  
  softdrop.set_reclustering(true,recluster);
  fastjet::PseudoJet finaljet = softdrop(fOutputJets[0]);
  // fastjet::PseudoJet finaljet_antikt = softdrop_antikt(fOutputJets[0]);
  //cout<< finaljet.structure_of<fastjet::contrib::SoftDrop>().symmetry()<<endl;
  //cout<< finaljet_antikt.structure_of<fastjet::contrib::SoftDrop>().symmetry()<<endl;


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
  fhDroppedBranches->Fill(finaljet.structure_of<fastjet::contrib::SoftDrop>().dropped_count());
  if(!fTruthJet) fJetInfoVar[2]=SymParam;
  else fJetInfoVar[3]=SymParam;
  if(!fTruthJet) fJetInfoVar[12] = DeltaR;
  else fJetInfoVar[13] = DeltaR;
  if(!fTruthJet) fJetInfoVar[14]=finaljet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
  else fJetInfoVar[15]=finaljet.structure_of<fastjet::contrib::SoftDrop>().dropped_count();
  if(!fTruthJet) fJetInfoVar[16]=finaljet.perp();
  else fJetInfoVar[17]=finaljet.perp();
  
  
  return;

  
}

//--------------------------------------------------------------------------
Int_t AliAnalysisTaskRecoilJetYield::SelectTriggerHadron(Float_t PtMin, Float_t PtMax){

  AliTrackContainer *PartCont = NULL;
  AliParticleContainer *PartContMC = NULL;


  if (fJetShapeSub==kConstSub){
    if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) PartContMC = GetParticleContainer(1);
    else PartCont = GetTrackContainer(1);
  }
  else{
    if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) PartContMC = GetParticleContainer(0);
    else PartCont = GetTrackContainer(0);
  }
  
  TClonesArray *TracksArray = NULL;
  TClonesArray *TracksArrayMC = NULL;
  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) TracksArrayMC = PartContMC->GetArray();
  else TracksArray = PartCont->GetArray();
 
  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly){
    if(!PartContMC || !TracksArrayMC) return -99999;
  }
  else {
    if(!PartCont || !TracksArray) return -99999;
  }
    
  AliAODTrack *Track = 0x0;
  Int_t Trigger_Index[100];
  for (Int_t i=0; i<100; i++) Trigger_Index[i] = 0;
  Int_t Trigger_Counter = 0;
  Int_t NTracks=0;
  if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly) NTracks = TracksArrayMC->GetEntriesFast();
  else NTracks = TracksArray->GetEntriesFast();
  for(Int_t i=0; i < NTracks; i++){
    if (fJetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly){
      if((Track = static_cast<AliAODTrack*>(PartContMC->GetAcceptParticle(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= PtMin) && (Track->Pt()< PtMax)) {
	  Trigger_Index[Trigger_Counter] = i;
	  Trigger_Counter++;
	}
      }
    }
    else{ 
      if((Track = static_cast<AliAODTrack*>(PartCont->GetAcceptTrack(i)))){
	if (!Track) continue;
	if(TMath::Abs(Track->Eta())>0.9) continue;
	if (Track->Pt()<0.15) continue;
	if ((Track->Pt() >= PtMin) && (Track->Pt()< PtMax)) {
	  Trigger_Index[Trigger_Counter] = i;
	  Trigger_Counter++;
	}
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
Double_t AliAnalysisTaskRecoilJetYield::GetFractionSharedPt_SubMatching(const AliEmcalJet *jet1, AliParticleContainer *cont2) const
{
  AliEmcalJet *jet2 = jet1->ClosestJet();
  if (!jet2) return -1;

  Double_t jetPt2 = jet2->Pt();
  if (jetPt2 <= 0) return -1;

  Int_t bgeom = kTRUE;
  if (!cont2) bgeom = kFALSE;
  Double_t sumPt = 0.;
  AliVParticle *vpf = 0x0;
  Int_t iFound = 0;
  for (Int_t icc = 0; icc < jet2->GetNumberOfTracks(); icc++) {
    Int_t idx = (Int_t)jet2->TrackAt(icc);
    //get particle
    AliVParticle *p2 = 0x0;
    if (bgeom) p2 = static_cast<AliVParticle*>(jet2->TrackAt(icc, cont2->GetArray()));
    iFound = 0;
    for (Int_t icf = 0; icf < jet1->GetNumberOfTracks(); icf++) {
      if (!bgeom && idx == jet1->TrackAt(icf) && iFound == 0 ) {
        iFound = 1;
        vpf = jet1->Track(icf);
        if (vpf) sumPt += vpf->Pt();
        continue;
      }
      if (bgeom){
        vpf = jet1->Track(icf);
        if (!vpf) continue;
        if (!SameParticle(vpf, p2, 1.e-4)) continue; //not the same particle
        sumPt += vpf->Pt();
      }
    }
  }

  Double_t fraction = sumPt / jetPt2;

  return fraction;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRecoilJetYield::SameParticle(const AliVParticle* part1, const AliVParticle* part2, Double_t dist)
{
  if(!part1) return kFALSE;
  if(!part2) return kFALSE;
  Double_t dPhi = TMath::Abs(part1->Phi() - part2->Phi());
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  if (dPhi > dist) return kFALSE;
  return kTRUE;
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
