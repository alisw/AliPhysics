/*_____________________________________________________________
 
 Class AliAnalysisHFCorrOnFlySim: On the fly Simulation class for
 heavy flavour correlations and general event/part properties 
 AUTHROS:
 Jitendra Kumar (jitendra.kumar@cern.ch)
 Andrea Rossi   (andrea.rossi@cern.ch)
 _____________________________________________________________*/

#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TChain.h"
#include <THnSparse.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include "TObjArray.h"
#include "TArrayI.h"
#include <TClonesArray.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliAnalysisHFCorrOnFlySim.h"
#include "AliGenEventHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
using namespace std;

ClassImp(AliAnalysisHFCorrOnFlySim)
//______________________________| Default Constructor
AliAnalysisHFCorrOnFlySim::AliAnalysisHFCorrOnFlySim():
  fMcEvent(0x0),
  fMcHandler(0x0),
  fHistEventsProcessed(0x0),
  fOutputQA(0),
  fOutputList(0),
  fEtaMin(-20),
  fEtaMax(20),
  fYMin(-20),
  fYMax(20),
  fPtMin(0),
  fPtMax(999),
  fMinMultiplicity(0),
  fMaxMultiplicity(2000),
  fCorrPart1(0x0),
  fCorrPart2(0x0),
  fChargeSel1(1),
  fChargeSel2(-1),
  fStack(0),
  fParticleArray(0),
  fIsEventProp(kTRUE),
  fIsPartProp(kTRUE),
  fIsCorrOfQQbar(kFALSE),
  fIsCorrOfHeavyFlavor(kFALSE),
  fIsCorrOfHadronHadron(kFALSE),
  flastdaugh(0)
{
  fArrayTrk        = new TArrayI(5000);
  fArraySkipDDaugh = new TArrayI(200);
}

//______________________________| Specific Constructor
AliAnalysisHFCorrOnFlySim::AliAnalysisHFCorrOnFlySim(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fMcEvent(0x0),
  fMcHandler(0x0),
  fHistEventsProcessed(0x0),
  fOutputQA(0),
  fOutputList(0),
  fEtaMin(-20),
  fEtaMax(20),
  fYMin(-20),
  fYMax(20),
  fPtMin(0),
  fPtMax(999.0),
  fMinMultiplicity(0),
  fMaxMultiplicity(2000),
  fCorrPart1(0x0),
  fCorrPart2(0x0),
  fChargeSel1(1),
  fChargeSel2(-1),
  fStack(0),
  fParticleArray(0X0),
  fIsEventProp(kTRUE),
  fIsPartProp(kTRUE),
  fIsCorrOfQQbar(kFALSE),
  fIsCorrOfHeavyFlavor(kFALSE),
  fIsCorrOfHadronHadron(kFALSE),
  flastdaugh(0)
{
  fArrayTrk        = new TArrayI(5000);
  fArraySkipDDaugh = new TArrayI(200);
  
  Info("AliAnalysisHFCorrOnFlySim","Calling Constructor");
  // Output slot #1 writes into a TList container (nevents histogram)
  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//______________________________| Copy Constructor
AliAnalysisHFCorrOnFlySim::AliAnalysisHFCorrOnFlySim(const AliAnalysisHFCorrOnFlySim& c) :
  AliAnalysisTaskSE(c),
  fMcEvent(c.fMcEvent),
  fMcHandler(c.fMcHandler),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fOutputQA(c.fOutputQA),
  fOutputList(c.fOutputList),
  fEtaMin(c.fEtaMin),
  fEtaMax(c.fEtaMax),
  fYMin(c.fYMin),
  fYMax(c.fYMax),
  fPtMin(c.fPtMin),
  fPtMax(c.fPtMax),
  fMinMultiplicity(c.fMinMultiplicity),
  fMaxMultiplicity(c.fMaxMultiplicity),
  fCorrPart1(c.fCorrPart1),
  fCorrPart2(c.fCorrPart2),
  fChargeSel1(c.fChargeSel1),
  fChargeSel2(c.fChargeSel2),
  fStack(c.fStack),
  fParticleArray(c.fParticleArray),
  fIsEventProp(c.fIsEventProp),
  fIsPartProp(c.fIsPartProp),
  fIsCorrOfQQbar(c.fIsCorrOfQQbar),
  fIsCorrOfHeavyFlavor(c.fIsCorrOfHeavyFlavor),
  fIsCorrOfHadronHadron(c.fIsCorrOfHadronHadron),
  fArraySkipDDaugh(c.fArraySkipDDaugh),
  fArrayTrk(c.fArrayTrk),
  flastdaugh(c.flastdaugh)
{
  // Copy Constructor
}

//______________________________| Destructor
AliAnalysisHFCorrOnFlySim::~AliAnalysisHFCorrOnFlySim()
{
  // Destructor
  Info("~AliAnalysisHFCorrOnFlySim","Calling Destructor");
  if (fHistEventsProcessed) delete fHistEventsProcessed;
  if (fOutputQA) delete fOutputQA;
  if (fOutputList) delete fOutputList;
  if (fArraySkipDDaugh) delete fArraySkipDDaugh;
  if (fArrayTrk) delete fArrayTrk;
}

//______________________________| Assignment Operator
AliAnalysisHFCorrOnFlySim& AliAnalysisHFCorrOnFlySim::operator=(const AliAnalysisHFCorrOnFlySim& c)
{
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fMcEvent = c.fMcEvent;
    fMcHandler = c.fMcHandler;
    fHistEventsProcessed = c.fHistEventsProcessed;
    fOutputQA = c.fOutputQA;
    fOutputList = c.fOutputList;
    fEtaMin = c.fEtaMin;
    fEtaMax = c.fEtaMax;
    fYMin = c.fYMin;
    fYMax = c.fYMax;
    fPtMin = c.fPtMin;
    fPtMax = c.fPtMax;
    fMinMultiplicity = c.fMinMultiplicity;
    fMaxMultiplicity = c.fMaxMultiplicity;
    fCorrPart1 = c.fCorrPart1;
    fCorrPart2 = c.fCorrPart2;
    fChargeSel1 = c.fChargeSel1;
    fChargeSel2 = c.fChargeSel2;
    fStack = c.fStack;
    fParticleArray = c.fParticleArray;
    fIsEventProp = c.fIsEventProp;
    fIsPartProp = c.fIsPartProp;
    fIsCorrOfQQbar = c.fIsCorrOfQQbar;
    fIsCorrOfHeavyFlavor = c.fIsCorrOfHeavyFlavor;
    fIsCorrOfHadronHadron = c.fIsCorrOfHadronHadron;
    fArraySkipDDaugh = c.fArraySkipDDaugh;
    fArrayTrk = c.fArrayTrk;
    flastdaugh = c.flastdaugh;
  }
    
  return *this;
}


//______________________________| User Output
void AliAnalysisHFCorrOnFlySim::UserCreateOutputObjects()
{  
  Info("AliAnalysisHFCorrOnFlySim","CreateOutputObjects of task %s", GetName());
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  fOutputQA = new TList();
  fOutputQA->SetOwner(kTRUE);
    
  fParticleArray = new TObjArray;
  fParticleArray->SetOwner(kTRUE);
    
  DefineHistoNames();

  if(!fArraySkipDDaugh)fArraySkipDDaugh=new TArrayI(200);
  if(!fArrayTrk)fArrayTrk=new TArrayI(5000);

  PostData(1, fOutputQA);
  PostData(2, fOutputList);
}

//______________________________| Init
void AliAnalysisHFCorrOnFlySim::Init()
{
  if(fDebug > 1) printf("AliAnalysisHFCorrOnFlySim::Init() \n");
  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}


//______________________________| User Exec
void AliAnalysisHFCorrOnFlySim::UserExec(Option_t *)
{

  if(fDebug > 1) printf("AliAnalysisHFCorrOnFlySim::UserExec \n");
  Init();
  
  if(fMcHandler){
    fMcEvent = fMcHandler->MCEvent();
  }else{
    if(fDebug > 1) printf("AliAnalysisHFCorrOnFlySim::Handler() fMcHandler=NULL\n");
    return;
  }
  if(!fMcEvent){
    if(fDebug > 1) printf("AliAnalysisHFCorrOnFlySim::UserExec()   fMcEvent=NULL \n");
    return;
  }
  
  fStack = ((AliMCEvent*)fMcEvent)->Stack();

  Bool_t IsEventMCSelected = IsMCEventSelected(fMcEvent);
  if(!IsEventMCSelected)return;
  if(fIsEventProp)CalculateEventProperties(fMcEvent);
  if(fIsPartProp)CalculateParticleProperties(fMcEvent);
      
  for(Int_t iPart = 0; iPart < fMcEvent->GetNumberOfTracks(); iPart++){
    
    AliVParticle* mcPart = (AliVParticle*)fMcEvent->GetTrack(iPart);
    if (!mcPart)continue;
    fHistEventsProcessed->Fill(3.5);
    
    Bool_t IsParticleMCSelected = IsMCParticleGenerated(mcPart);
    if(!IsParticleMCSelected)continue;
    fHistEventsProcessed->Fill(4.5);
    
    Bool_t IsParticleMCKineAccepted = IsMCParticleInKineCriteria(mcPart);
    if(!IsParticleMCKineAccepted)continue;
    fHistEventsProcessed->Fill(5.5);
    
    //Storing part array after event+part selections but right now its not used
    fParticleArray->Add(mcPart);  
  }
  
  if(fIsCorrOfHeavyFlavor)CalculateHFHadronCorrelations();
  if(fIsCorrOfQQbar)CalculateQQBarCorrelations();
  if(fIsCorrOfHadronHadron)CalculateHadronHadronCorrelations(fParticleArray);
  
  PostData(1, fOutputQA);
  PostData(2, fOutputList);
  
  return;
}

//______________________________| Terminate
void AliAnalysisHFCorrOnFlySim::Terminate(Option_t*)
{  
  Info("Terminate","Start and end of Method");
  
  fOutputQA = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputQA) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputList) {
    Printf("ERROR: Output list not available");
    return;
  }
  
  return;  
}

//______________________________| Event Cut
Bool_t AliAnalysisHFCorrOnFlySim::IsMCEventSelected(TObject* obj)
{
  if(!obj)return kFALSE;
  AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
  if(!event)return kFALSE;
    
  Bool_t isSelected = kTRUE;
  fHistEventsProcessed->Fill(0.5);

  const AliVVertex *vtxMC = event->GetPrimaryVertex();
  if(!vtxMC)return kFALSE;
  if(TMath::Abs(vtxMC->GetZ()) > 10.00)return kFALSE;
  else fHistEventsProcessed->Fill(1.5);

  Int_t PartMultiplicity = event->GetNumberOfTracks();
  if(PartMultiplicity < fMinMultiplicity || PartMultiplicity > fMaxMultiplicity)return kFALSE;
  else fHistEventsProcessed->Fill(2.5);

  return isSelected;
}

//______________________________| Event Properties
void AliAnalysisHFCorrOnFlySim::CalculateEventProperties(TObject* obj){
    
    if(!obj)return;
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
    if(!event) return;
    
    //Zvtx
    const AliVVertex *vtxMC = event->GetPrimaryVertex();
    if(!vtxMC) return;
    Float_t zVtx = vtxMC->GetZ();
    ((TH1F*)fOutputQA->FindObject(Form("fHistZvtx")))->Fill(zVtx);
    
    //Multiplicities Distributions
    //1. all
    Int_t MultTot = event->GetNumberOfTracks(); //fStack->GetNtrack();
    ((TH1F*)fOutputQA->FindObject(Form("fMultTot")))->Fill(MultTot);
    
    //2. all Primary
    Int_t MultPrimary = fStack->GetNprimary();
    ((TH1F*)fOutputQA->FindObject(Form("fMultPrimary")))->Fill(MultPrimary);
    
    //3. PhyPrimary and NOT PhyPrimary (eta<1.0)
    Int_t MultPrmyPartYes=0, MultPrmyPartYesCharge=0, MultPrmyPartNO=0;
   
    for(Int_t jAss=0; jAss<event->GetNumberOfTracks(); jAss++){

        AliVParticle *mcPartN=(AliVParticle*)event->GetTrack(jAss);
        if(!mcPartN) return;
        Int_t pdgAss = TMath::Abs(mcPartN->PdgCode());
        if(TMath::Abs(mcPartN->Eta())>1.0)continue; // Mult

        if(event->IsPhysicalPrimary(jAss)){
           ((TH1F*)fOutputQA->FindObject(Form("fHistNCharge")))->Fill(0.5);
           if(mcPartN->Charge()>0)((TH1F*)fOutputQA->FindObject(Form("fHistNCharge")))->Fill(1.5);
           else if(mcPartN->Charge()<0)((TH1F*)fOutputQA->FindObject(Form("fHistNCharge")))->Fill(2.5);
           else if(mcPartN->Charge()==0)((TH1F*)fOutputQA->FindObject(Form("fHistNCharge")))->Fill(3.5);
            
           MultPrmyPartYes++;
           if(pdgAss==11||pdgAss==13||pdgAss==211||pdgAss==321||pdgAss==2212)MultPrmyPartYesCharge++;
        }
        else MultPrmyPartNO++;
    }

    if(MultPrmyPartYes!=0)((TH1F*)fOutputQA->FindObject(Form("fMultPhyPriPart")))->Fill(MultPrmyPartYes);
    if(MultPrmyPartNO!=0)((TH1F*)fOutputQA->FindObject(Form("fMultNotPhyPriPart")))->Fill(MultPrmyPartNO);
    if(MultPrmyPartYesCharge!=0)((TH1F*)fOutputQA->FindObject(Form("fMultPhyPriPartCharge")))->Fill(MultPrmyPartYesCharge);

}


//______________________________| Generated Particle Cut
Bool_t AliAnalysisHFCorrOnFlySim::IsMCParticleGenerated(TObject* obj)
{
  // Check generated particles (primary, charged, pdgcode)
  // Add more related to issue..
  if(!obj) return kFALSE;
  if(!obj->InheritsFrom("AliVParticle")) {
    printf("object must derived from AliVParticle !");
    return kFALSE;
  }
  
  Bool_t isSelected = kTRUE;
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if(!particle) return kFALSE;
  
  //if(!fstack->IsPhysicalPrimary(particle->GetLabel())){
    //if(particle->PdgCode()!=11)
    //isSelected = kFALSE;
  //}
    
  return isSelected;
}

//______________________________| Kinematic Particle Cut
Bool_t AliAnalysisHFCorrOnFlySim::IsMCParticleInKineCriteria(TObject *obj)
{
  //Check generated particles (eta, y, pt)
  if(!obj) return  kFALSE;
  if(!obj->InheritsFrom("AliVParticle")){
    printf("object must derived from AliVParticle !");
    return kFALSE;   
  }
  
  Bool_t isSelected = kTRUE;
  AliVParticle* particle = dynamic_cast<AliVParticle *>(obj);
  if(!particle) return kFALSE;
  
  if(particle->Pt()<fPtMin || particle->Pt()>fPtMax) isSelected = kFALSE;
  if(particle->Eta()<fEtaMin || particle->Eta()>fEtaMax) isSelected = kFALSE;
  //if(particle->Y()<fYMin || particle->Y()>fYMax) isSelected = kFALSE; //Prblem now..
  //Double_t energy = particle->E(), pz = particle->Pz();
  //Double_t particleY = (TMath::Abs(energy-pz)>1e-10) ?  0.5*TMath::Log( (energy+pz)/(energy-pz) ) : 1e6;
    
  return isSelected;
}

//______________________________| Particle Properties
void AliAnalysisHFCorrOnFlySim::CalculateParticleProperties(TObject* obj){
    
    if(!obj)return;
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
    if(!event) return;
    Double_t MCPartProperties[4] = {-999};

    for (Int_t jPart = 0; jPart < event->GetNumberOfTracks(); jPart++){
       
        AliVParticle *mcPartN=(AliVParticle*)event->GetTrack(jPart);
        if(!mcPartN)continue;
        
        Bool_t IsParticleMCSelected = IsMCParticleGenerated(mcPartN);
        if(!IsParticleMCSelected)continue;
        
        Bool_t IsParticleMCKineAccepted = IsMCParticleInKineCriteria(mcPartN);
        if(!IsParticleMCKineAccepted)continue;
        
        MCPartProperties[1] = mcPartN->Pt();
        MCPartProperties[2] = mcPartN->Phi();
        MCPartProperties[3] = mcPartN->Eta();
        
        Int_t PDGv = TMath::Abs(mcPartN->PdgCode());
        
        if(PDGv==1||PDGv==2||PDGv==3)MCPartProperties[0] = 0.5;
        else if(PDGv==4)    MCPartProperties[0] = 1.5;
        else if(PDGv==5)    MCPartProperties[0] = 2.5;
        else if(PDGv==421||PDGv==411||PDGv==413){
             if(PDGv==421)   MCPartProperties[0] = 3.5;
             if(PDGv==411)   MCPartProperties[0] = 4.5;
             if(PDGv==413)   MCPartProperties[0] = 5.5;
        }
        else if(PDGv==511||PDGv==521){
            if(PDGv==511)   MCPartProperties[0] = 6.5;
            if(PDGv==521)   MCPartProperties[0] = 7.5;
        }
        else if(PDGv==11||PDGv==13||PDGv==211||PDGv==321||PDGv==2212){

            if(event->IsPhysicalPrimary(jPart)){
            if(PDGv==211)   MCPartProperties[0] = 8.5;
            if(PDGv==321)   MCPartProperties[0] = 9.5;
            if(PDGv==2212)  MCPartProperties[0] = 10.5;
            if(PDGv==11)    MCPartProperties[0] = 11.5;
            if(PDGv==13)    MCPartProperties[0] = 12.5;
            ((TH1F*)fOutputQA->FindObject(Form("fPhiDist")))->Fill(mcPartN->Phi());
            ((TH1F*)fOutputQA->FindObject(Form("fPhiEta")))->Fill(mcPartN->Eta());
            ((TH1F*)fOutputQA->FindObject(Form("fPhiY")))->Fill(mcPartN->Y());

            }
            else {
                if(PDGv==211)   MCPartProperties[0] = 13.5;
                if(PDGv==321)   MCPartProperties[0] = 14.5;
                if(PDGv==2212)  continue;//MCPartProperties[0] = 15.5; //beam particle
                if(PDGv==11)    MCPartProperties[0] = 16.5;
                if(PDGv==13)    MCPartProperties[0] = 17.5;
            }
        }
        else if(event->IsPhysicalPrimary(jPart))MCPartProperties[0] = 18.5;
        else  MCPartProperties[0] = 19.5;
        
        ((THnSparseD*)fOutputQA->FindObject(Form("fHistNParticle")))->Fill(MCPartProperties);
        //Comment --> Beam Proton giving peak in Phi at = Pi();
    }
}


//______________________________| HF-Hadron Correlations
void AliAnalysisHFCorrOnFlySim::CalculateHFHadronCorrelations(){
    
  for (Int_t jPart = 0; jPart < fMcEvent->GetNumberOfTracks(); jPart++){
        
      AliVParticle* mcPartHF = (AliVParticle*)fMcEvent->GetTrack(jPart);
      if (!mcPartHF)continue;
      Int_t pdgHF = TMath::Abs(mcPartHF->PdgCode());
      
      if(pdgHF==421 || pdgHF==411 || pdgHF==413){ ///D0, D+, D*
	
	if(mcPartHF->Pt() <2 || mcPartHF->Pt() >20)continue;
	if(TMath::Abs(mcPartHF->Y()) > 0.5)continue;
	fArraySkipDDaugh->Reset(0);
	fArraySkipDDaugh->AddAt(jPart,0);
	flastdaugh=1;
	RemoveNDaughterParticleArray(mcPartHF);
	HeavyFlavourCorrelations(mcPartHF);
      }
  }
}

//______________________________| Remove decay particles of trigger
void AliAnalysisHFCorrOnFlySim::RemoveNDaughterParticleArray(TObject* obj){
    
    if(!obj)return;
    AliVParticle* TrgPart = (AliVParticle*)obj;
    if(!TrgPart)return;
    
    Int_t DauPosF = TrgPart->GetFirstDaughter();
    Int_t DauPosL = TrgPart->GetLastDaughter();
    
    if(DauPosL<0)DauPosL=DauPosF;
    if(DauPosF > 0){
        for(Int_t jd = DauPosF; jd <= DauPosL; jd++){
            AliVParticle* Trg_DouPart = (AliVParticle*)fMcEvent->GetTrack(jd);
            if (!Trg_DouPart)continue;
            Int_t PDG_DouPart = TMath::Abs(Trg_DouPart->PdgCode());

            if(PDG_DouPart==11||PDG_DouPart==13||PDG_DouPart==211||PDG_DouPart==321||PDG_DouPart==2212){
                fArraySkipDDaugh->AddAt(jd,flastdaugh);
                flastdaugh++;
                RemoveNDaughterParticleArray(Trg_DouPart);
            }
            else {
                RemoveNDaughterParticleArray(Trg_DouPart);
            }
        }
    }
}


//______________________________| HF-Correlations Calculations
void AliAnalysisHFCorrOnFlySim::HeavyFlavourCorrelations(TObject *obj){
    
    if(!obj) return;
    AliVParticle* TrigPart = (AliVParticle*)obj;
    if(!TrigPart) return;
    
    Int_t PDG_TrigPart  = TMath::Abs(TrigPart->PdgCode());
         if(PDG_TrigPart==421)    PDG_TrigPart = 1; //D0
    else if(PDG_TrigPart==411)    PDG_TrigPart = 2; //D+
    else if(PDG_TrigPart==413)    PDG_TrigPart = 3; //D*
    else if(PDG_TrigPart==11)     PDG_TrigPart = 4; //e
    else                          PDG_TrigPart = 0; //NULL
    
    Bool_t hasToSkip = kFALSE;
    Int_t softpi = -1;
    Int_t TrgMomPos  = TrigPart->GetMother();
    
    if(TrgMomPos > 0){
        
        AliVParticle *MotherOfTrg = (AliVParticle*)fMcEvent->GetTrack(TrgMomPos);
        if(!MotherOfTrg) return;
        Int_t pdgOfMother = TMath::Abs(MotherOfTrg->PdgCode());
        
      
        if(PDG_TrigPart==4){ //Electron Stuff
            if(!(pdgOfMother==5||pdgOfMother==4||(400<pdgOfMother&&pdgOfMother<600)||(4000<pdgOfMother&&pdgOfMother<6000))){// NON HF electron
                if(pdgOfMother==22) PDG_TrigPart=5; //GAMMA CONV
                else if(pdgOfMother==111||pdgOfMother==113||pdgOfMother==221||pdgOfMother==223) PDG_TrigPart=6; //DALITZ DECAY
                else PDG_TrigPart=7; // Others NonHF
            }
        }
        
        if(pdgOfMother==413){ // Dstar --> D0: check soft pion
            if(PDG_TrigPart==1){
                for(Int_t isp = MotherOfTrg->GetFirstDaughter(); isp <= MotherOfTrg->GetLastDaughter(); isp++){
                    AliVParticle *sfp=(AliVParticle*)fMcEvent->GetTrack(isp);
                    Int_t pdgsp=TMath::Abs(sfp->PdgCode());
                    if(pdgsp==211)softpi=isp;
                }
            }
            TrgMomPos        =  MotherOfTrg->GetMother();
            MotherOfTrg      =  (AliVParticle*)fMcEvent->GetTrack(TrgMomPos);
            pdgOfMother      =  TMath::Abs(MotherOfTrg->PdgCode());
            
        }else if(pdgOfMother==423){// D*0 -> D0 or D+
            TrgMomPos        =  MotherOfTrg->GetMother();
            MotherOfTrg      =  (AliVParticle*)fMcEvent->GetTrack(TrgMomPos);
            pdgOfMother      =  TMath::Abs(MotherOfTrg->PdgCode());
        }
        
       if(pdgOfMother==5||(500<pdgOfMother&&pdgOfMother<600)||(5000<pdgOfMother&&pdgOfMother<6000)){
            PDG_TrigPart*=-1;
        }else if(pdgOfMother==4){
            TrgMomPos = MotherOfTrg->GetMother();
            if(TrgMomPos >= 0){
                MotherOfTrg   =  (AliVParticle*)fMcEvent->GetTrack(TrgMomPos);
                pdgOfMother   =  TMath::Abs(MotherOfTrg->PdgCode());
                if(pdgOfMother==5||(500<pdgOfMother&&pdgOfMother<600)||(5000<pdgOfMother&&pdgOfMother<6000)){
                    PDG_TrigPart*=-1;
                }
            }
        }
    }
  
   
    Double_t phiTrig = TrigPart->Phi();
    Double_t etaTrig = TrigPart->Eta();
    Double_t ptTrig  = TrigPart->Pt();
    
    Double_t PartProperties[8] = {static_cast<Double_t>(PDG_TrigPart),ptTrig,etaTrig,0,0,0,0,0};
    ((THnSparseD*)fOutputList->FindObject(Form("HFTrgiggerProp")))->Fill(PartProperties);
 
    Int_t nPartAss = 0;
    TArrayI* fArrayAssoPart = CalculateNPartType("Charge", nPartAss, 0);
    if(!fArrayAssoPart)return;
    
    for(Int_t ipartAsso = 0; ipartAsso < nPartAss; ipartAsso++){
    
        AliVParticle *partAsso  = (AliVParticle*)fMcEvent->GetTrack(fArrayAssoPart->At(ipartAsso));
        if(!partAsso)continue;
        hasToSkip=kFALSE;
        
        for(Int_t jdskip = 1; jdskip < flastdaugh; jdskip++){
            if(fArrayAssoPart->At(ipartAsso)==fArraySkipDDaugh->At(jdskip)){
                hasToSkip=kTRUE;
                break;
            }
        }
        
        if(hasToSkip)continue;
        
        PartProperties[3] = partAsso->Pt();
        PartProperties[4] = partAsso->Eta();
        PartProperties[5] = AssignCorrectPhiRange(phiTrig-partAsso->Phi());
        PartProperties[6] = etaTrig-partAsso->Eta();
        UInt_t pdgassmum  = TMath::Abs(partAsso->PdgCode());
        
        if(fArrayAssoPart->At(ipartAsso) == softpi)PartProperties[7]=-1;
        else if(pdgassmum==13)  PartProperties[7]=2;
        else if(pdgassmum==211) PartProperties[7]=3;
        else if(pdgassmum==321) PartProperties[7]=4;
        else if(pdgassmum==2212)PartProperties[7]=5;
        else if(pdgassmum==11){ //e-
            
            if(fStack->IsPhysicalPrimary(partAsso->GetLabel()))PartProperties[7]=1;//Prim e-

            Int_t eleMother = partAsso->GetMother();
            if(eleMother>0){
                
                AliVParticle *partEleMum=(AliVParticle*)fMcEvent->GetTrack(eleMother);
                if(!partEleMum)continue;
                UInt_t pdgeleMother = TMath::Abs(partEleMum->PdgCode());
                
                if(pdgeleMother==111||pdgeleMother==113||pdgeleMother==221||pdgeleMother==223)PartProperties[7]=6;// DALITZ DECAY
                else if(pdgeleMother==22)PartProperties[7]=7;// GAMMA CONV
                else if(pdgeleMother==4||pdgeleMother==5||(pdgeleMother>400&&pdgeleMother<600)||(pdgeleMother>4000&&pdgeleMother<6000))PartProperties[7]=8;
            }
            else PartProperties[7]=9;
        }
        else PartProperties[7]=0;
        
        Double_t ptLim_Sparse = ((THnSparseD*)fOutputList->FindObject(Form("2PCorrBtwn_HF-hadron")))->GetAxis(3)->GetXmax();
        if(PartProperties[3] > ptLim_Sparse) PartProperties[3] = ptLim_Sparse - 0.01;
        
        ((THnSparseD*)fOutputList->FindObject(Form("2PCorrBtwn_HF-hadron")))->Fill(PartProperties);
        
    }
}

//______________________________| HH Calculations
void AliAnalysisHFCorrOnFlySim::CalculateHadronHadronCorrelations(TObjArray *ParticleArray){
    
    //Under Development !
    //Bool_t hasToSkip=kFALSE;
    Int_t nPartCharge =0;
    TArrayI* fArrayChargePart = CalculateNPartType("Charge", nPartCharge, 0);
    if(!fArrayChargePart)return;
    
    for(Int_t ipart1 = 0; ipart1 < nPartCharge ; ipart1++){
        
        AliVParticle *partTrgCharge=(AliVParticle*)ParticleArray->At(fArrayChargePart->At(ipart1));
        if(!partTrgCharge)return;
        
        Int_t pdgTrigCharge  =TMath::Abs(partTrgCharge->PdgCode());
        if(pdgTrigCharge  == 11)        pdgTrigCharge = 0;
        else if(pdgTrigCharge == 13)    pdgTrigCharge = 1;
        else if(pdgTrigCharge == 211)   pdgTrigCharge = 2;
        else if(pdgTrigCharge == 321)   pdgTrigCharge = 3;
        else if(pdgTrigCharge == 2212)  pdgTrigCharge = 4;
        else pdgTrigCharge =-1;
        
        Double_t etaTrigCharge = partTrgCharge->Pt();
        Double_t phiTrigCharge = partTrgCharge->Phi();
        Double_t pTTrigCharge  = partTrgCharge->Eta();
        
        Double_t PropertiesOfTrigCharge[3] = {static_cast<Double_t>(pdgTrigCharge), pTTrigCharge, etaTrigCharge};//
        ((THnSparseD*)fOutputList->FindObject(Form("ChargeTrgiggerProp")))->Fill(PropertiesOfTrigCharge);
        
        for(Int_t ipart2 = 0; ipart2 < nPartCharge ; ipart2++){
            
            if(ipart2==ipart1)continue;
            
            AliVParticle *partAssoCharge=(AliVParticle*)ParticleArray->At(fArrayChargePart->At(ipart2));
            if(!partAssoCharge)return;
            
            Double_t etaAssoCharge = partAssoCharge->Pt();
            Double_t phiAssoCharge = partAssoCharge->Phi();
            Double_t pTAssoCharge  = partAssoCharge->Eta();
            
            Double_t DeltaPhi = AssignCorrectPhiRange(phiTrigCharge - phiAssoCharge);
            Double_t DeltaEta = etaTrigCharge - pTAssoCharge;
            
            Double_t HHCorr[8] = {static_cast<Double_t>(pdgTrigCharge), pTTrigCharge, etaTrigCharge, pTAssoCharge, etaAssoCharge, DeltaPhi, DeltaEta, 0};//
            
            Double_t ptLim_Sparse = ((THnSparseD*)fOutputList->FindObject(Form("2PCorrBtwn_%s-%s","TrgCharge","hadron")))->GetAxis(3)->GetXmax();
            if(pTAssoCharge > ptLim_Sparse) pTAssoCharge = ptLim_Sparse - 0.01;
            HHCorr[3] = pTAssoCharge;
            
            HHCorr[7]=TMath::Abs(partAssoCharge->PdgCode());
            if(HHCorr[7]==11)         HHCorr[7]=0;
            else if(HHCorr[7]==13)    HHCorr[7]=1;
            else if(HHCorr[7]==211)   HHCorr[7]=2;
            else if(HHCorr[7]==321)   HHCorr[7]=3;
            else if(HHCorr[7]==2212)  HHCorr[7]=4;
            else HHCorr[7]=-1;
            
            ((THnSparseD*)fOutputList->FindObject(Form("2PCorrBtwn_TrgCharge-hadron")))->Fill(HHCorr);
        }
    }
}


//______________________________| QQbar Calculations
void AliAnalysisHFCorrOnFlySim::CalculateQQBarCorrelations(){
    
    Int_t nPart1 =0, nPart2 = 0;
    
    TArrayI* fArray1Part = CalculateNPartType(fCorrPart1.Data(), nPart1, fChargeSel1);
    TArrayI* fArray1Part1 = new TArrayI(*fArray1Part);
    
    if(nPart1>0){
    
    TArrayI* fArray2Part = CalculateNPartType(fCorrPart2.Data(), nPart2, fChargeSel2);
    TArrayI* fArray2Part2 = new TArrayI(*fArray2Part);

    if(nPart2>0){
    
    for(Int_t ipart1 = 0; ipart1 < nPart1 ; ipart1++){
                
            AliVParticle *part1=(AliVParticle*)fMcEvent->GetTrack(fArray1Part1->At(ipart1));
            ((TH1F*)fOutputList->FindObject(Form("fHistPhiPart1")))->Fill(part1->Phi());
        
            for(Int_t ipart2 = 0; ipart2 < nPart2 ; ipart2++){
                AliVParticle *part2=(AliVParticle*)fMcEvent->GetTrack(fArray2Part2->At(ipart2));
                if(ipart2==0){((TH1F*)fOutputList->FindObject(Form("fHistPhiPart2")))->Fill(part2->Phi());
                }
                Double_t DeltaPhi = AssignCorrectPhiRange(part1->Phi() - part2->Phi());
                Double_t DeltaEta = part1->Eta() - part2->Eta();
                Double_t CorrVar[2] = {DeltaEta, DeltaPhi};
                ((THnSparseD*)fOutputList->FindObject(Form("2PCorrBtwn_%s-%s",fCorrPart1.Data(),fCorrPart2.Data())))->Fill(CorrVar);
            }
          }
       }
        delete fArray2Part2;
    }
    delete fArray1Part1;
}


//______________________________| Particle Type Calculations
TArrayI* AliAnalysisHFCorrOnFlySim::CalculateNPartType(TString pname, Int_t &count, Int_t ChargeSel){
    
    count=0;
    fArrayTrk->Reset(0);
    
    for(Int_t jAss=0; jAss<fMcEvent->GetNumberOfTracks(); jAss++){
        
        AliVParticle *partAss=(AliVParticle*)fMcEvent->GetTrack(jAss);
        
        Bool_t IsParticleMCSelected = IsMCParticleGenerated(partAss);
        if(!IsParticleMCSelected)continue;
        
        Bool_t IsParticleMCKineAccepted = IsMCParticleInKineCriteria(partAss);
        if(!IsParticleMCKineAccepted)continue;
        
        Int_t pdgAss=partAss->PdgCode();
        
        if(ChargeSel == -1 && pdgAss>0)continue;
        else if(ChargeSel == 1 && pdgAss<0)continue;
        
        if(pname == "c" && (TMath::Abs(pdgAss)==4)){
            fArrayTrk->AddAt(jAss,count);
            count++;
        }else if(pname == "b" && (TMath::Abs(pdgAss)==5)){
            fArrayTrk->AddAt(jAss,count);
            count++;
        }else if(pname == "D" && (TMath::Abs(pdgAss)==411||TMath::Abs(pdgAss)==413||TMath::Abs(pdgAss)==421)){
            fArrayTrk->AddAt(jAss,count);
            count++;
        }else if(pname == "B" && (TMath::Abs(pdgAss)==511||TMath::Abs(pdgAss)==521)){
            fArrayTrk->AddAt(jAss,count);
            count++;
        }else if(pname == "Charge" && (TMath::Abs(pdgAss)==11||TMath::Abs(pdgAss)==13||TMath::Abs(pdgAss)==211||TMath::Abs(pdgAss)==321||TMath::Abs(pdgAss)==2212)){
            if(!fStack->IsPhysicalPrimary(partAss->GetLabel()))continue;
            fArrayTrk->AddAt(jAss,count);
            count++;
        }
    }
    return fArrayTrk;
}



//______________________________| Definations of histo etc.
void AliAnalysisHFCorrOnFlySim::DefineHistoNames(){
    
    Double_t Pi = TMath::Pi();
    
    //1. Event Properties
    fHistEventsProcessed = new TH1I("fHistNEvents","fHistEventsProcessed",6,0,6) ;
    fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"All events");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Event with Ztx");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Event with Mult");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(4,"All Particle");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(5,"Only Physical Particle");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(6,"After Kine Cuts");
    fOutputQA->Add(fHistEventsProcessed);
    
    TH1F        *HistogramQA[11];
    for ( Int_t i = 0; i < 11; i++)HistogramQA[i] = NULL;
    HistogramQA[0] = new TH1F("fHistZvtx", "ZVtx distribution", 30, -30., 30.);
    HistogramQA[0]->GetXaxis()->SetTitle("ZVtx");
    HistogramQA[0]->GetYaxis()->SetTitle("N_tot");
    HistogramQA[0]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[1] = new TH1F("fMultTot", "Mult Total", 100, 0., 2000.);
    HistogramQA[1]->GetXaxis()->SetTitle("Total Multiplicity");
    HistogramQA[1]->GetYaxis()->SetTitle("N_tot");
    HistogramQA[1]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[2] = new TH1F("fMultPrimary", "Mult Primary", 100, 0., 2000.);
    HistogramQA[2]->GetXaxis()->SetTitle("Primary Multiplicity");
    HistogramQA[2]->GetYaxis()->SetTitle("N_tot");
    HistogramQA[2]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[3] = new TH1F("fMultPhyPriPart", "Mult PhyPrimary", 200, 0, 400);
    HistogramQA[3]->GetXaxis()->SetTitle("PhyPrimary Multiplicity");
    HistogramQA[3]->GetYaxis()->SetTitle("N_Pri");
    HistogramQA[3]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[4] = new TH1F("fMultNotPhyPriPart", "Mult NotPhyPrimary", 200, 0, 400.);
    HistogramQA[4]->GetXaxis()->SetTitle("NotPhyPrimary Multiplicity");
    HistogramQA[4]->GetYaxis()->SetTitle("No_Pri");
    HistogramQA[4]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[5] = new TH1F("fMultPhyPriPartCharge", "Mult Charge PhyPrimary", 150, 0, 150);
    HistogramQA[5]->GetXaxis()->SetTitle("Chrg_PhyPrimaryMultiplicity");
    HistogramQA[5]->GetYaxis()->SetTitle("N_PhyPri_Charge");
    HistogramQA[5]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[6] = new TH1F("fHistNCharge","Particle Counts",4,0,4) ;
    HistogramQA[6]->GetXaxis()->SetBinLabel(1,"All Particle");
    HistogramQA[6]->GetXaxis()->SetBinLabel(2,"+ive Charge");
    HistogramQA[6]->GetXaxis()->SetBinLabel(3,"-ive Charge");
    HistogramQA[6]->GetXaxis()->SetBinLabel(4,"Neutral Charge");
    
    HistogramQA[7] = new TH1F("fRemoveDaugh", "Remove Daughters", 20, 0., 20);
    HistogramQA[7]->GetXaxis()->SetTitle("N Daughters");
    HistogramQA[7]->GetYaxis()->SetTitle("N_removed");
    HistogramQA[7]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[8] = new TH1F("fPhiDist", "Phi Dist", 45, -0.5*Pi, 2*Pi);
    HistogramQA[8]->GetXaxis()->SetTitle("#phi");
    HistogramQA[8]->GetYaxis()->SetTitle("N_pri");
    HistogramQA[8]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[9] = new TH1F("fPhiEta", "Eta Dist", 90, -9., 9.);
    HistogramQA[9]->GetXaxis()->SetTitle("#eta");
    HistogramQA[9]->GetYaxis()->SetTitle("N_pri");
    HistogramQA[9]->SetMarkerStyle(kFullCircle);
    
    HistogramQA[10] = new TH1F("fPhiY", "Y Dist", 90, -9., 9.);
    HistogramQA[10]->GetXaxis()->SetTitle("Y");
    HistogramQA[10]->GetYaxis()->SetTitle("N_pri");
    HistogramQA[10]->SetMarkerStyle(kFullCircle);
    
    
    for(Int_t i = 0; i <11; i++)
    {
        HistogramQA[i]->Sumw2();
        fOutputQA->Add(HistogramQA[i]);
    }
    
    //2. Particle Properties
    Int_t    fBins[4] = {20,  200,             18,    30 };
    Double_t  fMin[4] = {0,   0,    0*TMath::Pi(),   -15.};
    Double_t  fMax[4] = {20,  40,   2*TMath::Pi(),    15.};
    THnSparseD *HistogramQAPartProp   = new THnSparseD("fHistNParticle",  "Particle Properties",  4, fBins, fMin, fMax);
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(1, "LtQuarks");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(2, "Charm");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(3, "Beauty");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(4, "DZero");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(5, "DPlus");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(6, "DStar");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(7, "BZero");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(8, "BPlus");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(9, "PPions");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(10,"PKaons");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(11,"PProtons");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(12,"PElectoon");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(13,"PMuons");
    //Physics Primary..
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(14,"nPPions");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(15,"nPKaons");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(16,"nPProtons");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(17,"nPElectoon");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(18,"nPMuons");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(19,"Others PhyPri");
    HistogramQAPartProp->GetAxis(0)->SetBinLabel(20,"Nt PPryOthers");
    
    HistogramQAPartProp->GetAxis(1)->SetTitle("p_{T}");
    HistogramQAPartProp->GetAxis(2)->SetTitle("#phi");
    HistogramQAPartProp->GetAxis(3)->SetTitle("#eta");
    fOutputQA->Add(HistogramQAPartProp);
    
    
    //3. Correlations histograms
    //3a. Quark-AntiQuark Correlations
    if(fIsCorrOfQQbar){
        
        Int_t    fBinsCorr[2] = {20,     32}; //DEta, DPhi
        Double_t  fMinCorr[2] = {-2.0, -0.5*TMath::Pi()};
        Double_t  fMaxCorr[2] = {+2.0,  1.5*TMath::Pi()};
        
        TString Corr2PQQbar = "2PCorrBtwn_";
        Corr2PQQbar += fCorrPart1.Data();
        Corr2PQQbar += "-";
        Corr2PQQbar += fCorrPart2.Data();
        
        THnSparseD *Correlation2PQQbar    = new THnSparseD(Corr2PQQbar.Data(),"2PCorrQQbar", 2, fBinsCorr, fMinCorr, fMaxCorr);
        Correlation2PQQbar->Sumw2();
        fOutputList->Add(Correlation2PQQbar);
        
        TH1F *hPhiPart1 = new TH1F("fHistPhiPart1", "Phi_Dist_of_Part1", 45, -0.5*Pi, 2*Pi);
        TH1F *hPhiPart2 = new TH1F("fHistPhiPart2", "Phi_Dist_of_Part2", 45, -0.5*Pi, 2*Pi);
        hPhiPart1->Sumw2();
        hPhiPart2->Sumw2();
        fOutputList->Add(hPhiPart1);
        fOutputList->Add(hPhiPart2);
        
    }
    
    
    //3b. D-Hadron Correlations
    if(fIsCorrOfHeavyFlavor){
        
        Int_t     nbinsTrigHF[3] = {  15, 20,  20};
        Double_t binlowTrigHF[3] = {-7.5, 0., -2.};
        Double_t  binupTrigHF[3] = { 7.5, 20., 2.};
        
        Int_t     nbinsCorrHF[8] = {  15, 20,  20, 10,   30,               32,   20,    11};
        Double_t binlowCorrHF[8] = {-7.5, 0., -2., 0., -15., -0.5*TMath::Pi(),  -5.,  -1.5};
        Double_t  binupCorrHF[8] = { 7.5, 20., 2., 5.,  15.,  1.5*TMath::Pi(),   5.,   9.5};
        
        THnSparseD *trigDPartPr   = new THnSparseD("HFTrgiggerProp","fHFTrgiggerProp;pdg;ptTrig;etaTrig;",3,nbinsTrigHF,binlowTrigHF,binupTrigHF);
        THnSparseD *trigDPartCorr = new THnSparseD("2PCorrBtwn_HF-hadron","HFCorrelations;pdg;ptTrig;etaTrig;ptAss;etaAss;deltaPhi;deltaEta;pdgAss;",8,nbinsCorrHF,binlowCorrHF,binupCorrHF);
        
        trigDPartPr->Sumw2();
        trigDPartCorr->Sumw2();
        fOutputList->Add(trigDPartPr);
        fOutputList->Add(trigDPartCorr);
    }
    
    //3c. Hadron-Hadron Correlations
    if(fIsCorrOfHadronHadron){
        
        Int_t     nbinsTrigHH[3] = {   6, 30,  20 };
        Double_t binlowTrigHH[3] = {-1.5, 0., -2.};
        Double_t  binupTrigHH[3] = { 4.5, 30., 2.};
        
        Int_t     nbinsCorrHH[8] = {  6,  30,  20,  6,   30,               72,    20,    6};
        Double_t binlowCorrHH[8] = {-1.5, 0., -2., 0., -15., -0.5*TMath::Pi(),  -5,  -1.5};
        Double_t  binupCorrHH[8] = { 4.5, 30., 2., 3.,  15.,  1.5*TMath::Pi(),   5.,  4.5};
        
        THnSparseD *trigChargePartPr   = new THnSparseD("ChargeTrgiggerProp","fChargeTrgiggerProp;pdg;ptTrig;etaTrig;",3,nbinsTrigHH,binlowTrigHH,binupTrigHH);
        THnSparseD *trigHHPartCorr = new THnSparseD("2PCorrBtwn_TrgCharge-hadron","HHCorrelations;pdg;ptTrig;etaTrig;ptAss;etaAss;deltaPhi;deltaEta;pdgAss;",8,nbinsCorrHH,binlowCorrHH,binupCorrHH);
        trigChargePartPr->Sumw2();
        trigHHPartCorr->Sumw2();
        fOutputList->Add(trigChargePartPr);
        fOutputList->Add(trigHHPartCorr);
        
    }
    
    PostData(1, fOutputQA);
    PostData(2, fOutputList);
}

