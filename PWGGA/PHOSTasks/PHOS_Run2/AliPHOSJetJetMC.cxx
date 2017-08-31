#include <TObject.h>
#include <TClonesArray.h>
#include <TH2.h>
#include <AliLog.h>
#include <AliPHOSGeometry.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliVCaloTrigger.h>
#include <AliVCluster.h>
#include <AliPHOSJetJetMC.h>
#include <AliCaloPhoton.h>

#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"


#include "TParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAnalysisManager.h"

// Author: Daiki Sekihata (Hiroshima University)
ClassImp(AliPHOSJetJetMC)

//________________________________________________________________________
AliPHOSJetJetMC::AliPHOSJetJetMC():
  fPtHardBin(-1),
  fPtHard(0),
  fXsection(0),
  fNTrials(0),
  fPtHardAndJetPtFactor(3),
  fPtHardAndSinglePtFactor(1.5),
  fFirstJetIndex(-1),
  fLastJetIndex(-1),
  fGenJetID(-1)
{
  //Constructor
  

}
//________________________________________________________________________
AliPHOSJetJetMC::AliPHOSJetJetMC(Int_t pThardbin):
  fPtHardBin(-1),
  fPtHard(0),
  fXsection(0),
  fNTrials(0),
  fPtHardAndJetPtFactor(3),
  fPtHardAndSinglePtFactor(1.5),
  fFirstJetIndex(-1),
  fLastJetIndex(-1),
  fGenJetID(-1)
{
  //Constructor
  fPtHardBin = pThardbin; 
  AliInfo(Form("pT-hard-bin : %d",fPtHardBin));

}
//________________________________________________________________________
AliPHOSJetJetMC::~AliPHOSJetJetMC()
{



}
//________________________________________________________________________
void AliPHOSJetJetMC::ConfigureJetJetMC(AliVEvent *event)
{
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);

  AliMCEvent *mcevent = 0x0;
  AliGenPythiaEventHeader* pythiaGenHeader = 0x0;

  if(esd){
    AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    if(eventHandler){
      AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(eventHandler);
      if(mcEventHandler) mcevent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
    }
    pythiaGenHeader = GetPythiaEventHeader(mcevent);
  }
  else if(aod) pythiaGenHeader = GetPythiaEventHeader(aod);

  Int_t pdg=0;
  Double_t pT=0;
  Int_t firstindexJet = 0;
  Int_t lastindexJet = -1;
  Int_t genIDJet = -1;

  if(esd){
    AliStack *stack = (AliStack*)mcevent->Stack();
    if(!stack){
      AliError("Could not get MC Stack!");
      return;
    }

    TList *genHeaders = GetGenHeaderList(mcevent);
    Int_t Ngen = genHeaders->GetEntries();
    AliInfo(Form("Ngen = %d.",Ngen));
    for(Int_t igen=0;igen<Ngen;igen++){
      AliGenEventHeader *gh = (AliGenEventHeader*)genHeaders->At(igen);
      TString GeneratorName = gh->GetName();

      AliInfo(Form("GeneratorName = %s , NProduced = %d.",GeneratorName.Data(),gh->NProduced()));
      lastindexJet += gh->NProduced();

      if(GeneratorName.Contains("Jets",TString::kIgnoreCase) || GeneratorName.Contains("Pythia",TString::kIgnoreCase)){
        //this string supports LHC16h3, LHC16h4.
        genIDJet = igen;
        break;
      }
      firstindexJet = gh->NProduced();
    }

    AliInfo(Form("genID of jet = %d , firstindexJet = %d , lastindexJet = %d.",genIDJet,firstindexJet,lastindexJet));
    fFirstJetIndex = firstindexJet;
    fLastJetIndex = lastindexJet;
    fGenJetID = genIDJet;
 
    const Int_t Ntrack = stack->GetNtrack();
    AliInfo(Form("%d MC particles are generated in ESD.",Ntrack));

  }//end of ESD
  else if(aod){
    TClonesArray *MCArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!MCArray){
      AliError("Could not retrieve AOD event!");
      return;
    }

    TList *genHeaders = GetGenHeaderList(aod);
    Int_t Ngen = genHeaders->GetEntries();
    AliInfo(Form("Ngen = %d.",Ngen));
    for(Int_t igen=0;igen<Ngen;igen++){
      AliGenEventHeader *gh = (AliGenEventHeader*)genHeaders->At(igen);
      TString GeneratorName = gh->GetName();

      AliInfo(Form("GeneratorName = %s , NProduced = %d.",GeneratorName.Data(),gh->NProduced()));
      lastindexJet += gh->NProduced();

      if(GeneratorName.Contains("Jets",TString::kIgnoreCase) || GeneratorName.Contains("Pythia",TString::kIgnoreCase)){
        //this string supports LHC16h3, LHC16h4.
        //genIDJet = igen;
        genIDJet = igen;
        break;
      }
      firstindexJet = gh->NProduced();
    }

    AliInfo(Form("genID of jet = %d , firstindexJet = %d , lastindexJet = %d.",genIDJet,firstindexJet,lastindexJet));
    fFirstJetIndex = firstindexJet;
    fLastJetIndex = lastindexJet;
    fGenJetID = genIDJet;

    const Int_t Ntrack = MCArray->GetEntriesFast();
    AliInfo(Form("%d MC particles are generated in AOD.",Ntrack));

  }//end of AOD


}
//________________________________________________________________________
AliGenPythiaEventHeader* AliPHOSJetJetMC::GetPythiaEventHeader(AliVEvent *MCEvent)
{
  Int_t NTrials = -1;
  Float_t XSection = -1;

  AliGenCocktailEventHeader* cHeader = 0;
  AliAODMCHeader *cHeaderAOD = 0;
  Bool_t headerFound = kFALSE;

  if(MCEvent->IsA()==AliMCEvent::Class()){
    if(dynamic_cast<AliMCEvent*>(MCEvent)){
      cHeader                  = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
      if(cHeader) headerFound  = kTRUE;
    }
  }
  if(MCEvent->IsA()==AliAODEvent::Class()){ // MCEvent is a AODEvent in case of AOD
    cHeaderAOD                 = dynamic_cast<AliAODMCHeader*>(MCEvent->FindListObject(AliAODMCHeader::StdBranchName()));
    if(cHeaderAOD) headerFound = kTRUE;
  }

  if(headerFound){
    TList *genHeaders = 0x0;
    if(cHeader) genHeaders = cHeader->GetHeaders();
    if(cHeaderAOD){
      genHeaders = cHeaderAOD->GetCocktailHeaders();

    }
    AliGenEventHeader* gh = 0;
    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName   = gh->GetName();
      AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(gh);
      if(gPythia){
        //if dynamic_cast is successfully done, gPythia should not by 0x0.
        //this is a property of dynamic_cast.
        NTrials = gPythia->Trials();
        XSection = gPythia->GetXsection();

        return gPythia;
      }
    }

    //AliGenPythiaEventHeader are not found.
    NTrials = -1;
    XSection = -1;
    AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Pythia event header not found");
    return 0;

  }
  else {
    AliGenEventHeader * eventHeader = dynamic_cast<AliMCEvent*>(MCEvent)->GenEventHeader();
    TString eventHeaderName  = eventHeader->ClassName();
    if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0){
      AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader);
      NTrials = gPythia->Trials();
      XSection = gPythia->GetXsection();
      return gPythia;
    }
    else{
      NTrials = -1;
      XSection = -1;
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Pythia event header not found");
      return 0;
    }

  }

}
//_______________________________________________________________________________
Bool_t AliPHOSJetJetMC::ComparePtHardBin(AliVEvent *event)
{
  const Double_t pthardbin_loweredges[20]  = {5, 7,  9, 12, 16, 21, 28, 36, 45, 57, 70, 85,  99, 115, 132, 150, 169, 190, 212, 235};
  const Double_t pthardbin_higheredges[20] = {7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235,  -1};

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);

  AliMCEvent *mcevent = 0x0;
  AliGenPythiaEventHeader* pythiaGenHeader = 0x0;

  if(esd){
    AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    if(eventHandler){
      AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
      if(mcEventHandler) mcevent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
    }
    pythiaGenHeader = GetPythiaEventHeader(mcevent);
  }
  else if(aod) pythiaGenHeader = GetPythiaEventHeader(aod);

  Float_t xsection = pythiaGenHeader->GetXsection();
  Float_t pThard   = pythiaGenHeader->GetPtHard();

  if(pThard < pthardbin_loweredges[fPtHardBin-1]){
    //printf("pT-hard = %f GeV/c , this value is lower than edge of pT-hard-bin[%d] settings. reject!\n",pThard,fPtHardBin);
    AliInfo(Form("pT-hard = %f GeV/c , this value is lower than edge of pT-hard-bin[%d] settings. reject!",pThard,fPtHardBin));
    return kFALSE;
  }
  if(fPtHardBin < 20 && pthardbin_higheredges[fPtHardBin-1] < pThard){
    AliInfo(Form("pT-hard = %f GeV/c , this value is higher than edge of pT-hard-bin[%d] settings. reject!",pThard,fPtHardBin));
    return kFALSE;
  }

  return kTRUE;
}
//_______________________________________________________________________________
Bool_t AliPHOSJetJetMC::ComparePtHardWithJet(AliVEvent *event)
{
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);

  AliMCEvent *mcevent = 0x0;
  AliGenPythiaEventHeader* pythiaGenHeader = 0x0;

  if(esd){
    AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    if(eventHandler){
      AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
      if(mcEventHandler) mcevent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
    }
    pythiaGenHeader = GetPythiaEventHeader(mcevent);
  }
  else if(aod) pythiaGenHeader = GetPythiaEventHeader(aod);

  Float_t xsection = pythiaGenHeader->GetXsection();
  Float_t pThard   = pythiaGenHeader->GetPtHard();

  TParticle * jet = 0;
  Int_t nTriggerJets =  pythiaGenHeader->NTriggerJets();

  AliInfo(Form("Njets: %d, pT Hard %f",nTriggerJets, pThard));

  Float_t tmpjet[]={0,0,0,0};
  for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
  {
    pythiaGenHeader->TriggerJet(ijet, tmpjet);
    jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);

    AliInfo(Form("jet %d; pycell jet pT %f",ijet, jet->Pt()));

    //Compare jet pT and pt Hard
    if(jet->Pt() > fPtHardAndJetPtFactor * pThard)
    {
      AliInfo(Form("Reject jet event with : pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n", pThard, jet->Pt(), fPtHardAndJetPtFactor));
      if(jet) delete jet;
      return kFALSE;
    }
  }

  if(jet) delete jet;

  return kTRUE;

}
//_______________________________________________________________________________
Bool_t AliPHOSJetJetMC::ComparePtHardWithSingleParticle(AliVEvent *event)
{
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event);

  Int_t pdg=0;
  Double_t pT=0;
  Int_t firstindexJet = 0;
  Int_t lastindexJet = -1;
  Int_t genIDJet = -1;
  AliGenPythiaEventHeader* pythiaGenHeader = 0x0;

  if(esd){
    AliMCEvent *mcevent = 0x0;
    AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();

    if(eventHandler){
      AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(eventHandler);
      if(mcEventHandler) mcevent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
    }

    pythiaGenHeader = GetPythiaEventHeader(mcevent);
    Float_t xsection = pythiaGenHeader->GetXsection();
    Float_t pThard   = pythiaGenHeader->GetPtHard();

    AliStack *stack = (AliStack*)mcevent->Stack();
    if(!stack){
      AliError("Could not get MC Stack!");
      return kFALSE;
    }
 
    AliInfo(Form("genID of jet = %d , firstindexJet = %d , lastindexJet = %d.",fGenJetID,fFirstJetIndex,fLastJetIndex));
    const Int_t Ntrack = stack->GetNtrack();
    AliInfo(Form("%d MC particles are generated in ESD.",Ntrack));
    for(Int_t i=fFirstJetIndex;i<=fLastJetIndex;i++){
      TParticle *p = (TParticle*)stack->Particle(i);
      pdg = p->GetPdgCode();

      if(pdg==111 || pdg==221 || pdg==22){
        pT = p->Pt();
        if(pT > fPtHardAndSinglePtFactor*pThard){
          AliInfo(Form("Reject jet event with : pT Hard %2.2f, particle %d : pT %2.2f, rejection factor %1.1f\n", pThard, pdg, pT, fPtHardAndSinglePtFactor));
          return kFALSE;//reject this event
        }

      }//end of particle selection

    }
  }//end of ESD
  else if(aod){
    pythiaGenHeader = GetPythiaEventHeader(aod);
    Float_t xsection = pythiaGenHeader->GetXsection();
    Float_t pThard   = pythiaGenHeader->GetPtHard();

    TClonesArray *MCArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!MCArray){
      AliError("Could not retrieve AOD event!");
      return kFALSE;
    }
    AliInfo(Form("genID of jet = %d , firstindexJet = %d , lastindexJet = %d.",fGenJetID,fFirstJetIndex,fLastJetIndex));

    const Int_t Ntrack = MCArray->GetEntriesFast();
    AliInfo(Form("%d MC particles are generated in AOD.",Ntrack));
    for(Int_t i=fFirstJetIndex;i<=fLastJetIndex;i++){
      AliAODMCParticle *p = (AliAODMCParticle*)MCArray->At(i);
      pdg = p->PdgCode();
      //another way to select particles generated by JJMC is to use AliAODMCParticle::GetGeneratorIndex()

      if(pdg==111 || pdg==221 || pdg==22){
        pT = p->Pt();
        if(pT > fPtHardAndSinglePtFactor*pThard){
          AliInfo(Form("Reject jet event with : pT Hard %2.2f, particle %d : pT %2.2f, rejection factor %1.1f", pThard, pdg, pT, fPtHardAndSinglePtFactor));
          return kFALSE;//reject this event
        }

      }//end of particle selection

    }//end of mc particle loop

  }//end of AOD

  return kTRUE;//accept this event

}
//________________________________________________________________________
TList *AliPHOSJetJetMC::GetGenHeaderList(AliVEvent *MCEvent)
{
  AliGenCocktailEventHeader* cHeader = 0;
  AliAODMCHeader *cHeaderAOD = 0;
  Bool_t headerFound = kFALSE;

  if(MCEvent->IsA()==AliMCEvent::Class()){
    if(dynamic_cast<AliMCEvent*>(MCEvent)){
      cHeader                  = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
      if(cHeader) headerFound  = kTRUE;
    }
  }
  if(MCEvent->IsA()==AliAODEvent::Class()){ // MCEvent is a AODEvent in case of AOD
    cHeaderAOD                 = dynamic_cast<AliAODMCHeader*>(MCEvent->FindListObject(AliAODMCHeader::StdBranchName()));
    if(cHeaderAOD) headerFound = kTRUE;
  }

  TList *genHeaders = 0x0;
  if(cHeader){
    genHeaders = cHeader->GetHeaders();
    return genHeaders;
  }
  if(cHeaderAOD){
    genHeaders = cHeaderAOD->GetCocktailHeaders();
    return genHeaders;
  }
  return genHeaders;
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

