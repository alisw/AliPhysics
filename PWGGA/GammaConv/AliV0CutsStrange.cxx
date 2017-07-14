#include "AliV0CutsStrange.h"

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "AliMCEvent.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliVCaloCells.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEMCALTriggerPatchInfo.h"

class iostream;

using namespace std;

ClassImp(AliV0CutsStrange)


const char* AliV0CutsStrange::fgkCutNames[AliV0CutsStrange::kNCuts] = {
  "V0FinderType"            // 0
};


//________________________________________________________________________
AliV0CutsStrange::AliV0CutsStrange(const char *name,const char *title) : 
  AliAnalysisCuts(name, title),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fIsQA(kFALSE),
  fV0ReaderStrangeName("V0ReaderStrange"),
  fCutString(NULL),
  fCutStringRead(""),
  fPreSelCut(kFALSE),
  fUseOnFlyV0Finder(kTRUE),
  fUseOnFlyV0FinderSameSign(0),
  fMinClsTPC(0.),
  fMinClsTPCToF(0.),
  fUseCorrectedTPCClsInfo(kFALSE),
  fPIDTPCnSigmaProtonLow(0),
  fPIDTPCnSigmaProtonUp(0),
  fPIDTPCnSigmaPionLow(0),
  fPIDTPCnSigmaPionUp(0),
  fPIDTOFnSigmaProtonLow(0),
  fPIDTOFnSigmaProtonUp(0),
  fPIDTOFnSigmaPionLow(0),
  fPIDTOFnSigmaPionUp(0),
  fHistoCutIndex(NULL),
  fHistodEdxCutsProton(NULL),
  fHistoTPCdEdxProtonBefore(NULL),
  fHistoTPCdEdxSigmaProtonBefore(NULL),
  fHistoTPCdEdxProtonAfter(NULL),
  fHistoTPCdEdxSigmaProtonAfter(NULL),
  fHistoTOFdEdxProtonBefore(NULL),
  fHistoTOFdEdxSigmaProtonBefore(NULL),
  fHistoTOFdEdxProtonAfter(NULL),
  fHistoTOFdEdxSigmaProtonAfter(NULL),
  fHistodEdxCutsPion(NULL),
  fHistoTPCdEdxPionBefore(NULL),
  fHistoTPCdEdxSigmaPionBefore(NULL),
  fHistoTPCdEdxPionAfter(NULL),
  fHistoTPCdEdxSigmaPionAfter(NULL),
  fHistoTOFdEdxPionBefore(NULL),
  fHistoTOFdEdxSigmaPionBefore(NULL),
  fHistoTOFdEdxPionAfter(NULL),
  fHistoTOFdEdxSigmaPionAfter(NULL)
{
  InitPIDResponse();
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());
}



//________________________________________________________________________
AliV0CutsStrange::AliV0CutsStrange(const AliV0CutsStrange &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fIsQA(kFALSE),
  fV0ReaderStrangeName("V0ReaderStrange"),
  fCutString(ref.fCutString),
  fCutStringRead(ref.fCutStringRead),
  fPreSelCut(ref.fPreSelCut),
  fUseOnFlyV0Finder(ref.fUseOnFlyV0Finder),
  fUseOnFlyV0FinderSameSign(ref.fUseOnFlyV0FinderSameSign),
  fMinClsTPC(ref.fMinClsTPC),
  fMinClsTPCToF(ref.fMinClsTPCToF),
  fUseCorrectedTPCClsInfo(ref.fUseCorrectedTPCClsInfo),
  fPIDTPCnSigmaProtonLow(ref.fPIDTPCnSigmaProtonLow),
  fPIDTPCnSigmaProtonUp(ref.fPIDTPCnSigmaProtonUp),
  fPIDTPCnSigmaPionLow(ref.fPIDTPCnSigmaPionLow),
  fPIDTPCnSigmaPionUp(ref.fPIDTPCnSigmaPionUp),
  fPIDTOFnSigmaProtonLow(ref.fPIDTOFnSigmaProtonLow),
  fPIDTOFnSigmaProtonUp(ref.fPIDTOFnSigmaProtonUp),
  fPIDTOFnSigmaPionLow(ref.fPIDTOFnSigmaPionLow),
  fPIDTOFnSigmaPionUp(ref.fPIDTOFnSigmaPionUp),
  fHistoCutIndex(NULL),
  fHistodEdxCutsProton(NULL),
  fHistoTPCdEdxProtonBefore(NULL),
  fHistoTPCdEdxSigmaProtonBefore(NULL),
  fHistoTPCdEdxProtonAfter(NULL),
  fHistoTPCdEdxSigmaProtonAfter(NULL),
  fHistoTOFdEdxProtonBefore(NULL),
  fHistoTOFdEdxSigmaProtonBefore(NULL),
  fHistoTOFdEdxProtonAfter(NULL),
  fHistoTOFdEdxSigmaProtonAfter(NULL),
  fHistodEdxCutsPion(NULL),
  fHistoTPCdEdxPionBefore(NULL),
  fHistoTPCdEdxSigmaPionBefore(NULL),
  fHistoTPCdEdxPionAfter(NULL),
  fHistoTPCdEdxSigmaPionAfter(NULL),
  fHistoTOFdEdxPionBefore(NULL),
  fHistoTOFdEdxSigmaPionBefore(NULL),
  fHistoTOFdEdxPionAfter(NULL),
  fHistoTOFdEdxSigmaPionAfter(NULL)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){
    fCuts[jj]=ref.fCuts[jj];
  }
  fCutString=new TObjString((GetCutNumber()).Data());
}



//________________________________________________________________________
AliV0CutsStrange::~AliV0CutsStrange() {
  // Destructor
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }
}


//________________________________________________________________________
void AliV0CutsStrange::InitCutHistograms(TString name, Bool_t preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);
  
  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("V0Cuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }
  
    // IsV0Selected
  fHistoCutIndex=new TH1F(Form("IsV0Selected %s",GetCutNumber().Data()),"IsV0Selected",12,-0.5,11.5);
  fHistoCutIndex->GetXaxis()->SetBinLabel(kV0In+1,"in");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoV0+1,"miss. V0 in AOD");
//   if (!fSwitchToKappa)fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"PID");
//   else fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"Kappa+[TOF,ITS,TRD] PID");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kSameSign+1,"Same sign");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kEventPlane+1,"EventPlane");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kV0Out+1,"out");
  fHistograms->Add(fHistoCutIndex);
  
  // dEdx Cuts proton
  fHistodEdxCutsProton=new TH2F(Form("dEdxCuts proton %s",GetCutNumber().Data()),"dEdxCuts proton vs p_{T,e}",11,-0.5,10.5,250,0,50);
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(1,"in");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(2,"TPCproton");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(3,"TOF proton");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(4,"out");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(5,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(6,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(7,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(8,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(9,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(10,"empty");
  fHistodEdxCutsProton->GetXaxis()->SetBinLabel(11,"empty");
  fHistograms->Add(fHistodEdxCutsProton);
    
  if(fIsQA==kTRUE){
    fHistoTPCdEdxProtonBefore=new TH2F(Form("Proton_dEdx_TPC_before %s",GetCutNumber().Data()),"dEdx Proton before; p_{T}; TPC signal" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxProtonBefore);
    fHistoTPCdEdxSigmaProtonBefore=new TH2F(Form("Proton_dEdxSigma_TPC_before %s",GetCutNumber().Data()),"dEdx Sigma Proton before; p_{T}; TPC sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTPCdEdxSigmaProtonBefore);
    fHistoTPCdEdxProtonAfter=new TH2F(Form("Proton_dEdx_TPC_after%s",GetCutNumber().Data()),"dEdx Proton after; p_{T}; TPC signal" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxProtonAfter);
    fHistoTPCdEdxSigmaProtonAfter=new TH2F(Form("Proton_dEdxSigma_TPC_after %s",GetCutNumber().Data()),"dEdx Sigma Proton after; p_{T}; TPC sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTPCdEdxSigmaProtonAfter);
    
    fHistoTOFdEdxProtonBefore=new TH2F(Form("Proton_dEdx_TOF_before %s",GetCutNumber().Data()),"dEdx Proton before; p_{T}; TOF signal" ,150,0.03,20,800,0,30000);
    fHistograms->Add(fHistoTOFdEdxProtonBefore);
    fHistoTOFdEdxSigmaProtonBefore=new TH2F(Form("Proton_dEdxSigma_TOF_before %s",GetCutNumber().Data()),"dEdx Sigma Proton before; p_{T}; TOF sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTOFdEdxSigmaProtonBefore);
    fHistoTOFdEdxProtonAfter=new TH2F(Form("Proton_dEdx_TOF_after%s",GetCutNumber().Data()),"dEdx Proton after; p_{T}; TOF signal" ,150,0.03,20,800,0,30000);
    fHistograms->Add(fHistoTOFdEdxProtonAfter);
    fHistoTOFdEdxSigmaProtonAfter=new TH2F(Form("Proton_dEdxSigma_TOF_after %s",GetCutNumber().Data()),"dEdx Sigma Proton after; p_{T}; TOF sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTOFdEdxSigmaProtonAfter);
  }
  

  // dEdx Cuts pion
  fHistodEdxCutsPion=new TH2F(Form("dEdxCuts pion %s",GetCutNumber().Data()),"dEdxCuts pion vs p_{T,e}",11,-0.5,10.5,250,0,50);
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(1,"in");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(2,"TPCpion");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(3,"TPCpion");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(4,"TPCpionhighp");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(5,"TPCkaonlowprej");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(6,"TPCprotonlowprej");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(7,"TPCpionlowprej");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(8,"TOFelectron");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(9,"ITSelectron");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(10,"TRDelectron");
  fHistodEdxCutsPion->GetXaxis()->SetBinLabel(11,"out");
  fHistograms->Add(fHistodEdxCutsPion);
  
  if(fIsQA==kTRUE){
    fHistoTPCdEdxPionBefore=new TH2F(Form("Pion_dEdx_TPC_before %s",GetCutNumber().Data()),"dEdx Pion before; p_{T}; TPC signal" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxPionBefore);
    fHistoTPCdEdxSigmaPionBefore=new TH2F(Form("Pion_dEdxSigma_TPC_before %s",GetCutNumber().Data()),"dEdx Sigma Pion before; p_{T}; TPC sigma pion" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTPCdEdxSigmaPionBefore);
    fHistoTPCdEdxPionAfter=new TH2F(Form("Pion_dEdx_TPC_after%s",GetCutNumber().Data()),"dEdx Pion after; p_{T}; TPC signal" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxPionAfter);
    fHistoTPCdEdxSigmaPionAfter=new TH2F(Form("Pion_dEdxSigma_TPC_after %s",GetCutNumber().Data()),"dEdx Sigma Pion after; p_{T}; TPC sigma pion" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTPCdEdxSigmaPionAfter); 
    
    fHistoTOFdEdxPionBefore=new TH2F(Form("Pion_dEdx_TOF_before %s",GetCutNumber().Data()),"dEdx Pion before; p_{T}; TOF signal" ,150,0.03,20,800,0,30000);
    fHistograms->Add(fHistoTOFdEdxPionBefore);
    fHistoTOFdEdxSigmaPionBefore=new TH2F(Form("Pion_dEdxSigma_TOF_before %s",GetCutNumber().Data()),"dEdx Sigma Pion before; p_{T}; TOF sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTOFdEdxSigmaPionBefore);
    fHistoTOFdEdxPionAfter=new TH2F(Form("Pion_dEdx_TOF_after%s",GetCutNumber().Data()),"dEdx Pion after; p_{T}; TOF signal" ,150,0.03,20,800,0,30000);
    fHistograms->Add(fHistoTOFdEdxPionAfter);
    fHistoTOFdEdxSigmaPionAfter=new TH2F(Form("Pion_dEdxSigma_TOF_after %s",GetCutNumber().Data()),"dEdx Sigma Pion after; p_{T}; TOF sigma proton" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoTOFdEdxSigmaPionAfter);
  }
  
  TH1::AddDirectory(kTRUE);
}


///________________________________________________________________________
TString AliV0CutsStrange::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}


//________________________________________________________________________
Bool_t AliV0CutsStrange::InitPIDResponse(){
  // Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;
  }
  return kFALSE;
}
  

///________________________________________________________________________
Bool_t AliV0CutsStrange::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  switch (cutID) {

    case kv0FinderType:
      if( SetV0Finder(value)) {
        fCuts[kv0FinderType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kclsTPCCut:
      if( SetTPCClusterCut(value)) {
        fCuts[kclsTPCCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;    
    case kProtonPIDcut:
      if( SetProtonPIDCut(value)) {
        fCuts[kProtonPIDcut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;    
    case kPionPIDcut:
      if( SetPionPIDCut(value)) {
        fCuts[kPionPIDcut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;
}


///________________________________________________________________________
Bool_t AliV0CutsStrange::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());
  
  // Initialize Cuts from a given Cut string

  AliInfo(Form("Set V0 Cut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

//   if (fV0ReaderStrangeName.CompareTo("") == 0){
//     fV0ReaderStrangeName = "V0ReaderV1";
//   }
  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  PrintCutsWithValues();

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliV0CutsStrange::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
void AliV0CutsStrange::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

///________________________________________________________________________
void AliV0CutsStrange::PrintCutsWithValues() {
  printf("\nConversion cutnumber \n");
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%d",fCuts[ic]);
  }
  
  if (!fUseCorrectedTPCClsInfo) printf("\t # TPC clusters > %3.2f \n", fMinClsTPC);
  if (fUseCorrectedTPCClsInfo) printf("\t #cluster TPC/ #findable clusters TPC (corrected for radius) > %3.2f\n", fMinClsTPCToF );
  
}



///________________________________________________________________________
Bool_t AliV0CutsStrange::SetV0Finder(Int_t v0FinderType){   // Set Cut
  switch (v0FinderType){
  case 0:  // on fly V0 finder
    cout << "have chosen onfly V0" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 1:  // offline V0 finder
    cout << "have chosen offline V0" << endl;
    fUseOnFlyV0Finder=kFALSE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 2:  // on fly V0 finder with same signs
    cout << "have chosen onfly V0 same sign pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=1;
    break;
  case 3:  // on fly V0 finder with unlike signs and same signs
    cout << "have chosen onfly V0 unlike sign and same signs pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=2;
    break;
  default:
    AliError(Form(" v0FinderType not defined %d",v0FinderType));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliV0CutsStrange::GetTrack(AliVEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
  if(esdEvent) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = esdEvent->GetTrack(label);
    return track;

  } else {
    if(label == -999999) return NULL; // if AOD relabelling goes wrong, immediately return NULL
    AliVTrack * track = 0x0;
    if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderStrangeName.Data()))->AreAODsRelabeled()){
      if(event->GetTrack(label)) track = dynamic_cast<AliVTrack*>(event->GetTrack(label));
      return track;
    }
    else{
      for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
        if(event->GetTrack(ii)) track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
        if(track){
          if(track->GetID() == label) {
            return track;
          }
        }
      }
    }
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}

///________________________________________________________________________
AliESDtrack *AliV0CutsStrange::GetESDTrack(AliESDEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  if(event) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = event->GetTrack(label);
    return track;
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}



///________________________________________________________________________
Bool_t AliV0CutsStrange::PhotonIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent,Bool_t checkForConvertedGamma)
{
  // MonteCarlo Photon Selection

  if(!mcEvent)return kFALSE;

//   if (particle->GetPdgCode() == 22){

//many many cuts to be included
    //     if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
    //       return kFALSE;
    //     if(fEtaCutMin>-0.1){
    //       if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
    //         return kFALSE;
    //     }
    
    //     if(particle->GetMother(0) >-1 && mcEvent->Particle(particle->GetMother(0))->GetPdgCode() == 22){
    //       return kFALSE; // no photon as mothers!
    //     }
    /*
    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    TParticle* ePos = NULL;
    TParticle* eNeg = NULL;
    
    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
        if(daughterIndex<0) continue;
        TParticle *tmpDaughter = mcEvent->Particle(daughterIndex);
        
        if(tmpDaughter->GetUniqueID() == 5){
          if(tmpDaughter->GetPdgCode() == 11){
            eNeg = tmpDaughter;
          } else if(tmpDaughter->GetPdgCode() == -11){
            ePos = tmpDaughter;
          }
        }
      }
    }
    
    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }
    
    if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
      return kFALSE; // no reconstruction below the Pt cut
    }

    if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
      eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
      return kFALSE;

    if(fEtaCutMin > -0.1){
      if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
        (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
        return kFALSE;
    }

    if(ePos->R()>fMaxR){
      return kFALSE; // cuts on distance from collision point
    }

    if(fabs(ePos->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }
    if(fabs(eNeg->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }

    if( ePos->R() <= ((fabs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((fabs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    if( eNeg->R() <= ((fabs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((fabs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
  }*/
  return kFALSE;
}



///________________________________________________________________________
Bool_t AliV0CutsStrange::GetPIDpion(AliVTrack *fCurrentTrack){
  // Pion Identification Cuts for V0 reconstruction
  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error
  
    // use TOF signal if available, for tracks with p > 0.75
  AliPIDResponse::EDetPidStatus statusPosTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fCurrentTrack);
  Bool_t isPIDTOF = kFALSE;
  if (AliPIDResponse::kDetPidOk == statusPosTOF && fCurrentTrack->GetP()>0.75) {
    isPIDTOF = kTRUE;
  }
  
  Int_t cutIndex=0;
  if(fHistodEdxCutsPion)fHistodEdxCutsPion->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigmaPionBefore)fHistoTPCdEdxSigmaPionBefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kPion));
  if(fHistoTPCdEdxPionBefore)fHistoTPCdEdxPionBefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  cutIndex++;
  
  //select protons
  if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDTPCnSigmaPionLow || fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)>fPIDTPCnSigmaPionUp){
    return kFALSE;
  }
  
  if(fHistodEdxCutsPion)fHistodEdxCutsPion->Fill(cutIndex,fCurrentTrack->Pt());
  cutIndex++;
    
  //select protons with TOF
  if(isPIDTOF){
    if(fHistoTOFdEdxSigmaPionBefore)fHistoTOFdEdxSigmaPionBefore->Fill(fCurrentTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kPion));
    if(fHistoTOFdEdxPionBefore)fHistoTOFdEdxPionBefore->Fill(fCurrentTrack->Pt(),fCurrentTrack->GetTOFsignal());
    if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kPion)<fPIDTOFnSigmaPionLow || fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kPion)>fPIDTOFnSigmaPionUp){
      return kFALSE;
    }
  }
  
  if(fHistodEdxCutsProton)fHistodEdxCutsPion->Fill(cutIndex,fCurrentTrack->Pt());
  cutIndex++;
  
  if(fHistodEdxCutsPion)fHistodEdxCutsPion->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigmaPionAfter)fHistoTPCdEdxSigmaPionAfter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kPion));
  if(fHistoTPCdEdxPionAfter)fHistoTPCdEdxPionAfter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  if(isPIDTOF){
    if(fHistoTOFdEdxSigmaPionAfter)fHistoTOFdEdxSigmaPionAfter->Fill(fCurrentTrack->Pt(), fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kPion));
    if(fHistoTOFdEdxPionAfter)fHistoTOFdEdxPionAfter->Fill(fCurrentTrack->Pt(),fCurrentTrack->GetTOFsignal());
  }  
  
  return kTRUE;
}



///________________________________________________________________________
Bool_t AliV0CutsStrange::GetPIDproton(AliVTrack *fCurrentTrack){
  // Proton Identification Cuts for V0 reconstruction
  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error
  
  // use TOF signal if available, for tracks with p > 0.75
  AliPIDResponse::EDetPidStatus statusPosTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fCurrentTrack);
  Bool_t isPIDTOF = kFALSE;
  if (AliPIDResponse::kDetPidOk == statusPosTOF && fCurrentTrack->GetP()>0.75) {
    isPIDTOF = kTRUE;
  }
  
  Int_t cutIndex=0;
  if(fHistodEdxCutsProton)fHistodEdxCutsProton->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigmaProtonBefore)fHistoTPCdEdxSigmaProtonBefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kProton));
  if(fHistoTPCdEdxProtonBefore)fHistoTPCdEdxProtonBefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  cutIndex++;
  
  //select protons with TPC
  if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton)<fPIDTPCnSigmaProtonLow || fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton)>fPIDTPCnSigmaProtonUp){
    return kFALSE;
  }
  
  if(fHistodEdxCutsProton)fHistodEdxCutsProton->Fill(cutIndex,fCurrentTrack->Pt());
  cutIndex++;
    
  //select protons with TOF
  if(isPIDTOF){
    if(fHistoTOFdEdxSigmaProtonBefore)fHistoTOFdEdxSigmaProtonBefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kProton));
    if(fHistoTOFdEdxProtonBefore)fHistoTOFdEdxProtonBefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTOFsignal());
    if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kProton)<fPIDTOFnSigmaProtonLow || fPIDResponse->NumberOfSigmasTOF(fCurrentTrack,AliPID::kProton)>fPIDTOFnSigmaProtonUp){
      return kFALSE;
    }
  }
  
  if(fHistodEdxCutsProton)fHistodEdxCutsProton->Fill(cutIndex,fCurrentTrack->Pt());
  cutIndex++;
  
  if(fHistodEdxCutsProton)fHistodEdxCutsProton->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigmaProtonAfter)fHistoTPCdEdxSigmaProtonAfter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kProton));
  if(fHistoTPCdEdxProtonAfter)fHistoTPCdEdxProtonAfter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  
  if(isPIDTOF){
    if(fHistoTOFdEdxSigmaProtonAfter)fHistoTOFdEdxSigmaProtonAfter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kProton));
    if(fHistoTOFdEdxProtonAfter)fHistoTOFdEdxProtonAfter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTOFsignal());
  }

  return kTRUE;
}


///________________________________________________________________________
Bool_t AliV0CutsStrange::SetTPCClusterCut(Int_t clsTPCCut){   // Set Cut
  switch(clsTPCCut){
  case 0: // 0
    fMinClsTPC= 0.;
    break;
  case 1:  // 60
    fMinClsTPC= 60.;
    break;
  case 2:  // 80
    fMinClsTPC= 80.;
    break;
  case 3:  // 100
    fMinClsTPC= 100.;
    break;
  case 4:  // 95% of findable clusters
    fMinClsTPCToF= 0.95;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 5:  // 0% of findable clusters
    fMinClsTPCToF= 0.0;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 6:  // 70% of findable clusters
    fMinClsTPCToF= 0.7;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 7:  // 0% of findable clusters
    fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=0;
    break;
  case 8:
    fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 9:
    fMinClsTPCToF= 0.6;
    fUseCorrectedTPCClsInfo=1;
    break;
  default:
    AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliV0CutsStrange::SetProtonPIDCut(Int_t pPIDcut){   // Set Cut
  switch(pPIDcut){
    case 0:
      fPIDTPCnSigmaProtonLow=-3;
      fPIDTPCnSigmaProtonUp = 3;
      fPIDTOFnSigmaProtonLow=-3;
      fPIDTOFnSigmaProtonUp = 3;
      break;
    case 1:
      fPIDTPCnSigmaProtonLow=-4;
      fPIDTPCnSigmaProtonUp = 4;
      fPIDTOFnSigmaProtonLow=-4;
      fPIDTOFnSigmaProtonUp = 4;
    default:
    AliError(Form("Warning: ProtonPIDcut not defined %d",pPIDcut));
    return kFALSE;
  }
  return kTRUE;
  
  
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliV0CutsStrange::SetPionPIDCut(Int_t pPIDcut){   // Set Cut
  switch(pPIDcut){
    case 0:
      fPIDTPCnSigmaPionLow=-3;
      fPIDTPCnSigmaPionUp = 3;
      fPIDTOFnSigmaPionLow=-3;
      fPIDTOFnSigmaPionUp = 3;
      break;
    case 1:
      fPIDTPCnSigmaPionLow=-4;
      fPIDTPCnSigmaPionUp = 4;
      fPIDTOFnSigmaPionLow=-4;
      fPIDTOFnSigmaPionUp = 4;
    default:
    AliError(Form("Warning: Pion PIDcut not defined %d",pPIDcut));
    return kFALSE;
  }
  return kTRUE;
}
