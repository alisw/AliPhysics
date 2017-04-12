#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "TParticle.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"


#include "AliESDv0.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronV0Cuts.h"

#include "TSystem.h"
#include <string>
#include "AliAnalysisTaskSimpleTreeMaker.h"

/*************** Tree Maker Class **********************
*                                                      *
* Created: 05.10.2016                                  *
* Authors: Aaron Capon      (aaron.capon@cern.ch)      *
*          Sebastian Lehner (sebastian.lehner@cern.ch) *
*                                                      *
*******************************************************/


ClassImp(AliAnalysisTaskSimpleTreeMaker)

Int_t eventNum = 0;

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker():
  AliAnalysisTaskSE(),
  fTree(0),
  fStream(0),
  fQAhist(0),
  esdEvent(0), 
  mcEvent(0), 
  fESDtrackCuts(0),
  fPIDResponse(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(80), 
  fPtMin(0.2),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigITSMin(-3.),
  fESigITSMax(3.),
  fESigTPCMin(-5.),
  fESigTPCMax(5.),
  fESigTOFMin(-3.),
  fESigTOFMax(3.),
  fPSigTPCMin(-99.),
  fPSigTPCMax(-3.),
  fPIDcutITS(kFALSE),
  fPIDcutTOF(kFALSE),
  fPionPIDcutTPC(kFALSE),
  fIsIonColl(kFALSE),
  fIsMC(kTRUE),
  fHasSDD(kTRUE),
  fIsGRIDanalysis(kTRUE),
  fIsV0tree(kFALSE),
  fArmPlot(0),
  fGridPID(-1)
{

}

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker(const char *name) :
  AliAnalysisTaskSE(name),
  fTree(0),
  fStream(0),
  fQAhist(0),
  esdEvent(0), 
  mcEvent(0), 
  fESDtrackCuts(0),
  fPIDResponse(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(80), 
  fPtMin(0.2),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigITSMin(-3.),
  fESigITSMax(3.),
  fESigTPCMin(-4.),
  fESigTPCMax(4.),
  fESigTOFMin(-3.),
  fESigTOFMax(3.),
  fPSigTPCMin(-99.),
  fPSigTPCMax(-3.),
  fPIDcutITS(kFALSE),
  fPIDcutTOF(kFALSE),
  fPionPIDcutTPC(kFALSE),
  fIsIonColl(kFALSE),
  fIsMC(kTRUE),
  fHasSDD(kTRUE),
  fIsGRIDanalysis(kTRUE),
  fIsV0tree(kFALSE),
  fArmPlot(0),
  fGridPID(-1)
{
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class()); //will be connected to fTree
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH2F::Class());
}

//________________________________________________________________________

//~ AliAnalysisTaskSimpleTreeMaker::~AliAnalysisTaskSimpleTreeMaker() {

  //~ // Destructor

  //~ // ... not implemented

//~ }


//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserCreateOutputObjects() {
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  inputHandler->SetNeedField();
     
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse){
    AliFatal("This task needs the PID response attached to the inputHandler");
	  return;
  } 

  fStream = new TTreeStream("tracks");
  fTree   = (TTree*)fStream->GetTree();

  //Get grid PID
  if( fIsGRIDanalysis ){
    const char* gridIDchar = gSystem->Getenv("ALIEN_PROC_ID");
    std::string str(gridIDchar);
    SetGridPID(str);
  }
  else{ 
    fGridPID = -1;
  }
  if(fIsV0tree){
    fArmPlot = new TH2F("ArmPlot", "Armenteros Plot", 100, -1, 1, 100, 0, 0.4);
    fArmPlot->GetXaxis()->SetTitle("#alpha = (p^{+}-p^{-})/(p^{+}+p^{-})");
    fArmPlot->GetYaxis()->SetTitle("p_{T}");
  }  
  fQAhist = new TH1F("h1", "Event and track QA", 6, 0, 1);
  PostData(1, fTree);
  PostData(2, fQAhist);
  PostData(3, fArmPlot);
}

//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  esdEvent = (AliESDEvent*)InputEvent();
  if(!esdEvent) {
    AliError("No AliESDEvent");
    return;
  } 
  fQAhist->Fill("Events_ESDcheck",1);

  // Process also MC truth  
  mcEvent = MCEvent();
  if( fIsMC ){
    if( !mcEvent ){
      AliError("Could not retrieve MC event");
      return;
    }
    fQAhist->Fill("Events_MCcheck",1);
  }
  eventNum += 1;

  // check event cuts
  if( IsEventAccepted(esdEvent) == 0){ 
    return;
  }
  fQAhist->Fill("Events_accepted",1);
  
  // PID Response task active?
  fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();

  if (!fPIDResponse){ AliFatal("This task needs the PID response attached to the inputHandler"); }
 
  AliESDVertex* vertex = (AliESDVertex*)esdEvent->GetPrimaryVertex();

  Double_t primaryVertex[3];
  primaryVertex[0] = vertex->GetX();
  primaryVertex[1] = vertex->GetY();
  primaryVertex[2] = vertex->GetZ();

  Double_t impactParameter = -999;
  if(fIsMC){
    AliGenHijingEventHeader* hHijing = 0;
    AliGenEventHeader*  mcGenH  = mcEvent->GenEventHeader();

    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
    //Option 1: Just Hijing
      hHijing = (AliGenHijingEventHeader*)mcGenH;
    } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
      //Option 2: cocktail involving Hijing
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("hijing_0"));
    }   
    if(hHijing){
      impactParameter = hHijing->ImpactParameter();
    }
  }

  Int_t eventTracks = esdEvent->GetNumberOfTracks();
  Int_t runNumber = esdEvent->GetRunNumber();
  Int_t numV0s = esdEvent->GetNumberOfV0s();

  //Loop over tracks for event
  if( !fIsV0tree ){
    for (Int_t iTrack = 0; iTrack < eventTracks; iTrack++){ 
      AliESDtrack* track = esdEvent->GetTrack(iTrack);
      if (!track) {
        AliError(Form("Could not receive track %d", iTrack));
        continue;
      }

      fQAhist->Fill("Tracks_all",1);
      //Apply global track filter
      if(!fESDtrackCuts->AcceptTrack(track)){ continue; }

      //Apply momentum and eta cuts
      Double_t pt   = track->Pt();
      if( pt < fPtMin || pt > fPtMax ){ continue; }
      Double_t eta  = track->Eta();
      if( eta < fEtaMin || eta > fEtaMax ){ continue; } 
      fQAhist->Fill("Tracks_KineCuts", 1);

      Double_t phi  = track->Phi();

      //Get electron nSigma in TPC for cut (inclusive cut)
      Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
      if( EnSigmaTPC > fESigTPCMax || EnSigmaTPC < fESigTPCMin) { continue; }
      
      Double_t EnSigmaITS = -999;
      if(fHasSDD){
        //Get rest of electron nSigma values and apply cuts if requested (inclusive cuts)
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron);
        if(fPIDcutITS){
          if(EnSigmaITS < fESigITSMin || EnSigmaITS > fESigITSMax){ continue; }
        }
      }

      Double_t EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron);
      if(fPIDcutTOF){
        if(EnSigmaTOF < fESigTOFMin || EnSigmaTOF > fESigTOFMax){ continue; }
      }

      //Get pion nSigma for TPC and apply cut if requested (exclusive cut)
      Double_t PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion);
      if(fPionPIDcutTPC){
        if(PnSigmaTPC > fPSigTPCMin && PnSigmaTPC < fPSigTPCMax){ continue; }
      }
      
      fQAhist->Fill("Tracks_PIDcuts",1); 

      //Get rest of nSigma values for pion and kaon
      Double_t PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion);
      Double_t PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion);
      
      Double_t KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon);
      Double_t KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon);
      Double_t KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon);
      

      //Get ITS and TPC signals
      Double_t ITSsignal = track->GetITSsignal();
      Double_t TPCsignal = track->GetTPCsignal();
      Double_t TOFsignal = track->GetTOFsignal();
      
      Double_t nTPCclusters = track->GetTPCNcls(); 
      Double_t nTPCcrossed = track->GetTPCClusterInfo(2,1);
      Double_t TPCcrossOverFind = nTPCcrossed/track->GetTPCNclsF();
      Double_t nTPCshared = track->GetTPCnclsS();
      Double_t chi2TPC = track->GetTPCchi2();

      //DCA values
      Float_t DCAxy = 0.;
      Float_t DCAz = 0.;
      track->GetImpactParameters( &DCAxy, &DCAz);
      
      Int_t nITS = 0;
      Double_t fITS_shared = 0;
      Double_t chi2ITS = 0;
      if(fHasSDD){
        //ITS clusters and shared clusters
        nITS = track->GetNumberOfITSClusters();
        fITS_shared = 0.;
        for(Int_t d = 0; d < 6; d++){
                fITS_shared += static_cast<Double_t>(track->HasSharedPointOnITSLayer(d));
          }
        fITS_shared /= nITS;
     
        chi2ITS = track->GetITSchi2();
      }

      Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDtrackCuts->kVertexTracks | fESDtrackCuts->kVertexSPD;

      const AliESDVertex* vertex = 0;

      if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTracks){
        vertex = track->GetESDEvent()->GetPrimaryVertexTracks();}

      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexSPD){
        vertex = track->GetESDEvent()->GetPrimaryVertexSPD();}

      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTPC){
        vertex = track->GetESDEvent()->GetPrimaryVertexTPC();}

      Double_t goldenChi2 = track->GetChi2TPCConstrainedVsGlobal(vertex);

      
      Int_t charge = -10;
      charge = track->Charge();
      Int_t label = -999;
      label = track->GetLabel();

      //Declare MC variables
      Double_t mcEta   = -99;
      Double_t mcPhi   = -99;
      Double_t mcPt    = -99;
      Int_t iPdg       = 0;
      Int_t iPdgMother = 0;
      if(fIsMC){
        AliMCParticle* mcTrack = (AliMCParticle*) mcEvent->GetTrack(TMath::Abs(label));
        TParticle* mcParticle = (TParticle*)mcTrack->Particle();
        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt = mcTrack->Pt();

        iPdg = mcTrack->PdgCode();

        Int_t gMotherIndex = mcTrack->GetMother();
        AliMCParticle* motherTrack = (AliMCParticle*)(mcEvent->GetTrack(gMotherIndex));
        iPdgMother = motherTrack->PdgCode();
      }

      //Stream values into tree
      if( fIsMC ){
        (*fStream)    << "tracks" <<
        "pt="         << pt << 
        "eta="        << eta << 
        "phi="        << phi << 
        "EsigITS="    << EnSigmaITS <<
        "EsigTPC="    << EnSigmaTPC <<
        "EsigTOF="    << EnSigmaTOF <<
        "PsigITS="    << PnSigmaITS <<
        "PsigTPC="    << PnSigmaTPC <<
        "PsigTOF="    << PnSigmaTOF <<
        "KsigITS="    << KnSigmaITS <<
        "KsigTPC="    << KnSigmaTPC <<
        "KsigTOF="    << KnSigmaTOF <<

        "charge="     << charge <<
        "ITSsignal="  << ITSsignal <<
        "TPCsignal="  << TPCsignal << 
        "TOFsignal="  << TOFsignal <<
        "vertexX="    << primaryVertex[0] <<
        "vertexY="    << primaryVertex[1] <<
        "vertexZ="    << primaryVertex[2] <<

        "nTPCclusters=" << nTPCclusters <<
        "nTPCcrossed="<< nTPCcrossed <<
        "TPCcrossFind="<< TPCcrossOverFind <<
        "nTPCshared=" << nTPCshared <<
        "chi2TPC="    << chi2TPC <<

        "nITS="       << nITS <<
        "fITSshared=" << fITS_shared << 
        "chi2ITS="    << chi2ITS <<

        "DCAxy="      << DCAxy <<
        "DCAz="       << DCAz <<
        "goldenChi2=" << goldenChi2 <<

        "mcEta="      << mcEta <<
        "mcPhi="      << mcPhi <<
        "mcPt="       << mcPt <<
        "pdg="        << iPdg <<
        "pdgMother="  << iPdgMother <<
        "runNumber="  << runNumber << 
        "eventNum="   << eventNum <<
        "gridPID="    << fGridPID <<
        "\n";
      }
      else{
        (*fStream)    << "tracks" <<
        "pt="         << pt << 
        "eta="        << eta << 
        "phi="        << phi << 
        "EsigITS="    << EnSigmaITS <<
        "EsigTPC="    << EnSigmaTPC <<
        "EsigTOF="    << EnSigmaTOF <<
        "PsigITS="    << PnSigmaITS <<
        "PsigTPC="    << PnSigmaTPC <<
        "PsigTOF="    << PnSigmaTOF <<
        "KsigITS="    << KnSigmaITS <<
        "KsigTPC="    << KnSigmaTPC <<
        "KsigTOF="    << KnSigmaTOF <<

        "charge="     << charge <<
        "ITSsignal="  << ITSsignal <<
        "TPCsignal="  << TPCsignal << 
        "TOFsignal="  << TOFsignal <<
        "vertexX="    << primaryVertex[0] <<
        "vertexY="    << primaryVertex[1] <<
        "vertexZ="    << primaryVertex[2] <<

        "nTPCclusters=" << nTPCclusters <<
        "nTPCcrossed=" << nTPCcrossed <<
        "TPCcrossFind="<< TPCcrossOverFind <<
        "nTPCshared=" << nTPCshared <<
        "chi2TPC="    << chi2TPC <<

        "nITS="       << nITS <<
        "fITSshared=" << fITS_shared << 
        "chi2ITS="    << chi2ITS <<

        "DCAxy="      << DCAxy <<
        "DCAz="       << DCAz <<
        "goldenChi2=" << goldenChi2 <<
        "runNumber="  << runNumber << 
        "eventNum="   << eventNum <<
        "gridPID="    << fGridPID <<
        "\n";
      }
    } //End loop over tracks
  }
  else{
    for(Int_t iV0 = 0; iV0 < numV0s; iV0++){
     
      AliESDv0* V0vertex = esdEvent->GetV0(iV0);
        
      if(!V0vertex){
          AliError(Form("Could not receive V0 track %d", iV0));
          continue;
      }

      fQAhist->Fill("Tracks_all",1);

      //Get V0 daughter tracks
      AliESDtrack* negTrack = esdEvent->GetTrack(V0vertex->GetIndex(0));
      AliESDtrack* posTrack = esdEvent->GetTrack(V0vertex->GetIndex(1));
      if(!negTrack | !posTrack){
        Printf("Daughter track of v0 not found: %p - %p \n",negTrack, posTrack);
        continue;
      }
      //Check for like-sign V0 candidates
      if(negTrack->Charge() == posTrack->Charge()){ continue; }

      //Apply kinematic and PID cuts to both legs 
      if(isV0daughterAccepted(negTrack) != kTRUE){ continue; }
      if(isV0daughterAccepted(posTrack) != kTRUE){ continue; }

      Double_t pointingAngle = V0vertex->GetV0CosineOfPointingAngle();
      Double_t daughtersDCA = V0vertex->GetDcaV0Daughters();
      Double_t decayLength = V0vertex->GetRr();
      Double_t v0mass = V0vertex->M();
      //Super loose cuts on V0 topological qualities(stored in Tree to be cut on later)
      if( pointingAngle < 0.8 | daughtersDCA < 0.05 | decayLength < 0.01 ){ continue; }

      Double_t ptArm = V0vertex->PtArmV0();
      Double_t alpha = V0vertex->AlphaV0();
      fArmPlot->Fill(alpha, ptArm);
      
      fQAhist->Fill("Arm. cuts",1);

      //Get positive particle obsevables
      Double_t pt = posTrack->Pt();
      Double_t eta = posTrack->Eta();
      Double_t phi = posTrack->Phi();

      Double_t EnSigmaITS = -999;
      Double_t PnSigmaITS = -999;
      Double_t KnSigmaITS = -999;
      if(fHasSDD){
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kElectron);
        PnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kPion);
        KnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kKaon);
      }
      Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kElectron);
      Double_t EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kElectron);
      Double_t PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kPion);
      Double_t PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kPion);
      Double_t KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kKaon);
      Double_t KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kKaon);
      //DCA values
      Float_t ImpParamXY = 0.;
      Float_t ImpParamZ = 0.;
      posTrack->GetImpactParameters( &ImpParamXY, &ImpParamZ);
      
      Int_t nITS = 0;
      Double_t fITS_shared = 0;
      if(fHasSDD){
      //ITS clusters and shared clusters
        nITS = posTrack->GetNumberOfITSClusters();
        for(Int_t d = 0; d < 6; d++){
          fITS_shared += static_cast<Double_t>(posTrack->HasSharedPointOnITSLayer(d));
        }
        fITS_shared /= nITS;
      }
      Int_t daughtCharge = posTrack->Charge();
      
      //Declare MC variables
      Double_t mcEta   = -99;
      Double_t mcPhi   = -99;
      Double_t mcPt    = -99;
      Int_t iPdg       = 0;
      Int_t iPdgMother = 0;
      Int_t label = -999;

      if(fIsMC){
        label = posTrack->GetLabel();
        AliMCParticle* mcTrack = (AliMCParticle*) mcEvent->GetTrack(TMath::Abs(label));
        TParticle* mcParticle = (TParticle*)mcTrack->Particle();
        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt = mcTrack->Pt();

        iPdg = mcTrack->PdgCode();

        Int_t gMotherIndex = mcTrack->GetMother();
        AliMCParticle* motherTrack = (AliMCParticle*)(mcEvent->GetTrack(gMotherIndex));
        iPdgMother = motherTrack->PdgCode();
        (*fStream)    << "tracks" <<
          //Positive particle obsevables
          "pt="         << pt << 
          "eta="        << eta << 
          "phi="        << phi << 
          "EsigITS="    << EnSigmaITS <<
          "EsigTPC="    << EnSigmaTPC <<
          "EsigTOF="    << EnSigmaTOF <<
          "PsigITS="    << PnSigmaITS <<
          "PsigTPC="    << PnSigmaTPC <<
          "PsigTOF="    << PnSigmaTOF <<
          "KsigITS="    << KnSigmaITS <<
          "KsigTPC="    << KnSigmaTPC <<
          "KsigTOF="    << KnSigmaTOF <<
          "nITS="       << nITS <<
          "fITSshared=" << fITS_shared << 
          "impParamXY=" << ImpParamXY <<
          "impParamZ="  << ImpParamZ <<
          "charge="     << daughtCharge <<
          "DCA="        << daughtersDCA <<
          "mcEta="      << mcEta <<
          "mcPhi="      << mcPhi <<
          "mcPt="       << mcPt <<
          "pdg="        << iPdg <<
          "pdgMother="  << iPdgMother <<
          //V0 particle observables 
          "v0effMass="  << v0mass <<
          "pointing="   << pointingAngle << 
          "Rlength="    << decayLength <<
          "ptArm="      << ptArm <<
          "alpha="      << alpha <<
          "\n";
      }
      else{
        (*fStream)    << "tracks" <<
          //Positive particle obsevables
          "pt="         << pt << 
          "eta="        << eta << 
          "phi="        << phi << 
          "EsigITS="    << EnSigmaITS <<
          "EsigTPC="    << EnSigmaTPC <<
          "EsigTOF="    << EnSigmaTOF <<
          "PsigITS="    << PnSigmaITS <<
          "PsigTPC="    << PnSigmaTPC <<
          "PsigTOF="    << PnSigmaTOF <<
          "KsigITS="    << KnSigmaITS <<
          "KsigTPC="    << KnSigmaTPC <<
          "KsigTOF="    << KnSigmaTOF <<
          "nITS="       << nITS <<
          "nITSshared=" << fITS_shared << 
          "impParamXY=" << ImpParamXY <<
          "impParamZ="  << ImpParamZ <<
          "charge="     << daughtCharge <<
          "DCA="        << daughtersDCA <<
          //V0 particle observables 
          "v0effMass="  << v0mass <<
          "pointing="   << pointingAngle << 
          "Rlength="    << decayLength <<
          "ptArm="      << ptArm <<
          "alpha="      << alpha <<
          "\n";
      }
      //Get negative particle obsevables
      pt = negTrack->Pt();
      eta = negTrack->Eta();
      phi = negTrack->Phi();

      if(fHasSDD){
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kElectron);
        PnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kPion);
        KnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kKaon);
      }
      EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kElectron);
      EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kElectron);
      PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kPion);
      PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kPion);
      KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kKaon);
      KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kKaon);
      //DCA values
      ImpParamXY = 0.;
      ImpParamZ = 0.;
      negTrack->GetImpactParameters( &ImpParamXY, &ImpParamZ);
      
      if(fHasSDD){
        //ITS clusters andared clusters
        nITS = negTrack->GetNumberOfITSClusters();
        fITS_shared = 0.;
        for(Int_t d = 0; d < 6; d++){
                fITS_shared += static_cast<Double_t>(negTrack->HasSharedPointOnITSLayer(d));
          }
        fITS_shared /= nITS;
      }
      daughtCharge = negTrack->Charge(); 
      //Write negative observales to tree (v0 information written twice. Filter by looking at only pos or neg charge)
      if(fIsMC){
        label = negTrack->GetLabel();
        AliMCParticle* mcTrack = (AliMCParticle*) mcEvent->GetTrack(TMath::Abs(label));
        TParticle* mcParticle = (TParticle*)mcTrack->Particle();
        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt = mcTrack->Pt();

        iPdg = mcTrack->PdgCode();

        Int_t gMotherIndex = mcTrack->GetMother();
        AliMCParticle* motherTrack = (AliMCParticle*)(mcEvent->GetTrack(gMotherIndex));
        iPdgMother = motherTrack->PdgCode();
        (*fStream)    << "tracks" <<
          //Positive particle obsevables
          "pt="         << pt << 
          "eta="        << eta << 
          "phi="        << phi << 
          "EsigITS="    << EnSigmaITS <<
          "EsigTPC="    << EnSigmaTPC <<
          "EsigTOF="    << EnSigmaTOF <<
          "PsigITS="    << PnSigmaITS <<
          "PsigTPC="    << PnSigmaTPC <<
          "PsigTOF="    << PnSigmaTOF <<
          "KsigITS="    << KnSigmaITS <<
          "KsigTPC="    << KnSigmaTPC <<
          "KsigTOF="    << KnSigmaTOF <<
          "nITS="       << nITS <<
          "nITSshared=" << fITS_shared << 
          "impParamXY=" << ImpParamXY <<
          "impParamZ="  << ImpParamZ <<
          "charge="     << daughtCharge <<
          "DCA="        << daughtersDCA <<
          "mcEta="      << mcEta <<
          "mcPhi="      << mcPhi <<
          "mcPt="       << mcPt <<
          "pdg="        << iPdg <<
          "pdgMother="  << iPdgMother <<
          //V0 particle observables 
          "v0effMass="  << v0mass <<
          "pointing="   << pointingAngle << 
          "Rlength="    << decayLength <<
          "ptArm="      << ptArm <<
          "alpha="      << alpha <<
          "\n";
      }
      else{
        (*fStream)    << "tracks" <<
          //Positive particle obsevables
          "pt="         << pt << 
          "eta="        << eta << 
          "phi="        << phi << 
          "EsigITS="    << EnSigmaITS <<
          "EsigTPC="    << EnSigmaTPC <<
          "EsigTOF="    << EnSigmaTOF <<
          "PsigITS="    << PnSigmaITS <<
          "PsigTPC="    << PnSigmaTPC <<
          "PsigTOF="    << PnSigmaTOF <<
          "KsigITS="    << KnSigmaITS <<
          "KsigTPC="    << KnSigmaTPC <<
          "KsigTOF="    << KnSigmaTOF <<
          "nITS="       << nITS <<
          "nITSshared=" << fITS_shared << 
          "impParamXY=" << ImpParamXY <<
          "impParamZ="  << ImpParamZ <<
          "charge="     << daughtCharge <<
          "DCA="        << daughtersDCA <<
          //V0 particle observables 
          "v0effMass="  << v0mass <<
          "pointing="   << pointingAngle << 
          "Rlength="    << decayLength <<
          "ptArm="      << ptArm <<
          "alpha="      << alpha <<
          "\n";
      }

     }//End loop over v0's for this event
  }

}

//~ //________________________________________________________________________

void  AliAnalysisTaskSimpleTreeMaker::FinishTaskOutput(){
  // Finish task output

  // not implemented ...

}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::Terminate(Option_t *) {
  // Draw result to the screen

  // Called once at the end of the query

  // not implemented ...

}
//~ 


//________________________________________________________________________

Int_t AliAnalysisTaskSimpleTreeMaker::IsEventAccepted(AliESDEvent *event){


  if (TMath::Abs(event->GetVertex()->GetZ()) < 10  &&  event->GetPrimaryVertexSPD() ){
    if (event->GetPrimaryVertexSPD()->GetNContributors() >0){ return 1; }
    else{ 
      return 0;
    }
  }
  return 0;
}

Bool_t AliAnalysisTaskSimpleTreeMaker::isV0daughterAccepted(AliESDtrack* track){
    
  Bool_t answer = kFALSE;
  //Kinematic cuts
  Double_t pt   = track->Pt();
  if( pt < fPtMin || pt > fPtMax ){ return answer; }
  Double_t eta  = track->Eta();
  if( eta < fEtaMin || eta > fEtaMax ){ return answer; } 
  
  fESDtrackCuts->AcceptTrack(track);
  fQAhist->Fill("Tracks_KineCuts", 1);
  //PID cuts
  Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
  if( EnSigmaTPC > fESigTPCMax || EnSigmaTPC < fESigTPCMin) { return answer; }
  fQAhist->Fill("Tracks_PIDcuts",1); 

  answer = kTRUE;
  return answer;
}


