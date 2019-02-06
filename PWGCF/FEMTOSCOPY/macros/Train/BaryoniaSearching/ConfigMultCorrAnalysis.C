#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoBaryoniaAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoTrioMinvFctn.h"
#include "AliFemtoTrioCut.h"
#include "AliFemtoTrio.h"
#include "AliESDtrack.h"
#endif

enum EParticle { kPionPlus, kPionMinus, kKaonPlus, kKaonMinus, kProt, kAntiProton, kLambda, kAntiLambda, nParticles };
enum ESys  { kAPL, kPAL, nSys };

const char *sysNames[nSys]      = { "APL", "PAL"};
const bool runSys[nSys]         = {   1  ,   0   };

bool ppCollisions;

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis);

AliFemtoMultCorrAnalysis* GetAnalysis();
AliFemtoBasicEventCut* GetEventCut();
AliFemtoESDTrackCut* GetTrackCut(EParticle particle);
AliFemtoV0TrackCut* GetV0TrackCut(EParticle particle);

void GetParticlesForSystem(ESys system,EParticle &mother0,EParticle &mother1,EParticle &daughter0,EParticle &daughter1,EParticle &daughter2);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, int year=2015, bool ppAnalysis=false)
{
  ppCollisions = ppAnalysis;
  
  // create analysis managers
  AliFemtoManager* Manager = new AliFemtoManager();
  AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
  
  // add event reader
  if(ppAnalysis){
    AliFemtoEventReaderAODChain* ReaderPP = GetReaderPP(mcAnalysis);
    Manager->SetEventReader(ReaderPP);
  }
  else if(year==2015){
    AliFemtoEventReaderAODMultSelection* Reader2015 = GetReader2015(mcAnalysis);
    Manager->SetEventReader(Reader2015);
  }
  else if(year==2011){
    AliFemtoEventReaderAODChain* Reader2011 = GetReader2011(mcAnalysis);
    Manager->SetEventReader(Reader2011);
  }
  
  // declare necessary objects
  AliFemtoMultCorrAnalysis  *multCorrAnalysis[nSys];
  AliFemtoBasicEventCut     *eventCut[nSys];
  
  // setup analysis
  int anIter = 0;
  
  EParticle motherPart[2];
  EParticle daughterPart[3];
  
  for (int iSys=0; iSys<nSys; iSys++){
    if (!runSys[iSys]) continue;
    
    // get particle cuts
    GetParticlesForSystem((ESys)iSys,motherPart[0],motherPart[1],daughterPart[0],daughterPart[1],daughterPart[2]);
    
    AliFemtoESDTrackCut *daughterCut0  = GetTrackCut(daughterPart[0]);
    AliFemtoESDTrackCut *daughterCut1  = GetTrackCut(daughterPart[1]);
    AliFemtoESDTrackCut *daughterCut2  = GetTrackCut(daughterPart[2]);
    
    AliFemtoESDTrackCut *motherCut0  = GetTrackCut(motherPart[0]);
    AliFemtoV0TrackCut  *motherCut1  = GetV0TrackCut(motherPart[1]);
    
    // create new analysis
    multCorrAnalysis[anIter] = GetAnalysis();
        
    // get event cut
    eventCut[anIter] = GetEventCut();
        
    // setup anallysis cuts
    multCorrAnalysis[anIter]->SetEventCut(eventCut[anIter]);
    
    multCorrAnalysis[anIter]->SetDaughterCut(daughterCut0,0);
    multCorrAnalysis[anIter]->SetDaughterCut(daughterCut1,1);
    multCorrAnalysis[anIter]->SetDaughterCut(daughterCut2,2);
    
    multCorrAnalysis[anIter]->SetMotherCut(motherCut0,0);
    multCorrAnalysis[anIter]->SetMotherCut(motherCut1,1);
    
    // add analysis to the manager
    Manager->AddAnalysis(multCorrAnalysis[anIter]);
    anIter++;
  }
  
  return Manager;
}

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis)
{
  AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterBit(7);
  Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetReadMC(mcAnalysis);
  
  return Reader;
}

AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis)
{
  AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(7);
  Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetReadMC(mcAnalysis);
  
  return Reader;
}

AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis)
{
  AliFemtoEventReaderAODChain* Reader = new AliFemtoEventReaderAODChain();
  Reader->SetFilterBit(128);
  Reader->SetReadV0(true);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kReference);
  Reader->SetMinPlpContribSPD(3);
  Reader->SetIsPileUpEvent(true);
  Reader->SetReadMC(mcAnalysis);
  
  return Reader;
}

AliFemtoMultCorrAnalysis* GetAnalysis()
{
  AliFemtoMultCorrAnalysis *analysis = new AliFemtoMultCorrAnalysis();
  // here one can put some additional analysis settings in the future
  
  return analysis;
}

AliFemtoBasicEventCut* GetEventCut()
{
  AliFemtoBasicEventCut *eventCut = new AliFemtoBasicEventCut();
  eventCut->SetEventMult(0,10000);
  eventCut->SetVertZPos(-10,10);
  eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
  return eventCut;
}

AliFemtoESDTrackCut* GetTrackCut(EParticle particle)
{
  AliFemtoESDTrackCut *particleCut = new AliFemtoESDTrackCut();
  
  if(particle == kKaonPlus || particle==kKaonMinus){
    particleCut->SetMostProbableKaon();
    particleCut->SetMass(0.493677);
  }
  else if(particle == kPionPlus || particle==kPionMinus){
    particleCut->SetMostProbablePion();
    particleCut->SetMass(0.13956995);
  }
  else if(particle == kProt || particle == kAntiProton){
    particleCut->SetMostProbableProton();
    particleCut->SetMass(0.938272013);
  }
  
  if(particle == kKaonPlus || particle == kPionPlus || particle == kProt){ particleCut->SetCharge( 1.0); }
  else                                                                     { particleCut->SetCharge(-1.0); }
  
  if(particle == kProt || particle == kAntiProton)  particleCut->SetPt(0.7, 5.0);
  else                                                particleCut->SetPt(0.14,1.5);
  particleCut->SetEta(-0.8, 0.8);
  particleCut->SetMaxImpactZ( ppCollisions ? 9999999 : 0.25);
  particleCut->SetMaxImpactXY(ppCollisions ? 9999999 : 0.20);
  particleCut->SetStatus(AliESDtrack::kTPCin);
  particleCut->SetminTPCncls(80);
  particleCut->SetRemoveKinks(kTRUE);
  particleCut->SetLabel(kFALSE);
  particleCut->SetMaxTPCChiNdof(4.0);
  
  //  particleCut->SetNsigma(3.0);
  //  particleCut->SetNsigmaTPCTOF(kTRUE);
  return particleCut;
}

AliFemtoV0TrackCut* GetV0TrackCut(EParticle particle)
{
  if(particle != kLambda && particle != kAntiLambda) return 0;
  double LambdaMass = 1.115683;
  
  AliFemtoV0TrackCut *particleCut = new AliFemtoV0TrackCut();
  particleCut->SetMass(LambdaMass);
  particleCut->SetEta(0.8);
  particleCut->SetPt(0.5, 5.0);
  particleCut->SetEtaDaughters(0.8);
  particleCut->SetTPCnclsDaughters(80);
  particleCut->SetNdofDaughters(4.0);
  particleCut->SetStatusDaughters(AliESDtrack::kTPCrefit);
  particleCut->SetOnFlyStatus(kFALSE);
  particleCut->SetMaxDcaV0Daughters(0.4);
  particleCut->SetMaxDcaV0(0.5);
  particleCut->SetMaxV0DecayLength(60.0);
  particleCut->SetMaxCosPointingAngle(0.9993);
  
  if(particle == kLambda)
  {
    particleCut->SetParticleType(0);
    particleCut->SetPtPosDaughter(0.5, 4.0);
    particleCut->SetPtNegDaughter(0.16, 4.0);
    particleCut->SetMinDaughtersToPrimVertex(0.1, 0.3);
    particleCut->SetInvariantMassLambda(LambdaMass-0.0038, LambdaMass+0.0043);
  }
  else if(particle == kAntiLambda)
  {
    particleCut->SetParticleType(1);
    particleCut->SetPtPosDaughter(0.16, 4.0);
    particleCut->SetPtNegDaughter(0.3, 4.0);
    particleCut->SetMinDaughtersToPrimVertex(0.3, 0.1);
    particleCut->SetInvariantMassLambda(LambdaMass-0.0036, LambdaMass+0.0041);
  }
  return particleCut;
}

void GetParticlesForSystem(ESys system,EParticle &mother0,EParticle &mother1,EParticle &daughter0,EParticle &daughter1,EParticle &daughter2)
{
  if(system == kAPL){
    mother0 = kAntiProton;
    mother1 = kLambda;
    daughter0  = kKaonMinus;
    daughter1  = kPionPlus;
    daughter2  = kPionMinus;
  }
  if(system == kPAL){
    mother0 = kProt;
    mother1 = kAntiLambda;
    daughter0  = kKaonPlus;
    daughter1  = kPionPlus;
    daughter2  = kPionMinus;
  }
}










