#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoEventAnalysis.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoSpatialSeparationFunction.h"
#include "AliESDtrack.h"
#endif

enum EPart { kELambda , kEAntiLambda , kEProton , kEAntiProton};

enum ESys { kLL , kALAL , kLAL , kPL , kAPL , kPAL , kAPAL , kPP , kPAP , kAPAP, nSys };
const char *sysNames[nSys] = { "LL", "ALAL", "LAL", "PL", "APL", "PAL", "APAL","PP","PAP","APAP" };
int runSys[nSys] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

const int nMult = 10;
int runMult[nMult] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
int multBins[nMult+1] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoV0TrackCut* GetV0TrackCut(EPart particle);
AliFemtoESDTrackCut* GetESDTrackCut(EPart particle);
void GetParticlesForSystem(ESys system, EPart &firstParticle, EPart &secondParticle);
bool AreIdentical(ESys system);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, int year=2015)
{
  // create analysis managers
  AliFemtoManager* Manager=new AliFemtoManager();
  
  // add event reader
  if(year==2015)
  {
    AliFemtoEventReaderAODMultSelection* Reader2015 = GetReader2015(mcAnalysis);
    Manager->SetEventReader(Reader2015);
  }
  else if(year==2011)
  {
    AliFemtoEventReaderAODChain* Reader2011 = GetReader2011(mcAnalysis);
    Manager->SetEventReader(Reader2011);
  }
  
  AliFemtoEventAnalysis *femtoAnalysis[nSys*nMult];
  AliFemtoSpatialSeparationFunction *separationFunction[nSys*nMult];
  
  AliFemtoBasicEventCut *eventCut = GetEventCut();
  
  int anIter = 0;
  for (int imult=0; imult<nMult; imult++)
  {
    if (!runMult[imult]) continue;
    
    for(int iSys=0;iSys<nSys;iSys++)
    {
      if (!runSys[iSys]) continue;
      
      anIter = imult * nSys + iSys;
      
      femtoAnalysis[iSys] = new AliFemtoEventAnalysis(multBins[imult], multBins[imult+1]);
      femtoAnalysis[iSys]->SetIdenticalParticles(AreIdentical(iSys));
      separationFunction[iSys] = new AliFemtoSpatialSeparationFunction(Form("%s_M%i",sysNames[iSys],imult));
      
      EPart firstParticle, secondParticle;
      GetParticlesForSystem((ESys)iSys,firstParticle,secondParticle);
      
      AliFemtoV0TrackCut  *firstV0TrackCut    = GetV0TrackCut(firstParticle);
      AliFemtoV0TrackCut  *secondV0TrackCut   = GetV0TrackCut(secondParticle);
      AliFemtoESDTrackCut *firstESDTrackCut   = GetESDTrackCut(firstParticle);
      AliFemtoESDTrackCut *secondESDTrackCut  = GetESDTrackCut(secondParticle);
      
    
      femtoAnalysis[iSys]->SetEventCut(eventCut);
      femtoAnalysis[iSys]->SetNumEventsToMix(10);
      femtoAnalysis[iSys]->SetV0SharedDaughterCut(true);
      
      if(firstV0TrackCut)
        femtoAnalysis[iSys]->SetFirstParticleCut(firstV0TrackCut);
      else
        femtoAnalysis[iSys]->SetFirstParticleCut(firstESDTrackCut);
      if(secondV0TrackCut)
        femtoAnalysis[iSys]->SetSecondParticleCut(secondV0TrackCut);
      else
        femtoAnalysis[iSys]->SetSecondParticleCut(secondESDTrackCut);
      
      
      femtoAnalysis[iSys]->AddCorrFctn(separationFunction[iSys]);
      
      Manager->AddAnalysis(femtoAnalysis[iSys]);
    }
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
  if(mcAnalysis) Reader->SetReadMC(kTRUE);
  
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
  if(mcAnalysis) Reader->SetReadMC(kTRUE);
  
  return Reader;
}

AliFemtoBasicEventCut* GetEventCut()
{
  
  AliFemtoBasicEventCut *eventCut = new AliFemtoBasicEventCut();
  eventCut->SetEventMult(0,100000);
  eventCut->SetVertZPos(-8,8);
  eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
  
  return eventCut;
}

AliFemtoV0TrackCut* GetV0TrackCut(EPart particle)
{
  if(particle != kELambda && particle != kEAntiLambda) return 0;
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
  
  if(particle == kELambda)
  {
    particleCut->SetParticleType(0);
    particleCut->SetPtPosDaughter(0.5, 4.0);
    particleCut->SetPtNegDaughter(0.16, 4.0);
    particleCut->SetMinDaughtersToPrimVertex(0.1, 0.3);
    particleCut->SetInvariantMassLambda(LambdaMass-0.0038, LambdaMass+0.0043);
  }
  else if(particle == kEAntiLambda)
  {
    particleCut->SetParticleType(1);
    particleCut->SetPtPosDaughter(0.16, 4.0);
    particleCut->SetPtNegDaughter(0.3, 4.0);
    particleCut->SetMinDaughtersToPrimVertex(0.3, 0.1);
    particleCut->SetInvariantMassLambda(LambdaMass-0.0036, LambdaMass+0.0041);
  }
  return particleCut;
}

AliFemtoESDTrackCut* GetESDTrackCut(EPart particle)
{
  if(particle != kEProton && particle != kEAntiProton) return 0;
  
  AliFemtoESDTrackCut *particleCut = new AliFemtoESDTrackCut();
  
  particleCut->SetMostProbableProton();
  particleCut->SetMass(0.938272013);
  particleCut->SetEta(-0.8, 0.8);
  particleCut->SetStatus(AliESDtrack::kTPCin);
  particleCut->SetminTPCncls(80);
  particleCut->SetRemoveKinks(kTRUE);
  particleCut->SetLabel(kFALSE);
  particleCut->SetMaxTPCChiNdof(4.0);
  particleCut->SetMaxImpactXY(2.8);
  particleCut->SetMaxImpactZ(3.2);
  particleCut->SetNsigma(3.0);
  particleCut->SetNsigmaTPCTOF(kTRUE);
  particleCut->SetCharge(particle == kEProton ? 1.0 : -1.0);
  particleCut->SetPt(0.7, particle == kEProton ? 4.0 : 5.0);
  
  return particleCut;
}

void GetParticlesForSystem(ESys system, EPart &firstParticle, EPart &secondParticle)
{
  if(system == kLL)   {firstParticle = kELambda;       secondParticle = kELambda;}
  if(system == kLAL)  {firstParticle = kELambda;       secondParticle = kEAntiLambda;}
  if(system == kALAL) {firstParticle = kEAntiLambda;   secondParticle = kEAntiLambda;}
  if(system == kPL)   {firstParticle = kELambda;       secondParticle = kEProton;}
  if(system == kAPL)  {firstParticle = kELambda;       secondParticle = kEAntiProton;}
  if(system == kPAL)  {firstParticle = kEAntiLambda;   secondParticle = kEProton;}
  if(system == kAPAL) {firstParticle = kEAntiLambda;   secondParticle = kEAntiProton;}
  if(system == kPP)   {firstParticle = kEProton;       secondParticle = kEProton;}
  if(system == kPAP)  {firstParticle = kEProton;       secondParticle = kEAntiProton;}
  if(system == kAPAP) {firstParticle = kEAntiProton;   secondParticle = kEAntiProton;}
}

bool AreIdentical(ESys system)
{
  if(system == kLL || system == kALAL || system == kPP || system == kAPAP) return true;
  else return false;
}







