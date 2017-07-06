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

enum EPart { kELambda , kEAntiLambda , kEProton , kEAntiProton };

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoV0TrackCut* GetV0TrackCut(EPart particle);
AliFemtoESDTrackCut* GetESDTrackCut(EPart particle);

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
  
  AliFemtoEventAnalysis *femtoAnalysis1 = new AliFemtoEventAnalysis();
  AliFemtoEventAnalysis *femtoAnalysis2 = new AliFemtoEventAnalysis();
  AliFemtoEventAnalysis *femtoAnalysis3 = new AliFemtoEventAnalysis();
  AliFemtoEventAnalysis *femtoAnalysis4 = new AliFemtoEventAnalysis();
  
  AliFemtoBasicEventCut         *eventCut = GetEventCut();

  AliFemtoSpatialSeparationFunction *separationFunction1;
  AliFemtoSpatialSeparationFunction *separationFunction2;
  AliFemtoSpatialSeparationFunction *separationFunction3;
  AliFemtoSpatialSeparationFunction *separationFunction4;

  // get particle cuts
  AliFemtoESDTrackCut *protonCut      = GetESDTrackCut(kEProton);
  AliFemtoV0TrackCut  *lambdaCut      = GetV0TrackCut(kELambda);
  AliFemtoESDTrackCut *antiprotonCut  = GetESDTrackCut(kEAntiProton);
  AliFemtoV0TrackCut  *antilambdaCut  = GetV0TrackCut(kEAntiLambda);
  
  // setup anallysis cuts
  femtoAnalysis1->SetEventCut(eventCut);
  femtoAnalysis1->SetV0SharedDaughterCut(true);
  femtoAnalysis1->SetFirstParticleCut(protonCut);
  femtoAnalysis1->SetSecondParticleCut(lambdaCut);
  
  femtoAnalysis2->SetEventCut(eventCut);
  femtoAnalysis2->SetV0SharedDaughterCut(true);
  femtoAnalysis2->SetFirstParticleCut(antiprotonCut);
  femtoAnalysis2->SetSecondParticleCut(lambdaCut);
  
  femtoAnalysis3->SetEventCut(eventCut);
  femtoAnalysis3->SetV0SharedDaughterCut(true);
  femtoAnalysis3->SetFirstParticleCut(antiprotonCut);
  femtoAnalysis3->SetSecondParticleCut(antilambdaCut);
  
  femtoAnalysis4->SetEventCut(eventCut);
  femtoAnalysis4->SetV0SharedDaughterCut(true);
  femtoAnalysis4->SetFirstParticleCut(protonCut);
  femtoAnalysis4->SetSecondParticleCut(antilambdaCut);
  
  
  // add Average Separation correlation function
  separationFunction1 = new AliFemtoSpatialSeparationFunction("PL");
  separationFunction2 = new AliFemtoSpatialSeparationFunction("APL");
  separationFunction3 = new AliFemtoSpatialSeparationFunction("APAL");
  separationFunction4 = new AliFemtoSpatialSeparationFunction("PAL");
  
  femtoAnalysis1->AddCorrFctn(separationFunction1);
  femtoAnalysis2->AddCorrFctn(separationFunction2);
  femtoAnalysis3->AddCorrFctn(separationFunction3);
  femtoAnalysis4->AddCorrFctn(separationFunction4);
  
  Manager->AddAnalysis(femtoAnalysis1);
  Manager->AddAnalysis(femtoAnalysis2);
  Manager->AddAnalysis(femtoAnalysis3);
  Manager->AddAnalysis(femtoAnalysis4);

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










