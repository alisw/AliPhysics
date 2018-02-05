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


enum ESys  { kAPL, kPAL, kKKpi, nSys };

const char *sysNames[nSys]      = { "APL", "PAL", "KKpi"};
const bool runSys[nSys]         = {   1  ,   1  ,   1  };
const double distMin[nSys]      = {  1.5 ,  1.5 ,  1.0 }; // everything in GeV here
const double distMax[nSys]      = {  2.5 ,  2.5 ,  2.0 };
const double distBinWidth[nSys] = { 0.001, 0.001, 0.001};

bool separationCuts;
bool ppCollisions;

enum EDecaysPionPion { kRho_770, kOmega_783, kF0_980, kF2_1270, kRho_1450, kF0_1500, kRho3_1690, kF0_1710, kF4_2050, kK0s, nDecaysPionPion };

const char*  decayPionPionName[nDecaysPionPion]  = {"rho(770)","omega(783)","f0(980)","f2(1270)","rho(1450)","f0(1500)","rho3(1690)","f0(1710)","f4(2050)","K0s"};
const double decayPionPionMass[nDecaysPionPion]  = { 0.775    , 0.783      , 0.990   , 1.276    , 1.465     , 1.504    , 1.689      , 1.723    , 2.018    , 0.497};
const double decayPionPionGamma[nDecaysPionPion] = { 0.145    , 0.008      , 0.010   , 0.187    , 0.400     , 0.109    , 0.161      , 0.139    , 0.237    , 0.005};

enum EDecaysPionKaon { kKstar_892, kKstar_1410, kKstar2_1430, kKstar_1680, kKstar3_1780, kKstar4_2045, kD0, nDecaysPionKaon };
const char*  decayPionKaonName[nDecaysPionKaon]  = {"K*(892)","K*(1410)","K*2(1430)","K*(1680)","K*3(1780)","K*4(2045)","D0"};
const double decayPionKaonMass[nDecaysPionKaon]  = { 0.294   , 1.421    , 1.425     , 1.718    , 1.776     , 2.045     , 1.865};
const double decayPionKaonGamma[nDecaysPionKaon] = { 0.050   , 0.236    , 0.099     , 0.322    , 0.159     , 0.198     , 0.010};

double decayCutWidth;

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis);
AliFemtoBaryoniaAnalysis* GetAnalysis(bool doEventMixing);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoESDTrackCut* GetTrackCut(AliFemtoTrio::EPart particle);
AliFemtoTrioCut* GetTrioCutAllDecays(ESys system);
AliFemtoTrioCut* GetTrioCutPionPionDecay(ESys system, EDecaysPionPion decay);
AliFemtoTrioCut* GetTrioCutPionKaonDecay(ESys system, EDecaysPionKaon decay);
void GetParticlesForSystem(ESys system, AliFemtoTrio::EPart &firstParticle, AliFemtoTrio::EPart &secondParticle, AliFemtoTrio::EPart &thirdParticle);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false, int year=2015, bool ppAnalysis=false, bool eventMixing=false, double cutWidth=0.8)
{
  separationCuts = sepCuts;
  ppCollisions = ppAnalysis;
  decayCutWidth = cutWidth;
  
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
  AliFemtoBaryoniaAnalysis  *baryoniaAnalysis[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  AliFemtoBasicEventCut     *eventCut[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  
  AliFemtoTrioMinvFctn *distribution[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  
  // setup analysis
  int anIter = 0;
  
  for (int iSys=0; iSys<nSys; iSys++){
    if (!runSys[iSys]) continue;
    
    //-----------------------------------------------------------------------------------
    // main distribution with all cuts applied
    //-----------------------------------------------------------------------------------
    
    // create new analysis
    baryoniaAnalysis[anIter] = GetAnalysis(eventMixing);
    
    // get event cut
    eventCut[anIter] = GetEventCut();
    
    // get particle cuts
    AliFemtoTrio::EPart firstParticle, secondParticle, thirdParticle;
    GetParticlesForSystem((ESys)iSys,firstParticle,secondParticle, thirdParticle);
    
    AliFemtoESDTrackCut *firstTrackCut  = GetTrackCut(firstParticle);
    AliFemtoESDTrackCut *secondTrackCut = GetTrackCut(secondParticle);
    AliFemtoESDTrackCut *thirdTrackCut  = GetTrackCut(thirdParticle);
    
    // get trio cut
    AliFemtoTrioCut *trioCut = GetTrioCutAllDecays((ESys)iSys);
    
    // setup anallysis cuts
    baryoniaAnalysis[anIter]->SetEventCut(eventCut[anIter]);
    baryoniaAnalysis[anIter]->SetV0SharedDaughterCut(true);
    
    baryoniaAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
    baryoniaAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
    baryoniaAnalysis[anIter]->SetThirdParticleCut(thirdTrackCut);
    
    baryoniaAnalysis[anIter]->SetCollection1type(firstParticle);
    baryoniaAnalysis[anIter]->SetCollection2type(secondParticle);
    baryoniaAnalysis[anIter]->SetCollection3type(thirdParticle);
    
    // create m_inv distribution and add to the analysis
    distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutWidth_%.2f",sysNames[iSys],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys]);
    distribution[anIter]->SetTrioCut(trioCut);
    
    baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
    
    // add analysis to the manager
    Manager->AddAnalysis(baryoniaAnalysis[anIter]);
    anIter++;
    
    //-----------------------------------------------------------------------------------
    // distributions with only selected pion-pion cuts applied
    //-----------------------------------------------------------------------------------
    
    for(int iDec=0;iDec<nDecaysPionPion;iDec++){
      baryoniaAnalysis[anIter] = GetAnalysis(eventMixing);
      eventCut[anIter] = GetEventCut();
      AliFemtoTrioCut *trioCut = GetTrioCutPionPionDecay((ESys)iSys, (EDecaysPionPion)iDec);
      baryoniaAnalysis[anIter]->SetEventCut(eventCut[anIter]);
      baryoniaAnalysis[anIter]->SetV0SharedDaughterCut(true);
      baryoniaAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
      baryoniaAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
      baryoniaAnalysis[anIter]->SetThirdParticleCut(thirdTrackCut);
      baryoniaAnalysis[anIter]->SetCollection1type(firstParticle);
      baryoniaAnalysis[anIter]->SetCollection2type(secondParticle);
      baryoniaAnalysis[anIter]->SetCollection3type(thirdParticle);
      distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutOn_%s_cutWidth_%.2f",sysNames[iSys],decayPionPionName[iDec],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys]);
      distribution[anIter]->SetTrioCut(trioCut);
      baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
      Manager->AddAnalysis(baryoniaAnalysis[anIter]);
      anIter++;
    }
    
    //-----------------------------------------------------------------------------------
    // distributions with only selected pion-kaon cuts applied
    //-----------------------------------------------------------------------------------
    
    for(int iDec=0;iDec<nDecaysPionKaon;iDec++){
      baryoniaAnalysis[anIter] = GetAnalysis(eventMixing);
      eventCut[anIter] = GetEventCut();
      AliFemtoTrioCut *trioCut = GetTrioCutPionKaonDecay((ESys)iSys, (EDecaysPionKaon)iDec);
      baryoniaAnalysis[anIter]->SetEventCut(eventCut[anIter]);
      baryoniaAnalysis[anIter]->SetV0SharedDaughterCut(true);
      baryoniaAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
      baryoniaAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
      baryoniaAnalysis[anIter]->SetThirdParticleCut(thirdTrackCut);
      baryoniaAnalysis[anIter]->SetCollection1type(firstParticle);
      baryoniaAnalysis[anIter]->SetCollection2type(secondParticle);
      baryoniaAnalysis[anIter]->SetCollection3type(thirdParticle);
      distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutOn_%s_cutWidth_%.2f",sysNames[iSys],decayPionKaonName[iDec],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys]);
      distribution[anIter]->SetTrioCut(trioCut);
      baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
      Manager->AddAnalysis(baryoniaAnalysis[anIter]);
      anIter++;
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

AliFemtoBaryoniaAnalysis* GetAnalysis(bool doEventMixing)
{
  AliFemtoBaryoniaAnalysis *analysis = new AliFemtoBaryoniaAnalysis();
  analysis->SetDoEventMixing(doEventMixing);
  // here one can put some additional analysis settings in the future
  
  return analysis;
}

AliFemtoBasicEventCut* GetEventCut()
{
  AliFemtoBasicEventCut *eventCut = new AliFemtoBasicEventCut();
  eventCut->SetEventMult(3,100000);
  eventCut->SetVertZPos(-8,8);
  eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
  
  return eventCut;
}

AliFemtoESDTrackCut* GetTrackCut(AliFemtoTrio::EPart particle)
{
  AliFemtoESDTrackCut *particleCut = new AliFemtoESDTrackCut();
  
  if(particle == AliFemtoTrio::kKaonPlus || particle==AliFemtoTrio::kKaonMinus){
    particleCut->SetMostProbableKaon();
    particleCut->SetMass(0.493677);
  }
  else if(particle == AliFemtoTrio::kPionPlus || particle==AliFemtoTrio::kPionMinus){
    particleCut->SetMostProbablePion();
    particleCut->SetMass(0.13956995);
  }
  
  if(particle == AliFemtoTrio::kKaonPlus || particle == AliFemtoTrio::kPionPlus){ particleCut->SetCharge( 1.0); }
  else                                              { particleCut->SetCharge(-1.0); }
  
  particleCut->SetPt(0.14,1.5);
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

AliFemtoTrioCut* GetTrioCutAllDecays(ESys system)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
  
  
  // ππ cuts
  for(int iDec=0;iDec<nDecaysPionPion;iDec++){
    trioCut->SetExcludePair(decayPionPionMass[iDec],decayCutWidth*decayPionPionGamma[iDec],
                            AliFemtoTrio::kPionPlus,AliFemtoTrio::kPionMinus);
  }

  // Kπ cuts:
  for(int iDec=0;iDec<nDecaysPionKaon;iDec++){
    trioCut->SetExcludePair(decayPionKaonMass[iDec],decayCutWidth*decayPionKaonGamma[iDec],
                            AliFemtoTrio::kKaonPlus,AliFemtoTrio::kPionMinus);
    trioCut->SetExcludePair(decayPionKaonMass[iDec],decayCutWidth*decayPionKaonGamma[iDec],
                            AliFemtoTrio::kKaonMinus,AliFemtoTrio::kPionPlus);
  }
  
//  pairCut->SetPhiStarDifferenceMinimum(0.012);
//  pairCut->SetEtaDifferenceMinimum(0.017);
//  pairCut->SetShareQualityMax(1.0);
//  pairCut->SetShareFractionMax(0.05);
//  pairCut->SetRemoveSameLabel(kFALSE);
//  pairCut->SetMaxEEMinv(0.002);
//  pairCut->SetMaxThetaDiff(0.008);
//  pairCut->SetDataType(AliFemtoPairCut::kAOD);
//  if(separationCuts)
//  {
//    pairCut->SetTPCEntranceSepMinimum(0.00001);
//    pairCut->SetAvgsepMinimum(5.0);
//  }
  
  return trioCut;
}

AliFemtoTrioCut* GetTrioCutPionPionDecay(ESys system, EDecaysPionPion decay)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
  trioCut->SetExcludePair(decayPionPionMass[decay],decayCutWidth*decayPionPionGamma[decay],
                            AliFemtoTrio::kPionPlus,AliFemtoTrio::kPionMinus);
  return trioCut;
}

AliFemtoTrioCut* GetTrioCutPionKaonDecay(ESys system, EDecaysPionKaon decay)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
  trioCut->SetExcludePair(decayPionKaonMass[decay],decayCutWidth*decayPionKaonGamma[decay],
                          AliFemtoTrio::kKaonPlus,AliFemtoTrio::kPionMinus);
  trioCut->SetExcludePair(decayPionKaonMass[decay],decayCutWidth*decayPionKaonGamma[decay],
                          AliFemtoTrio::kKaonMinus,AliFemtoTrio::kPionPlus);
  return trioCut;
}

void GetParticlesForSystem(ESys system, AliFemtoTrio::EPart &firstParticle, AliFemtoTrio::EPart &secondParticle, AliFemtoTrio::EPart &thirdParticle)
{
  if(system == kAPL){
    firstParticle  = AliFemtoTrio::kKaonMinus;
    secondParticle = AliFemtoTrio::kPionPlus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kPAL){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kPionPlus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kKKpi){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kKaonMinus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
}










