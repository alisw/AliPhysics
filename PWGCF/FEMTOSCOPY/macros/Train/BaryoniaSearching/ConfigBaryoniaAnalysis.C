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


enum ESys  { kAPL, kPAL, kApXz, kPrAxz, kKKpi, kDplus, nSys };

const char *sysNames[nSys]      = { "APL", "PAL", "ApXz", "PrAxz",  "D+" , "KKpi"};
const bool runSys[nSys]         = {   1  ,   1  ,   0   ,    0   ,    1  ,   0   };
const double distMin[nSys]      = {  1.5 ,  1.5 ,  1.5  ,   1.5  ,   1.0 ,  1.0  }; // everything in GeV here
const double distMax[nSys]      = {  2.5 ,  2.5 ,  3.0  ,   3.0  ,   2.8 ,  2.0  };
const double distBinWidth[nSys] = { 0.001, 0.001, 0.001 ,  0.001 , 0.001 , 0.001 };

const double dalitzCutMass[nSys]= { 2.148, 2.148, 2.330 ,  2.330 ,1.86959,   0   };
const double dalitzCutGamma[nSys]={ 0.134, 0.134, 0.100 ,  0.100 ,0.00063,   0   };

const double dalitzMin       = 0.0;
const double dalitzMax       = 6.0;
const double dalitzBinWidth  = 0.005;

bool separationCuts;
bool ppCollisions;
double decayCutWidth;
bool doAngles;

enum EDecaysPionPion { kRho_770, kOmega_783, kF0_980, kF2_1270, kRho_1450, kF0_1500, kRho3_1690, kF0_1710, kF4_2050, kK0s, nDecaysPionPion };

const bool   doDecayPionPion[nDecaysPionPion]    = {   0      ,    1       ,    0    ,    1     ,    1      ,    1     ,    1       ,     1    ,    1     ,  1  };
const char*  decayPionPionName[nDecaysPionPion]  = {"rho(770)","omega(783)","f0(980)","f2(1270)","rho(1450)","f0(1500)","rho3(1690)","f0(1710)","f4(2050)","K0s"};
const double decayPionPionMass[nDecaysPionPion]  = { 0.775    , 0.783      , 0.990   , 1.276    , 1.465     , 1.504    , 1.689      , 1.723    , 2.018    , 0.497};
const double decayPionPionGamma[nDecaysPionPion] = { 0.145    , 0.008      , 0.010   , 0.187    , 0.400     , 0.109    , 0.161      , 0.139    , 0.237    , 0.005};

enum EDecaysPionKaon { kKstar_892, kKstar_1410, kKstar2_1430, kKstar_1680, kKstar3_1780, kKstar4_2045, kD0, nDecaysPionKaon };

const bool   doDecayPionKaon[nDecaysPionKaon]    = {   0     ,    1     ,     1     ,    1     ,    1      ,    1     ,   0 };
const char*  decayPionKaonName[nDecaysPionKaon]  = {"K*(892)","K*(1410)","K*2(1430)","K*(1680)","K*3(1780)","K*4(2045)","D0"};
const double decayPionKaonMass[nDecaysPionKaon]  = { 0.294   , 1.421    , 1.425     , 1.718    , 1.776     , 2.045     , 1.865};
const double decayPionKaonGamma[nDecaysPionKaon] = { 0.050   , 0.236    , 0.099     , 0.322    , 0.159     , 0.198     , 0.010};




AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReaderPP(bool mcAnalysis);
AliFemtoBaryoniaAnalysis* GetAnalysis(bool doEventMixing);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoESDTrackCut* GetTrackCut(AliFemtoTrio::EPart particle);
AliFemtoTrioCut* GetTrioCutAllDecays(ESys system, double cutWidth);
AliFemtoTrioCut* GetTrioCutDalitz(ESys system, bool twoBody=false, double widthFraction=1.0);
AliFemtoTrioCut* GetTrioCutPionPionDecay(ESys system, EDecaysPionPion decay);
AliFemtoTrioCut* GetTrioCutPionKaonDecay(ESys system, EDecaysPionKaon decay);
void GetParticlesForSystem(ESys system, AliFemtoTrio::EPart &firstParticle, AliFemtoTrio::EPart &secondParticle, AliFemtoTrio::EPart &thirdParticle);

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false, int year=2015, bool ppAnalysis=false, bool eventMixing=false, double cutWidth=0.8, bool pionPionDecays=false, bool pionKaonDecays=false, bool _doAngles=false)
{
  separationCuts = sepCuts;
  ppCollisions = ppAnalysis;
  decayCutWidth = cutWidth;
  doAngles = _doAngles;
  
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
  
  AliFemtoTrioMinvFctn    *distribution[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  AliFemtoTrioMinvFctn    *distributionCutUp[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  AliFemtoTrioMinvFctn    *distributionCutDown[nSys*(nDecaysPionPion+nDecaysPionKaon+1)];
  
  AliFemtoTrioMinvFctn    *dalitz[nSys];
  AliFemtoTrioMinvFctn    *dalitzCutUp[nSys];
  AliFemtoTrioMinvFctn    *dalitzCutDown[nSys];
  
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
    
    // setup anallysis cuts
    baryoniaAnalysis[anIter]->SetEventCut(eventCut[anIter]);
    baryoniaAnalysis[anIter]->SetV0SharedDaughterCut(true);
    
    baryoniaAnalysis[anIter]->SetFirstParticleCut(firstTrackCut);
    baryoniaAnalysis[anIter]->SetSecondParticleCut(secondTrackCut);
    baryoniaAnalysis[anIter]->SetThirdParticleCut(thirdTrackCut);
    
    baryoniaAnalysis[anIter]->SetCollection1type(firstParticle);
    baryoniaAnalysis[anIter]->SetCollection2type(secondParticle);
    baryoniaAnalysis[anIter]->SetCollection3type(thirdParticle);
    
    bool doMinv = true;
    bool doDalitz = false;
    
    // create m_inv distribution and add to the analysis
    
    // get trio cut
    AliFemtoTrioCut *trioCut = GetTrioCutAllDecays((ESys)iSys,decayCutWidth);
    distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutWidth_%.2f",sysNames[iSys],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys],doMinv,doDalitz);
    distribution[anIter]->SetTrioCut(trioCut);
    
    AliFemtoTrioCut *trioCutUp = GetTrioCutAllDecays((ESys)iSys,decayCutWidth+0.10);
    distributionCutUp[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutWidth_%.2f",sysNames[iSys],decayCutWidth+0.1),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys],doMinv,doDalitz);
    distributionCutUp[anIter]->SetTrioCut(trioCutUp);
    
    AliFemtoTrioCut *trioCutDown = GetTrioCutAllDecays((ESys)iSys,decayCutWidth-0.10);
    distributionCutDown[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutWidth_%.2f",sysNames[iSys],decayCutWidth-0.1),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys],doMinv,doDalitz);
    distributionCutDown[anIter]->SetTrioCut(trioCutDown);
    
    doMinv = false;
    doDalitz = true;
    bool twoBodyCuts = false;
    
    AliFemtoTrioCut *dalitzTrioCut = GetTrioCutDalitz((ESys)iSys,twoBodyCuts,1.0);
    dalitz[iSys] = new AliFemtoTrioMinvFctn(Form("%s",sysNames[iSys]),
                                            (dalitzMax-dalitzMin)/dalitzBinWidth,dalitzMin,dalitzMax,
                                            doMinv,doDalitz,doAngles);
    dalitz[iSys]->SetTrioCut(dalitzTrioCut);
    
    AliFemtoTrioCut *dalitzTrioCutUp = GetTrioCutDalitz((ESys)iSys,twoBodyCuts,1.5);
    dalitzCutUp[iSys] = new AliFemtoTrioMinvFctn(Form("%s_cutUp",sysNames[iSys]),
                                                 (dalitzMax-dalitzMin)/dalitzBinWidth,dalitzMin,dalitzMax,
                                                 doMinv,doDalitz,doAngles);
    dalitzCutUp[iSys]->SetTrioCut(dalitzTrioCutUp);
    
    AliFemtoTrioCut *dalitzTrioCutDown = GetTrioCutDalitz((ESys)iSys,twoBodyCuts,0.5);
    dalitzCutDown[iSys] = new AliFemtoTrioMinvFctn(Form("%s_cutDown",sysNames[iSys]),
                                                   (dalitzMax-dalitzMin)/dalitzBinWidth,dalitzMin,dalitzMax,
                                                   doMinv,doDalitz,doAngles);
    dalitzCutDown[iSys]->SetTrioCut(dalitzTrioCutDown);
    
    baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
    baryoniaAnalysis[anIter]->AddDistribution(distributionCutUp[anIter]);
    baryoniaAnalysis[anIter]->AddDistribution(distributionCutDown[anIter]);
    baryoniaAnalysis[anIter]->AddDistribution(dalitz[iSys]);
    baryoniaAnalysis[anIter]->AddDistribution(dalitzCutUp[iSys]);
    baryoniaAnalysis[anIter]->AddDistribution(dalitzCutDown[iSys]);
    
    // add analysis to the manager
    Manager->AddAnalysis(baryoniaAnalysis[anIter]);
    anIter++;
    
    //-----------------------------------------------------------------------------------
    // distributions with only selected pion-pion cuts applied
    //-----------------------------------------------------------------------------------
    
    doMinv = true;
    doDalitz = false;
    
    if(pionPionDecays){
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
        distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutOn_%s_cutWidth_%.2f",sysNames[iSys],decayPionPionName[iDec],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys],doMinv,doDalitz);
        distribution[anIter]->SetTrioCut(trioCut);
        baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
        Manager->AddAnalysis(baryoniaAnalysis[anIter]);
        anIter++;
      }
    }
      
    //-----------------------------------------------------------------------------------
    // distributions with only selected pion-kaon cuts applied
    //-----------------------------------------------------------------------------------
    if(pionKaonDecays){
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
        distribution[anIter] = new AliFemtoTrioMinvFctn(Form("%s_cutOn_%s_cutWidth_%.2f",sysNames[iSys],decayPionKaonName[iDec],decayCutWidth),(distMax[iSys]-distMin[iSys])/distBinWidth[iSys],distMin[iSys],distMax[iSys],doMinv,doDalitz);
        distribution[anIter]->SetTrioCut(trioCut);
        baryoniaAnalysis[anIter]->AddDistribution(distribution[anIter]);
        Manager->AddAnalysis(baryoniaAnalysis[anIter]);
        anIter++;
      }
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

AliFemtoTrioCut* GetTrioCutAllDecays(ESys system, double cutWidth)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
  
  
  // ππ cuts
  for(int iDec=0;iDec<nDecaysPionPion;iDec++){
    if(!doDecayPionPion[iDec]) continue;
    
    trioCut->SetExcludePair(decayPionPionMass[iDec],cutWidth*decayPionPionGamma[iDec],
                            AliFemtoTrio::kPionPlus,AliFemtoTrio::kPionMinus);
  }

  // Kπ cuts:
  for(int iDec=0;iDec<nDecaysPionKaon;iDec++){
    if(!doDecayPionKaon[iDec]) continue;
    
    trioCut->SetExcludePair(decayPionKaonMass[iDec],cutWidth*decayPionKaonGamma[iDec],
                            AliFemtoTrio::kKaonPlus,AliFemtoTrio::kPionMinus);
    trioCut->SetExcludePair(decayPionKaonMass[iDec],cutWidth*decayPionKaonGamma[iDec],
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

AliFemtoTrioCut* GetTrioCutDalitz(ESys system, bool twoBody, double widthFraction)
{
  AliFemtoTrioCut *trioCut = new AliFemtoTrioCut();
    trioCut->SetIncludeTrioOnly(dalitzCutMass[system],widthFraction*dalitzCutGamma[system]);
  
  if(twoBody){
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
  }
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
  if(system == kApXz){
    firstParticle  = AliFemtoTrio::kKaonMinus;
    secondParticle = AliFemtoTrio::kKaonMinus;
    thirdParticle  = AliFemtoTrio::kPionPlus;
  }
  if(system == kPrAxz){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kKaonPlus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kKKpi){
    firstParticle  = AliFemtoTrio::kKaonPlus;
    secondParticle = AliFemtoTrio::kKaonMinus;
    thirdParticle  = AliFemtoTrio::kPionMinus;
  }
  if(system == kDplus){
    firstParticle  = AliFemtoTrio::kKaonMinus;
    secondParticle = AliFemtoTrio::kPionPlus;
    thirdParticle  = AliFemtoTrio::kPionPlus;
  }
}










