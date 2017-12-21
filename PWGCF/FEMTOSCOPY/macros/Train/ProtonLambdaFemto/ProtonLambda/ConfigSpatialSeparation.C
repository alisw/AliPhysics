///
/// \file ConfigSpatialSeparation.C
///
/// \brief Configuration macro for the spatial separation analysis.
///
/// \author Jeremi Niedziela
/// \email jeremi.niedziela@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoEventAnalysis.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoSpatialSeparationFunction.h"
#include "AliFemtoAngularSpatialSeparationFunction.h"
#include "AliESDtrack.h"
#endif

enum EPart { kELambda , kEAntiLambda , kEProton , kEAntiProton};

// system
enum ESys { kLL , kALAL , kLAL , kPL , kAPL , kPAL , kAPAL , kPP , kPAP , kAPAP, nSys };
const char *sysNames[nSys] = { "LL", "ALAL", "LAL", "PL", "APL", "PAL", "APAL","PP","PAP","APAP" };
int runSys[nSys] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

// centrality/multiplicity
const int nMult = 10;
int multBins[nMult+1] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};
int runMult[nMult] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};

// event plane
const int nEventPlane = 7;
int runEventPlane[nEventPlane] = {1, 1, 1, 1, 1, 1, 1};

// pT binning
const int nPt = 6;
int runPt[nPt] = {1, 1, 1, 1, 1, 1};
double pTbins[nPt] = {0.0, 1.32, 1.68, 2.056, 2.624, 5.0};

// declaration of functions
AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoV0TrackCut* GetV0TrackCut(EPart particle);
AliFemtoESDTrackCut* GetESDTrackCut(EPart particle);
void GetParticlesForSystem(ESys system, EPart &firstParticle, EPart &secondParticle);
bool AreIdentical(ESys system);

// main function
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, int year=2015, bool doPtBinning = false, bool doAngular = false)
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
  
  // prepare objets for analysis and cuts
  AliFemtoEventAnalysis                  *femtoAnalysis[nSys*nMult*nEventPlane*nPt];
  AliFemtoBasicEventCut                       *eventCut[nSys*nMult*nEventPlane*nPt];
  AliFemtoSpatialSeparationFunction *separationFunction[nSys*nMult*nEventPlane*nPt];
  AliFemtoAngularSpatialSeparationFunction *angularSeparationFunction[nSys*nMult*nEventPlane*nPt];
  
  AliFemtoV0TrackCut    *firstV0TrackCut[nSys*nMult*nEventPlane*nPt];
  AliFemtoV0TrackCut   *secondV0TrackCut[nSys*nMult*nEventPlane*nPt];
  AliFemtoESDTrackCut  *firstESDTrackCut[nSys*nMult*nEventPlane*nPt];
  AliFemtoESDTrackCut *secondESDTrackCut[nSys*nMult*nEventPlane*nPt];
  
  int anIter = 0;
  
  for(int iSys=0;iSys<nSys;iSys++) // iterate on systems
  {
    if (!runSys[iSys]) continue;
    
    EPart firstParticle, secondParticle;
    GetParticlesForSystem((ESys)iSys,firstParticle,secondParticle);
    
    for (int iMult=0; iMult<nMult; iMult++) // iterate on centrality bins
    {
      if (!runMult[iMult]) continue;
     
      for (int iEventPlane=0; iEventPlane<nEventPlane; iEventPlane++) // iterate on event plane bins
      {
        if (!runEventPlane[iEventPlane]) continue;
      
        for (int iPt=0; iPt<(nPt-1); iPt++) // iterate on pT bins
        {
          if (!runPt[iPt] && doPtBinning) continue;
          if(!doPtBinning && iPt!=0) continue;
          
          // create analysis
          femtoAnalysis[anIter] = new AliFemtoEventAnalysis(multBins[iMult], multBins[iMult+1]);
          femtoAnalysis[anIter]->SetIdenticalParticles(AreIdentical((ESys)iSys));
          femtoAnalysis[anIter]->SetNumEventsToMix(10);
          femtoAnalysis[anIter]->SetV0SharedDaughterCut(true);
          
          // create a separation function
          
          separationFunction[anIter]        = new AliFemtoSpatialSeparationFunction(Form("%s_M%i_EP%i_pT%i",sysNames[iSys],iMult,iEventPlane,iPt));
          angularSeparationFunction[anIter] = new AliFemtoAngularSpatialSeparationFunction(Form("%s_M%i_EP%i_pT%i",sysNames[iSys],iMult,iEventPlane,iPt));
          
          // setup event cuts
          eventCut[anIter] = GetEventCut();
          
          if (iEventPlane == (nEventPlane-1))
            eventCut[anIter]->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
          else
            eventCut[anIter]->SetEPVZERO(-TMath::Pi()/2. +  iEventPlane    * TMath::Pi()/6.,
                                         -TMath::Pi()/2. + (iEventPlane+1) * TMath::Pi()/6.);
          
          femtoAnalysis[anIter]->SetEventCut(eventCut[anIter]);
          
          
          firstV0TrackCut[anIter]    = GetV0TrackCut(firstParticle);
          secondV0TrackCut[anIter]   = GetV0TrackCut(secondParticle);
          firstESDTrackCut[anIter]   = GetESDTrackCut(firstParticle);
          secondESDTrackCut[anIter]  = GetESDTrackCut(secondParticle);
          
          // setup particle cuts
          if(firstV0TrackCut[anIter])
          {
            if(doPtBinning) firstV0TrackCut[anIter]->SetPt(pTbins[iPt],pTbins[iPt+1]);
            femtoAnalysis[anIter]->SetFirstParticleCut(firstV0TrackCut[anIter]);
          }
          else
          {
            if(doPtBinning) firstESDTrackCut[anIter]->SetPt(pTbins[iPt],pTbins[iPt+1]);
            femtoAnalysis[anIter]->SetFirstParticleCut(firstESDTrackCut[anIter]);
          }
          if(secondV0TrackCut[anIter])
          {
            if(doPtBinning) secondV0TrackCut[anIter]->SetPt(pTbins[iPt],pTbins[iPt+1]);
            femtoAnalysis[anIter]->SetSecondParticleCut(secondV0TrackCut[anIter]);
          }
          else
          {
            if(doPtBinning) secondESDTrackCut[anIter]->SetPt(pTbins[iPt],pTbins[iPt+1]);
            femtoAnalysis[anIter]->SetSecondParticleCut(secondESDTrackCut[anIter]);
          }
          // add correlation function to analysis and analysis to the manager
          if(doAngular) femtoAnalysis[anIter]->AddCorrFctn(angularSeparationFunction[anIter]);
          else        femtoAnalysis[anIter]->AddCorrFctn(separationFunction[anIter]);
          Manager->AddAnalysis(femtoAnalysis[anIter]);
          anIter++;
        }
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







