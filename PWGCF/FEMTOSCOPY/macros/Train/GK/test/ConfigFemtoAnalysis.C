#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoMJTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoXiTrackPairCut.h"
#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoXiTrackCut.h"
#include "AliFemtoPairCutRadialDistance.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliESDtrack.h"
#include "AliFemtoCutMonitorEventMult.h"
#endif

enum ESys { kPXim, kAPXim ,kPXip, kAPXip , kLL , kALAL , kLAL , kPL , kAPL , kPAL , kAPAL , kPP , kPAP , kAPAP, nSys };
enum EPart {kEXi, kEAntiXi, kELambda , kEAntiLambda , kEProton , kEAntiProton };

const char *sysNames[nSys] = {"PXim","APXim" ,"PXip","APXip", "V0LL", "V0ALAL", "V0LAL", "V0PL", "V0APL", "V0PAL", "V0APAL","PP","PAP","APAP" };

const int nMult = 10;
int runMult[nMult] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int multBins[nMult+1] = {0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

int runSys[nSys] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

bool separationCuts;

AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis);
AliFemtoEventReaderAODChain* GetReader2011(bool mcAnalysis);
AliFemtoVertexMultAnalysis* GetAnalysis(double multMin, double multMax);
AliFemtoBasicEventCut* GetEventCut();
AliFemtoXiTrackCut* GetXiTrackCut(EPart particle);
AliFemtoV0TrackCut* GetV0TrackCut(EPart particle);
AliFemtoMJTrackCut* GetMJTrackCut(EPart particle);
AliFemtoV0PairCut* GetV0PairCut(ESys system);
AliFemtoV0PairCut* GetV0PairCut(ESys system);
AliFemtoV0TrackPairCut* GetV0TrackPairCut(ESys system);
AliFemtoPairCutRadialDistance* GetTracksPairCut(ESys system);
void GetParticlesForSystem(ESys system, EPart &firstParticle, EPart &secondParticle);




AliFemtoEventReaderAODMultSelection* GetReader2015(bool mcAnalysis)
{
  AliFemtoEventReaderAODMultSelection* Reader = new AliFemtoEventReaderAODMultSelection();
  Reader->SetFilterMask(96);
  Reader->SetReadV0(1);
  Reader->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  Reader->SetEPVZERO(kTRUE);
  Reader->SetCentralityFlattening(kTRUE);
  Reader->SetReadCascade(kTRUE);
  Reader->SetPrimaryVertexCorrectionTPCPoints(kTRUE);

  Reader->SetUseAliEventCuts(kTRUE);
  Reader->SetTrackPileUpRemoval(kTRUE);
  
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
  Reader->SetReadCascade(kTRUE);
  Reader->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  if(mcAnalysis) Reader->SetReadMC(kTRUE);
  
  return Reader;
}

AliFemtoVertexMultAnalysis* GetAnalysis(double multMin, double multMax)
{
  AliFemtoVertexMultAnalysis *analysis = new AliFemtoVertexMultAnalysis(8, -8.0, 8.0, 4, multMin, multMax);
  analysis->SetNumEventsToMix(10);
  analysis->SetMinSizePartCollection(1);
  analysis->SetVerboseMode(kFALSE);
  
  return analysis;
}

AliFemtoBasicEventCut* GetEventCut()
{
  AliFemtoBasicEventCut *eventCut = new AliFemtoBasicEventCut();
  eventCut->SetEventMult(0,100000);
  eventCut->SetVertZPos(-8,8);
  eventCut->SetEPVZERO(-TMath::Pi()/2.,TMath::Pi()/2.);
  
  return eventCut;
}

AliFemtoXiTrackCut* GetXiTrackCut(EPart particle)
{
  if(particle != kEXi && particle != kEAntiXi) return 0;
  double XiMass = 1.32171;
  double LambdaMass = 1.115683;
  
  //xi cut
  //NOTE: the SetMass call actually is important
  //      This should be set to the mass of the particle of interest, here the Xi
  //      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)
  

  AliFemtoXiTrackCut *tXiCut = new AliFemtoXiTrackCut();
  // %%%%%%%%%%%%%%%%%%%%%%%% Version 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(particle == kEXi)
    {
      tXiCut->SetParticleTypeXi(AliFemtoXiTrackCut::kXiMinus);  //kXiMinus = 0 //sprawdza Nsigma<3 czy bachelor jest pionem
      tXiCut->SetChargeXi(-1);
    }
  else if(particle == kEAntiXi)
    {
      tXiCut->SetParticleTypeXi(AliFemtoXiTrackCut::kXiPlus); //kXiPlus = 1 //sprawdza Nsigma<3 czy bachelor jest pionem
      tXiCut->SetChargeXi(1);
    }
  

  tXiCut->SetPtXi(0.5,100);
  tXiCut->SetEtaXi(0.8);
  tXiCut->SetMass(XiMass);
  
  tXiCut->SetInvariantMassXi(XiMass-0.005,XiMass+0.005); //++ bylo 006
  
  tXiCut->SetMinCosPointingAngleXi(0.97); //++ bylo 0.99
  tXiCut->SetMaxDecayLengthXi(100.);
  tXiCut->SetMaxDcaXi(100);
  tXiCut->SetInvariantMassRejectOmega(1.667,1.677);//++ NEW: omega rejection od 1.667 do 1.677 !

  tXiCut->SetIgnoreOnFlyStatus(kTRUE); 
  
  //XiDaughters
  tXiCut->SetMaxDcaXiDaughters(1.6);//++ bylo 0.3 
  tXiCut->SetRadiusXiMin(0.8); //++ NEW!
  tXiCut->SetRadiusXiMax(200); //++ NEW!



  
  if(particle == kEXi)
    {
      //Xi -> Lam Pi-									
      //Bachelor cuts (here = PiM)
      tXiCut->SetMinDcaXiBac(0.05); //++ bylo 0.03
      tXiCut->SetEtaBac(0.8);
      tXiCut->SetTPCnclsBac(70); //++
      tXiCut->SetPtBac(0.3,100.);//++
      tXiCut->SetStatusBac(AliESDtrack::kTPCrefit);  //yes or no?
      //++ brakuje cut-u z bachelor do V0 (<1.3)!!!!++++++++++
			        
					
      //Lambda cuts (regular V0)
      tXiCut->SetParticleType(AliFemtoV0TrackCut::kLambda); //0=lambda
      tXiCut->SetMinDcaV0(0.07); //++ bylo 0.1
      tXiCut->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
      tXiCut->SetMinCosPointingAngle(0.97); //++ bylo 0.998
      tXiCut->SetEta(0.8);
      tXiCut->SetPt(0.0,100);
      //tXiCut->SetOnFlyStatus(kFALSE);
      tXiCut->SetMaxV0DecayLength(100.);
      tXiCut->SetRadiusV0Min(1.4); //++ NEW!
      tXiCut->SetRadiusV0Max(200); //++ NEW!
					
      //Lambda daughter cuts
      tXiCut->SetMinDaughtersToPrimVertex(0.04,0.04); //++   pierwsza pos, druga neg, bylio (0.1,0.1); albo  pion 0.04, proton 0.03
      tXiCut->SetMaxDcaV0Daughters(1.6); //++ bylo 0.8
      tXiCut->SetEtaDaughters(0.8);    //++
      tXiCut->SetPtPosDaughter(0.3,99); //++
      tXiCut->SetPtNegDaughter(0.3,99); //++

      tXiCut->SetTPCnclsDaughters(70); //++
      tXiCut->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?
      tXiCut->SetNsigmaPosDaughter(5.0); //++
      tXiCut->SetNsigmaNegDaughter(5.0); //++
      tXiCut->SetRequireTOFProton(kFALSE); //++
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      tXiCut->SetMinvPurityAidHistoXi("XiPurityAid","XiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
      tXiCut->SetMinvPurityAidHistoV0("LambdaPurityAid","LambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);	
    }
  else if(particle == kEAntiXi)
    {
      //AXi -> ALam Pi+
      //Bachelor cuts (here = PiP)
      tXiCut->SetMinDcaXiBac(0.05);
      tXiCut->SetEtaBac(0.8);
      tXiCut->SetTPCnclsBac(70);
      tXiCut->SetPtBac(0.3,100);
      tXiCut->SetStatusBac(AliESDtrack::kTPCrefit);  //yes or no?
      
      //AntiLambda cuts (regular V0)
      tXiCut->SetParticleType(AliFemtoV0TrackCut::kAntiLambda); //1=anti-lambda
      tXiCut->SetMinDcaV0(0.07);
      tXiCut->SetInvariantMassLambda(LambdaMass-0.005,LambdaMass+0.005);
      tXiCut->SetMinCosPointingAngle(0.97);
      tXiCut->SetEta(0.8);
      tXiCut->SetPt(0.,100);
      //tXiCut->SetOnFlyStatus(kFALSE);  //CHECK kTRUE STATUS AS WELL?
      tXiCut->SetMaxV0DecayLength(100.);
      tXiCut->SetRadiusV0Min(1.4); //++ NEW!
      tXiCut->SetRadiusV0Max(200); //++ NEW!
      //Lambda daughter cuts
      tXiCut->SetMinDaughtersToPrimVertex(0.04,0.04);
      tXiCut->SetMaxDcaV0Daughters(1.6);
      tXiCut->SetEtaDaughters(0.8);
      tXiCut->SetPtPosDaughter(0.3,99); //0.16 for pions
      tXiCut->SetPtNegDaughter(0.3,99); //0.5 for anti-protons
      tXiCut->SetTPCnclsDaughters(70);
      tXiCut->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?
      tXiCut->SetNsigmaPosDaughter(5.0); //++
      tXiCut->SetNsigmaNegDaughter(5.0); //++
      tXiCut->SetRequireTOFProton(kFALSE); //++
      
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      tXiCut->SetMinvPurityAidHistoXi("AXiPurityAid","AXiMinvBeforeFinalCut",100,XiMass-0.035,XiMass+0.035);
      tXiCut->SetMinvPurityAidHistoV0("AntiLambdaPurityAid","AntiLambdaMinvBeforeFinalCut",100,LambdaMass-0.035,LambdaMass+0.035);
    }


  return tXiCut;
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
  //particleCut->SetTPCnclsDaughters(80);
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

AliFemtoMJTrackCut* GetMJTrackCut(EPart particle)
{
  if(particle != kEProton && particle != kEAntiProton) return 0;
  
  AliFemtoMJTrackCut *particleCut = new AliFemtoMJTrackCut();
  double ProtonMass = 0.938272013;
  
  particleCut->SetMostProbable(18);
  particleCut->SetMass(ProtonMass);
  particleCut->SetEta(-0.8, 0.8);
  //particleCut->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  //particleCut->SetminTPCncls(80);
  //particleCut->SetRemoveKinks(kTRUE);
  //particleCut->SetLabel(kFALSE);
  //particleCut->SetMaxTPCChiNdof(4.0);
  //particleCut->SetMaxImpactXY(2.8);
  //particleCut->SetMaxImpactZ(3.2);
  particleCut->SetNsigma(3.0);
  particleCut->SetNsigma2(3.0);
  particleCut->SetNsigmaTPCTOF(kTRUE);
  particleCut->SetElectronRejection(kTRUE);
  particleCut->SetCharge(particle == kEProton ? 1.0 : -1.0);
  particleCut->SetPt(0.7, particle == kEProton ? 4.0 : 5.0);
  
  return particleCut;
}

AliFemtoXiTrackPairCut* GetXiTrackPairCut(ESys system)
{
  AliFemtoXiTrackPairCut *pairCut = 0;
  
  if(system == kPXim || system == kAPXim || system == kPXip || system == kAPXip)
  {

    pairCut = new AliFemtoXiTrackPairCut();
  }
  return pairCut;
}


AliFemtoV0PairCut* GetV0PairCut(ESys system)
{
  AliFemtoV0PairCut *pairCut = 0;
  
  if(system == kLL || system == kLAL || system == kALAL)
  {
    pairCut = new AliFemtoV0PairCut();
    pairCut->SetDataType(AliFemtoPairCut::kAOD);
    pairCut->SetTPCEntranceSepMinimum(0.00001);
    pairCut->SetTPCExitSepMinimum(-1.);
    if(separationCuts)
    {
      pairCut->SetMinAvgSeparation(0, 5.5); //proton-pion+
      pairCut->SetMinAvgSeparation(1, 5.5); //proton-antiproton
      pairCut->SetMinAvgSeparation(2, 0); //pion- - pion+
      pairCut->SetMinAvgSeparation(3, 5.5); //antiproton - pion-
    }
  }
  return pairCut;
}

AliFemtoV0TrackPairCut* GetV0TrackPairCut(ESys system)
{
  AliFemtoV0TrackPairCut *pairCut = 0;
  
  if(system == kPL || system == kAPAL)
  {
    pairCut = new AliFemtoV0TrackPairCut(); //lambda-proton
    pairCut->SetShareQualityMax(1.0); //between V0 daughter and track
    pairCut->SetShareFractionMax(0.05);
    pairCut->SetTPCOnly(kTRUE);
    pairCut->SetDataType(AliFemtoPairCut::kAOD);
    
    if(system == kPL)
    {
      pairCut->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
    }
    else if(system == kAPAL)
    {
      pairCut->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
    }
    
    if(separationCuts)
    {
      pairCut->SetTPCEntranceSepMinimum(0.00001);
      pairCut->SetTPCExitSepMinimum(-1.);
      pairCut->SetMinAvgSeparation(0, 5); //0 - track-pos, 1 - track-neg
      pairCut->SetMinAvgSeparation(1, 5);
    }
  }
  else if(system == kAPL || system == kPAL)
  {
    pairCut = new AliFemtoV0TrackPairCut(); //lambda-antiproton, antilambda-proton
    pairCut->SetShareQualityMax(1.0); //between V0 daughter and track
    pairCut->SetShareFractionMax(0.05);
    pairCut->SetTPCOnly(kTRUE);
    pairCut->SetDataType(AliFemtoPairCut::kAOD);
    if(separationCuts)
    {
      pairCut->SetTPCEntranceSepMinimum(0.00001);
      pairCut->SetTPCExitSepMinimum(-1.);
      if(system == kAPL)
      {
        pairCut->SetMinAvgSeparation(0, 5);  // antiproton - proton
        pairCut->SetMinAvgSeparation(1, 8); // antiproton - pion-
      }
      else if(system == kPAL)
      {
        pairCut->SetMinAvgSeparation(0, 8); // proton - pion+
        pairCut->SetMinAvgSeparation(1, 5); // proton - antiproton-
      }
    }
  }
  return pairCut;
}

AliFemtoPairCutRadialDistance* GetTracksPairCut(ESys system)
{
  AliFemtoPairCutRadialDistance *pairCut = 0;
  
  if(system == kPP || system == kPAP || system == kAPAP)
  {
    pairCut = new AliFemtoPairCutRadialDistance();
    pairCut->SetPhiStarDifferenceMinimum(0.012);
    pairCut->SetEtaDifferenceMinimum(0.017);
    pairCut->SetShareQualityMax(1.0);
    pairCut->SetShareFractionMax(0.05);
    pairCut->SetRemoveSameLabel(kFALSE);
    pairCut->SetMaxEEMinv(0.002);
    pairCut->SetMaxThetaDiff(0.008);
    pairCut->SetDataType(AliFemtoPairCut::kAOD);
    if(separationCuts)
    {
      pairCut->SetTPCEntranceSepMinimum(0.00001);
      pairCut->SetAvgsepMinimum(5.0);
    }
    
  }
  return pairCut;
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
  if(system == kPXim) {firstParticle = kEXi; secondParticle = kEProton;}
  if(system == kAPXim) {firstParticle = kEXi; secondParticle = kEAntiProton;}
  if(system == kPXip) {firstParticle = kEAntiXi; secondParticle = kEProton;}
  if(system == kAPXip) {firstParticle = kEAntiXi; secondParticle = kEAntiProton;}

}


//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false, int year=2015)
{
  separationCuts = sepCuts;
  
  // create analysis managers
  AliFemtoManager* Manager=new AliFemtoManager();
  AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
  
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
  
  // declare necessary objects
  AliFemtoVertexMultAnalysis    *femtoAnalysis[nSys*nMult];
  AliFemtoBasicEventCut         *eventCut[nSys*nMult];
  
  AliFemtoAvgSepCorrFctn        *avgSepCF[nSys*nMult];
  AliFemtoCorrFctnNonIdDR       *nonIdenticalCF[nSys*nMult];
  AliFemtoQinvCorrFctn          *identicalCF[nSys*nMult];
  AliFemtoModelCorrFctn         *modelCF[nSys*nMult];
  AliFemtoPairOriginMonitor     *pairOriginPass[nSys*nMult];
  AliFemtoPairOriginMonitor     *pairOriginFail[nSys*nMult];
  
  AliFemtoCutMonitorEventMult   *eventMultMonitorPass[nSys*nMult];
  AliFemtoCutMonitorEventMult   *eventMultMonitorFail[nSys*nMult];
  
  // setup analysis
  int anIter = 0;
  for (int imult=0; imult<nMult; imult++)
  {
    if (!runMult[imult]) continue;
    
    for (int iSys=0; iSys<nSys; iSys++)
    {
      if (!runSys[iSys]) continue;
      
      anIter = imult * nSys + iSys;
      
      // create new analysis
      femtoAnalysis[anIter] = GetAnalysis(multBins[imult], multBins[imult+1]);
      
      // get event cut
      eventCut[anIter] = GetEventCut();
      
      // get particle cuts
      EPart firstParticle, secondParticle;
      GetParticlesForSystem((ESys)iSys,firstParticle,secondParticle);

      AliFemtoXiTrackCut  *firstXiTrackCut    = GetXiTrackCut(firstParticle);
      AliFemtoXiTrackCut  *secondXiTrackCut   = GetXiTrackCut(secondParticle);  
      AliFemtoV0TrackCut  *firstV0TrackCut    = GetV0TrackCut(firstParticle);
      AliFemtoV0TrackCut  *secondV0TrackCut   = GetV0TrackCut(secondParticle);
      AliFemtoMJTrackCut *firstMJTrackCut   = GetMJTrackCut(firstParticle);
      AliFemtoMJTrackCut *secondMJTrackCut  = GetMJTrackCut(secondParticle);
      
      // get pair cut
      AliFemtoV0PairCut               *V0pairCut      = GetV0PairCut((ESys)iSys);
      AliFemtoXiTrackPairCut          *XitrackPairCut = GetXiTrackPairCut((ESys)iSys);
      AliFemtoV0TrackPairCut          *V0trackPairCut = GetV0TrackPairCut((ESys)iSys);
      AliFemtoPairCutRadialDistance   *tracksPairCut  = GetTracksPairCut((ESys)iSys);

      // create monitors
      if(mcAnalysis)
      {
        pairOriginPass[anIter] = new AliFemtoPairOriginMonitor(Form("Pass%stpcM%i", sysNames[iSys], imult));
        pairOriginFail[anIter] = new AliFemtoPairOriginMonitor(Form("Fail%stpcM%i", sysNames[iSys], imult));
      }
      eventMultMonitorPass[anIter] = new AliFemtoCutMonitorEventMult(Form("Pass%sM%i",sysNames[iSys],imult));
      eventMultMonitorFail[anIter] = new AliFemtoCutMonitorEventMult(Form("Fail%sM%i",sysNames[iSys],imult));
      eventCut[anIter]->AddCutMonitor(eventMultMonitorPass[anIter],eventMultMonitorFail[anIter]);
      
      // setup anallysis cuts
      femtoAnalysis[anIter]->SetEventCut(eventCut[anIter]);
      femtoAnalysis[anIter]->SetV0SharedDaughterCut(true);
      femtoAnalysis[anIter]->SetEnablePairMonitors(false);

      if(firstXiTrackCut)
	 femtoAnalysis[anIter]->SetFirstParticleCut(firstXiTrackCut);
      else if(firstV0TrackCut)
        femtoAnalysis[anIter]->SetFirstParticleCut(firstV0TrackCut);
      else
        femtoAnalysis[anIter]->SetFirstParticleCut(firstMJTrackCut);
      if(secondXiTrackCut)
        femtoAnalysis[anIter]->SetSecondParticleCut(secondXiTrackCut);
      else if(secondV0TrackCut)
        femtoAnalysis[anIter]->SetSecondParticleCut(secondV0TrackCut);
      else
        femtoAnalysis[anIter]->SetSecondParticleCut(secondMJTrackCut);



      if(XitrackPairCut)
	{
	  if(mcAnalysis) XitrackPairCut->AddCutMonitor(pairOriginPass[anIter], pairOriginFail[anIter]);

	  femtoAnalysis[anIter]->SetPairCut(XitrackPairCut);
	}
      else if(V0pairCut){
        if(mcAnalysis) V0pairCut->AddCutMonitor(pairOriginPass[anIter], pairOriginFail[anIter]);
        femtoAnalysis[anIter]->SetPairCut(V0pairCut);
        
      }
      else if(V0trackPairCut){
        if(mcAnalysis) V0trackPairCut->AddCutMonitor(pairOriginPass[anIter], pairOriginFail[anIter]);
        femtoAnalysis[anIter]->SetPairCut(V0trackPairCut);
      }
      else
      {
        if(mcAnalysis) tracksPairCut->AddCutMonitor(pairOriginPass[anIter], pairOriginFail[anIter]);
        femtoAnalysis[anIter]->SetPairCut(tracksPairCut);
      }
      
      if(iSys==kLL || iSys == kALAL || iSys == kLAL || iSys == kPL || iSys == kAPL || iSys == kPAL || iSys == kAPAL)
      {
        femtoAnalysis[anIter]->SetV0SharedDaughterCut(true);
      }
            
      // add Average Separation correlation function
      /*
      avgSepCF[anIter] = new AliFemtoAvgSepCorrFctn(Form("Avgsep%stpcM%iPsi6", sysNames[iSys], imult),5000,0,500);
      
      if(iSys == kLL || iSys == kALAL || iSys == kLAL)
        avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
      else if(iSys == kPL || iSys == kAPL || iSys == kPAL || iSys == kAPAL)
        avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
      else if(iSys ==kPP || iSys ==kPAP || iSys == kAPAP)
        avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
      
      femtoAnalysis[anIter]->AddCorrFctn(avgSepCF[anIter]);
      */
      // add femtoscopic correlation function (identical or non-identical masses)
      if(iSys == kPXim || iSys == kAPXim || iSys == kPXip || iSys == kAPXip || iSys==kPL || iSys==kAPL || iSys==kPAL || iSys==kAPAL)
      {
        // this is in k*
        nonIdenticalCF[anIter] = new AliFemtoCorrFctnNonIdDR(Form("CF_kstar_%sM%i", sysNames[iSys], imult),500,0.0,5.0);
        femtoAnalysis[anIter]->AddCorrFctn(nonIdenticalCF[anIter]);
      }
      else // pairs with the same mass of particles
      {
        // this is in q_inv
        identicalCF[anIter] = new AliFemtoQinvCorrFctn(Form("CF_qinv_%sM%i", sysNames[iSys], imult),500,0.0,10.0);
        femtoAnalysis[anIter]->AddCorrFctn(identicalCF[anIter]);
      }
      
      // add correlation function on model data
      if(mcAnalysis)
      {
        modelCF[anIter] = new AliFemtoModelCorrFctn(Form("CF_qinv_Model_%sM%i", sysNames[iSys],imult),500,0.0,10.0);
        modelCF[anIter]->ConnectToManager(modelMgr);
        femtoAnalysis[anIter]->AddCorrFctn(modelCF[anIter]);
        femtoAnalysis[anIter]->SetEnablePairMonitors(true);
      }
      
      Manager->AddAnalysis(femtoAnalysis[anIter]);
    }
    
  }
  return Manager;
}










