#if !defined(__CINT__) || defined(__MAKECINT_)
#include "AliFemtoManager.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoV0PairCut.h"
#include "AliFemtoV0TrackPairCut.h"
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoCorrFctnNonIdDR.h"
#include "AliFemtoAvgSepCorrFctn.h"
#endif

enum ESys { kLL , kALAL , kLAL , kPL , kAPL , kPAL , kAPAL , kPP , kPAP , kAPAP, nSys };
enum EPart { kLambda , kAntiLambda , kProton , kAntiProton };

const char *sysNames[nSys] = { "V0LL", "V0ALAL", "V0LAL", "V0PL", "V0APL", "V0PAL", "V0APAL","PP","PAP","APAP" };

const int nMult = 10;
int runMult[nMult] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
int multBins[nMult+1] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

int runSys[nSys] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

bool separationCuts;

//________________________________________________________________________
AliFemtoManager* ConfigFemtoAnalysis(bool mcAnalysis=false, bool sepCuts=false)
{
    separationCuts = sepCuts;
    
    // create analysis managers
    AliFemtoManager* Manager=new AliFemtoManager();
    AliFemtoModelManager *modelMgr = new AliFemtoModelManager();
    
    // add event reader
    AliFemtoEventReaderAODMultSelection* Reader = GetReader(mcAnalysis);
    Manager->SetEventReader(Reader);
    
    // declare necessary objects
    AliFemtoVertexMultAnalysis    *femtoAnalysis[nSys*nMult];
    AliFemtoBasicEventCut         *eventCut[nSys*nMult];

    AliFemtoAvgSepCorrFctn        *avgSepCF[nSys*nMult];
    AliFemtoCorrFctnNonIdDR       *nonIdenticalCF[nSys*nMult];
    AliFemtoQinvCorrFctn          *identicalCF[nSys*nMult];
    AliFemtoModelCorrFctn         *modelCF[nSys*nMult];
    
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
            GetParticlesForSystem(iSys,firstParticle,secondParticle);
            
            AliFemtoV0TrackCut  *firstV0TrackCut    = GetV0TrackCut(firstParticle);
            AliFemtoV0TrackCut  *secondV0TrackCut   = GetV0TrackCut(secondParticle);
            AliFemtoESDTrackCut *firstESDTrackCut   = GetESDTrackCut(firstParticle);
            AliFemtoESDTrackCut *secondESDTrackCut  = GetESDTrackCut(secondParticle);
            
            // get pair cut
            AliFemtoV0PairCut               *V0pairCut      = GetV0PairCut(iSys);
            AliFemtoV0TrackPairCut          *V0trackPairCut = GetV0TrackPairCut(iSys);
            AliFemtoPairCutRadialDistance   *tracksPairCut  = GetTracksPairCut(iSys);
            
            // setup anallysis cuts
            femtoAnalysis[anIter]->SetEventCut(eventCut[anIter]);
            femtoAnalysis[anIter]->SetV0SharedDaughterCut(true);
            femtoAnalysis[anIter]->SetEnablePairMonitors(false);
            femtoAnalysis[anIter]->SetFirstParticleCut( firstV0TrackCut  ? firstV0TrackCut  : firstESDTrackCut);
            femtoAnalysis[anIter]->SetSecondParticleCut(secondV0TrackCut ? secondV0TrackCut : secondESDTrackCut);
            femtoAnalysis[anIter]->SetPairCut(V0pairCut ? V0pairCut : (V0trackPairCut ? V0trackPairCut : tracksPairCut));
            
            if(iSys==kLL || iSys == kALAL || iSys == kLAL || iSys == kPL || iSys == kAPL || iSys == kPAL || iSys == kAPAL)
            {
                femtoAnalysis[anIter]->SetV0SharedDaughterCut(true);
            }
            
            // add Average Separation correlation function
            avgSepCF[anIter] = new AliFemtoAvgSepCorrFctn(Form("Avgsep%stpcM%iPsi6", sysNames[iSys], imult),5000,0,500);
            
            if(iSys == kLL || iSys == kALAL || iSys == kLAL)
                avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
            else if(iSys == kPL || iSys == kAPL || iSys == kPAL || iSys == kAPAL)
                avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
            else if(iSys ==kPP || iSys ==kPAP || iSys == kAPAP)
                avgSepCF[anIter]->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
            
            femtoAnalysis[anIter]->AddCorrFctn(avgSepCF[anIter]);
            
            // add femtoscopic correlation function (identical or non-identical)
            if(iSys==kPL || iSys==kAPL || iSys==kPAL || iSys==kAPAL || iSys==kPAP || iSys==kLAL)
            {
                nonIdenticalCF[anIter] = new AliFemtoCorrFctnNonIdDR(Form("cnonid%stpcM%iPsi6", sysNames[iSys], imult), 100, 0.0,1.0);
                femtoAnalysis[anIter]->AddCorrFctn(nonIdenticalCF[anIter]);
            }
            else
            {
                identicalCF[anIter] = new AliFemtoQinvCorrFctn(Form("cqinv%stpcM%iPsi6", sysNames[iSys], imult),100,0.0,1.0);
                femtoAnalysis[anIter]->AddCorrFctn(identicalCF[anIter]);
            }
        
            // add correlation function on model data
            if(mcAnalysis)
            {
                modelCF[anIter] = new AliFemtoModelCorrFctn(Form("cQinv_Model_%s_M%i", sysNames[iSys],imult), 400, 0, 2);
                modelCF[anIter]->ConnectToManager(modelMgr);
                femtoAnalysis[anIter]->AddCorrFctn(modelCF[anIter]);
            }
            
            Manager->AddAnalysis(femtoAnalysis[anIter]);
        }
        
    }
    return Manager;
}

AliFemtoEventReaderAODMultSelection* GetReader(bool mcAnalysis)
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

AliFemtoV0TrackCut* GetV0TrackCut(EPart particle)
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

AliFemtoESDTrackCut* GetESDTrackCut(EPart particle)
{
    if(particle != kProton && particle != kAntiProton) return 0;
    
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
    particleCut->SetCharge(particle == kProton ? 1.0 : -1.0);
    particleCut->SetPt(0.7, particle == kProton ? 4.0 : 5.0);
    
    return particleCut;
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
    
    if(system == kPL)
    {
        pairCut = new AliFemtoV0TrackPairCut(); //lambda-proton
        pairCut->SetShareQualityMax(1.0); //between V0 daughter and track
        pairCut->SetShareFractionMax(0.05);
        pairCut->SetTPCOnly(kTRUE);
        pairCut->SetDataType(AliFemtoPairCut::kAOD);
        pairCut->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kLambda,AliFemtoV0TrackPairCut::kProton); //0 - lambda, 2 - proton
        if(separationCuts)
        {
            pairCut->SetTPCEntranceSepMinimum(0.00001);
            pairCut->SetTPCExitSepMinimum(-1.);
            pairCut->SetMinAvgSeparation(0, 5); //0 - track-pos, 1 - track-neg
            pairCut->SetMinAvgSeparation(1, 5);
        }
    }
    else if(system == kAPAL)
    {
        pairCut = new AliFemtoV0TrackPairCut(); //antilambda-antiproton
        pairCut->SetShareQualityMax(1.0); //between V0 daughter and track
        pairCut->SetShareFractionMax(0.05);
        pairCut->SetTPCOnly(kTRUE);
        pairCut->SetDataType(AliFemtoPairCut::kAOD);
        pairCut->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
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
    if(system == kLL)   {firstParticle = kLambda;       secondParticle = kLambda;}
    if(system == kLAL)  {firstParticle = kLambda;       secondParticle = kAntiLambda;}
    if(system == kALAL) {firstParticle = kAntiLambda;   secondParticle = kAntiLambda;}
    if(system == kPL)   {firstParticle = kLambda;       secondParticle = kProton;}
    if(system == kAPL)  {firstParticle = kLambda;       secondParticle = kAntiProton;}
    if(system == kPAL)  {firstParticle = kAntiLambda;   secondParticle = kProton;}
    if(system == kAPAL) {firstParticle = kAntiLambda;   secondParticle = kAntiProton;}
    if(system == kPP)   {firstParticle = kProton;       secondParticle = kProton;}
    if(system == kPAP)  {firstParticle = kProton;       secondParticle = kAntiProton;}
    if(system == kAPAP) {firstParticle = kAntiProton;   secondParticle = kAntiProton;}
}










