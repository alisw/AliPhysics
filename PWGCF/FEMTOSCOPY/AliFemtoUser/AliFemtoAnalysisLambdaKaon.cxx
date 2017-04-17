///
/// \file AliFemtoAnalysisLambdaKaon.cxx
///

#include "AliFemtoAnalysisLambdaKaon.h"
#include "TObjArray.h"
#include "AliESDtrack.h"
#ifdef __ROOT__
ClassImp(AliFemtoAnalysisLambdaKaon)
#endif

static const double PionMass = 0.13957018,
                    KchMass = 0.493677,
                    K0ShortMass = 0.497614,
                    ProtonMass = 0.938272013,
                    LambdaMass = 1.115683,
		    XiMass     = 1.32171;

//____________________________
const char* const AliFemtoAnalysisLambdaKaon::fAnalysisTags[] = {"LamK0", "ALamK0", "LamKchP", "ALamKchP", "LamKchM", "ALamKchM", "LamLam", "ALamALam", "LamALam", "LamPiP", "ALamPiP", "LamPiM", "ALamPiM", "XiKchP", "AXiKchP", "XiKchM", "AXiKchM", "ProtPiM", "AProtPiP", "PiPPiM"};



//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AliFemtoAnalysisLambdaKaon::AnalysisType aAnalysisType, 
                                                       unsigned int binsVertex, double minVertex, double maxVertex, 
                                                       unsigned int binsMult, double minMult, double maxMult, 
                                                       bool aIsMCRun, bool aImplementAvgSepCuts, bool aWritePairKinematics, TString aDirNameModifier) :

  AliFemtoVertexMultAnalysis(binsVertex,minVertex,maxVertex,binsMult,minMult,maxMult),
  fAnalysisParams(DefaultAnalysisParams()),
  fAnalysisType(aAnalysisType),
  fGeneralAnalysisType(),
  fParticlePDGType1(),
  fParticlePDGType2(),
  fGeneralParticleType1(),
  fGeneralParticleType2(),

  fOutputName(fAnalysisTags[aAnalysisType]),
  fMultHist(NULL),
  fImplementAvgSepCuts(aImplementAvgSepCuts),
  fWritePairKinematics(aWritePairKinematics),
  fIsMCRun(aIsMCRun),
  fIsMBAnalysis(kFALSE),
  fBuildMultHist(kFALSE),
  fMinCent(-1000),
  fMaxCent(1000),

  fCollectionOfCfs(NULL),

  KStarCf(NULL),
  AvgSepCf(NULL),

  KStarModelCfs(NULL)
{
  fAnalysisParams.nBinsVertex = binsVertex;
  fAnalysisParams.minVertex = minVertex;
  fAnalysisParams.maxVertex = maxVertex;

  fAnalysisParams.nBinsMult = binsMult;
  fAnalysisParams.minMult = minMult;
  fAnalysisParams.maxMult = maxMult;

  fAnalysisParams.analysisType = aAnalysisType;

  fAnalysisParams.isMCRun = aIsMCRun;
  fAnalysisParams.implementVertexCorrections = aImplementAvgSepCuts;
  fAnalysisParams.writePairKinematics = aWritePairKinematics;

  SetParticleTypes(fAnalysisType);
  SetVerboseMode(kFALSE);
  SetNumEventsToMix(5);
  SetMinSizePartCollection(1);
  SetV0SharedDaughterCut(kTRUE);
  SetEnablePairMonitors(fIsMCRun);
  SetMultHist(fAnalysisTags[aAnalysisType]);

  fMinCent = minMult/10.;
  fMaxCent = maxMult/10.;

  fOutputName += TString::Format("_%0.0f%0.0f",fMinCent,fMaxCent);
  if(!aDirNameModifier.IsNull())
  {
    fOutputName += "_";
    fOutputName += aDirNameModifier;
  }


  fCollectionOfCfs = new AliFemtoCorrFctnCollection;

  if(fWritePairKinematics) KStarCf = CreateCorrFctnKStar(fAnalysisTags[aAnalysisType],62,0.,0.31); //TNtuple is huge, and I don't need data out to 1 GeV
  else KStarCf = CreateCorrFctnKStar(fAnalysisTags[aAnalysisType],200,0.,1.0);

  AvgSepCf = CreateAvgSepCorrFctn(fAnalysisTags[aAnalysisType],200,0.,20.);

  KStarModelCfs = CreateModelCorrFctnKStarFull(fAnalysisTags[aAnalysisType],200,0.,1.0);

  if(fWritePairKinematics) fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarCf);
  else
  {
    fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarCf);
    if(fAnalysisType!=kXiKchP && fAnalysisType!=kAXiKchP && fAnalysisType!=kXiKchM && fAnalysisType!=kAXiKchM) fCollectionOfCfs->push_back((AliFemtoCorrFctn*)AvgSepCf);
    if(fAnalysisType==kProtPiM || fAnalysisType==kAProtPiP || fAnalysisType==kPiPPiM) fCollectionOfCfs->push_back((AliFemtoV0PurityBgdEstimator*)CreateV0PurityBgdEstimator());
  }

  if(fIsMCRun) fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarModelCfs);

}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, TString aDirNameModifier) :
  AliFemtoVertexMultAnalysis(aAnParams.nBinsVertex,aAnParams.minVertex,aAnParams.maxVertex,
                             aAnParams.nBinsMult,aAnParams.minMult,aAnParams.maxMult),
  fAnalysisParams(aAnParams),
  fAnalysisType(aAnParams.analysisType),
  fGeneralAnalysisType(),
  fParticlePDGType1(),
  fParticlePDGType2(),
  fGeneralParticleType1(),
  fGeneralParticleType2(),
  fOutputName(fAnalysisTags[fAnalysisType]),
  fMultHist(NULL),
  fImplementAvgSepCuts(aAnParams.implementAvgSepCuts),
  fWritePairKinematics(aAnParams.writePairKinematics),
  fIsMCRun(aAnParams.isMCRun),
  fIsMBAnalysis(aAnParams.isMBAnalysis),
  fBuildMultHist(aAnParams.buildMultHist),
  fMinCent(-1000),
  fMaxCent(1000),

  fCollectionOfCfs(NULL),

  KStarCf(NULL),
  AvgSepCf(NULL),

  KStarModelCfs(NULL)
{
  SetParticleTypes(fAnalysisType);
  SetVerboseMode(aAnParams.verbose);
  SetNumEventsToMix(aAnParams.nEventsToMix);
  SetMinSizePartCollection(aAnParams.minCollectionSize);
  SetV0SharedDaughterCut(aAnParams.setV0SharedDaughterCut);
  SetEnablePairMonitors(aAnParams.isMCRun);
  SetMultHist(fAnalysisTags[fAnalysisType]);


  fMinCent = aAnParams.minMult/10.;
  fMaxCent = aAnParams.maxMult/10.;

  fOutputName += TString::Format("_%0.0f%0.0f",fMinCent,fMaxCent);
  if(!aDirNameModifier.IsNull())
  {
    fOutputName += "_";
    fOutputName += aDirNameModifier;
  }

  fCollectionOfCfs = new AliFemtoCorrFctnCollection;


  if(fWritePairKinematics)
  {
    KStarCf = CreateCorrFctnKStar(fAnalysisTags[fAnalysisType],62,0.,0.31); //TNtuple is huge, and I don't need data out to 1 GeV
    fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarCf);
  }
  else
  {
    KStarCf = CreateCorrFctnKStar(fAnalysisTags[fAnalysisType],200,0.,1.0);
    AvgSepCf = CreateAvgSepCorrFctn(fAnalysisTags[fAnalysisType],200,0.,20.);
    fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarCf);
    if(fAnalysisType!=kXiKchP && fAnalysisType!=kAXiKchP && fAnalysisType!=kXiKchM && fAnalysisType!=kAXiKchM) fCollectionOfCfs->push_back((AliFemtoCorrFctn*)AvgSepCf);
    if(fAnalysisType==kProtPiM || fAnalysisType==kAProtPiP || fAnalysisType==kPiPPiM) fCollectionOfCfs->push_back((AliFemtoV0PurityBgdEstimator*)CreateV0PurityBgdEstimator());
  }

  if(fIsMCRun)
  {
    KStarModelCfs = CreateModelCorrFctnKStarFull(fAnalysisTags[fAnalysisType],200,0.,1.0);
    KStarModelCfs->SetRemoveMisidentified(fAnalysisParams.removeMisidentifiedMCParticles);
    fCollectionOfCfs->push_back((AliFemtoCorrFctn*)KStarModelCfs);
  }

}


//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, V0CutParams &aV0CutParams1, V0CutParams &aV0CutParams2, TString aDirNameModifier) :
  AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aDirNameModifier)
{
  AliFemtoV0TrackCutNSigmaFilter* tV0Cut1 = CreateV0Cut(aV0CutParams1);
  AliFemtoV0TrackCutNSigmaFilter* tV0Cut2 = CreateV0Cut(aV0CutParams2);
  AliFemtoV0PairCut* tV0PairCut = CreateV0PairCut(aPairCutParams);

  if(aAnParams.isMBAnalysis)
  {
    AliFemtoBasicEventCut* tBasicEvCut = CreateBasicEventCut(aEvCutParams);
    SetAnalysis(tBasicEvCut,tV0Cut1,tV0Cut2,tV0PairCut);
  }
  else
  {
    AliFemtoEventCutEstimators* tEvCutEst = CreateEventCutEstimators(aEvCutParams);
    SetAnalysis(tEvCutEst,tV0Cut1,tV0Cut2,tV0PairCut);
  }

}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, V0CutParams &aV0CutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier) :
  AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aDirNameModifier)
{
  AliFemtoV0TrackCutNSigmaFilter* tV0Cut1 = CreateV0Cut(aV0CutParams1);
  AliFemtoESDTrackCutNSigmaFilter* tESDCut2 = CreateESDCut(aESDCutParams2);
  AliFemtoV0TrackPairCut* tV0TrackPairCut = CreateV0TrackPairCut(aPairCutParams);

  if(aAnParams.isMBAnalysis)
  {
    AliFemtoBasicEventCut* tBasicEvCut = CreateBasicEventCut(aEvCutParams);
    SetAnalysis(tBasicEvCut,tV0Cut1,tESDCut2,tV0TrackPairCut);
  }
  else
  {
    AliFemtoEventCutEstimators* tEvCutEst = CreateEventCutEstimators(aEvCutParams);
    SetAnalysis(tEvCutEst,tV0Cut1,tESDCut2,tV0TrackPairCut);
  }
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, XiCutParams &aXiCutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier) :
  AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aDirNameModifier)
{
  AliFemtoXiTrackCutNSigmaFilter* tXiCut1 = CreateXiCut(aXiCutParams1);
  AliFemtoESDTrackCutNSigmaFilter* tESDCut2 = CreateESDCut(aESDCutParams2);
  AliFemtoXiTrackPairCut* tXiTrackPairCut = CreateXiTrackPairCut(aPairCutParams);

  if(aAnParams.isMBAnalysis)
  {
    AliFemtoBasicEventCut* tBasicEvCut = CreateBasicEventCut(aEvCutParams);
    SetAnalysis(tBasicEvCut,tXiCut1,tESDCut2,tXiTrackPairCut);
  }
  else
  {
    AliFemtoEventCutEstimators* tEvCutEst = CreateEventCutEstimators(aEvCutParams);
    SetAnalysis(tEvCutEst,tXiCut1,tESDCut2,tXiTrackPairCut);
  }
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(AnalysisParams &aAnParams, EventCutParams &aEvCutParams, PairCutParams &aPairCutParams, ESDCutParams &aESDCutParams1, ESDCutParams &aESDCutParams2, TString aDirNameModifier) :
  AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aDirNameModifier)
{
  AliFemtoESDTrackCutNSigmaFilter* tESDCut1 = CreateESDCut(aESDCutParams1);
  AliFemtoESDTrackCutNSigmaFilter* tESDCut2 = CreateESDCut(aESDCutParams2);
  AliFemtoDummyPairCut *tDummyPairCut = new AliFemtoDummyPairCut();

  if(aAnParams.isMBAnalysis)
  {
    AliFemtoBasicEventCut* tBasicEvCut = CreateBasicEventCut(aEvCutParams);
    SetAnalysis(tBasicEvCut,tESDCut1,tESDCut2,tDummyPairCut);
  }
  else
  {
    AliFemtoEventCutEstimators* tEvCutEst = CreateEventCutEstimators(aEvCutParams);
    SetAnalysis(tEvCutEst,tESDCut1,tESDCut2,tDummyPairCut);
  }
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AliFemtoAnalysisLambdaKaon(const AliFemtoAnalysisLambdaKaon& a) :
  AliFemtoVertexMultAnalysis(a),

  fAnalysisParams(a.fAnalysisParams),
  fAnalysisType(a.fAnalysisType),
  fGeneralAnalysisType(a.fGeneralAnalysisType),
  fParticlePDGType1(a.fParticlePDGType1),
  fParticlePDGType2(a.fParticlePDGType2),
  fGeneralParticleType1(a.fGeneralParticleType1),
  fGeneralParticleType2(a.fGeneralParticleType2),
  fOutputName(a.fOutputName),
  fMultHist(NULL),
  fImplementAvgSepCuts(a.fImplementAvgSepCuts),
  fWritePairKinematics(a.fWritePairKinematics),
  fIsMCRun(a.fIsMCRun),
  fIsMBAnalysis(a.fIsMBAnalysis),
  fBuildMultHist(a.fBuildMultHist),
  fMinCent(a.fMinCent),
  fMaxCent(a.fMaxCent),

  fCollectionOfCfs(NULL),

  KStarCf(NULL),
  AvgSepCf(NULL),

  KStarModelCfs(NULL)
{

  if(a.fMultHist) fMultHist = new TH1F(*(a.fMultHist));
  else fMultHist = NULL;

  fCollectionOfCfs = new AliFemtoCorrFctnCollection;
  AliFemtoCorrFctnIterator iter;
  for(iter=a.fCollectionOfCfs->begin(); iter!=a.fCollectionOfCfs->end(); iter++)
  {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if(fctn) {AddCorrFctn(fctn);}
    else {cout << " myTrainAnalysisConstructor::myTrainAnalysisConstructor(const myTrainAnalysisConstructor& a) - correlation function not found " << endl;}
  }

  if(a.KStarCf) KStarCf = new AliFemtoCorrFctnKStar(*(a.KStarCf));
  else KStarCf = NULL;

  if(a.AvgSepCf) AvgSepCf = new AliFemtoAvgSepCorrFctn(*(a.AvgSepCf));
  else AvgSepCf = NULL;

  if(a.KStarModelCfs) KStarModelCfs = new AliFemtoModelCorrFctnKStarFull(*(a.KStarModelCfs));
  else KStarModelCfs = NULL;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon& AliFemtoAnalysisLambdaKaon::operator=(const AliFemtoAnalysisLambdaKaon& a)
{
  if(this == &a) {return *this;}

  AliFemtoVertexMultAnalysis::operator=(a);


  fAnalysisParams = a.fAnalysisParams;
  fAnalysisType = a.fAnalysisType;
  fGeneralAnalysisType = a.fGeneralAnalysisType;
  fParticlePDGType1 = a.fParticlePDGType1;
  fParticlePDGType2 = a.fParticlePDGType2;
  fGeneralParticleType1 = a.fGeneralParticleType1;
  fGeneralParticleType2 = a.fGeneralParticleType2;
  fOutputName = a.fOutputName;
  fMultHist = NULL;
  fImplementAvgSepCuts = a.fImplementAvgSepCuts;
  fWritePairKinematics = a.fWritePairKinematics;
  fIsMCRun = a.fIsMCRun;
  fIsMBAnalysis = a.fIsMBAnalysis;
  fBuildMultHist = a.fBuildMultHist;
  fMinCent = a.fMinCent;
  fMaxCent = a.fMaxCent;
  fCollectionOfCfs = NULL;
  KStarCf = NULL;
  AvgSepCf = NULL;

  KStarModelCfs = NULL;

  if(a.fMultHist) fMultHist = new TH1F(*(a.fMultHist));
  else fMultHist = NULL;

  fCollectionOfCfs = new AliFemtoCorrFctnCollection;
  AliFemtoCorrFctnIterator iter;
  for(iter=a.fCollectionOfCfs->begin(); iter!=a.fCollectionOfCfs->end(); iter++)
  {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if(fctn) {AddCorrFctn(fctn);}
    else {cout << " myTrainAnalysisConstructor::myTrainAnalysisConstructor(const myTrainAnalysisConstructor& a) - correlation function not found " << endl;}
  }

  if(a.KStarCf) KStarCf = new AliFemtoCorrFctnKStar(*(a.KStarCf));
  else KStarCf = NULL;

  if(a.AvgSepCf) AvgSepCf = new AliFemtoAvgSepCorrFctn(*(a.AvgSepCf));
  else AvgSepCf = NULL;

  if(a.KStarModelCfs) KStarModelCfs = new AliFemtoModelCorrFctnKStarFull(*(a.KStarModelCfs));
  else KStarModelCfs = NULL;

  return *this;
}


//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::~AliFemtoAnalysisLambdaKaon()
{
//no-op
}

//____________________________
void AliFemtoAnalysisLambdaKaon::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  double multiplicity = hbtEvent->UncorrectedNumberOfPrimaries();
  if(fBuildMultHist) fMultHist->Fill(multiplicity);
  AliFemtoVertexMultAnalysis::ProcessEvent(hbtEvent);
}

//____________________________
TList* AliFemtoAnalysisLambdaKaon::GetOutputList()
{
  TList *olist = new TList();
  TObjArray *temp = new TObjArray();
  olist->SetName(fOutputName.Data());
  temp->SetName(fOutputName.Data());

  TList *tOutputList = AliFemtoSimpleAnalysis::GetOutputList(); 

  if(fBuildMultHist) tOutputList->Add(fMultHist);

  TListIter next(tOutputList);
  while (TObject *obj = next())
  {
    temp->Add(obj);
  }

  olist->Add(temp);    
  return olist;
}


//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::SetParticleTypes(AliFemtoAnalysisLambdaKaon::AnalysisType aAnType)
{
  switch(aAnType) {
  case AliFemtoAnalysisLambdaKaon::kLamK0:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 =AliFemtoAnalysisLambdaKaon:: kPDGK0;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamK0:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGK0;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamKchP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchP;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamKchP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchP;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamKchM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchM;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamKchM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchM;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamLam:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamALam:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamALam:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamPiP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiP;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamPiP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiP;
    break;

  case AliFemtoAnalysisLambdaKaon::kLamPiM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGLam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiM;
    break;

  case AliFemtoAnalysisLambdaKaon::kALamPiM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kV0Track;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGALam;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiM;
    break;

  case AliFemtoAnalysisLambdaKaon::kXiKchP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kXiTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGXiC;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchP;
    break;

  case AliFemtoAnalysisLambdaKaon::kAXiKchP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kXiTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGAXiC;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchP;
    break;

  case AliFemtoAnalysisLambdaKaon::kXiKchM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kXiTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGXiC;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchM;
    break;

  case AliFemtoAnalysisLambdaKaon::kAXiKchM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kXiTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGAXiC;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGKchM;
    break;

  case AliFemtoAnalysisLambdaKaon::kProtPiM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kTrackTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGProt;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiM;
    break;

  case AliFemtoAnalysisLambdaKaon::kAProtPiP:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kTrackTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGAntiProt;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiP;
    break;

  case AliFemtoAnalysisLambdaKaon::kPiPPiM:
    fGeneralAnalysisType = AliFemtoAnalysisLambdaKaon::kTrackTrack;
    fParticlePDGType1 = AliFemtoAnalysisLambdaKaon::kPDGPiP;
    fParticlePDGType2 = AliFemtoAnalysisLambdaKaon::kPDGPiM;
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::SetParticleTypes: Invalid AnalysisType"
            "selection '" << aAnType << endl;
  }

  switch(fGeneralAnalysisType) {
  case AliFemtoAnalysisLambdaKaon::kV0V0:
    fGeneralParticleType1 = AliFemtoAnalysisLambdaKaon::kV0;
    fGeneralParticleType2 = AliFemtoAnalysisLambdaKaon::kV0;
    break;

  case AliFemtoAnalysisLambdaKaon::kV0Track:
    fGeneralParticleType1 = AliFemtoAnalysisLambdaKaon::kV0;
    fGeneralParticleType2 = AliFemtoAnalysisLambdaKaon::kTrack;
    break;

  case AliFemtoAnalysisLambdaKaon::kXiTrack:
    fGeneralParticleType1 = AliFemtoAnalysisLambdaKaon::kCascade;
    fGeneralParticleType2 = AliFemtoAnalysisLambdaKaon::kTrack;
    break;

  case AliFemtoAnalysisLambdaKaon::kTrackTrack:
    fGeneralParticleType1 = AliFemtoAnalysisLambdaKaon::kTrack;
    fGeneralParticleType2 = AliFemtoAnalysisLambdaKaon::kTrack;
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::SetParticleTypes" << endl;
  }

}


//___________________________________________________________________
AliFemtoBasicEventCut* AliFemtoAnalysisLambdaKaon::CreateBasicEventCut(EventCutParams &aEvCutParams)
{
  AliFemtoBasicEventCut* mec = new AliFemtoBasicEventCut();
    //Accept events with the given multiplicity
    mec->SetEventMult(aEvCutParams.minMult,aEvCutParams.maxMult);
    //and z-vertex distance to the center of the TPC
    mec->SetVertZPos(aEvCutParams.minVertexZ,aEvCutParams.maxVertexZ);
  return mec;
}



//___________________________________________________________________
AliFemtoEventCutEstimators* AliFemtoAnalysisLambdaKaon::CreateEventCutEstimators(EventCutParams &aEvCutParams)
{
  AliFemtoEventCutEstimators* EvCutEst = new AliFemtoEventCutEstimators();
    EvCutEst->SetCentEst1Range(aEvCutParams.minCentrality,aEvCutParams.maxCentrality);
    EvCutEst->SetVertZPos(aEvCutParams.minVertexZ,aEvCutParams.maxVertexZ);
    EvCutEst->SetVerboseMode(aEvCutParams.verboseMode);
  return EvCutEst;
}




//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomV0SelectionFilters(ParticlePDGType aV0Type, AliFemtoV0TrackCutNSigmaFilter* aCut)
{
  switch(aV0Type) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    //--Proton(+) daughter selection filter
    aCut->CreateCustomProtonNSigmaFilter();
    aCut->AddProtonTPCNSigmaCut(0.,0.8,3.);
    aCut->AddProtonTPCAndTOFNSigmaCut(0.8,1000.,3.,3.);
    aCut->AddProtonTPCNSigmaCut(0.8,1000.,3.);

    //--Pion(-) daughter selection filter
/*
    //the standard cuts in AliFemtoV0TrackCut
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,1000.,3.);
*/

    //RequireTOFPion
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,0.5,3.);
    aCut->AddPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    aCut->AddPionTPCNSigmaCut(0.5,1000.,3.);

    break;

  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    //--(Anti)Proton(-) daughter selection filter
    //for now, the custom filters will match the standard cuts in AliFemtoV0TrackCut
    //these also match my personal (proton) cuts in myAliFemtoV0TrackCut
    aCut->CreateCustomProtonNSigmaFilter();
    aCut->AddProtonTPCNSigmaCut(0.,0.8,3.);
    aCut->AddProtonTPCAndTOFNSigmaCut(0.8,1000.,3.,3.);
    aCut->AddProtonTPCNSigmaCut(0.8,1000.,3.);

/*
    //the standard cuts in AliFemtoV0TrackCut
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,1000.,3.);
*/


    //RequireTOFPion
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,0.5,3.);
    aCut->AddPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    aCut->AddPionTPCNSigmaCut(0.5,1000.,3.);

    break;

  case AliFemtoAnalysisLambdaKaon::kPDGK0:
    //--Pion(+) daughter selection filter
/*
    //the standard cuts in AliFemtoV0TrackCut
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,1000.,3.);
*/


    //RequireTOFPion
    aCut->CreateCustomPionNSigmaFilter();
    aCut->AddPionTPCNSigmaCut(0.,0.5,3.);
    aCut->AddPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    aCut->AddPionTPCNSigmaCut(0.5,1000.,3.);

    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomV0SelectionFilters" << endl;
  }
}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomV0RejectionFilters(ParticlePDGType aV0Type, AliFemtoV0TrackCutNSigmaFilter* aCut)
{
  switch(aV0Type) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    aCut->CreateCustomV0Rejection(AliFemtoV0TrackCut::kK0s);
    aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.,0.5,3.,  //positive daughter
                                         0.,0.5,3.); //negative daughter
    aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                               0.5,1000.,3.,3.,  //positive daughter
                                               0.5,1000.,3.,3.); //negative daughter
    aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.5,1000.,3.,  //positive daughter
                                         0.5,1000.,3.); //negative daughter
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    aCut->CreateCustomV0Rejection(AliFemtoV0TrackCut::kK0s);
    aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.,0.5,3.,  //positive daughter
                                         0.,0.5,3.); //negative daughter
    aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                               0.5,1000.,3.,3.,  //positive daughter
                                               0.5,1000.,3.,3.); //negative daughter
    aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.5,1000.,3.,  //positive daughter
                                         0.5,1000.,3.); //negative daughter
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGK0:
    //Lambda rejection
    aCut->CreateCustomV0Rejection(AliFemtoV0TrackCut::kLambda);
      //Positive daughter (Proton)
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,1,0.,0.8,3.);
      aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,1,0.8,1000.,3.,3.);
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,1,0.8,1000.,3.);
      //Negative daughter (Pion)
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,-1,0.,0.5,3.);
      aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,-1,0.5,1000.,3.,3.);
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kLambda,-1,0.5,1000.,3.);

    //AntiLambda rejection
    aCut->CreateCustomV0Rejection(AliFemtoV0TrackCut::kAntiLambda);
      //Positive daughter (Pion)
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,1,0.,0.5,3.);
      aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,1,0.5,1000.,3.,3.);
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,1,0.5,1000.,3.);
      //Negative daughter (AntiProton)
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,-1,0.,0.8,3.);
      aCut->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,-1,0.8,1000.,3.,3.);
      aCut->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kAntiLambda,-1,0.8,1000.,3.);
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomV0RejectionFilters" << endl;
  }
}

//___________________________________________________________________
AliFemtoV0TrackCutNSigmaFilter* AliFemtoAnalysisLambdaKaon::CreateV0Cut(V0CutParams &aCutParams)
{
  AliFemtoV0TrackCutNSigmaFilter* tV0Cut = new AliFemtoV0TrackCutNSigmaFilter();

  tV0Cut->SetParticleType(aCutParams.v0Type);
  tV0Cut->SetMass(aCutParams.mass);
  tV0Cut->SetInvariantMassLambda(aCutParams.minInvariantMass,aCutParams.maxInvariantMass);
  tV0Cut->SetLooseInvMassCut(aCutParams.useLooseInvMassCut, aCutParams.minLooseInvMass, aCutParams.maxLooseInvMass);

  tV0Cut->SetEta(aCutParams.eta);
  tV0Cut->SetPt(aCutParams.minPt, aCutParams.maxPt);
  tV0Cut->SetOnFlyStatus(aCutParams.onFlyStatus);
  tV0Cut->SetMaxDcaV0(aCutParams.maxDcaV0);
  tV0Cut->SetMinCosPointingAngle(aCutParams.minCosPointingAngle);
  tV0Cut->SetMaxV0DecayLength(aCutParams.maxV0DecayLength);
  //-----
  tV0Cut->SetEtaDaughters(aCutParams.etaDaughters);
  tV0Cut->SetPtPosDaughter(aCutParams.minPtPosDaughter,aCutParams.maxPtPosDaughter);
  tV0Cut->SetPtNegDaughter(aCutParams.minPtNegDaughter,aCutParams.maxPtNegDaughter);
  tV0Cut->SetTPCnclsDaughters(aCutParams.minTPCnclsDaughters);
  //tV0Cut->SetNdofDaughters(4.0); //4.0
  tV0Cut->SetStatusDaughters(AliESDtrack::kTPCrefit/* | AliESDtrack::kITSrefit*/);
  if(fAnalysisType==kProtPiM || fAnalysisType==kAProtPiP || fAnalysisType==kPiPPiM) tV0Cut->SetStatusDaughters(AliESDtrack::kTPCpid);  //TODO currently a work around
  tV0Cut->SetMaxDcaV0Daughters(aCutParams.maxDcaV0Daughters);
  tV0Cut->SetMinDaughtersToPrimVertex(aCutParams.minPosDaughterToPrimVertex,aCutParams.minNegDaughterToPrimVertex);

  if(aCutParams.useCustomFilter) AddCustomV0SelectionFilters(aCutParams.particlePDGType,tV0Cut);

  //Misidentification cuts -----*****-----*****-----*****-----*****-----*****-----*****
  tV0Cut->SetRemoveMisidentified(aCutParams.removeMisID);
  tV0Cut->SetUseSimpleMisIDCut(aCutParams.useSimpleMisID);
  if(!aCutParams.useSimpleMisID && aCutParams.useCustomMisID) AddCustomV0RejectionFilters(aCutParams.particlePDGType,tV0Cut);
  tV0Cut->SetBuildMisIDHistograms(aCutParams.buildMisIDHistograms);

  TString tTitle, tName;
  switch(aCutParams.particlePDGType) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    tName = TString("LambdaPurityAid");
    tTitle = TString("LambdaMinvBeforeFinalCut");

    tV0Cut->SetInvMassReject(AliFemtoV0TrackCut::kK0s, aCutParams.minInvMassReject,aCutParams.maxInvMassReject, aCutParams.removeMisID);

    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kLambda,100,LambdaMass-0.035,LambdaMass+0.035);
    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kK0s,100,K0ShortMass-0.070,K0ShortMass+0.070);
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    tName = TString("AntiLambdaPurityAid");
    tTitle = TString("AntiLambdaMinvBeforeFinalCut");

    tV0Cut->SetInvMassReject(AliFemtoV0TrackCut::kK0s, aCutParams.minInvMassReject,aCutParams.maxInvMassReject, aCutParams.removeMisID);

    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kAntiLambda,100,LambdaMass-0.035,LambdaMass+0.035);
    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kK0s,100,K0ShortMass-0.070,K0ShortMass+0.070);
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGK0:
    tName = TString("K0ShortPurityAid");
    tTitle = TString("K0ShortMinvBeforeFinalCut");

    tV0Cut->SetInvMassReject(AliFemtoV0TrackCut::kLambda, aCutParams.minInvMassReject,aCutParams.maxInvMassReject, aCutParams.removeMisID);
    tV0Cut->SetInvMassReject(AliFemtoV0TrackCut::kAntiLambda, aCutParams.minInvMassReject,aCutParams.maxInvMassReject, aCutParams.removeMisID); 
 
    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kLambda,100,LambdaMass-0.035,LambdaMass+0.035);
    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kAntiLambda,100,LambdaMass-0.035,LambdaMass+0.035);
    tV0Cut->SetMisIDHisto(AliFemtoV0TrackCut::kK0s,100,K0ShortMass-0.070,K0ShortMass+0.070);
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::CreateV0Cut" << endl;
  }

  //-----*****-----*****-----*****-----*****-----*****-----*****-----*****-----*****-----*****

  tV0Cut->SetMinvPurityAidHistoV0(tName,tTitle,aCutParams.nBinsPurity,aCutParams.minPurityMass,aCutParams.maxPurityMass);

  return tV0Cut;
}


//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomESDSelectionFilters(ParticlePDGType aESDType, AliFemtoESDTrackCutNSigmaFilter* aCut)
{
  switch(aESDType) {
  case AliFemtoAnalysisLambdaKaon::kPDGKchP:
  case AliFemtoAnalysisLambdaKaon::kPDGKchM:
    //Kaon filter
/*
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kKaon);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.0,0.5,2.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.5,0.8,3.0,2.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.8,1.0,3.0,1.5);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,1.0,1.5,3.0,1.0);
*/

    //New Kaon filter (Konstantin Mikhaylov)
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kKaon);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.0,0.4,2.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.4,0.45,1.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.45,0.8,3.0,2.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,0.8,1.0,3.0,1.5);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kKaon,1.0,99.,3.0,1.0);

    break;

  case AliFemtoAnalysisLambdaKaon::kPDGPiP:
  case AliFemtoAnalysisLambdaKaon::kPDGPiM:
    //Pion filter
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kPion);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.0,0.65,3.0,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.0,0.35,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.35,0.5,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.5,0.65,2.0);

      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.65,1.5,5.0,3.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,1.5,99.,5.0,2.0);
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGProt:
  case AliFemtoAnalysisLambdaKaon::kPDGAntiProt:
    //Proton filter
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kProton);
    aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kProton,0.,0.8,3.);
    aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kProton,0.8,1000.,3.,3.);
    aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kProton,0.8,1000.,3.);
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomESDSelectionFilters" << endl;
  }
}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomESDRejectionFilters(ParticlePDGType aESDType, AliFemtoESDTrackCutNSigmaFilter* aCut)
{
  switch(aESDType) {
  case AliFemtoAnalysisLambdaKaon::kPDGKchP:
  case AliFemtoAnalysisLambdaKaon::kPDGKchM:
    //Electron filter for removing misidentified
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kElectron);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kElectron,0.,99.,3.0);

    //Pion filter for removing misidentified
    aCut->CreateCustomNSigmaFilter(AliFemtoESDTrackCutNSigmaFilter::kPion);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.0,0.65,3.0,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.0,0.35,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.35,0.5,3.0);
      aCut->AddTPCNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.5,0.65,2.0);

      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,0.65,1.5,5.0,3.0);
      aCut->AddTPCAndTOFNSigmaCut(AliFemtoESDTrackCutNSigmaFilter::kPion,1.5,99.,5.0,2.0);
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomESDSelectionFilters" << endl;
  }
}

//___________________________________________________________________
AliFemtoESDTrackCutNSigmaFilter* AliFemtoAnalysisLambdaKaon::CreateESDCut(ESDCutParams &aCutParams)
{
  AliFemtoESDTrackCutNSigmaFilter* tESDCut = new AliFemtoESDTrackCutNSigmaFilter();

  tESDCut->SetPidProbPion(aCutParams.minPidProbPion,aCutParams.maxPidProbPion);
  tESDCut->SetPidProbMuon(aCutParams.minPidProbMuon,aCutParams.maxPidProbMuon);
  tESDCut->SetPidProbKaon(aCutParams.minPidProbKaon,aCutParams.maxPidProbKaon);
  tESDCut->SetPidProbProton(aCutParams.minPidProbProton,aCutParams.maxPidProbProton);
  tESDCut->SetMostProbable(aCutParams.mostProbable);  //this uses P().Mag() as first argument to IsKaonNSigma()
  tESDCut->SetCharge(aCutParams.charge);
  // so we set the correct mass
  tESDCut->SetMass(aCutParams.mass);
  // we select low pt
  tESDCut->SetPt(aCutParams.minPt,aCutParams.maxPt);
  tESDCut->SetEta(-aCutParams.eta,aCutParams.eta);
//  tESDCut->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);  //This cuts out all particles when used in conjunction with SetFilterBit(7)
  tESDCut->SetminTPCncls(aCutParams.minTPCncls);
  tESDCut->SetRemoveKinks(aCutParams.removeKinks);
  tESDCut->SetLabel(aCutParams.setLabel);
  tESDCut->SetMaxITSChiNdof(aCutParams.maxITSChiNdof);
  tESDCut->SetMaxTPCChiNdof(aCutParams.maxTPCChiNdof);
  tESDCut->SetMaxSigmaToVertex(aCutParams.maxSigmaToVertex);
  tESDCut->SetMinImpactXY(aCutParams.minImpactXY);
  tESDCut->SetMaxImpactXY(aCutParams.maxImpactXY);
  tESDCut->SetMaxImpactZ(aCutParams.maxImpactZ);

  tESDCut->SetElectronRejection(aCutParams.useElectronRejection);
  tESDCut->SetPionRejection(aCutParams.usePionRejection);

  if(aCutParams.useCustomFilter) AddCustomESDSelectionFilters(aCutParams.particlePDGType,tESDCut);
  if(aCutParams.useCustomMisID) AddCustomESDRejectionFilters(aCutParams.particlePDGType,tESDCut);

  tESDCut->SetUseCustomElectronRejection(aCutParams.useCustomElectronRejection);

  return tESDCut;
}



//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomXiSelectionFilters(ParticlePDGType aXiType, AliFemtoXiTrackCutNSigmaFilter* aCut)
{
  bool tDaughterV0FilterExists = false;
  AliFemtoV0TrackCutNSigmaFilter *tDaughterV0Filter;
  if(aCut->GetDaughterV0Filter() != NULL) {tDaughterV0FilterExists = true; tDaughterV0Filter = aCut->GetDaughterV0Filter();}
  else tDaughterV0Filter = new AliFemtoV0TrackCutNSigmaFilter();

  switch(aXiType) {
  case AliFemtoAnalysisLambdaKaon::kPDGXiC:
    //--Proton(+) daughter of Lambda daughter selection filter
    tDaughterV0Filter->CreateCustomProtonNSigmaFilter();
    tDaughterV0Filter->AddProtonTPCNSigmaCut(0.,0.8,3.);
    tDaughterV0Filter->AddProtonTPCAndTOFNSigmaCut(0.8,1000.,3.,3.);
    tDaughterV0Filter->AddProtonTPCNSigmaCut(0.8,1000.,3.);

    //--Pion(-) daughter of Lambda daughter selection filter
    //RequireTOFPion
    tDaughterV0Filter->CreateCustomPionNSigmaFilter();
    tDaughterV0Filter->AddPionTPCNSigmaCut(0.,0.5,3.);
    tDaughterV0Filter->AddPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    tDaughterV0Filter->AddPionTPCNSigmaCut(0.5,1000.,3.);

    //--Bachelor pion daughter
    aCut->CreateCustomBacPionNSigmaFilter();
    aCut->AddBacPionTPCNSigmaCut(0.,0.5,3.);
    aCut->AddBacPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    aCut->AddBacPionTPCNSigmaCut(0.5,1000.,3.);

    break;

  case AliFemtoAnalysisLambdaKaon::kPDGAXiC:
    //--(Anti)Proton(-) daughter of AntiLambda daughter selection filter
    tDaughterV0Filter->CreateCustomProtonNSigmaFilter();
    tDaughterV0Filter->AddProtonTPCNSigmaCut(0.,0.8,3.);
    tDaughterV0Filter->AddProtonTPCAndTOFNSigmaCut(0.8,1000.,3.,3.);
    tDaughterV0Filter->AddProtonTPCNSigmaCut(0.8,1000.,3.);

    //--Pion(-) daughter of AntiLambda daughter selection filter
    //RequireTOFPion
    tDaughterV0Filter->CreateCustomPionNSigmaFilter();
    tDaughterV0Filter->AddPionTPCNSigmaCut(0.,0.5,3.);
    tDaughterV0Filter->AddPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    tDaughterV0Filter->AddPionTPCNSigmaCut(0.5,1000.,3.);

    //--Bachelor pion daughter
    aCut->CreateCustomBacPionNSigmaFilter();
    aCut->AddBacPionTPCNSigmaCut(0.,0.5,3.);
    aCut->AddBacPionTPCAndTOFNSigmaCut(0.5,1000.,3.,3.);
    aCut->AddBacPionTPCNSigmaCut(0.5,1000.,3.);

    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomXiSelectionFilters" << endl;
  }
  if(!tDaughterV0FilterExists) aCut->SetDaughterV0Filter(new AliFemtoV0TrackCutNSigmaFilter(*tDaughterV0Filter));

}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCustomXiV0RejectionFilters(ParticlePDGType aXiType, AliFemtoXiTrackCutNSigmaFilter* aCut)
{
  bool tDaughterV0FilterExists = false;
  AliFemtoV0TrackCutNSigmaFilter *tDaughterV0Filter;
  if(aCut->GetDaughterV0Filter() != NULL) {tDaughterV0FilterExists = true; tDaughterV0Filter = aCut->GetDaughterV0Filter();}
  else tDaughterV0Filter = new AliFemtoV0TrackCutNSigmaFilter();

  ParticlePDGType tV0Type;
  if(aXiType==AliFemtoAnalysisLambdaKaon::kPDGXiC) tV0Type = AliFemtoAnalysisLambdaKaon::kPDGLam;
  else if(aXiType==AliFemtoAnalysisLambdaKaon::kPDGAXiC) tV0Type = AliFemtoAnalysisLambdaKaon::kPDGALam;
  else tV0Type = kPDGNull;

  switch(tV0Type) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    tDaughterV0Filter->CreateCustomV0Rejection(AliFemtoV0TrackCut::kK0s);
    tDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.,0.5,3.,  //positive daughter
                                         0.,0.5,3.); //negative daughter
    tDaughterV0Filter->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                               0.5,1000.,3.,3.,  //positive daughter
                                               0.5,1000.,3.,3.); //negative daughter
    tDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.5,1000.,3.,  //positive daughter
                                         0.5,1000.,3.); //negative daughter
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    tDaughterV0Filter->CreateCustomV0Rejection(AliFemtoV0TrackCut::kK0s);
    tDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.,0.5,3.,  //positive daughter
                                         0.,0.5,3.); //negative daughter
    tDaughterV0Filter->AddTPCAndTOFNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                               0.5,1000.,3.,3.,  //positive daughter
                                               0.5,1000.,3.,3.); //negative daughter
    tDaughterV0Filter->AddTPCNSigmaCutToV0Rejection(AliFemtoV0TrackCut::kK0s,
                                         0.5,1000.,3.,  //positive daughter
                                         0.5,1000.,3.); //negative daughter
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCustomV0RejectionFilters" << endl;
  }
  if(!tDaughterV0FilterExists) aCut->SetDaughterV0Filter(new AliFemtoV0TrackCutNSigmaFilter(*tDaughterV0Filter));
}



//___________________________________________________________________
AliFemtoXiTrackCutNSigmaFilter* AliFemtoAnalysisLambdaKaon::CreateXiCut(XiCutParams &aCutParams)
{
  //NOTE: the SetMass call actually is important
  //      This should be set to the mass of the particle of interest, here the Xi
  //      Be sure to not accidentally set it again in the Lambda cuts (for instance, when copy/pasting the lambda cuts from above!)

  AliFemtoXiTrackCutNSigmaFilter* tXiCut = new AliFemtoXiTrackCutNSigmaFilter();

  //Xi Cuts
  tXiCut->SetChargeXi(aCutParams.charge);
  tXiCut->SetParticleTypeXi(aCutParams.xiType);  //kXiMinus = 0, kXiPlus = 1
  tXiCut->SetPtXi(aCutParams.minPt,aCutParams.maxPt);
  tXiCut->SetEtaXi(aCutParams.eta);
  tXiCut->SetMass(aCutParams.mass);
  tXiCut->SetInvariantMassXi(aCutParams.minInvariantMass,aCutParams.maxInvariantMass);
  tXiCut->SetMaxDecayLengthXi(aCutParams.maxDecayLengthXi);
  tXiCut->SetMinCosPointingAngleXi(aCutParams.minCosPointingAngleXi);
  tXiCut->SetMinCosPointingAngleV0toXi(aCutParams.minCosPointingAngleV0toXi);
  tXiCut->SetMaxDcaXi(aCutParams.maxDcaXi);
    //XiDaughters
    tXiCut->SetMaxDcaXiDaughters(aCutParams.maxDcaXiDaughters);


  //Bachelor cuts (here = PiM)
  tXiCut->SetMinDcaXiBac(aCutParams.minDcaXiBac);
  tXiCut->SetEtaBac(aCutParams.etaBac);
  tXiCut->SetTPCnclsBac(aCutParams.minTPCnclsBac);
  tXiCut->SetPtBac(aCutParams.minPtBac,aCutParams.maxPtBac);
  tXiCut->SetStatusBac(AliESDtrack::kTPCrefit);  //yes or no?


  //Lambda cuts (regular V0)
  tXiCut->SetParticleType(aCutParams.v0Type); //0=lambda
  tXiCut->SetMinDcaV0(aCutParams.minDcaV0);
  tXiCut->SetInvariantMassLambda(aCutParams.minInvMassV0,aCutParams.maxInvMassV0);
  tXiCut->SetMinCosPointingAngle(aCutParams.minCosPointingAngleV0);
  tXiCut->SetEta(aCutParams.etaV0);
  tXiCut->SetPt(aCutParams.minPtV0,aCutParams.maxPtV0);
  tXiCut->SetOnFlyStatus(aCutParams.onFlyStatusV0);
  tXiCut->SetMaxV0DecayLength(aCutParams.maxV0DecayLength);
    //Lambda daughter cuts
    tXiCut->SetMinDaughtersToPrimVertex(aCutParams.minV0PosDaughterToPrimVertex,aCutParams.minV0NegDaughterToPrimVertex);
    tXiCut->SetMaxDcaV0Daughters(aCutParams.maxDcaV0Daughters);
    tXiCut->SetEtaDaughters(aCutParams.etaV0Daughters);
    tXiCut->SetPtPosDaughter(aCutParams.minPtPosV0Daughter,aCutParams.maxPtPosV0Daughter); //0.5 for protons
    tXiCut->SetPtNegDaughter(aCutParams.minPtNegV0Daughter,aCutParams.maxPtNegV0Daughter); //0.16 for pions
    tXiCut->SetTPCnclsDaughters(aCutParams.minTPCnclsV0Daughters);
    tXiCut->SetStatusDaughters(AliESDtrack::kTPCrefit);  //yes or no?

  if(aCutParams.useCustomV0Filter) AddCustomXiSelectionFilters(aCutParams.particlePDGType,tXiCut);
  if(aCutParams.useCustomV0MisID) AddCustomXiV0RejectionFilters(aCutParams.particlePDGType,tXiCut);

  TString tTitleXi, tNameXi;
  TString tTitleLam, tNameLam;
  switch(aCutParams.particlePDGType) {
  case AliFemtoAnalysisLambdaKaon::kPDGXiC:
    tNameXi = TString("XiPurityAid");
    tTitleXi = TString("XiMinvBeforeFinalCut");

    tNameLam = TString("LambdaPurityAid");
    tTitleLam = TString("LambdaMinvBeforeFinalCut");
    break;

  case AliFemtoAnalysisLambdaKaon::kPDGAXiC:
    tNameXi = TString("AXiPurityAid");
    tTitleXi = TString("AXiMinvBeforeFinalCut");

    tNameLam = TString("AntiLambdaPurityAid");
    tTitleLam = TString("AntiLambdaMinvBeforeFinalCut");
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::CreateXiCut" << endl;
  }

  tXiCut->SetMinvPurityAidHistoXi(tNameXi,tTitleXi,100,XiMass-0.035,XiMass+0.035);
  tXiCut->SetMinvPurityAidHistoV0(tNameLam,tTitleLam,100,LambdaMass-0.035,LambdaMass+0.035);

  return tXiCut;
}


//___________________________________________________________________
AliFemtoV0PairCut* AliFemtoAnalysisLambdaKaon::CreateV0PairCut(PairCutParams &aPairCutParams)
{
  AliFemtoV0PairCut *v0pc1 = new AliFemtoV0PairCut();  //K0Short-lambda
//  v0pc1->SetV0Max(0.25);
//  v0pc1->SetShareFractionMax(0.05)  //how do I implement this in AliFemtoV0PairCut?
  v0pc1->SetRemoveSameLabel(aPairCutParams.removeSameLabel);
  v0pc1->SetTPCExitSepMinimum(aPairCutParams.tpcExitSepMinimum);  //Default is 0, but for some reason distExitPos(Neg) always end up as 0?

  double tPosPos=0., tPosNeg=0., tNegPos=0., tNegNeg=0.;
  if(fImplementAvgSepCuts)
  {
    switch(fAnalysisType) {
    case AliFemtoAnalysisLambdaKaon::kLamK0:
    case AliFemtoAnalysisLambdaKaon::kLamLam:
      tPosPos = aPairCutParams.minAvgSepPosPos;
      tPosNeg = aPairCutParams.minAvgSepPosNeg;
      tNegPos = aPairCutParams.minAvgSepNegPos;
      tNegNeg = aPairCutParams.minAvgSepNegNeg;
      break;

    case AliFemtoAnalysisLambdaKaon::kALamK0:
    case AliFemtoAnalysisLambdaKaon::kALamALam:
      tPosPos = aPairCutParams.minAvgSepNegNeg;
      tPosNeg = aPairCutParams.minAvgSepNegPos;
      tNegPos = aPairCutParams.minAvgSepPosNeg;
      tNegNeg = aPairCutParams.minAvgSepPosPos;
      break;

    case AliFemtoAnalysisLambdaKaon::kLamALam:  //Jai (typically PosNeg=NegPos=0)
      tPosPos = 0.5*(aPairCutParams.minAvgSepPosPos + aPairCutParams.minAvgSepNegNeg);
      tPosNeg = 0.5*(aPairCutParams.minAvgSepPosNeg + aPairCutParams.minAvgSepNegPos);
      tNegPos = 0.5*(aPairCutParams.minAvgSepPosNeg + aPairCutParams.minAvgSepNegPos);
      tNegNeg = 0.5*(aPairCutParams.minAvgSepPosPos + aPairCutParams.minAvgSepNegNeg);
      break;

    default:
      break;
    }
  }

  v0pc1->SetMinAvgSeparation(0,tPosPos);
  v0pc1->SetMinAvgSeparation(1,tPosNeg);
  v0pc1->SetMinAvgSeparation(2,tNegPos);
  v0pc1->SetMinAvgSeparation(3,tNegNeg);

  return v0pc1;
}

//___________________________________________________________________
AliFemtoV0TrackPairCut* AliFemtoAnalysisLambdaKaon::CreateV0TrackPairCut(PairCutParams &aPairCutParams)
{
  AliFemtoV0TrackPairCut *v0TrackPairCut1 = new AliFemtoV0TrackPairCut();
    v0TrackPairCut1->SetShareQualityMax(aPairCutParams.shareQualityMax);
    v0TrackPairCut1->SetShareFractionMax(aPairCutParams.shareFractionMax);
    v0TrackPairCut1->SetTPCOnly(aPairCutParams.tpcOnly);
    v0TrackPairCut1->SetDataType(AliFemtoPairCut::kAOD);
    v0TrackPairCut1->SetTPCEntranceSepMinimum(aPairCutParams.tpcEntranceSepMinimum);
    v0TrackPairCut1->SetTPCExitSepMinimum(aPairCutParams.tpcExitSepMinimum);
//    v0TrackPairCut1->SetKstarCut(0.04,AliFemtoV0TrackPairCut::kAntiLambda,AliFemtoV0TrackPairCut::kAntiProton); //1 - antilambda, 3 - antiproton
//    v0TrackPairCut1->SetMinAvgSeparation(0,0); //0 - track-pos, 1 - track-neg
//    v0TrackPairCut1->SetMinAvgSeparation(1,11);
    v0TrackPairCut1->SetRemoveSameLabel(aPairCutParams.removeSameLabel);


  double tTrackPos=0., tTrackNeg=0.;
  if(fImplementAvgSepCuts)
  {
    switch(fAnalysisType) {
    case AliFemtoAnalysisLambdaKaon::kLamKchP:
    case AliFemtoAnalysisLambdaKaon::kALamKchP:
    case AliFemtoAnalysisLambdaKaon::kLamPiP:
    case AliFemtoAnalysisLambdaKaon::kALamPiP:
      tTrackPos = aPairCutParams.minAvgSepTrackPos;
      tTrackNeg = aPairCutParams.minAvgSepTrackNeg;
      break;

    case AliFemtoAnalysisLambdaKaon::kLamKchM:
    case AliFemtoAnalysisLambdaKaon::kALamKchM:
    case AliFemtoAnalysisLambdaKaon::kLamPiM:
    case AliFemtoAnalysisLambdaKaon::kALamPiM:
      tTrackPos = aPairCutParams.minAvgSepTrackNeg;
      tTrackNeg = aPairCutParams.minAvgSepTrackPos;
      break;

    default:
      break;
    }
  }

  v0TrackPairCut1->SetMinAvgSeparation(0,tTrackPos);
  v0TrackPairCut1->SetMinAvgSeparation(1,tTrackNeg);

  return v0TrackPairCut1;
}


//___________________________________________________________________
AliFemtoXiTrackPairCut* AliFemtoAnalysisLambdaKaon::CreateXiTrackPairCut(PairCutParams &aPairCutParams)
{
  AliFemtoXiTrackPairCut *tXiTrackPairCut = new AliFemtoXiTrackPairCut();

  return tXiTrackPairCut;
}






//___________________________________________________________________
AliFemtoCorrFctnKStar* AliFemtoAnalysisLambdaKaon::CreateCorrFctnKStar(const char* name, unsigned int bins, double min, double max)
{
  AliFemtoCorrFctnKStar *cf = new AliFemtoCorrFctnKStar(TString::Format("KStarCf_%s",name),bins,min,max);
    cf->SetCalculatePairKinematics(fWritePairKinematics);
    cf->SetBuildkTBinned(false);
    cf->SetBuildmTBinned(true);
    cf->SetBuild3d(false);
  return cf;
}

//___________________________________________________________________
AliFemtoAvgSepCorrFctn* AliFemtoAnalysisLambdaKaon::CreateAvgSepCorrFctn(const char* name, unsigned int bins, double min, double max)
{
  AliFemtoAvgSepCorrFctn *cf = new AliFemtoAvgSepCorrFctn(TString::Format("AvgSepCf_%s", name),bins,min,max);

  switch(fGeneralAnalysisType) {
  case AliFemtoAnalysisLambdaKaon::kV0V0:
    cf->SetPairType(AliFemtoAvgSepCorrFctn::kV0s);
    break;

  case AliFemtoAnalysisLambdaKaon::kV0Track:
    cf->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
    break;

  case AliFemtoAnalysisLambdaKaon::kXiTrack:
    cf->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0); //TODO
    break;

  case AliFemtoAnalysisLambdaKaon::kTrackTrack:
    cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
    break;

  default:
    cerr << "E-myTrainAnalysisConstructor::CreateAvgSepCorrFctn" << endl;
  }
  return cf;
}


//___________________________________________________________________
AliFemtoModelCorrFctnKStarFull* AliFemtoAnalysisLambdaKaon::CreateModelCorrFctnKStarFull(const char* name, unsigned int bins, double min, double max)
{
  bool tUseWeightGenerator = true;

  AliFemtoModelCorrFctnKStarFull *cf = new AliFemtoModelCorrFctnKStarFull(TString::Format("KStarModelCf_%s",name),bins,min,max);
    cf->SetRemoveMisidentified(fAnalysisParams.removeMisidentifiedMCParticles);
    cf->SetExpectedPDGCodes((int)fParticlePDGType1,(int)fParticlePDGType2);

  if(tUseWeightGenerator)
  {
    AliFemtoModelWeightGeneratorBasicLednicky *tGenerator = new AliFemtoModelWeightGeneratorBasicLednicky();
    tGenerator->SetIdenticalParticles(false);
    tGenerator->SetParamAlpha(0.);
    if(fAnalysisType == AliFemtoAnalysisLambdaKaon::kLamKchP || fAnalysisType == AliFemtoAnalysisLambdaKaon::kALamKchM)
    {
      tGenerator->SetParamLambda(0.1403);
      tGenerator->SetParamRadius(4.241);
      tGenerator->SetParamRef0(-1.981);
      tGenerator->SetParamImf0(0.8138);
      tGenerator->SetParamd0(2.621);
      tGenerator->SetParamNorm(1.);    
    }
    else if(fAnalysisType == AliFemtoAnalysisLambdaKaon::kLamKchM || fAnalysisType == AliFemtoAnalysisLambdaKaon::kALamKchP)
    {
      tGenerator->SetParamLambda(0.331);
      tGenerator->SetParamRadius(4.107);
      tGenerator->SetParamRef0(0.1362);
      tGenerator->SetParamImf0(0.4482);
      tGenerator->SetParamd0(6.666);
      tGenerator->SetParamNorm(1.);    
    }
    else if(fAnalysisType == AliFemtoAnalysisLambdaKaon::kLamK0 || fAnalysisType == AliFemtoAnalysisLambdaKaon::kALamK0)
    {
      tGenerator->SetParamLambda(0.232);
      tGenerator->SetParamRadius(3.08);
      tGenerator->SetParamRef0(-0.3319);
      tGenerator->SetParamImf0(0.3922);
      tGenerator->SetParamd0(-9.093);
      tGenerator->SetParamNorm(1.);    
    }

    AliFemtoModelManager *tManager = new AliFemtoModelManager();
    tManager->AcceptWeightGenerator((AliFemtoModelWeightGenerator*)tGenerator);

    cf->ConnectToManager(tManager);
  }

  return cf;
}


//___________________________________________________________________
AliFemtoV0PurityBgdEstimator* AliFemtoAnalysisLambdaKaon::CreateV0PurityBgdEstimator()
{
  TString tName = "V0PurityBgdEstimator_";
  V0CutParams tV0CutParams;
  unsigned int tNbins = 100;
  double tMinvMin, tMinvMax;

  switch(fAnalysisType) {
  case AliFemtoAnalysisLambdaKaon::kProtPiM:
    tName += TString("Lambda");
    tMinvMin = LambdaMass-0.035;
    tMinvMax = LambdaMass+0.035;
    tV0CutParams = DefaultLambdaCutParams();
    break;

  case AliFemtoAnalysisLambdaKaon::kAProtPiP:
    tName += TString("AntiLambda");
    tMinvMin = LambdaMass-0.035;
    tMinvMax = LambdaMass+0.035;
    tV0CutParams = DefaultAntiLambdaCutParams();
    break;

  case AliFemtoAnalysisLambdaKaon::kPiPPiM:
    tName += TString("K0Short");
    tMinvMin = K0ShortMass-0.070;
    tMinvMax = K0ShortMass+0.070;
    tV0CutParams = DefaultK0ShortCutParams();
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::CreateV0PurityBgdEstimator" << endl;
  }

  AliFemtoV0PurityBgdEstimator *tV0PurityBgdEstimator = new AliFemtoV0PurityBgdEstimator(tName,tNbins,tMinvMin,tMinvMax);
  AliFemtoV0TrackCutNSigmaFilter* tV0TrackCut = CreateV0Cut(tV0CutParams);
  if(fAnalysisType==kProtPiM || fAnalysisType==kAProtPiP) tV0TrackCut->SetInvariantMassLambda(tMinvMin,tMinvMax);
  else tV0TrackCut->SetInvariantMassK0s(tMinvMin,tMinvMax);
  tV0PurityBgdEstimator->SetV0TrackCut(tV0TrackCut);

  return tV0PurityBgdEstimator;
}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::AddCutMonitors(AliFemtoEventCut* aEventCut, AliFemtoParticleCut* aPartCut1, AliFemtoParticleCut* aPartCut2, AliFemtoPairCut* aPairCut)
{
  aEventCut->AddCutMonitorPass(new AliFemtoCutMonitorEventMult("_EvPass"));
  aEventCut->AddCutMonitorPass(new AliFemtoCutMonitorEventPartCollSize("_Part1",100,0,100,"_Part2",100,0,100));
  if(!fAnalysisParams.monitorEvCutPassOnly) aEventCut->AddCutMonitorFail(new AliFemtoCutMonitorEventMult("_EvFail"));

  TString tPartName1, tPartName2;

  int tPartType1 = -1;  //Only used if ESD track
  double tPartMass1 = 0;  //Only used if ESD track

  switch(fParticlePDGType1) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    tPartName1 = TString("_Lam");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    tPartName1 = TString("_ALam");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGK0:
    tPartName1 = TString("_K0");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGXiC:
    tPartName1 = TString("_Xi");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGAXiC:
    tPartName1 = TString("_AXi");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGProt:
    tPartName1 = TString("_Prot");
    tPartType1 = 2;
    tPartMass1 = ProtonMass;
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGAntiProt:
    tPartName1 = TString("_AProt");
    tPartType1 = 2;
    tPartMass1 = ProtonMass;
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGPiP:
    tPartName1 = TString("_PiP");
    tPartType1 = 0;
    tPartMass1 = PionMass;
    break;
  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCutMonitors" << endl;
  }

  int tPartType2 = -1;  //Only used if ESD track
  double tPartMass2 = 0;  //Only used if ESD track

  switch(fParticlePDGType2) {
  case AliFemtoAnalysisLambdaKaon::kPDGLam:
    tPartName2 = TString("_Lam");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGALam:
    tPartName2 = TString("_ALam");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGK0:
    tPartName2 = TString("_K0");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGXiC:
    tPartName2 = TString("_Xi");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGAXiC:
    tPartName2 = TString("_AXi");
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGKchP:
    tPartName2 = TString("_KchP");
    tPartType2 = 1;
    tPartMass2 = KchMass;
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGKchM:
    tPartName2 = TString("_KchM");
    tPartType2 = 1;
    tPartMass2 = KchMass;
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGPiP:
    tPartName2 = TString("_PiP");
    tPartType2 = 0;
    tPartMass2 = PionMass;
    break;
  case AliFemtoAnalysisLambdaKaon::kPDGPiM:
    tPartName2 = TString("_PiM");
    tPartType2 = 0;
    tPartMass2 = PionMass;
    break;
  default:
    cerr << "AliFemtoAnalysisLambdaKaon::AddCutMonitors" << endl;
  }

  //------------------------------------------------

  TString tNamePass1 = tPartName1 + TString("_Pass");
  TString tNameFail1 = tPartName1 + TString("_Fail");

  TString tNamePass2 = tPartName2 + TString("_Pass");
  TString tNameFail2 = tPartName2 + TString("_Fail");

  switch(fGeneralAnalysisType) {
  case AliFemtoAnalysisLambdaKaon::kV0V0:
    aPartCut1->AddCutMonitorPass(new AliFemtoCutMonitorV0(tNamePass1));
    if(!fAnalysisParams.monitorPart1CutPassOnly) aPartCut1->AddCutMonitorFail(new AliFemtoCutMonitorV0(tNameFail1));

    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorV0(tNamePass2));
    if(!fAnalysisParams.monitorPart2CutPassOnly) aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorV0(tNameFail2));
    break;

  case AliFemtoAnalysisLambdaKaon::kV0Track:
    aPartCut1->AddCutMonitorPass(new AliFemtoCutMonitorV0(tNamePass1));
    if(!fAnalysisParams.monitorPart1CutPassOnly) aPartCut1->AddCutMonitorFail(new AliFemtoCutMonitorV0(tNameFail1));

    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticleYPt(tNamePass2,tPartMass2));
    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticlePID(tNamePass2,tPartType2));
    if(!fAnalysisParams.monitorPart2CutPassOnly)
    {
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticleYPt(tNameFail2,tPartMass2));
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticlePID(tNameFail2,tPartType2));
    }
    break;

  case AliFemtoAnalysisLambdaKaon::kXiTrack:
    aPartCut1->AddCutMonitorPass(new AliFemtoCutMonitorXi(tNamePass1));
    if(!fAnalysisParams.monitorPart1CutPassOnly) aPartCut1->AddCutMonitorFail(new AliFemtoCutMonitorXi(tNameFail1));

    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticleYPt(tNamePass2,tPartMass2));
    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticlePID(tNamePass2,tPartType2));
    if(!fAnalysisParams.monitorPart2CutPassOnly)
    {
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticleYPt(tNameFail2,tPartMass2));
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticlePID(tNameFail2,tPartType2));
    }
    break;

  case AliFemtoAnalysisLambdaKaon::kTrackTrack:
    aPartCut1->AddCutMonitorPass(new AliFemtoCutMonitorParticleYPt(tNamePass1,tPartMass1));
    aPartCut1->AddCutMonitorPass(new AliFemtoCutMonitorParticlePID(tNamePass1,tPartType1));
    if(!fAnalysisParams.monitorPart1CutPassOnly)
    {
      aPartCut1->AddCutMonitorFail(new AliFemtoCutMonitorParticleYPt(tNameFail1,tPartMass1));
      aPartCut1->AddCutMonitorFail(new AliFemtoCutMonitorParticlePID(tNameFail1,tPartType1));
    }

    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticleYPt(tNamePass2,tPartMass2));
    aPartCut2->AddCutMonitorPass(new AliFemtoCutMonitorParticlePID(tNamePass2,tPartType2));
    if(!fAnalysisParams.monitorPart2CutPassOnly)
    {
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticleYPt(tNameFail2,tPartMass2));
      aPartCut2->AddCutMonitorFail(new AliFemtoCutMonitorParticlePID(tNameFail2,tPartType2));
    }
    break;

  default:
    cerr << "E-AliFemtoAnalysisLambdaKaon::AddCutMonitors" << endl;
  }

  //------------------------------------------------
  if(fIsMCRun) aPairCut->AddCutMonitorPass(new AliFemtoPairOriginMonitor("Pass"));

}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::SetAnalysis(AliFemtoEventCut* aEventCut, AliFemtoParticleCut* aPartCut1, AliFemtoParticleCut* aPartCut2, AliFemtoPairCut* aPairCut)
{
  //AliFemtoCutMonitorHandler does not implement deep-copying on the AliFemtoCutMonitorCollection objects
  //Therefore, it is necessary to create the cut monitors uniquely for each analysis, making it impossible to simply throw pre-made particle cuts, etc. into
  //new analysis objects.  The cut monitors would be shared thorughout all analyses
  AliFemtoEventCut* tEventCut = aEventCut->Clone();
  AliFemtoParticleCut* tPartCut1 = aPartCut1->Clone(); 
  AliFemtoParticleCut* tPartCut2 = aPartCut2->Clone(); 
  AliFemtoPairCut* tPairCut = aPairCut->Clone();

  AddCutMonitors(tEventCut,tPartCut1,tPartCut2,tPairCut);

  SetEventCut(tEventCut);
  SetFirstParticleCut(tPartCut1);
  SetSecondParticleCut(tPartCut2);
  SetPairCut(tPairCut);

  AliFemtoCorrFctnIterator iter;
  for(iter=fCollectionOfCfs->begin(); iter!=fCollectionOfCfs->end(); iter++)
  {
    AddCorrFctn(*iter);
  }
}

//___________________________________________________________________
void AliFemtoAnalysisLambdaKaon::SetMultHist(const char* name, int aNbins, double aMin, double aMax)
{
  fMultHist = new TH1F(
                       TString::Format("MultHist_%s", name), "Multiplicity",
                       aNbins, aMin, aMax);
  fBuildMultHist = true;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::AnalysisParams 
AliFemtoAnalysisLambdaKaon::DefaultAnalysisParams()
{
  AliFemtoAnalysisLambdaKaon::AnalysisParams tReturnParams;

  tReturnParams.nBinsVertex = 10;
  tReturnParams.minVertex = -10.;
  tReturnParams.maxVertex = 10.;

  tReturnParams.nBinsMult = 20;
  tReturnParams.minMult = 0.;
  tReturnParams.maxMult = 1000.;

  tReturnParams.analysisType = AliFemtoAnalysisLambdaKaon::kLamK0;
  tReturnParams.generalAnalysisType = AliFemtoAnalysisLambdaKaon::kV0V0;

  tReturnParams.nEventsToMix = 5;
  tReturnParams.minCollectionSize = 1;


  tReturnParams.verbose = false;
  tReturnParams.implementAvgSepCuts = true;
  tReturnParams.writePairKinematics = false;
  tReturnParams.isMCRun = false;
  tReturnParams.isMBAnalysis = false;
  tReturnParams.buildMultHist = false;
  tReturnParams.implementVertexCorrections = true;
  tReturnParams.removeMisidentifiedMCParticles = false;
  tReturnParams.setV0SharedDaughterCut = true;

  tReturnParams.monitorEvCutPassOnly = false;
  tReturnParams.monitorPart1CutPassOnly = false;
  tReturnParams.monitorPart2CutPassOnly = false;
  tReturnParams.monitorPairCutPassOnly = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::EventCutParams 
AliFemtoAnalysisLambdaKaon::DefaultEventCutParams()
{
  AliFemtoAnalysisLambdaKaon::EventCutParams tReturnParams;

  tReturnParams.minCentrality = 0.;
  tReturnParams.maxCentrality = 100.;

  tReturnParams.minMult = 0.;
  tReturnParams.maxMult = 100000.;

  tReturnParams.minVertexZ = -8.;
  tReturnParams.maxVertexZ = 8.;

  tReturnParams.verboseMode = false;

  return tReturnParams;
}


//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::V0CutParams 
AliFemtoAnalysisLambdaKaon::DefaultLambdaCutParams()
{
  AliFemtoAnalysisLambdaKaon::V0CutParams tReturnParams;

  tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGLam;
  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kV0;

  tReturnParams.v0Type = 0;

  tReturnParams.mass = LambdaMass;
  tReturnParams.minInvariantMass = LambdaMass-0.0038;
  tReturnParams.maxInvariantMass = LambdaMass+0.0038;

  tReturnParams.useLooseInvMassCut = true;
  tReturnParams.minLooseInvMass = LambdaMass-0.035;
  tReturnParams.maxLooseInvMass = LambdaMass+0.035;

  tReturnParams.nBinsPurity = 100;
  tReturnParams.minPurityMass = LambdaMass-0.035;
  tReturnParams.maxPurityMass = LambdaMass+0.035;

  tReturnParams.useCustomFilter = true;

  tReturnParams.removeMisID = true;
  tReturnParams.minInvMassReject = K0ShortMass-0.009;
  tReturnParams.maxInvMassReject = K0ShortMass+0.009;

  tReturnParams.useSimpleMisID = false;
  tReturnParams.buildMisIDHistograms = true;
  tReturnParams.useCustomMisID = true;

  tReturnParams.eta = 0.8;
  tReturnParams.minPt = 0.4;
  tReturnParams.maxPt = 100.;
  tReturnParams.onFlyStatus = false;
  tReturnParams.maxDcaV0 = 0.5;
  tReturnParams.minCosPointingAngle = 0.9993;
  tReturnParams.maxV0DecayLength = 60.;

  tReturnParams.etaDaughters = 0.8;
  tReturnParams.minPtPosDaughter = 0.5;
  tReturnParams.maxPtPosDaughter = 99.;
  tReturnParams.minPtNegDaughter = 0.16;
  tReturnParams.maxPtNegDaughter = 99.;
  tReturnParams.minTPCnclsDaughters = 80;
  tReturnParams.maxDcaV0Daughters = 0.4;
  tReturnParams.minPosDaughterToPrimVertex = 0.1;
  tReturnParams.minNegDaughterToPrimVertex = 0.3;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::V0CutParams 
AliFemtoAnalysisLambdaKaon::DefaultAntiLambdaCutParams()
{
  AliFemtoAnalysisLambdaKaon::V0CutParams tReturnParams;

  tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGALam;
  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kV0;

  tReturnParams.v0Type = 1;

  tReturnParams.mass = LambdaMass;
  tReturnParams.minInvariantMass = LambdaMass-0.0038;
  tReturnParams.maxInvariantMass = LambdaMass+0.0038;

  tReturnParams.useLooseInvMassCut = true;
  tReturnParams.minLooseInvMass = LambdaMass-0.035;
  tReturnParams.maxLooseInvMass = LambdaMass+0.035;

  tReturnParams.nBinsPurity = 100;
  tReturnParams.minPurityMass = LambdaMass-0.035;
  tReturnParams.maxPurityMass = LambdaMass+0.035;

  tReturnParams.useCustomFilter = true;

  tReturnParams.removeMisID = true;
  tReturnParams.minInvMassReject = K0ShortMass-0.009;
  tReturnParams.maxInvMassReject = K0ShortMass+0.009;

  tReturnParams.useSimpleMisID = false;
  tReturnParams.buildMisIDHistograms = true;
  tReturnParams.useCustomMisID = true;

  tReturnParams.eta = 0.8;
  tReturnParams.minPt = 0.4;
  tReturnParams.maxPt = 100.;
  tReturnParams.onFlyStatus = false;
  tReturnParams.maxDcaV0 = 0.5;
  tReturnParams.minCosPointingAngle = 0.9993;
  tReturnParams.maxV0DecayLength = 60.;

  tReturnParams.etaDaughters = 0.8;
  tReturnParams.minPtPosDaughter = 0.16;
  tReturnParams.maxPtPosDaughter = 99.;
  tReturnParams.minPtNegDaughter = 0.3;
  tReturnParams.maxPtNegDaughter = 99.;
  tReturnParams.minTPCnclsDaughters = 80;
  tReturnParams.maxDcaV0Daughters = 0.4;
  tReturnParams.minPosDaughterToPrimVertex = 0.3;
  tReturnParams.minNegDaughterToPrimVertex = 0.1;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::V0CutParams 
AliFemtoAnalysisLambdaKaon::DefaultK0ShortCutParams()
{
  AliFemtoAnalysisLambdaKaon::V0CutParams tReturnParams;

  tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGK0;
  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kV0;

  tReturnParams.v0Type = 2;

  tReturnParams.mass = K0ShortMass;
  tReturnParams.minInvariantMass = K0ShortMass-0.013677;
  tReturnParams.maxInvariantMass = K0ShortMass+0.020323;

  tReturnParams.useLooseInvMassCut = true;
  tReturnParams.minLooseInvMass = K0ShortMass-0.070;
  tReturnParams.maxLooseInvMass = K0ShortMass+0.070;

  tReturnParams.nBinsPurity = 100;
  tReturnParams.minPurityMass = K0ShortMass-0.070;
  tReturnParams.maxPurityMass = K0ShortMass+0.070;

  tReturnParams.useCustomFilter = true;

  tReturnParams.removeMisID = true;
  tReturnParams.minInvMassReject = LambdaMass-0.009;
  tReturnParams.maxInvMassReject = LambdaMass+0.009;

  tReturnParams.useSimpleMisID = false;
  tReturnParams.buildMisIDHistograms = true;
  tReturnParams.useCustomMisID = true;

  tReturnParams.eta = 0.8;
  tReturnParams.minPt = 0.2;
  tReturnParams.maxPt = 100.;
  tReturnParams.onFlyStatus = false;
  tReturnParams.maxDcaV0 = 0.3;
  tReturnParams.minCosPointingAngle = 0.9993;
  tReturnParams.maxV0DecayLength = 30.;

  tReturnParams.etaDaughters = 0.8;
  tReturnParams.minPtPosDaughter = 0.15;
  tReturnParams.maxPtPosDaughter = 99.;
  tReturnParams.minPtNegDaughter = 0.15;
  tReturnParams.maxPtNegDaughter = 99.;
  tReturnParams.minTPCnclsDaughters = 80;
  tReturnParams.maxDcaV0Daughters = 0.3;
  tReturnParams.minPosDaughterToPrimVertex = 0.3;
  tReturnParams.minNegDaughterToPrimVertex = 0.3;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::ESDCutParams 
AliFemtoAnalysisLambdaKaon::DefaultKchCutParams(int aCharge)
{
  AliFemtoAnalysisLambdaKaon::ESDCutParams tReturnParams;

  if(aCharge>0) tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGKchP;
  else tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGKchM;

  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kTrack;

  tReturnParams.minPidProbPion = 0.;
  tReturnParams.maxPidProbPion = 0.1;
  tReturnParams.minPidProbMuon = 0.;
  tReturnParams.maxPidProbMuon = 0.8;
  tReturnParams.minPidProbKaon = 0.2;
  tReturnParams.maxPidProbKaon = 1.001;
  tReturnParams.minPidProbProton = 0.;
  tReturnParams.maxPidProbProton = 0.1;
  tReturnParams.mostProbable = 3;    //this uses P().Mag() as first argument to IsKaonNSigma()
//  tReturnParams.mostProbable = 11; //this looks for Kaons, and uses Pt() as first argument to IsKaonNSigma
  tReturnParams.charge = aCharge;
  tReturnParams.mass = KchMass;

  tReturnParams.minPt = 0.14;
  tReturnParams.maxPt = 1.5;
  tReturnParams.eta = 0.8;
  tReturnParams.minTPCncls = 80;

  tReturnParams.removeKinks = true;
  tReturnParams.setLabel = false;
  tReturnParams.maxITSChiNdof = 3.0;
  tReturnParams.maxTPCChiNdof = 4.0;
  tReturnParams.maxSigmaToVertex = 3.0;
  tReturnParams.minImpactXY = -1000.0;
  tReturnParams.maxImpactXY = 2.4;
  tReturnParams.maxImpactZ = 3.0;

  tReturnParams.useCustomFilter = true;
  tReturnParams.useCustomMisID = true;
  tReturnParams.useElectronRejection = true;
  tReturnParams.useCustomElectronRejection = true;
  tReturnParams.usePionRejection = true;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::ESDCutParams 
AliFemtoAnalysisLambdaKaon::DefaultPiCutParams(int aCharge)
{
  AliFemtoAnalysisLambdaKaon::ESDCutParams tReturnParams;

  if(aCharge>0) tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiP;
  else tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiM;

  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kTrack;

  tReturnParams.minPidProbPion = 0.2;
  tReturnParams.maxPidProbPion = 1.001;
  tReturnParams.minPidProbMuon = 0.;
  tReturnParams.maxPidProbMuon = 0.8;
  tReturnParams.minPidProbKaon = 0.0;
  tReturnParams.maxPidProbKaon = 0.1;
  tReturnParams.minPidProbProton = 0.;
  tReturnParams.maxPidProbProton = 0.1;
  tReturnParams.mostProbable = 2;    //this uses P().Mag() as first argument to IsPionNSigma()
//  tReturnParams.mostProbable = 10; //this looks for Pions, and uses Pt() as first argument to IsPionNSigma
  tReturnParams.charge = aCharge;
  tReturnParams.mass = PionMass;

  tReturnParams.minPt = 0.1;
  tReturnParams.maxPt = 2.0;
  tReturnParams.eta = 0.8;
  tReturnParams.minTPCncls = 80;

  tReturnParams.removeKinks = true;
  tReturnParams.setLabel = false;
  tReturnParams.maxITSChiNdof = 3.0;
  tReturnParams.maxTPCChiNdof = 2.0;
  tReturnParams.maxSigmaToVertex = 3.0;
  tReturnParams.minImpactXY = -1000.0;
  tReturnParams.maxImpactXY = 2.4;
  tReturnParams.maxImpactZ = 3.2;

  tReturnParams.useCustomFilter = false;
  tReturnParams.useCustomMisID = false;
  tReturnParams.useElectronRejection = false;
  tReturnParams.useCustomElectronRejection = false;
  tReturnParams.usePionRejection = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::XiCutParams 
AliFemtoAnalysisLambdaKaon::DefaultXiCutParams()
{
  AliFemtoAnalysisLambdaKaon::XiCutParams tReturnParams;

  tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGXiC;
  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kCascade;

  tReturnParams.charge = -1;
  tReturnParams.xiType = 0;
  tReturnParams.minPt = 0.8;
  tReturnParams.maxPt = 100.;
  tReturnParams.eta = 0.8;
  tReturnParams.mass = XiMass;
  tReturnParams.minInvariantMass = XiMass-0.003;
  tReturnParams.maxInvariantMass = XiMass+0.003;

  tReturnParams.maxDecayLengthXi = 100.;
  tReturnParams.minCosPointingAngleXi = 0.9992;
  tReturnParams.minCosPointingAngleV0toXi = 0.9993;
  tReturnParams.maxDcaXi = 100.;
  tReturnParams.maxDcaXiDaughters = 0.3;

  tReturnParams.minDcaXiBac = 0.03;
  tReturnParams.etaBac = 0.8;
  tReturnParams.minTPCnclsBac = 70;
  tReturnParams.minPtBac = 0.;
  tReturnParams.maxPtBac = 100.;

  tReturnParams.v0Type = 0;
  tReturnParams.minDcaV0 = 0.1;
  tReturnParams.minInvMassV0 = LambdaMass-0.005;
  tReturnParams.maxInvMassV0 = LambdaMass+0.005;
  tReturnParams.minCosPointingAngleV0 = 0.;  //TODO was 0.998, might need to revert back
                                             //changed because of new minCosPointingAngleV0toXi
  tReturnParams.etaV0 = 0.8;
  tReturnParams.minPtV0 = 0.4;
  tReturnParams.maxPtV0 = 100.;
  tReturnParams.onFlyStatusV0 = false;
  tReturnParams.maxV0DecayLength = 60.;
  tReturnParams.minV0PosDaughterToPrimVertex = 0.1;
  tReturnParams.minV0NegDaughterToPrimVertex = 0.3;
  tReturnParams.maxDcaV0Daughters = 0.4;
  tReturnParams.etaV0Daughters = 0.8;
  tReturnParams.minPtPosV0Daughter = 0.5;
  tReturnParams.maxPtPosV0Daughter = 99.;
  tReturnParams.minPtNegV0Daughter = 0.16;
  tReturnParams.maxPtNegV0Daughter = 99.;

  tReturnParams.minTPCnclsV0Daughters = 70;

  tReturnParams.useCustomV0Filter = true;
  tReturnParams.useCustomV0MisID = false;
  tReturnParams.useCustomBacPionFilter = true;
  tReturnParams.useCustomBacPionMisID = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::XiCutParams 
AliFemtoAnalysisLambdaKaon::DefaultAXiCutParams()
{
  AliFemtoAnalysisLambdaKaon::XiCutParams tReturnParams;

  tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGAXiC;
  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kCascade;

  tReturnParams.charge = 1;
  tReturnParams.xiType = 1;
  tReturnParams.minPt = 0.8;
  tReturnParams.maxPt = 100.;
  tReturnParams.eta = 0.8;
  tReturnParams.mass = XiMass;
  tReturnParams.minInvariantMass = XiMass-0.003;
  tReturnParams.maxInvariantMass = XiMass+0.003;

  tReturnParams.maxDecayLengthXi = 100.;
  tReturnParams.minCosPointingAngleXi = 0.9992;
  tReturnParams.minCosPointingAngleV0toXi = 0.9993;
  tReturnParams.maxDcaXi = 100.;
  tReturnParams.maxDcaXiDaughters = 0.3;

  tReturnParams.minDcaXiBac = 0.03;
  tReturnParams.etaBac = 0.8;
  tReturnParams.minTPCnclsBac = 70;
  tReturnParams.minPtBac = 0.;
  tReturnParams.maxPtBac = 100.;

  tReturnParams.v0Type = 1;
  tReturnParams.minDcaV0 = 0.1;
  tReturnParams.minInvMassV0 = LambdaMass-0.005;
  tReturnParams.maxInvMassV0 = LambdaMass+0.005;
  tReturnParams.minCosPointingAngleV0 = 0.;  //TODO was 0.998, might need to revert back
                                             //changed because of new minCosPointingAngleV0toXi
  tReturnParams.etaV0 = 0.8;
  tReturnParams.minPtV0 = 0.4;
  tReturnParams.maxPtV0 = 100.;
  tReturnParams.onFlyStatusV0 = true;
  tReturnParams.maxV0DecayLength = 60.;
  tReturnParams.minV0PosDaughterToPrimVertex = 0.3;
  tReturnParams.minV0NegDaughterToPrimVertex = 0.1;
  tReturnParams.maxDcaV0Daughters = 0.4;
  tReturnParams.etaV0Daughters = 0.8;
  tReturnParams.minPtPosV0Daughter = 0.16;
  tReturnParams.maxPtPosV0Daughter = 99.;
  tReturnParams.minPtNegV0Daughter = 0.3;
  tReturnParams.maxPtNegV0Daughter = 99.;

  tReturnParams.minTPCnclsV0Daughters = 70;

  tReturnParams.useCustomV0Filter = true;
  tReturnParams.useCustomV0MisID = false;
  tReturnParams.useCustomBacPionFilter = false;
  tReturnParams.useCustomBacPionMisID = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::PairCutParams 
AliFemtoAnalysisLambdaKaon::DefaultPairParams()
{
  AliFemtoAnalysisLambdaKaon::PairCutParams tReturnParams;

  tReturnParams.removeSameLabel = true;
  tReturnParams.shareQualityMax = 1.0;
  tReturnParams.shareFractionMax = 1.0;
  tReturnParams.tpcOnly = true;

  tReturnParams.tpcExitSepMinimum = -1.0;
  tReturnParams.tpcEntranceSepMinimum = 0.00001;

  //Default for LamK0 and ALamK0
  tReturnParams.minAvgSepPosPos = 6.;
  tReturnParams.minAvgSepPosNeg = 0.;
  tReturnParams.minAvgSepNegPos = 0.;
  tReturnParams.minAvgSepNegNeg = 6.;

  //Default for LamKchP and ALamKchP
  tReturnParams.minAvgSepTrackPos = 8.;
  tReturnParams.minAvgSepTrackNeg = 0.;

  return tReturnParams;
}


//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::ESDCutParams 
AliFemtoAnalysisLambdaKaon::LambdaPurityPiCutParams(int aCharge)
{
  //Used with AliFemtoV0PurityBgdEstimator to estimate background
  //in V0 Minv plot used to calculate purity

  // aCharge = +1 ==> AntiLambda is V0
  // aCharge = -1 ==> Lambda is V0

  AliFemtoAnalysisLambdaKaon::ESDCutParams tReturnParams;

  if(aCharge>0) tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiP;
  else tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiM;

  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kTrack;
/*
  tReturnParams.minPidProbPion = 0.2;
  tReturnParams.maxPidProbPion = 1.001;
  tReturnParams.minPidProbMuon = 0.;
  tReturnParams.maxPidProbMuon = 0.8;
  tReturnParams.minPidProbKaon = 0.0;
  tReturnParams.maxPidProbKaon = 0.1;
  tReturnParams.minPidProbProton = 0.;
  tReturnParams.maxPidProbProton = 0.1;
*/
//TODO temporarily turn off these cuts
  tReturnParams.minPidProbPion = -1;
  tReturnParams.maxPidProbPion = 2;
  tReturnParams.minPidProbMuon = -1;
  tReturnParams.maxPidProbMuon = 2;
  tReturnParams.minPidProbKaon = -1;
  tReturnParams.maxPidProbKaon = 2;
  tReturnParams.minPidProbProton = -1;
  tReturnParams.maxPidProbProton = 2;

  tReturnParams.mostProbable = 2;    //this uses P().Mag() as first argument to IsPionNSigma()
//  tReturnParams.mostProbable = 10; //this looks for Kaons, and uses Pt() as first argument to IsPionNSigma
  tReturnParams.charge = aCharge;
  tReturnParams.mass = PionMass;

  tReturnParams.minPt = 0.16;
  tReturnParams.maxPt = 99.;
  tReturnParams.eta = 0.8;
  tReturnParams.minTPCncls = 80;

  tReturnParams.removeKinks = true;
  tReturnParams.setLabel = false;
  tReturnParams.maxITSChiNdof = 1000;
  tReturnParams.maxTPCChiNdof = 1000;
  tReturnParams.maxSigmaToVertex = 1000;
  tReturnParams.minImpactXY = 0.3;
  tReturnParams.maxImpactXY = 1000;  //TODO may need adjusted
  tReturnParams.maxImpactZ = 1000;  //TODO may need adjusted

  tReturnParams.useCustomFilter = true;
  tReturnParams.useCustomMisID = false;
  tReturnParams.useElectronRejection = false;
  tReturnParams.useCustomElectronRejection = false;
  tReturnParams.usePionRejection = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::ESDCutParams 
AliFemtoAnalysisLambdaKaon::K0ShortPurityPiCutParams(int aCharge)
{
  //Used with AliFemtoV0PurityBgdEstimator to estimate background
  //in V0 Minv plot used to calculate purity

  AliFemtoAnalysisLambdaKaon::ESDCutParams tReturnParams;

  if(aCharge>0) tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiP;
  else tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGPiM;

  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kTrack;
/*
  tReturnParams.minPidProbPion = 0.2;
  tReturnParams.maxPidProbPion = 1.001;
  tReturnParams.minPidProbMuon = 0.;
  tReturnParams.maxPidProbMuon = 0.8;
  tReturnParams.minPidProbKaon = 0.0;
  tReturnParams.maxPidProbKaon = 0.1;
  tReturnParams.minPidProbProton = 0.;
  tReturnParams.maxPidProbProton = 0.1;
*/
//TODO temporarily turn off these cuts
  tReturnParams.minPidProbPion = -1;
  tReturnParams.maxPidProbPion = 2;
  tReturnParams.minPidProbMuon = -1;
  tReturnParams.maxPidProbMuon = 2;
  tReturnParams.minPidProbKaon = -1;
  tReturnParams.maxPidProbKaon = 2;
  tReturnParams.minPidProbProton = -1;
  tReturnParams.maxPidProbProton = 2;

  tReturnParams.mostProbable = 2;    //this uses P().Mag() as first argument to IsPionNSigma()
//  tReturnParams.mostProbable = 10; //this looks for Kaons, and uses Pt() as first argument to IsPionNSigma
  tReturnParams.charge = aCharge;
  tReturnParams.mass = PionMass;

  tReturnParams.minPt = 0.15;
  tReturnParams.maxPt = 99.;
  tReturnParams.eta = 0.8;
  tReturnParams.minTPCncls = 80;

  tReturnParams.removeKinks = true;
  tReturnParams.setLabel = false;
  tReturnParams.maxITSChiNdof = 1000;
  tReturnParams.maxTPCChiNdof = 1000;
  tReturnParams.maxSigmaToVertex = 1000;
  tReturnParams.minImpactXY = 0.3;
  tReturnParams.maxImpactXY = 1000;  //TODO may need adjusted
  tReturnParams.maxImpactZ = 1000;  //TODO may need adjusted

  tReturnParams.useCustomFilter = true;
  tReturnParams.useCustomMisID = false;
  tReturnParams.useElectronRejection = false;
  tReturnParams.useCustomElectronRejection = false;
  tReturnParams.usePionRejection = false;

  return tReturnParams;
}

//___________________________________________________________________
AliFemtoAnalysisLambdaKaon::ESDCutParams 
AliFemtoAnalysisLambdaKaon::LambdaPurityProtonCutParams(int aCharge)
{
  //Used with AliFemtoV0PurityBgdEstimator to estimate background
  //in V0 Minv plot used to calculate purity

  // aCharge = +1 ==> Lambda is V0
  // aCharge = -1 ==> AntiLambda is V0

  AliFemtoAnalysisLambdaKaon::ESDCutParams tReturnParams;

  if(aCharge>0) tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGProt;
  else tReturnParams.particlePDGType = AliFemtoAnalysisLambdaKaon::kPDGAntiProt;

  tReturnParams.generalParticleType = AliFemtoAnalysisLambdaKaon::kTrack;
/*
  tReturnParams.minPidProbPion = 0.;
  tReturnParams.maxPidProbPion = 0.1;
  tReturnParams.minPidProbMuon = 0.;
  tReturnParams.maxPidProbMuon = 0.8;
  tReturnParams.minPidProbKaon = 0.0;
  tReturnParams.maxPidProbKaon = 0.1;
  tReturnParams.minPidProbProton = 0.2;
  tReturnParams.maxPidProbProton = 1.001;
*/
//TODO temporarily turn off these cuts
  tReturnParams.minPidProbPion = -1;
  tReturnParams.maxPidProbPion = 2;
  tReturnParams.minPidProbMuon = -1;
  tReturnParams.maxPidProbMuon = 2;
  tReturnParams.minPidProbKaon = -1;
  tReturnParams.maxPidProbKaon = 2;
  tReturnParams.minPidProbProton = -1;
  tReturnParams.maxPidProbProton = 2;

  tReturnParams.mostProbable = 4;    //this uses P().Mag() as first argument to IsProtonNSigma()
//  tReturnParams.mostProbable = 12; //this looks for Kaons, and uses Pt() as first argument to IsProtonNSigma
  tReturnParams.charge = aCharge;
  tReturnParams.mass = ProtonMass;

  if(aCharge>0) tReturnParams.minPt = 0.5;
  else tReturnParams.minPt = 0.3;
  tReturnParams.maxPt = 99.;
  tReturnParams.eta = 0.8;
  tReturnParams.minTPCncls = 80;

  tReturnParams.removeKinks = true;
  tReturnParams.setLabel = false;
  tReturnParams.maxITSChiNdof = 1000;
  tReturnParams.maxTPCChiNdof = 1000;
  tReturnParams.maxSigmaToVertex = 1000;
  tReturnParams.minImpactXY = 0.1;
  tReturnParams.maxImpactXY = 1000;  //TODO may need adjusted
  tReturnParams.maxImpactZ = 1000;  //TODO may need adjusted

  tReturnParams.useCustomFilter = true;
  tReturnParams.useCustomMisID = false;
  tReturnParams.useElectronRejection = false;
  tReturnParams.useCustomElectronRejection = false;
  tReturnParams.usePionRejection = false;

  return tReturnParams;
}








