///
/// \file PWGCF/FEMTOSCOPY/macros/Train/LambdaKaonFemto/ConfigFemtoAnalysisSystematics.C
/// \brief The configuration macro which sets up Lambda-Kaon analyses
/// \ author Jesse Buxton, Ohio State University,
///
/// \ based largely off of work done by Andrew Kubera in
/// \ file PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/ConfigFemtoAnalysis.C
/// \   Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliFemtoAnalysisLambdaKaon.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderAODChain.h"

#include <TROOT.h>
#endif

typedef AliFemtoAnalysisLambdaKaon AFALK;

bool DEFAULT_DO_KT = kFALSE;

struct MacroParams {
  std::vector<int> centrality_ranges;
  std::vector<int> pair_codes;
  float qinv_bin_size_MeV;
  float qinv_max_GeV;
  bool do_qinv_cf;
  bool do_q3d_cf;
  bool do_deltaeta_deltaphi_cf;
  bool do_avg_sep_cf;
  bool do_kt_q3d;
  bool do_kt_qinv;
  bool do_ylm_cf; // not implemented yet
  int filter_bit;
  AliFemtoEventReaderAOD::EventMult multiplicity;
  int dca_global_track;
};

struct CutVariations {
  std::vector<double> valuesToVary;  //ex ?[0.2:0.3:0.4:0.5:0.6]
  TString cutName;  //ex ?$|Lam|maxDcaV0Daughters
};

void SetPairCodes(AFALK::AnalysisParams &aAnConfig, MacroParams &aMacroConfig);

AliFemtoAnalysisLambdaKaon* 
CreateCorrectAnalysis(
  const TString &aText,
  AFALK::AnalysisType aAnType,
  AFALK::AnalysisParams &aAnParams,
  AFALK::EventCutParams &aEvCutParams,
  AFALK::PairCutParams &aPairCutParams,
  const TString &aDirNameModifier
);

void
BuildConfiguration(
  const TString &aText,
  AFALK::AnalysisParams &aAnParams,
  AFALK::EventCutParams &aEvCutParams,
  AFALK::PairCutParams &aPairCutParams,
  MacroParams &aMac,
  CutVariations &aCutVariations
);

void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::V0CutParams &aV0CutParams
);

void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::ESDCutParams &aESDCutParams
);


void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::XiCutParams &aXiCutParams
);

TString GetParticleCuts(const TString &aText);

AliFemtoManager* ConfigFemtoAnalysis(const TString& aParamString="") 
{
  std::cout << "[ConfigFemtoAnalysisSystematics (LambdaKaon)]\n";

  // Get the default configurations
  AFALK::AnalysisParams tAnalysisConfig = AFALK::DefaultAnalysisParams();
  AFALK::EventCutParams tEventCutConfig = AFALK::DefaultEventCutParams();
  AFALK::PairCutParams tPairCutConfig = AFALK::DefaultPairParams();

  // default Macro config
  MacroParams tMacroConfig;
  tMacroConfig.do_qinv_cf = true;
  tMacroConfig.do_q3d_cf = true;
  tMacroConfig.do_deltaeta_deltaphi_cf = false;
  tMacroConfig.do_avg_sep_cf = false;
  tMacroConfig.do_kt_q3d = tMacroConfig.do_kt_qinv = DEFAULT_DO_KT;
  tMacroConfig.do_ylm_cf = false;
  tMacroConfig.qinv_bin_size_MeV = 5.0f;
  tMacroConfig.qinv_max_GeV = 1.0f;
  tMacroConfig.filter_bit = 7;
  tMacroConfig.multiplicity = AliFemtoEventReaderAOD::kCentrality;
  tMacroConfig.dca_global_track = 0;

  //
  CutVariations tCutVariations;

  // Read parameter string and update configurations
  //Initial call to BuildConfiguration
  BuildConfiguration(aParamString,tAnalysisConfig,tEventCutConfig,tPairCutConfig,tMacroConfig,tCutVariations);
  TString tParticleCuts = GetParticleCuts(aParamString);

  if(tAnalysisConfig.analysisType==AFALK::kProtPiM || tAnalysisConfig.analysisType==AFALK::kAProtPiP ||
            tAnalysisConfig.analysisType==AFALK::kPiPPiM)
  {
    tMacroConfig.dca_global_track = 1;
    tMacroConfig.filter_bit = 0;
  }

  // Begin to build the manager and analyses
  AliFemtoManager *tManager = new AliFemtoManager();

  //Setup the event reader for ALICE AOD
  AliFemtoEventReaderAODChain *rdr = new AliFemtoEventReaderAODChain();
    rdr->SetUseMultiplicity(tMacroConfig.multiplicity);  //Sets the type of the event multiplicity estimator
    if(tAnalysisConfig.analysisType!=AFALK::kProtPiM && tAnalysisConfig.analysisType!=AFALK::kAProtPiP &&
       tAnalysisConfig.analysisType!=AFALK::kPiPPiM) rdr->SetFilterBit(tMacroConfig.filter_bit);
    //rdr->SetCentralityPreSelection(0, 900);
    if(tAnalysisConfig.analysisType==AFALK::kXiKchP || tAnalysisConfig.analysisType==AFALK::kAXiKchP ||
       tAnalysisConfig.analysisType==AFALK::kXiKchM || tAnalysisConfig.analysisType==AFALK::kAXiKchM) rdr->SetReadCascade(1);
    else if(tAnalysisConfig.analysisType==AFALK::kXiK0   || tAnalysisConfig.analysisType==AFALK::kAXiK0) {rdr->SetReadV0(1); rdr->SetReadCascade(1);}
    else if(tAnalysisConfig.analysisType==AFALK::kProtPiM || tAnalysisConfig.analysisType==AFALK::kAProtPiP ||
            tAnalysisConfig.analysisType==AFALK::kPiPPiM) rdr->SetReadV0(0); 
    else rdr->SetReadV0(1);  //Read V0 information from the AOD and put it into V0Collection
    rdr->SetEPVZERO(kTRUE);  //to get event plane angle from VZERO
    rdr->SetCentralityFlattening(kFALSE);
    rdr->SetPrimaryVertexCorrectionTPCPoints(tAnalysisConfig.implementVertexCorrections);
    rdr->SetReadMC(tAnalysisConfig.isMCRun);
    rdr->SetDCAglobalTrack(tMacroConfig.dca_global_track);
  tManager->SetEventReader(rdr);


  if(tMacroConfig.centrality_ranges.empty())
  {
    tMacroConfig.centrality_ranges.push_back(0.);
    tMacroConfig.centrality_ranges.push_back(10.);
  }

  // Identify all sister analyses that go along with tAnalysisConfig.analysisType
  SetPairCodes(tAnalysisConfig, tMacroConfig);

  //loop over the different cut values for the cut being varied
  double tNCutVals = tCutVariations.valuesToVary.size();
  for(unsigned int iCutVal=0; iCutVal<tNCutVals; iCutVal++)
  {
    //change the configuration
    TString tNewParamString = tParticleCuts;
    tNewParamString += tCutVariations.cutName + TString(" = ") + TString::Format("%f",tCutVariations.valuesToVary[iCutVal]) + TString(";");

    int tBeginDirNameModifier = 1;
    if(tCutVariations.cutName[1] == TString("|")) tBeginDirNameModifier = 2;
    TString tDirNameModifier = tCutVariations.cutName(tBeginDirNameModifier, tCutVariations.cutName.Length() - 1) + TString("|") + TString::Format("%0.6f",tCutVariations.valuesToVary[iCutVal]);
    tDirNameModifier.ReplaceAll("|","_");

    // loop over centrality ranges
    for(unsigned int iCent = 0; iCent+1 < tMacroConfig.centrality_ranges.size(); iCent += 2)
    {
      const int tMultLow  = tMacroConfig.centrality_ranges[iCent],
                tMultHigh = tMacroConfig.centrality_ranges[iCent+1];

      tAnalysisConfig.minMult = 10.*tMultLow;
      tAnalysisConfig.maxMult = 10.*tMultHigh;
      tAnalysisConfig.nBinsMult = (tMultHigh-tMultLow)/5.0;

      tEventCutConfig.minCentrality = tMultLow;
      tEventCutConfig.maxCentrality = tMultHigh;

      //loop over pair types
      for(unsigned int iPair = 0; iPair < tMacroConfig.pair_codes.size(); iPair++)
      {
        
        //Update the configuration
        BuildConfiguration(tNewParamString,tAnalysisConfig,tEventCutConfig,tPairCutConfig,tMacroConfig,tCutVariations);

        //Build unique analysis for each pair type in each centrality bin
        //and update the particle cuts configurations
        AliFemtoAnalysisLambdaKaon *tAnalysis = CreateCorrectAnalysis(tNewParamString,tMacroConfig.pair_codes[iPair],tAnalysisConfig,tEventCutConfig,tPairCutConfig,tDirNameModifier);
        //TODO get pair cut to change for LamKchP and LamKchM

        tManager->AddAnalysis(tAnalysis);
      }

    }

  }

  return tManager;
}



void SetPairCodes(AFALK::AnalysisParams &aAnConfig, MacroParams &aMacroConfig)
{
  aMacroConfig.pair_codes.clear();

  switch(aAnConfig.analysisType) {
  case AFALK::kLamK0:
  case AFALK::kALamK0:
    aMacroConfig.pair_codes.push_back(AFALK::kLamK0);
    aMacroConfig.pair_codes.push_back(AFALK::kALamK0);
    aAnConfig.generalAnalysisType = AFALK::kV0V0;
    break;

  case AFALK::kLamKchP:
  case AFALK::kALamKchP:
  case AFALK::kLamKchM:
  case AFALK::kALamKchM:
    aMacroConfig.pair_codes.push_back(AFALK::kLamKchP);
    aMacroConfig.pair_codes.push_back(AFALK::kALamKchP);
    aMacroConfig.pair_codes.push_back(AFALK::kLamKchM);
    aMacroConfig.pair_codes.push_back(AFALK::kALamKchM);
    aAnConfig.generalAnalysisType = AFALK::kV0Track;
    break;

  case AFALK::kLamLam:
  case AFALK::kALamALam:
  case AFALK::kLamALam:
    aMacroConfig.pair_codes.push_back(AFALK::kLamLam);
    aMacroConfig.pair_codes.push_back(AFALK::kALamALam);
    aMacroConfig.pair_codes.push_back(AFALK::kLamALam);
    aAnConfig.generalAnalysisType = AFALK::kV0V0;
    break;

  case AFALK::kLamPiP:
  case AFALK::kALamPiP:
  case AFALK::kLamPiM:
  case AFALK::kALamPiM:
    aMacroConfig.pair_codes.push_back(AFALK::kLamPiP);
    aMacroConfig.pair_codes.push_back(AFALK::kALamPiP);
    aMacroConfig.pair_codes.push_back(AFALK::kLamPiM);
    aMacroConfig.pair_codes.push_back(AFALK::kALamPiM);
    aAnConfig.generalAnalysisType = AFALK::kV0Track;
    break;

  case AFALK::kXiKchP:
  case AFALK::kAXiKchP:
  case AFALK::kXiKchM:
  case AFALK::kAXiKchM:
    aMacroConfig.pair_codes.push_back(AFALK::kXiKchP);
    aMacroConfig.pair_codes.push_back(AFALK::kAXiKchP);
    aMacroConfig.pair_codes.push_back(AFALK::kXiKchM);
    aMacroConfig.pair_codes.push_back(AFALK::kAXiKchM);
    aAnConfig.generalAnalysisType = AFALK::kXiTrack;
    break;

  case AFALK::kXiK0:
  case AFALK::kAXiK0:
    aMacroConfig.pair_codes.push_back(AFALK::kXiK0);
    aMacroConfig.pair_codes.push_back(AFALK::kAXiK0);
    aAnConfig.generalAnalysisType = AFALK::kXiV0;
    break;

  case AFALK::kProtPiM:
    aMacroConfig.pair_codes.push_back(AFALK::kProtPiM);
    aAnConfig.generalAnalysisType = AFALK::kTrackTrack;
  break;

  case AFALK::kAProtPiP:
    aMacroConfig.pair_codes.push_back(AFALK::kAProtPiP);
    aAnConfig.generalAnalysisType = AFALK::kTrackTrack;
  break;

  case AFALK::kPiPPiM:
    aMacroConfig.pair_codes.push_back(AFALK::kPiPPiM);
    aAnConfig.generalAnalysisType = AFALK::kTrackTrack;
  break;

  default:
    break;
  }
}

AliFemtoAnalysisLambdaKaon* 
CreateCorrectAnalysis(
  const TString &aText,
  AFALK::AnalysisType aAnType,
  AFALK::AnalysisParams &aAnParams,
  AFALK::EventCutParams &aEvCutParams,
  AFALK::PairCutParams &aPairCutParams,
  const TString &aDirNameModifier
)
{
  AliFemtoAnalysisLambdaKaon *tAnalysis;

  aAnParams.analysisType = aAnType;

  AFALK::V0CutParams tV0CutConfig1,
                     tV0CutConfig2;
  AFALK::ESDCutParams tESDCutConfig1;
  AFALK::ESDCutParams tESDCutConfig2;
  AFALK::XiCutParams tXiCutConfig;

  switch(aAnParams.generalAnalysisType) {

  case AFALK::kV0V0:
    switch(aAnType) {
    case AFALK::kLamK0:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tV0CutConfig2 = AFALK::DefaultK0ShortCutParams();
      break;

    case AFALK::kALamK0:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tV0CutConfig2 = AFALK::DefaultK0ShortCutParams();
      break;

    case AFALK::kLamLam:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tV0CutConfig2 = AFALK::DefaultLambdaCutParams();
      break;

    case AFALK::kALamALam:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tV0CutConfig2 = AFALK::DefaultAntiLambdaCutParams();
      break;

    case AFALK::kLamALam:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tV0CutConfig2 = AFALK::DefaultAntiLambdaCutParams();
      break;
    }
    BuildParticleConfiguration(aText,tV0CutConfig1);
    BuildParticleConfiguration(aText,tV0CutConfig2);
    tAnalysis = new AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aPairCutParams,tV0CutConfig1,tV0CutConfig2,aDirNameModifier);
    break;


  case AFALK::kV0Track:
    switch(aAnType) {
    case AFALK::kLamKchP:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(1);
      break;

    case AFALK::kALamKchP:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(1);
      break;

    case AFALK::kLamKchM:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(-1);
      break;

    case AFALK::kALamKchM:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(-1);
      break;

    case AFALK::kLamPiP:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultPiCutParams(1);
      break;

    case AFALK::kALamPiP:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultPiCutParams(1);
      break;

    case AFALK::kLamPiM:
      tV0CutConfig1 = AFALK::DefaultLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultPiCutParams(-1);
      break;

    case AFALK::kALamPiM:
      tV0CutConfig1 = AFALK::DefaultAntiLambdaCutParams();
      tESDCutConfig1 = AFALK::DefaultPiCutParams(-1);
      break;
    }
    BuildParticleConfiguration(aText,tV0CutConfig1);
    BuildParticleConfiguration(aText,tESDCutConfig1);
    tAnalysis = new AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aPairCutParams,tV0CutConfig1,tESDCutConfig1,aDirNameModifier);
    break;


  case AFALK::kXiTrack:
    switch(aAnType) {
    case AFALK::kXiKchP:
      tXiCutConfig = AFALK::DefaultXiCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(1);
      break;

    case AFALK::kAXiKchP:
      tXiCutConfig = AFALK::DefaultAXiCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(1);
      break;

    case AFALK::kXiKchM:
      tXiCutConfig = AFALK::DefaultXiCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(-1);
      break;

    case AFALK::kAXiKchM:
      tXiCutConfig = AFALK::DefaultAXiCutParams();
      tESDCutConfig1 = AFALK::DefaultKchCutParams(-1);
      break;
    }
    BuildParticleConfiguration(aText,tXiCutConfig);
    BuildParticleConfiguration(aText,tESDCutConfig1);
    tAnalysis = new AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aPairCutParams,tXiCutConfig,tESDCutConfig1,aDirNameModifier);
    break;

  case AFALK::kXiV0:
    switch(aAnType) {
    case AFALK::kXiK0:
      tXiCutConfig = AFALK::DefaultXiCutParams();
      tV0CutConfig1 = AFALK::DefaultK0ShortCutParams();
      break;

    case AFALK::kAXiK0:
      tXiCutConfig = AFALK::DefaultAXiCutParams();
      tV0CutConfig1 = AFALK::DefaultK0ShortCutParams();
      break;
    }
    BuildParticleConfiguration(aText,tXiCutConfig);
    BuildParticleConfiguration(aText,tV0CutConfig1);
    tAnalysis = new AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aPairCutParams,tXiCutConfig,tV0CutConfig1);
    break;

  case AFALK::kTrackTrack:
    switch(aAnType) {
    case AFALK::kProtPiM:
      tESDCutConfig1 = AFALK::LambdaPurityProtonCutParams(1);
      tESDCutConfig2 = AFALK::LambdaPurityPiCutParams(-1);
      break;

    case AFALK::kAProtPiP:
      tESDCutConfig1 = AFALK::LambdaPurityProtonCutParams(-1);
      tESDCutConfig2 = AFALK::LambdaPurityPiCutParams(1);
      break;

    case AFALK::kPiPPiM:
      tESDCutConfig1 = AFALK::K0ShortPurityPiCutParams(1);
      tESDCutConfig2 = AFALK::K0ShortPurityPiCutParams(-1);
      break;
    }
    BuildParticleConfiguration(aText,tESDCutConfig1);
    BuildParticleConfiguration(aText,tESDCutConfig2);
    tAnalysis = new AliFemtoAnalysisLambdaKaon(aAnParams,aEvCutParams,aPairCutParams,tESDCutConfig1,tESDCutConfig2);
    break;
  }

  return tAnalysis;
}


void
BuildConfiguration(
  const TString &aText,
  AFALK::AnalysisParams &aAnParams,
  AFALK::EventCutParams &aEvCutParams,
  AFALK::PairCutParams &aPairCutParams,
  MacroParams &aMac,
  CutVariations &aCutVariations)
{
//  std::cout << "I-BuildConfiguration:" << TBase64::Encode(text) << " \n";

  const TString tAnalysisVarName = "aAnParams",//
                tEvCutVarName    = "aEvCutParams",
                tPairCutVarName  = "aPairCutParams",
                tMacroVarName    = "aMac",
                tCutVariationsName = "aCutVariations";

  TObjArray* tLines = aText.Tokenize("\n;");
  TIter tNextLine(tLines);
  TObject *tLineObj = NULL;

  while(tLineObj = tNextLine())
  {
    const TString tLine = ((TObjString*)tLineObj)->String().Strip(TString::kBoth, ' ');
    TString tCmd("");

    switch(tLine[0]) {
    case '@':  //Analysis Params
      tCmd = tAnalysisVarName + "." + tLine(1, tLine.Length() - 1);
      break;

    case '&':  //Event Cut Params
      tCmd = tEvCutVarName + "." + tLine(1, tLine.Length() - 1);
      break;

    case '%':  //Pair Cut Params
      tCmd = tPairCutVarName + "." + tLine(1, tLine.Length() - 1);
      break;

    case '~':  //Macro Params
      tCmd = tMacroVarName + "." + tLine(1, tLine.Length() - 1);
      break;

    case '[':  //Centrality in Macro Params
    {
      unsigned int tRangeEnd = tLine.Index("]");
      if(tRangeEnd == -1) tRangeEnd = tLine.Length();

      TString tCentralityRanges = tLine(1, tRangeEnd - 1);
      TObjArray *tRangeGroups = tCentralityRanges.Tokenize(",");
      TIter tNextRangeGroup(tRangeGroups);
      TObjString *tRangeGroup = NULL;

      while(tRangeGroup = (TObjString*)tNextRangeGroup())
      {
        TObjArray *tSubRange = tRangeGroup->String().Tokenize(":");
        TIter tNextSubRange(tSubRange);
        TObjString *tSubRange_it = (TObjString*)tNextSubRange();
        TString tPrev = TString::Format("%0.2d", tSubRange_it->String().Atoi());
        while(tSubRange_it = (TObjString*)tNextSubRange())
        {
          TString tNext = TString::Format("%0.2d", tSubRange_it->String().Atoi());

          tCmd = tMacroVarName + ".centrality_ranges.push_back(" + tPrev + ");";
          cout << "I-BuildConfiguration: `" << tCmd << "`\n";
          gROOT->ProcessLineFast(tCmd);

          tCmd = tMacroVarName + ".centrality_ranges.push_back(" + tNext + ");";
          cout << "I-BuildConfiguration: `" << tCmd << "`\n";
          gROOT->ProcessLineFast(tCmd);
          tPrev = tNext;
        }
      }
    }
      continue;

    case '?':  //Centrality in Macro Params
    {
      if(!(tLine[1] == '[')) 
      {
        tCmd = tCutVariationsName + ".cutName = '" + tLine(1, tLine.Length() - 1) + "'";
        tCmd.ReplaceAll("'", '"');
      }
      else
      {
        unsigned int tRangeEnd = tLine.Index("]");
        if(tRangeEnd == -1) tRangeEnd = tLine.Length();

        TString tCentralityRanges = tLine(2, tRangeEnd - 1);
        TObjArray *tRangeGroups = tCentralityRanges.Tokenize(",");
        TIter tNextRangeGroup(tRangeGroups);
        TObjString *tRangeGroup = NULL;

        while(tRangeGroup = (TObjString*)tNextRangeGroup())
        {
          TObjArray *tSubRange = tRangeGroup->String().Tokenize(":");
          TIter tNextSubRange(tSubRange);
          TObjString *tSubRange_it = NULL;
          while(tSubRange_it = (TObjString*)tNextSubRange())
          {
            TString tNext = TString::Format("%f", tSubRange_it->String().Atof());

            tCmd = tCutVariationsName + ".valuesToVary.push_back(" + tNext + ");";
            cout << "I-BuildConfiguration: `" << tCmd << "`\n";
            gROOT->ProcessLineFast(tCmd);

          }
        }
        continue;
      }
    }
      break;

    default:
      continue;
    }

    if(!tCmd.IsNull())
    {
      tCmd += ";";
      cout << "I-BuildConfiguration: `" << tCmd << "`\n";
      gROOT->ProcessLineFast(tCmd);
    }

  }
}

void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::V0CutParams &aV0CutParams
)
{
//  std::cout << "I-BuildParticleConfiguration:" << TBase64::Encode(text) << " \n";

  const TString tV0CutVarName = "aV0CutParams";

  TString tDesiredName;

  switch(aV0CutParams.particlePDGType) {
  case AFALK::kPDGLam:
    tDesiredName = TString("Lam");
    break;

  case AFALK::kPDGALam:
    tDesiredName = TString("ALam");
    break;

  case AFALK::kPDGK0:
    tDesiredName = TString("K0s");
    break;

  default:
    tDesiredName = TString("");
  }

  TObjArray* tLines = aText.Tokenize("\n;");
  TIter tNextLine(tLines);
  TObject *tLineObj = NULL;

  while(tLineObj = tNextLine())
  {
    const TString tLine = ((TObjString*)tLineObj)->String().Strip(TString::kBoth, ' ');
    TString tCmd("");

    if(tLine[0] == '$')
    {
      TObjArray* tCutFullLine = tLine.Tokenize("|");
      const TString tParticleType = ((TObjString*)tCutFullLine->At(1))->String().Strip(TString::kBoth, ' ');
      const TString tParticleCut = ((TObjString*)tCutFullLine->At(2))->String().Strip(TString::kBoth, ' ');

      if(tParticleType.EqualTo(tDesiredName) || tParticleType.EqualTo("ALL") || tParticleType.EqualTo("ALLV0S")) tCmd = tV0CutVarName + "." + tParticleCut(0, tParticleCut.Length());

      if(tParticleType.EqualTo("CLAM"))  //do for both Lam and ALam
      {
        if(aV0CutParams.particlePDGType==AFALK::kPDGLam || aV0CutParams.particlePDGType==AFALK::kPDGALam) tCmd = tV0CutVarName + "." + tParticleCut(0, tParticleCut.Length());
      }
    }

    if(!tCmd.IsNull())
    {
      tCmd += ";";
      cout << "I-BuildParticleConfiguration: `" << tCmd << "`\n";
      gROOT->ProcessLineFast(tCmd);
    }

  }
}


void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::ESDCutParams &aESDCutParams
)
{
//  std::cout << "I-BuildParticleConfiguration:" << TBase64::Encode(text) << " \n";

  const TString tESDCutVarName = "aESDCutParams";

  TString tDesiredName;

  switch(aESDCutParams.particlePDGType) {
  case AFALK::kPDGKchP:
    tDesiredName = TString("KchP");
    break;

  case AFALK::kPDGKchM:
    tDesiredName = TString("KchP");
    break;

  case AFALK::kPDGPiP:
    tDesiredName = TString("PiP");
    break;

  case AFALK::kPDGPiM:
    tDesiredName = TString("PiM");
    break;

  case AFALK::kPDGProt:
    tDesiredName = TString("Prot");
    break;

  case AFALK::kPDGAntiProt:
    tDesiredName = TString("AProt");
    break;

  default:
    tDesiredName = TString("");
  }

  TObjArray* tLines = aText.Tokenize("\n;");
  TIter tNextLine(tLines);
  TObject *tLineObj = NULL;

  while(tLineObj = tNextLine())
  {
    const TString tLine = ((TObjString*)tLineObj)->String().Strip(TString::kBoth, ' ');
    TString tCmd("");

    if(tLine[0] == '$')
    {
      TObjArray* tCutFullLine = tLine.Tokenize("|");
      const TString tParticleType = ((TObjString*)tCutFullLine->At(1))->String().Strip(TString::kBoth, ' ');
      const TString tParticleCut = ((TObjString*)tCutFullLine->At(2))->String().Strip(TString::kBoth, ' ');

      if(tParticleType.EqualTo(tDesiredName) || tParticleType.EqualTo("ALL") || tParticleType.EqualTo("ALLTRACKS")) tCmd = tESDCutVarName + "." + tParticleCut(0, tParticleCut.Length());
    }

    if(!tCmd.IsNull())
    {   
      tCmd += ";";
      cout << "I-BuildParticleConfiguration: `" << tCmd << "`\n";
      gROOT->ProcessLineFast(tCmd);
    }

  }
}


void
BuildParticleConfiguration(
  const TString &aText,
  AFALK::XiCutParams &aXiCutParams
)
{
//  std::cout << "I-BuildParticleConfiguration:" << TBase64::Encode(text) << " \n";

  const TString tXiCutVarName = "aXiCutParams";

  TString tDesiredName;

  switch(aXiCutParams.particlePDGType) {
  case AFALK::kPDGXiC:
    tDesiredName = TString("Xi");
    break;

  case AFALK::kPDGAXiC:
    tDesiredName = TString("AXi");
    break;

  default:
    tDesiredName = TString("");
  }

  TObjArray* tLines = aText.Tokenize("\n;");
  TIter tNextLine(tLines);
  TObject *tLineObj = NULL;

  while(tLineObj = tNextLine())
  {
    const TString tLine = ((TObjString*)tLineObj)->String().Strip(TString::kBoth, ' ');
    TString tCmd("");

    if(tLine[0] == '$')
    {
      TObjArray* tCutFullLine = tLine.Tokenize("|");
      const TString tParticleType = ((TObjString*)tCutFullLine->At(1))->String().Strip(TString::kBoth, ' ');
      const TString tParticleCut = ((TObjString*)tCutFullLine->At(2))->String().Strip(TString::kBoth, ' ');

      if(tParticleType.EqualTo(tDesiredName) || tParticleType.EqualTo("ALL") || tParticleType.EqualTo("ALLXIS")) tCmd = tXiCutVarName + "." + tParticleCut(0, tParticleCut.Length());
    }

    if(!tCmd.IsNull())
    {
      tCmd += ";";
      cout << "I-BuildParticleConfiguration: `" << tCmd << "`\n";
      gROOT->ProcessLineFast(tCmd);
    }

  }
}


TString GetParticleCuts(const TString &aText)
{
  TObjArray* tLines = aText.Tokenize("\n;");
  TIter tNextLine(tLines);
  TObject *tLineObj = NULL;

  TString tReturnString("");

  while(tLineObj = tNextLine())
  {
    const TString tLine = ((TObjString*)tLineObj)->String().Strip(TString::kBoth, ' ');
    TString tCmd("");

    if(tLine[0] == '$')
    {
      tReturnString += ((TObjString*)tLineObj)->String();
      tReturnString += TString("; ");
    }
  }
  return tReturnString;
}

