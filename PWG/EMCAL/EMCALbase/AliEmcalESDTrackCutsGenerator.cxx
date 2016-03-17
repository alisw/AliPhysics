/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TString.h>
#include <TFormula.h>

#include "AliEmcalTrackSelection.h"
#include "AliEmcalESDTrackCutsGenerator.h"

#include "AliESDtrackCuts.h"

#include "TError.h"

const Int_t AliEmcalESDTrackCutsGenerator::fgkAddCutFactor = 10000;

/**
 * Function to create track cuts for PWG Jet analysis
 * User can select a specific set by indicating cutMode
 *
 * \param cutMode has 8 digits: first 4 digits additional cuts, last 4 digits standard cuts
 *                 additional cuts are variations of standard cuts (used for hybrid track selection and QA)
 *                 Numbering starts from 1000 For standard and additional cut numbers
 *
 * \return AliESDtrackCuts object
 */
AliESDtrackCuts* AliEmcalESDTrackCutsGenerator::CreateTrackCutsPWGJE(Int_t cutMode)
{
  //Get standard cuts: last 4 digits of cutMode
  Int_t stdCutMode = cutMode % fgkAddCutFactor;

  //Get additional cut mode: first 4 digits of cutMode
  Int_t addCutMode = (int)((float)cutMode/(float)fgkAddCutFactor);

  return CreateTrackCutsPWGJE(stdCutMode, addCutMode);
}

/**
 * Function to create track cuts for PWG Jet analysis
 * User can select a specific set by indicating cutMode
 *
 * \param stdCutMode standard cuts, should be one of the EStdCutMode_t enum values
 * \param addCutMode additional cuts, should be one of the EAddCutMode_t enum values
 *
 * \return AliESDtrackCuts object
 */
AliESDtrackCuts* AliEmcalESDTrackCutsGenerator::CreateTrackCutsPWGJE(Int_t stdCutMode, Int_t addCutMode)
{
  AliESDtrackCuts* trackCuts = 0;

  TString tag;

  tag = SetStandardCuts(trackCuts, stdCutMode);
  tag += SetAdditionalCuts(trackCuts, addCutMode);

  ::Info("AliEmcalESDTrackCutsGenerator::CreateTrackCutsPWGJE", "Created track cuts for: %s", tag.Data());

  return trackCuts;
}

/**
 * Function to create track cuts for PWG Jet analysis
 * User can select a specific set by indicating cutMode
 *
 * \param stdCutMode standard cuts, should be one of the EStdCutMode_t enum values
 * \param addCutMode1 additional cuts, should be one of the EAddCutMode_t enum values
 * \param addCutMode2 additional cuts, should be one of the EAddCutMode_t enum values
 *
 * \return AliESDtrackCuts object
 */
AliESDtrackCuts* AliEmcalESDTrackCutsGenerator::CreateTrackCutsPWGJE(Int_t stdCutMode, Int_t addCutMode1, Int_t addCutMode2)
{
  AliESDtrackCuts* trackCuts = 0;

  TString tag;

  tag = SetStandardCuts(trackCuts, stdCutMode);
  tag += SetAdditionalCuts(trackCuts, addCutMode1);
  tag += SetAdditionalCuts(trackCuts, addCutMode2);

  ::Info("AliEmcalESDTrackCutsGenerator::CreateTrackCutsPWGJE", "Created track cuts for: %s", tag.Data());

  return trackCuts;
}

/**
 * Function to set standard track cuts
 * User can select a specific set by indicating stdCutMode
 *
 * \param trackCuts AliESDtrackCuts to be set
 * \param stdCutMode has 4 digits, >= 1000
 *
 * \return string with track cut label
 */
TString AliEmcalESDTrackCutsGenerator::SetStandardCuts(AliESDtrackCuts*& trackCuts, Int_t stdCutMode)
{
  TString tag;

  if (trackCuts) {
    delete trackCuts;
    trackCuts = 0;
  }

  switch (stdCutMode) {
  case kRAA2011:
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    trackCuts->SetMaxChi2TPCConstrainedGlobal(36);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1e10);

    tag = "Global track RAA analysis QM2011 + Chi2ITS<36";

    break;
  }

  case kGlobalTracksNCls90NoSPD:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    trackCuts->SetMinNClustersTPC(90);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=90, noSPD requirement";

    break;
  }

  case kGlobalTracksNCls80NoSPD:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    trackCuts->SetMinNClustersTPC(80);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

    tag = "Global tracks jet analysis with ITSrefit and Ncls=80, noSPD requirement";

    break;
  }

  case kGlobalTracks2010NCrossRows120:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // tight global tracks
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,1);
    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially switches it off
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);

    tag = "Global tracks ITSTPC2010 + NCrossedRows + loose ITS";

    break;
  }

  case kGlobalTracksNCls70NoSPD:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=70, noSPD requirement";

    break;
  }

  case kGlobalTracksNCls70NoSPDNoPtCut:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1E+15);

    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=70, noSPD requirement, no upper pt cut";

    break;
  }

  case kGlobalTracksNClsPtDepNoSPDNoPtCut:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
    trackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
    trackCuts->SetMinNClustersTPC(70);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackCuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackCuts->SetMaxChi2PerClusterITS(36);
    trackCuts->SetMaxChi2TPCConstrainedGlobal(36);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1E+15);

    tag = "Global tracks jet analysis with ITSrefit and NclsIter1=PtDep, noSPD requirement, no upper pt cut, golden chi2";

    break;
  }

  case kGlobalTracks2011:
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);

    trackCuts->SetMaxChi2TPCConstrainedGlobal(36);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 1E+15);

    tag = "Global tracks with AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE)";

    break;
  }

  case kGlobalTracks2011NoSPD:
  {
    trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);

    trackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

    trackCuts->SetMaxFractionSharedTPCClusters(0.4);

    tag = "Global tracks 2011 with AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE) and no SPD requirement";

    break;
  }

  case kGlobalTracksNCls90NoITS:
  {
    trackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    // TPC
    trackCuts->SetMinNClustersTPC(90);
    trackCuts->SetMaxChi2PerClusterTPC(4);
    trackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackCuts->SetAcceptKinkDaughters(kFALSE);
    trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMaxFractionSharedTPCClusters(0.4);
    //accept secondaries
    trackCuts->SetMaxDCAToVertexXY(2.4);
    trackCuts->SetMaxDCAToVertexZ(3.2);
    trackCuts->SetDCAToVertex2D(kTRUE);

    trackCuts->SetRequireSigmaToVertex(kFALSE);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

    tag = "Global tracks jet analysis, loose cuts, NClsIter1=90, no ITS requirements";

    break;
  }

  case kTPCOnlyTracksNCls70:
  {
    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    // trackCuts->SetRequireTPCRefit(kTRUE);
    trackCuts->SetMinNClustersTPC(70);

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);


    tag = "TPConly track cuts, loose cuts, NCls=70, no ITS requirements";

    break;
  }

  case kTPCOnlyTracksNCrossRows120:
  {
    trackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(120);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.1);// essentially switches it off

    trackCuts->SetEtaRange(-0.9,0.9);
    trackCuts->SetPtRange(0.15, 100.);

    tag = "TPConly track cuts, loose cuts, NCrossRows=120, no ITS requirements";

    break;
  }

  default:
  {
    Printf("AliEmcalESDTrackCutsGenerator: standard cuts not recognized.");
    break;
  }
  }

  return tag;
}

/**
 * Function to set additional track cuts on top of the standard ones
 * User can select a specific set by indicating addCutMode
 *
 * \param trackCuts AliESDtrackCuts to be set
 * \param addCutMode has 4 digits, >= 1000
 *
 * \return string with track cut label
 */
TString AliEmcalESDTrackCutsGenerator::SetAdditionalCuts(AliESDtrackCuts*& trackCuts, Int_t addCutMode)
{
  TString tag;

  switch (addCutMode) {
  case kSPDAny:
  {
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

    tag += " + additonal: SPD any requirement";

    break;
  }

  case kSPDNone:
  {
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

    tag += " + additional: w/o hits in SPD";

    break;
  }

  case kNoITSChi2:
  {
    trackCuts->SetMaxChi2PerClusterITS(1E10);

    tag += " + additional: maxITSChi2=1e10";

    break;
  }

  case kNoMinTPCCls:
  {
    trackCuts->SetMinNClustersTPC(0);
    trackCuts->SetMinNCrossedRowsTPC(0);
    trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.);

    tag += " + additional: minClusters=0 minCrossedRows=0 minCrossedRowsOverFindable=0";

    break;
  }

  case kNoITSRefit:
  {
    trackCuts->SetRequireITSRefit(kFALSE);

    tag += " + additional: ITSrefit=kFALSE";

    break;
  }

  case kSPDOff:
  {
    trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

    tag += " + additional: no SPD requirement (kOff)";

    break;
  }
  }

  return tag;
}

/**
 * Helper function to steer period string into an enum type
 *
 * \param period String identifying a data or MC period
 *
 * \return enum type identifying a data or MC period
 */
AliEmcalESDTrackCutsGenerator::EDataSet_t AliEmcalESDTrackCutsGenerator::SteerDataSetFromString(TString period)
{
  EDataSet_t dataSet = kUnknown;

  TString strPeriod(period);
  strPeriod.ToLower();

  if (strPeriod == "lhc10h") {
    dataSet = kLHC10h;
  } else if (strPeriod == "lhc11a" || strPeriod == "lhc12a15a") {
    dataSet = kLHC11a;
  } else if (strPeriod == "lhc10b" || strPeriod == "lhc10c" ||
      strPeriod == "lhc10d" || strPeriod == "lhc10e") {
    dataSet = kLHC10bcde;
  } else if (strPeriod == "lhc11a1a" ||  strPeriod == "lhc11a1b" ||
      strPeriod == "lhc11a1c" ||  strPeriod == "lhc11a1d" ||
      strPeriod == "lhc11a1e" ||  strPeriod == "lhc11a1f" ||
      strPeriod == "lhc11a1g" ||  strPeriod == "lhc11a1h" ||
      strPeriod == "lhc11a1i" ||  strPeriod == "lhc11a1j") {
    dataSet = kLHC11a;
  } else if (strPeriod == "lhc11c") {
    dataSet = kLHC11c;
  } else if (strPeriod == "lhc11d") {
    dataSet = kLHC11d;
  } else if (strPeriod == "lhc11h" || strPeriod == "lhc12a15e") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc12g") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc12") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13b") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13c") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13d") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13e") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13f") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13g") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc12a15f") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc13b4") {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc12a15g") {
    dataSet = kLHC11d;
  } else if (strPeriod == "lhc12f2a") {
    dataSet = kLHC11d;
  } else if (strPeriod.BeginsWith("lhc12a17")) {
    dataSet = kLHC11h;
  } else if (strPeriod == "lhc14a1") {
    dataSet = kLHC11h;
  } else if (strPeriod.BeginsWith("lhc15g6")) {
    dataSet = kLHC10bcde;
  } else {
    ::Error("AliEmcalESDTrackCutsGenerator::SteerDataSetFromString", "Dataset %s not recognized!", period.Data());
  }

  return dataSet;
}

/**
 * Adds hybrid track cuts to an AliEmcalTrackSelection object
 *
 * \param trkSel AliEmcalTrackSelection object
 * \param period enum type identifying a data or MC period
 */
void AliEmcalESDTrackCutsGenerator::AddHybridTrackCuts(AliEmcalTrackSelection* trkSel, EDataSet_t period)
{
  switch (period) {
  case kLHC11c:
  case kLHC11d:
  case kLHC11h:
  {
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(kGlobalTracks2011NoSPD, kSPDAny);
    trkSel->AddTrackCuts(cutsp);
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(kGlobalTracks2011NoSPD, kSPDOff);
    trkSel->AddTrackCuts(hybsp);

    break;
  }
  case kLHC10h:
  case kLHC11a:
  case kLHC10bcde:
  {
    /* hybrid track cuts*/
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(kGlobalTracksNClsPtDepNoSPDNoPtCut, kSPDAny);
    trkSel->AddTrackCuts(cutsp);
    AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(kGlobalTracksNClsPtDepNoSPDNoPtCut, kSPDOff, kNoITSRefit);
    trkSel->AddTrackCuts(hybsp);
    break;
  }
  default:
  {
    ::Error("AliEmcalESDTrackCutsGenerator::AddHybridTrackCuts", "Hybrid track cuts not available for dataset %d", period);
    break;
  }
  }
}

/**
 * Adds TPC only track cuts to an AliEmcalTrackSelection object
 *
 * \param trkSel AliEmcalTrackSelection object
 * \param period enum type identifying a data or MC period
 */
void AliEmcalESDTrackCutsGenerator::AddTPCOnlyTrackCuts(AliEmcalTrackSelection* trkSel, EDataSet_t period)
{
  switch (period) {
  case kLHC11c:
  case kLHC11d:
  case kLHC11h:
  {
    AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(kTPCOnlyTracksNCls70);
    trkSel->AddTrackCuts(cutsp);
    break;
  }
  default:
  {
    Printf("AliEmcalESDTrackCutsGenerator::AddTPCOnlyTrackCuts: TPC only track cuts not available for dataset %d", period);
    break;
  }
  }
}
