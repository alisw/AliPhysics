/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliRDHFCutsBPlustoD0Pi.cxx $ */

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed BPlus->D0Pi->KPiPi
//
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//     Several AliPhysics classes have been used as a basis for this code
//
//
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsBPlustoD0Pi.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsBPlustoD0Pi);
/// \endcond



//--------------------------------------------------------------------------
AliRDHFCutsBPlustoD0Pi::AliRDHFCutsBPlustoD0Pi(const char* name) :
  AliRDHFCuts(name),
  fTrackCutsSoftPi(0),
  fMaxPtPid(9999.),
  fTPCflag(999.),
  fCircRadius(0.),

  fIsCutUsed(0x0),

  fnVarsD0forD0ptbin(0),
  fnPtBinsD0forD0ptbin(1),
  fGlobalIndexD0forD0ptbin(0),
  fCutsRDD0forD0ptbin(0x0),
  fnPtBinLimitsD0forD0ptbin(1),
  fPtBinLimitsD0forD0ptbin(0x0),
  fIsUpperCutD0forD0ptbin(0x0),
  fIsCutUsedD0forD0ptbin(0x0),
  fVarNamesD0forD0ptbin(0x0),

  fMinITSNclsD0FirstDaughter(0),
  fMinTPCNclsD0FirstDaughter(0),
  fUseITSRefitD0FirstDaughter(0),
  fUseTPCRefitD0FirstDaughter(0),
  fUseFilterBitD0FirstDaughter(0),
  fFilterBitD0FirstDaughter(0),
  fMinPtD0FirstDaughter(0),
  fMaxAbsEtaD0FirstDaughter(999.),
  fHardSelectionArrayITSD0FirstDaughter(),
  fSoftSelectionArrayITSD0FirstDaughter(),
  fNSoftITSCutD0FirstDaughter(0),

  fMinITSNclsD0SecondDaughter(0),
  fMinTPCNclsD0SecondDaughter(0),
  fUseITSRefitD0SecondDaughter(0),
  fUseTPCRefitD0SecondDaughter(0),
  fUseFilterBitD0SecondDaughter(0),
  fFilterBitD0SecondDaughter(0),
  fMinPtD0SecondDaughter(0),
  fMaxAbsEtaD0SecondDaughter(999.),
  fHardSelectionArrayITSD0SecondDaughter(),
  fSoftSelectionArrayITSD0SecondDaughter(),
  fNSoftITSCutD0SecondDaughter(0),

  fMinITSNclsBPlusPion(0),
  fMinTPCNclsBPlusPion(0),
  fUseITSRefitBPlusPion(0),
  fUseTPCRefitBPlusPion(0),
  fUseFilterBitBPlusPion(0),
  fFilterBitBPlusPion(0),
  fMinPtBPlusPion(0),
  fMaxAbsEtaBPlusPion(999.),
  fHardSelectionArrayITSBPlusPion(),
  fSoftSelectionArrayITSBPlusPion(),
  fNSoftITSCutBPlusPion(0),

  fMind0D0FirstDaughter(0),
  fMind0D0SecondDaughter(0),
  fMind0BPlusPion(0.),
  fFiducialYCut(0.8),

  fnVariablesForCutOptimization(0),
  fnCutsForOptimization(0),
  fGlobalIndexCutOptimization(0),
  fCutsRDForCutOptimization(0x0),
  fIsUpperCutForCutOptimization(0x0),
  fCutIndexForCutOptimization(0x0),
  fSigmaForCutOptimization(0x0),
  fNumberOfSigmaBinsForCutOptimization(0)
{
  //
  // Default Constructor
  //

  // Main cut setup as function of BPlus pt bins
  const Int_t nvars = 75;
  SetNVars(nvars);

  TString varNames[nvars];
  Int_t iterator = 0;

  //D0 cut variables
  varNames[iterator++] =   /*-00-*/ "inv. mass width[GeV]";
  varNames[iterator++] =   /*-01-*/ "delta mass width  [GeV]";   //not used for D0
  varNames[iterator++] =   /*-02-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++] =   /*-03-*/ "dca [cm]";
  varNames[iterator++] =   /*-04-*/ "Pt D0 [GeV/c]";
  varNames[iterator++] =   /*-05-*/ "Pt first daughter [GeV/c]";
  varNames[iterator++] =   /*-06-*/ "Pt second daughter [GeV/c]";
  varNames[iterator++] =   /*-07-*/ "d0 D0 [cm]";
  varNames[iterator++] =   /*-08-*/ "d0 first daughter [cm]";
  varNames[iterator++] =   /*-09-*/ "d0 second daughter [cm]";
  varNames[iterator++] =   /*-10-*/ "d0d0 [cm^2]";
  varNames[iterator++] =   /*-11-*/ "d0d0 XY [cm^2]";

  varNames[iterator++] =   /*-12-*/ "angle between both daughters";
  varNames[iterator++] =   /*-13-*/ "angle mother with first daughter";
  varNames[iterator++] =   /*-14-*/ "angle mother with second daughter";
  varNames[iterator++] =   /*-15-*/ "cosThetaStar";
  varNames[iterator++] =   /*-16-*/ "vertexDistance";
  varNames[iterator++] =   /*-17-*/ "pseudoProperDecayTime";
  varNames[iterator++] =   /*-18-*/ "DecayTime";
  varNames[iterator++] =   /*-19-*/ "normalizedDecayTime";
  varNames[iterator++] =   /*-20-*/ "normDecayLength";
  varNames[iterator++] =   /*-21-*/ "topomatic first daughter";
  varNames[iterator++] =   /*-22-*/ "topomatic second daughter";
  varNames[iterator++] =   /*-23-*/ "topomatic max";
  varNames[iterator++] =   /*-24-*/ "topomatic min";

  varNames[iterator++] =   /*-25-*/ "pointing angle XY [Cos(theta)]";
  varNames[iterator++] =   /*-26-*/ "vertex distance XY [cm]";
  varNames[iterator++] =   /*-27-*/ "normDecayLength XY";
  varNames[iterator++] =   /*-28-*/ "Chi2 per NDF vertex";

  varNames[iterator++] =   /*-29-*/ "pointingAngleToBPlus";
  varNames[iterator++] =   /*-30-*/ "d0MotherToBPlus";
  varNames[iterator++] =   /*-31-*/ "d0FirstDaughterToBPlus";
  varNames[iterator++] =   /*-32-*/ "d0SecondDaughterToBPlus";
  varNames[iterator++] =   /*-33-*/ "impactProductToBPlus";
  varNames[iterator++] =   /*-34-*/ "impactProductXYToBPlus";
  varNames[iterator++] =   /*-35-*/ "normDecayLengthToBPlus";
  varNames[iterator++] =   /*-36-*/ "pseudoProperDecayTimeToBPlus";
  varNames[iterator++] =   /*-37-*/ "DecayTimeToBPlus";
  varNames[iterator++] =   /*-38-*/ "normalizedDecayTimeToBPlus";

  //BPlus cut variables
  varNames[iterator++] =   /*-39-*/ "inv. mass width[GeV]";
  varNames[iterator++] =   /*-40-*/ "delta mass width  [GeV]";
  varNames[iterator++] =   /*-41-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++] =   /*-42-*/ "dca [cm]";
  varNames[iterator++] =   /*-43-*/ "Pt BPlus [GeV/c]";
  varNames[iterator++] =   /*-44-*/ "Pt D0 [GeV/c]";
  varNames[iterator++] =   /*-45-*/ "Pt Pion [GeV/c]";
  varNames[iterator++] =   /*-46-*/ "d0 BPlus [cm]";
  varNames[iterator++] =   /*-47-*/ "d0 D0 [cm]";
  varNames[iterator++] =   /*-48-*/ "d0 Pion [cm]";
  varNames[iterator++] =   /*-49-*/ "d0d0 [cm^2]";
  varNames[iterator++] =   /*-50-*/ "d0d0 XY [cm^2]";

  varNames[iterator++] =   /*-51-*/ "angle between both daughters";
  varNames[iterator++] =   /*-52-*/ "angle mother with first daughter";
  varNames[iterator++] =   /*-53-*/ "angle mother with second daughter";
  varNames[iterator++] =   /*-54-*/ "cosThetaStar";
  varNames[iterator++] =   /*-55-*/ "vertexDistance";
  varNames[iterator++] =   /*-56-*/ "pseudoProperDecayTime";
  varNames[iterator++] =   /*-57-*/ "DecayTime";
  varNames[iterator++] =   /*-58-*/ "normalizedDecayTime";
  varNames[iterator++] =   /*-59-*/ "normDecayLength";
  varNames[iterator++] =   /*-60-*/ "topomatic first daughter";
  varNames[iterator++] =   /*-61-*/ "topomatic second daughter";
  varNames[iterator++] =   /*-62-*/ "topomatic max";
  varNames[iterator++] =   /*-63-*/ "topomatic min";

  varNames[iterator++] =   /*-64-*/ "pointing angle XY [Cos(theta)]";
  varNames[iterator++] =   /*-65-*/ "vertex distance XY [cm]";
  varNames[iterator++] =   /*-66-*/ "normDecayLength XY";
  varNames[iterator++] =   /*-67-*/ "Chi2 per NDF vertex";

  varNames[iterator++] =   /*-68-*/ "Normalized d0 D0 first daughter [cm]";
  varNames[iterator++] =   /*-69-*/ "Normalized d0 D0 second daughter [cm]";
  varNames[iterator++] =   /*-70-*/ "Normalized d0 D0 [cm]";
  varNames[iterator++] =   /*-71-*/ "Normalized d0 BPlus pion [cm]";
  varNames[iterator++] =   /*-72-*/ "Normalized d0 BPlus [cm]";
  varNames[iterator++] =   /*-73-*/ "Normalized impact product D0 [cm]";
  varNames[iterator++] =   /*-74-*/ "Normalized impact product BPlus [cm]";

  Bool_t isUpperCut[nvars] = {0};

  SetVarNames(nvars, varNames, isUpperCut);

  Float_t limits[2] = {0, 999999999.};
  SetPtBins(2, limits);


  //
  // Initialization of D0 cuts for D0 pt bins
  //

  const Int_t nvarsD0forD0ptbin = 29;
  SetNVarsD0forD0ptbin(nvarsD0forD0ptbin);

  TString varNamesD0forD0ptbin[nvarsD0forD0ptbin];
  iterator = 0;

  //D0 cut variables
  varNamesD0forD0ptbin[iterator++] =   /*-00-*/ "inv. mass width[GeV]";
  varNamesD0forD0ptbin[iterator++] =   /*-01-*/ "delta mass width  [GeV]";   //not used for D0
  varNamesD0forD0ptbin[iterator++] =   /*-02-*/ "pointing angle [Cos(theta)]";
  varNamesD0forD0ptbin[iterator++] =   /*-03-*/ "dca [cm]";
  varNamesD0forD0ptbin[iterator++] =   /*-04-*/ "Pt D0 [GeV/c]";
  varNamesD0forD0ptbin[iterator++] =   /*-05-*/ "Pt first daughter [GeV/c]";
  varNamesD0forD0ptbin[iterator++] =   /*-06-*/ "Pt second daughter [GeV/c]";
  varNamesD0forD0ptbin[iterator++] =   /*-07-*/ "d0 D0 [cm]";
  varNamesD0forD0ptbin[iterator++] =   /*-08-*/ "d0 first daughter [cm]";
  varNamesD0forD0ptbin[iterator++] =   /*-09-*/ "d0 second daughter [cm]";
  varNamesD0forD0ptbin[iterator++] =   /*-10-*/ "d0d0 [cm^2]";
  varNamesD0forD0ptbin[iterator++] =   /*-11-*/ "d0d0 XY [cm^2]";

  varNamesD0forD0ptbin[iterator++] =   /*-12-*/ "angle between both daughters";
  varNamesD0forD0ptbin[iterator++] =   /*-13-*/ "angle mother with first daughter";
  varNamesD0forD0ptbin[iterator++] =   /*-14-*/ "angle mother with second daughter";
  varNamesD0forD0ptbin[iterator++] =   /*-15-*/ "cosThetaStar";
  varNamesD0forD0ptbin[iterator++] =   /*-16-*/ "vertexDistance";
  varNamesD0forD0ptbin[iterator++] =   /*-17-*/ "pseudoProperDecayTime";
  varNamesD0forD0ptbin[iterator++] =   /*-18-*/ "DecayTime";
  varNamesD0forD0ptbin[iterator++] =   /*-19-*/ "normalizedDecayTime";
  varNamesD0forD0ptbin[iterator++] =   /*-20-*/ "normDecayLength";
  varNamesD0forD0ptbin[iterator++] =   /*-21-*/ "topomatic first daughter";
  varNamesD0forD0ptbin[iterator++] =   /*-22-*/ "topomatic second daughter";
  varNamesD0forD0ptbin[iterator++] =   /*-23-*/ "topomatic max";
  varNamesD0forD0ptbin[iterator++] =   /*-24-*/ "topomatic min";

  varNamesD0forD0ptbin[iterator++] =   /*-25-*/ "pointing angle XY [Cos(theta)]";
  varNamesD0forD0ptbin[iterator++] =   /*-26-*/ "vertex distance XY [cm]";
  varNamesD0forD0ptbin[iterator++] =   /*-27-*/ "normDecayLength XY";
  varNamesD0forD0ptbin[iterator++] =   /*-28-*/ "Chi2 per NDF vertex";

  Bool_t isUpperCutD0forD0ptbin[nvarsD0forD0ptbin] = {0};

  SetVarNamesD0forD0ptbin(nvarsD0forD0ptbin, varNamesD0forD0ptbin, isUpperCutD0forD0ptbin);

  Float_t limitsD0forD0ptbin[2] = {0, 999999999.};
  SetPtBinsD0forD0ptbin(2, limitsD0forD0ptbin);

  Bool_t forOpt[16] = {0}; //not yet used for BPlus analysis
  SetVarsForOpt(16, forOpt);

}
//--------------------------------------------------------------------------
AliRDHFCutsBPlustoD0Pi::AliRDHFCutsBPlustoD0Pi(const AliRDHFCutsBPlustoD0Pi &source) :
  AliRDHFCuts(source),
  fTrackCutsSoftPi(0),
  fMaxPtPid(source.fMaxPtPid),
  fTPCflag(source.fTPCflag),
  fCircRadius(source.fCircRadius),

  fIsCutUsed(0x0),

  fnVarsD0forD0ptbin(source.fnVarsD0forD0ptbin),
  fnPtBinsD0forD0ptbin(source.fnPtBinsD0forD0ptbin),
  fGlobalIndexD0forD0ptbin(source.fGlobalIndexD0forD0ptbin),
  fCutsRDD0forD0ptbin(0x0),
  fnPtBinLimitsD0forD0ptbin(source.fnPtBinLimitsD0forD0ptbin),
  fPtBinLimitsD0forD0ptbin(0x0),
  fIsUpperCutD0forD0ptbin(0x0),
  fIsCutUsedD0forD0ptbin(0x0),
  fVarNamesD0forD0ptbin(0x0),

  fMinITSNclsD0FirstDaughter(source.fMinITSNclsD0FirstDaughter),
  fMinTPCNclsD0FirstDaughter(source.fMinTPCNclsD0FirstDaughter),
  fUseITSRefitD0FirstDaughter(source.fUseITSRefitD0FirstDaughter),
  fUseTPCRefitD0FirstDaughter(source.fUseTPCRefitD0FirstDaughter),
  fUseFilterBitD0FirstDaughter(source.fUseFilterBitD0FirstDaughter),
  fFilterBitD0FirstDaughter(source.fFilterBitD0FirstDaughter),
  fMinPtD0FirstDaughter(source.fMinPtD0FirstDaughter),
  fMaxAbsEtaD0FirstDaughter(source.fMaxAbsEtaD0FirstDaughter),
  fHardSelectionArrayITSD0FirstDaughter(),
  fSoftSelectionArrayITSD0FirstDaughter(),
  fNSoftITSCutD0FirstDaughter(source.fNSoftITSCutD0FirstDaughter),

  fMinITSNclsD0SecondDaughter(source.fMinITSNclsD0SecondDaughter),
  fMinTPCNclsD0SecondDaughter(source.fMinTPCNclsD0SecondDaughter),
  fUseITSRefitD0SecondDaughter(source.fUseITSRefitD0SecondDaughter),
  fUseTPCRefitD0SecondDaughter(source.fUseTPCRefitD0SecondDaughter),
  fUseFilterBitD0SecondDaughter(source.fUseFilterBitD0SecondDaughter),
  fFilterBitD0SecondDaughter(source.fFilterBitD0SecondDaughter),
  fMinPtD0SecondDaughter(source.fMinPtD0SecondDaughter),
  fMaxAbsEtaD0SecondDaughter(source.fMaxAbsEtaD0SecondDaughter),
  fHardSelectionArrayITSD0SecondDaughter(),
  fSoftSelectionArrayITSD0SecondDaughter(),
  fNSoftITSCutD0SecondDaughter(source.fNSoftITSCutD0SecondDaughter),

  fMinITSNclsBPlusPion(source.fMinITSNclsBPlusPion),
  fMinTPCNclsBPlusPion(source.fMinTPCNclsBPlusPion),
  fUseITSRefitBPlusPion(source.fUseITSRefitBPlusPion),
  fUseTPCRefitBPlusPion(source.fUseTPCRefitBPlusPion),
  fUseFilterBitBPlusPion(source.fUseFilterBitBPlusPion),
  fFilterBitBPlusPion(source.fFilterBitBPlusPion),
  fMinPtBPlusPion(source.fMinPtBPlusPion),
  fMaxAbsEtaBPlusPion(source.fMaxAbsEtaBPlusPion),
  fHardSelectionArrayITSBPlusPion(),
  fSoftSelectionArrayITSBPlusPion(),
  fNSoftITSCutBPlusPion(source.fNSoftITSCutBPlusPion),

  fMind0D0FirstDaughter(source.fMind0D0FirstDaughter),
  fMind0D0SecondDaughter(source.fMind0D0SecondDaughter),
  fMind0BPlusPion(source.fMind0BPlusPion),
  fFiducialYCut(source.fFiducialYCut),

  fnVariablesForCutOptimization(source.fnVariablesForCutOptimization),
  fnCutsForOptimization(source.fnCutsForOptimization),
  fGlobalIndexCutOptimization(source.fGlobalIndexCutOptimization),
  fCutsRDForCutOptimization(0x0),
  fIsUpperCutForCutOptimization(0x0),
  fCutIndexForCutOptimization(0x0),
  fSigmaForCutOptimization(0x0),
  fNumberOfSigmaBinsForCutOptimization(source.fNumberOfSigmaBinsForCutOptimization)
{
  //
  // Copy constructor
  //

  if(source.GetTrackCutsSoftPi()) AddTrackCutsSoftPi(source.GetTrackCutsSoftPi());

  if (source.fPtBinLimitsD0forD0ptbin) SetPtBinsD0forD0ptbin(source.fnPtBinLimitsD0forD0ptbin, source.fPtBinLimitsD0forD0ptbin);
  if (source.fVarNamesD0forD0ptbin) SetVarNamesD0forD0ptbin(source.fnVarsD0forD0ptbin, source.fVarNamesD0forD0ptbin, source.fIsUpperCut);
  if (source.fIsCutUsed)
  {
    if (fIsCutUsed) {
      delete [] fIsCutUsed;
      fIsCutUsed = nullptr;
    }
    fIsCutUsed = new Bool_t[(source.GetNPtBins()) * (source.GetNVars())];

    for (Int_t i = 0; i < source.fnVars; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBins; j++)
      {
        Bool_t bUse = source.GetIsCutUsed(i, j);
        SetIsCutUsed(i, j, bUse);
      }
    }
  }
  if (source.fIsCutUsedD0forD0ptbin)
  {
    if (fIsCutUsedD0forD0ptbin) {
      delete [] fIsCutUsedD0forD0ptbin;
      fIsCutUsedD0forD0ptbin = nullptr;
    }
    fIsCutUsedD0forD0ptbin = new Bool_t[(source.GetNPtBinsD0forD0ptbin()) * (source.GetNVarsD0forD0ptbin())];
    for (Int_t i = 0; i < source.fnVarsD0forD0ptbin; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBinsD0forD0ptbin; j++)
      {
        Bool_t bUse = source.GetIsCutUsedD0forD0ptbin(i, j);
        SetIsCutUsedD0forD0ptbin(i, j, bUse);
      }
    }
  }
  if (source.fCutsRDForCutOptimization) InitializeCutsForCutOptimization(source.fnCutsForOptimization, source.fnVariablesForCutOptimization);
  if (source.fIsUpperCutForCutOptimization)
  {
    if (fIsUpperCutForCutOptimization) {
      delete [] fIsUpperCutForCutOptimization;
      fIsUpperCutForCutOptimization = nullptr;
    }
    fIsUpperCutForCutOptimization = new Bool_t[source.fnVariablesForCutOptimization];
    for (Int_t i = 0; i < source.fnVariablesForCutOptimization; i++)
    {
      Bool_t bUpperCut = source.GetIsUpperCutForCutOptimization(i);
      SetIsUpperCutForCutOptimization(i, bUpperCut);
    }
  }
  if (source.fCutIndexForCutOptimization)
  {
    if (fCutIndexForCutOptimization) {
      delete [] fCutIndexForCutOptimization;
      fCutIndexForCutOptimization = nullptr;
    }
    fCutIndexForCutOptimization = new Int_t[source.fnVariablesForCutOptimization];
    for (Int_t i = 0; i < source.fnVariablesForCutOptimization; i++)
    {
      Int_t nCutIndex = source.GetCutIndexForCutOptimization(i);
      SetCutIndexForCutOptimization(i, nCutIndex);
    }
  }
  if (source.fSigmaForCutOptimization)
  {
    if (fSigmaForCutOptimization) {
      delete [] fSigmaForCutOptimization;
      fSigmaForCutOptimization = nullptr;
    }
    fSigmaForCutOptimization = new Float_t[source.fnPtBins];
    for (Int_t i = 0; i < source.fnPtBins; i++)
    {
      Double_t binSigma = source.GetSigmaForCutOptimization(i);
      SetSigmaForCutOptimization(binSigma, i);
    }
  }
  if (source.fCutsRDD0forD0ptbin) SetCutsD0forD0ptbin(source.fGlobalIndexD0forD0ptbin, source.fCutsRDD0forD0ptbin);
  if (source.fCutsRDForCutOptimization) SetCutsForCutOptimization(source.fGlobalIndexCutOptimization, source.fCutsRDForCutOptimization);

  SetHardSelectionArrayITSD0FirstDaughter(source.fHardSelectionArrayITSD0FirstDaughter);
  SetSoftSelectionArrayITSD0FirstDaughter(source.fSoftSelectionArrayITSD0FirstDaughter);
  SetHardSelectionArrayITSD0SecondDaughter(source.fHardSelectionArrayITSD0SecondDaughter);
  SetSoftSelectionArrayITSD0SecondDaughter(source.fSoftSelectionArrayITSD0SecondDaughter);
  SetHardSelectionArrayITSBPlusPion(source.fHardSelectionArrayITSBPlusPion);
  SetSoftSelectionArrayITSBPlusPion(source.fSoftSelectionArrayITSBPlusPion);
}
//--------------------------------------------------------------------------
AliRDHFCutsBPlustoD0Pi::~AliRDHFCutsBPlustoD0Pi() {
  //
  // Default Destructor
  //
  if (fTrackCutsSoftPi) { delete fTrackCutsSoftPi; fTrackCutsSoftPi = nullptr;}
  if (fIsCutUsed) { delete [] fIsCutUsed; fIsCutUsed = nullptr; }
  if (fCutsRDD0forD0ptbin) { delete [] fCutsRDD0forD0ptbin; fCutsRDD0forD0ptbin = nullptr; }
  if (fPtBinLimitsD0forD0ptbin) { delete [] fPtBinLimitsD0forD0ptbin; fPtBinLimitsD0forD0ptbin = nullptr; }
  if (fIsUpperCutD0forD0ptbin) { delete [] fIsUpperCutD0forD0ptbin; fIsUpperCutD0forD0ptbin = nullptr; }
  if (fIsCutUsedD0forD0ptbin) { delete [] fIsCutUsedD0forD0ptbin; fIsCutUsedD0forD0ptbin = nullptr; }
  if (fVarNamesD0forD0ptbin) { delete [] fVarNamesD0forD0ptbin; fVarNamesD0forD0ptbin = nullptr; }
  if (fCutsRDForCutOptimization) { delete [] fCutsRDForCutOptimization; fCutsRDForCutOptimization = nullptr; }
  if (fIsUpperCutForCutOptimization) { delete [] fIsUpperCutForCutOptimization; fIsUpperCutForCutOptimization = nullptr; }
  if (fCutIndexForCutOptimization) { delete [] fCutIndexForCutOptimization; fCutIndexForCutOptimization = nullptr; }
  if (fSigmaForCutOptimization) { delete [] fSigmaForCutOptimization; fSigmaForCutOptimization = nullptr; }
}
//--------------------------------------------------------------------------
AliRDHFCutsBPlustoD0Pi &AliRDHFCutsBPlustoD0Pi::operator=(const AliRDHFCutsBPlustoD0Pi &source)
{
  //
  // assignment operator
  //

  if (&source == this) return *this;

  AliRDHFCuts::operator=(source);

  if(source.GetTrackCutsSoftPi()) {
    delete fTrackCutsSoftPi;
    fTrackCutsSoftPi = new AliESDtrackCuts(*(source.GetTrackCutsSoftPi()));
  }
  fMaxPtPid = source.fMaxPtPid;
  fTPCflag = source.fTPCflag;
  fCircRadius = source.fCircRadius;
  fnVarsD0forD0ptbin = source.fnVarsD0forD0ptbin;
  fnPtBinsD0forD0ptbin = source.fnPtBinsD0forD0ptbin;
  fGlobalIndexD0forD0ptbin = source.fGlobalIndexD0forD0ptbin;
  fnPtBinLimitsD0forD0ptbin = source.fnPtBinLimitsD0forD0ptbin;
  fMinITSNclsD0FirstDaughter = source.fMinITSNclsD0FirstDaughter;
  fMinTPCNclsD0FirstDaughter = source.fMinTPCNclsD0FirstDaughter;
  fUseITSRefitD0FirstDaughter = source.fUseITSRefitD0FirstDaughter;
  fUseTPCRefitD0FirstDaughter = source.fUseTPCRefitD0FirstDaughter;
  fUseFilterBitD0FirstDaughter = source.fUseFilterBitD0FirstDaughter;
  fFilterBitD0FirstDaughter = source.fFilterBitD0FirstDaughter;
  fMinPtD0FirstDaughter = source.fMinPtD0FirstDaughter;
  fMinITSNclsD0SecondDaughter = source.fMinITSNclsD0SecondDaughter;
  fMinTPCNclsD0SecondDaughter = source.fMinTPCNclsD0SecondDaughter;
  fUseITSRefitD0SecondDaughter = source.fUseITSRefitD0SecondDaughter;
  fUseTPCRefitD0SecondDaughter = source.fUseTPCRefitD0SecondDaughter;
  fUseFilterBitD0SecondDaughter = source.fUseFilterBitD0SecondDaughter;
  fFilterBitD0SecondDaughter = source.fFilterBitD0SecondDaughter;
  fMinPtD0SecondDaughter = source.fMinPtD0SecondDaughter;
  fMinITSNclsBPlusPion = source.fMinITSNclsBPlusPion;
  fMinTPCNclsBPlusPion = source.fMinTPCNclsBPlusPion;
  fUseITSRefitBPlusPion = source.fUseITSRefitBPlusPion;
  fUseTPCRefitBPlusPion = source.fUseTPCRefitBPlusPion;
  fUseFilterBitBPlusPion = source.fUseFilterBitBPlusPion;
  fFilterBitBPlusPion = source.fFilterBitBPlusPion;
  fMinPtBPlusPion = source.fMinPtBPlusPion;
  fMaxAbsEtaD0FirstDaughter = source.fMaxAbsEtaD0FirstDaughter;
  fNSoftITSCutD0FirstDaughter = source.fNSoftITSCutD0FirstDaughter;
  fMaxAbsEtaD0SecondDaughter = source.fMaxAbsEtaD0SecondDaughter;
  fNSoftITSCutD0SecondDaughter = source.fNSoftITSCutD0SecondDaughter;
  fMaxAbsEtaBPlusPion = source.fMaxAbsEtaBPlusPion;
  fNSoftITSCutBPlusPion = source.fNSoftITSCutBPlusPion;
  fMind0D0FirstDaughter = source.fMind0D0FirstDaughter;
  fMind0D0SecondDaughter = source.fMind0D0SecondDaughter;
  fMind0BPlusPion = source.fMind0BPlusPion;
  fFiducialYCut = source.fFiducialYCut;
  fnVariablesForCutOptimization = source.fnVariablesForCutOptimization;
  fnCutsForOptimization = source.fnCutsForOptimization;
  fGlobalIndexCutOptimization = source.fGlobalIndexCutOptimization;
  fNumberOfSigmaBinsForCutOptimization = source.fNumberOfSigmaBinsForCutOptimization;

  if (source.fPtBinLimitsD0forD0ptbin) SetPtBinsD0forD0ptbin(source.fnPtBinLimitsD0forD0ptbin, source.fPtBinLimitsD0forD0ptbin);
  if (source.fVarNamesD0forD0ptbin) SetVarNamesD0forD0ptbin(source.fnVarsD0forD0ptbin, source.fVarNamesD0forD0ptbin, source.fIsUpperCut);
  if (source.fIsCutUsed)
  {
    if (fIsCutUsed) {
      delete [] fIsCutUsed;
      fIsCutUsed = nullptr;
    }
    fIsCutUsed = new Bool_t[(source.GetNPtBins()) * (source.GetNVars())];

    for (Int_t i = 0; i < source.fnVars; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBins; j++)
      {
        Bool_t bUse = source.GetIsCutUsed(i, j);
        SetIsCutUsed(i, j, bUse);
      }
    }
  }
  if (source.fIsCutUsedD0forD0ptbin)
  {
    if (fIsCutUsedD0forD0ptbin) {
      delete [] fIsCutUsedD0forD0ptbin;
      fIsCutUsedD0forD0ptbin = nullptr;
    }
    fIsCutUsedD0forD0ptbin = new Bool_t[(source.GetNPtBinsD0forD0ptbin()) * (source.GetNVarsD0forD0ptbin())];
    for (Int_t i = 0; i < source.fnVarsD0forD0ptbin; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBinsD0forD0ptbin; j++)
      {
        Bool_t bUse = source.GetIsCutUsedD0forD0ptbin(i, j);
        SetIsCutUsedD0forD0ptbin(i, j, bUse);
      }
    }
  }
  if (source.fCutsRDForCutOptimization) InitializeCutsForCutOptimization(source.fnCutsForOptimization, source.fnVariablesForCutOptimization);
  if (source.fIsUpperCutForCutOptimization)
  {
    if (fIsUpperCutForCutOptimization) {
      delete [] fIsUpperCutForCutOptimization;
      fIsUpperCutForCutOptimization = nullptr;
    }
    fIsUpperCutForCutOptimization = new Bool_t[source.fnVariablesForCutOptimization];
    for (Int_t i = 0; i < source.fnVariablesForCutOptimization; i++)
    {
      Bool_t bUpperCut = source.GetIsUpperCutForCutOptimization(i);
      SetIsUpperCutForCutOptimization(i, bUpperCut);
    }
  }
  if (source.fCutIndexForCutOptimization)
  {
    if (fCutIndexForCutOptimization) {
      delete [] fCutIndexForCutOptimization;
      fCutIndexForCutOptimization = nullptr;
    }
    fCutIndexForCutOptimization = new Int_t[source.fnVariablesForCutOptimization];
    for (Int_t i = 0; i < source.fnVariablesForCutOptimization; i++)
    {
      Int_t nCutIndex = source.GetCutIndexForCutOptimization(i);
      SetCutIndexForCutOptimization(i, nCutIndex);
    }
  }
  if (source.fSigmaForCutOptimization)
  {
    if (fSigmaForCutOptimization) {
      delete [] fSigmaForCutOptimization;
      fSigmaForCutOptimization = nullptr;
    }
    fSigmaForCutOptimization = new Float_t[source.fnPtBins];
    for (Int_t i = 0; i < source.fnPtBins; i++)
    {
      Double_t binSigma = source.GetSigmaForCutOptimization(i);
      SetSigmaForCutOptimization(binSigma, i);
    }
  }
  if (source.fCutsRDD0forD0ptbin) SetCutsD0forD0ptbin(source.fGlobalIndexD0forD0ptbin, source.fCutsRDD0forD0ptbin);
  if (source.fCutsRDForCutOptimization) SetCutsForCutOptimization(source.fGlobalIndexCutOptimization, source.fCutsRDForCutOptimization);

  SetHardSelectionArrayITSD0FirstDaughter(source.fHardSelectionArrayITSD0FirstDaughter);
  SetSoftSelectionArrayITSD0FirstDaughter(source.fSoftSelectionArrayITSD0FirstDaughter);
  SetHardSelectionArrayITSD0SecondDaughter(source.fHardSelectionArrayITSD0SecondDaughter);
  SetSoftSelectionArrayITSD0SecondDaughter(source.fSoftSelectionArrayITSD0SecondDaughter);
  SetHardSelectionArrayITSBPlusPion(source.fHardSelectionArrayITSBPlusPion);
  SetSoftSelectionArrayITSBPlusPion(source.fSoftSelectionArrayITSBPlusPion);

  return *this;
}
//--------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::GetCutVarsForOpt(AliAODRecoDecayHF* /*d*/, Float_t* /*vars*/, Int_t /*nvars*/, Int_t* /*pdgdaughters*/) {
  // not yet used

  return;
}
//--------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[75]) {
  //
  // In this function we apply the selection cuts on the BPlus candidate and its daughters.
  // The function returns 0 if the candidate is cut and is able to return information on which cuts the candidate passes.
  //

  fIsSelectedCuts = 0;
  fIsSelectedPID = 0;

  // The cuts used in this class have to be set in the maketfile.
  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateBPlus = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateBPlus) {
    cout << "candidateBPlus null" << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateBPlus->GetDaughter(1);
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateBPlus->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 1" << endl;
    return 0;
  }

  AliAODVertex * vertexBPlus = candidateBPlus->GetSecondaryVtx();
  if (!vertexBPlus) {
    cout << "vertexBPlus null" << endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  Int_t returnvalue = 1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();

  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    // We check to which pt bin the candidate belongs
    Int_t ptbin = PtBin(candidateBPlus->Pt());
    if (ptbin == -1) return -1;

    // We obtain the variable values in the section below
    // D0Mass and BPlusmass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mBPlusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mBPlusPDG - mD0PDG;

    // Half width BPlus mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 421;
    Double_t invMassBPlus = candidateBPlus->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mBPlusPDG - invMassBPlus);
    Double_t invMassDelta = TMath::Abs(deltaPDG - (DeltaInvMassBPlusKpipi(candidateBPlus)));

    Double_t pointingAngle = candidateBPlus->CosPointingAngle();
    Double_t dcaMother = candidateBPlus->GetDCA();
    Double_t ptMother = candidateBPlus->Pt();
    Double_t momentumMother = candidateBPlus->P();
    Double_t ptD0 = candidateD0->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateBPlus);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateBPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateBPlus->Getd0Prong(1));

    Double_t impactProduct = candidateBPlus->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateBPlus->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateBPlus->Px() * candidatePion->Px() + candidateBPlus->Py() * candidatePion->Py() + candidateBPlus->Pz() * candidatePion->Pz()) / (candidateBPlus->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateBPlus->Px() * candidateD0->Px() + candidateBPlus->Py() * candidateD0->Py() + candidateBPlus->Pz() * candidateD0->Pz()) / (candidateBPlus->P() * candidateD0->P());

    Double_t cosThetaStar = candidateBPlus->CosThetaStar(0, 521, 211, 421);
    Double_t vertexDistance = vertexBPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateBPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    Double_t pseudoProperDecayLength = ((vertexBPlus->GetX() - primaryVertex->GetX()) * candidateBPlus->Px() / TMath::Abs(candidateBPlus->Pt())) + ((vertexBPlus->GetY() - primaryVertex->GetY()) * candidateBPlus->Py() / TMath::Abs(candidateBPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateBPlus->Phi();
    Double_t theta = candidateBPlus->Theta();
    Double_t covMatrix[21];
    candidateBPlus->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateBPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateBPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexBPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateBPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexBPlus->GetChi2perNDF();

    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateBPlus->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}

    Double_t Normalizedd0D0firstdaughter = TMath::Abs(candidateD0->Getd0Prong(0)/candidateD0->Getd0errProng(0));
    Double_t Normalizedd0D0seconddaughter = TMath::Abs(candidateD0->Getd0Prong(1)/candidateD0->Getd0errProng(1));
    Double_t Normalizedd0D0 = TMath::Abs(candidateBPlus->Getd0Prong(1)/candidateBPlus->Getd0errProng(1));
    Double_t Normalizedd0BPluspion = TMath::Abs(candidateBPlus->Getd0Prong(0)/candidateBPlus->Getd0errProng(0));
    Double_t Normalizedd0BPlus = TMath::Abs(d0[0]/covd0z0[0]);
    Double_t NormalizedimpactproductD0 = (candidateD0->Getd0Prong(0)/candidateD0->Getd0errProng(0)) * (candidateD0->Getd0Prong(1)/candidateD0->Getd0errProng(1));
    Double_t NormalizedimpactproductBPlus = (candidateBPlus->Getd0Prong(0)/candidateBPlus->Getd0errProng(0)) * (candidateBPlus->Getd0Prong(1)/candidateBPlus->Getd0errProng(1));

    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 39;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 40;
    cutVariableValue = invMassDelta;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 41;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 42;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt BPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 43;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" -------------------------------------------------
    nCutIndex = 44;
    cutVariableValue = ptD0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 45;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 BPlus [cm]" -------------------------------------------------------
    nCutIndex = 46;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 D0 [cm]"-----------------------------------------------------
    nCutIndex = 47;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 48;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 49;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 50;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 51;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 52;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 53;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 54;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 55;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 56;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 57;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 58;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 59;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 60;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 61;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 62;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 63;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 64;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 65;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 66;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 67;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0D0firstdaughter" ----------------------------------------------------
    nCutIndex = 68;
    cutVariableValue = Normalizedd0D0firstdaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0D0seconddaughter" ----------------------------------------------------
    nCutIndex = 69;
    cutVariableValue = Normalizedd0D0seconddaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0D0" ----------------------------------------------------
    nCutIndex = 70;
    cutVariableValue = Normalizedd0D0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0BPluspion" ----------------------------------------------------
    nCutIndex = 71;
    cutVariableValue = Normalizedd0BPluspion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0BPlus" ----------------------------------------------------
    nCutIndex = 72;
    cutVariableValue = Normalizedd0BPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "NormalizedimpactproductD0" ----------------------------------------------------
    nCutIndex = 73;
    cutVariableValue = NormalizedimpactproductD0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "NormalizedimpactproductBPlus" ----------------------------------------------------
    nCutIndex = 74;
    cutVariableValue = NormalizedimpactproductBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------


    // select D0
    bPassedCut = IsD0FromBPlusSelected(ptMother, candidateBPlus, selectionLevel, aod, bCutArray);
  }

  if (bPassedCut == kFALSE)
  {
    returnvalue = 0;
  } else
  {
    for (Int_t i = 39; i < 75; ++i)
    {
      if (bCutArray[i] == kTRUE) {
        returnvalue = 0;
        break;
      }
    }
  }

  fIsSelectedCuts = returnvalue;

  return returnvalue;
}
//_________________________________________________________________________________________________
Int_t AliRDHFCutsBPlustoD0Pi::IsD0FromBPlusSelected(Double_t ptBPlus, TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[75]) {
  //
  // Apply selection on D0 candidate from BPlus candidate. We have to pass the BPlus candidate to this function to get variables w.r.t. BPlus vertex.
  //

  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateBPlus = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateBPlus) {
    cout << "candidateBPlus null" << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateBPlus->GetDaughter(1);
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 2" << endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if (!candidateKaon) {
    cout << "candidateKaon null" << endl;
    return 0;
  }

  AliAODVertex * vertexBPlus = candidateBPlus->GetSecondaryVtx();
  if (!vertexBPlus) {
    cout << "vertexBPlus null" << endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if (!vertexD0) {
    cout << "vertexD0 null" << endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  Int_t returnvalue = 1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();


  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBin(ptBPlus);

    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // D0 window - invariant mass
    Int_t chargeBPlus = candidateBPlus->Charge();
    UInt_t prongs[2];
    if (chargeBPlus == -1)
    {
      prongs[0] = 211;
      prongs[1] = 321;
    }
    else if (chargeBPlus == 1)
    {
      prongs[0] = 321;
      prongs[1] = 211;
    }
    else
    {
      cout << "Wrong charge BPlus." << endl;
      return 0;
    }
    Double_t invMassD0 = candidateD0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));

    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) / (candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) / (candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0, 421, 211, 321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateD0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexD0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateD0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexD0->GetChi2perNDF();


    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 25;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------



    AliAODRecoDecay* candidateD0toBPlus = (AliAODRecoDecay*)candidateD0;
    AliExternalTrackParam pionD0Track;
    AliExternalTrackParam kaonD0Track;

    Double_t d0z0DSVert[2], covd0z0DSVert[3], d0DSVert[2];

    pionD0Track.CopyFromVTrack(candidatePion);
    pionD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    kaonD0Track.CopyFromVTrack(candidateKaon);
    kaonD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(candidateD0);
    Double_t d0z0D0DSVert[2], covd0z0D0DSVert[3], d0D0DSVert;
    motherTrack.PropagateToDCA(vertexBPlus, bz, 100., d0z0D0DSVert, covd0z0D0DSVert);
    d0D0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToBPlus = d0DSVert[0] * d0DSVert[1];
    Double_t impactProductXYToBPlus = candidateD0toBPlus->ImpParXY(vertexBPlus);

    Double_t pointingAngleToBPlus = candidateD0toBPlus->CosPointingAngle(vertexBPlus);
    Double_t d0FirstDaughterToBPlus = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToBPlus = TMath::Abs(d0DSVert[1]);
    Double_t normDecayLengthToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus);

    Double_t pseudoProperDecayLengthDSVert = ((vertexD0->GetX() - vertexBPlus->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - vertexBPlus->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTimeToBPlus = pseudoProperDecayLengthDSVert * pdgMassMother / ptMother;
    Double_t DecayTimeToBPlus = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phiDSVert = candidateD0->Phi();
    Double_t thetaDSVert = candidateD0->Theta();
    Double_t covMatrixDSVert[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                    + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                    + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                    + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                    + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                    + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTimeToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    // "pointingAngleToBPlus" ---------------------------------------------
    nCutIndex = 29;
    cutVariableValue = pointingAngleToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0MotherToBPlus" --------------------------------------------------
    nCutIndex = 30;
    cutVariableValue = d0D0DSVert;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0FirstDaughterToBPlus" -------------------------------------------
    nCutIndex = 31;
    cutVariableValue = d0FirstDaughterToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0SecondDaughterToBPlus" ------------------------------------------
    nCutIndex = 32;
    cutVariableValue = d0SecondDaughterToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "impactProductToBPlus" ---------------------------------------------
    nCutIndex = 33;
    cutVariableValue = impactProductToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "impactProductXYToBPlus" -------------------------------------------
    nCutIndex = 34;
    cutVariableValue = impactProductXYToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLengthToBPlus" -------------------------------------------
    nCutIndex = 35;
    cutVariableValue = normDecayLengthToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToBPlus" -------------------------------------
    nCutIndex = 36;
    cutVariableValue = pseudoProperDecayTimeToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTimeToBPlus" -------------------------------------------------
    nCutIndex = 37;
    cutVariableValue = DecayTimeToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToBPlus" ---------------------------------------------
    nCutIndex = 38;
    cutVariableValue = normalizedDecayTimeToBPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------
  }

  if (bPassedCut == kFALSE)
  {
    returnvalue = 0;
  } else
  for (Int_t i = 0; i < 39; ++i)
  {
    if (bCutArray[i] == kTRUE) {
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsD0forD0ptbinSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[29]) {
//
  // Apply selection on D0 candidate.
  //

  if (!fCutsRDD0forD0ptbin) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 3" << endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if (!candidateKaon) {
    cout << "candidateKaon null" << endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if (!vertexD0) {
    cout << "vertexD0 null" << endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  Int_t returnvalue = 1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();


  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBinD0forD0ptbin(candidateD0->Pt());
    if (ptbin == -1) return -1;

    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // D0 window - invariant mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 321;
    Double_t invMassD0 = candidateD0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    UInt_t prongs2[2];
    prongs2[1] = 211; prongs2[0] = 321;
    Double_t invMassD02 = candidateD0->InvMass(2, prongs2);
    Double_t invMassDifference2 = TMath::Abs(mD0PDG - invMassD02);

    if (invMassDifference > invMassDifference2) invMassDifference = invMassDifference2;

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));

    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) / (candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) / (candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0, 421, 211, 321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateD0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexD0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateD0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexD0->GetChi2perNDF();

    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 25;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

  }

  if (bPassedCut == kFALSE)
  {
    returnvalue = 0;
  } else
  for (Int_t i = 0; i < 29; ++i)
  {
    if (bCutArray[i] == kTRUE) {
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//--------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod) {
  //
  // In this function we apply the selection cuts on the BPlus candidate and its daughters.
  // The function returns 0 if the candidate is cut.
  //

  // The cuts used in this class have to be set in the maketfile.
  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateBPlus = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateBPlus) {
    cout << "candidateBPlus null" << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateBPlus->GetDaughter(1);
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateBPlus->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 1" << endl;
    return 0;
  }

  AliAODVertex * vertexBPlus = candidateBPlus->GetSecondaryVtx();
  if (!vertexBPlus) {
    cout << "vertexBPlus null" << endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();

  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    // We check to which pt bin the candidate belongs
    Int_t ptbin = PtBin(candidateBPlus->Pt());
    if (ptbin == -1) return -1;

    // We obtain the variable values in the section below
    // D0Mass and BPlusmass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mBPlusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mBPlusPDG - mD0PDG;

    // Half width BPlus mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 421;
    Double_t invMassBPlus = candidateBPlus->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mBPlusPDG - invMassBPlus);
    Double_t invMassDelta = TMath::Abs(deltaPDG - (DeltaInvMassBPlusKpipi(candidateBPlus)));

    Double_t pointingAngle = candidateBPlus->CosPointingAngle();
    Double_t dcaMother = candidateBPlus->GetDCA();
    Double_t ptMother = candidateBPlus->Pt();
    Double_t momentumMother = candidateBPlus->P();
    Double_t ptD0 = candidateD0->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateBPlus);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateBPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateBPlus->Getd0Prong(1));

    Double_t impactProduct = candidateBPlus->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateBPlus->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateBPlus->Px() * candidatePion->Px() + candidateBPlus->Py() * candidatePion->Py() + candidateBPlus->Pz() * candidatePion->Pz()) / (candidateBPlus->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateBPlus->Px() * candidateD0->Px() + candidateBPlus->Py() * candidateD0->Py() + candidateBPlus->Pz() * candidateD0->Pz()) / (candidateBPlus->P() * candidateD0->P());

    Double_t cosThetaStar = candidateBPlus->CosThetaStar(0, 521, 211, 421);
    Double_t vertexDistance = vertexBPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateBPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    Double_t pseudoProperDecayLength = ((vertexBPlus->GetX() - primaryVertex->GetX()) * candidateBPlus->Px() / TMath::Abs(candidateBPlus->Pt())) + ((vertexBPlus->GetY() - primaryVertex->GetY()) * candidateBPlus->Py() / TMath::Abs(candidateBPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateBPlus->Phi();
    Double_t theta = candidateBPlus->Theta();
    Double_t covMatrix[21];
    candidateBPlus->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateBPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateBPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexBPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateBPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexBPlus->GetChi2perNDF();

    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateBPlus->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 39;
    cutVariableValue = invMassDifference;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 40;
    cutVariableValue = invMassDelta;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 41;
    cutVariableValue = pointingAngle;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 42;
    cutVariableValue = dcaMother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt BPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 43;
    cutVariableValue = ptMother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" -------------------------------------------------
    nCutIndex = 44;
    cutVariableValue = ptD0;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 45;
    cutVariableValue = ptPion;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 BPlus [cm]" -------------------------------------------------------
    nCutIndex = 46;
    cutVariableValue = d0Mother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]"-----------------------------------------------------
    nCutIndex = 47;
    cutVariableValue = d0firstTrack;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 48;
    cutVariableValue = d0secondTrack;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 49;
    cutVariableValue = impactProduct;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 50;
    cutVariableValue = impactProductXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 51;
    cutVariableValue = angleBetweenBothDaughters;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 52;
    cutVariableValue = angleMotherFirstDaughter;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 53;
    cutVariableValue = angleMotherSecondDaughter;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 54;
    cutVariableValue = cosThetaStar;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 55;
    cutVariableValue = vertexDistance;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 56;
    cutVariableValue = pseudoProperDecayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 57;
    cutVariableValue = decayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 58;
    cutVariableValue = normalizedDecayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 59;
    cutVariableValue = normDecayLength;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 60;
    cutVariableValue = dd0pr1;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 61;
    cutVariableValue = dd0pr2;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 62;
    cutVariableValue = dd0max;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 63;
    cutVariableValue = dd0min;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 64;
    cutVariableValue = cosPointingAngleXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 65;
    cutVariableValue = distanceXYToVertex;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 66;
    cutVariableValue = normalizedDecayLengthXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 67;
    cutVariableValue = chi2Vertex;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // select D0
    if (!IsD0FromBPlusSelectedMVA(ptMother, candidateBPlus, selectionLevel, aod, primaryVertex, bz)) return 0;
  }

  return 1;
}

Bool_t AliRDHFCutsBPlustoD0Pi::IsThisDaughterSelected(AliAODTrack *track, AliAODVertex *primary, const AliAODEvent* aod) {
  //
  // Daughter track selection
  //

  if(!fTrackCuts) return kTRUE;
  if(!track) {return kFALSE;}
  if(track->Charge()==0) return kFALSE; // it's not a track, but a V0

  Double_t pos[3],cov[6];
  primary->GetXYZ(pos);
  primary->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  Bool_t retval=kTRUE;

  // SetUseTPCtrackCutsOnThisDaughter(kTRUE);
  if(!IsDaughterSelected(track,&vESD,fTrackCuts,aod)) retval = kFALSE;

  return retval;
}

//_________________________________________________________________________________________________
Int_t AliRDHFCutsBPlustoD0Pi::IsD0FromBPlusSelectedMVA(Double_t ptBPlus, TObject* obj, Int_t selectionLevel, AliAODEvent* /*aod*/, AliAODVertex *primaryVertex, Double_t bz) {
  //
  // Apply selection on D0 candidate from BPlus candidate. We have to pass the BPlus candidate to this function to get variables w.r.t. BPlus vertex.
  //

  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateBPlus = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateBPlus) {
    cout << "candidateBPlus null" << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateBPlus->GetDaughter(1);
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 2" << endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if (!candidateKaon) {
    cout << "candidateKaon null" << endl;
    return 0;
  }

  AliAODVertex * vertexBPlus = candidateBPlus->GetSecondaryVtx();
  if (!vertexBPlus) {
    cout << "vertexBPlus null" << endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if (!vertexD0) {
    cout << "vertexD0 null" << endl;
    return 0;
  }

  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBin(ptBPlus);

    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // D0 window - invariant mass
    Int_t chargeBPlus = candidateBPlus->Charge();
    UInt_t prongs[2];
    if (chargeBPlus == -1)
    {
      prongs[0] = 211;
      prongs[1] = 321;
    }
    else if (chargeBPlus == 1)
    {
      prongs[0] = 321;
      prongs[1] = 211;
    }
    else
    {
      cout << "Wrong charge BPlus." << endl;
      return 0;
    }
    Double_t invMassD0 = candidateD0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));

    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) / (candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) / (candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0, 421, 211, 321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateD0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexD0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateD0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexD0->GetChi2perNDF();


    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 25;
    cutVariableValue = cosPointingAngleXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = distanceXYToVertex;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayLengthXY;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = chi2Vertex;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------



    AliAODRecoDecay* candidateD0toBPlus = (AliAODRecoDecay*)candidateD0;
    AliExternalTrackParam pionD0Track;
    AliExternalTrackParam kaonD0Track;

    Double_t d0z0DSVert[2], covd0z0DSVert[3], d0DSVert[2];

    pionD0Track.CopyFromVTrack(candidatePion);
    pionD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    kaonD0Track.CopyFromVTrack(candidateKaon);
    kaonD0Track.PropagateToDCA(vertexBPlus, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(candidateD0);
    Double_t d0z0D0DSVert[2], covd0z0D0DSVert[3], d0D0DSVert;
    motherTrack.PropagateToDCA(vertexBPlus, bz, 100., d0z0D0DSVert, covd0z0D0DSVert);
    d0D0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToBPlus = d0DSVert[0] * d0DSVert[1];
    Double_t impactProductXYToBPlus = candidateD0toBPlus->ImpParXY(vertexBPlus);

    Double_t pointingAngleToBPlus = candidateD0toBPlus->CosPointingAngle(vertexBPlus);
    Double_t d0FirstDaughterToBPlus = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToBPlus = TMath::Abs(d0DSVert[1]);
    Double_t normDecayLengthToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus);

    Double_t pseudoProperDecayLengthDSVert = ((vertexD0->GetX() - vertexBPlus->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - vertexBPlus->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTimeToBPlus = pseudoProperDecayLengthDSVert * pdgMassMother / ptMother;
    Double_t DecayTimeToBPlus = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phiDSVert = candidateD0->Phi();
    Double_t thetaDSVert = candidateD0->Theta();
    Double_t covMatrixDSVert[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                    + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                    + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                    + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                    + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                    + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTimeToBPlus = candidateD0toBPlus->NormalizedDecayLength(vertexBPlus) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    // "pointingAngleToBPlus" ---------------------------------------------
    nCutIndex = 29;
    cutVariableValue = pointingAngleToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0MotherToBPlus" --------------------------------------------------
    nCutIndex = 30;
    cutVariableValue = d0D0DSVert;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0FirstDaughterToBPlus" -------------------------------------------
    nCutIndex = 31;
    cutVariableValue = d0FirstDaughterToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0SecondDaughterToBPlus" ------------------------------------------
    nCutIndex = 32;
    cutVariableValue = d0SecondDaughterToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "impactProductToBPlus" ---------------------------------------------
    nCutIndex = 33;
    cutVariableValue = impactProductToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "impactProductXYToBPlus" -------------------------------------------
    nCutIndex = 34;
    cutVariableValue = impactProductXYToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normDecayLengthToBPlus" -------------------------------------------
    nCutIndex = 35;
    cutVariableValue = normDecayLengthToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToBPlus" -------------------------------------
    nCutIndex = 36;
    cutVariableValue = pseudoProperDecayTimeToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "DecayTimeToBPlus" -------------------------------------------------
    nCutIndex = 37;
    cutVariableValue = DecayTimeToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToBPlus" ---------------------------------------------
    nCutIndex = 38;
    cutVariableValue = normalizedDecayTimeToBPlus;
    if(!ApplyCutOnVariableMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------
  }

  return 1;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsD0forD0ptbinSelectedMVA(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, AliAODVertex *primaryVertex, Double_t bz) {
  //
  // Apply selection on D0 candidate (for MVA, need to be merged with function for standard analysis)
  //

  if (!fCutsRDD0forD0ptbin) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateD0) {
    cout << "candidateD0 null" << endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if (!candidatePion) {
    cout << "candidatePion null 3" << endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if (!candidateKaon) {
    cout << "candidateKaon null" << endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if (!vertexD0) {
    cout << "vertexD0 null" << endl;
    return 0;
  }

  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return 0;
  }

  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBinD0forD0ptbin(candidateD0->Pt());
    if (ptbin == -1) return -1;

    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // D0 window - invariant mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 321;
    Double_t invMassD0 = candidateD0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    UInt_t prongs2[2];
    prongs2[1] = 211; prongs2[0] = 321;
    Double_t invMassD02 = candidateD0->InvMass(2, prongs2);
    Double_t invMassDifference2 = TMath::Abs(mD0PDG - invMassD02);

    if (invMassDifference > invMassDifference2) invMassDifference = invMassDifference2;

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));

    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) / (candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) / (candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) / (candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0, 421, 211, 321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9] * cp * cp * ct * ct // GetCovPxPx
                             + covMatrix[13] * 2.*cp * sp * ct * ct // GetCovPxPy
                             + covMatrix[18] * 2.*cp * ct * st // GetCovPxPz
                             + covMatrix[14] * sp * sp * ct * ct // GetCovPyPy
                             + covMatrix[19] * 2.*sp * ct * st // GetCovPyPz
                             + covMatrix[20] * st * st; // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateD0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexD0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateD0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexD0->GetChi2perNDF();

    //Topomatic
    Double_t dd0pr1 = 0.;
    Double_t dd0pr2 = 0.;
    Double_t dd0max = 0.;
    Double_t dd0min = 0.;
    for (Int_t ipr = 0; ipr < 2; ipr++)
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr, bz, diffIP, errdiffIP);
      Double_t normdd0 = 0.;
      if (errdiffIP > 0.) normdd0 = diffIP / errdiffIP;
      if (ipr == 0) dd0pr1 = normdd0;
      if (ipr == 1) dd0pr2 = normdd0;
    }
    if (TMath::Abs(dd0pr1) > TMath::Abs(dd0pr2)) {dd0max = dd0pr1; dd0min = dd0pr2;}
    else {dd0max = dd0pr2; dd0min = dd0pr1;}


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 25;
    cutVariableValue = cosPointingAngleXY;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = distanceXYToVertex;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayLengthXY;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = chi2Vertex;
    if(!ApplyCutOnVariableD0forD0ptbinMVA(nCutIndex, ptbin, cutVariableValue)) return 0;
    //---------------------------------------------------------------------
  }
    
  if(selectionLevel == AliRDHFCuts::kAll || selectionLevel == AliRDHFCuts::kTracks){

    if(!AreDaughtersSelected(candidateD0,aod)) return 0;

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(candidatePion);
    Double_t d0[2], covd0[3];
    particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

    if (TMath::Abs(d0[0]) < GetMind0D0FirstDaughter()) return 0;

    AliExternalTrackParam particleTrack2;
    particleTrack2.CopyFromVTrack(candidateKaon);
    Double_t d02[2], covd02[3];
    particleTrack2.PropagateToDCA(primaryVertex, bz, 100., d02, covd02);

    if (TMath::Abs(d0[0]) < GetMind0D0SecondDaughter()) return 0;
      
    if (UseFilterBitD0FirstDaughter() == kTRUE) {
      if (!(candidatePion->TestFilterMask(BIT(GetFilterBitD0FirstDaughter())))) return 0;
      if (!(candidateKaon->TestFilterMask(BIT(GetFilterBitD0SecondDaughter())))) return 0;
    }
  }

  return 1;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsBplusPionSelectedMVA(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, AliAODVertex *primaryVertex, Double_t bz) {
  //
  // Apply selection on D0 candidate.
  //

  AliAODTrack* candidatePion = (AliAODTrack*)obj;
  if (!candidatePion) {
    cout << "candidatePion null" << endl;
    return 0;
  }

  if(selectionLevel == AliRDHFCuts::kAll || selectionLevel == AliRDHFCuts::kTracks){

    //quick track quality pre selection to speed things up.
    //The first 4 cuts are also in AliESDtrackCuts object
    if (candidatePion->GetITSNcls() < 1) return 0;
    if (candidatePion->GetTPCNcls() < 1) return 0;
    if (candidatePion->GetStatus()&AliESDtrack::kITSpureSA) return 0;
    if (!(candidatePion->GetStatus()&AliESDtrack::kITSin)) return 0;
    if (candidatePion->GetID() < 0) return 0;
    Double_t covtest[21];
    if (!candidatePion->GetCovarianceXYZPxPyPz(covtest)) return 0;

    Double_t pos[3],cov[6];
    primaryVertex->GetXYZ(pos);
    primaryVertex->GetCovarianceMatrix(cov);
    const AliESDVertex vESD(pos,cov,100.,100);

    if(!IsDaughterSelected(candidatePion,&vESD,fTrackCutsSoftPi,aod)) return 0;

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(candidatePion);
    Double_t d0[2], covd0[3];
    particleTrack.PropagateToDCA(primaryVertex, bz, 100., d0, covd0);

    if (TMath::Abs(d0[0]) < GetMind0BPlusPion()) return 0;

    if (UseFilterBitBPlusPion() == kTRUE) {
      if (!(candidatePion->TestFilterMask(BIT(GetFilterBitBPlusPion())))) return 0;
    }

  }

  if(selectionLevel == AliRDHFCuts::kAll || selectionLevel == AliRDHFCuts::kPID){
    if (!(SelectPID(candidatePion, 2))) return 0;
  }

  return 1;

}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::IsD0SelectedPreRecVtxMVA(AliAODRecoDecayHF2Prong* d, AliAODTrack* pion, AliAODVertex *primaryVertex, Double_t bz, Int_t selLevel = 0)
{
  //
  // Preselection before reconstruction Bplus vertex to save time.
  //
  // selLevel =
  //  0: Everything
  //  1: PID of D0 daughters (using charge pion track. Always true when PID is off)
  //  2: Invariant mass cut D0 (Update earlier D0forD0ptbin cut with correct D0/D0bar mass)
  //  3: Bplus invariant mass cut without reconstruction vertex
  //  4: Everything except Bplus invariant mass cut
  //
  // returnvalue:
  //  1: selected
  //  0: selected (but Wrong-Sign pair)
  // -1: all other (rejected)

  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return -1;
  }

  if (!d) {
    cout << "candidate D0 null" << endl;
    return -1;
  }

  if (!pion) {
    cout << "candidate Pion null" << endl;
    return -1;
  }

  if (!primaryVertex) {
    cout << "primaryVertex null" << endl;
    return -1;
  }
    
  if(selLevel > 4 || selLevel < 0){
    AliWarning(Form("selLevel %d is not supported. return -1.",selLevel));
    return -1;
  }

  Int_t returnvaluePID = 1;
  Int_t returnvalueD0 = 1;
  Int_t returnvalueBplus = 1;

  // PID D0 daughters (check if the pions have the opposite charge)
  if(selLevel == 0 || selLevel == 1 || selLevel == 4){
    AliAODTrack* track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack* track1 = (AliAODTrack*)d->GetDaughter(1);

    if (pion->Charge() == -1) {
      if      ((SelectPID(track0, 2)) && (SelectPID(track1, 3))) returnvaluePID = 1;
      else if ((SelectPID(track0, 3)) && (SelectPID(track1, 2))) returnvaluePID = 0;
      else returnvaluePID = -1;
    } else if (pion->Charge() == 1) {
      if      ((SelectPID(track0, 3)) && (SelectPID(track1, 2))) returnvaluePID = 1;
      else if ((SelectPID(track0, 2)) && (SelectPID(track1, 3))) returnvaluePID = 0;
      else returnvaluePID = -1;
    } else {
      returnvaluePID = -1;
    }
  }

  // D0 window - invariant mass cut
  if(selLevel == 0 || selLevel == 2 || selLevel == 4){
    Double_t pdgMassD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Int_t ptbin = PtBinD0forD0ptbin(d->Pt());
    if(ptbin == -1) returnvalueD0 = -1;

    if (pion->Charge() == -1) {
      Double_t invMassDifference = TMath::Abs(d->InvMassD0() - pdgMassD0);
      if(!ApplyCutOnVariableD0forD0ptbinMVA(0, ptbin, invMassDifference)) returnvalueD0 = -1;
    } else if (pion->Charge() == 1) {
      Double_t invMassDifference = TMath::Abs(d->InvMassD0bar() - pdgMassD0);
      if(!ApplyCutOnVariableD0forD0ptbinMVA(0, ptbin, invMassDifference)) returnvalueD0 = -1;
    } else {
      returnvalueD0 = -1;
    }
  }

  // Bplus invariant mass window (using PV instead of SV, to speed thing up)
  if(selLevel == 0 || selLevel == 3){
    AliExternalTrackParam firstTrack;
    firstTrack.CopyFromVTrack(pion);
    AliExternalTrackParam secondTrack;
    secondTrack.CopyFromVTrack(d);
    
    Double_t d0z0[2], covd0z0[3], d0[2], d0err[2];
    firstTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0]; d0err[0] = TMath::Sqrt(covd0z0[0]);
    secondTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[1] = d0z0[0]; d0err[1] = TMath::Sqrt(covd0z0[0]);
    
    Double_t px[2],py[2],pz[2],momentum[3];
    firstTrack.GetPxPyPz(momentum);
    px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
    secondTrack.GetPxPyPz(momentum);
    px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];
    
    Double_t xdummy = 0., ydummy = 0.;
    Double_t dca = secondTrack.GetDCA(&firstTrack, bz, xdummy, ydummy);
    
    AliAODRecoDecayHF2Prong trackBPlus(primaryVertex, px, py, pz, d0, d0err, dca);
    
    UInt_t pdg2[2] = {211,421};
    Double_t invMassBPlus = trackBPlus.InvMass(2,pdg2);
    Double_t mBPlusPDG = TDatabasePDG::Instance()->GetParticle(521)->Mass();
    Double_t invMassDifference = TMath::Abs(mBPlusPDG - invMassBPlus) - 0.1; //-0.1 to be on the safe side as invmass is not completely correct here, will be recomputed in IsSelected()

    Int_t ptbin = PtBin(trackBPlus.Pt());
    if (ptbin == -1) returnvalueBplus = -1;
    else{
      if(!ApplyCutOnVariableMVA(39, ptbin, invMassDifference)) returnvalueBplus = -1;
    }
  }

  if(returnvaluePID == 1 && returnvalueD0 == 1 && returnvalueBplus == 1) return 1;
  else if(returnvaluePID == 0 && returnvalueD0 == 1 && returnvalueBplus == 1) return 0;
  else return -1;

}
//----------------------------------------------------------------------------------
Bool_t AliRDHFCutsBPlustoD0Pi::IsInFiducialAcceptance(Double_t /*pt*/, Double_t y) const
{
  //
  // Should be improved (for example pT dependence)
  //

  if (y > GetFiducialYCut()) return kFALSE;

  return kTRUE;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsBPlustoD0Pi::IsSelectedPID(AliAODRecoDecayHF* /*obj*/)
{
  //
  // PID method, n sigma approach default // not used for BPlus, done seperately for each daughter
  //

  // AliAODRecoDecayHF2Prong* dstar = (AliAODRecoDecayHF2Prong*)obj;
  // if(!dstar){
  //   cout<<"AliAODRecoDecayHF2Prong null"<<endl;
  //   return 0;
  // }

  // if(!fUsePID || dstar->Pt() > fMaxPtPid) return 3;

  // AliAODRecoDecayHF2Prong* d0 = (AliAODRecoDecayHF2Prong*)dstar->Get2Prong();
  // if(!d0){
  //   cout<<"AliAODRecoDecayHF2Prong null"<<endl;
  //   return 0;
  // }

  // //  here the PID
  // AliAODTrack *pos = (AliAODTrack*)dstar->Get2Prong()->GetDaughter(0);
  // AliAODTrack *neg = (AliAODTrack*)dstar->Get2Prong()->GetDaughter(1);

  // if (dstar->Charge()>0){
  //   if(!SelectPID(pos,2)) return 0;//pion+
  //   if(!SelectPID(neg,3)) return 0;//kaon-
  // }else{
  //   if(!SelectPID(pos,3)) return 0;//kaon+
  //   if(!SelectPID(neg,2)) return 0;//pion-
  // }

  // if ((fPidHF->GetMatch() == 10 || fPidHF->GetMatch() == 11) && fPidHF->GetITS()) { //ITS n sigma band
  //   AliAODTrack *softPion = (AliAODTrack*) dstar->GetBachelor();

  //   if (fPidHF->CheckBands(AliPID::kPion, AliPIDResponse::kITS, softPion) == -1) {
  //     return 0;
  //   }
  // }

  return 3;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsBPlustoD0Pi::SelectPID(AliAODTrack *track, Int_t type)
{
  //
  //  here the PID

  Bool_t isParticle = kTRUE;
  if (!fUsePID) return isParticle;
  Int_t match = fPidHF->GetMatch();

  if (match == 1) { //n-sigma
    Bool_t TPCon = TMath::Abs(2) > 1e-4 ? kTRUE : kFALSE;
    Bool_t TOFon = TMath::Abs(3) > 1e-4 ? kTRUE : kFALSE;

    Bool_t isTPC = kTRUE;
    Bool_t isTOF = kTRUE;

    if (TPCon) { //TPC
      if (fPidHF->CheckStatus(track, "TPC")) {
        if (type == 2) isTPC = fPidHF->IsPionRaw(track, "TPC");
        if (type == 3) isTPC = fPidHF->IsKaonRaw(track, "TPC");
      }
    }
    if (TOFon) { //TOF
      if (fPidHF->CheckStatus(track, "TOF")) {
        if (type == 2) isTOF = fPidHF->IsPionRaw(track, "TOF");
        if (type == 3) isTOF = fPidHF->IsKaonRaw(track, "TOF");
      }
    }

    //--------------------------------
    // cut on high momentum in the TPC
    //--------------------------------
    Double_t pPIDcut = track->P();
    if (pPIDcut > fTPCflag) isTPC = 1;

    isParticle = isTPC && isTOF;
  }

  if (match == 2) { //bayesian
    //Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
    Double_t prob[5] = {1., 1., 1., 1., 1.};

    //fPidHF->SetPriors(priors,5);
    //    fPidHF->BayesianProbability(track,prob);

    Double_t max = 0.;
    Int_t k = -1;
    for (Int_t i = 0; i < 5; i++) {
      if (prob[i] > max) {k = i; max = prob[i];}
    }
    isParticle = Bool_t(k == type);
  }

  if (match == 10 || match == 11) { //Assymetric PID using histograms
    Int_t checkTPC = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTPC, track);
    Int_t checkTOF = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTOF, track);

    isParticle = checkTPC >= 0 && checkTOF >= 0 ? kTRUE : kFALSE; //Standard: both at least compatible
    if (match == 11) { //Extra requirement: at least one identified
      isParticle = isParticle && checkTPC + checkTOF >= 1 ? kTRUE : kFALSE;
    }
  }

  if (match == 12) { //Circular cut
    Double_t nSigmaTOF = 0;
    Double_t nSigmaTPC = 0;

    Double_t radius = fCircRadius;

    isParticle = kTRUE;
    if (radius > 0) {
      Int_t TPCok = fPidHF->GetnSigmaTPC(track, type, nSigmaTPC);
      Int_t TOFok = fPidHF->GetnSigmaTOF(track, type, nSigmaTOF);
      if (TPCok != -1 && TOFok != -1) {
        //Both detectors gave PID information
        isParticle = TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF) <= radius ? kTRUE : kFALSE;
      }
      else {
        //Maximum one detector gave PID information
        if (TPCok != -1) {
          isParticle = nSigmaTPC <= radius ? kTRUE : kFALSE;
        }
        if (TOFok != -1) {
          isParticle = nSigmaTOF <= radius ? kTRUE : kFALSE;
        }
      }
    }
  }

  return isParticle;

}
//-------------------------------------------------------------------------------------
Double_t AliRDHFCutsBPlustoD0Pi::DeltaInvMassBPlusKpipi(AliAODRecoDecayHF2Prong * BPlus) const
{
  ///
  /// 3 prong invariant mass of the D0 daughters, the soft pion, and the BPlus pion
  ///

  AliAODRecoDecayHF2Prong * D0 = (AliAODRecoDecayHF2Prong*)BPlus->GetDaughter(1);

  Double_t e[3] = {0};
  if (BPlus->Charge() == -1)
  {
    e[0] = D0->EProng(0, 211);
    e[1] = D0->EProng(1, 321);
  } else if (BPlus->Charge() == 1) {
    e[0] = D0->EProng(0, 321);
    e[1] = D0->EProng(1, 211);
  }
  e[2] = BPlus->EProng(0, 211);

  Double_t esum = e[0] + e[1] + e[2];
  Double_t invMassBPlus = TMath::Sqrt(esum * esum - BPlus->P2());

  Double_t invMassD0 = -1;

  if (BPlus->Charge() == -1) {invMassD0 = D0->InvMassD0();}
  else {invMassD0 = D0->InvMassD0bar();}
  if (invMassD0 == -1) {cout << "wrong invmass delta D0 BPlus" << endl;}

  return invMassBPlus - invMassD0;
}


//---------------------------------------------------------------------------
//
//  DO for D0 pt bin functions
//
//---------------------------------------------------------------------------

void AliRDHFCutsBPlustoD0Pi::SetCutsD0forD0ptbin(Int_t nVars, Int_t nPtBins, Float_t **cutsRDD0forD0ptbin) {
  //
  // store the cuts
  //

  if (nVars != fnVarsD0forD0ptbin) {
    printf("Wrong number of variables: it has to be %d\n", fnVarsD0forD0ptbin);
    AliFatal("exiting");
  }
  if (nPtBins != fnPtBinsD0forD0ptbin) {
    printf("Wrong number of pt bins: it has to be %d\n", fnPtBinsD0forD0ptbin);
    AliFatal("exiting");
  }

  if (!fCutsRDD0forD0ptbin)  fCutsRDD0forD0ptbin = new Float_t[fGlobalIndexD0forD0ptbin];


  for (Int_t iv = 0; iv < fnVarsD0forD0ptbin; iv++)
  {
    for (Int_t ib = 0; ib < fnPtBinsD0forD0ptbin; ib++)
    {
      //check

      if (GetGlobalIndexD0forD0ptbin(iv, ib) >= fGlobalIndexD0forD0ptbin)
      {
        cout << "Overflow, exit..." << endl;
        return;
      }

      fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(iv, ib)] = cutsRDD0forD0ptbin[iv][ib];

    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetCutsD0forD0ptbin(Int_t glIndex, Float_t *cutsRDD0forD0ptbin) {
  //
  // store the cuts
  //
  if (glIndex != fGlobalIndexD0forD0ptbin) {
    cout << "Wrong array size: it has to be " << fGlobalIndexD0forD0ptbin << endl;
    AliFatal("exiting");
  }
  if (!fCutsRDD0forD0ptbin)  fCutsRDD0forD0ptbin = new Float_t[fGlobalIndexD0forD0ptbin];

  for (Int_t iGl = 0; iGl < fGlobalIndexD0forD0ptbin; iGl++) {
    fCutsRDD0forD0ptbin[iGl] = cutsRDD0forD0ptbin[iGl];
  }
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::PtBinD0forD0ptbin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin = -1;
  if (pt < fPtBinLimitsD0forD0ptbin[0])return ptbin;
  for (Int_t i = 0; i < fnPtBinsD0forD0ptbin; i++) {
    if (pt < fPtBinLimitsD0forD0ptbin[i + 1]) {
      ptbin = i;
      break;
    }
  }
  return ptbin;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetPtBinsD0forD0ptbin(Int_t nPtBinLimits, Float_t *ptBinLimits) {
  // Set the pt bins

  if (fPtBinLimitsD0forD0ptbin) {
    delete [] fPtBinLimitsD0forD0ptbin;
    fPtBinLimitsD0forD0ptbin = NULL;
    printf("Changing the pt bins\n");
  }

  if (nPtBinLimits != fnPtBinsD0forD0ptbin + 1) {
    cout << "Warning: ptBinLimits dimension " << nPtBinLimits << " != nPtBins+1 (" << fnPtBinsD0forD0ptbin + 1 << ")\nSetting nPtBins to " << nPtBinLimits - 1 << endl;
    SetNPtBinsD0forD0ptbin(nPtBinLimits - 1);
  }

  fnPtBinLimitsD0forD0ptbin = nPtBinLimits;
  SetGlobalIndexD0forD0ptbin();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimitsD0forD0ptbin = new Float_t[fnPtBinLimitsD0forD0ptbin];
  for (Int_t ib = 0; ib < nPtBinLimits; ib++) fPtBinLimitsD0forD0ptbin[ib] = ptBinLimits[ib];

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[75]) {

  if (GetIsCutUsed(nCutIndex, ptbin) == kTRUE)
  {
    Bool_t bCut = kFALSE;
    if (GetIsUpperCut(nCutIndex) == kTRUE)
    {
      if (cutVariableValue > fCutsRD[GetGlobalIndex(nCutIndex, ptbin)]) bCut = kTRUE;
    } else
    {
      if (cutVariableValue < fCutsRD[GetGlobalIndex(nCutIndex, ptbin)]) bCut = kTRUE;
    }
    if (bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::ApplyCutOnVariableD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[29]) {

  // std::cout << "index: " << nCutIndex << ", ptbin: " << ptbin << ", used: " << GetIsCutUsedD0forD0ptbin(nCutIndex,ptbin) << ", upper: " << GetIsUpperCutD0forD0ptbin(nCutIndex) << ", value: " << cutVariableValue << ", cutvalue: " << fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)] << std::endl;
  if (GetIsCutUsedD0forD0ptbin(nCutIndex, ptbin) == kTRUE)
  {
    Bool_t bCut = kFALSE;
    if (GetIsUpperCutD0forD0ptbin(nCutIndex) == kTRUE)
    {
      if (cutVariableValue > fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptbin)]) bCut = kTRUE;
    } else
    {
      if (cutVariableValue < fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptbin)]) bCut = kTRUE;
    }
    if (bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::ApplyCutOnVariableMVA(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue) {

  if (GetIsCutUsed(nCutIndex, ptbin) == kTRUE)
  {
    if (GetIsUpperCut(nCutIndex) == kTRUE)
    {
      if (cutVariableValue > fCutsRD[GetGlobalIndex(nCutIndex, ptbin)]) return 0;
    } else
    {
      if (cutVariableValue < fCutsRD[GetGlobalIndex(nCutIndex, ptbin)]) return 0;
    }
  }
  return 1;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsBPlustoD0Pi::ApplyCutOnVariableD0forD0ptbinMVA(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue) {

  if (GetIsCutUsedD0forD0ptbin(nCutIndex, ptbin) == kTRUE)
  {
    if (GetIsUpperCutD0forD0ptbin(nCutIndex) == kTRUE)
    {
      if (cutVariableValue > fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptbin)]) return 0;
    } else
    {
      if (cutVariableValue < fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptbin)]) return 0;
    }
  }
  return 1;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetVarNamesD0forD0ptbin(Int_t nVars, TString *varNames, Bool_t *isUpperCut) {
  // Set the variable names

  if (fVarNamesD0forD0ptbin) {
    delete [] fVarNamesD0forD0ptbin;
    fVarNamesD0forD0ptbin = NULL;
    //printf("Changing the variable names\n");
  }
  if (nVars != fnVarsD0forD0ptbin) {
    printf("Wrong number of variables: it has to be %d\n", fnVarsD0forD0ptbin);
    return;
  }
  //fnVars=nVars;
  fVarNamesD0forD0ptbin = new TString[nVars];
  fIsUpperCutD0forD0ptbin = new Bool_t[nVars];
  for (Int_t iv = 0; iv < nVars; iv++) {
    fVarNamesD0forD0ptbin[iv] = varNames[iv];
    fIsUpperCutD0forD0ptbin[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::InitializeCuts() {
  // Set the use of all cuts to kFALSE and the cut value to zero. This function has to be called after setting the pt bins.

  if (fIsCutUsed) {
    delete [] fIsCutUsed;
    fIsCutUsed = NULL;
  }

  fIsCutUsed = new Bool_t[(GetNPtBins()) * (GetNVars())];
  for (Int_t iv = 0; iv < (GetNPtBins()) * (GetNVars()); iv++) {
    fIsCutUsed[iv] = kFALSE;
  }

  if (!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];


  for (Int_t iv = 0; iv < fnVars; iv++) {

    for (Int_t ib = 0; ib < fnPtBins; ib++) {

      //check
      if (GetGlobalIndex(iv, ib) >= fGlobalIndex) {
        cout << "Overflow, exit..." << endl;
        return;
      }

      fCutsRD[GetGlobalIndex(iv, ib)] = 0;

    }
  }

  // D0 for D0 pt bin

  if (fIsCutUsedD0forD0ptbin) {
    delete [] fIsCutUsedD0forD0ptbin;
    fIsCutUsedD0forD0ptbin = NULL;
  }

  fIsCutUsedD0forD0ptbin = new Bool_t[(GetNPtBinsD0forD0ptbin()) * (GetNVarsD0forD0ptbin())];
  for (Int_t iv = 0; iv < (GetNPtBinsD0forD0ptbin()) * (GetNVarsD0forD0ptbin()); iv++) {
    fIsCutUsedD0forD0ptbin[iv] = kFALSE;

  }

  if (!fCutsRDD0forD0ptbin)  fCutsRDD0forD0ptbin = new Float_t[fGlobalIndexD0forD0ptbin];

  for (Int_t iv = 0; iv < fnVarsD0forD0ptbin; iv++)
  {
    for (Int_t ib = 0; ib < fnPtBinsD0forD0ptbin; ib++)
    {
      //check

      if (GetGlobalIndexD0forD0ptbin(iv, ib) >= fGlobalIndexD0forD0ptbin)
      {
        cout << "Overflow, exit..." << endl;
        return;
      }

      fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(iv, ib)] = 0;
    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t cutValue) {
  // Set the cut value and direction

  fIsCutUsed[GetGlobalIndex(nCutIndex, ptBin)] = kTRUE;
  fIsUpperCut[nCutIndex] = cutDirection;
  fCutsRD[GetGlobalIndex(nCutIndex, ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetCutD0forD0ptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t cutValue) {
  // Set the cut value and direction

  fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptBin)] = kTRUE;
  fIsUpperCutD0forD0ptbin[nCutIndex] = cutDirection;
  fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex, ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::InitializeCutsForCutOptimization(Int_t nCutsForOptimization, Int_t nVariables) {
  //
  // Initialize the cuts for optimization
  //

  fnVariablesForCutOptimization = nVariables;
  fnCutsForOptimization = nCutsForOptimization;
  SetGlobalIndexForCutOptimization();


  if (fIsUpperCutForCutOptimization) 
  {
    delete [] fIsUpperCutForCutOptimization;
    fIsUpperCutForCutOptimization = nullptr;
  }
  fIsUpperCutForCutOptimization = new Bool_t[fnVariablesForCutOptimization];

  if (fCutIndexForCutOptimization) 
  {
    delete [] fCutIndexForCutOptimization;
    fCutIndexForCutOptimization = nullptr;
  }
  fCutIndexForCutOptimization = new Int_t[fnVariablesForCutOptimization];

  if (fSigmaForCutOptimization) 
  {
    delete [] fSigmaForCutOptimization;
    fSigmaForCutOptimization = nullptr;
  }
  fSigmaForCutOptimization = new Float_t[fnPtBins];

  if (!fCutsRDForCutOptimization)  fCutsRDForCutOptimization = new Float_t[fGlobalIndexCutOptimization];


  for (Int_t iv = 0; iv < fnVariablesForCutOptimization; iv++)
  {
    for (Int_t ib = 0; ib < fnPtBins; ib++)
    {
      for (Int_t ic = 0; ic < fnCutsForOptimization; ++ic)
      {
        if (GetGlobalIndexForCutOptimization(ic, iv, ib) >= fGlobalIndexCutOptimization)
        {
          cout << "Overflow, exit..." << endl;
          return;
        }

        fCutsRDForCutOptimization[GetGlobalIndexForCutOptimization(ic, iv, ib)] = 0;
      }
    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetCutsForCutOptimization(Int_t glIndex, Float_t *cutsRDForCutOptimization) {
  //
  // store the cuts
  //
  if (glIndex != fGlobalIndexCutOptimization) {
    cout << "Wrong array size: it has to be " << fGlobalIndexCutOptimization << endl;
    AliFatal("exiting");
  }
  if (!fCutsRDForCutOptimization)  fCutsRDForCutOptimization = new Float_t[fGlobalIndexCutOptimization];

  for (Int_t iGl = 0; iGl < fGlobalIndexCutOptimization; iGl++) {
    fCutsRDForCutOptimization[iGl] = cutsRDForCutOptimization[iGl];
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsBPlustoD0Pi::SetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t * cutValues) {
  // Set the cut value and direction
  if (nVariable > GetnVariablesForCutOptimization() || nVariable < 0) cout << "wrong variable number for GetNCutsForOptimization" << endl;

  fIsUpperCutForCutOptimization[nVariable] = cutDirection;
  fCutIndexForCutOptimization[nVariable] = nCutIndex;

  for (Int_t iCut = 0; iCut < fnCutsForOptimization; ++iCut)
  {
    fCutsRDForCutOptimization[GetGlobalIndexForCutOptimization(iCut, nVariable, ptBin)] = cutValues[iCut];
  }

  return;
}
