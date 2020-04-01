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

/* $Id: AliRDHFCutsB0toDPi.cxx $ */

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed B0->DPlusPi->KPiPi
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
#include "AliAODRecoDecayHF3Prong.h"
#include "AliRDHFCutsB0toDPi.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsB0toDPi);
/// \endcond



//--------------------------------------------------------------------------
AliRDHFCutsB0toDPi::AliRDHFCutsB0toDPi(const char* name) :
  AliRDHFCuts(name),
  fMaxPtPid(9999.),
  fTPCflag(999.),
  fCircRadius(0.),

  fIsCutUsed(0x0),

  fnVarsDPlusforDPlusptbin(0),
  fnPtBinsDPlusforDPlusptbin(1),
  fGlobalIndexDPlusforDPlusptbin(0),
  fCutsRDDPlusforDPlusptbin(0x0),
  fnPtBinLimitsDPlusforDPlusptbin(1),
  fPtBinLimitsDPlusforDPlusptbin(0x0),
  fIsUpperCutDPlusforDPlusptbin(0x0),
  fIsCutUsedDPlusforDPlusptbin(0x0),
  fVarNamesDPlusforDPlusptbin(0x0),

  fMinITSNclsDaughterType(),
  fMinTPCNclsDaughterType(),
  fUseITSRefitDaughterType(),
  fUseTPCRefitDaughterType(),
  fUseFilterBitDaughterType(),
  fFilterBitDaughterType(),
  fMinPtDaughterType(),
  fMaxAbsEtaDaughterType(),
  fHardSelectionArrayITSDaughterType(),
  fSoftSelectionArrayITSDaughterType(),
  fNSoftITSCutDaughterType(),
  fMind0DaughterType(),
  fMinNormd0DaughterType(),

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

  // Main cut setup as function of B0 pt bins
  const Int_t nvars = 78;
  SetNVars(nvars);

  TString varNames[nvars];
  Int_t iterator = 0;

  //DPlus cut variables
  varNames[iterator++] =   /*-00-*/ "inv. mass width[GeV]";
  varNames[iterator++] =   /*-01-*/ "delta mass width  [GeV]";
  varNames[iterator++] =   /*-02-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++] =   /*-03-*/ "pointing angle XY [Cos(theta)]";
  varNames[iterator++] =   /*-04-*/ "dca12 [cm]";
  varNames[iterator++] =   /*-05-*/ "dca13 [cm]";
  varNames[iterator++] =   /*-06-*/ "dca23 [cm]";
  varNames[iterator++] =   /*-07-*/ "Pt DPlus [GeV/c]";
  varNames[iterator++] =   /*-08-*/ "Pt first daughter [GeV/c]";
  varNames[iterator++] =   /*-09-*/ "Pt second daughter [GeV/c]";
  varNames[iterator++] =   /*-10-*/ "Pt third daughter [GeV/c]";
  varNames[iterator++] =   /*-11-*/ "d0 DPlus [cm]";
  varNames[iterator++] =   /*-12-*/ "d0 first daughter [cm]";
  varNames[iterator++] =   /*-13-*/ "d0 second daughter [cm]";
  varNames[iterator++] =   /*-14-*/ "d0 third daughter [cm]";
  varNames[iterator++] =   /*-15-*/ "d0d0d0 [cm^3]";
  varNames[iterator++] =   /*-16-*/ "d0d012 [cm^2]";
  varNames[iterator++] =   /*-17-*/ "d0d013 [cm^2]";
  varNames[iterator++] =   /*-18-*/ "d0d023 [cm^2]";
  varNames[iterator++] =   /*-19-*/ "d0d0 XY [cm^2]";

  varNames[iterator++] =   /*-20-*/ "smallest Angle Mother Daughters";
  varNames[iterator++] =   /*-21-*/ "largest Angle Mother Daughter";
  varNames[iterator++] =   /*-22-*/ "angle difference";
  varNames[iterator++] =   /*-23-*/ "vertexDistance [cm]";
  varNames[iterator++] =   /*-24-*/ "vertex distance XY [cm]";
  varNames[iterator++] =   /*-25-*/ "pseudoProperDecayTime";
  varNames[iterator++] =   /*-26-*/ "DecayTime";
  varNames[iterator++] =   /*-27-*/ "normalizedDecayTime";
  varNames[iterator++] =   /*-28-*/ "normDecayLength";
  varNames[iterator++] =   /*-29-*/ "normDecayLength XY";
  varNames[iterator++] =   /*-30-*/ "Chi2 per NDF vertex";

  varNames[iterator++] =   /*-31-*/ "Normalized d0 DPlus First pion [cm]";
  varNames[iterator++] =   /*-32-*/ "Normalized d0 DPlus kaon [cm]";
  varNames[iterator++] =   /*-33-*/ "Normalized d0 DPlus Second pion [cm]";
  varNames[iterator++] =   /*-34-*/ "Normalized d0 DPlus [cm]";
  varNames[iterator++] =   /*-35-*/ "Normalized impact product DPlus [cm^3]";
  varNames[iterator++] =   /*-36-*/ "Dist12 [cm]";
  varNames[iterator++] =   /*-37-*/ "Dist23 [cm]";
  varNames[iterator++] =   /*-38-*/ "Sigma vertex";

  varNames[iterator++] =   /*-39-*/ "pointingAngleToB0";
  varNames[iterator++] =   /*-40-*/ "d0MotherToB0";
  varNames[iterator++] =   /*-41-*/ "d0FirstDaughterToB0";
  varNames[iterator++] =   /*-42-*/ "d0SecondDaughterToB0";
  varNames[iterator++] =   /*-43-*/ "d0ThirdDaughterToB0";
  varNames[iterator++] =   /*-44-*/ "impactProductToB0";
  varNames[iterator++] =   /*-45-*/ "impactProductXYToB0";
  varNames[iterator++] =   /*-46-*/ "normDecayLengthToB0";
  varNames[iterator++] =   /*-47-*/ "pseudoProperDecayTimeToB0";
  varNames[iterator++] =   /*-48-*/ "DecayTimeToB0";
  varNames[iterator++] =   /*-49-*/ "normalizedDecayTimeToB0";

  //B0 cut variables
  varNames[iterator++] =   /*-50-*/ "inv. mass width[GeV]";
  varNames[iterator++] =   /*-51-*/ "delta mass width  [GeV]";
  varNames[iterator++] =   /*-52-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++] =   /*-53-*/ "pointing angle XY [Cos(theta)]";
  varNames[iterator++] =   /*-54-*/ "dca [cm]";
  varNames[iterator++] =   /*-55-*/ "Pt B0 [GeV/c]";
  varNames[iterator++] =   /*-56-*/ "Pt DPlus [GeV/c]";
  varNames[iterator++] =   /*-57-*/ "Pt Pion [GeV/c]";
  varNames[iterator++] =   /*-58-*/ "d0 B0 [cm]";
  varNames[iterator++] =   /*-59-*/ "d0 DPlus [cm]";
  varNames[iterator++] =   /*-60-*/ "d0 Pion [cm]";
  varNames[iterator++] =   /*-61-*/ "d0d0 [cm^2]";
  varNames[iterator++] =   /*-62-*/ "d0d0 XY [cm^2]";
  varNames[iterator++] =   /*-63-*/ "angleMotherFirstDaughter";
  varNames[iterator++] =   /*-64-*/ "angleMotherSecondDaughter";
  varNames[iterator++] =   /*-65-*/ "angleBetweenBothDaughters";
  varNames[iterator++] =   /*-66-*/ "cosThetaStar";

  varNames[iterator++] =   /*-67-*/ "vertexDistance [cm]";
  varNames[iterator++] =   /*-68-*/ "vertex distance XY [cm]";
  varNames[iterator++] =   /*-69-*/ "pseudoProperDecayTime";
  varNames[iterator++] =   /*-70-*/ "DecayTime";
  varNames[iterator++] =   /*-71-*/ "normalizedDecayTime";
  varNames[iterator++] =   /*-72-*/ "normDecayLength";
  varNames[iterator++] =   /*-73-*/ "normDecayLength XY";
  varNames[iterator++] =   /*-74-*/ "Chi2 per NDF vertex";

  varNames[iterator++] =   /*-75-*/ "Normalized d0 B0 pion [cm]";
  varNames[iterator++] =   /*-76-*/ "Normalized d0 B0 [cm]";
  varNames[iterator++] =   /*-77-*/ "Normalized impact product B0 [cm^2]";

  Bool_t isUpperCut[nvars] = {0};

  SetVarNames(nvars, varNames, isUpperCut);

  Float_t limits[2] = {0, 999999999.};
  SetPtBins(2, limits);


  //
  // Initialization of DPlus cuts for DPlus pt bins
  //

  const Int_t nvarsDPlusforDPlusptbin = 39;
  SetNVarsDPlusforDPlusptbin(nvarsDPlusforDPlusptbin);

  TString varNamesDPlusforDPlusptbin[nvarsDPlusforDPlusptbin];
  iterator = 0;

  //DPlus cut variables
  varNamesDPlusforDPlusptbin[iterator++] =   /*-00-*/ "inv. mass width[GeV]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-01-*/ "delta mass width  [GeV]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-02-*/ "pointing angle [Cos(theta)]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-03-*/ "pointing angle XY [Cos(theta)]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-04-*/ "dca12 [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-05-*/ "dca13 [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-06-*/ "dca23 [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-07-*/ "Pt DPlus [GeV/c]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-08-*/ "Pt first daughter [GeV/c]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-09-*/ "Pt second daughter [GeV/c]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-10-*/ "Pt third daughter [GeV/c]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-11-*/ "d0 DPlus [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-12-*/ "d0 first daughter [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-13-*/ "d0 second daughter [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-14-*/ "d0 third daughter [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-15-*/ "d0d0d0 [cm^3]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-16-*/ "d0d012 [cm^2]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-17-*/ "d0d013 [cm^2]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-18-*/ "d0d023 [cm^2]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-19-*/ "d0d0 XY [cm^2]";

  varNamesDPlusforDPlusptbin[iterator++] =   /*-20-*/ "smallest Angle Mother Daughters";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-21-*/ "largest Angle Mother Daughter";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-22-*/ "angle difference";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-23-*/ "vertexDistance [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-24-*/ "vertex distance XY [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-25-*/ "pseudoProperDecayTime";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-26-*/ "DecayTime";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-27-*/ "normalizedDecayTime";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-28-*/ "normDecayLength";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-29-*/ "normDecayLength XY";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-30-*/ "Chi2 per NDF vertex";

  varNamesDPlusforDPlusptbin[iterator++] =   /*-31-*/ "Normalized d0 DPlus First pion [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-32-*/ "Normalized d0 DPlus kaon [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-33-*/ "Normalized d0 DPlus Second pion [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-34-*/ "Normalized d0 DPlus [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-35-*/ "Normalized impact product DPlus [cm^3]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-36-*/ "Dist12 [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-37-*/ "Dist23 [cm]";
  varNamesDPlusforDPlusptbin[iterator++] =   /*-38-*/ "Sigma vertex";

  Bool_t isUpperCutDPlusforDPlusptbin[nvarsDPlusforDPlusptbin] = {0};

  SetVarNamesDPlusforDPlusptbin(nvarsDPlusforDPlusptbin, varNamesDPlusforDPlusptbin, isUpperCutDPlusforDPlusptbin);

  Float_t limitsDPlusforDPlusptbin[2] = {0, 999999999.};
  SetPtBinsDPlusforDPlusptbin(2, limitsDPlusforDPlusptbin);

  Bool_t forOpt[16] = {0}; //not used for B0 analysis
  SetVarsForOpt(16, forOpt);

}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDPi::AliRDHFCutsB0toDPi(const AliRDHFCutsB0toDPi &source) :
  AliRDHFCuts(source),
  fMaxPtPid(source.fMaxPtPid),
  fTPCflag(source.fTPCflag),
  fCircRadius(source.fCircRadius),
  fIsCutUsed(0x0),

  fnVarsDPlusforDPlusptbin(source.fnVarsDPlusforDPlusptbin),
  fnPtBinsDPlusforDPlusptbin(source.fnPtBinsDPlusforDPlusptbin),
  fGlobalIndexDPlusforDPlusptbin(source.fGlobalIndexDPlusforDPlusptbin),
  fCutsRDDPlusforDPlusptbin(0x0),
  fnPtBinLimitsDPlusforDPlusptbin(source.fnPtBinLimitsDPlusforDPlusptbin),
  fPtBinLimitsDPlusforDPlusptbin(0x0),
  fIsUpperCutDPlusforDPlusptbin(0x0),
  fIsCutUsedDPlusforDPlusptbin(0x0),
  fVarNamesDPlusforDPlusptbin(0x0),

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

  if (source.fPtBinLimitsDPlusforDPlusptbin) SetPtBinsDPlusforDPlusptbin(source.fnPtBinLimitsDPlusforDPlusptbin, source.fPtBinLimitsDPlusforDPlusptbin);
  if (source.fVarNamesDPlusforDPlusptbin) SetVarNamesDPlusforDPlusptbin(source.fnVarsDPlusforDPlusptbin, source.fVarNamesDPlusforDPlusptbin, source.fIsUpperCut);
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
  if (source.fIsCutUsedDPlusforDPlusptbin)
  {
    if (fIsCutUsedDPlusforDPlusptbin) {
      delete [] fIsCutUsedDPlusforDPlusptbin;
      fIsCutUsedDPlusforDPlusptbin = nullptr;
    }
    fIsCutUsedDPlusforDPlusptbin = new Bool_t[(source.GetNPtBinsDPlusforDPlusptbin()) * (source.GetNVarsDPlusforDPlusptbin())];
    for (Int_t i = 0; i < source.fnVarsDPlusforDPlusptbin; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBinsDPlusforDPlusptbin; j++)
      {
        Bool_t bUse = source.GetIsCutUsedDPlusforDPlusptbin(i, j);
        SetIsCutUsedDPlusforDPlusptbin(i, j, bUse);
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
  if (source.fCutsRDDPlusforDPlusptbin) SetCutsDPlusforDPlusptbin(source.fGlobalIndexDPlusforDPlusptbin, source.fCutsRDDPlusforDPlusptbin);
  if (source.fCutsRDForCutOptimization) SetCutsForCutOptimization(source.fGlobalIndexCutOptimization, source.fCutsRDForCutOptimization);

  Int_t nDaughterTypes = 3;
  for (int i = 0; i < nDaughterTypes; ++i)
  {
    SetMinITSNclsDaughterType(i, source.GetMinITSNclsDaughterType(i));
    SetMinTPCNclsDaughterType(i, source.GetMinTPCNclsDaughterType(i));
    SetUseITSRefitDaughterType(i, source.UseITSRefitDaughterType(i));
    SetUseTPCRefitDaughterType(i, source.UseTPCRefitDaughterType(i));
    SetUseFilterBitDaughterType(i, source.UseFilterBitDaughterType(i));
    SetFilterBitDaughterType(i, source.GetFilterBitDaughterType(i));
    SetMinPtDaughterType(i, source.GetMinPtDaughterType(i));
    SetMaxAbsEtaDaughterType(i, source.GetMaxAbsEtaDaughterType(i));
    SetNSoftITSCutDaughterType(i, source.GetNSoftITSCutDaughterType(i));
    SetMind0DaughterType(i, source.GetMind0DaughterType(i));
    SetMinNormd0DaughterType(i, source.GetMinNormd0DaughterType(i));

    Bool_t bHardArray[7] = {0}; Bool_t bSoftArray[7] = {0};
    source.GetHardSelectionArrayITSDaughterType(i, bHardArray);
    source.GetSoftSelectionArrayITSDaughterType(i, bSoftArray);
    SetHardSelectionArrayITSDaughterType(i, bHardArray);
    SetSoftSelectionArrayITSDaughterType(i, bSoftArray);
  }

}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDPi::~AliRDHFCutsB0toDPi() {
  //
  // Default Destructor
  //

  if (fIsCutUsed) { delete [] fIsCutUsed; fIsCutUsed = nullptr; }
  if (fCutsRDDPlusforDPlusptbin) { delete [] fCutsRDDPlusforDPlusptbin; fCutsRDDPlusforDPlusptbin = nullptr; }
  if (fPtBinLimitsDPlusforDPlusptbin) { delete [] fPtBinLimitsDPlusforDPlusptbin; fPtBinLimitsDPlusforDPlusptbin = nullptr; }
  if (fIsUpperCutDPlusforDPlusptbin) { delete [] fIsUpperCutDPlusforDPlusptbin; fIsUpperCutDPlusforDPlusptbin = nullptr; }
  if (fIsCutUsedDPlusforDPlusptbin) { delete [] fIsCutUsedDPlusforDPlusptbin; fIsCutUsedDPlusforDPlusptbin = nullptr; }
  if (fVarNamesDPlusforDPlusptbin) { delete [] fVarNamesDPlusforDPlusptbin; fVarNamesDPlusforDPlusptbin = nullptr; }
  if (fCutsRDForCutOptimization) { delete [] fCutsRDForCutOptimization; fCutsRDForCutOptimization = nullptr; }
  if (fIsUpperCutForCutOptimization) { delete [] fIsUpperCutForCutOptimization; fIsUpperCutForCutOptimization = nullptr; }
  if (fCutIndexForCutOptimization) { delete [] fCutIndexForCutOptimization; fCutIndexForCutOptimization = nullptr; }
  if (fSigmaForCutOptimization) { delete [] fSigmaForCutOptimization; fSigmaForCutOptimization = nullptr; }
}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDPi &AliRDHFCutsB0toDPi::operator=(const AliRDHFCutsB0toDPi &source)
{
  //
  // assignment operator
  //

  if (&source == this) return *this;

  AliRDHFCuts::operator=(source);

  fMaxPtPid = source.fMaxPtPid;
  fTPCflag = source.fTPCflag;
  fCircRadius = source.fCircRadius;
  fnVarsDPlusforDPlusptbin = source.fnVarsDPlusforDPlusptbin;
  fnPtBinsDPlusforDPlusptbin = source.fnPtBinsDPlusforDPlusptbin;
  fGlobalIndexDPlusforDPlusptbin = source.fGlobalIndexDPlusforDPlusptbin;
  fnPtBinLimitsDPlusforDPlusptbin = source.fnPtBinLimitsDPlusforDPlusptbin;

  fFiducialYCut = source.fFiducialYCut;
  fnVariablesForCutOptimization = source.fnVariablesForCutOptimization;
  fnCutsForOptimization = source.fnCutsForOptimization;
  fGlobalIndexCutOptimization = source.fGlobalIndexCutOptimization;
  fNumberOfSigmaBinsForCutOptimization = source.fNumberOfSigmaBinsForCutOptimization;

  if (source.fPtBinLimitsDPlusforDPlusptbin) SetPtBinsDPlusforDPlusptbin(source.fnPtBinLimitsDPlusforDPlusptbin, source.fPtBinLimitsDPlusforDPlusptbin);
  if (source.fVarNamesDPlusforDPlusptbin) SetVarNamesDPlusforDPlusptbin(source.fnVarsDPlusforDPlusptbin, source.fVarNamesDPlusforDPlusptbin, source.fIsUpperCut);
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
  if (source.fIsCutUsedDPlusforDPlusptbin)
  {
    if (fIsCutUsedDPlusforDPlusptbin) {
      delete [] fIsCutUsedDPlusforDPlusptbin;
      fIsCutUsedDPlusforDPlusptbin = nullptr;
    }
    fIsCutUsedDPlusforDPlusptbin = new Bool_t[(source.GetNPtBinsDPlusforDPlusptbin()) * (source.GetNVarsDPlusforDPlusptbin())];
    for (Int_t i = 0; i < source.fnVarsDPlusforDPlusptbin; ++i)
    {
      for (Int_t j = 0; j < source.fnPtBinsDPlusforDPlusptbin; j++)
      {
        Bool_t bUse = source.GetIsCutUsedDPlusforDPlusptbin(i, j);
        SetIsCutUsedDPlusforDPlusptbin(i, j, bUse);
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
  if (source.fCutsRDDPlusforDPlusptbin) SetCutsDPlusforDPlusptbin(source.fGlobalIndexDPlusforDPlusptbin, source.fCutsRDDPlusforDPlusptbin);
  if (source.fCutsRDForCutOptimization) SetCutsForCutOptimization(source.fGlobalIndexCutOptimization, source.fCutsRDForCutOptimization);

  Int_t nDaughterTypes = 3;
  for (int i = 0; i < nDaughterTypes; ++i)
  {
    SetMinITSNclsDaughterType(i, source.GetMinITSNclsDaughterType(i));
    SetMinTPCNclsDaughterType(i, source.GetMinTPCNclsDaughterType(i));
    SetUseITSRefitDaughterType(i, source.UseITSRefitDaughterType(i));
    SetUseTPCRefitDaughterType(i, source.UseTPCRefitDaughterType(i));
    SetUseFilterBitDaughterType(i, source.UseFilterBitDaughterType(i));
    SetFilterBitDaughterType(i, source.GetFilterBitDaughterType(i));
    SetMinPtDaughterType(i, source.GetMinPtDaughterType(i));
    SetMaxAbsEtaDaughterType(i, source.GetMaxAbsEtaDaughterType(i));
    SetNSoftITSCutDaughterType(i, source.GetNSoftITSCutDaughterType(i));
    SetMind0DaughterType(i, source.GetMind0DaughterType(i));
    SetMinNormd0DaughterType(i, source.GetMinNormd0DaughterType(i));

    Bool_t bHardArray[7] = {0}; Bool_t bSoftArray[7] = {0};
    source.GetHardSelectionArrayITSDaughterType(i, bHardArray);
    source.GetSoftSelectionArrayITSDaughterType(i, bSoftArray);
    SetHardSelectionArrayITSDaughterType(i, bHardArray);
    SetSoftSelectionArrayITSDaughterType(i, bSoftArray);
  }

  return *this;
}
//--------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::GetCutVarsForOpt(AliAODRecoDecayHF* /*d*/, Float_t* /*vars*/, Int_t /*nvars*/, Int_t* /*pdgdaughters*/) {
  // not used - alternative cut optimization method below

  return;
}
//--------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDPi::IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[78]) {
  //
  // In this function we apply the selection cuts on the B0 candidate and its daughters.
  // The function returns 0 if the candidate is cut and is able to return information on which cuts the candidate passes.
  //

  fIsSelectedCuts = 0;
  fIsSelectedPID = 0;

  // The cuts used in this class have to be set in the maketfile.
  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateB0 = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateB0) {
    cout << "candidateB0 null" << endl;
    return 0;
  }

  AliAODRecoDecayHF3Prong* candidateDPlus = (AliAODRecoDecayHF3Prong*)candidateB0->GetDaughter(0);
  if (!candidateDPlus) {
    std::cout << "candidateDPlus null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateB0Pion = (AliAODTrack*)candidateB0->GetDaughter(1);
  if (!candidateB0Pion) {
    std::cout << "candidateB0Pion null" << std::endl;
    return 0;
  }

  // AliAODTrack *candidateFirstDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(0);
  // if (!candidateFirstDaughter) {
  //   std::cout << "candidatePion null" << std::endl;
  //   return 0;
  // }

  // AliAODTrack *candidateSecondDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(1);
  // if (!candidateSecondDaughter) {
  //   std::cout << "candidateKaon null" << std::endl;
  //   return 0;
  // }

  // AliAODTrack *candidateThirdDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(2);
  // if (!candidateSecondDaughter) {
  //   std::cout << "candidatePion null" << std::endl;
  //   return 0;
  // }

  AliAODVertex * vertexB0 = candidateB0->GetSecondaryVtx();
  if (!vertexB0) {
    std::cout << "vertexB0 null" << std::endl;
    return 0;
  }

  // AliAODVertex * vertexDPlus = candidateDPlus->GetSecondaryVtx();
  // if (!vertexDPlus) {
  //   std::cout << "vertexDPlus null" << std::endl;
  //   return 0;
  // }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    std::cout << "primaryVertex null" << std::endl;
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
    Int_t ptbin = PtBin(candidateB0->Pt());
    if (ptbin == -1) return -1;

    // We obtain the variable values in the section below
    // DPlusMass and B0mass
    Double_t mDPlusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mB0PDG = TDatabasePDG::Instance()->GetParticle(511)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mB0PDG - mDPlusPDG;

    // Half width B0 mass
    UInt_t prongs[2];
    prongs[0] = 411; prongs[1] = 211;
    Double_t invMassB0 = candidateB0->InvMass(2, prongs);
    Double_t invMassDifference = TMath::Abs(mB0PDG - invMassB0);
    Double_t invMassDelta = TMath::Abs(deltaPDG - (DeltaInvMassB0Kpipipi(candidateB0)));

    Double_t pointingAngle = candidateB0->CosPointingAngle();
    Double_t dcaMother = candidateB0->GetDCA();
    Double_t ptMother = candidateB0->Pt();
    Double_t momentumMother = candidateB0->P();
    Double_t ptDPlus = candidateDPlus->Pt();
    Double_t ptPion = candidateB0Pion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateB0);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateB0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateB0->Getd0Prong(1));

    Double_t impactProduct = candidateB0->Getd0Prong(0) * candidateB0->Getd0Prong(1);
    Double_t impactProductXY = TMath::Abs(candidateB0->ImpParXY());

    Double_t angleMotherFirstDaughter = (candidateB0->Px() * candidateDPlus->Px() + candidateB0->Py() * candidateDPlus->Py() + candidateB0->Pz() * candidateDPlus->Pz()) / (candidateB0->P() * candidateDPlus->P());
    Double_t angleMotherSecondDaughter = (candidateB0->Px() * candidateB0Pion->Px() + candidateB0->Py() * candidateB0Pion->Py() + candidateB0->Pz() * candidateB0Pion->Pz()) / (candidateB0->P() * candidateB0Pion->P());

    Double_t angleBetweenBothDaughters  = (candidateDPlus->Px() * candidateB0Pion->Px() + candidateDPlus->Py() * candidateB0Pion->Py() + candidateDPlus->Pz() * candidateB0Pion->Pz()) / (candidateDPlus->P() * candidateB0Pion->P());

    Double_t cosThetaStar = candidateB0->CosThetaStar(0, 511, 411, 211);

    Double_t vertexDistance = vertexB0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateB0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();
    Double_t pseudoProperDecayLength = ((vertexB0->GetX() - primaryVertex->GetX()) * candidateB0->Px() / TMath::Abs(candidateB0->Pt())) + ((vertexB0->GetY() - primaryVertex->GetY()) * candidateB0->Py() / TMath::Abs(candidateB0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateB0->Phi();
    Double_t theta = candidateB0->Theta();
    Double_t covMatrix[21];
    candidateB0->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normalizedDecayTime = candidateB0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateB0->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexB0->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateB0->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexB0->GetChi2perNDF();

    Double_t Normalizedd0B0pion = TMath::Abs(candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));
    Double_t Normalizedd0B0 = TMath::Abs(d0[0] / TMath::Sqrt(covd0z0[0]));
    Double_t NormalizedimpactproductB0 = (candidateB0->Getd0Prong(0) / candidateB0->Getd0errProng(0)) * (candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));

    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 50;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 51;
    cutVariableValue = invMassDelta;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 52;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 53;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca12 (cm)" ---------------------------------------------------------
    nCutIndex = 54;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt B0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 55;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt DPlus [GeV/c]" -------------------------------------------------
    nCutIndex = 56;
    cutVariableValue = ptDPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 57;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 B0 (cm)" -------------------------------------------------------
    nCutIndex = 58;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 DPlus (cm)"-----------------------------------------------------
    nCutIndex = 59;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 First Pion (cm)" -----------------------------------------------------
    nCutIndex = 60;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 61;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 62;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angleMotherFirstDaughter" -------------------------------------
    nCutIndex = 63;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angleMotherSecondDaughter" ---------------------------------
    nCutIndex = 64;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angleBetweenBothDaughters" --------------------------------
    nCutIndex = 65;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "cosThetaStar" --------------------------------
    nCutIndex = 66;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 67;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 68;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 69;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 70;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 71;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 72;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 73;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 74;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);

    // "Normalizedd0B0Pion" ----------------------------------------------------
    nCutIndex = 75;
    cutVariableValue = Normalizedd0B0pion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0B0" ----------------------------------------------------
    nCutIndex = 76;
    cutVariableValue = Normalizedd0B0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "NormalizedimpactproductB0" ----------------------------------------------------
    nCutIndex = 77;
    cutVariableValue = NormalizedimpactproductB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------


    // select DPlus
    bPassedCut = IsDPlusFromB0Selected(ptMother, candidateB0, selectionLevel, aod, bCutArray);
  }

  if (bPassedCut == kFALSE)
  {
    returnvalue = 0;
  } else
  {
    for (Int_t i = 50; i < 78; ++i) 
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
Int_t AliRDHFCutsB0toDPi::IsDPlusFromB0Selected(Double_t ptB0, TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[78]) {
  //
  // Apply selection on DPlus candidate from B0 candidate. We have to pass the B0 candidate to this function to get variables w.r.t. B0 vertex.
  //

  if (!fCutsRD) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateB0 = (AliAODRecoDecayHF2Prong*)obj;
  if (!candidateB0) {
    cout << "candidateB0 null" << endl;
    return 0;
  }

  AliAODRecoDecayHF3Prong* candidateDPlus = (AliAODRecoDecayHF3Prong*)candidateB0->GetDaughter(0);
  if (!candidateDPlus) {
    std::cout << "candidateDPlus null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateB0Pion = (AliAODTrack*)candidateB0->GetDaughter(1);
  if (!candidateB0Pion) {
    std::cout << "candidateB0Pion null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateFirstDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(0);
  if (!candidateFirstDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateSecondDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(1);
  if (!candidateSecondDaughter) {
    std::cout << "candidateKaon null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateThirdDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(2);
  if (!candidateSecondDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return 0;
  }

  AliAODVertex * vertexB0 = candidateB0->GetSecondaryVtx();
  if (!vertexB0) {
    std::cout << "vertexB0 null" << std::endl;
    return 0;
  }

  AliAODVertex * vertexDPlus = candidateDPlus->GetSecondaryVtx();
  if (!vertexDPlus) {
    std::cout << "vertexDPlus null" << std::endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    std::cout << "primaryVertex null" << std::endl;
    return 0;
  }

  Int_t returnvalue = 1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();


  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBin(ptB0);

    // DPlusmass
    Double_t mDPlusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mPionPDG = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    // DPlus window - invariant mass
    UInt_t prongs[3];
    prongs[0] = 211; 
    prongs[1] = 321;
    prongs[2] = 211;

    // delta mass PDG
    Double_t deltaPDG = mDPlusPDG - 2 * mPionPDG;

    Double_t invMassDPlus = candidateDPlus->InvMass(3, prongs);
    Double_t invMassDifference = TMath::Abs(mDPlusPDG - invMassDPlus);
    Double_t eKaon1 = candidateDPlus->EProng(1, 321);
    Double_t ePion2 = candidateDPlus->EProng(2, 211);
    Double_t eSum = eKaon1 + ePion2;
    Double_t pxPions = (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2)) * (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2));
    Double_t pyPions = (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2)) * (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2));
    Double_t pzPions = (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2)) * (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2));
    Double_t pSum = pxPions + pyPions + pzPions;
    Double_t invMassPions = TMath::Sqrt(eSum * eSum - pSum);
    Double_t invMassDelta = invMassDPlus - invMassPions;
    Double_t invMassDeltaDifference = TMath::Abs(deltaPDG - invMassDelta);

    Double_t pointingAngle = candidateDPlus->CosPointingAngle();
    Double_t dcaMother12 = candidateDPlus->GetDCA(0);
    Double_t dcaMother13 = candidateDPlus->GetDCA(1);
    Double_t dcaMother23 = candidateDPlus->GetDCA(2);
    Double_t ptMother = candidateDPlus->Pt();
    Double_t momentumMother = candidateDPlus->P();
    Double_t ptKaon = candidateSecondDaughter->Pt();
    Double_t ptFirstPion = candidateFirstDaughter->Pt();
    Double_t ptSecondPion = candidateThirdDaughter->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateDPlus);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateDPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateDPlus->Getd0Prong(1));
    Double_t d0thirdTrack = TMath::Abs(candidateDPlus->Getd0Prong(2));
    
    Double_t impactProduct123 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct12 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1);
    Double_t impactProduct13 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct23 = candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProductXY = TMath::Abs(candidateDPlus->ImpParXY());

    Double_t angleMotherFirstDaughter = (candidateDPlus->Px() * candidateFirstDaughter->Px() + candidateDPlus->Py() * candidateFirstDaughter->Py() + candidateDPlus->Pz() * candidateFirstDaughter->Pz()) / (candidateDPlus->P() * candidateFirstDaughter->P());
    Double_t angleMotherSecondDaughter = (candidateDPlus->Px() * candidateSecondDaughter->Px() + candidateDPlus->Py() * candidateSecondDaughter->Py() + candidateDPlus->Pz() * candidateSecondDaughter->Pz()) / (candidateDPlus->P() * candidateSecondDaughter->P());
    Double_t angleMotherThirdDaughter = (candidateDPlus->Px() * candidateThirdDaughter->Px() + candidateDPlus->Py() * candidateThirdDaughter->Py() + candidateDPlus->Pz() * candidateThirdDaughter->Pz()) / (candidateDPlus->P() * candidateThirdDaughter->P());

    Double_t smallestAngleMotherDaughter = 0.0;
    Double_t largestAngleMotherDaughter = 0.0;

    if(angleMotherFirstDaughter > angleMotherSecondDaughter) {smallestAngleMotherDaughter = angleMotherSecondDaughter;} else {smallestAngleMotherDaughter = angleMotherFirstDaughter;}
    if(angleMotherThirdDaughter < smallestAngleMotherDaughter) smallestAngleMotherDaughter = angleMotherThirdDaughter;
    if(angleMotherFirstDaughter > angleMotherSecondDaughter) {largestAngleMotherDaughter = angleMotherFirstDaughter;} else {largestAngleMotherDaughter = angleMotherSecondDaughter;}
    if(angleMotherThirdDaughter > largestAngleMotherDaughter) largestAngleMotherDaughter = angleMotherThirdDaughter;

    Double_t vertexDistance = vertexDPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateDPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t pseudoProperDecayLength = ((vertexDPlus->GetX() - primaryVertex->GetX()) * candidateDPlus->Px() / TMath::Abs(candidateDPlus->Pt())) + ((vertexDPlus->GetY() - primaryVertex->GetY()) * candidateDPlus->Py() / TMath::Abs(candidateDPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateDPlus->Phi();
    Double_t theta = candidateDPlus->Theta();
    Double_t covMatrix[21];
    candidateDPlus->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normalizedDecayTime = candidateDPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateDPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexDPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateDPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexDPlus->GetChi2perNDF();


    Double_t Normalizedd0DPlusfirstdaughter = TMath::Abs(candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0));
    Double_t Normalizedd0DPlusseconddaughter = TMath::Abs(candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1));
    Double_t Normalizedd0DPlusthirddaughter = TMath::Abs(candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));
    Double_t Normalizedd0DPlus = TMath::Abs(candidateB0->Getd0Prong(1) / candidateB0->Getd0errProng(1));
    Double_t NormalizedimpactproductDPlus = (candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0)) * (candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1))* (candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));

    Double_t Dist12 = candidateDPlus->GetDist12toPrim();
    Double_t Dist23 = candidateDPlus->GetDist23toPrim();
    Double_t SigmaVertex = candidateDPlus->GetSigmaVert();


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 1;
    cutVariableValue = invMassDeltaDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca12 [cm]" ---------------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = dcaMother12;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca13 [cm]" ---------------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = dcaMother13;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca23 [cm]" ---------------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = dcaMother23;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt DPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt First Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = ptFirstPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Second Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = ptSecondPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 DPlus [cm]" -------------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 12;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 First Pion [cm]" -----------------------------------------------------
    nCutIndex = 13;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Second Pion [cm]" -----------------------------------------------------
    nCutIndex = 14;
    cutVariableValue = d0thirdTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0d0 [cm^3]" ------------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = impactProduct123;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = impactProduct12;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 17;
    cutVariableValue = impactProduct13;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = impactProduct23;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "smallestAngleMotherDaughter" -------------------------------------
    nCutIndex = 20;
    cutVariableValue = smallestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "largestAngleMotherDaughter" ---------------------------------
    nCutIndex = 21;
    cutVariableValue = largestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle difference" --------------------------------
    nCutIndex = 22;
    cutVariableValue = largestAngleMotherDaughter - smallestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 25;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 29;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 30;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusfirstdaughter" ----------------------------------------------------
    nCutIndex = 31;
    cutVariableValue = Normalizedd0DPlusfirstdaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusseconddaughter" ----------------------------------------------------
    nCutIndex = 32;
    cutVariableValue = Normalizedd0DPlusseconddaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusthirddaughter" ----------------------------------------------------
    nCutIndex = 33;
    cutVariableValue = Normalizedd0DPlusthirddaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlus" ----------------------------------------------------
    nCutIndex = 34;
    cutVariableValue = Normalizedd0DPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "NormalizedimpactproductDPlus" ----------------------------------------------------
    nCutIndex = 35;
    cutVariableValue = NormalizedimpactproductDPlus;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Dist12" ----------------------------------------------------
    nCutIndex = 36;
    cutVariableValue = Dist12;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Dist23" ----------------------------------------------------
    nCutIndex = 37;
    cutVariableValue = Dist23;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Sigma vertex" ----------------------------------------------------
    nCutIndex = 38;
    cutVariableValue = SigmaVertex;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------


    AliAODRecoDecay* candidateDPlustoB0 = (AliAODRecoDecay*)candidateDPlus;
    AliExternalTrackParam firstDaughterDPlusTrack;
    AliExternalTrackParam secondDaughterDPlusTrack;
    AliExternalTrackParam thirdDaughterDPlusTrack;

    Double_t d0z0DSVert[2], covd0z0DSVert[3], d0DSVert[3];

    firstDaughterDPlusTrack.CopyFromVTrack(candidateFirstDaughter);
    firstDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    secondDaughterDPlusTrack.CopyFromVTrack(candidateSecondDaughter);
    secondDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    thirdDaughterDPlusTrack.CopyFromVTrack(candidateThirdDaughter);
    thirdDaughterDPlusTrack.PropagateToDCA(vertexB0, bz, 100., d0z0DSVert, covd0z0DSVert);
    d0DSVert[2] = d0z0DSVert[0];

    AliExternalTrackParam DPlusTrack;
    DPlusTrack.CopyFromVTrack(candidateDPlus);
    Double_t d0z0D0DSVert[2], covd0z0D0DSVert[3], d0d0DSVert;
    motherTrack.PropagateToDCA(vertexB0, bz, 100., d0z0D0DSVert, covd0z0D0DSVert);
    d0d0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToB0 = d0DSVert[0] * d0DSVert[1] * d0DSVert[2];
    Double_t impactProductXYToB0 = candidateDPlustoB0->ImpParXY(vertexB0);

    Double_t pointingAngleToB0 = candidateDPlustoB0->CosPointingAngle(vertexB0);
    Double_t d0FirstDaughterToB0 = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToB0 = TMath::Abs(d0DSVert[1]);
    Double_t d0ThirdDaughterToB0 = TMath::Abs(d0DSVert[2]);
    Double_t normDecayLengthToB0 = candidateDPlustoB0->NormalizedDecayLength(vertexB0);

    Double_t pseudoProperDecayLengthDSVert = ((vertexDPlus->GetX() - vertexB0->GetX()) * candidateDPlus->Px() / TMath::Abs(candidateDPlus->Pt())) + ((vertexDPlus->GetY() - vertexB0->GetY()) * candidateDPlus->Py() / TMath::Abs(candidateDPlus->Pt()));
    Double_t pseudoProperDecayTimeToB0 = pseudoProperDecayLengthDSVert * pdgMassMother / ptMother;
    Double_t DecayTimeToB0 = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phiDSVert = candidateDPlus->Phi();
    Double_t thetaDSVert = candidateDPlus->Theta();
    Double_t covMatrixDSVert[21];
    candidateDPlus->GetCovarianceXYZPxPyPz(covMatrixDSVert);

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
    Double_t normalizedDecayTimeToB0 = candidateDPlustoB0->NormalizedDecayLength(vertexB0) / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    // "pointingAngleToB0" ---------------------------------------------
    nCutIndex = 39;
    cutVariableValue = pointingAngleToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0MotherToB0" --------------------------------------------------
    nCutIndex = 40;
    cutVariableValue = d0d0DSVert;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0FirstDaughterToB0" -------------------------------------------
    nCutIndex = 41;
    cutVariableValue = d0FirstDaughterToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0SecondDaughterToB0" ------------------------------------------
    nCutIndex = 42;
    cutVariableValue = d0SecondDaughterToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0ThirdDaughterToB0" ------------------------------------------
    nCutIndex = 43;
    cutVariableValue = d0ThirdDaughterToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "impactProductToB0" ---------------------------------------------
    nCutIndex = 44;
    cutVariableValue = impactProductToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "impactProductXYToB0" -------------------------------------------
    nCutIndex = 45;
    cutVariableValue = impactProductXYToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLengthToB0" -------------------------------------------
    nCutIndex = 46;
    cutVariableValue = normDecayLengthToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToB0" -------------------------------------
    nCutIndex = 47;
    cutVariableValue = pseudoProperDecayTimeToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTimeToB0" -------------------------------------------------
    nCutIndex = 48;
    cutVariableValue = DecayTimeToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToB0" ---------------------------------------------
    nCutIndex = 49;
    cutVariableValue = normalizedDecayTimeToB0;
    bPassedCut = ApplyCutOnVariable(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------
  }

  if (bPassedCut == kFALSE)
  {
    returnvalue = 0;
  } else
  for (Int_t i = 0; i < 50; ++i)
  {
    if (bCutArray[i] == kTRUE) {
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDPi::IsDPlusforDPlusptbinSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[39]) {
//
  // Apply selection on DPlus candidate.
  //

  if (!fCutsRDDPlusforDPlusptbin) {
    cout << "Cut matrice not inizialized. Exit..." << endl;
    return 0;
  }

  AliAODRecoDecayHF3Prong* candidateDPlus = (AliAODRecoDecayHF3Prong*)obj;
  if (!candidateDPlus) {
    cout << "candidateDPlus null" << endl;
    return 0;
  }

  AliAODTrack *candidateFirstDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(0);
  if (!candidateFirstDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateSecondDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(1);
  if (!candidateSecondDaughter) {
    std::cout << "candidateKaon null" << std::endl;
    return 0;
  }

  AliAODTrack *candidateThirdDaughter = (AliAODTrack*)candidateDPlus->GetDaughter(2);
  if (!candidateSecondDaughter) {
    std::cout << "candidatePion null" << std::endl;
    return 0;
  }

  AliAODVertex * vertexDPlus = candidateDPlus->GetSecondaryVtx();
  if (!vertexDPlus) {
    std::cout << "vertexDPlus null" << std::endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if (!primaryVertex) {
    std::cout << "primaryVertex null" << std::endl;
    return 0;
  }
  Int_t returnvalue = 1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField();


  // selection on candidate
  if (selectionLevel == AliRDHFCuts::kAll ||
      selectionLevel == AliRDHFCuts::kCandidate) {

    Int_t ptbin = PtBinDPlusforDPlusptbin(candidateDPlus->Pt());
    if (ptbin == -1) return -1;

    // DPlusmass
    Double_t mDPlusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mPionPDG = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    // DPlus window - invariant mass
    UInt_t prongs[3];
    prongs[0] = 211; 
    prongs[1] = 321;
    prongs[2] = 211;

    // delta mass PDG
    Double_t deltaPDG = mDPlusPDG - 2 * mPionPDG;

    Double_t invMassDPlus = candidateDPlus->InvMass(3, prongs);
    Double_t invMassDifference = TMath::Abs(mDPlusPDG - invMassDPlus);
    Double_t eKaon1 = candidateDPlus->EProng(1, 321);
    Double_t ePion2 = candidateDPlus->EProng(2, 211);
    Double_t eSum = eKaon1 + ePion2;
    Double_t pxPions = (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2)) * (candidateDPlus->PxProng(1) + candidateDPlus->PxProng(2));
    Double_t pyPions = (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2)) * (candidateDPlus->PyProng(1) + candidateDPlus->PyProng(2));
    Double_t pzPions = (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2)) * (candidateDPlus->PzProng(1) + candidateDPlus->PzProng(2));
    Double_t pSum = pxPions + pyPions + pzPions;
    Double_t invMassPions = TMath::Sqrt(eSum * eSum - pSum);
    Double_t invMassDelta = invMassDPlus - invMassPions;
    Double_t invMassDeltaDifference = TMath::Abs(deltaPDG - invMassDelta);

    Double_t pointingAngle = candidateDPlus->CosPointingAngle();
    Double_t dcaMother12 = candidateDPlus->GetDCA(0);
    Double_t dcaMother13 = candidateDPlus->GetDCA(1);
    Double_t dcaMother23 = candidateDPlus->GetDCA(2);
    Double_t ptMother = candidateDPlus->Pt();
    Double_t momentumMother = candidateDPlus->P();
    Double_t ptKaon = candidateSecondDaughter->Pt();
    Double_t ptFirstPion = candidateFirstDaughter->Pt();
    Double_t ptSecondPion = candidateThirdDaughter->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateDPlus);
    Double_t d0z0[2], covd0z0[3], d0[2];
    motherTrack.PropagateToDCA(primaryVertex, bz, 100., d0z0, covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateDPlus->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateDPlus->Getd0Prong(1));
    Double_t d0thirdTrack = TMath::Abs(candidateDPlus->Getd0Prong(2));
    
    Double_t impactProduct123 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct12 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(1);
    Double_t impactProduct13 = candidateDPlus->Getd0Prong(0) * candidateDPlus->Getd0Prong(2);
    Double_t impactProduct23 = candidateDPlus->Getd0Prong(1) * candidateDPlus->Getd0Prong(2);
    Double_t impactProductXY = TMath::Abs(candidateDPlus->ImpParXY());

    Double_t angleMotherFirstDaughter = (candidateDPlus->Px() * candidateFirstDaughter->Px() + candidateDPlus->Py() * candidateFirstDaughter->Py() + candidateDPlus->Pz() * candidateFirstDaughter->Pz()) / (candidateDPlus->P() * candidateFirstDaughter->P());
    Double_t angleMotherSecondDaughter = (candidateDPlus->Px() * candidateSecondDaughter->Px() + candidateDPlus->Py() * candidateSecondDaughter->Py() + candidateDPlus->Pz() * candidateSecondDaughter->Pz()) / (candidateDPlus->P() * candidateSecondDaughter->P());
    Double_t angleMotherThirdDaughter = (candidateDPlus->Px() * candidateThirdDaughter->Px() + candidateDPlus->Py() * candidateThirdDaughter->Py() + candidateDPlus->Pz() * candidateThirdDaughter->Pz()) / (candidateDPlus->P() * candidateThirdDaughter->P());

    Double_t smallestAngleMotherDaughter = 0.0;
    Double_t largestAngleMotherDaughter = 0.0;

    if(angleMotherFirstDaughter > angleMotherSecondDaughter) {smallestAngleMotherDaughter = angleMotherSecondDaughter;} else {smallestAngleMotherDaughter = angleMotherFirstDaughter;}
    if(angleMotherThirdDaughter < smallestAngleMotherDaughter) smallestAngleMotherDaughter = angleMotherThirdDaughter;
    if(angleMotherFirstDaughter > angleMotherSecondDaughter) {largestAngleMotherDaughter = angleMotherFirstDaughter;} else {largestAngleMotherDaughter = angleMotherSecondDaughter;}
    if(angleMotherThirdDaughter > largestAngleMotherDaughter) largestAngleMotherDaughter = angleMotherThirdDaughter;

    Double_t vertexDistance = vertexDPlus->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateDPlus->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t pseudoProperDecayLength = ((vertexDPlus->GetX() - primaryVertex->GetX()) * candidateDPlus->Px() / TMath::Abs(candidateDPlus->Pt())) + ((vertexDPlus->GetY() - primaryVertex->GetY()) * candidateDPlus->Py() / TMath::Abs(candidateDPlus->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother / ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother / (momentumMother * momentumMother)) + 1)));

    Double_t phi = candidateDPlus->Phi();
    Double_t theta = candidateDPlus->Theta();
    Double_t covMatrix[21];
    candidateDPlus->GetCovarianceXYZPxPyPz(covMatrix);

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
    Double_t normalizedDecayTime = candidateDPlus->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1 / ((pdgMassMother * pdgMassMother * errorMomentum * errorMomentum / (momentumMother * momentumMother)) + 1)));

    Double_t cosPointingAngleXY = candidateDPlus->CosPointingAngleXY();
    Double_t distanceXYToVertex = vertexDPlus->DistanceXYToVertex(primaryVertex);
    Double_t normalizedDecayLengthXY = candidateDPlus->NormalizedDecayLengthXY();
    Double_t chi2Vertex = vertexDPlus->GetChi2perNDF();


    Double_t Normalizedd0DPlusfirstdaughter = TMath::Abs(candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0));
    Double_t Normalizedd0DPlusseconddaughter = TMath::Abs(candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1));
    Double_t Normalizedd0DPlusthirddaughter = TMath::Abs(candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));
    Double_t Normalizedd0DPlus = TMath::Abs(candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1));
    Double_t NormalizedimpactproductDPlus = (candidateDPlus->Getd0Prong(0) / candidateDPlus->Getd0errProng(0)) * (candidateDPlus->Getd0Prong(1) / candidateDPlus->Getd0errProng(1))* (candidateDPlus->Getd0Prong(2) / candidateDPlus->Getd0errProng(2));

    Double_t Dist12 = candidateDPlus->GetDist12toPrim();
    Double_t Dist23 = candidateDPlus->GetDist23toPrim();
    Double_t SigmaVertex = candidateDPlus->GetSigmaVert();


    // We apply the cuts
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 1;
    cutVariableValue = invMassDeltaDifference;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pointing angle XY" ----------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = cosPointingAngleXY;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca12 [cm]" ---------------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = dcaMother12;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca13 [cm]" ---------------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = dcaMother13;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "dca23 [cm]" ---------------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = dcaMother23;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt DPlus [GeV/c]" ----------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt First Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = ptFirstPion;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Pt Second Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = ptSecondPion;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 DPlus [cm]" -------------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 12;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 First Pion [cm]" -----------------------------------------------------
    nCutIndex = 13;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 Second Pion [cm]" -----------------------------------------------------
    nCutIndex = 14;
    cutVariableValue = d0thirdTrack;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0d0 [cm^3]" ------------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = impactProduct123;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = impactProduct12;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 17;
    cutVariableValue = impactProduct13;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = impactProduct23;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "smallestAngleMotherDaughter" -------------------------------------
    nCutIndex = 20;
    cutVariableValue = smallestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "largestAngleMotherDaughter" ---------------------------------
    nCutIndex = 21;
    cutVariableValue = largestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "angle difference" --------------------------------
    nCutIndex = 22;
    cutVariableValue = largestAngleMotherDaughter - smallestAngleMotherDaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "vertex distance XY" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = distanceXYToVertex;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 25;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 27;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 28;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "normalized decay length XY" ----------------------------------------------------
    nCutIndex = 29;
    cutVariableValue = normalizedDecayLengthXY;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "chi squared per NDF" ----------------------------------------------------
    nCutIndex = 30;
    cutVariableValue = chi2Vertex;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusfirstdaughter" ----------------------------------------------------
    nCutIndex = 31;
    cutVariableValue = Normalizedd0DPlusfirstdaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusseconddaughter" ----------------------------------------------------
    nCutIndex = 32;
    cutVariableValue = Normalizedd0DPlusseconddaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlusthirddaughter" ----------------------------------------------------
    nCutIndex = 33;
    cutVariableValue = Normalizedd0DPlusthirddaughter;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Normalizedd0DPlus" ----------------------------------------------------
    nCutIndex = 34;
    cutVariableValue = Normalizedd0DPlus;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "NormalizedimpactproductDPlus" ----------------------------------------------------
    nCutIndex = 35;
    cutVariableValue = NormalizedimpactproductDPlus;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Dist12" ----------------------------------------------------
    nCutIndex = 36;
    cutVariableValue = Dist12;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Dist23" ----------------------------------------------------
    nCutIndex = 37;
    cutVariableValue = Dist23;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
    //---------------------------------------------------------------------

    // "Sigma vertex" ----------------------------------------------------
    nCutIndex = 38;
    cutVariableValue = SigmaVertex;
    bPassedCut = ApplyCutOnVariableDPlusforDPlusptbin(nCutIndex, ptbin, cutVariableValue, bCutArray);
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
//--------------------------------------------------------------------------
Bool_t AliRDHFCutsB0toDPi::IsInFiducialAcceptance(Double_t /*pt*/, Double_t y) const
{
  //
  // Chould be improved (for example pT dependence)
  //

  if (y > GetFiducialYCut()) return kFALSE;

  return kTRUE;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsB0toDPi::IsSelectedPID(AliAODRecoDecayHF* /*obj*/)
{
  //
  // PID method, n sigma approach default // not used for B0, done seperately for each daughter
  //

  // AliAODRecoDecayHF2Prong* dstar = (AliAODRecoDecayHF2Prong*)obj;
  // if(!dstar){
  //   cout<<"AliAODRecoDecayHF2Prong null"<<endl;
  //   return 0;
  // }

  // if(!fUsePID || dstar->Pt() > fMaxPtPid) return 3;

  // AliAODRecoDecayHF2Prong* dPlus = (AliAODRecoDecayHF2Prong*)dstar->Get2Prong();
  // if(!dPlus){
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
Int_t AliRDHFCutsB0toDPi::SelectPID(AliAODTrack *track, Int_t type)
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
Double_t AliRDHFCutsB0toDPi::DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const
{
  ///
  /// delta invariant mass
  ///

  AliAODRecoDecayHF3Prong * DPlus = (AliAODRecoDecayHF3Prong*)Bzero->GetDaughter(0);

  Double_t e[4] = {0};
  e[0] = DPlus->EProng(0, 211); // to be done: Check effect of taking energies of DPlus daughters at B vertex by propagating D to B first ////////////////////////////////////////////////////////////////////////////////////////////////////////
  e[1] = DPlus->EProng(1, 321);
  e[2] = DPlus->EProng(2, 211);
  e[3] = Bzero->EProng(1, 211);

  Double_t esum = e[0] + e[1] + e[2] + e[3];
  Double_t invMassB0 = TMath::Sqrt(esum * esum - Bzero->P2());

  Double_t eD[3] = {0};
  eD[0] = DPlus->EProng(0, 211);
  eD[1] = DPlus->EProng(1, 321);
  eD[2] = DPlus->EProng(2, 211);

  Double_t esumD = eD[0] + eD[1] + eD[2];
  Double_t invMassDPlus = TMath::Sqrt(esumD * esumD - DPlus->P2());

  return invMassB0 - invMassDPlus;
}


//---------------------------------------------------------------------------
//
//  DO for DPlus pt bin functions
//
//---------------------------------------------------------------------------

void AliRDHFCutsB0toDPi::SetCutsDPlusforDPlusptbin(Int_t nVars, Int_t nPtBins, Float_t **cutsRDDPlusforDPlusptbin) {
  //
  // store the cuts
  //

  if (nVars != fnVarsDPlusforDPlusptbin) {
    printf("Wrong number of variables: it has to be %d\n", fnVarsDPlusforDPlusptbin);
    AliFatal("exiting");
  }
  if (nPtBins != fnPtBinsDPlusforDPlusptbin) {
    printf("Wrong number of pt bins: it has to be %d\n", fnPtBinsDPlusforDPlusptbin);
    AliFatal("exiting");
  }

  if (!fCutsRDDPlusforDPlusptbin)  fCutsRDDPlusforDPlusptbin = new Float_t[fGlobalIndexDPlusforDPlusptbin];


  for (Int_t iv = 0; iv < fnVarsDPlusforDPlusptbin; iv++)
  {
    for (Int_t ib = 0; ib < fnPtBinsDPlusforDPlusptbin; ib++)
    {
      //check

      if (GetGlobalIndexDPlusforDPlusptbin(iv, ib) >= fGlobalIndexDPlusforDPlusptbin)
      {
        cout << "Overflow, exit..." << endl;
        return;
      }

      fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(iv, ib)] = cutsRDDPlusforDPlusptbin[iv][ib];

    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::SetCutsDPlusforDPlusptbin(Int_t glIndex, Float_t *cutsRDDPlusforDPlusptbin) {
  //
  // store the cuts
  //
  if (glIndex != fGlobalIndexDPlusforDPlusptbin) {
    cout << "Wrong array size: it has to be " << fGlobalIndexDPlusforDPlusptbin << endl;
    AliFatal("exiting");
  }
  if (!fCutsRDDPlusforDPlusptbin)  fCutsRDDPlusforDPlusptbin = new Float_t[fGlobalIndexDPlusforDPlusptbin];

  for (Int_t iGl = 0; iGl < fGlobalIndexDPlusforDPlusptbin; iGl++) {
    fCutsRDDPlusforDPlusptbin[iGl] = cutsRDDPlusforDPlusptbin[iGl];
  }
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDPi::PtBinDPlusforDPlusptbin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin = -1;
  if (pt < fPtBinLimitsDPlusforDPlusptbin[0])return ptbin;
  for (Int_t i = 0; i < fnPtBinsDPlusforDPlusptbin; i++) {
    if (pt < fPtBinLimitsDPlusforDPlusptbin[i + 1]) {
      ptbin = i;
      break;
    }
  }
  return ptbin;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::SetPtBinsDPlusforDPlusptbin(Int_t nPtBinLimits, Float_t *ptBinLimits) {
  // Set the pt bins

  if (fPtBinLimitsDPlusforDPlusptbin) {
    delete [] fPtBinLimitsDPlusforDPlusptbin;
    fPtBinLimitsDPlusforDPlusptbin = NULL;
    printf("Changing the pt bins\n");
  }

  if (nPtBinLimits != fnPtBinsDPlusforDPlusptbin + 1) {
    cout << "Warning: ptBinLimits dimension " << nPtBinLimits << " != nPtBins+1 (" << fnPtBinsDPlusforDPlusptbin + 1 << ")\nSetting nPtBins to " << nPtBinLimits - 1 << endl;
    SetNPtBinsDPlusforDPlusptbin(nPtBinLimits - 1);
  }

  fnPtBinLimitsDPlusforDPlusptbin = nPtBinLimits;
  SetGlobalIndexDPlusforDPlusptbin();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimitsDPlusforDPlusptbin = new Float_t[fnPtBinLimitsDPlusforDPlusptbin];
  for (Int_t ib = 0; ib < nPtBinLimits; ib++) fPtBinLimitsDPlusforDPlusptbin[ib] = ptBinLimits[ib];

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDPi::ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[78]) {

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
Int_t AliRDHFCutsB0toDPi::ApplyCutOnVariableDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[39]) {

  // std::cout << "index: " << nCutIndex << ", ptbin: " << ptbin << ", used: " << GetIsCutUsedDPlusforDPlusptbin(nCutIndex,ptbin) << ", upper: " << GetIsUpperCutDPlusforDPlusptbin(nCutIndex) << ", value: " << cutVariableValue << ", cutvalue: " << fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex,ptbin)] << std::endl;
  if (GetIsCutUsedDPlusforDPlusptbin(nCutIndex, ptbin) == kTRUE)
  {
    Bool_t bCut = kFALSE;
    if (GetIsUpperCutDPlusforDPlusptbin(nCutIndex) == kTRUE)
    {
      if (cutVariableValue > fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex, ptbin)]) bCut = kTRUE;
    } else
    {
      if (cutVariableValue < fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex, ptbin)]) bCut = kTRUE;
    }
    if (bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::SetVarNamesDPlusforDPlusptbin(Int_t nVars, TString *varNames, Bool_t *isUpperCut) {
  // Set the variable names

  if (fVarNamesDPlusforDPlusptbin) {
    delete [] fVarNamesDPlusforDPlusptbin;
    fVarNamesDPlusforDPlusptbin = NULL;
    //printf("Changing the variable names\n");
  }
  if (nVars != fnVarsDPlusforDPlusptbin) {
    printf("Wrong number of variables: it has to be %d\n", fnVarsDPlusforDPlusptbin);
    return;
  }
  //fnVars=nVars;
  fVarNamesDPlusforDPlusptbin = new TString[nVars];
  fIsUpperCutDPlusforDPlusptbin = new Bool_t[nVars];
  for (Int_t iv = 0; iv < nVars; iv++) {
    fVarNamesDPlusforDPlusptbin[iv] = varNames[iv];
    fIsUpperCutDPlusforDPlusptbin[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::InitializeCuts() {
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

  // DPlus for DPlus pt bin

  if (fIsCutUsedDPlusforDPlusptbin) {
    delete [] fIsCutUsedDPlusforDPlusptbin;
    fIsCutUsedDPlusforDPlusptbin = NULL;
  }

  fIsCutUsedDPlusforDPlusptbin = new Bool_t[(GetNPtBinsDPlusforDPlusptbin()) * (GetNVarsDPlusforDPlusptbin())];
  for (Int_t iv = 0; iv < (GetNPtBinsDPlusforDPlusptbin()) * (GetNVarsDPlusforDPlusptbin()); iv++) {
    fIsCutUsedDPlusforDPlusptbin[iv] = kFALSE;

  }

  if (!fCutsRDDPlusforDPlusptbin)  fCutsRDDPlusforDPlusptbin = new Float_t[fGlobalIndexDPlusforDPlusptbin];

  for (Int_t iv = 0; iv < fnVarsDPlusforDPlusptbin; iv++)
  {
    for (Int_t ib = 0; ib < fnPtBinsDPlusforDPlusptbin; ib++)
    {
      //check

      if (GetGlobalIndexDPlusforDPlusptbin(iv, ib) >= fGlobalIndexDPlusforDPlusptbin)
      {
        cout << "Overflow, exit..." << endl;
        return;
      }

      fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(iv, ib)] = 0;
    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t cutValue) {
  // Set the cut value and direction

  fIsCutUsed[GetGlobalIndex(nCutIndex, ptBin)] = kTRUE;
  fIsUpperCut[nCutIndex] = cutDirection;
  fCutsRD[GetGlobalIndex(nCutIndex, ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::SetCutDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t cutValue) {
  // Set the cut value and direction

  fIsCutUsedDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex, ptBin)] = kTRUE;
  fIsUpperCutDPlusforDPlusptbin[nCutIndex] = cutDirection;
  fCutsRDDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex, ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDPi::InitializeCutsForCutOptimization(Int_t nCutsForOptimization, Int_t nVariables) {
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
void AliRDHFCutsB0toDPi::SetCutsForCutOptimization(Int_t glIndex, Float_t *cutsRDForCutOptimization) {
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
void AliRDHFCutsB0toDPi::SetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t * cutValues) {
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
