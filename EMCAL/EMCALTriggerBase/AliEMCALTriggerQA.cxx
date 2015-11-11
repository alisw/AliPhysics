/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
/**
 * @file AliEMCALTriggerQA.cxx
 * @date Nov. 1, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */

#include <THashList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliLog.h"

#include "AliEMCALTriggerQA.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerQA)
/// \endcond

//const Int_t AliEMCALTriggerQA::fgkNTriggerTypes = 6;
const TString AliEMCALTriggerQA::fgkTriggerTypeNames[AliEMCALTriggerQA::fgkNTriggerTypes] = {"EJE", "EGA", "EL0", "REJE", "REGA", "EBKG"};

/**
 * Dummy constructor
 */
AliEMCALTriggerQA::AliEMCALTriggerQA():
  TNamed(),
  fBkgOfflineAmp(kFALSE),
  fBkgPatchType(3),
  fDebugLevel(0),
  fHistos(0)
{
  for (int itype = 0; itype < fgkNTriggerTypes; itype++) {
    fADCAmpEMCal[itype].Set(100);
    fNPatchesEMCal[itype] = 0;
    fMaxPatchEMCal[itype] = 0;

    fADCAmpDCal[itype].Set(100);
    fNPatchesDCal[itype] = 0;
    fMaxPatchDCal[itype] = 0;
  }
}

/**
 * Constructor
 */
AliEMCALTriggerQA::AliEMCALTriggerQA(const char* name):
  TNamed(name,name),
  fBkgOfflineAmp(kFALSE),
  fBkgPatchType(3),
  fDebugLevel(0),
  fHistos(0)
{
  for (int itype = 0; itype < fgkNTriggerTypes; itype++) {
    fADCAmpEMCal[itype].Set(100);
    fNPatchesEMCal[itype] = 0;
    fMaxPatchEMCal[itype] = 0;

    fADCAmpDCal[itype].Set(100);
    fNPatchesDCal[itype] = 0;
    fMaxPatchDCal[itype] = 0;
  }
}

/**
 * Destructor
 */
AliEMCALTriggerQA::~AliEMCALTriggerQA()
{
}

/**
 * Calculate the patch size
 * \param patch Pointer to a valid trigger patch
 * \return the patch size in number of FastORs.
 */
Int_t AliEMCALTriggerQA::GetPatchType(AliEMCALTriggerPatchInfo* patch) 
{
  Int_t type = -1;

  if (patch->IsJetLow() || patch->IsJetHigh()) { 
    type = kTMEMCalJet;
  }
  else if (patch->IsGammaLow() || patch->IsGammaHigh()) {
    type = kTMEMCalGamma;
  }
  else if (patch->IsLevel0()) {
    type = kTMEMCalLevel0;
  }
  else if (patch->IsRecalcJet() || patch->GetPatchSize() == 16) {
    type = kTMEMCalRecalcJet;
  }
  else if (patch->IsRecalcGamma() || patch->GetPatchSize() == 2) {
    type = kTMEMCalRecalcGamma;
  }
  else if (patch->GetPatchSize() == 8) {
    type = kTMEMCalBackground;
  }
  
  return type;
}

/**
 * Initialize the class, i.e. allocate histograms.
 */
void AliEMCALTriggerQA::Init()
{
  TString hname;

  const char *patchtypes[2] = {"Online", "Offline"};

  fHistos = new THashList();
  fHistos->SetName(Form("histos%s", GetName()));
  fHistos->SetOwner(kTRUE);

  for (int itype = 0; itype < fgkNTriggerTypes; itype++) {
    for (const char **patchtype = patchtypes; patchtype < patchtypes + 2; ++patchtype) {
      hname = Form("histEMCalPatchAmp%s%s", fgkTriggerTypeNames[itype].Data(), *patchtype);
      CreateTH1(hname, hname, 1000, 0, 1000);

      hname = Form("histDCalPatchAmp%s%s", fgkTriggerTypeNames[itype].Data(), *patchtype);
      CreateTH1(hname, hname, 1000, 0, 1000);
    }

    hname = Form("histEdgePos%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 200, 0, 200, 200, 0, 200);

    hname = Form("histCMPos%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 200, 0, TMath::TwoPi(), 60, -1, 1);

    hname = Form("histGeoPos%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 200, 0, TMath::TwoPi(), 60, -1, 1);

    hname = Form("histEMCalPatchEnergy%s", fgkTriggerTypeNames[itype].Data());
    CreateTH1(hname, hname, 200, 0, 200);

    hname = Form("histDCalPatchEnergy%s", fgkTriggerTypeNames[itype].Data());
    CreateTH1(hname, hname, 200, 0, 200);

    hname = Form("histEMCalMedianVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);

    hname = Form("histDCalMedianVsEMCalMax%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);

    hname = Form("histEMCalMedianVsEMCalMax%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);

    hname = Form("histDCalMedianVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);

    hname = Form("histEMCalMedianVsDCalMedian%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);

    hname = Form("histEMCalMaxVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    CreateTH2(hname, hname, 1000, 0, 1000, 1000, 0, 1000);
  }
}

/**
 * Process a patch, filling relevant histograms.
 * \param patch Pointer to a valid trigger patch
 */
void AliEMCALTriggerQA::ProcessPatch(AliEMCALTriggerPatchInfo* patch)
{
  TString hname;
  
  Int_t type = GetPatchType(patch);
  
  if (type < 0) return;

  hname = Form("histEdgePos%s", fgkTriggerTypeNames[type].Data());
  FillTH2(hname, patch->GetRowStart(), patch->GetColStart());

  hname = Form("histCMPos%s", fgkTriggerTypeNames[type].Data());
  FillTH2(hname, patch->GetPhiCM(), patch->GetEtaCM());

  hname = Form("histGeoPos%s", fgkTriggerTypeNames[type].Data());
  FillTH2(hname, patch->GetPhiGeo(), patch->GetEtaGeo());

  TString det;
  
  Int_t amp = fBkgOfflineAmp ? patch->GetADCOfflineAmp() : patch->GetADCAmp();
  
  if (patch->IsEMCal()) {
    det = "EMCal";
    if (fNPatchesEMCal[type] >= fADCAmpEMCal[type].GetSize()) {
      fADCAmpEMCal[type].Set((fNPatchesEMCal[type]+1)*2);
    }
    fADCAmpEMCal[type].AddAt(amp, fNPatchesEMCal[type]);
    if (fMaxPatchEMCal[type] < amp) fMaxPatchEMCal[type] = amp;
    fNPatchesEMCal[type]++;
  }
  else if (patch->IsDCalPHOS()) {
    det = "DCal";
    if (fNPatchesDCal[type] >= fADCAmpDCal[type].GetSize()) {
      fADCAmpDCal[type].Set((fNPatchesDCal[type]+1)*2);
    }
    fADCAmpDCal[type].AddAt(amp, fNPatchesDCal[type]);
    if (fMaxPatchDCal[type] < amp) fMaxPatchDCal[type] = amp;
    fNPatchesDCal[type]++;
  }
  else {
    AliWarning(Form("Patch is not EMCal nor DCal/PHOS (pos: %d, %d)", patch->GetRowStart(), patch->GetColStart()));
  }

  if (patch->GetADCAmp() > 0) {
    hname = Form("hist%sPatchAmp%sOnline", det.Data(), fgkTriggerTypeNames[type].Data());
    FillTH1(hname, patch->GetADCAmp());
  }

  if (patch->GetADCOfflineAmp() > 0) {
    hname = Form("hist%sPatchAmp%sOffline", det.Data(), fgkTriggerTypeNames[type].Data());
    FillTH1(hname, patch->GetADCOfflineAmp());
  }

  if (patch->GetPatchE() > 0) {
    hname = Form("hist%sPatchEnergy%s", det.Data(), fgkTriggerTypeNames[type].Data());
    FillTH1(hname, patch->GetPatchE());
  }

  if (fDebugLevel >= 2) {
    Printf("Type = %s; global pos = (%d, %d); Amp (online) = %d; Amp (offline) = %d; Patch energy = %.3f\n"
	   "Position (CM): Eta=%.3f, Phi=%.3f\n"
	   "Position (Geo): Eta=%.3f, Phi=%.3f\n",
	   fgkTriggerTypeNames[type].Data(), patch->GetRowStart(), patch->GetColStart(), patch->GetADCAmp(), patch->GetADCOfflineAmp(), patch->GetPatchE(),
	   patch->GetEtaCM(), patch->GetPhiCM(),
	   patch->GetEtaGeo(), patch->GetPhiGeo());
  }
}

/**
 * This method should be called at the end of each event.
 */
void AliEMCALTriggerQA::EventCompleted()
{
  TString hname;

  Double_t medianEMCal = TMath::Median(fNPatchesEMCal[fBkgPatchType], fADCAmpEMCal[fBkgPatchType].GetArray());
  Double_t medianDCal = TMath::Median(fNPatchesDCal[fBkgPatchType], fADCAmpDCal[fBkgPatchType].GetArray());

  for (int itype = 0; itype < fgkNTriggerTypes; itype++) {
    if (fNPatchesEMCal[itype] == 0 && fNPatchesDCal[itype] == 0) continue;

    if (fDebugLevel >= 2) Printf("Type %s: fNPatchesEMCal = %d, fNPatchesDCal = %d", fgkTriggerTypeNames[itype].Data(), fNPatchesEMCal[itype], fNPatchesDCal[itype]);
    
    hname = Form("histEMCalMedianVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, medianEMCal, fMaxPatchDCal[itype]);

    hname = Form("histDCalMedianVsEMCalMax%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, medianDCal, fMaxPatchEMCal[itype]);

    hname = Form("histDCalMedianVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, medianDCal, fMaxPatchDCal[itype]);

    hname = Form("histEMCalMedianVsEMCalMax%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, medianEMCal, fMaxPatchEMCal[itype]);

    hname = Form("histEMCalMaxVsDCalMax%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, fMaxPatchEMCal[itype], fMaxPatchDCal[itype]);

    hname = Form("histEMCalMedianVsDCalMedian%s", fgkTriggerTypeNames[itype].Data());
    FillTH2(hname, medianEMCal, medianDCal);

    fADCAmpEMCal[itype].Reset();
    fNPatchesEMCal[itype] = 0;
    fMaxPatchEMCal[itype] = 0;

    fADCAmpDCal[itype].Reset();
    fNPatchesDCal[itype] = 0;
    fMaxPatchDCal[itype] = 0;
  }
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax)
{
  /*
   * Create a new TH1 within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param title: Title of the histogram
   * @param nbins: number of bins
   * @param xmin: min. value of the range
   * @param xmax: max. value of the range
   * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEMCALTriggerQA::CreateTH1", "Object %s already exists", name);
    return;
  }
  TH1* hist = new TH1D(name, title, nbins, xmin, xmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH2(const char *name, const char *title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax)
{
  /*
   * Create a new TH2 within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param title: Title of the histogram
   * @param nbinsx: number of bins in x-direction
   * @param xmin: min. value of the range in x-direction
   * @param xmax: max. value of the range in x-direction
   * @param nbinsy: number of bins in y-direction
   * @param ymin: min. value of the range in y-direction
   * @param ymax: max. value of the range in y-direction
   * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEMCALTriggerQA::CreateTH2", "Object %s already exists", name);
    return;
  }
  TH2* hist = new TH2D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::CreateTH3(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, int nbinsz, double zmin, double zmax)
{
  /*
   * Create a new TH3 within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param nbinsx: number of bins in x-direction
   * @param xmin: min. value of the range in x-direction
   * @param xmax: max. value of the range in x-direction
   * @param nbinsy: number of bins in y-direction
   * @param ymin: min. value of the range in y-direction
   * @param ymax: max. value of the range in y-direction
   * @param nbinsz: number of bins in z-direction
   * @param zmin: min. value of the range in z-direction
   * @param zmax: max. value of the range in z-direction
   * Raises fatals in case the parent group does not exist or the object is attempted to be duplicated within the group
   */
  if (fHistos->FindObject(name)) {
    Fatal("AliEMCALTriggerQA::CreateTH3", "Object %s already exists", name);
    return;
  }
  TH3* hist = new TH3D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
  fHistos->Add(hist);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH1(const char *name, double x, double weight)
{
  /*
   * Fill a 1D histogram within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
   */
  TH1* hist = dynamic_cast<TH1*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEMCALTriggerQA::FillTH1", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, weight);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH2(const char *name, double x, double y, double weight)
{
  /*
   * Fill a 2D histogram within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param y: y-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
   */
  TH2* hist = dynamic_cast<TH2*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEMCALTriggerQA::FillTH2", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, weight);
}

//______________________________________________________________________________
void AliEMCALTriggerQA::FillTH3(const char* name, double x, double y, double z, double weight)
{
  /*
   * Fill a 3D histogram within the container. The histogram name also contains the parent group(s) according to the common
   * group notation.
   *
   * @param name: Name of the histogram
   * @param x: x-coordinate
   * @param y: y-coordinate
   * @param z: z-coordinate
   * @param weight (@default 1): optional weight of the entry
   * Raises fatals in case the parent group is not found or the histogram is not found in the parent group
   */

  TH3* hist = dynamic_cast<TH3*>(fHistos->FindObject(name));
  if (!hist) {
    Fatal("AliEMCALTriggerQA::FillTH3", "Histogram %s not found", name);
    return;
  }
  hist->Fill(x, y, z, weight);
}

//______________________________________________________________________________
TObject *AliEMCALTriggerQA::FindObject(const char *name) const
{
  /*
   * Find an object inside the container. The object can also be within a
   * histogram group. For this the name has to follow the common notation
   *
   * @param name: Name of the object to find inside the container
   * @return: pointer to the object (NULL if not found)
   */

  return fHistos->FindObject(name);
}

//______________________________________________________________________________
TObject* AliEMCALTriggerQA::FindObject(const TObject* obj) const
{
  /*
   * Find and object inside the container. The object name is expected to contain the
   * full path of the histogram object, including parent groups
   *
   * @param obj: the object to find
   * @return: pointer to the object (NULL if not found)
   */
  TString hname(obj->GetName());
  return fHistos->FindObject(hname);
}
