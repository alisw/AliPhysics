#ifndef ALIHFETRACKFILTER_H
#define ALIHFETRACKFILTER_H

/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */ 

//
// Track filter class 
// Apply cut steps to all tracks in one event and returns a list of
// filtered tracks
//
#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_TObjArray
#include <TObjArray.h>
#endif

class AliCFContainer;
class AliESDEvent;
class AliHFEcontainer;
class AliHFEcutStep;
class AliMCEvent;

class AliHFEtrackFilter : public TNamed{
  public:
    AliHFEtrackFilter(const Char_t *name);
    AliHFEtrackFilter(const AliHFEtrackFilter &o);
    AliHFEtrackFilter &operator=(const AliHFEtrackFilter &o);
    virtual void Copy(TObject &o) const;
    ~AliHFEtrackFilter();
  
    // Base Functionality
    void AddCutStep(AliHFEcutStep *cuts);
    void GenerateCutSteps();
    void SetRecEvent(AliVEvent *rec);
    AliHFEcutStep *GetCutStep(Int_t istep);
    AliHFEcutStep *GetCutStep(const Char_t *name);
    void FilterTracks(AliESDEvent *const esd);
    TObjArray *GetFilteredTracks() const { return fFilteredTracks; }
    void Flush();

    // Usage of the Correction Framework
    void InitCF();
    void InitCF(AliHFEcontainer *cont);
    AliCFContainer *GetEfficiencyContainer(Int_t icont); 
    void ReleaseContainers() { SetBit(kOwnCFContainers, kFALSE); }
    void OwnContainers() { SetBit(kOwnCFContainers, kTRUE); }
    void SetPtBins(Int_t nBins, Double_t *binning);
    void SetEtaBins(Int_t nBins, Double_t *binning);
    void SetPhiBins(Int_t nBins, Double_t *binning);

    // Add Possibility to check Monte-Carlo Information
    void SetMC(AliMCEvent * const mc);
    void SetMCSignalStep(AliHFEcutStep * const sig) { fMCsignal = sig; };
    AliHFEcutStep *GetMCSignalCuts() const { return fMCsignal; }

    // Creators for Standard Cut Steps
    AliHFEcutStep *MakeCutStepRecKineITSTPC();
    AliHFEcutStep *MakeCutStepPrimary();
    AliHFEcutStep *MakeCutStepHFEITS();
    AliHFEcutStep *MakeCutStepHFETRD();
    AliHFEcutStep *MakeMCSignalCuts();

  private:
    enum{
      kOwnCFContainers = BIT(14)
    };
    TObjArray *fFilteredTracks;             //! List of Tracks which survived the filter
    TObjArray *fCutSteps;                   //! Container of Cut Objects
    TObjArray *fEfficiencyContainers;       //! Efficiency Container
    AliMCEvent * fMC;                       //! Monte Carlo Event
    AliHFEcutStep *fMCsignal;               //! Cuts for MC Signal

    // Helpers for the correction framework
    Int_t fPtBins;                          // Number of Pt Bins
    Int_t fEtaBins;                         // Number of Eta Bins
    Int_t fPhiBins;                         // Number of Phi Bins
    Double_t *fPtBinning;                   // Pt binning
    Double_t *fEtaBinning;                  // Eta binning
    Double_t *fPhiBinning;                  // Phi binning

    ClassDef(AliHFEtrackFilter, 1)
};
#endif

