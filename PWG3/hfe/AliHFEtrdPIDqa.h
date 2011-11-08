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
//
// QA class for TRD PID
// Evaluate TRD PID using well identified reference tracks
// For more information see implementation file
//
#ifndef ALIHFETRDPIDQA_H
#define ALIHFETRDPIDQA_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

#ifndef ROOT_TH2
#include <TH2.h>
#endif

#ifndef ALIHFECOLLECTION_H
#include "AliHFEcollection.h"
#endif

class TArrayD;
class TBrowser;
class TCollection;
class TF1;
class TGraph;
class TGraphErrors;
class TH1;
//class TH2;
class TList;
class TObjArray;

class AliAODTrack;
class AliESDtrack;
class AliVTrack;
class AliHFEpidTRD;

class AliHFEtrdPIDqa : public TNamed{
  public:
    AliHFEtrdPIDqa();
    AliHFEtrdPIDqa(const Char_t *name);
    AliHFEtrdPIDqa(const AliHFEtrdPIDqa &ref);
    AliHFEtrdPIDqa &operator=(const AliHFEtrdPIDqa &ref);
    virtual void Copy(TObject &o) const;
    virtual Long64_t Merge(TCollection *coll);
    virtual ~AliHFEtrdPIDqa();
   
    virtual void Browse(TBrowser *b);
    virtual Bool_t IsFolder() const { return kTRUE; }

    void ProcessTracks(const TObjArray * const  l, Int_t species);
    void ProcessTrack(const AliVTrack *const track, Int_t species);

    void Init();
    void FinishAnalysis();
    void ShowMessages() { fShowMessage = kTRUE; }
    void SetTotalChargeInSlice0() { fTotalChargeInSlice0 = kTRUE; }
    void StoreResults(const Char_t *filename = "HFEtrdPIDqa.root");
    void SaveThresholdParameters(const Char_t * filename = "TRD.Thresholds.root", Double_t lowerLimit = 0.5, Double_t upperLimit = 3.5);

    TGraphErrors *GetPionEfficiency(Int_t ntracklets, Int_t eleffpercent);
    TGraphErrors *GetProtonEfficiency(Int_t ntracklets, Int_t eleffpercent);
    TGraph *GetThreshold(Int_t ntracklets, Int_t eleffpercent);

    void DrawTracklet(Int_t tracklet, Double_t pmin = -1., Double_t pmax = -1., Bool_t doFit = kFALSE);
    void ClearLists();

    Double_t EvalPionEfficiency(Int_t ntls, Int_t eEff, Double_t p);
    Double_t EvalProtonEfficiency(Int_t ntls, Int_t eEff, Double_t p);
    Double_t EvalThreshold(Int_t ntls, Int_t eEff, Double_t p);

    Double_t CalculateIntegratedPionEfficiency(UInt_t nTracklets, Double_t electronEff, Double_t pmin, Double_t pmax, Double_t *error = NULL);

    //---------------------------------------------------
    // Getters for Histograms
    THnSparseF *GetLikelihoodHistogram() const { return fHistos ? dynamic_cast<THnSparseF *>(fHistos->Get("fLikeTRD")) : NULL; }
    THnSparseF *GetQAHistogram() const { return fHistos ? dynamic_cast<THnSparseF *>(fHistos->Get("fQAtrack")) : NULL; }
    THnSparseF *GetdEdxHistogram() const { return fHistos ? dynamic_cast<THnSparseF *>(fHistos->Get("fQAdEdx")) : NULL; }
    THnSparseF *GetHistoTruncMean() const { return fHistos ? dynamic_cast<THnSparseF *>(fHistos->Get("fTRDtruncMean")) : NULL; }
    TH2 *GetSliceChargePions(Bool_t pions = kTRUE) const { 
      const Char_t * species = pions ? "Pions" : "Electrons";
      return fHistos ? dynamic_cast<TH2 *>(fHistos->Get(Form("fTRDslices%s", species))) : NULL; 
    }
    AliHFEcollection *GetHistos() const { return fHistos; }
    //---------------------------------------------------
  protected:
    // Description of the containers we use to store basic information
    // we access in the Post Processing. For all containers we have a 
    // common part containing species, momentum and number of tracklets,
    // and a specific part. For both containers we define the number of 
    // variables too
    enum QuantitiesCommon_t{
      kSpecies = 0,
      kP = 1,
      kNTracklets = 2,
      kQuantitiesCommon = 3
    };
    enum QuantitiesLike_t{
      kElectronLike = 3,
      kNClustersLike = 4,
      kQuantitiesLike = 5
    };
    enum QuantitiesQAtrack_t{
      kNonZeroTrackletCharge = 3,
      kNClusters = 4,
      kQuantitiesQA = 5
    };
    enum QuantitiesdEdx_t{
      kNclusters = 3,
      kNonZeroSlices = 4,
      kdEdx = 5,
      kQuantitiesdEdx = 6
    };
    enum QuantitiesTruncMean_t{
      kTPCdEdx = 3,
      kTRDdEdxMethod1 = 4,
      kTRDdEdxMethod2 = 5,
      kQuantitiesTruncMean = 6
    };

    void ProcessTrackESD(const AliESDtrack * const track, Int_t species);
    void ProcessTrackAOD(const AliAODTrack * const track, Int_t species);

    void FillTRDLikelihoods(const AliESDtrack * const track, Int_t species);
    void FillTRDQAplots(const AliESDtrack *const track, Int_t species);

    void AnalyseNTracklets(Int_t nTracklets);
    Double_t CalculateHadronEfficiency(const TH1 * const hadron, const TH1 *const electron, Double_t eff, Int_t &threshbin, Double_t &error);
    TF1 *MakeThresholds(TGraph *input, Double_t lowerLimit, Double_t upperLimit);

    void CreateLikelihoodHistogram();
    void CreateQAHistogram();
    void CreatedEdxHistogram();
    void CreateHistoTruncatedMean();

  private:
    enum{
      kNElectronEffs = 6
    };
    static const Double_t fgkElectronEff[kNElectronEffs];       // Electron efficiency bins
    static const Int_t    fgkNBinsCommon[kQuantitiesCommon];    // Number of bins for common quantities
    static const Double_t fgkMinBinCommon[kQuantitiesCommon];   // Bin Limits for common quantities (lower limit)
    static const Double_t fgkMaxBinCommon[kQuantitiesCommon];   // Bin Limits for common quantities (upper limit)
    AliHFEpidTRD *fTRDpid;        // HFE PID for TRD
    AliHFEcollection *fHistos;    // Histogram collection

    // List for Histograms:
    TList *fPionEfficiencies;     //! List for Pion efficiencies
    TList *fProtonEfficiencies;   //! List for Proton efficiencies
    TList *fKaonEfficiencies;     //! List for Kaon efficiencies

    TList *fThresholds;           //! List for Threshold Graphs

    Bool_t fShowMessage;         // Display debug messages
    Bool_t fTotalChargeInSlice0;  // Flag for Foreward/Backward compatibility in TRD total charge calculation
  
  ClassDef(AliHFEtrdPIDqa, 3)     // QA class for TRD PID 
};
#endif

