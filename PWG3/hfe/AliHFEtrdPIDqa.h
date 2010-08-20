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

class TCollection;
class TF1;
class TGraph;
class TH1;
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

    void ProcessTracks(TObjArray * const  l, Int_t species);
    void ProcessTrack(AliVTrack *track, Int_t species);

    void Init();
    void FinishAnalysis();
    void StoreResults(const Char_t *filename = "HFEtrdPIDqa.root");
    void SaveThresholdParameters(const Char_t * filename = "TRD.Thresholds.root");

    void DrawTracklet(Int_t tracklet);
    void ClearLists();

    //---------------------------------------------------
    // Getters for Histograms
    THnSparseF *GetLikelihoodHistogram() const { return fLikeTRD; }
    THnSparseF *GetQAHistogram() const { return fQAtrack; }
    THnSparseF *GetdEdxHistogram() const { return fQAdEdx; }
    THnSparseF *GetHistoTruncMean() const { return fTRDtruncMean; }
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
      kQuantitiesLike = 4
    };
    enum QuantitiesQAtrack_t{
      kNonZeroTrackletCharge = 3,
      kNClusters = 4,
      kQuantitiesQA = 5
    };
    enum QuantitiesdEdx_t{
      kdEdx = 3,
      kQuantitiesdEdx = 4
    };
    enum QuantitiesTruncMean_t{
      kTPCdEdx = 3,
      kTRDdEdxMethod1 = 4,
      kTRDdEdxMethod2 = 5,
      kQuantitiesTruncMean = 6
    };

    void ProcessTrackESD(AliESDtrack *track, Int_t species);
    void ProcessTrackAOD(AliAODTrack * const track, Int_t species);

    void FillTRDLikelihoods(AliESDtrack *track, Int_t species);
    void FillTRDQAplots(AliESDtrack *track, Int_t species);

    void AnalyseNTracklets(Int_t nTracklets);
    Int_t GetThresholdBin(TH1 * const input, Double_t efficiency);
    Bool_t CalculateEfficiency(TH1 * const input, Int_t threshbin, Double_t *params);
    TF1 *MakeThresholds(TGraph *input);

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
    THnSparseF *fLikeTRD;         // Histo for Likelihoods
    THnSparseF *fQAtrack;         // QA histo for quantities based on track level
    THnSparseF *fQAdEdx;          // QA for tracklet charge
    THnSparseF *fTRDtruncMean;    // QA for truncated mean

    // List for Histograms:
    TList *fPionEfficiencies;     //! List for Pion efficiencies
    TList *fProtonEfficiencies;   //! List for Proton efficiencies
    TList *fKaonEfficiencies;     //! List for Kaon efficiencies

    TList *fThresholds;           //! List for Threshold Graphs
  
  ClassDef(AliHFEtrdPIDqa, 2)     // QA class for TRD PID 
};
#endif

