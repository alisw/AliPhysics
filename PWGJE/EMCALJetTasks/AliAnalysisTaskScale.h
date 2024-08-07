#ifndef AliAnalysisTaskScale_h
#define AliAnalysisTaskScale_h

// $Id$

#include "AliAnalysisTaskEmcal.h"
#include "TH2.h"
#include "TF1.h"
#include "AliEmcalJet.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

class AliAnalysisTaskScale : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskScale();
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  void                   UserCreateOutputObjects();

  void                   SetScaleFunction(TF1* sf)  { fScaleFunction = sf   ; }

  static AliAnalysisTaskScale* AddTaskScale(
    const TString nTracks        = "Tracks",
    const TString nClusters      = "CaloClustersCorr",
    const Double_t trackptcut    = 0.150,
    const Double_t clusptcut     = 0.150,
    const TString taskname       = "Scale",
    const char *sfuncPath        = 0,
    const char *sfuncName        = 0
  );
  
 protected:
  void                   ExecOnce();
  Double_t               GetScaleFactor(Double_t cent);
  Bool_t                 FillHistograms();

#if !(defined(__CINT__) || defined(__MAKECINT__))
  // Handle mapping between index and containers
  AliEmcalContainerIndexMap <AliClusterContainer, AliVCluster> fClusterContainerIndexMap;    //!<! Mapping between index and cluster containers
  AliEmcalContainerIndexMap <AliParticleContainer, AliVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
#endif

 private:
  TF1                   *fScaleFunction;                 // scale factor as a function of centrality

  Double_t               fEmcalArea;                     //!Emcal area
  Double_t               fTpcArea;                       //!Tpc area

  TH2                   *fHistPtTPCvsCent;               //!output histogram
  TH2                   *fHistPtEMCALvsCent;             //!output histogram
  TH2                   *fHistEtvsCent;                  //!output histogram
  TH2                   *fHistScalevsCent;               //!output histogram
  TH2                   *fHistDeltaScalevsCent;          //!output histogram
  TH2                   *fHistScaleEmcalvsCent;          //!output histogram
  TH2                   *fHistScale2EmcalvsCent;         //!output histogram
  TH2                   *fHistDeltaScale2EmcalvsCent;    //!output histogram
  TH2                   *fHistScale3EmcalvsCent;         //!output histogram
  TH2                   *fHistDeltaScale3EmcalvsCent;    //!output histogram
  TH2                   *fHistChScalevsCent;             //!output histogram
  TH2                   *fHistChScale2EmcalvsCent;       //!output histogram
  TH2                   *fHistChScale3EmcalvsCent;       //!output histogram
  TH2                   *fHistPtTPCvsNtrack;             //!output histogram
  TH2                   *fHistPtEMCALvsNtrack;           //!output histogram
  TH2                   *fHistEtvsNtrack;                //!output histogram
  TH2                   *fHistScalevsNtrack;             //!output histogram
  TH2                   *fHistDeltaScalevsNtrack;        //!output histogram
  TH2                   *fHistScaleEmcalvsNtrack;        //!output histogram
  TH2                   *fHistScale2EmcalvsNtrack;       //!output histogram
  TH2                   *fHistScale3EmcalvsNtrack;       //!output histogram
  TH2                   *fHistChScalevsNtrack;           //!output histogram
  TH2                   *fHistChScale2EmcalvsNtrack;     //!output histogram
  TH2                   *fHistChScale3EmcalvsNtrack;     //!output histogram
  TH2                   *fHistTrackPtvsCent;             //!output histogram
  TH2                   *fHistClusterPtvsCent;           //!output histogram
  TH2                   *fHistTrackEtaPhi;               //!output histogram
  TH2                   *fHistClusterEtaPhi;             //!output histogram
  TH2                   *fHistScalevsScale2Emcal;        //!output histogram
  TH2                   *fHistScalevsScale3Emcal;        //!output histogram
  TH2                   *fHistScalevsScaleEmcal;         //!output histogram
  TH2                   *fHistScaleEmcalvsScale2Emcal;   //!output histogram
  TH2                   *fHistScaleEmcalvsScale3Emcal;   //!output histogram
  TH2                   *fHistScaleShift1EmcalvsCent;    //!output histogram
  TH2                   *fHistScaleShift2EmcalvsCent;    //!output histogram
  TH2                   *fHistScaleShiftMeanEmcalvsCent; //!output histogram
  TH2                   *fHistScaleEmcalvsPhi;           //!output histogram
  TH2                   *fHistPtEmcalvsPhi;              //!output histogram
  TH2                   *fHistPtvsPhi;                   //!output histogram
  TH2                   *fHistEtaPhiScaleEMCAL;          //!output histogram
  TH2                   *fHistEtaPhiScale2EMCAL;         //!output histogram
  TH2                   *fHistEtaPhiScale3EMCAL;         //!output histogram
  TH2                   *fHistEtaPhiScaleShift1Emcal;    //!output histogram
  TH2                   *fHistEtaPhiScaleShift2Emcal;    //!output histogram
  TH2                   *fHistEtaPhiScaleShiftMeanEmcal; //!output histogram
  TH2                   *fHistEtaPhiScaleEmcalvsPhi;     //!output histogram

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 11); // Scale task
};
#endif
