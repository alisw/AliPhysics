#ifndef AliAnalysisTaskScale_h
#define AliAnalysisTaskScale_h

// $Id$

class TH2;
class TF1;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskScale : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskScale();
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  void                   UserCreateOutputObjects();

  void                   SetScaleFunction(TF1* sf)  { fScaleFunction = sf   ; }
  
 protected:
  void                   ExecOnce();
  Double_t               GetScaleFactor(Double_t cent);
  Bool_t                 FillHistograms();

 private:
  TF1                   *fScaleFunction;               // scale factor as a function of centrality

  Double_t               fEmcalArea;                   //!Emcal area
  Double_t               fTpcArea;                     //!Tpc area

  TH2                   *fHistPtTPCvsCent;             //!output histogram
  TH2                   *fHistPtEMCALvsCent;           //!output histogram
  TH2                   *fHistEtvsCent;                //!output histogram
  TH2                   *fHistScalevsCent;             //!output histogram
  TH2                   *fHistDeltaScalevsCent;        //!output histogram
  TH2                   *fHistScaleEmcalvsCent;        //!output histogram
  TH2                   *fHistScale2EmcalvsCent;       //!output histogram
  TH2                   *fHistDeltaScale2EmcalvsCent;  //!output histogram
  TH2                   *fHistChScalevsCent;           //!output histogram
  TH2                   *fHistChScale2EmcalvsCent;     //!output histogram
  TH2                   *fHistPtTPCvsNtrack;           //!output histogram
  TH2                   *fHistPtEMCALvsNtrack;         //!output histogram
  TH2                   *fHistEtvsNtrack;              //!output histogram
  TH2                   *fHistScalevsNtrack;           //!output histogram
  TH2                   *fHistDeltaScalevsNtrack;      //!output histogram
  TH2                   *fHistScaleEmcalvsNtrack;      //!output histogram
  TH2                   *fHistScale2EmcalvsNtrack;     //!output histogram
  TH2                   *fHistChScalevsNtrack;         //!output histogram
  TH2                   *fHistChScale2EmcalvsNtrack;   //!output histogram
  TH2                   *fHistTrackPtvsCent;           //!output histogram
  TH2                   *fHistClusterPtvsCent;         //!output histogram
  TH2                   *fHistTrackEtaPhi;             //!output histogram
  TH2                   *fHistClusterEtaPhi;           //!output histogram
  TH2                   *fHistScalevsScale2Emcal;      //!output histogram
  TH2                   *fHistScalevsScaleEmcal;       //!output histogram
  TH2                   *fHistScaleEmcalvsScale2Emcal; //!output histogram

  AliParticleContainer  *fTracksCont;                  //!Tracks
  AliClusterContainer   *fCaloClustersCont;            //!Clusters 

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 11); // Scale task
};
#endif
