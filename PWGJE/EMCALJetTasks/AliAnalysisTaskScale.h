#ifndef AliAnalysisTaskScale_h
#define AliAnalysisTaskScale_h

// $Id$

class TH2F;
class TF1;

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

  TH2F                  *fHistPtTPCvsCent;             //!output histogram
  TH2F                  *fHistPtEMCALvsCent;           //!output histogram
  TH2F                  *fHistEtvsCent;                //!output histogram
  TH2F                  *fHistScalevsCent;             //!output histogram
  TH2F                  *fHistDeltaScalevsCent;        //!output histogram
  TH2F                  *fHistScaleEmcalvsCent;        //!output histogram
  TH2F                  *fHistScale2EmcalvsCent;       //!output histogram
  TH2F                  *fHistChScalevsCent;           //!output histogram
  TH2F                  *fHistChScale2EmcalvsCent;     //!output histogram
  TH2F                  *fHistPtTPCvsNtrack;           //!output histogram
  TH2F                  *fHistPtEMCALvsNtrack;         //!output histogram
  TH2F                  *fHistEtvsNtrack;              //!output histogram
  TH2F                  *fHistScalevsNtrack;           //!output histogram
  TH2F                  *fHistDeltaScalevsNtrack;      //!output histogram
  TH2F                  *fHistScaleEmcalvsNtrack;      //!output histogram
  TH2F                  *fHistScale2EmcalvsNtrack;     //!output histogram
  TH2F                  *fHistChScalevsNtrack;         //!output histogram
  TH2F                  *fHistChScale2EmcalvsNtrack;   //!output histogram
  TH2F                  *fHistTrackPtvsCent;           //!output histogram
  TH2F                  *fHistClusterPtvsCent;         //!output histogram
  TH2F                  *fHistTrackEtaPhi;             //!output histogram
  TH2F                  *fHistClusterEtaPhi;           //!output histogram
  TH2F                  *fHistScalevsScale2Emcal;      //!output histogram
  TH2F                  *fHistScalevsScaleEmcal;       //!output histogram
  TH2F                  *fHistScaleEmcalvsScale2Emcal; //!output histogram

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 10); // Scale task
};
#endif
