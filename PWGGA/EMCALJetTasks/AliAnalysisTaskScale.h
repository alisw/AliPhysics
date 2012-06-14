#ifndef AliAnalysisTaskScale_h
#define AliAnalysisTaskScale_h

// $Id$

class TList;
class TH1F;
class TH2F;
class TF1;
class AliEMCALGeometry;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskScale : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskScale();
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  void                   UserCreateOutputObjects();
  void                   Terminate(Option_t *);

  void                   SetScaleFunction(TF1* sf)  { fScaleFunction = sf   ; }
  
 protected:
  virtual Double_t       GetScaleFactor(Double_t cent);
  virtual Bool_t         FillHistograms();
  void                   Init();

 private:
  TF1                   *fScaleFunction;          // scale factor as a function of centrality

  AliEMCALGeometry      *fGeom;                   //!ptr to emcal geometry object
  TH1F                  *fHistCentrality;         //!output histogram
  TH2F                  *fHistPtTPCvsCent;        //!output histogram
  TH2F                  *fHistPtEMCALvsCent;      //!output histogram
  TH2F                  *fHistEtvsCent;           //!output histogram
  TH2F                  *fHistScalevsCent;        //!output histogram
  TH2F                  *fHistDeltaScalevsCent;   //!output histogram
  TH2F                  *fHistPtTPCvsNtrack;      //!output histogram
  TH2F                  *fHistPtEMCALvsNtrack;    //!output histogram
  TH2F                  *fHistEtvsNtrack;         //!output histogram
  TH2F                  *fHistScalevsNtrack;      //!output histogram
  TH2F                  *fHistDeltaScalevsNtrack; //!output histogram
  TH2F                  *fHistTrackPtvsCent;      //!output histogram
  TH2F                  *fHistClusterPtvsCent;    //!output histogram
  TH2F                  *fHistTrackEtaPhi;        //!output histogram
  TH2F                  *fHistClusterEtaPhi;      //!output histogram

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 7); // Scale task
};
#endif
