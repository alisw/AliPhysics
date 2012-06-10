#ifndef AliAnalysisTaskScale_h
#define AliAnalysisTaskScale_h

// $Id$

class TList;
class TH1F;
class TH2F;
class TF1;
class AliEMCALGeometry;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskScale : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskScale();
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);
  virtual void           Terminate(Option_t *);

  void                   SetClustersName(const char *n)                        { fClustersName  = n    ; }
  void                   SetMinClusterPt(Double_t min)                         { fMinClusterPt  = min  ; }
  void                   SetMinTrackPt(Double_t min)                           { fMinTrackPt    = min  ; }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction = sf   ; }
  void                   SetTracksName(const char *n)                          { fTracksName    = n    ; }
  
 protected:
  virtual Double_t       GetScaleFactor(Double_t cent);

 private:
  TString                fTracksName;             // name of track collection
  TString                fClustersName;           // name of clusters collection
  Double_t               fMinTrackPt;             // pt cut for scale factor calculation
  Double_t               fMinClusterPt;           // pt cut for scale factor calculation
  TF1                   *fScaleFunction;          // scale factor as a function of centrality
  AliEMCALGeometry      *fGeom;                   //!ptr to emcal geometry object
  TList                 *fOutputList;             //!output list
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
  
  ClassDef(AliAnalysisTaskScale, 6); // Scale task
};
#endif
