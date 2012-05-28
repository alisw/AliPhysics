#ifndef AliAnalysisTaskScale_cxx
#define AliAnalysisTaskScale_cxx

// $Id$

class TList;
class TH1F;
class TH2F;
class TF1;
class AliEMCALGeometry;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskScale : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskScale() : AliAnalysisTaskSE(), fTracksName(), fClustersName(), fScaleFunction(0), fGeom(0),
    fOutputList(0), fHistCentrality(0), fHistPtTPCvsCent(0), fHistPtEMCALvsCent(0), fHistEtvsCent(0),  
    fHistScalevsCent(0),  fHistDeltaScalevsCent(0), fHistPtTPCvsNtrack(0), fHistPtEMCALvsNtrack(0), 
    fHistEtvsNtrack(0), fHistScalevsNtrack(0), fHistDeltaScalevsNtrack(0), fHistTrackPtvsCent(0), 
    fHistClusterPtvsCent(0), fHistTrackEtaPhi(0), fHistClusterEtaPhi(0),  fMinTrackPt(0.15), 
    fMinClusterPt(0.15) {}
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  virtual void           UserCreateOutputObjects();
  virtual void           UserExec(Option_t *option);
  virtual void           Terminate(Option_t *);

  void                   SetTracksName(const char *n)                          { fTracksName    = n    ; }
  void                   SetClustersName(const char *n)                        { fClustersName  = n    ; }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction = sf   ; }
  void                   SetMinTrackPt(Double_t min)                           { fMinTrackPt    = min  ; }
  void                   SetMinClusterPt(Double_t min)                         { fMinClusterPt  = min  ; }
  
 protected:
  virtual Double_t       GetScaleFactor(Double_t cent);

 private:
  TString                fTracksName;             // name of track collection
  TString                fClustersName;           // name of clusters collection
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
  Double_t               fMinTrackPt;             //pt cut for scale factor calculation
  Double_t               fMinClusterPt;           //pt cut for scale factor calculation

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 4); // Scale task
};
#endif
