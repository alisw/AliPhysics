#ifndef AliAnalysisTaskPt_cxx
#define AliAnalysisTaskPt_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class AliESDEvent;
// class AliEMCALGeometry;
// class AliEMCALGeoUtils;

// #include "AliEMCALGeometry.h"
// #include "AliEMCALGeoUtils.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPt : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPt() : AliAnalysisTaskSE(), fESD(0),
  fOutputList(0), fMyTr(0), jev(0), iev(0) {} 

  AliAnalysisTaskPt(const char *name);
  virtual ~AliAnalysisTaskPt() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TTree       *fMyTr;
  int          jev, iev;
/*
  TH1F        *fHstPt;
  TH1F        *fEMCe;
  TH1F        *fEMCt;
  TH1F        *fEMCn;
  TH2F        *fEMCm;
  TH1F        *fCelle;
  TH1F        *fCellf;
  TH1F        *fCellt;
  TH2F        *fCellc;
  TH2F        *fCelld;
  TH2F        *fCellm;
//AliEMCALGeometry *fEMCALGeo;
  AliEMCALGeoUtils *fEMCALGeo;
*/
   
  AliAnalysisTaskPt(const AliAnalysisTaskPt&); // not implemented
  AliAnalysisTaskPt& operator=(const AliAnalysisTaskPt&); // not implemented
  
  ClassDef(AliAnalysisTaskPt, 1); // example of analysis
};

#endif
