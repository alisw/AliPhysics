#ifndef AliAnalysisTaskEpRatio_cxx
#define AliAnalysisTaskEpRatio_cxx

// E/p analysis task.
// Authors: Boris Polishchuk, Tsubasa Okubo

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliPHOSGeometry;
class AliAnalysisTaskSE;
class AliPIDResponse;


#include "TH2I.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEpRatio : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskEpRatio(const char *name = "AliAnalysisTaskEpRatio");
  virtual ~AliAnalysisTaskEpRatio() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
private:
  AliAnalysisTaskEpRatio(const AliAnalysisTaskEpRatio&); // not implemented
  AliAnalysisTaskEpRatio& operator=(const AliAnalysisTaskEpRatio&); // not implemented

  void SetGeometry();
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
 
private:

  Int_t fRunNumber;
  TList * fOutputContainer;       //final histogram container
  AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
  AliPIDResponse *fPIDResponse; // PID Response

  ClassDef(AliAnalysisTaskEpRatio, 2); // PHOS analysis task
};

#endif
