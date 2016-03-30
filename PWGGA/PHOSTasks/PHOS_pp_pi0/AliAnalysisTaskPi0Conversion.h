#ifndef AliAnalysisTaskPi0Conversion_cxx
#define AliAnalysisTaskPi0Conversion_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliPHOSCalibData ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0Conversion : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0Conversion(const char *name = "AliAnalysisTaskPi0Conversion");
  virtual ~AliAnalysisTaskPi0Conversion() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *){}
  
protected:
  AliAnalysisTaskPi0Conversion(const AliAnalysisTaskPi0Conversion&); // not implemented
  AliAnalysisTaskPi0Conversion& operator=(const AliAnalysisTaskPi0Conversion&); // not implemented
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  Int_t FindAODLabel(Int_t esdLabel)const ;
 
private:
  
  THashList * fOutputContainer;   //final histogram container
  TList * fPHOSEvents ;           //Container for PHOS photons
  TClonesArray * fPHOSEvent ;     //PHOS photons in current event
  TClonesArray *fStack;           //! Array of primary particles 
 
  AliPHOSGeometry  *fPHOSGeo;     //! PHOS geometry
 
  ClassDef(AliAnalysisTaskPi0Conversion, 1); // PHOS analysis task
};

#endif
