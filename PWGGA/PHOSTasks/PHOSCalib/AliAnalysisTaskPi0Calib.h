#ifndef AliAnalysisTaskPi0Calib_cxx
#define AliAnalysisTaskPi0Calib_cxx

// PHOS calibration with pi0
// Authors: DP, YK

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH2D;
class TH3F;
class TF1 ;
class AliPHOSGeometry;
class AliESDEvent ;
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0Calib : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0Calib(const char *name = "PHOS");
  virtual ~AliAnalysisTaskPi0Calib() ;
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetPHOSBadMap(Int_t mod,TH2I * h)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod]=new TH2I(*h) ;
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }

  void SetPHOSCalib(Int_t mod,TH2D * h)
  {
    if(fPHOSAmp[mod]) delete fPHOSAmp[mod] ;
    fPHOSAmp[mod]=new TH2D(*h) ;
    printf("Set %s \n",fPHOSAmp[mod]->GetName());
  }
  
  
private:

  AliAnalysisTaskPi0Calib(const AliAnalysisTaskPi0Calib&); // not implemented
  AliAnalysisTaskPi0Calib& operator=(const AliAnalysisTaskPi0Calib&); // not implemented
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  Bool_t TestTOF(Double_t e,Double_t tof) ;  //Evaluate TOF cut
  Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
  Int_t ConvertRunNumber(Int_t run) ; 

private:
  TList            * fOutputContainer; //! final histogram container
  TList            * fPHOSEvents;      //! Containers for events with PHOS photons
  TClonesArray     * fPHOSEvent ;      //! PHOS photons in current event
  AliPHOSCalibData * fPHOSCalibData;   //! PHOS calibration object
  TH2I             * fPHOSBadMap[6] ;  //Container for PHOS bad channels map
  TH2D             * fPHOSAmp[6] ;     //Container for PHOS bad channels map

  AliPHOSGeometry  * fPHOSGeo;         //! PHOS geometry
  Int_t              fEventCounter;    // number of analyzed events

  ClassDef(AliAnalysisTaskPi0Calib, 1); // PHOS analysis task
};

#endif
