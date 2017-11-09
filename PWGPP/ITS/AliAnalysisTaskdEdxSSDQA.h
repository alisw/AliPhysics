//SSD dEdX QA task
//Marek Chojnacki
//Marek.Chojnacki@cern.ch
#ifndef ALIANALYSISTASKDEDXSSDQA_H
#define ALIANALYSISTASKDEDXSSDQA_H

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
class AliITSCalibrationSSD;
class TList;

class AliAnalysisTaskdEdxSSDQA : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskdEdxSSDQA(const char *name = "AliAnalysisTaskdEdxSSDQA");
  virtual ~AliAnalysisTaskdEdxSSDQA() ;
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   LocalInit();
  
  
  void SetPcut(Float_t pcut){ fPcut=pcut;}
  Float_t GetPcut() const{return fPcut;}
  void SetDoChipCorretions(const char* filename);
 private:

  TH1F*   fHistBadMods;   // bad module map
  TH1F*   fHistBadPstr;   // bad p strip
  TH1F*   fHistBadNstr;   // bad n strip
  TH2F*   fHist1;         // CR for each module
  TH2F*   fHist2;         // landau distributions for each module	
  TH3F*   fHist3;         // CR as function of Charge for the AliTrackPoint 
  TH2F*   fHist4;         // Q on chips
  TH2F*   fHist5;         // Q on chips corrected
  TH2F*   fHist6;         // QNvQP not corrected for track inclinaition 
  TH2F*   fHist7;         // zeta-phi of points in lay 5 
  TH2F*   fHist8;         // zeta-phi of points in lay 6 
  TH2F*   fHist9;         // detec-lad of points in lay 5
  TH2F*   fHist10;         // detec-lad of points in lay 6


  TH2F*   fHist1sa;         // SA tracks: CR for each module
  TH2F*   fHist2sa;         // SA tracks: landau distributions for each module	
  TH3F*   fHist3sa;         // SA tracks: CR as function of Charge for the AliTrackPoint 
  TH2F*   fHist4sa;         // SA tracks: Q on chips
  TH2F*   fHist5sa;         // SA tracks: Q on chips corrected
  TH2F*   fHist6sa;         // SA tracks: QNvQP not corrected for track inclinaition 
  TH2F*   fHist7sa;         // zeta-phi of points in lay 5 
  TH2F*   fHist8sa;         // zeta-phi of points in lay 6 
  TH2F*   fHist9sa;         // detec-lad of points in lay 5
  TH2F*   fHist10sa;         // detec-lad of points in lay 6

  TList*  fListOfHistos;  // output list	
  Float_t fPcut;          // Momentum cut

  Bool_t fdothecorrection; //do the correction  	
  Float_t fcorrections[20376] ; //[20376]chip corrections
  Bool_t fInitCalib;        //flag for calib initialization
  AliITSCalibrationSSD* fSSDCalibration; //SSD calibration object from OCDB

 AliAnalysisTaskdEdxSSDQA(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 AliAnalysisTaskdEdxSSDQA& operator=(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 Int_t Pstrip5(Float_t x,Float_t z) const;
 Int_t Pstrip6(Float_t x,Float_t z) const;
 Int_t Nstrip5(Float_t x,Float_t z) const;
 Int_t Nstrip6(Float_t x,Float_t z) const;

 ClassDef(AliAnalysisTaskdEdxSSDQA, 4); // example of analysis
};

#endif
