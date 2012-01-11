//SSD dEdX QA task
//Marek Chojnacki
//Marek.Chojnacki@cern.ch
#ifndef ALIANALYSISTASKDEDXSSDQA_H
#define ALIANALYSISTASKDEDXSSDQA_H

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
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

  TH2F*   fHist1;         // CR for each module
  TH2F*   fHist2;         // landau distributions for each module	
  TH3F*   fHist3;         // CR as function of Charge for the AliTrackPoint 
  TH2F*   fHist4;         // Q on chips
  TH2F*   fHist5;         // Q on chips corrected
  TH2F*   fHist6;         // QNvQP not corrected for track inclinaition 
  TList*  fListOfHistos;  // output list	
  Float_t fPcut;          // Momentum cut

  Bool_t fdothecorrection; //do the correction  	
  Float_t fcorrections[20376] ; //[20376]chip corrections


 AliAnalysisTaskdEdxSSDQA(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 AliAnalysisTaskdEdxSSDQA& operator=(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 Int_t Pstrip5(Float_t x,Float_t z) const;
 Int_t Pstrip6(Float_t x,Float_t z) const;
 Int_t Nstrip5(Float_t x,Float_t z) const;
 Int_t Nstrip6(Float_t x,Float_t z) const;
 ClassDef(AliAnalysisTaskdEdxSSDQA, 3); // example of analysis
};

#endif
