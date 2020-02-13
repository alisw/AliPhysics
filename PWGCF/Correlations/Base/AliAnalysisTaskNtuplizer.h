// Ntuplizer of event data into TTree.

#ifndef AliAnalysisTaskNtuplizer_H
#define AliAnalysisTaskNtuplizer_H

class TH1;
class TH2D;
class TTree;
class TObjArray;
class AliEventCuts;

#include <AliAnalysisTaskSE.h>
#include <AliEventCuts.h>

class AliAnalysisTaskNtuplizer : public AliAnalysisTaskSE {

 public:

   AliAnalysisTaskNtuplizer();
   AliAnalysisTaskNtuplizer(const char* name);
   ~AliAnalysisTaskNtuplizer();

   void UserCreateOutputObjects();
   void UserExec(Option_t*);

   void WarnIncErr(const char* msg);

   void SetFilterBit(Int_t FilterBit) {fFilterBit=FilterBit;}
 void SetKinematicCuts_chtracks(Float_t minPt, Float_t maxPt,Float_t eta)
  {
    fminPt=minPt;
    fmaxPt=maxPt;
    feta=eta;
  }
 void SetFillchtracks(Bool_t Fillchtracks) {fFillchtracks=Fillchtracks;}
 void SetFilllambdas(Bool_t Filllambdas) {fFilllambdas=Filllambdas;}
 void SetRunPeriod(Int_t ds) {year = ds;}
 void SetFillSPD(Bool_t fillSPD) {fFillSPD = fillSPD;}

 private:

   // dummy (not implemented) copy constructions to avoid compilation warnings
   AliAnalysisTaskNtuplizer(const AliAnalysisTaskNtuplizer&);
   AliAnalysisTaskNtuplizer& operator=(const AliAnalysisTaskNtuplizer&);


  Bool_t fFillchtracks;
  Bool_t fFilllambdas;
  Int_t fFilterBit;
  Float_t fminPt;
  Float_t fmaxPt;
  Float_t feta;
  Int_t year;
  Bool_t fFillSPD;

     
   TH1*        fhStats;         //! histogram with event statistics
   TH2D* eta_phi_accpt;   //!
   TObjArray*  fHistos;         //! array with output histograms
   TTree*      fTree;           //! output tree
   AliEventCuts fEventCut_d;   //!  
   AliAODEvent* vEvent; //!

   ClassDef(AliAnalysisTaskNtuplizer, 1);
};

#endif

