// Ntuplizer of event data into TTree.

#ifndef ALIANALYSISTASKNTUPLIZER_H
#define ALIANALYSISTASKNTUPLIZER_H

class TH1;
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

 private:

   // dummy (not implemented) copy constructions to avoid compilation warnings
   AliAnalysisTaskNtuplizer(const AliAnalysisTaskNtuplizer&);
   AliAnalysisTaskNtuplizer& operator=(const AliAnalysisTaskNtuplizer&);

   TH1*        fhStats;         //! histogram with event statistics
   TObjArray*  fHistos;         //! array with output histograms
   TTree*      fTree;           //! output tree
   AliEventCuts fEventCut_d;   //!  
   AliAODEvent* vEvent; //!

   ClassDef(AliAnalysisTaskNtuplizer, 1);
};

#endif
