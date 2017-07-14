#ifndef ALIPRIMARYPIONSELECTOR_H
#define ALIPRIMARYPIONSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliPrimaryPionCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliPrimaryPionSelector : public AliAnalysisTaskSE {

 public:

   AliPrimaryPionSelector(const char *name="PionSelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliPrimaryPionSelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected electron/positron array
   vector <Int_t> GetReconstructedNegPionIndex(){ return fNegPionsIndex; }
   vector <Int_t> GetReconstructedPosPionIndex(){ return fPosPionsIndex; }
   AliPrimaryPionCuts *GetPrimaryPionCuts(){   return fPionCuts; }
   TList *GetCutHistograms(){ if(fPionCuts){return fPionCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetPrimaryPionCuts(const TString cut);
   void SetPrimaryPionCuts(AliPrimaryPionCuts *cuts){fPionCuts=cuts;}

 protected:
   //selected electron arrays
   
   Bool_t ProcessESDs();
   AliPrimaryPionCuts *fPionCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPosPionsIndex;
   vector<Int_t> fNegPionsIndex;
   Bool_t fEventIsSelected;

 private:
   AliPrimaryPionSelector (const AliPrimaryPionSelector&); // not implemented
   AliPrimaryPionSelector & operator=(const AliPrimaryPionSelector&); // not implemented


   
   ClassDef(AliPrimaryPionSelector,1)
      };

inline void AliPrimaryPionSelector::SetPrimaryPionCuts(const TString cut){
   if(fPionCuts != NULL){
      delete fPionCuts;
      fPionCuts=NULL;
   }
   if(fPionCuts == NULL){
      fPionCuts=new AliPrimaryPionCuts("ElectronCuts","ElectronCuts");
      fPionCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
