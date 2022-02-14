#ifndef ALIPRIMARYDEUTERONSELECTOR_H
#define ALIPRIMARYDEUTERONSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliPrimaryDeuteronCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliPrimaryDeuteronSelector : public AliAnalysisTaskSE {

 public:

   AliPrimaryDeuteronSelector(const char *name="DeuteronSelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliPrimaryDeuteronSelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected electron/positron array
   vector <Int_t> GetReconstructedNegDeuteronIndex(){ return fNegDeuteronsIndex; }
   vector <Int_t> GetReconstructedPosDeuteronIndex(){ return fPosDeuteronsIndex; }
   AliPrimaryDeuteronCuts *GetPrimaryDeuteronCuts(){   return fDeuteronCuts; }
   TList *GetCutHistograms(){ if(fDeuteronCuts){return fDeuteronCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetPrimaryDeuteronCuts(const TString cut);
   void SetPrimaryDeuteronCuts(AliPrimaryDeuteronCuts *cuts){fDeuteronCuts=cuts;}

 protected:
   //selected electron arrays
   
   Bool_t ProcessESDs();
   Bool_t ProcessAODs();
   AliPrimaryDeuteronCuts *fDeuteronCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPosDeuteronsIndex;
   vector<Int_t> fNegDeuteronsIndex;
   Bool_t fEventIsSelected;

 private:
   AliPrimaryDeuteronSelector (const AliPrimaryDeuteronSelector&); // not implemented
   AliPrimaryDeuteronSelector & operator=(const AliPrimaryDeuteronSelector&); // not implemented


   
   ClassDef(AliPrimaryDeuteronSelector,2)
      };

inline void AliPrimaryDeuteronSelector::SetPrimaryDeuteronCuts(const TString cut){
   if(fDeuteronCuts != NULL){
      delete fDeuteronCuts;
      fDeuteronCuts=NULL;
   }
   if(fDeuteronCuts == NULL){
      fDeuteronCuts=new AliPrimaryDeuteronCuts("DeuteronCuts","DeuteronCuts");
      fDeuteronCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
