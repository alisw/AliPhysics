#ifndef ALIPRIMARYKAONSELECTOR_H
#define ALIPRIMARYKAONSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliPrimaryKaonCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliPrimaryKaonSelector : public AliAnalysisTaskSE {

 public:

   AliPrimaryKaonSelector(const char *name="KaonSelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliPrimaryKaonSelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected NegKaon/PosKaon array
   vector <Int_t> GetReconstructedNegKaonIndex(){ return fNegKaonsIndex; }
   vector <Int_t> GetReconstructedPosKaonIndex(){ return fPosKaonsIndex; }
   AliPrimaryKaonCuts *GetPrimaryKaonCuts(){   return fKaonCuts; }
   TList *GetCutHistograms(){ if(fKaonCuts){return fKaonCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetPrimaryKaonCuts(const TString cut);
   void SetPrimaryKaonCuts(AliPrimaryKaonCuts *cuts){fKaonCuts=cuts;}

 protected:
   //selected kaon arrays
   
   Bool_t ProcessESDs();
   Bool_t ProcessAODs();
   AliPrimaryKaonCuts *fKaonCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPosKaonsIndex;
   vector<Int_t> fNegKaonsIndex;
   Bool_t fEventIsSelected;

 private:
   AliPrimaryKaonSelector (const AliPrimaryKaonSelector&); // not implemented
   AliPrimaryKaonSelector & operator=(const AliPrimaryKaonSelector&); // not implemented


   
   ClassDef(AliPrimaryKaonSelector,2)
      };

inline void AliPrimaryKaonSelector::SetPrimaryKaonCuts(const TString cut){
   if(fKaonCuts != NULL){
      delete fKaonCuts;
      fKaonCuts=NULL;
   }
   if(fKaonCuts == NULL){
      fKaonCuts=new AliPrimaryKaonCuts("KaonCuts","KaonCuts");
      fKaonCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
