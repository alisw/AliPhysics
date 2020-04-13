#ifndef ALIIDENTIFIEDPRIMARYSELECTOR_H
#define ALIIDENTIFIEDPRIMARYSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliIdentifiedPrimaryCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliIdentifiedPrimarySelector : public AliAnalysisTaskSE {

 public:

   AliIdentifiedPrimarySelector(const char *name="IdentifiedPrimarySelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliIdentifiedPrimarySelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected NegIdentified/PosIdentified array
   vector <Int_t> GetReconstructedNegIdentifiedIndex(){ return fNegIdentifiedIndex; }
   vector <Int_t> GetReconstructedPosIdentifiedIndex(){ return fPosIdentifiedIndex; }
   AliIdentifiedPrimaryCuts *GetIdentifiedPrimaryCuts(){   return fIdentifiedCuts; }
   TList *GetCutHistograms(){ if(fIdentifiedCuts){return fIdentifiedCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetIdentifiedPrimaryCuts(const TString cut);
   void SetIdentifiedPrimaryCuts(AliIdentifiedPrimaryCuts *cuts){fIdentifiedCuts=cuts;}

 protected:
   //selected Identified arrays
   
   Bool_t ProcessESDs();
   Bool_t ProcessAODs();
   AliIdentifiedPrimaryCuts *fIdentifiedCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPosIdentifiedIndex;
   vector<Int_t> fNegIdentifiedIndex;
   Bool_t fEventIsSelected;

 private:
   AliIdentifiedPrimarySelector (const AliIdentifiedPrimarySelector&); // not implemented
   AliIdentifiedPrimarySelector & operator=(const AliIdentifiedPrimarySelector&); // not implemented


   
   ClassDef(AliIdentifiedPrimarySelector,2)
      };

inline void AliIdentifiedPrimarySelector::SetIdentifiedPrimaryCuts(const TString cut){
   if(fIdentifiedCuts != NULL){
      delete fIdentifiedCuts;
      fIdentifiedCuts=NULL;
   }
   if(fIdentifiedCuts == NULL){
      fIdentifiedCuts=new AliIdentifiedPrimaryCuts("IdentifiedCuts","IdentifiedCuts");
      fIdentifiedCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
