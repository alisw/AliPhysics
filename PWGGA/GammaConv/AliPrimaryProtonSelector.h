#ifndef ALIPRIMARYPROTONSELECTOR_H
#define ALIPRIMARYPROTONSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliPrimaryProtonCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliPrimaryProtonSelector : public AliAnalysisTaskSE {

 public:

   AliPrimaryProtonSelector(const char *name="ProtonSelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliPrimaryProtonSelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected proton/antiproton array
   vector <Int_t> GetReconstructedNegProtonIndex(){ return fNegProtonsIndex; }
   vector <Int_t> GetReconstructedPosProtonIndex(){ return fPosProtonsIndex; }
   AliPrimaryProtonCuts *GetPrimaryProtonCuts(){   return fProtonCuts; }
   TList *GetCutHistograms(){ if(fProtonCuts){return fProtonCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetPrimaryProtonCuts(const TString cut);
   void SetPrimaryProtonCuts(AliPrimaryProtonCuts *cuts){fProtonCuts=cuts;}

 protected:
   //selected electron arrays
   
   Bool_t ProcessESDs();
   Bool_t ProcessAODs();
   AliPrimaryProtonCuts *fProtonCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPosProtonsIndex;
   vector<Int_t> fNegProtonsIndex;
   Bool_t fEventIsSelected;

 private:
   AliPrimaryProtonSelector (const AliPrimaryProtonSelector&); // not implemented
   AliPrimaryProtonSelector & operator=(const AliPrimaryProtonSelector&); // not implemented


   
   ClassDef(AliPrimaryProtonSelector,2)
      };

inline void AliPrimaryProtonSelector::SetPrimaryProtonCuts(const TString cut){
   if(fProtonCuts != NULL){
      delete fProtonCuts;
      fProtonCuts=NULL;
   }
   if(fProtonCuts == NULL){
      fProtonCuts=new AliPrimaryProtonCuts("ProtonCuts","ProtonCuts");
      fProtonCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
