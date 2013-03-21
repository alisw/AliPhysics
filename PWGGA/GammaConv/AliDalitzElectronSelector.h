#ifndef ALIDALITZELECTRONSELECTOR_H
#define ALIDALITZELECTRONSELECTOR_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliDalitzElectronCuts.h"
#include "AliExternalTrackParam.h"

class TRandom3;
class AliStack;
class TList;
class TString;
class TH1F;
class TH2F;

using namespace std;

class AliDalitzElectronSelector : public AliAnalysisTaskSE {

 public:

   AliDalitzElectronSelector(const char *name="ElectronSelector");
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliDalitzElectronSelector();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);
   virtual void Init();

   Bool_t ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent=NULL);
   Bool_t IsEventSelected(){return fEventIsSelected;}

   // Return selected electron/positron array
   vector <Int_t> GetReconstructedElectronsIndex(){ return fElectronsIndex; }
   vector <Int_t> GetReconstructedPositronsIndex(){ return fPositronsIndex; }
   AliDalitzElectronCuts *GetDalitzElectronCuts(){   return fElectronCuts; }
   TList *GetCutHistograms(){ if(fElectronCuts){return fElectronCuts->GetCutHistograms();} return NULL;}
   // Set Options

   void SetDalitzElectronCuts(const TString cut);
   void SetDalitzElectronCuts(AliDalitzElectronCuts *cuts){fElectronCuts=cuts;}

 protected:
   //selected electron arrays
   
   Bool_t ProcessESDs();
   AliDalitzElectronCuts *fElectronCuts; // Pointer to the ConversionCut Selection
   vector<Int_t> fPositronsIndex;
   vector<Int_t> fElectronsIndex;
   Bool_t fEventIsSelected;

 private:
   AliDalitzElectronSelector (const AliDalitzElectronSelector&); // not implemented
   AliDalitzElectronSelector & operator=(const AliDalitzElectronSelector&); // not implemented


   
   ClassDef(AliDalitzElectronSelector,1)
      };

inline void AliDalitzElectronSelector::SetDalitzElectronCuts(const TString cut){
   if(fElectronCuts != NULL){
      delete fElectronCuts;
      fElectronCuts=NULL;
   }
   if(fElectronCuts == NULL){
      fElectronCuts=new AliDalitzElectronCuts("ElectronCuts","ElectronCuts");
      fElectronCuts->InitializeCutsFromCutString(cut.Data());
   }
}


#endif
