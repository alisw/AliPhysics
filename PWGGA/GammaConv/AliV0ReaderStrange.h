#ifndef ALIV0READERSTRANGE_H
#define ALIV0READERSTRANGE_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliConversionPhotonCuts.h"
#include "AliV0CutsStrange.h"
#include "AliConvEventCuts.h"
#include "AliExternalTrackParam.h"
#include "TObject.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include <vector>
#include "AliESDpid.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliAnalysisManager.h"
#include "AliV0ParticleStrange.h"


class TRandom3;
class TList;
class AliKFConversionPhoton;
class TString;
class TClonesArray;
class TH1F;
class TH2F;
class AliAODConversionPhoton;

using namespace std;

class AliV0ReaderStrange : public AliAnalysisTaskSE {

  public:
    AliV0ReaderStrange(const char *name="V0ReaderStrange");
    virtual                    ~AliV0ReaderStrange();                            //virtual destructor

    void                      UserCreateOutputObjects();
    virtual void              UserExec(Option_t *option);
    virtual void              Init();

    Bool_t                    ProcessEvent (AliVEvent *inputEvent, AliMCEvent *mcEvent=NULL);
    
    AliV0ParticleStrange      *ReconstructV0(AliESDEvent *fESDEvent, AliESDv0 *fCurrentV0,Int_t currentV0Index);
    AliV0ParticleStrange      *ReconstructV0(AliAODEvent *fAODEvent, AliAODv0 *fCurrentV0,Int_t currentV0Index);
    
    // Return Reconstructed Gammas
    TClonesArray*             GetReconstructedV0s()              {return fConversionGammas;}
    Int_t                     GetNReconstructedV0s()             {if(fConversionGammas){return fConversionGammas->GetEntriesFast();} else{ return 0;}}

    AliV0CutsStrange*         GetV0Cuts()                           {return fV0Cuts;}
    AliConvEventCuts*         GetEventCuts()                        {return fEventCuts;}
//     TList*                    GetCutHistograms()                    {if(fConversionCuts) {return fConversionCuts->GetCutHistograms();}
//                                                                      return NULL;}
    TList*                    GetEventCutHistograms()               {if(fEventCuts) {return fEventCuts->GetCutHistograms();}
                                                                     return NULL;}

    // Set Options
    void               SetEventCuts(const TString cut);
    void               SetEventCuts(AliConvEventCuts *cuts)             {fEventCuts=cuts; return;}

    void               SetV0Cuts(const TString cut);
    void               SetV0Cuts(AliV0CutsStrange *cuts)                       {fV0Cuts=cuts; return;}
    
  protected:
    const AliExternalTrackParam*   GetExternalTrackParam(AliESDv0 *fCurrentV0, Int_t &tracklabel, Int_t charge);
    
    Bool_t                  ProcessESDV0s();
    Bool_t                  ProcessAODV0s();
    Bool_t                  GetAODConversionGammas();

    AliConvEventCuts         *fEventCuts;         // Pointer to the Event Cut Selection
    AliV0CutsStrange                *fV0Cuts;
    TClonesArray             *fConversionGammas;  // TClonesArray holding the reconstructed photons
   
   Bool_t         fEventIsSelected;
   
    vector<Int_t>  fVectorFoundGammas;            // vector with found MC labels of gammas

  private:
    AliV0ReaderStrange(AliV0ReaderStrange &original);
    AliV0ReaderStrange &operator=(const AliV0ReaderStrange &ref);

    ClassDef(AliV0ReaderStrange, 1)

};

inline void AliV0ReaderStrange::SetEventCuts(const TString cut){
  if(fEventCuts != NULL){
  delete fEventCuts;
  fEventCuts=NULL;
  }
  if(fEventCuts == NULL){
    fEventCuts=new AliConvEventCuts("V0ReaderEventCuts","V0ReaderEventCuts");
    fEventCuts->InitializeCutsFromCutString(cut.Data());
  }
}

inline void AliV0ReaderStrange::SetV0Cuts(const TString cut){
  if(fV0Cuts != NULL){
  delete fV0Cuts;
  fV0Cuts=NULL;
  }
  if(fV0Cuts == NULL){
    fV0Cuts=new AliV0CutsStrange("V0ReaderCutsStrange","V0ReaderCutsStrange");
    fV0Cuts->InitializeCutsFromCutString(cut.Data());
  }
}



#endif
