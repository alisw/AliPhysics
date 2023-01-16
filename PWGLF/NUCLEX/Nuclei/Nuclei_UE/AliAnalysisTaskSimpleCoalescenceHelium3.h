#ifndef AliAnalysisTaskSimpleCoalescenceHelium3_cxx
#define AliAnalysisTaskSimpleCoalescenceHelium3_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "TVector3.h"
#include "TList.h"
#include "TH1D.h"
#include "TF1.h"

//==================== SIMPLE COALESCENCE ====================//
//                                                            //
//                                                            //
//                                                            //
//                                                            //
//============================================================//

//__________________________________________________________________________________________________________________
class AliAnalysisTaskSimpleCoalescenceHelium3 : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskSimpleCoalescenceHelium3();
    AliAnalysisTaskSimpleCoalescenceHelium3(const char *name);
    virtual ~AliAnalysisTaskSimpleCoalescenceHelium3();
        
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //Proton Weights
    void SetProtonWeights (TF1 *fWeight) { fProtonWeights = fWeight; }

    //User Functions
    Bool_t    GetEvent ();
    Bool_t    IsINELgtZERO ();
    Bool_t    IsInjectedParticle    (AliMCParticle *particle);
    Double_t  GetProtonWeight       (Double_t pt);
    Double_t  GetHelium3Weight      (Double_t pt_prot1, Double_t pt_prot2, Double_t pt_neut);
    Bool_t    ThreeBodyCoalescence  (Double_t deltaP1, Double_t deltaP2, Double_t deltaP3, Double_t p0_pp, Double_t p0_pn);
    TLorentzVector LorentzTransform (TLorentzVector R, TVector3 beta_vect);

private:
    AliAODEvent *fAODevent;//!
    AliMCEvent  *fMCEvent;//!
    TList       *fOutputList;//!
    TList       *fQAList;//!
  
    //Proton Weights
    TF1 *fProtonWeights;//
    
    //Event Counter
    TH1D *hNumberOfEvents;//!

    //p_{T} Spectra: Protons
    TH1D *hProtons;//!
    TH1D *hProtons_reshaped;//!
    
    //p_{T} Spectra: Neutrons
    TH1D *hNeutrons;//!
    TH1D *hNeutrons_reshaped;//!
   
    //p_{T} Spectra: Helium3
    TH1D *hHelium3[100];//!
  
    //DeltaP
    TH1D *hDeltaP_pn;//!
    TH1D *hDeltaP_pp;//!

    
    AliAnalysisTaskSimpleCoalescenceHelium3(const AliAnalysisTaskSimpleCoalescenceHelium3&);
    AliAnalysisTaskSimpleCoalescenceHelium3& operator=(const AliAnalysisTaskSimpleCoalescenceHelium3&);
        
    ClassDef(AliAnalysisTaskSimpleCoalescenceHelium3, 1);
    
};
//__________________________________________________________________________________________________________________

#endif
