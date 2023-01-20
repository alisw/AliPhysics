#ifndef AliAnalysisTaskSimpleCoalescenceDeuteronInJets_cxx
#define AliAnalysisTaskSimpleCoalescenceDeuteronInJets_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "TVector3.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TF1.h"

//============================ DEUTERONS IN JET  ============================//
//                                                                           //
//    Deuteron production inside jets using anti-kt jet finder algorithm.    //
//    Coalescence parameter B2 vs. jet radius R and multiplicity.            //
//                                                                           //
//===========================================================================//

//__________________________________________________________________________________________________________________
class AliAnalysisTaskSimpleCoalescenceDeuteronInJets : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskSimpleCoalescenceDeuteronInJets();
    AliAnalysisTaskSimpleCoalescenceDeuteronInJets(const char *name);
    virtual ~AliAnalysisTaskSimpleCoalescenceDeuteronInJets();
        
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //External Settings
    void SetProtonWeights (TF1 *fWeight)     { fProtonWeights = fWeight; }
    void SetJetRadius     (Double_t Radius)  { fJetRadius     = Radius;  }
    void SetMaximumPt     (Double_t ptMax)   { fMaximumPt     = ptMax;   }

    //User Functions
    Bool_t    GetEvent ();
    Bool_t    IsINELgtZERO ();
    Bool_t    IsInjectedParticle    (AliMCParticle *particle);
    Double_t  GetProtonWeight       (Double_t pt);
    Double_t  GetDeuteronWeight     (Double_t pt_prot, Double_t pt_neut);
    Double_t  Minimum               (Double_t x1, Double_t x2);
    Bool_t    TwoBodyCoalescence    (Double_t deltaP, Double_t p0);
    TLorentzVector LorentzTransform (TLorentzVector R, TVector3 beta_vect);

private:
    AliAODEvent *fAODevent;//!
    AliMCEvent  *fMCEvent;//!
    TList       *fOutputList;//!
    TList       *fQAList;//!
    Double_t     fMaximumPt;//
    Double_t     fJetRadius;//

    //Proton Weights
    TF1 *fProtonWeights;//
    
    //Event Counter
    TH1D *hNumberOfEvents;//!

    //p_{T} Spectra
    TH1D *hProtons;//!
    TH1D *hNeutrons;//!
    TH1D *hDeuterons[4];//!
    TH1D *hDeuterons_ptoverA[4];//!

    //General Histograms
    TH1I *hNumberOfParticlesInJet;//!
    TH1D *hMultiplicityBin;//!
    
    
    AliAnalysisTaskSimpleCoalescenceDeuteronInJets(const AliAnalysisTaskSimpleCoalescenceDeuteronInJets&);
    AliAnalysisTaskSimpleCoalescenceDeuteronInJets& operator=(const AliAnalysisTaskSimpleCoalescenceDeuteronInJets&);
        
    ClassDef(AliAnalysisTaskSimpleCoalescenceDeuteronInJets, 1);
    
};
//__________________________________________________________________________________________________________________

#endif
