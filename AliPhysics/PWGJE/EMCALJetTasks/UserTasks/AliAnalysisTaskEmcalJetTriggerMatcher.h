#ifndef AliAnalysisTaskEmcalJetTriggerMatcher_h
#define AliAnalysisTaskEmcalJetTriggerMatcher_h

class TClonesArray;
class TH1;
class TH2;
class TH3;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TLorentzVector;
class TGraph;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrack;
class AliESDCaloCluster;
class TClonesArray;
class TArrayI;
class TProfile;
class AliParticleContainer;
class AliClusterContainer;

#include <TRef.h>
#include <TBits.h>
#include <TMath.h>

// this whole section of includes added 
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>
#include "AliAnalysisTaskEmcalJet.h"
#include <AliESDCaloCluster.h>
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisFilter.h"

class AliAnalysisTaskEmcalJetTriggerMatcher : public AliAnalysisTaskEmcalJet {
 public:
    
    enum MainPatchType {
        kManual = 0,    //just select highest energy patch in array
        kEmcalJet = 1,   //use functionality of AliAnalysisTaskEmcal
        kEmcalGamma = 2,
        kEmcalJetcal = 3,
        kEmcalGammacal = 4

    };
    
  AliAnalysisTaskEmcalJetTriggerMatcher();
  AliAnalysisTaskEmcalJetTriggerMatcher(const char *name);
  virtual ~AliAnalysisTaskEmcalJetTriggerMatcher();
  
  virtual void           UserCreateOutputObjects();
    
  void SetMainPatchType(MainPatchType t)    { fMainPatchType      = t;}
  void SetMainTriggerTypeCat(TriggerCategory cat, Bool_t b) {fMainTrigCat = cat; fMainTrigSimple = b;}
    
  //JetTrigger Sparse
  virtual THnSparse*      NewTHnSparseDJetTrigger(const char* name, UInt_t entries);
  virtual void            GetDimParamsJetTrigger(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  // setters
  void                          SetJetPt(Double_t jpt)                      { fJetHIpt    = jpt;   }  // jet threshold pt cut
  virtual void                  SetTrackEta(Double_t e)                     { fTrackEta   = e;     }  //eta range of the associated tracks
  virtual void                  SetFillHistograms(Bool_t hists)             { fFillHists  = hists; }  // Fill Histograms
  virtual void                  SetTrigClusterMatchDistance(Double_t dist)  { fMatchDist  = dist;  }  // Jet Cluster to Patch Matching Distance in units of towers (0.014 * fMatchDist)

  // eta and phi limits of jets - setters
  virtual void            SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void            SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }
  
  // event no.
  Int_t event;          // event number (processed)
 protected:
  Bool_t                          Run();
  virtual void                    Terminate(Option_t *);
  virtual Int_t                   AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  void 			                  ExecOnce();
  void                            ExtractMainPatch();
  
  // Trigger bit
  Bool_t                     IsEventEMCALL1Gamma1()            const { return fEventTrigEMCALL1Gamma1  ; }
  Bool_t                     IsEventEMCALL1Gamma2()            const { return fEventTrigEMCALL1Gamma2  ; }
  Bool_t                     IsEventEMCALL1Gamma()             const { return (fEventTrigEMCALL1Gamma1 || fEventTrigEMCALL1Gamma2) ; }
  Bool_t                     fEventTrigEMCALL1Gamma1 ;   // Event is L1-Gamma, threshold 1 on its name,  it should correspond kEMCEGA
  Bool_t                     fEventTrigEMCALL1Gamma2 ;   // Event is L1-Gamma, threshold 2 on its name,  it should correspond kEMCEGA
  
  // data type switch
  AliVEvent                  * fInputEvent;                //! pointer to esd or aod input
  
  // cuts
  Double_t                   fPhimin;
  Double_t                   fPhimax;
  Double_t                   fEtamin;
  Double_t                   fEtamax;
  Double_t                   fAreacut;                     //area cut
  Double_t                   fJetHIpt;                    // high jet pt
  Double_t                   fTrackEta;
  Double_t                   fTrkQAcut;                    //trkQA cut
  Bool_t                     fFillHists;
  Double_t                   fMatchDist;
    
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters

 private:
  AliESDEvent                *fESD;//!  // ESD object
  AliAODEvent                *fAOD;//!  // AOD Object
    
  AliEMCALTriggerPatchInfo  *fMaxPatch;             // main patch
  MainPatchType             fMainPatchType;         // method to select main patch
  TriggerCategory           fMainTrigCat;           // trigger category for main trigger from AliAnalysisTaskEmcal::GetMainTriggerPatch
  Bool_t                    fMainTrigSimple;        // use offline trigger instead of online from AliAnalysisTaskEmcal::GetMainTriggerPatch
  
  THnSparse                 *fhnTriggerInfo;//! patch energy and event observables
  TH1F                      *fHistMatchedClusJet;//!
  TH2F                      *fHistMatchdetadphi;//!
  TH1F                      *fHistTriggerBitInfo;//!
  TH1F                      *fHistMaxTriggerBitInfo;//!
  TH1F                      *fHistEventSelection;//!
  TH2F                      *fHistRecalcGASize;//!
  TH1F                      *fHistRecalcGAEnergy;//!
  TH2F                      *fHistClusEvPatchE;//!      Matched Cluster E vs Patch E
  TH2F                      *fHistdEtaPatchvdPhiPatch;//!
  TH2F                      *fHistRawJetPtvPatchE;//!
  TH2F                      *fHistdEtaPatchvdPhiPatchCMtoGeo;//!  Plot of Trigger Geo pos - Trigger COM pos
  TH1F                      *fHistEventClusSpect;//!
  THnSparse                 *fhnJetTrigger;//!             // Jet Trigger Sparse

  //Declare it private to avoid compilation warning
  AliAnalysisTaskEmcalJetTriggerMatcher(const AliAnalysisTaskEmcalJetTriggerMatcher & g) ; // cpy ctor
  AliAnalysisTaskEmcalJetTriggerMatcher& operator=(const AliAnalysisTaskEmcalJetTriggerMatcher&); // not implemented
  


  ClassDef(AliAnalysisTaskEmcalJetTriggerMatcher, 4); 
};
#endif
