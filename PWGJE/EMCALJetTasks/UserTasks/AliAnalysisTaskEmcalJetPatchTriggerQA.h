#ifndef AliAnalysisTaskEmcalJetPatchTriggerQA_h
#define AliAnalysisTaskEmcalJetPatchTriggerQA_h

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class AliLocalRhoParameter;

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

class AliAnalysisTaskEmcalJetPatchTriggerQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetPatchTriggerQA();
  AliAnalysisTaskEmcalJetPatchTriggerQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetPatchTriggerQA();  
  
  virtual void           UserCreateOutputObjects();
  virtual THnSparse*     NewTHnSparseF(const char* name, UInt_t entries);
  virtual void           GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual void           SetLocalRhoName(const char *n)        { fLocalRhoName  = n; }

  virtual void           SetAreaCut(Double_t a)                { fAreacut    = a; }
  virtual void           SetJetEta(Double_t emin, Double_t emax)  { fEtamin = emin; fEtamax = emax; }
  virtual void           SetJetPhi(Double_t pmin, Double_t pmax)  { fPhimin = pmin; fPhimax = pmax; }

  // getters
  TString		 GetLocalRhoName() const		 {return fLocalRhoName; } 

 protected:
  Bool_t                 Run();
  virtual void           Terminate(Option_t *);
  virtual Int_t          AcceptMyJet(AliEmcalJet *jet);   // applies basic jet tests/cuts before accepting
  virtual Int_t          GetCentBin(Double_t cent) const;
  Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const;

  // parameters of detector to cut on for event
  Double_t               fPhimin;                  // phi min
  Double_t               fPhimax;                  // phi max
  Double_t               fEtamin;                  // eta min
  Double_t               fEtamax;                  // eta max
  Double_t               fAreacut;                 // area cut

  void 			ExecOnce();
  Double_t		fLocalRhoVal;

 private:
  TH2F                  *fHistNjetvsCent;          //!number of jets versus Centrality
  THnSparse             *fhnJetTriggerQA;      //! jet sparse


  AliAnalysisTaskEmcalJetPatchTriggerQA(const AliAnalysisTaskEmcalJetPatchTriggerQA&);            // not implemented
  AliAnalysisTaskEmcalJetPatchTriggerQA& operator=(const AliAnalysisTaskEmcalJetPatchTriggerQA&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetPatchTriggerQA, 4); // ChristineQA
};
#endif
