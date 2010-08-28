#ifndef AliAnalysisTaskTotEt_cxx
#define AliAnalysisTaskTotEt_cxx

class AliAnalysisEt;
class TTree;
class AliVParticle;
class TH1F;
class TH2F;
class TNtuple;
class TObjArray;
class AliESDEvent;
class AliMCParticle;
class TDatabasePDG;

#include "AliAnalysisTaskSE.h"
#include "TObject.h"

//class ParticleVars : public TObject        // Inherit from TObject to put in TClonesArray
//{
//public:
//  
//  ParticleVars() : TObject(){}
//  Int_t fPdgCode; // from MC
//  Int_t fPid; //from ESDs
//  Int_t fMass;
//  Int_t fCharge;
//  Double_t fEt;
//  Double_t fPhi;
//  Double_t fEta;
//  
//  ClassDef(ParticleVars, 1);
//  
//};
//
class AliAnalysisTaskTotEt : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskTotEt(const char *name = "AliAnalysisTaskTotEt");
  virtual ~AliAnalysisTaskTotEt() {}
private:
  //Declare it private to avoid compilation warning
  AliAnalysisTaskTotEt & operator = (const AliAnalysisTaskTotEt & g) ;//cpy assignment
  AliAnalysisTaskTotEt(const AliAnalysisTaskTotEt & g) ; // cpy ctor
  
public:
  
  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  virtual void SetTriggerSelection(Bool_t v) {
    fTriggerSelection = v;
  }
  
  /* // Not yet implemented methods commented out for now..
   private:
   
   Float_t CorrectForCaloAcceptance(Float_t energy);
   bool CheckGoodVertex(AliVParticle *track);
   bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);
   bool ParticleInCalorimeter(AliMCParticle *part);
   */
  
private:
  
  AliESDEvent *fESD;    //ESD object
  
  TList *fOutputList;
  
  AliAnalysisEt *fRecAnalysis;
  AliAnalysisEt *fMCAnalysis;
  
  TH2F *fHistEtRecvsEtMC;
  
  Bool_t fTriggerSelection;
  
  Int_t fCount;
  
  const int fkPhotonPdg;
  
  const Float_t fkProtonMass;
  
  TDatabasePDG *fPdgDB;
  
  class EventVars
  {
  public:
    Double_t fTotEt;
    Double_t fTotEtAcc;
    Double_t fTotEnergy;
    
    Double_t fTotNeutralEt;
    Double_t fTotNeutralEtAcc;
    
    Double_t fTotChargedEt;
    Double_t fTotChargedEtAcc;
    
    Int_t fChargedMultiplicity;
    Int_t fNeutralMultiplicity;
    
  };
  
  EventVars *fRecEventVars;
  EventVars *fSimEventVars;
  
  
  TClonesArray *fRecParticleArray;
  TClonesArray *fSimParticleArray;
  
  ClassDef(AliAnalysisTaskTotEt, 1); // example of analysis
};

#endif
