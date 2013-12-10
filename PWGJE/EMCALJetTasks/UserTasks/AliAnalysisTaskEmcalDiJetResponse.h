#ifndef ALIANALYSISTASKEMCALDIJETRESPONSE_H
#define ALIANALYSISTASKEMCALDIJETRESPONSE_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisUtils;
class AliAnalysisManager;
class AliJetContainer;

#include "AliAnalysisTaskEmcalDiJetBase.h"

class AliAnalysisTaskEmcalDiJetResponse : public AliAnalysisTaskEmcalDiJetBase {
 public:
  AliAnalysisTaskEmcalDiJetResponse();
  AliAnalysisTaskEmcalDiJetResponse(const char *name);
  virtual ~AliAnalysisTaskEmcalDiJetResponse();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  //Setters
  void SetMatchFullCharged(Bool_t b)        { fDoMatchFullCharged = b; }
  void SetResponseVar(Int_t v)              { fnUsedResponseVar   = v; }    

  //Getters

 protected:
  Bool_t                      Run()              ;
  void                        CorrelateJets(const Int_t type);
  void                        CorrelateAllJets(const Int_t type);
  void                        CorrelateTwoJets(const Int_t type);

  void                        FillDiJetHistos(const AliEmcalJet *jet1 = 0, const AliEmcalJet *jet2 = 0, const Int_t mode = 0);
  void                        FillMatchHistos();
  Bool_t                      RetrieveEventObjects();

  void                        FillDiJetResponse(const AliEmcalJet *jetTrigMC = 0, const AliEmcalJet *jetAssocMC = 0, const AliEmcalJet *jetTrigDet = 0, const AliEmcalJet *jetAssocDet = 0, Int_t type = 0);

 private:
  Bool_t            fDoMatchFullCharged;          //  do full-charged matching histos
  THnSparse        *fhnDiJetResponseCharged;      //! sparse with di-jet properties (full-full)
  THnSparse        *fhnDiJetResponseFullCharged;  //! sparse with di-jet properties (full-full)
  TH1F             *fh1TriggersCharged[2];        //! charged jet triggers
  TH1F             *fh1TriggersFull[2];           //! full jet triggers
  TH1F             *fh1TriggersLostCharged;       //! lost charged jet triggers
  TH1F             *fh1TriggersLostFull;          //! lost full jet triggers
  TH3F             *fh3AssocLostPtDeltaPhiCharged;//! lost charged associated jet
  TH3F             *fh3AssocLostPtDeltaPhiFull;   //! lost full associated jet
  THnSparse        *fhnMatchingCharged;           //! sparse comparing matched particle and detector level charged jets
  THnSparse        *fhnMatchingFull;              //! sparse comparing matched particle and detector level charged jets
  Int_t             fnUsedResponseVar;            //  build response for kt (0) or dijet eta (1)


  AliAnalysisTaskEmcalDiJetResponse(const AliAnalysisTaskEmcalDiJetResponse&);            // not implemented
  AliAnalysisTaskEmcalDiJetResponse &operator=(const AliAnalysisTaskEmcalDiJetResponse&); // not implemented

  ClassDef(AliAnalysisTaskEmcalDiJetResponse, 2) // jet sample analysis task
};
#endif
