#ifndef ALIANALYSISTASKJETCLUSTERKINE_H
#define ALIANALYSISTASKJETCLUSTERKINE_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// *******************************************
// Jet Finding in the Kine train
// *******************************************

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h" // cannot forward declare ThnSparseF
#include "AliVParticle.h" //FK//
#ifndef __CINT__
# include "fastjet/ClusterSequenceArea.hh"
# include "fastjet/AreaDefinition.hh"
# include "fastjet/JetDefinition.hh"
#else
namespace fastjet {
  enum JetAlgorithm;
  enum Strategy;
  enum RecombinationScheme;
  enum AreaType;
}
#endif

////////////////
class AliJetHeader;
class AliESDEvent;
class AliAODEvent;
class AliAODExtension;
class AliAODJet;
class AliGenPythiaEventHeader;
class AliJetFinder;
class AliAODMCParticle;
class AliMCEvent;    //FK//
class AliMCEventHandler; //FK//
class AliVParticle; //FK//
class AliGenEventHeader; //FK//
class TList;
class TChain;
class TH2F;
class TH1F;
class TH3F;
class TRefArray;
class TClonesArray;
class TF1;
class TProfile;

class AliAnalysisTaskJetClusterKine : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetClusterKine();
    AliAnalysisTaskJetClusterKine(const char* name);
    virtual ~AliAnalysisTaskJetClusterKine();
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();  //FK//
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();
    
    
    virtual void SetTrackEtaWindow(Float_t f){fTrackEtaWindow = f;}
    virtual void SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
    virtual void SetTrackPtCut(Float_t x){fTrackPtCut = x;}
    virtual void SetVtxCuts(Float_t z){fVtxZCut = z;}

    virtual void SetJetOutputBranch(const char *c){fNonStdBranch = c;}
    virtual void SetJetOutputContainer(Int_t c){fOutContainer = c;} //FF//
    virtual const char* GetJetOutputBranch(){return fNonStdBranch.Data();}
    virtual void SetJetOutputFile(const char *c){fNonStdFile = c;}
    virtual const char* GetJetOutputFile(){return fNonStdFile.Data();}
    virtual void SetMaxTrackPtInJet(Float_t x){fMaxTrackPtInJet = x;}
    virtual void SetJetOutputMinPt(Float_t x){fJetOutputMinPt = x;}



    // for Fast Jet
    fastjet::JetAlgorithm        GetAlgorithm()         const {return fAlgorithm;}
    fastjet::Strategy            GetStrategy()          const {return fStrategy;}
    fastjet::RecombinationScheme GetRecombScheme()      const {return fRecombScheme;}
    fastjet::AreaType            GetAreaType()          const {return fAreaType;}
    // Setters
    void SetRparam(Double_t f)                           {fRparam = f;}
    // Temporary change to integer; problem with dictionary generation?
    //void SetAlgorithm(fastjet::JetAlgorithm f)           {fAlgorithm = f;}
    void SetAlgorithm(Int_t f)                           {fAlgorithm = (fastjet::JetAlgorithm) f;}
    void SetStrategy(fastjet::Strategy f)                {fStrategy = f;}
    void SetRecombScheme(fastjet::RecombinationScheme f) {fRecombScheme = f;}
    void SetAreaType(fastjet::AreaType f)                {fAreaType = f;}
    void SetGhostArea(Double_t f)                        {fGhostArea = f;}
    void SetActiveAreaRepeats(Int_t f)                   {fActiveAreaRepeats = f;}
    void SetGhostEtamax(Double_t f)                      {fGhostEtamax = f;}



    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackKineAll=0, kTrackKineCharged};
    enum { kNoOutput=0, kAODBranch=1, kExchCont=2 }; //FF//
    

 private:

    AliAnalysisTaskJetClusterKine(const AliAnalysisTaskJetClusterKine&);
    AliAnalysisTaskJetClusterKine& operator=(const AliAnalysisTaskJetClusterKine&);

    Int_t GetListOfTracks(TList *list,Int_t type);
	
    AliMCEvent*              fMcEvent;    //! MC event                       //FK//   FFF
    AliInputEventHandler*    fMcHandler;  //! MCEventHandler                 //FK//   FFF
    TRefArray       *fRef;                // ! trefarray for track references within the jet  FFF
    Int_t         fTrackTypeGen;          // type of tracks used for 0=full, 1=charged         FFF 
    Float_t       fAvgTrials;             // Average nimber of trials
    Float_t       fTrackEtaWindow;        // eta window used for corraltion plots between rec and gen  FFF 
    Float_t       fTrackPtCut;            // minimum track pt to be accepted    FFF
    Float_t       fJetOutputMinPt;        // minimum p_t for jets to be written out  FFF
    Float_t       fMaxTrackPtInJet;       // maximum track pt within a jet for flagging... FFF
    Float_t       fVtxZCut;               // zvtx cut

    // output configurartion
    TString       fNonStdBranch;      // the name of the non-std branch name, if empty no branch is filled  FFF
    Int_t         fOutContainer;    //FF//output container 1=AOD Branch   2=Exchange container 
    TString       fNonStdFile;        // The optional name of the output file the non-std branch is written to FFF



    // Fast jet
    Double_t fRparam;                  // fastjet distance parameter  FFF
    fastjet::JetAlgorithm fAlgorithm; //fastjet::kt_algorithm  FFF
    fastjet::Strategy fStrategy;  //= fastjet::Best;    FFF
    fastjet::RecombinationScheme fRecombScheme; // = fastjet::BIpt_scheme;   FFF
    fastjet::AreaType fAreaType;  // fastjet area type
    Double_t fGhostArea;          // fasjet ghost area               FFF
    Int_t fActiveAreaRepeats;     // fast jet active area repeats   FFF
    Double_t fGhostEtamax;        // fast jet ghost area     FFF

    TClonesArray  *fTCAJetsOut;   //! TCA of output jets     FFF

    TProfile*     fh1Xsec;   //! pythia cross section and trials
    TH1F*         fh1Trials; //! trials are added                         FFF
    TH1F*         fh1PtHard;  //! Pt har of the event...       
    TH1F*         fh1PtHardNoW;  //! Pt har of the event without weigt       
    TH1F*         fh1PtHardTrials;  //! Number of trials 

    TH1F*         fh1NJetsGen; //! number of generator level jets            FFF
    TH1F*         fh1NConstGen;//! number of constiutens in leading jet   FFF
    TH1F*         fh1NConstLeadingGen;//! number of constiutens in leading jet  FFF
    TH1F*         fh1PtJetsGenIn;  //! Jet pt for all jets               FFF
    TH1F*         fh1PtJetsLeadingGenIn;  //! Jet pt for the leading jets   FFF
    TH1F*         fh1PtJetConstGen;//! pt of constituents               FFF
    TH1F*         fh1PtJetConstLeadingGen;// pt of constituents    FFF
    TH1F*         fh1PtTracksGenIn;  //! track pt for all tracks  FFF

    TH1F*         fh1Nch;            //! (charged) particle mult      FFF
    TH1F*         fh1Z;                // ! centrality of anaylsed events   FFF 


    TH2F*         fh2NConstPt;           //! number of constituents vs. pt  FFF
    TH2F*         fh2NConstLeadingPt;           //! number of constituents vs. pt FFF
    TH2F*         fh2JetPhiEta;             //! jet phi eta FFF
    TH2F*         fh2LeadingJetPhiEta;      //! leading jet phi eta  FFF
    TH2F*         fh2JetEtaPt;              //! leading jet eta FFF
    TH2F*         fh2LeadingJetEtaPt;              //! leading jet eta  FFF
    TH2F*         fh2TrackEtaPt;              //! track eta for all tracks  FFF
    TH2F*         fh2JetsLeadingPhiEta;     //! jet delta phi delta eta w.r.t. leading jet  FFF
    TH2F*         fh2JetsLeadingPhiPt;      //! jet correlation with leading jet FFF
    TH2F*         fh2JetsLeadingPhiPtW;      //! jet correlation with leading jet FFF

    TList *fHistList; //!leading tracks to be skipped in the randomized event Output list
   

    ClassDef(AliAnalysisTaskJetClusterKine, 1) 
};
 
#endif
