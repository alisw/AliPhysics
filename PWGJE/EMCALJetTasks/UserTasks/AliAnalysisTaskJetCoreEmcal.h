#ifndef ALIANALYSISTASKJETCOREEMCAL_H
#define ALIANALYSISTASKJETCOREEMCAL_H

//
// Implementation of jet core task in Emcal framework
// Basic functionality copied from AliAnalysisTaskEmcalJetSample
//

class TH1F;
class TH1I;
class TH2F;
class TH3F;
class THnSparse;
class TRandom3;

#include "AliEventCuts.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"


class AliAnalysisTaskJetCoreEmcal : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetCoreEmcal()                                               ;
  AliAnalysisTaskJetCoreEmcal(const char *name)                               ;
  virtual ~AliAnalysisTaskJetCoreEmcal()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;


	// Setters
	virtual void		 SetJetShapeType(Int_t type){fJetShapeType=type;}
	virtual void		 SetCentrality(Float_t cmin, Float_t cmax){fCentMin=cmin; fCentMax=cmax;}
	virtual void     SetTTLowRef(Float_t ttlow){fTTLowRef=ttlow;}
	virtual void     SetTTUpRef(Float_t ttup){fTTUpRef=ttup;}
	virtual void     SetTTLowSig(Float_t ttlows){fTTLowSig=ttlows;}
	virtual void     SetTTUpSig(Float_t ttups){fTTUpSig=ttups;}
	virtual void		 SetNRPBins(Float_t nrpb){fNRPBins=nrpb;}
	virtual void		 SetSignalFraction(Float_t sfrac){fFrac=sfrac;}
	virtual void     SetJetHadronDeltaPhi(Float_t delta){fJetHadronDeltaPhi=delta;}
	virtual void     SetMinFractionSharedPt(Float_t min){fMinFractionSharedPt=min;}
	virtual void     SetMinEmbJetPt(Float_t min){fMinEmbJetPt=min;}
	virtual void		 SetJetContName(TString cont){fJetContName=cont;}
	virtual void		 SetJetContTrueName(TString cont){fJetContTrueName=cont;}
	virtual void		 SetJetContPartName(TString cont){fJetContPartName=cont;}
	virtual void		 SetFillTrackHistograms(Bool_t b){fFillTrackHistograms=b;}
	virtual void		 SetFillJetHistograms(Bool_t b){fFillJetHistograms=b;}
	virtual void		 SetFillRecoilTHnSparse(Bool_t b){fFillRecoilTHnSparse=b;}
	virtual void		 SetFillInclusiveTree(Bool_t b){fFillInclusiveTree=b;}
	virtual void		 SetFillRecoilTree(Bool_t b){fFillRecoilTree=b;}
	virtual void		 SetPtHardBin(Int_t bin){fPtHardBin=bin;}
	virtual void		 SetRejectionFactorInclusiveJets(Int_t f){fRejectionFactorInclusiveJets=f;}
	virtual void     SetMoreTreeVars(Bool_t more){fMoreTreeVars=more;} 

//  static AliAnalysisTaskJetCoreEmcal* AddTaskJetCoreEmcal( //      const char *ntracks            = "usedefault", //      const char *nclusters          = "usedefault", //      const char* ncells             = "usedefault",
//      const char *suffix             = "");

  enum JetShapeType {
    kMCTrue = 0,   // generated jets only
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
    kDetEmbPart = 3, // embedding
    kDetEmbPartCorr = 4, // embedding, do embedded h+jet correlation
    kDetPart = 5, // pp response
    kDetEmbDet = 6 // pp data embedding
//    kDetEmb = 3,  //detector embedded jets
//    kPythiaDef = 5,
//    kDetEmbPartPythia=6,
//    kGenOnTheFly = 7
  };

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateCellHistograms()                          ;
  void                        AllocateJetCoreHistograms()                          ;

  void                        DoJetLoop()                                       ;
  void                        DoTrackLoop()                                     ;
  void                        DoClusterLoop()                                   ;
  void                        DoCellLoop()                                      ;
  void                        DoJetCoreLoop()                                      ;
	void												DoMatchingLoop()
		;

	Int_t												SelectTrigger(TList *list,Double_t minT,Double_t maxT,Int_t &number);
	Double_t										RelativePhi(Double_t mphi,Double_t vphi);
	Int_t												GetPhiBin(Double_t phi);

  THistManager                fHistManager;///< Histogram manager

	// flags and selection
	AliEventCuts fEventCuts; ///< Event cuts
	Float_t fJetShapeType; ///<
	Float_t fCentMin; ///< minimum centrality
	Float_t fCentMax; ///< maximum centrality
	Float_t fTTLowRef; ///< minimum reference trigger track pt
	Float_t fTTUpRef; ///< maximum reference trigger track pt
	Float_t fTTLowSig; ///< minimum signal trigger track pt
	Float_t fTTUpSig; ///< maximum signal trigger track pt
	Int_t fNRPBins;	 ///< 
	Float_t fFrac; ///< fraction of events that are used to fill signal recoil jet population
	Float_t fJetHadronDeltaPhi; ///< max angle from pi (set <0 for no selection)
	Float_t fMinFractionSharedPt; ///< min fraction of pt between hybrid / detector jets
	Float_t fMinEmbJetPt; ///< min corrected jet pt to use in embedding
	TString fJetContName; ///< Base level jet container name
	TString fJetContTrueName; ///< True pp (detector) level jet container name
	TString fJetContPartName; ///< Particle(MC) level jet container name
	Bool_t fFillTrackHistograms; ///< switch to fill track histograms
	Bool_t fFillJetHistograms; ///< switch to fill jet histograms
	Bool_t fFillRecoilTHnSparse; ///< switch to fill recoil THnSparse for main analysis
	Bool_t fFillInclusiveTree; ///< switch to fill embedding tree with inclusive jet info
	Bool_t fFillRecoilTree; ///< switch to fill embedding tree with recoil jet info
	Bool_t fMoreTreeVars; ///< add more variables to the output tree
	Int_t fPtHardBin; ///< pt hard bin if running embedding
	Int_t fRejectionFactorInclusiveJets; ///< factor to reject inclusive jets, to reduce size of ttree
	//
	TRandom3 *fRandom; ///<
	Float_t fTreeVarsInclusive[9]; ///<
	Float_t fTreeVarsInclusiveMoreVars[13]; ///<
	Float_t fTreeVarsRecoil[8]; ///<
	Float_t fTreeVarsRecoilMoreVars[12]; ///<
	//histograms to fill
	TH1I *fHistEvtSelection; //!<!
	// recoil jet info contained in THnSparse
	THnSparse *fHJetSpec;  //!<!
	// recoil histograms
	TH1D *fh1TrigRef; //!<!
	TH1D *fh1TrigSig; //!<!
	TH2F *fh2Ntriggers; //!<!
	TH2F *fhRhoCentSig; //!<!
	TH2F *fhRhoCentRef; //!<!
	TH2F *fhDphiPtSigPi; //!<!
	TH2F *fhDphiPtSig; //!<!
	TH2F *fhDphiPtRefPi; //!<!
	TH2F *fhDphiPtRef; //!<!
	// embedding histograms
	// inclusive jets
	TH2F *fhPtDetPart; //!<!
	TH2F *fhPtHybrDet; //!<!
	TH2F *fhPtHybrPart; //!<!
	TH2F *fhPtHybrPartCor; //!<!
	TH2F *fhPhiHybrPartCor; //!<!
	TH1F *fhPtDet; //!<!
	TH1F *fhPtDetMatchedToPart; //!<!
	TH1F *fhPtPartMatched; //!<!
	TH2F *fhPtPartMatchedCent; //!<!
	TH2F *fhPtPartMatchedWrongCent; //!<!
	TH1F *fhResidual; //!<!
	TH2F *fhPtResidual; //!<!
	TH1F *fhPhiResidual; //!<!
	TH2F *fhPhiPhiResidual; //!<!
	//recoil jets
	TH2F *fhPtDetPartRecoil; //!<!
	TH2F *fhPtHybrDetRecoil; //!<!
	TH2F *fhPtHybrPartRecoil; //!<!
	TH2F *fhPtHybrPartCorRecoil; //!<!
	TH2F *fhPhiHybrPartCorRecoil; //!<!
	TH1F *fhPtDetRecoil; //!<!
	TH1F *fhPtDetMatchedToPartRecoil; //!<!
	TH1F *fhResidualRecoil; //!<!
	TH2F *fhPtResidualRecoil; //!<!
	TH1F *fhDphiResidualRecoil; //!<!
	TH2F *fhDphiphiResidualRecoil; //!<!
	TH2F *fhTTPtDetMatchedToPart; //!<!
	TH2F *fhTTPhiDetMatchedToPart; //!<!
	TH2F *fhDPhiHybrPartCorRecoil; //!<!
	TH2F *fhSelectedTrigger; //!<!
	TH2F *fhFractionSharedPtInclusive; //!<!
	TH2F *fhFractionSharedPtRecoil; //!<!
	// embedding trees
	TTree *fTreeEmbInclusive; //!<!
	TTree *fTreeEmbRecoil; //!<!

 private:
  AliAnalysisTaskJetCoreEmcal(const AliAnalysisTaskJetCoreEmcal&)           ; // not implemented
  AliAnalysisTaskJetCoreEmcal &operator=(const AliAnalysisTaskJetCoreEmcal&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetCoreEmcal, 11);
  /// \endcond
};
#endif
