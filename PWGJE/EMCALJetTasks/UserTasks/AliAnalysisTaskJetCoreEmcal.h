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
	virtual void		 SetCentrality(Float_t cmin, Float_t cmax){fCentMin=cmin; fCentMax=cmax;}
	virtual void     SetTTLowRef(Float_t ttlow){fTTLowRef=ttlow;}
	virtual void     SetTTUpRef(Float_t ttup){fTTUpRef=ttup;}
	virtual void     SetTTLowSig(Float_t ttlows){fTTLowSig=ttlows;}
	virtual void     SetTTUpSig(Float_t ttups){fTTUpSig=ttups;}
	virtual void		 SetNRPBins(Float_t nrpb){fNRPBins=nrpb;}
	virtual void		 SetSignalFraction(Float_t sfrac){fFrac=sfrac;}
	virtual void     SetJetEtaMin(Float_t eta){fJetEtaMin=eta;}
	virtual void     SetJetEtaMax(Float_t eta){fJetEtaMax=eta;}
	virtual void     SetJetHadronDeltaPhi(Float_t delta){fJetHadronDeltaPhi=delta;}
	virtual void		 SetJetContName(TString cont){fJetContName=cont;}
	virtual void		 SetRunAnaAzimuthalCorrelation(Bool_t b){fRunAnaAzimuthalCorrelation=b;}

//  static AliAnalysisTaskJetCoreEmcal* AddTaskJetCoreEmcal(
//      const char *ntracks            = "usedefault",
//      const char *nclusters          = "usedefault",
//      const char* ncells             = "usedefault",
//      const char *suffix             = "");

  enum JetShapeType {
    kMCTrue = 0,   // generated jets only
    kTrueDet =1,  // detector and generated jets  
    kData   = 2,  // raw data 
//    kDetEmb = 3,  //detector embedded jets
//    kDetEmbPart = 4,
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

	Int_t												SelectTrigger(TList *list,Double_t minT,Double_t maxT,Int_t &number);
	THnSparse*									NewTHnSparseF(const char* name, UInt_t entries);
	void												GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
	Double_t										RelativePhi(Double_t mphi,Double_t vphi);
	Int_t												GetPhiBin(Double_t phi);

  THistManager                fHistManager                                      ;///< Histogram manager

	// flags and selection
	AliEventCuts fEventCuts; ///< Event cuts
	Float_t fCentMin; ///<
	Float_t fCentMax; ///<
	Float_t fTTLowRef; ///<
	Float_t fTTUpRef; ///<
	Float_t fTTLowSig; ///<
	Float_t fTTUpSig; ///<
	Int_t fNRPBins;	 ///<
	Float_t fFrac; ///<
	Float_t fJetEtaMin; ///<
	Float_t fJetEtaMax; ///<
	Float_t fJetHadronDeltaPhi; ///<
	TString fJetContName; ///<
	Bool_t fRunAnaAzimuthalCorrelation; ///<
	//
	TRandom3 *fRandom; ///<
	//histograms to fill
	TH1I *fHistEvtSelection; //!<!
	THnSparse *fHJetSpec;  //!<!
	TH1D *fh1TrigRef; //!<!
	TH1D *fh1TrigSig; //!<!
	TH2F *fh2Ntriggers; //!<!
	TH2F *fh2RPJetsC10; //!<!
	TH2F *fh2RPJetsC20; //!<!
	TH2F *fh2RPTC10; //!<!
	TH2F *fh2RPTC20; //!<!

	THnSparse *fHJetPhiCorr; //!<!
	TH2F *fhDphiPtSig; //!<!
	TH2F *fhDphiPtRef; //!<!

 private:
  AliAnalysisTaskJetCoreEmcal(const AliAnalysisTaskJetCoreEmcal&)           ; // not implemented
  AliAnalysisTaskJetCoreEmcal &operator=(const AliAnalysisTaskJetCoreEmcal&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetCoreEmcal, 3);
  /// \endcond
};
#endif
