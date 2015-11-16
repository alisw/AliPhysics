/// \class AlidNdPtUnifiedAnalysisTask
/// \brief Unified Task to perform Spectra and Mean Pt analysis
///
/// The Task is designed to run on AOD's and ESD's both for Data and Monte-Carlo.
/// Also the task is supposed to work for all collision systems.
///
///
/// \author Michael Knichel <michael.linus.knichel@cern.ch>, University Heidelberg
/// \author Philipp Luettig <luettig@ikf.uni-frankfurt.de>, University Frankfurt
/// \author Patrick Huhn <phuhn@ikf.uni-frankfurt.de>, University Frankfurt
/// \author Federica Sozzi <f.sozzi@gsi.de>, GSI Darmstadt
/// \author Edgar Perez Lezama <e.perezlezama@gsi.de>, GSI Darmstadt
/// \author Tatiana Drozhzhova <tatiana.drozhzhova@cern.ch>, GSI Darmstadt
/// \author Julius Gronefeld <j.gronefeld@cern.ch>, GSI Darmstadt
/// \date May 29, 2015


#ifndef AlidNdPtUnifiedAnalysisTask_cxx
#define AlidNdPtUnifiedAnalysisTask_cxx

// class THnF;
class TParticle;

class AliESDEvent;
class AliVEvent;
class AliStack;
class AliESDtrackCuts;
class AlidNdPtEventCuts;

#include "THn.h"
#include "AliAnalysisTaskSE.h"

class AlidNdPtUnifiedAnalysisTask : public AliAnalysisTaskSE {
 public:
  enum ParticleType {kPrimary=0, kPion=211, kKaon=321, kProtons=2212, kRest=-1};
  AlidNdPtUnifiedAnalysisTask(const char *name = "AlidNdPtUnifiedAnalysisTask");
  virtual ~AlidNdPtUnifiedAnalysisTask() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // Getters
    
  // Binning
  TArrayD* GetBinsPt() {return fBinsPt;}
  TArrayD* GetBinsEta(){return fBinsEta;}
  TArrayD* GetBinsMultCent() {return fBinsMultCent;}
  TArrayD* GetBinsZv() {return fBinsZv;}
  
  // Setters
  /// Set the flag for the use of MC.
  void SetUseMC(Bool_t useMC=kTRUE){fIsMC=useMC;}
  /// Set the flag for the use of ESD \c (fIsESD=kTRUE)
  void SetUseESD(){fIsESD=kTRUE;}
  /// Set the flag for the use of AOD \c (fIsESD=kFALSE)
  void SetUseAOD(){fIsESD=kFALSE;}
  /// Set the flag for the particle composition
  void SetMCParticleType(AlidNdPtUnifiedAnalysisTask::ParticleType pLabel) {fParticleLabel = pLabel;}
  
  void SetUseMultiplicity(){fUseMultiplicity = kTRUE;}
  void SetUseCentrality(){fUseMultiplicity = kFALSE;}
  // Binning
  /// Set bins in Pt using a TArrayD
  void SetBinsPt(TArrayD *bins){fBinsPt = new TArrayD(*bins);}
  /// Set bins in Pt using number of bins and array of bin edges
  void SetBinsPt(Int_t nBins, Double_t *binEdges){fBinsPt = new TArrayD(nBins+1,binEdges);}
  /// Set bins in Eta using a TArrayD
  void SetBinsEta(TArrayD *bins){fBinsEta = new TArrayD(*bins);}
  /// Set bins in Eta using number of bins and array of bin edges
  void SetBinsEta(Int_t nBins, Double_t *binEdges){fBinsEta = new TArrayD(nBins+1,binEdges);}
  /// Set bins in Multiplicity/Centrality using a TArrayD
  void SetBinsMultCent(TArrayD *bins){fBinsMultCent = new TArrayD(*bins);}
  /// Set bins in Multiplicity/Centrality number of bins and array of bin edges
  void SetBinsMultCent(Int_t nBins, Double_t *binEdges){fBinsMultCent = new TArrayD(nBins+1,binEdges);}
  /// Set bins in Zv using a TArrayD
  void SetBinsZv(TArrayD *bins){fBinsZv = new TArrayD(*bins);}
  /// Set bins in Zv using number of bins and array of bin edges
  void SetBinsZv(Int_t nBins, Double_t *binEdges){fBinsZv = new TArrayD(nBins+1,binEdges);}
  
  // Acceptance cuts
  /// Set the minimum Eta cut
  void SetMinEta(Float_t minEta){fMinEta = minEta;}
  /// Set the maximum Eta cut
  void SetMaxEta(Float_t maxEta){fMaxEta = maxEta;}
  /// Set the minimum Pt cut
  void SetMinPt(Float_t minPt){fMinPt = minPt;}
  /// Set the maximum Pt cut
  void SetMaxPt(Float_t maxPt){fMaxPt = maxPt;}
  
  void SetSigmaMeanXYZv(Float_t sigmaXv, Float_t sigmaYv, Float_t sigmaZv){fSigmaMeanXYZv[0] = sigmaXv; fSigmaMeanXYZv[1] = sigmaYv; fSigmaMeanXYZv[2] = sigmaZv;}
  void SetMeanXYZv(Float_t meanXv, Float_t meanYv, Float_t meanZv){fMeanXYZv[0] = meanXv; fMeanXYZv[1] = meanYv; fMeanXYZv[2] = meanZv;}
  void SetEventTriggerRequired(Bool_t eventtrigger){fEventTriggerRequired = eventtrigger;}
  void SetZvtx(Float_t zvtx){fZvtx = zvtx;}
  void SetTPCRefit(Bool_t tpcrefit){fTPCRefit = tpcrefit;}
  void SetITSRefit(Bool_t itsrefit){fITSRefit = itsrefit;}
  void SetKinkDaughters(Bool_t kinkD){fAcceptKinks = kinkD;}
  void SetMinCrossedRowsTPC(Int_t minCRows){fminNCrossedRowsTPC = minCRows;}
  void SetRatioCrossedRowsOverFindableClustersTPC(Float_t crossedRoverFindClu){fminRatioCrossedRowsOverFindableClustersTPC = crossedRoverFindClu;}
  void SetFractionSharedClustersTPC(Float_t sharedclu){fmaxFractionSharedTPCCluster = sharedclu;}
  void SetMaxchi2perTPCclu(Float_t maxchi2TPCclu){fmaxchi2perTPCcl = maxchi2TPCclu;}
  void SetMaxchi2perITSclu(Float_t maxchi2ITSclu){fmaxchi2perITScl = maxchi2ITSclu;}
  void SetDCAtoVertex2D(Bool_t dcatovertex2d){fdcatoVertex2D = dcatovertex2d;}
  void SetSigmaToVertex(Bool_t sigtovertex){fsigmatovertex = sigtovertex;}
  void SetDCAtoVertexZ(Float_t dcatovertexz){fmaxdcazITSTPC = dcatovertexz;}
  void SetMaxChi2TPCConstrained(Float_t chi2TPCconstrained){fmaxchi2TPCconstrained = chi2TPCconstrained;}
  void SetMinLenghtInActiveZoneTPC(Int_t length){fminActiveLength = length;}
  void SetGeometricalCut(Bool_t usegeometricalCut){fUseGeomCut = usegeometricalCut;}
  
  Bool_t IsTrackAcceptedKinematics(AliVTrack *track);
  Bool_t IsTrackAcceptedKinematics(TParticle *mcTrack);
  Bool_t IsTrackAcceptedQuality(AliVTrack *track);
  Bool_t IsEventAcceptedGeometrics(AliVEvent *event);
  Bool_t IsEventAcceptedQuality(AliVEvent *event);
  
  Bool_t IsSelectedParticle(Int_t mcLabel);
  Bool_t IsSecondary(Int_t mcLabel);
  
  Double_t GetEventMultCent(AliVEvent *event);
  
  Bool_t IsVertexOK(AliVEvent *event);
  Bool_t SelectMCEventINEL0(AliMCEvent* mcEvent, Double_t ptmin, Double_t etarange);
  Bool_t IsTrackAcceptedGeometricalCut(AliVTrack *tr, Double_t bMagZ);
  
  void InitESDTrackCuts();
  void InitdNdPtEventCuts();

  
  

 private:
  AliVEvent   *fEvent;			//!<! Event object (AliVEvent)
  AliMCEvent  *fMCEvent;		//!<! MC event
  AliStack    *fMCStack;		//!<! MC stack
  TList       *fOutputList;		//!<! Output list
  
  AliESDtrackCuts   *fESDtrackCuts;
  Bool_t            fESDtrackCutsInit;
  AlidNdPtEventCuts *fEventCuts;
  Bool_t            fEventCutsInit;
  
  
  // Output Histograms
  
  THnF        	*fHistTrack;		///<  Histogram for tracks (pt,eta,Zv,mult/cent)
  THnF        	*fHistEvent;		///<  Histogram for events (Zv,mult/cent)
  
  THnF        	*fHistMCGenPrimTrack;	///<  Histogram for generated MC tracks (pt,eta,mult/cent)
  THnF        	*fHistMCRecTrack;		///<  Histogram for reconstructed MC tracks (pt,eta,mult/cent)
  THnF        	*fHistMCRecPrimTrack;	///<  Histogram for primary MC tracks (pt,eta,mult/cent) 
  THnF        	*fHistMCRecSecTrack;	///<  Histogram for secondary MC tracks (pt,eta,mult/cent)
  THnF        	*fHistMCRecEvent;		///<  Histogram for reconstructed MC events (Zv,mult/cent)
  THnF        	*fHistMCTrigEvent;	///<  Histogram for triggered MC events (Zv,mult/cent)
  THnF        	*fHistMCGenEvent;		///<  Histogram for generated MC events (Zv,mult/cent)
  THnF        	*fHistMCGenINEL0Event;	///<  Histogram for generated INEL>0 MC events (Zv,mult/cent)
  THnF        	*fHistMCTrigINEL0Event;   ///<  Histogram for triggered INEL>0 MC events (Zv,mult/cent)
  THnF        	*fHistMCRecINEL0Event;    ///<  Histogram for reconstructed INEL>0 MC events (Zv,mult/cent)
//   THnF	      *fHistMCGenTrackINEL0;    ///<  Histogram for generated MC tracks for INEL>0 events (pt,eta,mult/cent)
  
  // Binning
  
  TArrayD     	*fBinsPt;			///< Array of bins in pt
  TArrayD     	*fBinsEta;		///< Array of bins in eta
  TArrayD     	*fBinsMultCent;		///< Array of bins in multiplicity or centrality
  
  TArrayD     	*fBinsZv;			///< Array of bins in Zv (Z-position of primary vtx)
  
  Bool_t      	fIsMC;			///< Flag for MC usage
  Bool_t      	fIsESD;			///< Flag for ESD usage
  Bool_t      	fUseMultiplicity;         ///< Flag for Multiplicity or Centrality
  Bool_t      	fIsEventINEL0;            ///< Flag for INEL>0 event class
  
  // Acceptance cuts for tracks
  Float_t     	fMinEta;			///< Minimum eta cut
  Float_t     	fMaxEta;			///< Maximum eta cut
  Float_t     	fMinPt;			///< Minimum pT cut 
  Float_t     	fMaxPt;			///< Maximum pT cut 
  
  Float_t     	fSigmaMeanXYZv[3];	///<[3]
  Float_t     	fMeanXYZv[3];		///<[3]
  Bool_t      	fEventTriggerRequired;
  Float_t     	fZvtx;
  Bool_t      	fTPCRefit;		///< TPC refit
  Bool_t      	fITSRefit;		///< TPC refit
  Bool_t      	fAcceptKinks; 		///< Accept Kink Daughters
  Int_t       	fminNCrossedRowsTPC;	///< Minimum number of crossed Rows in the TPC
  Float_t     	fminRatioCrossedRowsOverFindableClustersTPC; 
  Float_t     	fmaxFractionSharedTPCCluster;
  Float_t     	fmaxchi2perTPCcl;
  Float_t     	fmaxchi2perITScl;
  Bool_t      	fdcatoVertex2D;
  Bool_t      	fsigmatovertex;
  Float_t     	fmaxdcazITSTPC;
  Float_t     	fmaxchi2TPCconstrained;
  Int_t       	fminActiveLength;
  Bool_t      	fUseGeomCut;
  ParticleType	fParticleLabel;		///< Particle type dependent MC
  
  AlidNdPtUnifiedAnalysisTask(const AlidNdPtUnifiedAnalysisTask&); // not implemented
  AlidNdPtUnifiedAnalysisTask& operator=(const AlidNdPtUnifiedAnalysisTask&); // not implemented
/// \cond CLASSIMP
  ClassDef(AlidNdPtUnifiedAnalysisTask, 1); // example of analysis
/// \endcond  
};

#endif
