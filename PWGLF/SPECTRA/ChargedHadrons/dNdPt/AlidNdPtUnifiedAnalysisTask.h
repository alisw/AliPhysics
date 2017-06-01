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
class AliAnalysisUtils;

#include "THn.h"
#include "TF1.h"
#include "AliAnalysisTaskSE.h"

class AlidNdPtUnifiedAnalysisTask : public AliAnalysisTaskSE {
  public:
    enum ParticleType {kPrimary=0, kPion=1, kKaon=2, kProtons=3, kSigmaPlus=4, kSigmaMinus=5, kOmegaMinus=6, kXiMinus=7, kElectron=8, kMuon=9,  kRest=10, kLambda=11};
    AlidNdPtUnifiedAnalysisTask(const char *name = "AlidNdPtUnifiedAnalysisTask");
    virtual ~AlidNdPtUnifiedAnalysisTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    // Getters
    
    // Binning
    TArrayD* GetBinsPt() {return fBinsPt;}
    TArrayD* GetBinsEta(){return fBinsEta;}
    TArrayD* GetBinsMultCent() {return fBinsMultCent;}
    TArrayD* GetBinsZv() {return fBinsZv;}


    void SetTriggerMask(UInt_t triggermask)  { fTriggerMask = triggermask; }
    UInt_t GetTriggerMask()  { return fTriggerMask; }

    // Setters
    /// Set the flag for the use of MC.
    void SetUseMC(Bool_t useMC=kTRUE){fIsMC=useMC;}
    /// Set the flag for the use of ESD \c (fIsESD=kTRUE)
    void SetUseESD(){fIsESD=kTRUE;}
    /// Set the flag for the use of AOD \c (fIsESD=kFALSE)
    void SetUseAOD(){fIsESD=kFALSE;}
    /// Set the flag for the particle composition

    void SetTrkEffParametrisation(TF1 *function){fFunTrkEff = function;}

    void SetUseMultiplicity(){fUseMultiplicity = kTRUE;}
    void SetUseCentrality(){fUseMultiplicity = kFALSE;}
    void SetUseCountedMult(){fUseCountedMult = kTRUE;}
    // Binning
    /// Set bins in Pt using a TArrayD
    void SetBinsPt(TArrayD *bins){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(*bins);}
    /// Set bins in Pt using number of bins and array of bin edges
    void SetBinsPt(Int_t nBins, Double_t *binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Eta using a TArrayD
    void SetBinsEta(TArrayD *bins){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(*bins);}
    /// Set bins in Eta using number of bins and array of bin edges
    void SetBinsEta(Int_t nBins, Double_t *binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Multiplicity/Centrality using a TArrayD
    void SetBinsMultCent(TArrayD *bins){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(*bins);}
    /// Set bins in Multiplicity/Centrality number of bins and array of bin edges
    void SetBinsMultCent(Int_t nBins, Double_t *binEdges){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Zv using a TArrayD
    void SetBinsZv(TArrayD *bins){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(*bins);}
    /// Set bins in Zv using number of bins and array of bin edges
    void SetBinsZv(Int_t nBins, Double_t *binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}

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
    void SetMinCrossedRowsTPC(Int_t minCRows){fMinNCrossedRowsTPC = minCRows;}
    void SetRatioCrossedRowsOverFindableClustersTPC(Float_t crossedRoverFindClu){fMinRatioCrossedRowsOverFindableClustersTPC = crossedRoverFindClu;}
    void SetFractionSharedClustersTPC(Float_t sharedclu){fMaxFractionSharedClustersTPC = sharedclu;}
    void SetMaxchi2perTPCclu(Float_t maxchi2TPCclu){fMaxChi2PerTPCCluster = maxchi2TPCclu;}
    void SetClusterReqITS(Bool_t cluReqITS){fRequiresClusterITS = cluReqITS;}
    void SetMaxchi2perITSclu(Float_t maxchi2ITSclu){fMaxChi2PerITSCluster = maxchi2ITSclu;}
    void SetDCAtoVertex2D(Bool_t dcatovertex2d){fDCAToVertex2D = dcatovertex2d;}
    void SetSigmaToVertex(Bool_t sigtovertex){fSigmaToVertex = sigtovertex;}
    void SetDCAtoVertexZ(Float_t dcatovertexz){fMaxDCAzITSTPC = dcatovertexz;}
    void SetDCAtoVertexXYPtDep(const char *dcaxypt){fDCAToVertexXYPtDep = dcaxypt;}
    void SetDCAtoVertexXY(Float_t dcatovertexxy){fDCAToVertexXY = dcatovertexxy;}
    void SetMaxChi2TPCConstrained(Float_t chi2TPCconstrained){fMaxChi2TPCConstrained = chi2TPCconstrained;}
    void SetMinLenghtInActiveZoneTPC(Int_t length){fMinActiveLength = length;}
    void SetGeometricalCut(Bool_t usegeometricalCut, Float_t deadzoneWidth, Float_t ncrnclgeomlength , Float_t ncrnclgeom1pt, Float_t fractionNcr, Float_t fractionNcl  ){fUseGeomCut = usegeometricalCut; fDeadZoneWidth = deadzoneWidth; fCutGeoNcrNclLenght = ncrnclgeomlength; fCutGeoNcrNclGeom1Pt = ncrnclgeom1pt; fCutGeoNcrNclFractionNcl = fractionNcr; fCutGeoNcrNclFractionNcl = fractionNcl;}

    /// TOF pileup -> only for Matching Efficiency calculations
    void SetTOFbunchCrossing(Bool_t isTOFbunch){fUseTOFBunchCrossing=isTOFbunch;}

    /// Event cuts for 2013 and 2015 data
    void Set2013pA(Bool_t is2013) { fIs2013pA = is2013; }
    void Set2015data(Bool_t is2015) {fIs2015data = is2015;}

    Bool_t IsTrackAcceptedKinematics(AliVTrack *track);
    Bool_t IsTrackAcceptedKinematics(TParticle *mcTrack, Bool_t useLowerPtCut = kTRUE);
    Bool_t IsTrackAcceptedQuality(AliVTrack *track);
    Bool_t IsEventAcceptedGeometrics(AliVEvent *event);
    Bool_t IsEventAcceptedQuality(AliVEvent *event);

    Bool_t IsEventAccepted2013pA(AliVEvent *event);
    Bool_t IsEventAccepted2015data(AliVEvent *event);

    Int_t  IdentifyMCParticle(Int_t mcLabel);

    Double_t GetEventMultCent(AliVEvent *event);

    Bool_t IsVertexOK(AliVEvent *event);
    Bool_t IsMCEventINEL0(AliMCEvent* mcEvent, Double_t ptmin, Double_t etarange);
    Bool_t IsTrackAcceptedGeometricalCut(AliVTrack *tr, Double_t bMagZ);
    Bool_t IsChargedPrimary(Int_t stackIndex);
    Bool_t IsChargedPrimaryOrLambda(Int_t stackIndex);

    void InitESDTrackCuts();
    void InitdNdPtEventCuts();
    
    
    void SetCentralityCut(Double_t lowerCut, Double_t upperCut){fUseCentralityCut = kTRUE; fLowerCentralityBound = lowerCut; fUpperCentralityBound = upperCut;}
    Bool_t IsSelectedCentrality();
    
    Bool_t fIncludeSigmas;
    void SetIncludeSigmas(Bool_t includeSigmas){fIncludeSigmas = includeSigmas;}

  private:
    AliVEvent   *fEvent;			//!<! Event object (AliVEvent)
    AliMCEvent  *fMCEvent;		//!<! MC event
    AliStack    *fMCStack;		//!<! MC stack
    TList       *fOutputList;		//!<! Output list

    AliESDtrackCuts   *fESDtrackCuts;
    AlidNdPtEventCuts *fEventCuts;
    UInt_t fTriggerMask;    // trigger mask

    TF1               *fFunTrkEff;

    TH1D* fHistV0Amp;
    THnF* fHistMCMultPt;
    Double_t fLowerCentralityBound;
    Double_t fUpperCentralityBound;
    Bool_t fUseCentralityCut;

    // Output Histograms

    THnF        	*fHistTrack;			///<  Histogram for tracks (pt,eta,Zv,mult/cent)
    THnF	       	*fHistTrackCharge;		///<  Control Histogram for track charge (pt,eta,mult/cent,charge)
    THnF        	*fHistEvent;			///<  Histogram for events (Zv,mult/cent)
    THnF        	*fHistMultEvent;		///<  Histogram for events (mult/cent,acc.mult,corr.mult)

    THnF        	*fHistMCGenPrimTrack;		///<  Histogram for generated MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecTrack;		///<  Histogram for reconstructed MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecPrimTrack;		///<  Histogram for primary MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecSecTrack;		///<  Histogram for secondary MC tracks (pt,eta,mult/cent)

    THnF        	*fHistMCGenPrimTrackParticle;	///<  Particle type histogram for generated MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecTrackParticle;	///<  Particle type histogram for reconstructed MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecPrimTrackParticle;	///<  Particle type histogram for primary MC tracks (pt,eta,mult/cent)
    THnF        	*fHistMCRecSecTrackParticle;	///<  Particle type histogram for secondary MC tracks (pt,eta,mult/cent)

    THnF        	*fHistMCRecEvent;		///<  Histogram for reconstructed MC events (Zv,mult/cent)
    THnF        	*fHistMultMCRecEvent;		///<  Histogram for reconstructed MCevents (mult/cent,acc.mult,corr.mult)
    THnF        	*fHistMultCorrelation;		///<  Histogram for multiplicity correlation (corr.mult,corr.mult)
    THnF        	*fHistMCTrigEvent;		///<  Histogram for triggered MC events (Zv,mult/cent)
    THnF        	*fHistMCGenEvent;		///<  Histogram for generated MC events (Zv,mult/cent)
    THnF        	*fHistMultMCGenEvent;		///<  Histogram for generated MCevents (mult/cent,acc.mult,corr.mult)
    THnF        	*fHistMCGenINEL0Event;		///<  Histogram for generated INEL>0 MC events (Zv,mult/cent)
    THnF        	*fHistMCTrigINEL0Event;   	///<  Histogram for triggered INEL>0 MC events (Zv,mult/cent)
    THnF        	*fHistMCRecINEL0Event;    	///<  Histogram for reconstructed INEL>0 MC events (Zv,mult/cent)
    THnF        	*fHistMCResponseMat;    	///<  Histogram for Detector Response N_ch vs. N_acc

    //   THnF	      *fHistMCGenTrackINEL0;    ///<  Histogram for generated MC tracks for INEL>0 events (pt,eta,mult/cent)

    // Binning

    TArrayD     	*fBinsPt;			///< Array of bins in pt
    TArrayD     	*fBinsEta;		///< Array of bins in eta
    TArrayD     	*fBinsMultCent;		///< Array of bins in multiplicity or centrality

    TArrayD     	*fBinsZv;			///< Array of bins in Zv (Z-position of primary vtx)

    Bool_t      	fUseCountedMult;		///< Flag to select wether to use GetMultiplicity() or N_acc in fHistEvent and fHistTrack

    Bool_t      	fIsMC;			///< Flag for MC usage
    Bool_t      	fIsESD;			///< Flag for ESD usage
    Bool_t      	fUseMultiplicity;         ///< Flag for Multiplicity or Centrality

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
    Int_t       	fMinNCrossedRowsTPC;	///< Minimum number of crossed Rows in the TPC
    Float_t     	fMinRatioCrossedRowsOverFindableClustersTPC;
    Float_t     	fMaxFractionSharedClustersTPC;
    Float_t     	fMaxChi2PerTPCCluster;
    Bool_t        fRequiresClusterITS;
    Float_t     	fMaxChi2PerITSCluster;
    Bool_t      	fDCAToVertex2D;
    Bool_t      	fSigmaToVertex;
    Float_t     	fMaxDCAzITSTPC;
    TString       fDCAToVertexXYPtDep;
    Float_t       fDCAToVertexXY;
    Float_t     	fMaxChi2TPCConstrained;
    Int_t       	fMinActiveLength;
    Bool_t      	fUseGeomCut;
    Float_t       fDeadZoneWidth;
    Float_t       fCutGeoNcrNclLenght;
    Float_t       fCutGeoNcrNclGeom1Pt;
    Float_t       fCutGeoNcrNclFractionNcr;
    Float_t       fCutGeoNcrNclFractionNcl;

    AliAnalysisUtils* fUtils;
    Bool_t fIs2013pA;
    Bool_t fIs2015data;
    Bool_t fUseTOFBunchCrossing;

    AlidNdPtUnifiedAnalysisTask(const AlidNdPtUnifiedAnalysisTask&); // not implemented
    AlidNdPtUnifiedAnalysisTask& operator=(const AlidNdPtUnifiedAnalysisTask&); // not implemented
    /// \cond CLASSIMP
    ClassDef(AlidNdPtUnifiedAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
