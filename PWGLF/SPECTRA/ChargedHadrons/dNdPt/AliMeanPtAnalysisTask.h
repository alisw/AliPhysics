/// \class AliMeanPtAnalysisTask
/// \brief Task to perform Mean Pt analysis
///
/// The Task is designed to run on AOD's and ESD's both for Data and Monte-Carlo.
/// Also the task is supposed to work for all collision systems.
///
///


#ifndef AliMeanPtAnalysisTask_cxx
#define AliMeanPtAnalysisTask_cxx

// class THnF;
class TParticle;

class AliESDEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliAnalysisUtils;

#include "THn.h"
#include "THnSparse.h"
#include "TF1.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskSE.h"

class AliMeanPtAnalysisTask : public AliAnalysisTaskSE {
  public:
    AliMeanPtAnalysisTask(const char *name = "AliMeanPtAnalysisTask");
    virtual ~AliMeanPtAnalysisTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t* option);
    virtual void   Terminate(Option_t*);

    // Getters

    // Binning
    TArrayD* GetBinsPt()    {return fBinsPt;}
    TArrayD* GetBinsEta()   {return fBinsEta;}
    TArrayD* GetBinsMult()  {return fBinsMult;}
    TArrayD* GetBinsCent()  {return fBinsCent;}
    TArrayD* GetBinsZv()    {return fBinsZv;}


    void SetTriggerMask(UInt_t triggermask)  {fTriggerMask = triggermask;}
    UInt_t GetTriggerMask()  {return fTriggerMask;}

    // Setters
    /// Set the flag for the use of MC.
    void SetUseMC(Bool_t useMC = kTRUE){fIsMC = useMC;}
    /// Set the flag for the use of ESD \c (fIsESD=kTRUE)
    void SetUseESD(){fIsESD = kTRUE;}
    /// Set the flag for the use of AOD \c (fIsESD=kFALSE)
    void SetUseAOD(){fIsESD = kFALSE;}
    /// Set the flag for additional Histograms
    void SetIncludeCrosscheckHistos(Bool_t includeHistos = kTRUE){fIncludeCrosscheckHistos = includeHistos;}

    // Binning
    /// Set bins in Pt using a TArrayD
    void SetBinsPt(TArrayD* bins){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(*bins);}
    /// Set bins in Pt using number of bins and array of bin edges
    void SetBinsPt(Int_t nBins, Double_t* binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Eta using a TArrayD
    void SetBinsEta(TArrayD* bins){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(*bins);}
    /// Set bins in Eta using number of bins and array of bin edges
    void SetBinsEta(Int_t nBins, Double_t* binEdges){if(fBinsEta) delete fBinsEta; fBinsEta = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Multiplicity using a TArrayD
    void SetBinsMult(TArrayD* bins){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(*bins);}
    /// Set bins in Multiplicity number of bins and array of bin edges
    void SetBinsMult(Int_t nBins, Double_t* binEdges){if(fBinsMult) delete fBinsMult; fBinsMult = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Centrality using a TArrayD
    void SetBinsCent(TArrayD* bins){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(*bins);}
    /// Set bins in Centrality number of bins and array of bin edges
    void SetBinsCent(Int_t nBins, Double_t* binEdges){if(fBinsCent) delete fBinsCent; fBinsCent = new TArrayD(nBins+1,binEdges);}
    /// Set bins in Zv using a TArrayD
    void SetBinsZv(TArrayD* bins){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(*bins);}
    /// Set bins in Zv using number of bins and array of bin edges
    void SetBinsZv(Int_t nBins, Double_t* binEdges){if(fBinsZv) delete fBinsZv; fBinsZv = new TArrayD(nBins+1,binEdges);}

    void SetBinsPtReso(Int_t nBins, Double_t* binEdges){if(fBinsPtReso) delete fBinsPtReso; fBinsPtReso = new TArrayD(nBins+1,binEdges);}
    void SetBins1Pt(Int_t nBins, Double_t* binEdges){if(fBins1Pt) delete fBins1Pt; fBins1Pt = new TArrayD(nBins+1,binEdges);}
    void SetBinsSigma1Pt(Int_t nBins, Double_t* binEdges){if(fBinsSigma1Pt) delete fBinsSigma1Pt; fBinsSigma1Pt = new TArrayD(nBins+1,binEdges);}

    // Acceptance cuts
    /// Set the minimum Eta cut
    void SetMinEta(Double_t minEta){fMinEta = minEta;}
    /// Set the maximum Eta cut
    void SetMaxEta(Double_t maxEta){fMaxEta = maxEta;}
    /// Set the minimum Pt cut
    void SetMinPt(Double_t minPt){fMinPt = minPt;}
    /// Set the maximum Pt cut
    void SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}

    void SetSigmaMeanXYZv(Float_t sigmaXv, Float_t sigmaYv, Float_t sigmaZv){fSigmaMeanXYZv[0] = sigmaXv; fSigmaMeanXYZv[1] = sigmaYv; fSigmaMeanXYZv[2] = sigmaZv;}
    void SetMeanXYZv(Float_t meanXv, Float_t meanYv, Float_t meanZv){fMeanXYZv[0] = meanXv; fMeanXYZv[1] = meanYv; fMeanXYZv[2] = meanZv;}
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
    void SetDCAtoVertexXYPtDep(const char* dcaxypt){fDCAToVertexXYPtDep = dcaxypt;}
    void SetDCAtoVertexXY(Float_t dcatovertexxy){fDCAToVertexXY = dcatovertexxy;}
    void SetMaxChi2TPCConstrained(Float_t chi2TPCconstrained){fMaxChi2TPCConstrained = chi2TPCconstrained;}
    void SetMinLenghtInActiveZoneTPC(Int_t length){fMinActiveLength = length;}
    void SetGeometricalCut(Bool_t usegeometricalCut, Float_t deadzoneWidth, Float_t ncrnclgeomlength , Float_t ncrnclgeom1pt, Float_t fractionNcr, Float_t fractionNcl  ){fUseGeomCut = usegeometricalCut; fDeadZoneWidth = deadzoneWidth; fCutGeoNcrNclLenght = ncrnclgeomlength; fCutGeoNcrNclGeom1Pt = ncrnclgeom1pt; fCutGeoNcrNclFractionNcl = fractionNcr; fCutGeoNcrNclFractionNcl = fractionNcl;}

    /// Event cuts for 2013 and 2015 data
    void Set2013pA(Bool_t is2013) { fIs2013pA = is2013; }
    void Set2015data(Bool_t is2015) {fIs2015data = is2015;}

    Bool_t IsTrackInKinematicRange(AliVTrack* track);
    Bool_t IsParticleInKinematicRange(AliMCParticle* mcParticle);

    Bool_t IsTrackAcceptedQuality(AliVTrack* track);
    Bool_t IsEventVertexOK(AliVEvent* event);

    Bool_t IsEventAccepted2013pA(AliVEvent* event);
    Bool_t IsEventAccepted2015data(AliVEvent* event);

    Double_t GetCentrality(AliVEvent* event);

    Bool_t IsVertexOK(AliVEvent* event);
    Bool_t IsChargedPrimary(Int_t mcLabel);

    void InitESDTrackCuts();
    void InitdNdPtEventCuts();
    void SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins);


  private:
    Double_t            fPRECISION;
    TList*              fOutputList;		//!<! Output list
    AliVEvent*          fEvent;			    //!<! Event object
    AliMCEvent*         fMCEvent;       //!<! MC event

    AliESDtrackCuts*    fESDtrackCuts;

    AliAnalysisUtils*   fUtils;
    Bool_t              fIsESD;			    ///< Flag for ESD usage
    Bool_t      	      fIsMC;			    ///< Flag for MC usage
    Bool_t	 			      fIs2013pA;
    Bool_t 				      fIs2015data;

    Bool_t      	fTPCRefit;		     ///< TPC refit
    Bool_t      	fITSRefit;		     ///< TPC refit
    Bool_t      	fAcceptKinks; 		 ///< Accept Kink Daughters
    Bool_t        fRequiresClusterITS;
    Bool_t      	fDCAToVertex2D;
    Bool_t      	fSigmaToVertex;
    Bool_t      	fUseGeomCut;

    Bool_t      	fIncludeCrosscheckHistos;


    // Acceptance cuts for tracks
    UInt_t        fTriggerMask;   // trigger mask
    Double_t     	fMinEta;			///< Minimum eta cut
    Double_t     	fMaxEta;			///< Maximum eta cut
    Double_t     	fMinPt;			  ///< Minimum pT cut
    Double_t     	fMaxPt;			  ///< Maximum pT cut

    Float_t     	fSigmaMeanXYZv[3];	///<[3]
    Float_t     	fMeanXYZv[3];		///<[3]
    Float_t     	fZvtx;

    Int_t       	fMinNCrossedRowsTPC;	///< Minimum number of crossed Rows in the TPC
    Float_t     	fMinRatioCrossedRowsOverFindableClustersTPC;
    Float_t     	fMaxFractionSharedClustersTPC;
    Float_t     	fMaxChi2PerTPCCluster;
    Float_t     	fMaxChi2PerITSCluster;
    Float_t     	fMaxDCAzITSTPC;
    TString      	fDCAToVertexXYPtDep;
    Float_t       fDCAToVertexXY;
    Float_t     	fMaxChi2TPCConstrained;
    Int_t       	fMinActiveLength;
    Float_t       fDeadZoneWidth;
    Float_t       fCutGeoNcrNclLenght;
    Float_t       fCutGeoNcrNclGeom1Pt;
    Float_t       fCutGeoNcrNclFractionNcr;
    Float_t       fCutGeoNcrNclFractionNcl;

    TArrayD*      fBinsMult;		///< Array of bins in multiplicity
    TArrayD*      fBinsCent;		///< Array of bins in centrality
    TArrayD*      fBinsPt;			///< Array of bins in pt
    TArrayD*      fBinsEta;		///< Array of bins in eta
    TArrayD*      fBinsZv;			///< Array of bins in Zv (Z-position of primary vtx)
    TArrayD*      fBinsPtReso;			   ///< Array of bins for relative pt resoulution
    TArrayD*      fBins1Pt;			       ///< Array of bins for 1/pt
    TArrayD*      fBinsSigma1Pt;			///< Array of bins for 1/pt resoulution


    // Output Histograms
    TH1F*         fEventCount;		        ///< Histogram for triggered events and events with vertex
    TH1F*         fHistMCTrackParticle;		///<

    THnF*         fHistEvent;			        ///<  Histogram for events

    THnSparseF*   fHistMCResponseMat;    	///<  Histogram for Detector Response N_ch vs. N_acc
    THnSparseF*   fHistMCResponseMatTracks;    	///<  Histogram for Detector Response N_ch vs. N_acc

    THnF*         fHistTrack;			///<  Histogram for tracks (pt,eta,Zv,mult/cent)
    THnF*         fHistRelPtResoFromCov;			///<  Histogram for relative pT resolution of tracks from covariance matrix

    THnF*         fHistMCRecTrack;		///<  Histogram for reconstructed MC tracks (pt,eta,mult/cent)
    THnF*         fHistMCGenPrimTrack;		///<  Histogram for generated MC tracks (pt,eta,mult/cent)
    THnF*         fHistMCRecPrimTrack;		///<  Histogram for primary MC tracks (pt,eta,mult/cent)
    THnF*         fHistMCRecSecTrack;		///<  Histogram for secondary MC tracks (pt,eta,mult/cent)

    THnF* 		    fHistMCMultPtGenerated;
    THnSparseF*   fHistMCTrackMultGen;		///<  Histogram for true tracks vs multiplicity (pt,Nacc,Nch)

    THnF*         fHistMCPtRes;                    ///<  Histogram for pT_gen vs pT_rec for resolution chrosschecks
    THnF*         fHistMCRelPtReso;                    ///<  Histogram for relative pt resolution vs pT_gen vs pT_rec vs cent
    THnF*         fHistMCEtaRes;                    ///<  Histogram for eta_gen vs eta_rec for resolution chrosschecks
    THnSparseF*   fHistMCMultRes;                    ///<  Histogram for Nacc vs Nrec for resolution chrosschecks

    THnF*         fHistMCParticle;			///<  Histogram for particles (pt,eta, mult, cent)


    AliMeanPtAnalysisTask(const AliMeanPtAnalysisTask&); // not implemented
    AliMeanPtAnalysisTask& operator=(const AliMeanPtAnalysisTask&); // not implemented
    /// \cond CLASSIMP
    ClassDef(AliMeanPtAnalysisTask, 1); // example of analysis
    /// \endcond
};

#endif
