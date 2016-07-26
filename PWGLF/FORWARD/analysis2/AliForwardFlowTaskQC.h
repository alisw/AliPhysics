//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDFLOWTASKQC_H
#define ALIFORWARDFLOWTASKQC_H
/**
 * @file AliForwardFlowTaskQC.h
 * @author Alexander Hansen
 * 
 * @brief
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include <TH2D.h>
class AliAODForwardMult;
class TH1I;
class TH1D;
class TH2I;
class TH2F;
class TH2D;
class TH3D;
class TAxis;
class AliAnalysisFilter;
class AliESDEvent;
/**
 * @defgroup pwglf_forward_tasks_flow Flow tasks 
 *
 * Code to do with flow 
 *
 * @ingroup pwglf_forward_tasks
 */
/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - forward_flow.root
 *
 * @ingroup pwglf_forward_flow
 *
 */
class AliForwardFlowTaskQC : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   */
  AliForwardFlowTaskQC();
  /** 
   * Constructor
   * 
   * @param name Name of task 
   */
  AliForwardFlowTaskQC(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardFlowTaskQC() {}
  /** 
   * @{ 
   * @name Task interface methods 
   */
  /** 
   * Create output objects 
   */
  virtual void UserCreateOutputObjects();
  /**
   * Initialize the task
   */
  virtual void Init() {} 
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
  virtual void UserExec(Option_t *option);
  /** 
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t *option);
  /* @} */
  /**
   * Returns the outputlist
   * 
   * @return TList* 
   */
  TList* GetOutputList() { return fOutputList; }
  /**
   * Set max flow moment to calculate.
   * 
   * @param  n Do v_2 to v_n
   */
  void SetMaxFlowMoment(Short_t n) { fMaxMoment = n; } 
  /**
   * Set vertex binning and range
   *
   * @param axis Use this vtx axis
   */
  void SetVertexAxis(TAxis* axis) { fVtxAxis = axis; }
  /**
   * Set centrality/multiplicity binning and range
   *
   * @param axis Use this binning
   */
  void SetCentralityAxis(TAxis* axis) { fCentAxis = axis; }
  /** 
   * Set the centrality axis to use based on a string.  The bin edges
   * are separated by colons.
   * 
   * @param bins String of bin edges
   */
  void SetCentralityAxis(const char* bins);
  /**
   * Set detector sigma cuts
   *
   * @param fmdCut FMD sigma cut
   * @param spdCut SPD sigma cut
   */
  void SetDetectorCuts(Double_t fmdCut, Double_t spdCut) { fFMDCut = fmdCut; fSPDCut = spdCut; }
  /**
   * Set flow flags, @f$\eta@f$-gap, sym. around @f$\eta=0@f$ or
   * sat. vtx. interactions, also sets which forward detector to use
   *
   * @param flags EFlowFlags 
   */
  void SetFlowFlags(UShort_t flags);
  /**
    * Get QC type
    *
    * @param flags EFlowFlags
    * @param prependUS prepend an underscore
    *
    * @return type
    */
  static const Char_t* GetQCType(UShort_t flags, Bool_t prependUS = kTRUE);
  /**
   *  Set @f$|\eta|@f$ value to make cut for @f$\eta@f$ gap at
   *
   * @param eg gap value
   */
  void SetEtaGapValue(Double_t eg) { fEtaGap = eg; }
  void SetTrackCuts(AliAnalysisFilter* trCuts) { fTrackCuts = trCuts; }
  /**
   * Enum for flow flags
   */
  enum EFlowFlags {
    kStdQC   = 0x0001, // Standard QC{2} and QC{4} calculations
    kEtaGap  = 0x0002, // QC{2} w/ an eta-gap
    k3Cor    = 0x0004, // 3 correlator method for QC{2} w/ an eta-gap
    kSymEta  = 0x0008, // Symmetrize ref flow in std. QC{2} and QC{4} around eta = 0
    kSatVtx  = 0x0010, // Do satellite vertex input (currently not implemented)
    kNUAcorr = 0x0020, // Apply full NUA correction
    kFMD     = 0x0040, // Use FMD for forward flow
    kVZERO   = 0x0080, // Use VZERO for forward flow
    kSPD     = 0x0100, // SPD object flag
    kMC      = 0x0200, // MC object flag
    kTracks  = 0x1000, // Use tracks for reference flow
    kTPC     = 0x3000, // Use TPC tracks
    kHybrid  = 0x5000  // Use hybrid tracks
  };
  /**
   * struct to handle cumulant calculations and control histograms.
   * Used internally and never streamed.
   */
  struct CumuHistos : public TObject
  {
  public:
    /*
     * Constructor
     */
    CumuHistos() : fMaxMoment(), fRefHists(), fDiffHists(), fNUA() {}
    /**
     * Constructor
     * 
     * @param n max flow moment contained
     * @param nua Make room for NUA corrected histograms
     */
    CumuHistos(Int_t n, UInt_t nua) : fMaxMoment(n), fRefHists(), fDiffHists(), fNUA(nua) {}
    /**
     * Copy constructor 
     * 
     * @param o Object to copy from
     * 
     * @return CumuHistos
     */
    CumuHistos(const CumuHistos& o)
      : TObject(o),
        fMaxMoment(o.fMaxMoment), // Max moment to compute
        fRefHists(o.fRefHists),   // List with ref hists
        fDiffHists(o.fDiffHists), // List with diff hists
        fNUA(o.fNUA)
    {}
    /**
     * Assignment operator 
     * 
     * @param o Object to assing from
     * 
     * @return reference to this 
     */
    CumuHistos& operator=(const CumuHistos& o) 
    { 
      if (&o == this) return *this;
      TObject::operator=(o);
      fMaxMoment = o.fMaxMoment;
      fRefHists  = o.fRefHists;
      fDiffHists = o.fDiffHists;
      fNUA       = o.fNUA;
      return *this;
    }
    /**
     * Destructor 
     */
    ~CumuHistos(){}
    /**
     * To access histograms
     * main function of this class
     *
     * @param t (r)eference or (d)iff
     * @param n    flow moment 
     * @param nua  nua type
     * 
     * @return requested histogram
     */
    TH1* Get(Char_t t, Int_t n, UInt_t nua = 0) const;
    /**
     * Connect internal lists to output
     *
     * @param name Name of VertexBin
     * @param l    Output list
     */
    void ConnectList(TString name, TList* l);
    /**
     * Make histograms to one of the lists
     *
     * @param h Hist to add
     */
     void Add(TH1* h) const;
    /** 
     * Check to see of lists are connected, 
     * needed for grid/proof
     *
     * @return is connected?
     */
    Bool_t IsConnected() { return (fRefHists && fDiffHists); }
    /**
     * enum for NUA histograms
     */
    enum {
      kNoNUA = 0, // No NUA applied
      kNUAOld,    // NUA correction from same moment applied
      kNUA        // Full cross-moment NUA correction applied
    };
  protected:
    /**
     * Get position of histogram in list
     *
     * @param n moment to get position of
     * @param nua nua type
     * 
     * @return position
     */
    Int_t GetPos(Int_t n, UInt_t nua) const;

    Int_t  fMaxMoment;  // Max flow moment contained
    TList* fRefHists;   // List of reference hists
    TList* fDiffHists;  // List of diff hists
    UInt_t fNUA;        // NUA tracker

    // ClassDef(CumuHistos, 1);
  }; // End of struct

protected:
  /**
   * Enum for filling flow histos
   */
  enum {
    kFillRef  = 0x1, // Fill only ref flow
    kFillDiff = 0x2, // Fill only diff flow
    kFillBoth = 0x3, // Fill both
    kReset    = 0x4, // Reset hists (used with one of the above)
  };
  /**
   * Enum for event diagnostics
   */
  enum {
    kNoEvent = 1, // No event found
    kNoForward,   // No forward object found
    kNoCentral,   // No central object found
    kNoTrigger,   // No (wrong) trigger
    kNoCent,      // No centrality
    kInvCent,     // Centrality outside range
    kNoVtx,       // No vertex
    kInvVtx,      // Vertex outside range
    kOK           // Event OK!
  };
  // ----------------- Being nested class ---------------------
  /**
   * Nested class to handle cumulant calculations in vertex bins.
   * Used internally and never streamed.
   */
  class VertexBin : public TNamed
  {
  public:
    /*
     * Constructor
     */
    VertexBin();
    /**
     * Constructor
     * 
     * @param vLow Min vertex z-coordinate
     * @param vHigh Max vertex z-coordinate
     * @param moment Flow moment
     * @param type Data type (FMD/VZERO/SPD/FMDTR/SPDTR/MC)
     * @param flags Flags
     * @param cut Cut value 
     * @param etaGap @f$\eta@f$ gap 
     */
    VertexBin(Int_t vLow, Int_t vHigh, 
              UShort_t moment, TString type,
              UShort_t flags, 
              Double_t cut = -1, Double_t etaGap = -1.);
    /**
     * Copy constructor 
     * 
     * @param o Object to copy from
     * 
     * @return VertexBin
     */
    VertexBin(const VertexBin& o){;}
    /**
     * Assignment operator 
     * 
     * @param v Object to assing from
     * 
     * @return reference to this 
     */
    VertexBin& operator=(const VertexBin& v);
    /**
     * Destructor 
     */
    ~VertexBin(){}
    /**
     * Add vertex bin output to list
     * 
     * @param list Histograms are added to this list
     * @param centAxis Axis to handle centrality binning
     * 
     * @return void 
     */
    void AddOutput(TList* list, TAxis* centAxis);
    /**
     * Fill reference and differential flow histograms for analysis
     * using histograms as input
     *
     * @param dNdetadphi 2D data histogram
     * @param cent Centrality
     * @param mode fill ref/diff or both
     *
     * @return false if bad event (det. hotspot)
     */
    Bool_t FillHists(TH2D& dNdetadphi, Double_t cent, UShort_t mode);
    /** 
     * Fill reference and differential flow histograms for analysis 
     * using tracks as input
     *
     * @param trList Array with tracks
     * @param esd  ESD event
     * @param trFilter analysis filter 
     * @param mode fill ref/diff or both
     *
     * @return false if bad event (det. hotspot)
     */
    Bool_t FillTracks(TObjArray* trList, AliESDEvent* esd,
		      AliAnalysisFilter* trFilter, UShort_t mode);
    /**
     * Do cumulants calculations for current event with 
     * centrality cent
     * 
     * @param cent Event centrality
     */
    void CumulantsAccumulate(Double_t cent);
    /**
     * Do 3 correlator cumulants calculations for current event with 
     * centrality cent
     * 
     * @param cent Event centrality
     */
    void CumulantsAccumulate3Cor(Double_t cent);
    /**
     * Get limits to do reference flow calculations for 3 correlator method
     *
     * @param bin Differential bin
     * @param aLow Lowest bin to be used for v_A
     * @param aHigh Highest bin to be used for v_A
     * @param bLow Lowest bin to be used for v_B
     * @param bHigh Highest bin to be used for v_B
     */
    void GetLimits(Int_t bin, Int_t& aLow, Int_t& aHigh, Int_t& bLow, Int_t& bHigh) const;
    /**
     * Finish cumulants calculations. Takes input and
     * output lists in case Terminate is called separately
     * 
     * @param inlist List with input histograms
     * @param outlist List with output histograms
     */
    void CumulantsTerminate(TList* inlist, TList* outlist);
    /*
     * Enumeration for cumulant histograms
     */
    enum { kW2Two = 1, // <w2*two>
	   kW2, // <w2>
	   kW4Four, // <w4*four>
	   kW4, // <w4>
	   kCosphi1phi2, // <cos(phi1+phi2)> 
	   kSinphi1phi2, // <sin(phi1+phi2)>
	   kCosphi1phi2phi3m, // <cos(phi1-phi2-phi3)>
	   kSinphi1phi2phi3m, // <sin(phi1-phi2-phi3)>
	   k3pWeight, // M(M-1)(M-1) or (mp*M-2mq)(M-1)
	   kCosphi1phi2phi3p, // <cos(phi1+phi2-phi3)>
	   kSinphi1phi2phi3p // <sin(phi1+phi2-phi3)>
	  };
  protected:
    /** 
     * Calculate reference flow
     *
     * @param cumu2h QC2 histos
     * @param cumu4h QC4 histos
     * @param quality QC Quality diag. histo
     * @param chist Centrality histogram
     * @param dNdetaRef dN/deta histogram
     */
    void CalculateReferenceFlow(CumuHistos& cumu2h, CumuHistos& cumu4h, TH2I* quality, TH1D* chist, TH2D* dNdetaRef) const;
    /** 
     * Calculate differential flow
     *
     * @param cumu2h QC2 histos
     * @param cumu4h QC4 histos
     * @param quality QC Quality diag. histo
     * @param dNdetaDiff dN/deta histogram
     */
    void CalculateDifferentialFlow(CumuHistos& cumu2h, CumuHistos& cumu4h, TH2I* quality, TH2D* dNdetaDiff) const;
    /** 
     * Calculate 3 correlator ref and fiff flow
     *
     * @param cumu2h QC2 histos
     * @param quality QC Quality diag. histo
     * @param chist Centrality histogram
     * @param dNdetaRef dN/deta histogram
     * @param dNdetaDiff dN/deta histogram
     */
    void Calculate3CorFlow(CumuHistos& cumu2h, TH2I* quality, TH1D* chist, TH2D* dNdetaRef, TH2D* dNdetaDiff) const;
    /**
     * Solve coupled eqs. to get v_n
     * 
     * @param cumu CumuHistos object with non-corrected flow results
     * @param type reference of differential flow ('r'/'d'/'a'/'b')
     */
    void SolveCoupledFlowEquations(CumuHistos& cumu, Char_t type) const;
    /**
     * Calculate NUA matrix elements to fill into the matrix
     * 
     * @param n row
     * @param m column
     * @param type reference of differential flow ('r'/'d'/'a'/'b')
     * @param binA Eta bin of phi1
     * @param cBin Centrality bin
     *
     * @return maxtrix element
     */
    Double_t CalculateNUAMatrixElement(Int_t n, Int_t m, Char_t type, Int_t binA, Int_t cBin) const;
    /**
     * Adds up the vertex bins to master profiles
     *
     * @param cumu QC histos
     * @param list output list
     * @param nNUA number of nua calculations
     */
    void AddVertexBins(CumuHistos& cumu, TList* list, UInt_t nNUA) const;
    /**
     * Get the bin number of <<cos(nphi)>>
     *
     * @param n moment
     *
     * @return bin number
     */
    Int_t GetBinNumberCos(Int_t n = 0) const;
    /**
     * Get the bin number of <<sin(nphi)>>
     *
     * @param n moment
     *
     * @return bin number
     */
    Int_t GetBinNumberSin(Int_t n = 0) const;
    /**
     * Setup NUA axis with labels
     *
     * @param a NUA axis
     */
    void SetupNUALabels(TAxis* a) const;
    /**
     * Make diagnostics hitogram
     *
     * @param name Name
     *
     * @return hist
     */
    TH2I* MakeQualityHist(const Char_t* name) const;
    /**
     * Make output histogram
     *
     * @param qc   # of particle correlations
     * @param n    flow moment
     * @param ctype  Type of flow
     * @param nua  For nua corrected hists
     *
     * @return hist
     */
    TH2D* MakeOutputHist(Int_t qc, Int_t n, const Char_t* ctype, UInt_t nua) const;

    UShort_t   fMaxMoment;     // Max flow moment 
    Int_t      fVzMin;         // z-vertex min must be in whole [cm]
    Int_t      fVzMax;         // z-vertex max must be in whole [cm]
    TString    fType;          // Data type
    UShort_t   fFlags;         // Flow flags, e.g., eta-gap or sat. vtx
    Double_t   fSigmaCut;      // Detector specific cut for outlier events
    Double_t   fEtaGap;        // Eta gap value
    Double_t   fEtaLims[6];    // Limits for binning in 3Cor method
    TH2D*      fCumuRef;       // Histogram for reference flow
    TH2D*      fCumuDiff;      // Histogram for differential flow
    CumuHistos fCumuHists;     // Array of histograms for cumulants calculations
    TH3D*      fCumuNUARef;    // histogram for NUA terms
    TH3D*      fCumuNUADiff;   // histogram for NUA terms
    TH2F*      fdNdedpRefAcc;  // Diagnostics histogram for acc. maps
    TH2F*      fdNdedpDiffAcc; // Diagnostics histogram for acc. maps
    TH2F*      fOutliers;      // Sigma <M> histogram 
    UShort_t   fDebug;         // Debug flag

    // ClassDef(VertexBin, 4); // object for eta dependent cumulants ananlysis
  };
  // ---------- End of nested class -------------
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o);
  /** 
   * Assignment operator 
   * 
   * @return Reference to this object 
   */
  AliForwardFlowTaskQC& operator=(const AliForwardFlowTaskQC&);
  /**
   * Initiate vertex bin objects
   */
  virtual void InitVertexBins();
  /**
   * Initiate diagnostics histograms
   */
  virtual void InitHists();
  /**
   * Analyze event
   *
   * @return true on success
   */
  virtual Bool_t Analyze();
  /**
   * Finalize analysis
   */
  virtual void Finalize();
  /**
   * Loops of vertex bins in list and runs analysis on those for current vertex
   *
   * @param list List of vertex bins
   * @param h1 dN/detadphi histogram
   * @param vtx Current vertex bin
   * @param flags Extra flags
   */
  void FillVtxBinList(const TList& list, TH2D& h1, Int_t vtx, UShort_t flags = 0x0) const;
  /**
   * Loops of vertex bins in list and runs analysis on those for current vertex
   *
   * @param list List of vertex bins
   * @param href dN/detadphi histogram for ref flow
   * @param hdiff dN/detadphi histogram for diff flow
   * @param vtx Current vertex bin
   * @param flags Extra flags
   */
  void FillVtxBinListEtaGap(const TList& list, TH2D& href, TH2D& hdiff, Int_t vtx, UShort_t flags = 0x0) const;
  /**
   * Loops of vertex bins in list and runs analysis on those for current vertex
   *
   * @param list List of vertex bins
   * @param hcent dN/detadphi histogram for central barrel
   * @param hfwd dN/detadphi histogram for fwd detectors
   * @param vtx Current vertex bin
   * @param flags Extra flags
   */
  void FillVtxBinList3Cor(const TList& list, TH2D& hcent, TH2D& hfwd, Int_t vtx, UShort_t flags = 0x0);
  /** 
   * Combine forward and central detector histograms to one histogram, to be used for 3 correlator method
   *
   * @param hcent Central data
   * @param hfwd Forward data
   *
   * @return combined hist
   */
  TH2D& CombineHists(TH2D& hcent, TH2D& hfwd);
  /**
   * Get Fill tracks from ESD or AOD input event
   *
   * @return true on success
   */
  Bool_t FillTracks(VertexBin* bin, UShort_t mode) const;
  /**
   * Loops over VertexBin list and calls terminate on each
   *
   * @param list VertexBin list
   */
  void EndVtxBinList(const TList& list) const;
  /**
   * Projects a list of TH2D's with flow
   * results to TH1's in centrality bins
   *
   * @param list List of flow results
   */
  void MakeCentralityHists(TList* list) const;
  /**
   * Check AODevent object for trigger, vertex and centrality
   * uses aod header if object is null
   * returns true if event is OK
   *
   * @param aodfm AliAODForwardMult object
   * 
   * @return Bool_t 
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult* aodfm);
  /**
   * Check trigger from AODForwardMult object
   * uses aod header if object is null
   * returns true if offline trigger is present
   *
   * @param aodfm AliAODForwardMult object
   * 
   * @return Bool_t 
   */
  virtual Bool_t CheckTrigger(const AliAODForwardMult* aodfm) const;
  /**
   * Check for centrality in AliAODForwardMult object, 
   * uses aod header if object is null
   * if present return true - also sets fCent value
   *
   * @param aodfm AliAODForwardMult object
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetCentrality(const AliAODForwardMult* aodfm);
  /**
   * Check for vertex in AliAODForwardMult
   * uses aod header if object is null
   * returns true if in range of fVtxAXis, also sets fVtx value
   *
   * @param aodfm AliAODForwardMult object
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetVertex(const AliAODForwardMult* aodfm);
  /**
   * Get VZERO Data
   *
   * @return VZERO data object
   */
  AliVVZERO* GetVZERO() const;
  /**
   * Fill VZERO d^2N/detadphi hist
   *
   * @param vzero AliAODVZERO object
   */
  void FillVZEROHist(AliVVZERO* vzero);
  /**
   * Print the setup of the task
   */
  void PrintFlowSetup() const;

  TAxis*             fVtxAxis;       //  Axis to control vertex binning
  TAxis*             fCentAxis;      //  Axis to control centrality/multiplicity binning
  Double_t           fFMDCut;        //  FMD sigma cut for outlier events
  Double_t           fSPDCut;        //  SPD sigma cut for outlier events
  UShort_t           fFlowFlags;     //  Flow flags, e.g., eta-gap, sat. vtx.
  Double_t           fEtaGap;        //  Eta gap value
  TList              fBinsForward;   //  List with forward VertexBin objects 
  TList              fBinsCentral;   //  List with central VertexBin objects
  TList*             fSumList;       //  Sum list
  TList*             fOutputList;    //  Output list
  AliAODEvent*       fAOD;           //  AOD event
  AliAnalysisFilter* fTrackCuts;     //  ESD track cuts
  Int_t              fMaxMoment;     //  Calculate v_{n} flag
  Float_t            fVtx;           //  Z vertex bin
  Double_t           fCent;          //  Centrality
  TH2D               fHistdNdedpV0;  //  VZERO d^2N/detadphi histogram
  TH2D               fHistdNdedp3Cor;//  3 correlator d^2N/detadphi histogram
  TH2D*              fHistFMDSPDCorr;//  Diagnostics hist for multiplicity correlations between FMD and SPD
  TH1D*              fHistCent;      //  Diagnostics hist for centrality
  TH1D*              fHistVertexSel; //  Diagnostics hist for selected vertices
  TH1I*              fHistEventSel;  //  Diagnostics hist for event selection

  ClassDef(AliForwardFlowTaskQC, 5); // Analysis task for flow analysis
};

#endif
// Local Variables:
//   mode: C++ 
// End:
