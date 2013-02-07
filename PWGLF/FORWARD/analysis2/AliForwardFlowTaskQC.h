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
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TArrayI.h"
class AliAODForwardMult;
class TH1D;
class TH2F;
class TH2D;
class TH3D;
class TAxis;
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
 *   - AnalysisResults.root
 *
 * @ingroup pwglf_forward_tasks_flow
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
   * Set which flow moments to calculate.
   * 
   * @param  n Do @f$ v_{n}@f$
   * 
   * @return void 
   */
  void AddFlowMoment(Short_t n); 
  /**
   * Set non-default vertex binning and range
   *
   * @param axis Use this vtx axis
   */
  void SetVertexAxis(TAxis* axis) { fVtxAxis = axis; }
  /**
   * Set detector sigma cuts
   *
   * @param fmdCut FMD sigma cut
   * @param spdCut SPD sigma cut
   */
  void SetDetectorCuts(Double_t fmdCut, Double_t spdCut) { fFMDCut = fmdCut; fSPDCut = spdCut; }
  /**
   * Set flow flags, @f$\eta@f$-gap, sym. around @f$\eta=0@f$ or
   * sat. vtx. interactions
   *
   * @param flags EFlowFlags 
   */
  void SetFlowFlags(UShort_t flags) { fFlowFlags = flags; }
  /**
   * Enum for flow flags
   */
  enum EFlowFlags {
    kEtaGap  = 0x1,
    kSymEta  = 0x2,
    kSatVtx  = 0x4
  };
  /**
   *  Set @f$|\eta|@f$ value to make cut for @f$\eta@f$ gap at
   *
   * @param eg gap value
   */
  void SetEtaGapValue(Double_t eg) { fEtaGap = eg; }
protected:
  /**
   * Enum for reference flow (eta-gap) mode
   */
  enum EFillFlow {
    kFillRef  = 0x1,
    kFillDiff = 0x2,
    kFillBoth = 0x3
  };
  // ----------------- Being nested class ---------------------
  /**
   * Nested class to handle cumulant calculations in vertex bins
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
     * @param type Data type (FMD/SPD/FMDTR/SPDTR/MC)
     * @param flags Flags
     * @param cut Cut value 
     * @param etaGap @f$\eta@f$ gap 
     */
    VertexBin(Int_t vLow, Int_t vHigh, 
              UShort_t moment, TString type,
              UShort_t flags = kSymEta, 
              Double_t cut = -1, Double_t etaGap = 2.);
    /**
     * Copy constructor 
     * 
     * @param o Object to copy from
     * 
     * @return VertexBin
     */
    VertexBin(const VertexBin& o);
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
     * 
     * @return void 
     */
    virtual void AddOutput(TList* list);
    /**
     * Fill reference and differential flow histograms for analysis
     *
     * @param dNdetadphi 2D data histogram
     * @param cent Centrality
     * @param mode fill ref/diff or both
     *
     * @return false if bad event (det. hotspot)
     */
    Bool_t FillHists(const TH2D& dNdetadphi, Double_t cent, EFillFlow mode = kFillBoth);
    /**
     * Do cumulants calculations for current event with 
     * centrality cent
     * 
     * @param cent Event centrality
     * @param skipFourP Skip ?
     *
     * @return void 
     */
    void CumulantsAccumulate(Double_t cent, Bool_t skipFourP = kFALSE);
    /**
     * Finish cumulants calculations. Takes input and
     * output lists in case Terminate is called separately
     * 
     * @param inlist List with input histograms
     * @param outlist List with output histograms
     * 
     * @return void 
     */
    void CumulantsTerminate(TList* inlist, TList* outlist);

  protected:
    /*
     * Enumeration for ref/diff histograms
     */
    enum { kHmultA = 1, kHmultB, kHQnReA, kHQnImA, kHQnReB, kHQnImB, kHQ2nRe, kHQ2nIm };
    /*
     * Enumeration for cumulant histograms
     */
    enum { kW2Two = 1, 
	   kW2, 
	   kW4Four, 
	   kW4, 
	   kQnReA, 
	   kQnImA, 
	   kMA,
	   kQnReB, 
	   kQnImB, 
	   kMB,
	   kCosphi1phi2, 
	   kSinphi1phi2, 
	   kCosphi1phi2phi3m, 
	   kSinphi1phi2phi3m, 
	   kMm1m2, 
	   kw2two, 
	   kw2, 
	   kw4four, 
	   kw4, 
	   kpnRe, 
	   kpnIm, 
	   kmp, 
	   kCospsi1phi2, 
	   kSinpsi1phi2, 
	   kCospsi1phi2phi3m, 
	   kSinpsi1phi2phi3m,
	   kmpmq, 
	   kCospsi1phi2phi3p, 
	   kSinpsi1phi2phi3p };
    /**
     * Set centrality axis
     *
     * @param axis Centrality axis
     *
     * @return void
     */
    void SetupCentAxis(TAxis* axis);

    const UShort_t fMoment;        // flow moment 
    const Int_t    fVzMin;         // z-vertex min must be in whole [cm]
    const Int_t    fVzMax;         // z-vertex max must be in whoe [cm]
    TString        fType;          // data type
    const UShort_t fFlags;         // Flow flags, e.g., eta-gap sat. vtx
    const Double_t fSigmaCut;      // Detector specific cut for outlier events
    const Double_t fEtaGap;        // Eta gap value
    TH2D*          fCumuRef;       // histogram for reference flow
    TH2D*          fCumuDiff;      // histogram for differential flow
    TH3D*          fCumuHist;      // histogram for cumulants calculations
    TH2F*          fdNdedpAcc;     // Diagnostics histogram to make acc. maps
    TH2F*          fOutliers;      // Sigma <M> histogram 
    UShort_t       fDebug;         // Debug flag

    ClassDef(VertexBin, 3); // object for cumulants ananlysis in FMD
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
   * @param h dN/detadphi histogram
   * @param vtx Current vertex bin
   *
   * @return true on success
   */
  Bool_t FillVtxBinList(const TList& list, const TH2D& h, Int_t vtx) const;
  /**
   * Loops of vertex bins in list and runs analysis on those for current vertex
   *
   * @param list List of vertex bins
   * @param href dN/detadphi histogram for ref flow
   * @param hdiff dN/detadphi histogram for diff flow
   * @param vtx Current vertex bin
   *
   * @return true on success
   */
  Bool_t FillVtxBinListEtaGap(const TList& list, const TH2D& href, const TH2D& hdiff, Int_t vtx) const;

  /**
   * Loops over VertexBin list and calls terminate on each
   *
   * @param list VertexBin list
   */
  void EndVtxBinList(const TList& list) const;
  /**
   * Projects a list of TProfile2D's with flow
   * results to TH1's in centrality bins
   *
   * @param list List of flow results
   *
   * @return void
   */
  void MakeCentralityHists(TList* list);
  /**
   * Check AODevent object for trigger, vertex and centrality
   * returns true if event is OK
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  Bool_t CheckEvent(const AliAODForwardMult* aodfm);
  /**
   * Check trigger from AODForwardMult object
   * returns true if offline trigger is present
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t CheckTrigger(const AliAODForwardMult* aodfm) const;
  /**
   * Check for centrality in AliAODForwardMult object, 
   * if present return true - also sets fCent value
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetCentrality(const AliAODForwardMult* aodfm);
  /**
   * Check for vertex in AliAODForwardMult
   * returns true if in range of fVtxAXis, also sets fVtx value
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetVertex(const AliAODForwardMult* aodfm);
  /**
   * Make diagnostics hitogram
   *
   * @return void
   */
  void MakeQualityHist(const Char_t* name) const;
  /**
   * Print the setup of the task
   *
   * @return void
   */
  virtual void PrintFlowSetup() const;

  TAxis*         fVtxAxis;          //  Axis to control vertex binning
  Double_t       fFMDCut;           //  FMD sigma cut for outlier events
  Double_t       fSPDCut;           //  SPD sigma cut for outlier events
  UShort_t       fFlowFlags;        //  Flow flags, e.g., eta-gap, sat. vtx.
  Double_t       fEtaGap;           //  Eta gap value
  TList          fBinsFMD;          //  list with FMD VertexBin objects 
  TList          fBinsSPD;          //  list with SPD VertexBin objects
  TList*         fSumList;          //  sum list
  TList*         fOutputList;       //  Output list
  AliAODEvent*   fAOD;              //  AOD event
  TArrayI        fV;                //  Calculate v_{n} flag
  Float_t  	 fVtx;              //  Z vertex bin
  Double_t       fCent;             //  Centrality
  TH1D*          fHistCent;         //  Diagnostics hist for centrality
  TH1D*          fHistVertexSel;    //  Diagnostics hist for selected vertices

  ClassDef(AliForwardFlowTaskQC, 3); // Analysis task for FMD analysis
};

#endif
// Local Variables:
//   mode: C++ 
// End:
