//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDFLOWTASKQC_H
#define ALIFORWARDFLOWTASKQC_H
/**
 * @file AliForwardFlowTaskQC.h
 * @author Alexander Hansen
 * @date   Tue Feb 14 2012
 * 
 * @brief
 * 
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
class AliAODForwardMult;
class TH1D;
class TH2D;
class TH3D;
class TAxis;

 /**
 * @defgroup pwglf_forward_tasks_flow Flow tasks 
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
   * Loops over VertexBin list and calls terminate on each
   *
   * @param list VertexBin list
   */
  void EndVtxBinList(const TList& list) const;
  /**
   * Returns the outputlist
   * 
   * @return TList* 
   */
  TList* GetOutputList() { return fOutputList; }
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
  /*
   * Check for vertex in AliAODForwardMult
   * returns true if in range of fVtxAXis, also sets fVtx value
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetVertex(const AliAODForwardMult* aodfm);
  /**
   * Set which harmonics to calculate. @f$ v_{1}@f$ to @f$ v_{4}@f$ is
   * available and calculated as default
   * 
   * @param  v1 Do @f$ v_{1}$f$
   * @param  v2 Do @f$ v_{2}$f$
   * @param  v3 Do @f$ v_{3}$f$
   * @param  v4 Do @f$ v_{4}$f$
   * @param  v5 Do @f$ v_{5}$f$
   * @param  v6 Do @f$ v_{6}$f$
   * 
   * @return void 
   */
  void SetDoHarmonics(Bool_t v1 = kTRUE, Bool_t v2 = kTRUE, 
		      Bool_t v3 = kTRUE, Bool_t v4 = kTRUE,
		      Bool_t v5 = kTRUE, Bool_t v6 = kTRUE) { 
    fv[1] = v1; fv[2] = v2; fv[3] = v3; fv[4] = v4; fv[5] = v5; fv[6] = v6;}
   /*
   * Set non-default vertex binning and range
   *
   * @param axis Use this vtx axis
   *
   * @return void
   */
  void SetVertexAxis(TAxis* axis) { fVtxAxis = axis; }
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
     * @parma sym Data is symmetric in eta
     */
    VertexBin(Int_t vLow, Int_t vHigh, 
              UShort_t moment, TString type,
              Bool_t sym = kTRUE);
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
     * @param VertexBin&
     * 
     * @return VertexBin& 
     */
    VertexBin& operator=(const VertexBin&);
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
     *
     * @return false if bad event (det. hotspot)
     */
    Bool_t FillHists(const TH2D& dNdetadphi);
    /**
     * Do cumulants calculations for current event with 
     * centrality cent
     * 
     * @param cent Event centrality
     * 
     * @return void 
     */
    void CumulantsAccumulate(Double_t cent);
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
    enum { kHmult = 1, kHQnRe, kHQnIm, kHQ2nRe, kHQ2nIm };
    /*
     * Enumeration for cumulant histograms
     */
    enum { kW2Two = 1, 
	   kW2, 
	   kW4Four, 
	   kW4, 
	   kQnRe, 
	   kQnIm, 
	   kM,
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

    const UShort_t fMoment;    // flow moment 
    const Int_t    fVzMin;     // z-vertex min must be in whole [cm]
    const Int_t    fVzMax;     // z-vertex max must be in whoe [cm]
    TString        fType;      // data type
    const Bool_t   fSymEta;    // Use forward-backward symmetry, if detector allows it
    TH2D*          fCumuRef;   // histogram for reference flow
    TH2D*          fCumuDiff;  // histogram for differential flow
    TH3D*          fCumuHist;  // histogram for cumulants calculations
    TH2D*          fdNdedpAcc; // Diagnostics histogram to make acc. maps
    UShort_t       fDebug;     // Debug flag

    ClassDef(VertexBin, 1); // object for cumulants ananlysis in FMD
  };

protected:
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
  /*
   * Initiate vertex bin objects
   */
  virtual void InitVertexBins();
  /*
   * Initiate diagnostics histograms
   */
  virtual void InitHists();
  /*
   * Analyze event
   */
  virtual Bool_t Analyze();
  /*
   * Finalize analysis
   */
  virtual void Finalize();

  TAxis*         fVtxAxis;        //  Axis to control vertex binning
  TList          fBinsFMD;        //  list with FMD VertexBin objects 
  TList          fBinsSPD;        //  list with SPD VertexBin objects
  TList*         fSumList;        //  sum list
  TList*         fOutputList;     //  Output list
  AliAODEvent*   fAOD;            //  AOD event
  Bool_t         fv[7];           //  Calculate v_{n} flag
  Float_t  	 fVtx;            //  Z vertex bin
  Double_t       fCent;           //  Centrality
  TH1D*          fHistCent;       //  Diagnostics hist for centrality
  TH1D*          fHistVertexSel;  //  Diagnostics hist for selected vertices
  TH1D*          fHistVertexAll;  //  Diagnostics hist for all vertices

  ClassDef(AliForwardFlowTaskQC, 1); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
