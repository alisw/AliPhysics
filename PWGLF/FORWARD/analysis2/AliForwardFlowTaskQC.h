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
class AliAODForwardMult;
class TH1D;
class TH2D;
class TH3D;

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
   * Returns the outputlist
   * 
   * @return TList* 
  */
  TList* GetOutputList() { return fOutputList; }
 /**
  * Check AODForwardMult object for trigger, vertex and centrality
  * returns true if event is OK
  * 
  * @param const aodfm
  * 
  * @return Bool_t 
  */
  Bool_t AODCheck(const AliAODForwardMult* aodfm);
   /*
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
    */
    VertexBin(const Int_t vLow, const Int_t vHigh, 
              const Int_t moment, const Char_t* type);
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
    * Check if vertex vZ is in range for this bin
    * 
    * @param vZ z-coordinate of vertex to check
    * 
    * @return Bool_t 
    */
    Bool_t CheckVertex(Double_t vZ);
  /**
    * Fill reference and differential flow histograms for analysis
    *
    * @param dNdetadphi 2D data histogram
    *
    * @return false if bad event (det. hotspot)
    */
    Bool_t FillHists(TH2D* dNdetadphi);
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
    * Enumeration for cumulant histogram
    */
    enum { kW2Two = 1, kW2, kW4Four, kW4, kQnRe, kQnIm, kM,
       kCosphi1phi2, kSinphi1phi2, kCosphi1phi2phi3m, kSinphi1phi2phi3m, kMm1m2, 
       kw2two, kw2, kw4four, kw4, kpnRe, kpnIm, kmp, 
       kCospsi1phi2, kSinpsi1phi2, kCospsi1phi2phi3m, kSinpsi1phi2phi3m,
       kmpmq, kCospsi1phi2phi3p, kSinpsi1phi2phi3p };

    const Int_t   fMoment;        // flow moment 
    const Int_t   fVzMin;         // z-vertex min
    const Int_t   fVzMax;         // z-vertex max
    const Char_t* fType;          // data type
    TH2D*         fCumuRef;       // histogram for reference flow
    TH2D*         fCumuDiff;      // histogram for differential flow
    TH3D*         fCumuHist;      // histogram for cumulants calculations
    TH2D*         fdNdedpAcc;     // Diagnostics histogram to make acc. maps

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

  TList          fBinsFMD;        //  list with FMD VertexBin objects 
  TList          fBinsSPD;        //  list with SPD VertexBin objects
  TList*         fSumList;        //  sum list
  TList*         fOutputList;     //  Output list
  AliAODEvent*   fAOD;            //  AOD event
  Bool_t         fv[7];           //  Calculate v_{n} flag
  Float_t  	 fZvertex;        //  Z vertex bin
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
