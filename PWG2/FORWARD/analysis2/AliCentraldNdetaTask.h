//
// Task to analyse the AOD for for dN/deta in the central regions 
//
#ifndef ALICENTRALDNDETATASK_H
#define ALICENTRALDNDETATASK_H
#include <AliAnalysisTaskSE.h>
// #include <AliAODCentralMult.h>
class TList;
class TH2D;
class TH1D;

/**
 * Task to determine the 
 */
class AliCentraldNdetaTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   */
  AliCentraldNdetaTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   * @param maxVtx  Set @f$v_z@f$ range
   */
  AliCentraldNdetaTask(const char* name);
  
  void SetVertexRange(Double_t min, Double_t max) { fVtxMin=min; fVtxMax=max; }
  void SetRebinning(Int_t rebin) { fRebin = rebin; }
  void SetTriggerMask(UShort_t mask) { fTriggerMask = mask; }
  void SetTriggerMask(const char* mask);
  /**
   * Destructor
   * 
   */
  virtual ~AliCentraldNdetaTask();
  /** 
   * Initialise on master - does nothing
   * 
   */
  virtual void   Init() {}
  /** 
   * Create output objects.  
   *
   * This is called once per slave process 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process a single event 
   * 
   * @param option Not used
   */
  virtual void UserExec(Option_t* option);
  /** 
   * Called at end of event processing.. 
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
protected:
  AliCentraldNdetaTask(const AliCentraldNdetaTask&);
  AliCentraldNdetaTask& operator=(const AliCentraldNdetaTask&) { return *this; }
  /** 
   * Clone a 2D histogram
   * 
   * @param in    Histogram to clone.
   * @param name  New name of clone.
   * 
   * @return The clone
   */
  TH2D* CloneHist(TH2D* in, const char* name);
  /** 
   * Make a copy of the input histogram and rebin that histogram
   * 
   * @param h  Histogram to rebin
   * 
   * @return New (rebinned) histogram
   */
  TH1D* Rebin(const TH1D* h) const;
  /** 
   * Set histogram graphical options, etc. 
   * 
   * @param h       Histogram to modify
   * @param colour  Marker color 
   * @param marker  Marker style
   * @param title   Title of histogram
   * @param ytitle  Title on y-axis. 
   */
  void  SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker, const char* title, 
			       const char* ytitle="#frac{1}{N} #frac{dN_{ch}}{d#eta}");
  /** 
   * Trigger histogram bins 
   */
  enum { 
    kAll        = 1, 
    kB          = 2, 
    kA          = 3, 
    kC          = 4, 
    kE          = 5,
    kMB         = 6, 
    kWithTrigger= 7,
    kWithVertex = 8, 
    kAccepted   = 9 
  };

  TH2D*           fSumCentral;    //  Sum of histograms 
  TH2D*           fSumCentralMC;  //  Sum of MC histograms (if any)

  TList*          fSums;          // Container of sums 
  TList*          fOutput;        // Container of outputs 

  TH1D*           fTriggers;      // Histogram of triggers 
  
  Double_t        fVtxMin;        // Minimum v_z
  Double_t        fVtxMax;        // Maximum v_z
  Int_t           fTriggerMask;   // Trigger mask 
  Int_t           fRebin;         // Rebinning factor 
  Bool_t          fCutEdges;      // Whether to cut edges when rebinning

  ClassDef(AliCentraldNdetaTask,1); // Determine multiplicity in central area
};

#endif
