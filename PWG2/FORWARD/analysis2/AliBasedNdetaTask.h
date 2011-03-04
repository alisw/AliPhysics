//
// Task to analyse the AOD for for dN/deta in the base regions 
//
#ifndef ALIBASEDNDETATASK_H
#define ALIBASEDNDETATASK_H
#include <AliAnalysisTaskSE.h>
// #include <AliAODBaseMult.h>
class TList;
class TH2D;
class TH1D;
class AliAODEvent;

/**
 * Task to determine the 
 */
class AliBasedNdetaTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   */
  AliBasedNdetaTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   * @param maxVtx  Set @f$v_z@f$ range
   */
  AliBasedNdetaTask(const char* name);

  /** 
   * Set the vertex range to use 
   * 
   * @param min Minimum (in centermeter)
   * @param max Maximum (in centermeter)
   */  
  void SetVertexRange(Double_t min, Double_t max) { fVtxMin=min; fVtxMax=max; }
  /** 
   * Set the rebinning factor 
   * 
   * @param rebin Rebinning factor 
   */
  void SetRebinning(Int_t rebin) { fRebin = rebin; }
  /** 
   * Set the trigger maskl 
   * 
   * @param mask Trigger mask
   */
  void SetTriggerMask(UShort_t mask) { fTriggerMask = mask; }
  /** 
   * Set the trigger mask 
   * 
   * @param mask trigger mask 
   */
  void SetTriggerMask(const char* mask);
  /** 
   * Trigger efficiency for selected trigger(s)
   * 
   * @param e Trigger efficiency 
   */
  void SetTriggerEff(Double_t e) { fTriggerEff = e; } 
  /** 
   * Set the shape correction (a.k.a., track correction) for selected
   * trigger(s)
   * 
   * @param h Correction
   */
  void SetShapeCorrection(const TH1* h);
  /**
   * Destructor
   * 
   */
  virtual ~AliBasedNdetaTask();
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
  AliBasedNdetaTask(const AliBasedNdetaTask&);
  AliBasedNdetaTask& operator=(const AliBasedNdetaTask&) { return *this; }

  /** 
   * Retrieve the histogram 
   * 
   * @param aod AOD event 
   * @param mc  Whether to get the MC histogram or not
   * 
   * @return Retrieved histogram or null
   */
  virtual TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc=false) = 0;
  /** 
   * Check the trigger and vertex 
   * 
   * @param aod 
   * 
   * @return 
   */
  Bool_t CheckEvent(const AliAODEvent* aod);
  /** 
   * Clone a 2D histogram
   * 
   * @param in    Histogram to clone.
   * @param name  New name of clone.
   * 
   * @return The clone
   */
  TH2D* CloneHist(const TH2D* in, const char* name);
  /** 
   * Make a copy of the input histogram and rebin that histogram
   * 
   * @param h  Histogram to rebin
   * 
   * @return New (rebinned) histogram
   */
  TH1D* Rebin(const TH1D* h) const;
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  TH1* Symmetrice(const TH1* h) const;
  /** 
   * Project onto the X axis 
   * 
   * @param h         2D histogram 
   * @param name      New name 
   * @param firstbin  First bin to use 
   * @param lastbin   Last bin to use
   * @param error     Whether to calculate errors
   * 
   * @return Newly created histogram or null
   */
  TH1D* ProjectX(const TH2D* h, 
		 const char* name,
		 Int_t firstbin, 
		 Int_t lastbin, 
		 bool  corr=true,
		 bool  error=true) const;
  /** 
   * Set histogram graphical options, etc. 
   * 
   * @param h       Histogram to modify
   * @param colour  Marker color 
   * @param marker  Marker style
   * @param title   Title of histogram
   * @param ytitle  Title on y-axis. 
   */
  void  SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker, 
			       const char* title, 
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
    kAccepted   = 9,
    kMCNSD      = 10
  };

  TH2D*           fSum;          // Sum of histograms 
  TH2D*           fSumMC;        // Sum of MC histograms (if any)

  TList*          fSums;         // Container of sums 
  TList*          fOutput;       // Container of outputs 

  TH1D*           fTriggers;     // Histogram of triggers 
   
  Double_t        fVtxMin;       // Minimum v_z
  Double_t        fVtxMax;       // Maximum v_z
  Int_t           fTriggerMask;  // Trigger mask 
  Int_t           fRebin;        // Rebinning factor 
  Bool_t          fCutEdges;     // Whether to cut edges when rebinning
  Bool_t          fSymmetrice;   // Whether to symmetrice data 
  Bool_t          fCorrEmpty;    // Correct for empty bins 
  Double_t        fTriggerEff;   // Trigger efficiency for selected trigger(s)
  TH1*            fShapeCorr;    // Shape correction 

  ClassDef(AliBasedNdetaTask,1); // Determine multiplicity in base area
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
