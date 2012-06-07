// This task finds the eventplane
// using the FMD
// 
#ifndef ALIFMDEVENTPLANEFINDER_H 
#define ALIFMDEVENTPLANEFINDER_H 
/**
 * @file AliFMDEventPlaneFinder.h
 * @author Alexander Hansen
 * @date   Tue Feb 14 2012
 * 
 * @brief
 * 
 * 
 * @ingroup pwglf_forward
 */
#include <TNamed.h>
#include <TVector2.h>
#include "AliForwardUtil.h"
class AliVEvent;
class TH1D;
class TH2F;
class TH2D;
class TString;
class AliOADBContainer;
class AliAODForwardEP;

/**
 * Find the event plane using the FMD
 * 
 */
class AliFMDEventPlaneFinder : public TNamed
{
public:
  /** 
   * Constructor 
   */
  AliFMDEventPlaneFinder();
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDEventPlaneFinder(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDEventPlaneFinder(const AliFMDEventPlaneFinder& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDEventPlaneFinder();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDEventPlaneFinder& operator=(const AliFMDEventPlaneFinder& o);
  /** 
   * Initialize this sub-algorithm
   * 
   * @param etaAxis  Eta axis to use 
   */
  virtual void Init(const TAxis& etaAxis);
  /** 
   * Do the calculations 
   * 
   * @param hists    Histogram cache
   * @param esd      Event 
   * @param aodEp    Output object 
   * @param h        Output histogram
   * 
   * @return true on successs 
   */
  Bool_t FindEventplane(AliVEvent* esd,
                        AliAODForwardEP& aodEp,
                        TH2D* h,
                        AliForwardUtil::Histos* hists);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  virtual void DefineOutput(TList* dir);
  /** 
   * Print information 
   * 
   * @param option Print options 
   *   - max  Print max weights 
   */
  void Print(Option_t* option="") const;
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /**
   * Calculate Q vectors
   *
   * @param h dN/detadphi histogram
   * @param eHist histogram for ep vs. eta
   */
  void CalcQVectors(TH2D* h, TH1D* eHist);
  /**
   * Calculate the eventplane from a vector
   *
   * @param v TVector2 of Q-vectors
   *
   * @return the eventplane as a double
   */
  Double_t CalcEventplane(const TVector2& v) const;
  /**
   * Set the run number, used for OADB object
   *
   * @param run Run number
   */
  void SetRunNumber(Int_t run);
  /**
   * Get the run number
   *
   * @return returns the run number
   */
  Int_t GetRunNumber() { return fRunNumber; }
  /**
   * Get the OADB phi distribution for flattening
   */
  void GetPhiDist();
  /**
   * Flag for setting the use of phi weights for flattening
   * 
   * @param use true or false
   */
  void SetUsePhiWeights(Bool_t use = kTRUE) { fUsePhiWeights = use; }
  /**
   * Fill diagnostics hists
   *
   * @param fmdEP Object containing results of FMD EP calculations
   */
  void FillHists(AliAODForwardEP* fmdEP);
  /**
   * Set the OADB path, for using a custom OADB path and file
   *
   * @param fname Name of the custom OADB file, including path
   */
  void SetOADBPath(Char_t* fname) { fOADBFileName = fname; }

protected:
  /**
   * Get the phi weight from OADB histogram for the ep flattening
   *
   * @param etaBin which eta bin
   * @param phiBin which phi bin
   *
   * @return phi weight for etaBin, phiBin as double
   */
  Double_t GetPhiWeight(Int_t etaBin, Int_t phiBin) const;
  /** 
   * Calculat the difference @f$a_1 - a_2@$f between two angles
   * @f$a_1, a_2@f$ and normalize to @f$[-\pi/2,\pi/2]@f$
   * 
   * @param a1 First angle @f$a_1@f$
   * @param a2 Second angle @f$a_2@f$
   * 
   * @return @f$a_1 - a_2 \in[-\pi/2,\pi/2]@f$
   */
  Double_t CalcDifference(Double_t a1, Double_t a2) const;
  /** 
   * Make a histogram of @f$\Psi_R@f$ values 
   * 
   * @param name   Name of histogram
   * @param title  Source 
   * @param color  Color of histogram 
   * 
   * @return Newly allocated histogram
   */
  TH1D* MakePsiRHist(const char* name, 
		     const char* title,
		     Int_t color);
  /** 
   * Make a difference histogram, and add to the output list
   * 
   * @param name    Name of histogram
   * @param first   First variable to correlate (X axis)
   * @param second  Second variable to correlate (Y axis)
   * @param color   Fill and line color of histogram
   * 
   * @return Newly allocated histogram
   */
  TH1D* MakeDiffHist(const char* name, 
		     const char* first,
		     const char* second,
		     Int_t color);
  /** 
   * Make a correlation histogram, and add it to the output list.
   * 
   * @param name    Name of histogram
   * @param first   First variable to correlate (X axis)
   * @param second  Second variable to correlate (Y axis)
   * 
   * @return Newly allocated histogram 
   */
  TH2F* MakeCorrHist(const char* name, 
		     const char* first,
		     const char* second);
  TList*            fList;               // List for diag. hists.
  AliVEvent*        fEvent;              // Current event
  TVector2          fQt;                 // Q vector for total ep
  TVector2          fQa;                 // Q vector for sub-ep A
  TVector2          fQc;                 // Q vector for sub-ep C
  TVector2          fQ1;                 // Q vector for sub-ep 1
  TVector2          fQ2;                 // Q vector for sub-ep 2
  TVector2          fQeta;               // Q vector for psi eta-dependence
  TH1D*             fHepFMD;             // Diagnostics histogram
  TH1D*             fHepFMDA;            // Diagnostics histogram
  TH1D*             fHepFMDC;            // Diagnostics histogram
  TH1D*             fHepFMDQC1;          // Diagnostics histogram
  TH1D*             fHepFMDQC2;          // Diagnostics histogram
  TH1D*             fHdiffFMDAC;         // Diagnostics histogram
  TH1D*             fHdiffFMDTPC;        // Diagnostics histogram
  TH1D*             fHdiffFMDVZERO;      // Diagnostics histogram
  TH2F*             fHcorrFMDAC;         // Diagnostics histogram
  TH2F*             fHcorrFMDTPC;        // Diagnostics histogram
  TH2F*             fHcorrFMDVZERO;      // Diagnostics histogram
  TH2D*             fHPhi;               // Diagnostics histogram
  Int_t             fDebug;              // Debug flag
  TString           fOADBFileName;       // Path to OADB container
  AliOADBContainer* fOADBContainer;      // OADBContainer object
  TH2D*             fPhiDist;            // Phi dist. for phi weights
  Int_t             fRunNumber;          // Run number supplied
  Bool_t            fUsePhiWeights;      // Flag for phi weights

  ClassDef(AliFMDEventPlaneFinder,1); //  
};

#endif
// Local Variables:
//   mode: C++
// End:

