#ifndef ALIFMDQACHECKER_H
#define ALIFMDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#include "AliQACheckerBase.h"
class TFile; 
class TH1F; 
class TH1I; 
class TFitResultPtr;

/** 
 * @class AliFMDQAChecker 
 * @brief Quality assurance checker for the FMD 
 */
class AliFMDQAChecker : public AliQACheckerBase 
{
public:
  /** Constructor */
  AliFMDQAChecker();
  /** Destructor */
  virtual ~AliFMDQAChecker() {}
  /** 
   * Member function called to do the actual checking
   * 
   * @param rv   Array of return values. 
   * @param what What to check 
   * @param list Array of arrays of histograms.  There's one arrat for
   *             each 'specie'
   * @param t    Reconstruction parameters - not used. 
   */
  void Check(Double_t* rv, AliQAv1::ALITASK_t what, 
	     TObjArray** list, const AliDetectorRecoParam* t);
  /** 
   * Make output images.  This is overridden relative to the base
   * class so that we can set the log(y) scale and put everything on
   * the same axis. 
   * 
   * @param list  List of specie array of histograms 
   * @param task  What to show 
   * @param mode  Mode 
   */
  void  MakeImage(TObjArray** list, 
		  AliQAv1::TASKINDEX_t task, 
		  AliQAv1::MODE_t mode);
  /** 
   * Set wether to scale the histograms to a common Y axis or not when
   * generating plots
   * 
   * @param on If true, do scale 
   */
  void SetDoScale(Bool_t on=true) { fDoScale = on; }
protected:
  // Return values 
  enum { 
    kOK,
    kProblem, 
    kBad, 
    kWhatTheFk
  };
  /** 
   * Check one histogram 
   * 
   * @param specie Event specie 
   * @param hist   Histogram to check 
   * 
   * @return 0 if all is good, increasing severity for increasingly
   * bad data
   */
  UShort_t CheckOne(AliQAv1::ALITASK_t          what,
		    AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  /** 
   * Check raw input.   If the reconstructor is enabled, then we try to fit 
   * 
   * @f[ 
   *   f(\Delta;\Delta_p,\xi,\sigma) = \int_{-\infty}^{\infty}d\Delta_p' 
   *     f_{L}(\Delta;\Delta_p',\xi) \times 
   *     e^{\frac{(\Delta-\Delta_p')^2}{\sigma^2}}
   * @f]
   * 
   * where @f$f_L@f$ is the Landau distribution, to the data.  The
   * quality is set according to the value of @f$\chi^2/\nu@f$.
   *
   * If no reconstructor is set, then simply check that the histogram
   * isn't empty
   * 
   * @param specie Event specie 
   * @param hist   Histogram to check 
   * 
   * @return 0 if all is good, increasing severity for increasingly
   * bad data
   */
  UShort_t CheckRaw(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  /** 
   * Check simulation output.  Does a simple test of whether the
   * histogram is empty or not.
   * 
   * @param specie Event specie 
   * @param hist   Histogram to check 
   * 
   * @return 0 if all is good, increasing severity for increasingly
   * bad data
   */
  UShort_t CheckSim(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  /** 
   * Check ESD.  Does a simple test of whether the histogram is empty
   * or not.
   * 
   * @param specie Event specie 
   * @param hist   Histogram to check 
   * 
   * @return 0 if all is good, increasing severity for increasingly
   * bad data
   */
  UShort_t CheckESD(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  /** 
   * Check reconstruction points.  Does a simple test of whether the
   * histogram is empty or not.
   * 
   * @param specie Event specie 
   * @param hist   Histogram to check 
   * 
   * @return 0 if all is good, increasing severity for increasingly
   * bad data
   */
  UShort_t CheckRec(AliRecoParam::EventSpecie_t specie, 
		    TH1*                        hist) const;
  /** 
   * Set the returned QA from this checker based on the values in the
   * array @a values.  Note, this by-passes the Low/High setting of
   * the base class (which are very confusing) 
   * 
   * @param index   Task index 
   * @param values  Array of values 
   */
  void SetQA(AliQAv1::ALITASK_t index, Double_t* values) const;
  /** 
   * Process external parameters 
   * 
   */
  void ProcessExternalParams();
  /** 
   * Process a single external parameter 
   * 
   * @param name Name of parameter 
   * @param v    On return, the value - as a double 
   */
  void ProcessExternalParam(const char* name, Double_t& v);
  /** 
   * Get the thresholds from OCDB
   * 
   */
  void GetThresholds();
  /** 
   * The basic check on a histogram 
   * 
   * @param hist Histogram 
   * 
   * @return 1 empty - 0 otherwise 
   */
  UShort_t BasicCheck(TH1* hist) const;
  /** 
   * Translate our internal quality measure to QA framework bit 
   * 
   * @param qual Internal quality 
   * 
   * @return QA framework quality bit 
   */
  AliQAv1::QABIT_t Quality2Bit(UShort_t qual) const;
  /** 
   * Add fit results to to plot 
   * 
   * @param hist  Histogram 
   * @param res   Fit result 
   * @param color Color to use for the text - depends on quality 
   * @param low   Lower bound on fit range 
   * @param high  Upper bound on fit range 
   */
  void AddFitResults(TH1* hist, const TFitResultPtr& res, Int_t color,
		     Double_t low, Double_t high) const;
  UShort_t CheckFit(TH1* hist, const TFitResultPtr& res, 
		    Double_t low, Double_t high, Int_t& color) const;
  Bool_t   fDoScale;           // Whether to scale all histograms 
  Bool_t   fDidExternal;       // Whether we've processed the external params 
  Bool_t   fShowFitResults;    // Whether to put the fit result on the plots
  Double_t fELossLowCut;       // Low cut on ELoss fits 
  Double_t fELossNRMS;         // Number of RMS to fit upward
  Double_t fELossBadChi2Nu;    // Cut on bad chi2/nu
  Double_t fELossFkupChi2Nu;   // Cut on F**ked up chi2/nu
  Int_t    fELossMinEntries;   // Least number of entries before fitting
  Double_t fELossGoodParError; // Least relative error
  Double_t fROErrorsBad;       // Cut on read-out errors 
  Double_t fROErrorsFkup;      // Cut on read-out errors 
private:
  /** 
   * Copy constructor - not implemented 
   * 
   * @param qac Object to copy from
   */
  AliFMDQAChecker(const AliFMDQAChecker& qac); 
  /** 
   * assignment operator - not implemented 
   * 
   * @param qac Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliFMDQAChecker& operator=(const AliFMDQAChecker& qac); 

  ClassDef(AliFMDQAChecker,0)  // Yves? what to do? 
};

#endif // AliFMDQAChecker_H
// Local Variables:
//  mode: c++
// End:
