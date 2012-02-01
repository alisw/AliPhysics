/**
 * @file   QARing.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 12:16:51 2011
 * 
 * @brief  Base class for QA rings
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
#ifndef QARING_H
#define QARING_H
#ifndef __CINT__
# include "QAStructs.h"
#else 
class RingQuantity;
class FitStatus;
class Merge;
class DataLoss;
class Correlation;
class TTree;
#endif

/**
 * Base class for QA rings.  This manages the per-ring entries in the
 * tree 
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct QARing 
{ 
  /** 
   * Constructor 
   * 
   * @param d  Detector
   * @param r  Ring
   * 
   * @return 
   */
  QARing(UShort_t d, Char_t r) 
  : fD(d), fR(r),
    fChi2(0), 
    fC(0), 
    fDelta(0), 
    fXi(0), 
    fSigma(0),
    fFitStatus(0),
    fMerge(0),
    fDataLoss(0),
    fOccupancy(0),
    fCorrelation(0)
  {}
  /** 
   * Destructor 
   */
  virtual ~QARing()
  {
    if (fChi2)		{ delete fChi2;		fChi2 		= 0; } 
    if (fC)		{ delete fC;		fC 		= 0; } 
    if (fDelta)		{ delete fDelta;	fDelta 		= 0; } 
    if (fXi)		{ delete fXi;		fXi 		= 0; } 
    if (fSigma)		{ delete fSigma;	fSigma 		= 0; }
    if (fFitStatus)	{ delete fFitStatus;	fFitStatus 	= 0; }
    if (fMerge)		{ delete fMerge;	fMerge 		= 0; }
    if (fDataLoss)      { delete fDataLoss;     fDataLoss       = 0; }
    if (fOccupancy)	{ delete fOccupancy;	fOccupancy 	= 0; }
    if (fCorrelation)	{ delete fCorrelation;	fCorrelation 	= 0; } 
  }
  /** 
   * Initialise 
   * 
   * @param tree Tree to read from/write to
   * @param read If true, read from tree, other wise write 
   * 
   * @return true on success
   */
  Bool_t Init(TTree* tree, bool read=false)
  {
    if (!read) {
      fChi2        = RingQuantity::MakeBranch(tree, fD, fR, "chi2");
      fC           = RingQuantity::MakeBranch(tree, fD, fR, "c");
      fDelta       = RingQuantity::MakeBranch(tree, fD, fR, "delta");
      fXi          = RingQuantity::MakeBranch(tree, fD, fR, "xi");
      fSigma       = RingQuantity::MakeBranch(tree, fD, fR, "sigma");
      fFitStatus   = FitStatus::MakeBranch(tree, fD, fR);
      fMerge       = Merge::MakeBranch(tree, fD, fR);
      fDataLoss    = DataLoss::MakeBranch(tree, fD, fR);
      fOccupancy   = RingQuantity::MakeBranch(tree, fD, fR, "occupancy");
      fCorrelation = Correlation::MakeBranch(tree, fD, fR);
      return true;
    }
    fChi2        = RingQuantity::SetBranch(tree, fD, fR, "chi2");
    fC           = RingQuantity::SetBranch(tree, fD, fR, "c");
    fDelta       = RingQuantity::SetBranch(tree, fD, fR, "delta");
    fXi          = RingQuantity::SetBranch(tree, fD, fR, "xi");
    fSigma       = RingQuantity::SetBranch(tree, fD, fR, "sigma");
    fFitStatus   = FitStatus::SetBranch(tree, fD, fR);
    fMerge       = Merge::SetBranch(tree, fD, fR);
    fDataLoss    = DataLoss::SetBranch(tree, fD, fR);
    fOccupancy   = RingQuantity::SetBranch(tree, fD, fR, "occupancy");
    fCorrelation = Correlation::SetBranch(tree, fD, fR);
    return true;
  }
  UShort_t      fD;              // Detector number
  Char_t        fR;              // Ring identifier 
  RingQuantity* fChi2;           // ELoss reduced chi-square
  RingQuantity* fC;              // ELoss constant 
  RingQuantity* fDelta;          // ELoss most probably value 
  RingQuantity* fXi;             // ELoss Landau width
  RingQuantity* fSigma;          // ELoss Gaus width
  FitStatus*    fFitStatus;      // ELoss fit status  
  Merge*        fMerge;          // Merging parameters
  DataLoss*     fDataLoss;       // Hit 'loss' summary
  RingQuantity* fOccupancy;      // Mean occupancy 
  Correlation*  fCorrelation;    // ELoss vs Poisson correlation 
};

#endif // QARING_H
// Local Variables:
//  mode: C++
// End:
