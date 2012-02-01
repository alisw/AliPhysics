/**
 * @file   QAStructs.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 12:20:49 2011
 * 
 * @brief  Data structures used in QA trending tree
 *
 * @ingroup pwg2_forward_qa_scripts
 */
#ifndef QASTRUCTS_H
#define QASTRUCTS_H
#ifndef __CINT__
# include <TTree.h>
# include <TString.h>
#else 
class TTree;
#endif

// --- A quantity ----------------------------------------------------
/**
 * A simple quantity with mean, variance, min, and max
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct Quantity
{
  Double_t mean;  // Mean
  Double_t var;   // Variance 
  Double_t min;   // Minimum
  Double_t max;   // Maximum
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Quantity* MakeBranch(TTree* tree, const char* name)
  {
    Quantity* q = new Quantity;
    if (tree) tree->Branch(name, q, "mean/D:var:min:max");
    return q;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Quantity* SetBranch(TTree* tree, const char* name)
  {
    Quantity* q = new Quantity;
    if (tree) tree->SetBranchAddress(name, q);
    return q;    
  }
};

// --- A quantity ----------------------------------------------------
/**
 * A per-ring quantity with mean, variance, min, and max
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct RingQuantity : public Quantity 
{
  UShort_t det;
  Char_t   ring;

  /** 
   * Return the branch name to use 
   * 
   * @param d Detector 
   * @param r Ring
   * 
   * @return Name of branch to use 
   */
  static const char* BranchName(UShort_t d, Char_t r, const char* name)
  {
    return Form("FMD%d%c_%s", d, r, name);
  }
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static RingQuantity* MakeBranch(TTree* tree, UShort_t d, Char_t r, 
				  const char* name)
  {
    RingQuantity* q = new RingQuantity;
    if (tree) 
      tree->Branch(BranchName(d, r, name), q, "mean/D:var:min:max:d/s:r/b");
    q->det = d;
    q->ring = r;
    return q;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static RingQuantity* SetBranch(TTree* tree, UShort_t d, Char_t r, 
				 const char* name)
  {
    RingQuantity* q = new RingQuantity;
    if (tree) tree->SetBranchAddress(BranchName(d, r, name), q);
    return q;    
  }
};

// --- A quantity ----------------------------------------------------
/**
 * Per run information
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct Global 
{
  UInt_t   runNo;
  UInt_t   nAccepted;
  Double_t meanVz;
  Double_t sigmaVz;
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Global* MakeBranch(TTree* tree)
  {
    Global* g = new Global;
    if (tree) tree->Branch("global", g, "run/i:accepted:meanVz/D:sigmaVz");
    return g;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Global* SetBranch(TTree* tree) 
  {
    Global* g = new Global;
    if (tree) tree->SetBranchAddress("global", g);
    return g;
  }
};

// --- A quantity ----------------------------------------------------
/**
 * Per-ring status information on ELoss fits
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct FitStatus
{
  UShort_t nLow;
  UShort_t nCandidates;
  UShort_t nFitted;
  UShort_t det;
  Char_t   ring;

  /** 
   * Return the branch name to use 
   * 
   * @param d Detector 
   * @param r Ring
   * 
   * @return Name of branch to use 
   */
  static const char* BranchName(UShort_t d, Char_t r)
  {
    return Form("FMD%d%c_fitstatus", d, r);
  }
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static FitStatus* MakeBranch(TTree* tree, UShort_t d, Char_t r)
  {
    FitStatus* s = new FitStatus;
    if (tree) 
      tree->Branch(BranchName(d, r), s, "low/s:candidates:fitted:d:r/b");
    s->det  = d;
    s->ring = r;
    return s;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static FitStatus* SetBranch(TTree* tree, UShort_t d, Char_t r)
  {
    FitStatus* s = new FitStatus;
    if (tree) tree->SetBranchAddress(BranchName(d, r), s);
    return s;
  }    
};


// --- A quantity ----------------------------------------------------
/**
 * Per-ring status information on merging (sharing filter)
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct Merge
{
  Double_t one;
  Double_t two;
  Double_t three;
  UShort_t det;
  Char_t   ring;

  /** 
   * Return the branch name to use 
   * 
   * @param d Detector 
   * @param r Ring
   * 
   * @return Name of branch to use 
   */
  static const char* BranchName(UShort_t d, Char_t r)
  {
    return Form("FMD%d%c_merge", d, r);
  }
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Merge* MakeBranch(TTree* tree, UShort_t d, Char_t r)
  {
    Merge* s = new Merge;
    if (tree) tree->Branch(BranchName(d, r), s, "one/D:two:three:d:r/b");
    s->det  = d;
    s->ring = r;
    return s;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Merge* SetBranch(TTree* tree, UShort_t d, Char_t r)
  {
    Merge* s = new Merge;
    if (tree) tree->SetBranchAddress(BranchName(d, r), s);
    return s;
  }    
};

// --- Cuts ----------------------------------------------------------
/**
 * Per-ring status information on hit 'loss'
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct DataLoss
{
  Double_t merge;
  Double_t density;
  Double_t full;
  UShort_t det;
  Char_t   ring;

  /** 
   * Return the branch name to use 
   * 
   * @param d Detector 
   * @param r Ring
   * 
   * @return Name of branch to use 
   */
  static const char* BranchName(UShort_t d, Char_t r)
  {
    return Form("FMD%d%c_dataloss", d, r);
  }
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static DataLoss* MakeBranch(TTree* tree, UShort_t d, Char_t r)
  {
    DataLoss* s = new DataLoss;
    if (tree) tree->Branch(BranchName(d, r), s, "merge/D:density:full:d:r/b");
    s->det  = d;
    s->ring = r;
    return s;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static DataLoss* SetBranch(TTree* tree, UShort_t d, Char_t r)
  {
    DataLoss* s = new DataLoss;
    if (tree) tree->SetBranchAddress(BranchName(d, r), s);
    return s;
  }    
};

// --- A quantity ----------------------------------------------------
/**
 * Per-ring status information on Poisson vs ELoss correlation
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
struct Correlation
{
  Double_t alpha;
  Double_t beta;
  Double_t a;
  Double_t ea;
  Double_t b;
  Double_t eb;
  Double_t chi2;
  UShort_t det;
  Char_t   ring;
  /** 
   * Return the branch name to use 
   * 
   * @param d Detector 
   * @param r Ring
   * 
   * @return Name of branch to use 
   */
  static const char* BranchName(UShort_t d, Char_t r)
  {
    return Form("FMD%d%c_correlation", d, r);
  }
  /** 
   * Make a branch to contain an object of this type
   * 
   * @param tree Tree to store in
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Correlation* MakeBranch(TTree* tree, UShort_t d, Char_t r)
  {
    Correlation* s = new Correlation;
    if (tree) tree->Branch(BranchName(d, r), s, 
			   "alpha/D:beta:a:ea:b:eb:chi2:d/s:r/b");
    s->det  = d;
    s->ring = r;
    return s;
  }
  /** 
   * Associate a branch with an object of this type
   * 
   * @param tree Tree to read from
   * @param name Name of branch
   * 
   * @return Newly allocated object of this type associated with the
   * branch 
   */  
  static Correlation* SetBranch(TTree* tree, UShort_t d, Char_t r)
  {
    Correlation* s = new Correlation;
    if (tree) tree->SetBranchAddress(BranchName(d, r), s);
    return s;
  }    
};
#endif // QASTRUCTS_H
// Local Variables:
//  mode: C++
// End:
