#ifndef ALIFORWARDMULTDIST_H
#define ALIFORWARDMULTDIST_H
#include <AliAnalysisTaskSE.h>
#include <TList.h>
#include <TString.h>
class TH1;
class TH2;
class AliAODForwardMult;

/**
 * Class to make raw @f$P(N_{ch})@f$ distributions 
 * 
 */
class AliForwardMultDists : public AliAnalysisTaskSE
{
public:
  enum { 
    kInvalidEta = 999
  };
  enum { 
    kAnalysis = 1, 
    kMC       = 2, 
    kTrigger  = 3, 
    kVertex   = 4, 
    kTriggerVertex = 5
  };
  /** 
   * Default constructor
   */
  AliForwardMultDists();
  /** 
   * User constructor 
   * 
   * @param name Name of the task 
   */
  AliForwardMultDists(const char* name);
  /** 
   * Copy constructor
   *
   * @param o object to copy fron
   */
  AliForwardMultDists(const AliForwardMultDists& o);
  /** 
   * Assignment operator
   * 
   * @param o object to assign from 
   *
   * @return Reference to this object
   */
  AliForwardMultDists& operator=(const AliForwardMultDists& o);
  /** 
   * Destructor
   */
  virtual ~AliForwardMultDists() {}
  /** 
   * Create output objects - called at start of job in slave 
   * 
   */
  void UserCreateOutputObjects();
  /** 
   * Analyse a single event 
   * 
   * @param option Not used
   */
  void UserExec(Option_t* option="");
  /** 
   * Called at the end of the final processing of the job on the
   * full data set (merged data)
   * 
   * 
   * @param option Not used
   */
  void Terminate(Option_t* option="");
  /** 
   * Set-up internal structures on first seen event 
   * 
   * @param hist Basic histogram template from AOD object 
   */
  void SetupForData(const TH2& hist, Bool_t useMC);
  void StoreInformation(const AliAODForwardMult* forward);
  /** 
   * Project a 2D histogram into a 1D histogram taking care to use
   * either the @f$\phi2f$ acceptance stored in the overflow bins, or
   * the @f$\eta@f$ coverage stored in the underflow bins.
   * 
   * @param input      2D histogram to project 
   * @param cache      1D histogram to project into 
   * @param usePhiAcc  If true, use the @f$\phi2f$ acceptance stored in
   * the overflow bins, or if false the @f$\eta@f$ coverage stored in
   * the underflow bins.
   */
  static void ProjectX(const TH2& input, TH1& cache, Bool_t usePhiAcc=true);
  /** 
   * Project on @f$\eta@f$ axis.  If any of the pointers passed is
   * zero, do nothing.
   * 
   * @param input 
   * @param cache 
   */
  static void ProjectX(const TH2* input, TH1* cache);
  /** 
   * Add an @f$\eta@f$ bin
   * 
   * @param etaLow Low cut on @f$\eta@f$
   * @param etaMax High cut on @f$\eta@f$
   */
  void AddBin(Double_t etaLow, Double_t etaMax=kInvalidEta); 
  /** 
   * Set the maximum @f$N_{ch}@f$ to consider 
   * 
   * @param maxN Maximum
   * @param nBins (Optionally) the number of bins per particle number 
   */
  void SetMaxN(UShort_t maxN, UShort_t nBins=-1);
  void SetNDivisions(UInt_t nDiv) { fNDivisions = nDiv; }
  /** 
   * Set the range of valid interaction points 
   * 
   * @param z1 Least Z coordinate 
   * @param z2 Largest Z coordinate  
   */
  void SetIpZRange(Double_t z1, Double_t z2) { fMinIpZ = z1; fMaxIpZ = z2; }
  /** 
   * Set the trigger mask 
   * 
   * @param mask Mask 
   */
  void SetTriggerMask(UInt_t mask) { fTriggerMask = mask; }
  /** 
   * Set the trigger mask 
   * 
   * @param mask Mask 
   */
  void SetTriggerMask(const char* mask); 
  /** 
   * Whether to use the stored phi acceptance 
   * 
   * @param use If true, use stored phi acceptance 
   */
  void SetUsePhiAcc(Bool_t use) { fUsePhiAcc = use; }
  /** 
   * Print this task 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * An @f$\eta@f$ bin 
   */
  struct EtaBin : public TObject
  {
    /** 
     * I/O constructor
     */
    EtaBin();
    /** 
     * User constructor 
     * 
     * @param minEta Least @f$\eta@f$ to consider 
     * @param maxEta Largest @f$\eta@f$ to consider 
     */
    EtaBin(Double_t minEta, Double_t maxEta); 
    /** 
     * Copy constructor
     *
     * @param o object to copy fron
     */
    EtaBin(const EtaBin& o);
    /** 
     * Assignment operator
     * 
     * @param o object to assign from 
     *
     * @return Reference to this object
     */
    EtaBin& operator=(const EtaBin& o);
    /** 
     * Destructor
     */
    virtual ~EtaBin() {}
    /** 
     * Get the name of the bin
     */
    const char* GetName() const { return fName.Data(); }
    /** 
     * Is this bin symmetric around 0?
     */
    Bool_t IsSymmetric() const;
    /** 
     * Is this bin positive only?
     */
    Bool_t IsNegative() const;
    /** 
     * Is this bin positive only?
     */
    Bool_t IsPositive() const;
    /** 
     * Get parent container name 
     * 
     * 
     * @return Parent container name
     */
    const char* ParentName() const;
    /** 
     * Find the parent container.  if not found, and @a create is
     * true, then make the container.
     * 
     * @param l       Top container 
     * @param create  If true, create container if not found
     * 
     * @return Container, or null
     */
    TList* FindParent(TList* l, Bool_t create=true) const;
    /** 
     * Set-up internal structures on first event. 
     * 
     * @param list  List to add information to
     * @param hist  Template histogram 
     * @param max   Maximum number of particles 
     * @param nDiv  Number of divisions per charged particle bin
     * @param useMC Whether to set-up for MC input 
     */
    void SetupForData(TList* list, const TH2& hist, UShort_t max, 
		      UShort_t nDiv, Bool_t useMC);
    /** 
     * Process a single event 
     * 
     * @param sumForward  Projection of forward data
     * @param sumCentral  Projection of the central data
     * @param forward     The original forward data 
     * @param central     The original central data
     * @param accepted    True if event is accepted for analysis
     * @param mc          Distribution of primary particles from MC
     */
    void Process(const TH1& sumForward, const TH1& sumCentral,
		 const TH2& forward,    const TH2& central,
		 Bool_t     accepted,   const TH1* mc);
    /** 
     * Called at the end of the final processing of the job on the
     * full data set (merged data)
     * 
     * @param in    Input list
     * @param out   Output list 
     * @param maxN  Maximum number of @f$N_{ch}@f$ to consider
     */
    void Terminate(TList* in, TList* out, UShort_t maxN);
      
    TString  fName;          // Name of this bin
    Double_t fMinEta;        // Least @f$\eta@f$ to consider
    Double_t fMaxEta;        // Largest @f$\eta@f$ to consider
    Int_t    fMinBin;        // Least @f$\eta@f$ bin to consider
    Int_t    fMaxBin;        // Largest @f$\eta@f$ bin to consider
    TH1*     fSum;           // Distribution 
    TH2*     fCorr;          // Correlation between forward and central
    TH2*     fResponse;      // Response matrix (for MC)
    TH1*     fTruth;         // `true' distribution 
    TH1*     fTruthAccepted; // `true' distribution for accepted events
    TH1*     fCoverage;      // How much was covered

    ClassDef(EtaBin,1);
  };
  TList    fBins;         // List of bins 
  TList*   fSymmetric;    // Bins symmetric around 0
  TList*   fNegative;     // Bins on negative side only 
  TList*   fPositive;     // Bins on the positive side only
  TList*   fList;         // Output 
  TH1*     fTriggers;     // Histogram of triggers
  TH1*     fStatus;       // Histogram of event selection status 
  TH1*     fVertex;       // Histogram of IpZ
  TH1*     fMCVertex;     // Histogram of MC IpZ
  TH2*     fDiag;         // Diagnostics
  UInt_t   fTriggerMask;  // Trigger mask
  Double_t fMinIpZ;       // Least @f$IP_{z}@f$ to consider 
  Double_t fMaxIpZ;       // Largest @f$IP_{z}@f$ to consider 
  Bool_t   fFirstEvent;   // First-event seen or not 
  TH1*     fForwardCache; // Projection cache 
  TH1*     fCentralCache; // Projection cache 
  TH1*     fMCCache;      // Projection cache 
  UShort_t fMaxN;         // Maximum of @f$N_{ch}@f$ 
  UShort_t fNDivisions;   // Number of particle number sub-divions
  Bool_t   fUsePhiAcc;    // If true, scale by phi acceptance 

  ClassDef(AliForwardMultDists,1);
};

#endif
// Local Variables:
//  mode: C++
// End:


