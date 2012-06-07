#ifndef ALIROOT_ALIFORWARDFLOWWEIGHTS
#define ALIROOT_ALIFORWARDFLOWWEIGHTS
#include <TObject.h>
class TGraph;
class TF1;
class TList;

/** 
 * Utility class to calculate flow weights 
 * 
 */
class AliForwardFlowWeights : public TObject
{
public:
  enum { 
    kEta   = 0x01,  // include eta effect 
    kPt    = 0x02,  // include pt effect
    kPID   = 0x04,  // include PID effect
    kCent  = 0x08,  // include centrality effect
    kB     = 0x10   // include impact parmaeter effect
  };
  /** 
   * Constructor 
   */
  AliForwardFlowWeights();
  /** 
   * copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardFlowWeights(const AliForwardFlowWeights& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   *
   * @return Reference to this object
   */
  AliForwardFlowWeights& operator=(const AliForwardFlowWeights& o);
  /** 
   * Destructor
   */
  virtual ~AliForwardFlowWeights();

  /** 
   * Initialize 
   * 
   * @param l List to add objects to 
   */
  virtual void Init(TList* l);
  /** 
   * @a what is a bit-wise or of 
   *
   * - kPt      Weight according to transverse momentum 
   * - kEta     Weight according to Pseudo-rapidity
   * - kPID     Weight according to particle type 
   * - kCent    Weight according to centrality 
   * - kB       Weight according to impact parameter. 
   *
   * Note, that kCent and kB are mutually exclusive 
   *
   * @a type can be one of 
   *
   * - 0  No weighting 
   *
   * - 1 Pt is weighted as mean of @f$v_2{2}@f$ and @f$v_2{4}@f$ from
   *     40-50% centrality, unity weight for p, other for @f$\pi@f$,
   *     other particles fixed, and the width of the @f$\eta@f$
   *     Gaussian is 9
   *
   * - 2 Pt is weighted by @f$v_2{2}@f$ from 40-50% centrality, fixed
   *     PID weight, and the width of the @f$\eta@f$ Gaussian is 3
   *
   * - 3 Pt is weighted by @f$v_2{4}@f$ from 30-40% centrality, unity
   *     weight for p, other for @f$\pi@f$, other particles fixed, and
   *     the width of the @f$\eta@f$ Gaussian is 15
   *
   * - 4 Pt is weighted by @f$v_2{4}@f$ from 40-50% centrality, unity
   *     weight for p, other for @f$\pi@f$, other particles fixed, and
   *     the width of the @f$\eta@f$ Gaussian is 9
   * 
   * 
   * @param eta   Pseudo-rapidity of particle (@f$\eta@f$)
   * @param pt    Transverse momentum of particle (@f$p_T@f$) in GeV
   * @param phi   Azimuthal angle of particle (@f$\phi@f$) in radians
   * @param id    Particle type of particle 
   * @param phiR  Event plane angle (@f$\phi_R@f$) in radians
   * @param bOrC  Impact paramter/Centrality of event (@f$b@f$) in fm
   * @param type  Type of flow to add 
   * @param order Order of flow weight 
   * @param what  Which effects to include
   * 
   * @return 
   */
  Double_t CalcWeight(Double_t eta, 
		      Double_t pt, 
		      Double_t phi,  
		      Int_t    id, 
		      Double_t phiR, 
		      Double_t bOrC, 
		      Int_t    type, 
		      UShort_t order,
		      UShort_t what) const;
  /** 
   * Calculate the weights 
   * 
   * @param eta  Pseudo-rapidity of particle (@f$\eta@f$)
   * @param pt   Transverse momentum of particle (@f$p_T@f$) in GeV
   * @param phi  Azimuthal angle of particle (@f$\phi@f$) in radians
   * @param id   Particle type of particle 
   * @param phiR Event plane angle (@f$\phi_R@f$) in radians
   * @param b    Impact paramter of event (@f$b@f$) in fm
   * 
   * @return 
   */
  virtual Double_t CalcWeight(Double_t eta, Double_t pt, 
			      Double_t phi, Int_t id, 
			      Double_t phiR, Double_t b) const;
  /** 
   * Construct an object from objects found in list, or null
   * 
   * @param l List to find objects in 
   * 
   * @return Newly created object, or null 
   */
  static AliForwardFlowWeights* FromList(TList* l);
protected:
  /** 
   * Calculate weight 
   * 
   * @param eta   Psuedo-rapidity 
   * @param type  Parameterization type
   * 
   * @return weight
   */
  Double_t CalcEtaWeight(Double_t eta, Int_t type) const;
  /** 
   * Calculate weight 
   * 
   * @param id    Particle (type) identifier 
   * @param type  Parameterization type
   * 
   * @return weight
   */
  Double_t CalcPidWeight(Int_t id, Int_t type) const;
  /** 
   * Calculate weight 
   * 
   * @param pt    Transverse momentum (GeV)
   * @param type  Parameterization type
   * 
   * @return weight
   */
  Double_t CalcPtWeight(Double_t pt, Int_t type) const;
  /** 
   * Calculate weight 
   * 
   * @param c Centrality
   * 
   * @return weight
   */
  Double_t CalcCentWeight(Double_t c) const;
  /** 
   * Calculate weight 
   * 
   * @param b Impact parameters (fm)
   * 
   * @return weight
   */
  Double_t CalcBWeight(Double_t b) const;

  TGraph* fV22Pt;    // Contribution from v2{2} as a function of pt
  TGraph* fV24Pt;    // Contribution from v2{4} as a function of pt
  TGraph* fV24AltPt; // Contribution from v2{4} as a function of pt
  TGraph* fV2B;      // Contribution from v2 as a function of b
  TGraph* fV2C;      // Contribution from v2 as a function of centrality
    
  ClassDef(AliForwardFlowWeights,1);
};

#endif
//
// Local Variables: 
//  mode: C++ 
// End:
//
