#ifndef ALIROOT_PWG2_FORWARD_ANALYSIS_ALIAODFORWARDMULT_H
#define ALIROOT_PWG2_FORWARD_ANALYSIS_ALIAODFORWARDMULT_H
#include <TObject.h>
#include <TH2D.h>
class TBrowser;
/**
 * Class that contains the forward multiplicity data per event 
 *
 * This class contains a histogram of 
 * @f[
 *   \frac{d^2N_{ch}}{d\eta d\phi}\quad,
 * @f]
 * as well as a trigger mask for each analysed event.  
 * 
 * The eta acceptance of the event is stored in the underflow bins of
 * the histogram.  So to build the final histogram, one needs to
 * correct for this acceptance (properly weighted by the events), and
 * the vertex efficiency.  This simply boils down to defining a 2D
 * histogram and summing the event histograms in that histogram.  One
 * should of course also do proper book-keeping of the accepted event.
 *
 * @code 
 * TTree* GetAODTree()
 * { 
 *    TFile* file = TFile::Open("AliAODs.root","READ");
 *    TTree* tree = static_cast<TTree*>(file->Get("aodTree"));
 *    return tree;
 * }
 * 
 * void Analyse()
 * { 
 *   TH2D*              sum        = 0;                  // Summed hist
 *   TTree*             tree       = GetAODTree();       // AOD tree
 *   AliAODForwardMult* mult       = 0;                  // AOD object
 *   Int_t              nTriggered = 0;                  // # of triggered ev.
 *   Int_t              nWithVertex= 0;                  // # of ev. w/vertex
 *   Int_t              nAccepted  = 0;                  // # of ev. used
 *   Int_t              nAvailable = tree->GetEntries(); // How many entries
 *   Float_t            vzLow      = -10;                // Lower ip cut
 *   Float_t            vzHigh     =  10;                // Upper ip cut
 *   Int_t              mask       = AliAODForwardMult::kInel;// Trigger mask
 *   tree->SetBranchAddress("forward", &forward);        // Set the address
 * 
 *   for (int i = 0; i < nAvailable; i++) { 
 *     // Create sum histogram on first event - to match binning to input
 *     if (!sum) sum = static_cast<TH2D*>(mult->Clone("d2ndetadphi"));
 * 
 *     tree->GetEntry(i);
 * 
 *     // Other trigger/event requirements could be defined 
 *     if (!mult->IsTriggerBits(mask)) continue; 
 *     nTriggered++;
 *
 *     // Check if we have vertex 
 *     if (!mult->HasIpZ()) continue;
 *     nWithVertex++;
 * 
 *     // Select vertex range (in centimeters) 
 *     if (!mult->InRange(vzLow, vzHigh) continue; 
 *     nAccepted++;
 * 
 *     // Add contribution from this event
 *     sum->Add(&(mult->GetHistogram()));
 *   }
 * 
 *   // Get acceptance normalisation from underflow bins 
 *   TH1D* norm   = sum->Projection("norm", 0, 1, "");
 *   // Project onto eta axis - _ignoring_underflow_bins_!
 *   TH1D* dndeta = sum->Projection("dndeta", 1, -1, "e");
 *   // Normalize to the acceptance 
 *   dndeta->Divide(norm);
 *   // Scale by the vertex efficiency 
 *   dndeta->Scale(Double_t(nWithVertex)/nTriggered, "width");
 *   // And draw the result
 *   dndeta->Draw();
 * }
 * @endcode   
 *     
 * The above code will draw the final @f$ dN_{ch}/d\eta@f$ for the
 * selected event class and vertex range
 *
 * The histogram can be used as input for other kinds of analysis too, 
 * like flow, event-plane, centrality, and so on. 
 *
 * @ingroup pwg2_forward 
 */
class AliAODForwardMult : public TObject
{
public:
  /** 
   * Bits of the trigger pattern
   */
  enum { 
    /** In-elastic collision */
    kInel     = 0x001, 
    /** In-elastic collision with at least one SPD tracklet */
    kInelGt0  = 0x002, 
    /** Non-single diffractive collision */
    kNSD      = 0x004, 
    /** Empty bunch crossing */
    kEmpty    = 0x008, 
    /** A-side trigger */
    kA        = 0x010, 
    /** B(arrel) trigger */
    kB        = 0x020, 
    /** C-side trigger */
    kC        = 0x080,  
    /** Empty trigger */
    kE        = 0x100
  };
  /** 
   * Default constructor 
   * 
   * Used by ROOT I/O sub-system - do not use
   */
  AliAODForwardMult();
  /** 
   * Constructor 
   * 
   */
  AliAODForwardMult(Bool_t);
  /** 
   * Destructor 
   */
  ~AliAODForwardMult() {}
  /** 
   * Initialize 
   * 
   * @param etaAxis  Pseudo-rapidity axis
   */
  void Init(const TAxis& etaAxis);
  /** 
   * Get the @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   *
   * @return @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   */  
  const TH2D& GetHistogram() const { return fHist; }
  /** 
   * Get the @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   *
   * @return @f$ d^2N_{ch}/d\eta d\phi@f$ histogram, 
   */  
  TH2D& GetHistogram() { return fHist; }
  /** 
   * Get the trigger mask 
   * 
   * @return Trigger mask 
   */
  UInt_t GetTriggerMask() const { return fTriggers; }
  /** 
   * Set the trigger mask 
   * 
   * @param trg Trigger mask
   */
  void SetTriggerMask(UInt_t trg) { fTriggers = trg; }
  /** 
   * Set bit(s) in trigger mask 
   * 
   * @param bits bit(s) to set 
   */
  void SetTriggerBits(UInt_t bits) { fTriggers |= bits; }
  /** 
   * Check if bit(s) are set in the trigger mask 
   * 
   * @param bits Bits to test for 
   * 
   * @return 
   */
  Bool_t IsTriggerBits(UInt_t bits) const;
  /** 
   * Whether we have any trigger bits 
   */
  Bool_t HasTrigger() const { return fTriggers != 0; }
  /** 
   * Clear all data 
   * 
   * @param option  Passed on to TH2::Reset verbatim
   */
  void Clear(Option_t* option="");
  /** 
   * browse this object 
   * 
   * @param b Browser 
   */
  void Browse(TBrowser* b);
  /** 
   * This is a folder 
   * 
   * @return Always true
   */
  Bool_t IsFolder() const { return kTRUE; }
  /** 
   * Print content 
   * 
   * @param option Passed verbatim to TH2::Print 
   */
  void Print(Option_t* option="") const;
  /** 
   * Set the z coordinate of the interaction point
   * 
   * @param ipZ Interaction point z coordinate
   */
  void SetIpZ(Float_t ipZ) { fIpZ = ipZ; }
  /** 
   * Set the z coordinate of the interaction point
   * 
   * @return Interaction point z coordinate
   */
  Float_t GetIpZ() const { return fIpZ; }
  /** 
   * Check if we have a valid z coordinate of the interaction point
   *
   * @return True if we have a valid interaction point z coordinate
   */
  Bool_t HasIpZ() const;
  /** 
   * Check if the z coordinate of the interaction point is within the
   * given limits.  Note that the convention used corresponds to the
   * convention used in ROOTs TAxis.
   * 
   * @param low  Lower cut (inclusive)
   * @param high Upper cut (exclusive)
   * 
   * @return true if @f$ low \ge ipz < high@f$ 
   */
  Bool_t InRange(Float_t low, Float_t high) const;
  /** 
   * Get the name of the object 
   * 
   * @return Name of object 
   */
  const Char_t* GetName() const { return "Forward"; }
  /** 
   * Get a string correspondig to the trigger mask
   * 
   * @param mask Trigger mask 
   * 
   * @return Static string (copy before use)
   */
  static const Char_t* GetTriggerString(UInt_t mask);
protected: 
  TH2D    fHist;     // Histogram of d^2N_{ch}/(deta dphi) for this event
  UInt_t  fTriggers; // Trigger bit mask 
  Float_t fIpZ;      // Z coordinate of the interaction point

  static const Float_t fgkInvalidIpZ; // Invalid IpZ value 
  ClassDef(AliAODForwardMult,1); // AOD forward multiplicity 
};

//____________________________________________________________________
inline Bool_t
AliAODForwardMult::InRange(Float_t low, Float_t high) const 
{
  return HasIpZ() && fIpZ >= low && fIpZ < high;
}

//____________________________________________________________________
inline Bool_t 
AliAODForwardMult::IsTriggerBits(UInt_t bits) const 
{ 
  return HasTrigger() && ((fTriggers & bits) != 0); 
}
  

#endif
// Local Variables:
//  mode: C++
// End:

