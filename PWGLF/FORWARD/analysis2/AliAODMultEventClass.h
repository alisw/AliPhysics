#ifndef ALIAODMULTEVENTCLASS_H
#define ALIAODMULTEVENTCLASS_H

#include <TObject.h>
class TAxis;

class AliAODMultEventClass : public TObject
{
public:
  /** 
   * Enum of centrality types 
   */
  enum {
    kV0M = 0x1,
    kV0A = 0x2,
    kV0C = 0x4,
    kCND = 0x8,
    kEq  = 0x100
  };
  /** 
   * Constructor 
   */
  AliAODMultEventClass()
    : fMult(-1),
      fUtilV0M(-1),
      fUtilV0A(-1),
      fUtilV0C(-1),
      fUtilV0MEq(-1),
      fUtilV0AEq(-1),
      fUtilV0CEq(-1),
      fSelCND(-1),
      fSelV0M(-1),
      fSelV0A(-1),
      fSelV0C(-1),
      fSelV0MEq(-1),
      fSelV0AEq(-1),
      fSelV0CEq(-1)
  {}
  /** 
   * Get the name of this object 
   * 
   * @return The string "MultClass"
   */
  virtual const char* GetName() const { return "MultClass"; }
  /** 
   * Clear this object 
   * 
   * @param option Not used
   */
  virtual void Clear(Option_t* option="");
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;

  /** 
   * @{ 
   * @name Stuff for getting reference multiplicity etc. 
   */
  Int_t GetMultBin() const;
  /** 
   * Get the reference multiplicity in @f$|\eta|<0.8@f$ 
   * 
   * @return Reference multiplicity 
   */
  Int_t GetMult() const { return fMult; }
  /** 
   * Get the "candle" centrality estimator from AliCentrality.  This
   * is based on the number of reference tracks
   * 
   * @return 
   */
  Float_t GetCND() const { return fSelCND; }
  /** 
   * Get the defined bins 
   * 
   * @return Pointer to defined bins 
   */
  static const Int_t* GetBins();
  /** 
   * Return standard axis 
   * 
   * @return Pointer to axis (static object)
   */
  static const TAxis* GetAxis();
  /* @} */
  /** 
   * @{ 
   * @name Getters 
   */
  /** 
   * General function to get the centrality 
   * 
   * @param which String selecting the type
   * 
   * @return The centrality 
   */
  Float_t GetCentrality(const TString& which) const;
  /** 
   * General function to get the centrality 
   * 
   * @param which Bit pattern to select centrality estimate 
   * @param util  If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   * 
   * @return 
   */
  Float_t GetCentrality(UShort_t which, Bool_t util) const;
  /* @} */
  /** 
   * @{ 
   * @name Setters 
   */
  void SetCentrality(UShort_t which, Bool_t util, Float_t c);
  /** 
   * Set the reference multiplicity 
   * 
   * @param m 
   */
  void SetMult(Int_t m) { fMult = m; }
  /** @} */
protected:
  Int_t   fMult;       // The reference multiplicity
  Float_t fUtilV0M;    // V0M centrality from AliPPVsMultUtils
  Float_t fUtilV0A;    // V0A centrality from AliPPVsMultUtils
  Float_t fUtilV0C;    // V0C centrality from AliPPVsMultUtils
  Float_t fUtilV0MEq;  // V0MEq centrality from AliPPVsMultUtils
  Float_t fUtilV0AEq;  // V0MEq centrality from AliPPVsMultUtils
  Float_t fUtilV0CEq;  // V0MEq centrality from AliPPVsMultUtils
  Float_t fSelCND;
  Float_t fSelV0M;     // V0M centrality from AliCentralitySelector
  Float_t fSelV0A;     // V0A centrality from AliCentralitySelector
  Float_t fSelV0C;     // V0C centrality from AliCentralitySelector
  Float_t fSelV0MEq;   // V0MEq centrality from AliCentralitySelector
  Float_t fSelV0AEq;   // V0MEq centrality from AliCentralitySelector
  Float_t fSelV0CEq;   // V0MEq centrality from AliCentralitySelector

  ClassDef(AliAODMultEventClass,1); // Multiplicity based pp event classes
};
      
  
#endif
// Local Variables:
//   mode: C++
// End:
