#ifndef ALIAODMULTEVENTCLASS_H
#define ALIAODMULTEVENTCLASS_H

#include <TObject.h>

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
  Float_t GetCND() const { return fSelCND; }
  /** 
   * Get the defined bins 
   * 
   * @return Pointer to defined bins 
   */
  static const Int_t* GetBins();
  /* @} */
  /** 
   * @{ 
   * @name Getters 
   */
  /**
   * Get V0M centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0M Centrality
   */
  Float_t GetV0M(Bool_t util=true) const { return util ? fUtilV0M : fSelV0M; }
  /**
   * Get V0A centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0A Centrality
   */
  Float_t GetV0A(Bool_t util=true) const { return util ? fUtilV0A : fSelV0A; }
  /**
   * Get V0C centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0C Centrality
   */
  Float_t GetV0C(Bool_t util=true) const
  {
    return util ? fUtilV0C : fSelV0C;
  }
  /**
   * Get V0MEq centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0MEq Centrality
   */
  Float_t GetV0MEq(Bool_t util=true) const
  {
    return util ? fUtilV0MEq : fSelV0MEq;
  }
  /**
   * Get V0AEq centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0AEq Centrality
   */
  Float_t GetV0AEq(Bool_t util=true) const
  {
    return util ? fUtilV0AEq : fSelV0AEq;
  }
  /**
   * Get V0CEq centrality 
   *
   * @param util If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   *
   * @return V0CEq Centrality
   */
  Float_t GetV0CEq(Bool_t util=true) const
  {
    return util ? fUtilV0CEq : fSelV0CEq;
  }
  /** 
   * General function to get the centrality 
   * 
   * @param which Bit pattern to select centrality estimate 
   * @param util  If true, from AliPPVsMultUtils, otherwise
   * AliCentralitySelector
   * 
   * @return 
   */
  Float_t GetCentrality(UShort_t which, Bool_t util) const
  {
    Bool_t isEq = (which & kEq);
    if      (which & kV0M) return isEq ? GetV0MEq(util) : GetV0M(util);
    else if (which & kV0A) return isEq ? GetV0AEq(util) : GetV0A(util);
    else if (which & kV0C) return isEq ? GetV0CEq(util) : GetV0C(util);
    else if (which & kCND) return fSelCND;
    return -1;
  }
  /* @} */
  /** 
   * @{ 
   * @name Setters 
   */
  void SetCentrality(UShort_t which, Bool_t util, Float_t c)
  {
    Bool_t isEq = (which & kEq);
    if      (which & kV0M) {
      if (isEq) 
	if (util) fUtilV0MEq = c; else fSelV0MEq = c;
      else 
	if (util) fUtilV0M   = c; else fSelV0M   = c;
    }
    else if (which & kV0A) {
      if (isEq) 
	if (util) fUtilV0AEq = c; else fSelV0AEq = c;
      else 
	if (util) fUtilV0A   = c; else fSelV0A   = c;
    }
    else if (which & kV0C) {
      if (isEq) 
	if (util) fUtilV0CEq = c; else fSelV0CEq = c;
      else 
	if (util) fUtilV0C   = c; else fSelV0C   = c;
    }
    else if (which & kCND && !util)
      fSelCND = c;
  }
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
