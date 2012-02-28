
// 
// Class that contains results from FMD eventplane calculations
//
#ifndef ALIAODFORWARDEP_H
#define ALIAODFORWARDEP_H
/**
 * @file AliAODForwardMult.h
 * @author Alexander Hansen
 * @date   Tue Feb 14 2012
 * 
 * @brief
 * 
 * 
 * @ingroup pwglf_forward
 */
#include <TObject.h>
#include <TH1D.h>
class TBrowser;
class TH1D;

class AliAODForwardEP : public TObject
{
public:
  /** 
   * Default constructor 
   * 
   * Used by ROOT I/O sub-system - do not use
   */
  AliAODForwardEP();
  /** 
   * Constructor 
   * 
   * @param isMC Whether this was from MC or not 
   */
  AliAODForwardEP(Bool_t isMC);
  /** 
   * Destructor 
   */
  virtual ~AliAODForwardEP() {} // Destructor 
  /** 
   * Initialize 
   * 
   * @param etaAxis  Pseudo-rapidity axis
   */
  void Init(const TAxis& etaAxis);
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
  Bool_t IsFolder() const { return kTRUE; } // Always true 
  /** 
   * Print content 
   * 
   * @param option Passed verbatim to TH2::Print 
   */
  void Print(Option_t* option="") const;
 /*
  * Get the name of the AOD object
  *
  * @return for now ForwardEP
  */
  const Char_t* GetName() const { return (fIsMC ? "ForwardEP" : "ForwardEP"); }
 /*
  * Set the overall FMD eventplane
  * 
  * @param ep FMD EP
  */
  void SetEventplane(Double_t ep) { fEpT = ep; }
 /* 
  * Get the overall FMD eventplane
  *
  * @return FMD EP
  */
  Double_t GetEventplane() { return fEpT; }
 /*
  * Set FMD eventplane from A side
  *
  * @param epA FMD A side EP
  */
  void SetEventplaneA(Double_t epA) { fEpA = epA; }
 /*
  * Get FMD eventplane from A side
  *
  * @return FMD A side EP
  */
  Double_t GetEventplaneA() { return fEpA; }
 /*
  * Set FMD eventplane from C side
  *
  * @param epC FMD C side EP
  */
  void SetEventplaneC(Double_t epC) { fEpC = epC; }
 /*
  * Get FMD eventplane from C side
  *
  * @return FMD C side EP
  */
  Double_t GetEventplaneC() { return fEpC; }
 /*
  * Set FMD eventplane 1 using random particles
  *
  * @param ep1 FMD EP 1 from random selection
  */
  void SetEventplane1(Double_t ep1) { fEp1 = ep1; }
 /* 
  * Get FMD eventplane 1 from random particles
  *
  * @return FMD EP 1 from random selection
  */
  Double_t GetEventplane1() { return fEp1; }
 /*
  * Set FMD eventplane 2 using random particles
  *
  * @param ep2 FMD EP 2 from random selection
  */
  void SetEventplane2(Double_t ep2) { fEp2 = ep2; }
 /* 
  * Get FMD eventplane 2 from random particles
  *
  * @return FMD EP 2 from random selection
  */
  Double_t GetEventplane2() { return fEp2; }
 /*
  * Get eta histogram from eta vs. Psi_2 calculatins
  *
  * @return Returns eta vs. Psi_2 histogram
  */
  const TH1D& GetHistogram() const { return fHist; } // Get histogram 
 /*
  * Get eta histogram from eta vs. Psi_2 calculatins
  *
  * @return Returns eta vs. Psi_2 histogram
  */
  TH1D& GetHistogram() { return fHist; } // Get histogram 

protected: 
  Bool_t   fIsMC;       // Whether this is from MC 
  Double_t fEpT;        // Total FMD event plane
  Double_t fEpA;        // FMD1+2 subevent plane
  Double_t fEpC;        // FMD3 subevent plane
  Double_t fEp1;        // Random FMD subevent plane 1
  Double_t fEp2;        // Random FMD subevent plane 2
  TH1D     fHist;       // Histogram of subevent planes in eta for this event

  ClassDef(AliAODForwardEP,1); // AOD forward multiplicity 
};

#endif
// Local Variables:
//  mode: C++
// End:
