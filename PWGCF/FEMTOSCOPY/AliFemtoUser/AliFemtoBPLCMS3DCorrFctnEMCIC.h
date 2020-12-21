/***************************************************************************
 *
 * $Id: AliFemtoBPLCMS3DCorrFctnEMCIC.h  $
 *
 * Author: Nicolas Bock, Ohio State University, bock@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: Calculates of the 3D Correlation Function, and also
 *              produces histograms to calculate Energy Momentum Conservation
 *              Induced Correlations  (EMCICs)
 *
 * This Class produces the following histograms as function of Qinv
 * (for both real and mixed pairs):
 *        1)   E1 + E2
 *        2)   E1 * E2
 *        3)   Pt1*Pt2
 *        4)   Pz1*Pz2
 *
 * The class is derived from AliFemtoBPLCMS3DCorrFctn, therefore it produces
 * also the histograms in that class.
 *
 * NOTE: The EMCIC histograms are not averaged in this class, to obtain
 * the average, the user needs to divide the real pair histograms by
 * the numerator, and the mixed pair histograms by the denominator
 *
 ***************************************************************************
 *
 **************************************************************************/


#ifndef ALIFEMTOBPLCMS3DCORRFCTNEMCIC_H
#define ALIFEMTOBPLCMS3DCORRFCTNEMCIC_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoBPLCMS3DCorrFctnEMCIC : public AliFemtoCorrFctn{
public:
  AliFemtoBPLCMS3DCorrFctnEMCIC(const char* title, const int& nbins, const float& QLo, const float& QHi);
  // Variable bin size constructor :
  //qBins array of low-edges for each bin. This is an array of size nbins+1
  AliFemtoBPLCMS3DCorrFctnEMCIC(const char* title, const int& nbins, const float* qBins);

  AliFemtoBPLCMS3DCorrFctnEMCIC(const AliFemtoBPLCMS3DCorrFctnEMCIC& aCorrFctn);
  virtual ~AliFemtoBPLCMS3DCorrFctnEMCIC();

  AliFemtoBPLCMS3DCorrFctnEMCIC& operator = (const AliFemtoBPLCMS3DCorrFctnEMCIC& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void Finish();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  void SetUseRPSelection(unsigned short aRPSel);


  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoBPLCMS3DCorrFctnEMCIC(*this); }

 private:

  TH3D* fNumerator;         // numerator
  TH3D* fDenominator;       // denominator
  //EMCIC histograms
  //TH3D* fEnergyTotalReal;       // E1+E2 from real pairs
  //TH3D* fEnergyMultReal;        // E1*E2
  //TH3D* fPzMultReal;            // Pz1*Pz2
  //TH3D* fPtMultReal;            // Pt1*Pt2
  TH3D* fEnergyTotalMix;       // E1+E2 from mixed pairs
  TH3D* fEnergyMultMix;        // E1*E2
  TH3D* fPzMultMix;            // Pz1*Pz2
  TH3D* fPtMultMix;            // Pt1*Pt2

 protected:
  unsigned short fUseRPSelection;  // The pair cut uses RP selection



#ifdef __ROOT__
  ClassDef(AliFemtoBPLCMS3DCorrFctnEMCIC, 1)
#endif
};

inline AliFemtoString AliFemtoBPLCMS3DCorrFctnEMCIC::Report(){AliFemtoString r="";return r;}
inline void AliFemtoBPLCMS3DCorrFctnEMCIC::Finish(){}


#endif
