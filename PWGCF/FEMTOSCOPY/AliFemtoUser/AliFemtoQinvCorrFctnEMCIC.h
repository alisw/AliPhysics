/***************************************************************************
 *
 * $Id: AliFemtoQinvCorrFctnEMCIC.h  $
 *
 * Author: Nicolas Bock, Ohio State University, bock@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: Calculates of the Qinv Correlation Function, and also
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
 * The class is derived from AliFemtoQinvCorrFctn, therefore it produces
 * also the histograms in that class.
 *
 * NOTE: The EMCIC histograms are not averaged in this class, to obtain
 * the average, the user needs to divide the real pair histograms by
 * the numerator, and the mixed pairs by denominator
 *
 ***************************************************************************
 *
 **************************************************************************/

#ifndef ALIFEMTOQINVCORRFCTNEMCIC_H
#define ALIFEMTOQINVCORRFCTNEMCIC_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoQinvCorrFctn.h"


class AliFemtoQinvCorrFctnEMCIC : public AliFemtoQinvCorrFctn
{
 public:
  AliFemtoQinvCorrFctnEMCIC(const char* title, const int& nbins,
			    const float& QinvLo, const float& QinvHi);
  AliFemtoQinvCorrFctnEMCIC(const AliFemtoQinvCorrFctnEMCIC& aCorrFctn);
  virtual ~AliFemtoQinvCorrFctnEMCIC();

  AliFemtoQinvCorrFctnEMCIC& operator=(const AliFemtoQinvCorrFctnEMCIC& aCorrFctn);

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);



  virtual TList* GetOutputList();
  void Write();

 private:
  //Emcic histograms:
  /*TH1D* fESumReal;   //  <E1+E2>   from real Pairs
  TH1D* fEMultReal;  //  <E1*E2>   from real Pairs
  TH1D* fPtMultReal; //  <Pt1*Pt2> from real Pairs
  TH1D* fPzMultReal; //  <Pz1*Pz2> from real Pairs */
  TH1D* fESumMix;    //  <E1+E2>   from mixed Pairs
  TH1D* fEMultMix;   //  <E1*E2>   from mixed Pairs
  TH1D* fPtMultMix;  //  <PT1*Pt2> from mixed Pairs
  TH1D* fPzMultMix;  //  <Pz1*Pz2> from mixed Pairs




#ifdef __ROOT__
  ClassDef(AliFemtoQinvCorrFctnEMCIC, 1)
#endif
};



#endif

