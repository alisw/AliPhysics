/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   This one does 3D Bertsch-Pratt decomposition in the LCMS frame
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.5  2002/06/07 22:51:39  lisa
 * Widely used AliFemtoBPLCMS3DCorrFctn class now accumulates UNcorrected denominator and has a WriteOutHistos method
 *
 * Revision 1.4  2001/05/23 00:19:04  lisa
 * Add in Smearing classes and methods needed for momentum resolution studies and correction
 *
 * Revision 1.3  2000/10/26 19:48:50  rcwells
 * Added functionality for Coulomb correction of <qInv> in 3D correltions
 *
 * Revision 1.2  2000/09/14 18:36:53  lisa
 * Added Qinv and ExitSep pair cuts and AliFemtoBPLCMS3DCorrFctn_SIM CorrFctn
 *
 * Revision 1.1  2000/08/17 20:48:39  lisa
 * Adding correlationfunction in LCMS frame
 *
 *
 *
 **************************************************************************/

#ifndef AliFemtoBPLCMS3DCorrFctn_hh
#define AliFemtoBPLCMS3DCorrFctn_hh

#include "Base/AliFemtoCorrFctn.h"
//#include "Infrastructure/AliFemtoCoulomb.h"
#include "Base/AliFemtoPairCut.h"
//#include "Infrastructure/AliFemtoHisto.h"
#include "TH3D.h"
//#include "Infrastructure/AliFemtoSmearPair.h"

class AliFemtoBPLCMS3DCorrFctn : public AliFemtoCorrFctn {
public:
  AliFemtoBPLCMS3DCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi);
  AliFemtoBPLCMS3DCorrFctn(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn);
  virtual ~AliFemtoBPLCMS3DCorrFctn();

  AliFemtoBPLCMS3DCorrFctn& operator=(const AliFemtoBPLCMS3DCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(const AliFemtoPair*);
  virtual void AddMixedPair(const AliFemtoPair*);

  virtual void Finish();

  TH3D* Numerator();
  TH3D* Denominator();
  TH3D* Ratio();
  TH3D* QinvHisto();

  // here are get and set for the range over which the correlation function 
  // is normalized (in Qinv).  The range is set to 0.15..0.18 in the constuctor
  // by default, but the Set's below override this
  void SetNormRangeLo(float qLo);
  void SetNormRangeHi(float qHi);
  float GetNormRangeLo();
  float GetNormRangeHi();

  void WriteOutHistos();

  //  void SetCoulombCorrection(AliFemtoCoulomb* Correction);

  void SetSpecificPairCut(AliFemtoPairCut*);

  //  void SetSmearPair(AliFemtoSmearPair*);
  void SetRout(double guess);
  void SetRside(double guess);
  void SetRlong(double guess);
  void SetLambda(double guess);


  // here are a whole bunch of histos that get filled if we do resolution correction
  TH3D* fIDNumHisto;
  TH3D* fIDDenHisto;
  TH3D* fIDRatHisto;
  //
  TH3D* fSMNumHisto;
  TH3D* fSMDenHisto;
  TH3D* fSMRatHisto;
  //
  TH3D* fCorrectionHisto;
  TH3D* fCorrCFHisto;




private:
  TH3D* fNumerator;
  TH3D* fDenominator;
  //  TH3D* fUncorrectedDenominator;
  TH3D* fRatio;
  TH3D* fQinvHisto;

  // for resolution correction
  //  AliFemtoSmearPair* fSmearPair; //!
  double fLambda;
  double fRout2;
  double fRside2;
  double fRlong2;

  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis

  // upper and lower bounds of Qinv region where to do normalization
  float fQinvNormLo;
  float fQinvNormHi;

  // and here are the number of pairs in that region...
  unsigned long int fNumRealsNorm;
  unsigned long int fNumMixedNorm;

  //  AliFemtoCoulomb* fCorrection; //!


#ifdef __ROOT__
  ClassDef(AliFemtoBPLCMS3DCorrFctn, 1)
#endif
};

inline  TH3D* AliFemtoBPLCMS3DCorrFctn::Numerator(){return fNumerator;}
inline  TH3D* AliFemtoBPLCMS3DCorrFctn::Denominator(){return fDenominator;}
//inline  TH3D* AliFemtoBPLCMS3DCorrFctn::UncorrectedDenominator(){return fUncorrectedDenominator;}
inline  TH3D* AliFemtoBPLCMS3DCorrFctn::Ratio(){return fRatio;}
inline  TH3D* AliFemtoBPLCMS3DCorrFctn::QinvHisto(){return fQinvHisto;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetNormRangeLo(float qLo){fQinvNormLo = qLo;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetNormRangeHi(float qHi){fQinvNormHi = qHi;}
inline  float AliFemtoBPLCMS3DCorrFctn::GetNormRangeLo(){return fQinvNormLo;}
inline  float AliFemtoBPLCMS3DCorrFctn::GetNormRangeHi(){return fQinvNormHi;}
//inline  void AliFemtoBPLCMS3DCorrFctn::SetCoulombCorrection(AliFemtoCoulomb* Correction){fCorrection = Correction;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetSpecificPairCut(AliFemtoPairCut* pc){fPairCut=pc;}
//inline  void AliFemtoBPLCMS3DCorrFctn::SetSmearPair(AliFemtoSmearPair* sp){fSmearPair = sp;}

inline  void AliFemtoBPLCMS3DCorrFctn::SetRout(double r){fRout2 = r*r;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetRside(double r){fRside2 = r*r;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetRlong(double r){fRlong2 = r*r;}
inline  void AliFemtoBPLCMS3DCorrFctn::SetLambda(double l){fLambda = l;}

#endif

