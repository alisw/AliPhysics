//
// Helper class to calculate Q cumulant in forward & central regions
//
#ifndef AliForwardGenericFramework_H
#define AliForwardGenericFramework_H
/**
 * @file AliForwardGenericFramework.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include <TObject.h>
#include <TH2D.h>
#include <TH3D.h>
#include "TString.h"
#include "TNtuple.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"
#include <THn.h>
#include "TString.h"
#include "AliForwardSettings.h"
#include "TComplex.h"
#include "AliForwardFlowUtil.h"
/**
 * Class to handle cumulant calculations.
 */
class AliForwardGenericFramework
{
public:
  /*
  * Constructor
  */
  AliForwardGenericFramework();

  /**
   * Destructor
   */
  virtual ~AliForwardGenericFramework(){}


  AliForwardSettings fSettings;

    // Utility class for filling histograms
  //AliForwardFlowUtil fUtil;
  /**
   * Do cumulants calculations for current event with
   * centrality cent
   *
   * @param cent Event centrality
   */
  void CumulantsAccumulate(TH2D*& dNdetadphi, double cent,double vertexpos,Bool_t useFMD,Bool_t doRefFlow, Bool_t doDiffFlow);

  void saveEvent(TList* outputList, double cent, double vertexpos,UInt_t r, Int_t ptn);

  /**
   * Constants
   */
  enum {
    ktpcOnly = 128,        // TPC only tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };

  /**
   * Reset histograms
   */
  void reset();

  THnD* fQvector;//!     // Accumulated reference particles
  THnD* fpvector;//!    // Accumulated differential particles
  THnD* fqvector;//!    // Accumulated differential particles

  TComplex Q(Int_t n, Int_t p, Int_t etaBin);
  TComplex p(Int_t n, Int_t p, Int_t etaBin);
  TComplex q(Int_t n, Int_t p, Int_t etaBin);

  // TH1F fAutoRef;
  // TH1F fAutoDiff;

  TComplex Two(Int_t n1, Int_t n2, Int_t eta1, Int_t eta2);
  TComplex TwoDiff(Int_t n1, Int_t n2, Int_t refetabin, Int_t diffetabin);
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4,Int_t eta1, Int_t eta2);
  TComplex FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t refetabinA, Int_t refetabinB, Int_t diffetabin,Int_t qetabin);


  ClassDef(AliForwardGenericFramework, 1); // object for eta dependent cumulant ananlysis
};

#endif
