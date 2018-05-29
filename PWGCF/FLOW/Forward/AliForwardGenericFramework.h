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
#include "AliForwardFlowRun2Settings.h"
#include "TComplex.h"
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
   
    
  AliForwardFlowRun2Settings fSettings;
  /**
   * Do cumulants calculations for current event with 
   * centrality cent
   * 
   * @param cent Event centrality
   */
  void CumulantsAccumulate(TH2D& dNdetadphi, TList* outputList, double cent,double vertexpos,TString detType,Bool_t doRefFlow, Bool_t doDiffFlow);

  void saveEvent(TList* outputList, double cent, double vertexpos,UInt_t r);

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


  THnD* fQvector;     // Accumulated reference particles
  THnD* fpvector;    // Accumulated differential particles
  THnD* fqvector;    // Accumulated differential particles


  TComplex Q(int n, int p, int etaBin);
  TComplex p(int n, int p, int etaBin);
  TComplex q(int n, int p, int etaBin);

  TH1F fAutoRef;
  TH1F fAutoDiff;
  bool useEvent;
  bool doNUA;

  TComplex Two(int n1, int n2, int eta1, int eta2);
  TComplex TwoDiff(int n1, int n2, int refetabin, int diffetabin);
  TComplex Four(int n1, int n2, int n3, int n4,int eta1, int eta2);
  TComplex FourDiff(int n1, int n2, int n3, int n4, int refetabin, int diffetabin,int qetabin);


  ClassDef(AliForwardGenericFramework, 1); // object for eta dependent cumulant ananlysis
};

#endif