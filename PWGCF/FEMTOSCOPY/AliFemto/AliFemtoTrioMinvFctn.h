///
/// \file AliFemtoTrioMinvFctn.h
/// \author Jeremi Niedziela

#ifndef AliFemtoTrioMinvFctn_hh
#define AliFemtoTrioMinvFctn_hh

//#include "AliFemtoAnalysis.h"
#include "AliFemtoEvent.h"
#include "AliFemtoTrio.h"
#include "AliFemtoTrioCut.h"

#include <TH1D.h>

/// \class AliFemtoTrioMinvFctn
/// \brief 3-body M_inv distribution
///
///
/// This class calculated and stores distributions of 3-particle invariant mass
///

class AliFemtoTrioMinvFctn
{
public:
  AliFemtoTrioMinvFctn(const char* name="", int nBins=1000, double min=0.0, double max=10.0);
  ~AliFemtoTrioMinvFctn();
  
  void AddRealTrio(AliFemtoTrio* trio);   // add real trio (same event)
  void AddMixedTrio(AliFemtoTrio* trio);  // add background trio (each particle from different event)
  
  inline void SetTrioCut(AliFemtoTrioCut* cut){fTrioCut = cut;}
  
  void Write();
  TList* GetOutputList();
  
  AliFemtoString Report(){return "";}
  void Finish(){}
  
private:
  AliFemtoTrioCut* fTrioCut;    //! this is a trio selection criteria for this distribution
  
  TH1D *fRealDistribution;  // real distribution of invariant mass
  TH1D *fMixedDistribution; // mixed distribution of invariant mass
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrioMinvFctn, 1);
  /// \endcond
#endif
};

#endif
