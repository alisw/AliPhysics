///
/// \file AliFemtoTrioMinvFctn.h
/// \author Jeremi Niedziela

#ifndef AliFemtoTrioMinvFctn_hh
#define AliFemtoTrioMinvFctn_hh

//#include "AliFemtoAnalysis.h"
#include "AliFemtoEvent.h"
#include "AliFemtoTrio.h"
#include "AliFemtoTrioCut.h"
#include "AliFemtoTrioFctn.h"

#include <TH1D.h>
#include <TH2D.h>

/// \class AliFemtoTrioMinvFctn
/// \brief 3-body M_inv distribution
///
///
/// This class calculated and stores distributions of 3-particle invariant mass
///

class AliFemtoTrioMinvFctn : public AliFemtoTrioFctn
{
public:
  AliFemtoTrioMinvFctn(const char* name="", int nBins=1000, double min=0.0, double max=10.0,
                       bool doMinv=true, bool doDalitz=false,bool doAngles=false);
  AliFemtoTrioMinvFctn(const AliFemtoTrioMinvFctn&);
  virtual ~AliFemtoTrioMinvFctn();

  AliFemtoTrioMinvFctn& operator=(const AliFemtoTrioMinvFctn&);
  
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
  TH2D *fDalitzPlot12_23;        // Dalitz plot for pairs 12 and 23
  TH2D *fDalitzPlot23_31;        // Dalitz plot for pairs 23 and 31
  TH2D *fDalitzPlot12_31;        // Dalitz plot for pairs 12 and 31
  
  TH1D *fAngle12;           // angle distribution between 12-particles plane and mother
  TH1D *fAngle23;           // angle distribution between 23-particles plane and mother
  TH1D *fAngle31;           // angle distribution between 31-particles plane and mother
  
  TH1D *fAngle1;           // angle distribution between 1 particle and mother
  TH1D *fAngle2;           // angle distribution between 2 particle and mother
  TH1D *fAngle3;           // angle distribution between 3 particle and mother
  
  TH1D *fCosAngle12;        // angle distribution between 12-particles plane and third particle
  TH1D *fCosAngle23;        // angle distribution between 23-particles plane and third particle
  TH1D *fCosAngle31;        // angle distribution between 31-particles plane and third particle
  
  TH1D *fCosAngle1;        // angle distribution between 1 particle and mother
  TH1D *fCosAngle2;        // angle distribution between 2 particle and mother
  TH1D *fCosAngle3;        // angle distribution between 3 particle and mother
  
  bool fDoMinv;
  bool fDoDalitz;
  bool fDoAngles;
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrioMinvFctn, 1);
  /// \endcond
#endif
};

#endif
