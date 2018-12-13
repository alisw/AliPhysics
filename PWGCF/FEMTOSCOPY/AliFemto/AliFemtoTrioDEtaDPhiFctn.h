///
/// \file AliFemtoTrioDEtaDPhiFctn.h
/// \author Lukasz Graczykowski

#ifndef AliFemtoTrioDEtaDPhiFctn_hh
#define AliFemtoTrioDEtaDPhiFctn_hh

//#include "AliFemtoAnalysis.h"
#include "AliFemtoEvent.h"
#include "AliFemtoTrio.h"
#include "AliFemtoTrioCut.h"
#include "AliFemtoTrioFctn.h"

#include <TH1D.h>
#include <TH2D.h>

/// \class AliFemtoTrioDEtaDPhiFctn
/// \3-particle azimuthal correlations class
///
///
/// This class calculates and stores distributions of 3-particle azimuthal correlations
///

class AliFemtoTrioDEtaDPhiFctn : public AliFemtoTrioFctn
{
public:
  AliFemtoTrioDEtaDPhiFctn(const char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoTrioDEtaDPhiFctn(const AliFemtoTrioDEtaDPhiFctn&);
  virtual ~AliFemtoTrioDEtaDPhiFctn();

  AliFemtoTrioDEtaDPhiFctn& operator=(const AliFemtoTrioDEtaDPhiFctn&);
  
  void AddRealTrio(AliFemtoTrio* trio);   // add real trio (same event)
  void AddMixedTrio(AliFemtoTrio* trio);  // add background trio (each particle from different event)
  
  inline void SetTrioCut(AliFemtoTrioCut* cut){fTrioCut = cut;}
  
  void Write();
  TList* GetOutputList();
  
  AliFemtoString Report(){return "";}
  void Finish(){}
  
private:
  AliFemtoTrioCut* fTrioCut;    //! this is a trio selection criteria for this distribution
  
  TH2D *fRealDPhi;
  TH2D *fMixedDPhi;

  TH2D *fRealDEta;
  TH2D *fMixedDEta;

  double fphiL;
  double fphiT;
  
  int fEtaBins;
  int fPhiBins;
  
  TString fTitle;

  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoTrioDEtaDPhiFctn, 1);
  /// \endcond
#endif
};

#endif
