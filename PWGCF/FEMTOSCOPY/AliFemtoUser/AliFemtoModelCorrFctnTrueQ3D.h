///
/// \file AliFemtoModelCorrFctnTrueQ3D.h
///

#pragma once

#include "AliFemtoCorrFctn.h"

class AliFemtoPair;
class AliFemtoModelManager;
class TH3D;
class TRandom;

/// \class AliFemtoModelCorrFctnTrueQ3D
/// \brief Correlation function storing true momentum in LCMS frame.
///
/// \author Andrew Kubera, The Ohio State University <andrew.michael.kubera@cern.ch>
///
class AliFemtoModelCorrFctnTrueQ3D : public AliFemtoCorrFctn {
public:
  
  struct Parameters {
    UInt_t bin_count;
    Double_t qmin;
    Double_t qmax;
    TString title;
    AliFemtoModelManager *mc_manager;
    
    static Parameters Default() { return {56, -0.14, 0.14, "CF_TrueQ3D"}; }
    Parameters Title(const TString& title) {
      Parameters p(*this);
      p.title = title;
      return p;
    }
    Parameters Manager(AliFemtoModelManager *manager) {
      Parameters p(*this);
      p.mc_manager = manager;
      return p;
    }
    operator AliFemtoModelCorrFctnTrueQ3D*() {
      return new AliFemtoModelCorrFctnTrueQ3D(*this);
    }
  };
  
  /// Deafult parameters
  ///
  /// - Name: "CF_TrueQ3D"
  /// - Binning Paramters: (56, -0.14, 0.14)
  ///
  AliFemtoModelCorrFctnTrueQ3D();
  
  /// Custom title
  ///
  /// Use default binning parameters
  AliFemtoModelCorrFctnTrueQ3D(const char *title);
  
  /// Symmetric constructor
  ///
  /// Construct with nbins from -qmax to qmax in both directions
  AliFemtoModelCorrFctnTrueQ3D(const char *title, UInt_t nbins, Double_t qmax);
  
  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ3D(const char *title, UInt_t nbins, Double_t aQinvLo, Double_t aQinvHi);
  
  /// Construct from parameter object
  ///
  AliFemtoModelCorrFctnTrueQ3D(const Parameters &params);
  
  /// copy constructor - Structure preserved - no histogram contents copied
  AliFemtoModelCorrFctnTrueQ3D(const AliFemtoModelCorrFctnTrueQ3D& aCorrFctn);
  
  /// Assignment Operator - unused
  //   AliFemtoModelCorrFctnTrueQ3D& operator=(const AliFemtoModelCorrFctnTrueQ3D& aCorrFctn);
  
  /// Destructor - histograms destroyed, ModelManager is NOT
  virtual ~AliFemtoModelCorrFctnTrueQ3D();
  
  /// Set the MC model manager
  virtual void SetManager(AliFemtoModelManager *);
  
  virtual AliFemtoString Report();
  
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);
  
  /// Finish Data
  virtual void Finish();
  
  virtual TList* GetOutputList();
  virtual TList* AppendOutputList(TList &);
  
  virtual AliFemtoCorrFctn* Clone();
  
  Double_t GetQinvTrue(AliFemtoPair*);
  
  //Special MC analysis for K selected by PDG code -->
  void SetKaonPDG(Bool_t aSetKaonAna);
  
protected:
  AliFemtoModelManager *fManager; //!<! Link back to the manager to retrieve weights
  
  /// Numerator made with pairs using monte-carlo generated momentum
  TH3D *fNumeratorGenerated;
  
  /// Numerator made with pairs from the same event
  TH3D *fNumeratorReconstructed;
  
  /// Denominator with the monte-carlo generated momentum
  TH3D *fDenominatorGenerated;
  
  /// Denominator with reconstructed data
  TH3D *fDenominatorReconstructed;

  /// Random number generator used for randomizing order of pair momentums
  TRandom *fRng;
};

inline
void
AliFemtoModelCorrFctnTrueQ3D::SetManager(AliFemtoModelManager *manager)
{
  fManager = manager;
  std::cout << "fManager set to " << fManager << "\n";
}

inline
void
AliFemtoModelCorrFctnTrueQ3D::Finish()
{
  // no-op
}

inline
AliFemtoCorrFctn*
AliFemtoModelCorrFctnTrueQ3D::Clone()
{
  return new AliFemtoModelCorrFctnTrueQ3D(*this);
}
