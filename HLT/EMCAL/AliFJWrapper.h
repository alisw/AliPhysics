#ifndef AliFJWrapper_HH
#define AliFJWrapper_HH

#include <vector>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/SISConePlugin.hh>

#include "TNamed.h"

class AliFastJetHeader;

class AliFJWrapper : public TNamed
{
 public:

  AliFJWrapper(const char *_name, const char *_title);

  virtual ~AliFJWrapper();
  
  virtual void CopySettingsFrom (const AliFJWrapper&     _wrapper);
  virtual void SetupFromHeader  (const AliFastJetHeader* _header);

  virtual void AddInputVector (double _px, double _py, double _pz, double _E, int _index = -99999);

  virtual void AddInputVector (const fastjet::PseudoJet& _vec,                int _index = -99999);
  virtual void AddInputVectors(const std::vector<fastjet::PseudoJet>& _vecs,  int _offsetIndex = -99999);

  virtual Int_t  Run();
  virtual void   Clear(const Option_t* _opt = "");

  void SetStrategy     (const fastjet::Strategy _strat)             { fStrategy_ = _strat;  }
  void SetAlgorithm    (const fastjet::JetAlgorithm _algor)         { fAlgor_    = _algor;  }
  void SetRecombScheme (const fastjet::RecombinationScheme _scheme) { fScheme_   = _scheme; }
  void SetAreaType     (const fastjet::AreaType _atype)             { fAreaType_ = _atype;  }

  void SetNRepeats      (const int    _nrepeat) { fNGhostRepeats_ = _nrepeat; }
  void SetGhostArea     (const double _gharea)  { fGhostArea_     = _gharea;  }
  void SetMaxRap        (const double _maxrap)  { fMaxRap_        = _maxrap;  }
  void SetR             (const double _r)       { fR_             = _r;       }
  void SetGridScatter   (const double _gridSc)  { fGridScatter_   = _gridSc;  }
  void SetKtScatter     (const double _ktSc)    { fKtScatter_     = _ktSc;    }
  void SetMeanGhostKt   (const double _meankt)  { fMeanGhostKt_   = _meankt;  }
  void SetPluginAlgor   (const int    _plugin)  { fPluginAlgor_   = _plugin;  }
  void SetUseArea4Vector(const bool   _useA4v)  { fUseArea4Vector_ = _useA4v; }

  virtual void   GetMedianAndSigma(double& _median, double& _sigma, int _remove = 0);
  virtual double GetMedianUsedForBgSubtraction() {return fMedianUsedForBgSubtraction_;}

  std::vector<fastjet::PseudoJet> GetInputVectors()    { return fInputVectors_;}
  std::vector<fastjet::PseudoJet> GetInclusiveJets()   { return fInclusiveJets_; }

  std::vector<fastjet::PseudoJet> GetJetConstituents(const unsigned int _idx);

  fastjet::ClusterSequenceArea*   GetClusterSequence() { return fClustSeq_;      }

  std::vector<double>             GetSubtractedJetsPts(const double _median_pt = -1, const bool _sorted = false);

  double                          GetJetArea         (const unsigned int _idx);
  double                          GetJetSubtractedPt (const unsigned int _idx);

 protected:

  AliFJWrapper();
  AliFJWrapper(const AliFJWrapper& _wrapper);
  AliFJWrapper& operator = (const AliFJWrapper& _wrapper);// {return *this;}

  std::vector<fastjet::PseudoJet> fInputVectors_;
  std::vector<fastjet::PseudoJet> fInclusiveJets_;

  std::vector<double> fSubtractedJetsPt_;
  
  fastjet::AreaDefinition        *fAreaDef_;//!
  fastjet::VoronoiAreaSpec       *fVorAreaSpec_;//!
  fastjet::GhostedAreaSpec       *fGhostedAreaSpec_;//!
  fastjet::JetDefinition         *fJetDef_;//!
  fastjet::JetDefinition::Plugin *fPlugin_;//!
  fastjet::RangeDefinition       *fRange_;//!
  fastjet::ClusterSequenceArea   *fClustSeq_;//!

  // same as in the alifjheader
  fastjet::Strategy               fStrategy_;
  fastjet::JetAlgorithm           fAlgor_;
  fastjet::RecombinationScheme    fScheme_;
  fastjet::AreaType               fAreaType_;

  int                             fNGhostRepeats_;
  double                          fGhostArea_;
  double                          fMaxRap_;
  double                          fR_;

  // not included in the header (within AliFastJetHeader)
  // no setters for the moment - used default values in the constructor
  double                          fGridScatter_;
  double                          fKtScatter_;
  double                          fMeanGhostKt_;
  int                             fPluginAlgor_;

  // extra parameters
  double                          fMedianUsedForBgSubtraction_;
  bool                            fUseArea4Vector_;

  virtual void   SubtractBackground(const double _median_pt = -1);

  ClassDef(AliFJWrapper, 0)
};

#endif
