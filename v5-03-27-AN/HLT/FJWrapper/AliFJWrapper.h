#ifndef AliFJWrapper_HH
#define AliFJWrapper_HH
/*
@Comments:
@
@
@
@
@
@
@
*/
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

  AliFJWrapper(const char *name, const char *title);

  virtual ~AliFJWrapper();
  
  virtual void CopySettingsFrom (const AliFJWrapper&     wrapper);
  virtual void SetupFromHeader  (const AliFastJetHeader* header);

  virtual void AddInputVector (double px, double py, double pz, double E, int index = -99999);

  virtual void AddInputVector (const fastjet::PseudoJet& vec,                int index = -99999);
  virtual void AddInputVectors(const std::vector<fastjet::PseudoJet>& vecs,  int offsetIndex = -99999);

  virtual Int_t  Run();
  virtual void   Clear(const Option_t* opt = "");

  void SetStrategy     (const fastjet::Strategy strat)             { fStrategy = strat;  }
  void SetAlgorithm    (const fastjet::JetAlgorithm algor)         { fAlgor    = algor;  }
  void SetRecombScheme (const fastjet::RecombinationScheme scheme) { fScheme   = scheme; }
  void SetAreaType     (const fastjet::AreaType atype)             { fAreaType = atype;  }

  void SetupAlgorithmfromOpt (const char *option);
  void SetupStrategyfromOpt  (const char *option);
  void SetupSchemefromOpt    (const char *option);
  void SetupAreaTypefromOpt  (const char *option);

  void SetNRepeats      (const int    nrepeat) { fNGhostRepeats  = nrepeat; }
  void SetGhostArea     (const double gharea)  { fGhostArea      = gharea;  }
  void SetMaxRap        (const double maxrap)  { fMaxRap         = maxrap;  }
  void SetR             (const double r)       { fR              = r;       }
  void SetGridScatter   (const double gridSc)  { fGridScatter    = gridSc;  }
  void SetKtScatter     (const double ktSc)    { fKtScatter      = ktSc;    }
  void SetMeanGhostKt   (const double meankt)  { fMeanGhostKt    = meankt;  }
  void SetPluginAlgor   (const int    plugin)  { fPluginAlgor    = plugin;  }
  void SetUseArea4Vector(const bool   useA4v)  { fUseArea4Vector = useA4v; }

  virtual void   GetMedianAndSigma(double& median, double& sigma, int remove = 0);
  virtual double GetMedianUsedForBgSubtraction() {return fMedianUsedForBgSubtraction;}

  const std::vector<fastjet::PseudoJet> GetInputVectors()    { return fInputVectors;}
  const std::vector<fastjet::PseudoJet> GetInclusiveJets()   { return fInclusiveJets; }

  std::vector<fastjet::PseudoJet> GetJetConstituents(const unsigned int idx);

  fastjet::ClusterSequenceArea*   GetClusterSequence() { return fClustSeq;      }

  std::vector<double>             GetSubtractedJetsPts(const double median_pt = -1, const bool sorted = false);

  double                          GetJetArea         (const unsigned int idx);
  double                          GetJetSubtractedPt (const unsigned int idx);

 protected:

  AliFJWrapper();
  AliFJWrapper(const AliFJWrapper& wrapper);
  AliFJWrapper& operator = (const AliFJWrapper& wrapper);// {return *this;}

  std::vector<fastjet::PseudoJet> fInputVectors; //!
  std::vector<fastjet::PseudoJet> fInclusiveJets; //!

  std::vector<double> fSubtractedJetsPt; //!
  
  fastjet::AreaDefinition        *fAreaDef;//!
  fastjet::VoronoiAreaSpec       *fVorAreaSpec;//!
  fastjet::GhostedAreaSpec       *fGhostedAreaSpec;//!
  fastjet::JetDefinition         *fJetDef;//!
  fastjet::JetDefinition::Plugin *fPlugin;//!
  fastjet::RangeDefinition       *fRange;//!
  fastjet::ClusterSequenceArea   *fClustSeq;//!

  // same as in the alifjheader
  fastjet::Strategy               fStrategy; //!
  fastjet::JetAlgorithm           fAlgor;    //!
  fastjet::RecombinationScheme    fScheme;   //!
  fastjet::AreaType               fAreaType; //!

  int                             fNGhostRepeats; //!
  double                          fGhostArea;	  //!
  double                          fMaxRap;	  //!
  double                          fR;             //!

  // not included in the header (within AliFastJetHeader)
  // no setters for the moment - used default values in the constructor
  double                          fGridScatter; //!
  double                          fKtScatter;	//!
  double                          fMeanGhostKt;	//!
  int                             fPluginAlgor; //!

  // extra parameters
  double                          fMedianUsedForBgSubtraction; //!
  bool                            fUseArea4Vector; //!

  virtual void   SubtractBackground(const double median_pt = -1);

  ClassDef(AliFJWrapper, 0)
};

#endif
