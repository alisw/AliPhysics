#include "AliFJWrapper.h"

#include "TString.h"

#include "AliLog.h"
#include "AliFastJetHeader.h"

#include <cmath>
#include <iostream>

using namespace std;
namespace fj = fastjet;

AliFJWrapper::AliFJWrapper(const char *_name, const char *_title)
  : TNamed(_name, _title)
    
  , fInputVectors_      ( )
  , fInclusiveJets_     ( )

  , fSubtractedJetsPt_  ( )

  , fAreaDef_           (0)
  , fVorAreaSpec_       (0)
  , fGhostedAreaSpec_   (0)
  , fJetDef_            (0)
  , fPlugin_            (0)
  , fRange_             (0)
  , fClustSeq_          (0)

  , fStrategy_          (fj::Best)
  , fAlgor_             (fj::kt_algorithm)
  , fScheme_            (fj::BIpt_scheme)
  , fAreaType_          (fj::active_area)

  , fNGhostRepeats_     (1)
  , fGhostArea_         (0.01)
  , fMaxRap_            (1.)
  , fR_                 (0.4)

  , fGridScatter_       (1.0)
  , fKtScatter_         (1.0)
  , fMeanGhostKt_       (1e-100)
  , fPluginAlgor_       (0)

  , fMedianUsedForBgSubtraction_ (0)
  , fUseArea4Vector_    (false)
{
  ;
}

void AliFJWrapper::CopySettingsFrom (const AliFJWrapper& _wrapper)
{
  //
  // copy some settings - avoid double typing...
  // you very often want to keep most of the settings 
  // but change only the algorithm or R - do it after call to this function
  //

  fStrategy_       = _wrapper.fStrategy_;
  fAlgor_          = _wrapper.fAlgor_;
  fScheme_         = _wrapper.fScheme_;
  fAreaType_       = _wrapper.fAreaType_;
  
  fNGhostRepeats_  = _wrapper.fNGhostRepeats_;
  fGhostArea_      = _wrapper.fGhostArea_;
  fMaxRap_         = _wrapper.fMaxRap_;
  fR_              = _wrapper.fR_;
  
  fGridScatter_    = _wrapper.fGridScatter_;
  fKtScatter_      = _wrapper.fKtScatter_;
  fMeanGhostKt_    = _wrapper.fMeanGhostKt_;
  fPluginAlgor_    = _wrapper.fPluginAlgor_;
  
  fUseArea4Vector_ = _wrapper.fUseArea4Vector_;
  
}

AliFJWrapper::~AliFJWrapper()
{  
  delete fAreaDef_;
  delete fVorAreaSpec_;
  delete fGhostedAreaSpec_;
  delete fJetDef_;
  delete fPlugin_;
  delete fRange_;
  delete fClustSeq_;  
}

void AliFJWrapper::Clear(const Option_t *_opt)
{
  //
  // simply clear the input vectors
  // make sure done on every event if the instance is reused
  // reset the median to zero
  //

  fInputVectors_.clear();
  fMedianUsedForBgSubtraction_ = 0;

  // for the moment brute force delete everything
  delete fAreaDef_;
  delete fVorAreaSpec_;
  delete fGhostedAreaSpec_;
  delete fJetDef_;
  delete fPlugin_;
  delete fRange_;
  delete fClustSeq_;  

  TNamed::Clear(_opt);
}

void AliFJWrapper::SetupFromHeader(const AliFastJetHeader *_header)
{
  //
  // feed the wrapper from the parameters of the header
  //
  
  fR_             = _header->GetRparam();
  fStrategy_      = _header->GetStrategy();
  fScheme_        = _header->GetRecombScheme();
  fAlgor_         = _header->GetAlgorithm(); 
  fAreaType_      = _header->GetAreaType();
  fMaxRap_        = _header->GetGhostEtaMax(); 
  fGhostArea_     = _header->GetGhostArea(); 
  fNGhostRepeats_ = _header->GetActiveAreaRepeats(); 

}

void AliFJWrapper::AddInputVector(double _px, double _py, double _pz, double _E, int _index)
{
  // make the input pseudojet 
  fastjet::PseudoJet inVec(_px, _py, _pz, _E);
  
  if (_index > -99999)
    {
      inVec.set_user_index(_index);
    }
  else
    {
      inVec.set_user_index(fInputVectors_.size());
    }

  // add to the fj container of input vectors
  fInputVectors_.push_back(inVec);
}

void AliFJWrapper::AddInputVector(const fj::PseudoJet& _vec, int _index)
{
  // add an input pseudojet 
  fj::PseudoJet inVec = _vec;
  
  if (_index > -99999)
    {
      inVec.set_user_index(_index);
    }
  else
    {
      inVec.set_user_index(fInputVectors_.size());
    }

  // add to the fj container of input vectors
  fInputVectors_.push_back(inVec);
}

void AliFJWrapper::AddInputVectors(const vector<fj::PseudoJet>& _vecs, int _offsetIndex)
{
  //
  // add the input from vector of pseudojets
  //

  for (unsigned int i = 0; i < _vecs.size(); i++)
    {
      fj::PseudoJet inVec = _vecs[i];
      
      if (_offsetIndex > -99999)
	{
	  inVec.set_user_index(fInputVectors_.size() + _offsetIndex);
	}
      else
	{
	  inVec.set_user_index(fInputVectors_.size());
	}
      
      // add to the fj container of input vectors
      fInputVectors_.push_back(inVec);
    }
}

void AliFJWrapper::GetMedianAndSigma(double &_median, double &_sigma, int _remove)
{
  //
  // get the median and sigma from fastjet
  // user can also do it on his own because the cluster sequence is exposed (via a getter)
  //

  if (!fClustSeq_)
    {
      AliError("[e] Run the jfinder first.");
    }

  double mean_area = 0;
  try 
    {
      if(0 == _remove)
	{
	  fClustSeq_->get_median_rho_and_sigma(*fRange_, fUseArea4Vector_, _median, _sigma, mean_area);
	}
      else
	{
	  std::vector<fastjet::PseudoJet> input_jets = sorted_by_pt(fClustSeq_->inclusive_jets());
	  input_jets.erase(input_jets.begin(), input_jets.begin() + _remove);
	  fClustSeq_->get_median_rho_and_sigma(input_jets, *fRange_, false, _median, _sigma, mean_area);
	  input_jets.clear();
	}
    }
  
  catch (fj::Error)
    {
      AliError(" [w] FJ Exception caught.");
      _median = -1.;
      _sigma = -1;
    }
}

int AliFJWrapper::Run()
{
  //
  // run the actual jet finder
  //

  if (fAreaType_ == fj::voronoi_area)
    {
      // Rfact - check dependence - default is 1.
      // NOTE: hardcoded variable!
      fVorAreaSpec_ = new fj::VoronoiAreaSpec(1.); 
      fAreaDef_     = new fj::AreaDefinition(*fVorAreaSpec_);      
    }
  else
    {
      fGhostedAreaSpec_ = new fj::GhostedAreaSpec(fMaxRap_,
						  fNGhostRepeats_, 
						  fGhostArea_,
						  fGridScatter_,
						  fKtScatter_,
						  fMeanGhostKt_);

      fAreaDef_ = new fj::AreaDefinition(*fGhostedAreaSpec_, fAreaType_);
    }

  // this is acceptable by fastjet:
  fRange_ = new fj::RangeDefinition(fMaxRap_ - 0.95 * fR_);

  if (fAlgor_ == fj::plugin_algorithm) 
    {
      if (fPluginAlgor_ == 0)
	{
	  // SIS CONE ALGOR
	  // NOTE: hardcoded split parameter
	  double overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
	  fPlugin_ = new fj::SISConePlugin(fR_, 
					   overlap_threshold,
					   0, //search of stable cones - zero = until no more
					   1.0); // this should be seed effectively for proto jets
	  //__INFO("Define siscone");
	  fJetDef_ = new fastjet::JetDefinition(fPlugin_);
	  //__INFO("SIS cone initialized");
	}
      else
	{
	  AliError("[e] Unrecognized plugin number!");
	}
    }
  else
    {
      //fJetDef_ = new fj::JetDefinition(fj::kt_algorithm, 0.4, fj::pt_scheme, fj::Best);
      fJetDef_ = new fj::JetDefinition(fAlgor_, fR_, fScheme_, fStrategy_);
    }
  
  try 
    {
      fClustSeq_ = new fj::ClusterSequenceArea(fInputVectors_, *fJetDef_, *fAreaDef_);
    }
  catch (fj::Error)
    {
      AliError(" [w] FJ Exception caught.");
      return -1;
    }

  // inclusive jets:
  fInclusiveJets_.clear();
  fInclusiveJets_ = fClustSeq_->inclusive_jets(0.0); //sorted_by_pt(fClustSeq_->inclusive_jets(0.0));

  return 0;
}

void AliFJWrapper::SubtractBackground(const double _median_pt)
{
  //
  // subtract the background (specify the value - see below the meaning)
  // negative argument means the bg will be determined with the current algorithm
  // this is the default behavior
  // zero is allowed
  // Note: user may set the switch for area4vector based subtraction
  //

  double median = 0;
  double sigma = 0;
  double mean_area = 0;

  // clear the subtracted jet pt's vector<double>
  fSubtractedJetsPt_.clear();
      
  // check what was specified (default is -1)
  if (_median_pt < 0)
    {
      try 
	{
	  fClustSeq_->get_median_rho_and_sigma(*fRange_, false, median, sigma, mean_area);
	}
      
      catch (fj::Error)
	{
	  AliError(" [w] FJ Exception caught.");
	  median = -9999.;
	  sigma = -1;
	  fMedianUsedForBgSubtraction_ = median;
	  return;
	}
      fMedianUsedForBgSubtraction_ = median;
    }
  else
    {
      // we do not know the sigma in this case
      sigma = -1;
      if (0.0 == _median_pt)
	{
	  // AliWarning(Form(" [w] Median specified for bg subtraction is ZERO: %f", median_pt )) << endl;
	  fMedianUsedForBgSubtraction_ = 0.;
	}
      else
	{
	  fMedianUsedForBgSubtraction_ = _median_pt;
	}
    }

  // subtract:
  for (unsigned i = 0; i < fInclusiveJets_.size(); i++) 
    {
      //fj::PseudoJet jet_sub = fClustSeq_->subtracted_jet(fInclusiveJets_[i], fMedianUsedForBgSubtraction_);
      //cout << fInclusiveJets_[i].perp() << " (" << fMedianUsedForBgSubtraction_ << ") -> " << jet_sub.perp() << endl;
      
      if (true == fUseArea4Vector_)
	{
	  // sutract the background using the area4vector
	  fj::PseudoJet area4v = fClustSeq_->area_4vector(fInclusiveJets_[i]);
	  fj::PseudoJet jet_sub = fInclusiveJets_[i] - area4v * fMedianUsedForBgSubtraction_;
	  fSubtractedJetsPt_.push_back(jet_sub.perp()); // here we put only the pt of the jet - note: this can be negative
	}
      else
	{
	  // subtract the background using scalars
	  // fj::PseudoJet jet_sub = fInclusiveJets_[i] - area * fMedianUsedForBgSubtraction_;

	  double area = fClustSeq_->area(fInclusiveJets_[i]);

	  // standard subtraction
	  double pt_sub = fInclusiveJets_[i].perp() - fMedianUsedForBgSubtraction_ * area;
	  fSubtractedJetsPt_.push_back(pt_sub); // here we put only the pt of the jet - note: this can be negative
	  //cout << "phi incl : " << fInclusiveJets_[i].phi() << " subtracted : " << jet_sub.phi() << endl;
	}
    }
}

double AliFJWrapper::GetJetArea(const unsigned int _idx)
{
  double retval = -1; // really wrong area..
  if ( _idx < fInclusiveJets_.size() )
    {
      retval = fClustSeq_->area(fInclusiveJets_[_idx]);
    }
  else
    {
      cerr << "[e] ::GetJetArea wrong index: " << _idx << endl;
    }
  return retval;
}

std::vector<double> AliFJWrapper::GetSubtractedJetsPts(const double _median_pt, const bool _sorted)
{ 
  SubtractBackground(_median_pt);
  
  if (true == _sorted)
    {
      std::sort(fSubtractedJetsPt_.begin(), fSubtractedJetsPt_.begin());
    }
  return fSubtractedJetsPt_;
}

double AliFJWrapper::GetJetSubtractedPt(const unsigned int _idx)
{
  double retval = -99999.; // really wrong pt..
  if ( _idx < fSubtractedJetsPt_.size() )
    {
      retval = fSubtractedJetsPt_[_idx];
    }
  return retval;
}

std::vector<fastjet::PseudoJet> 
AliFJWrapper::GetJetConstituents(const unsigned int _idx)
{
  std::vector<fastjet::PseudoJet> retval;
  
  if ( _idx < fInclusiveJets_.size() )
    {
      retval = fClustSeq_->constituents(fInclusiveJets_[_idx]);
    }  
  else
    {
      cerr << "[e] ::GetJetConstituents wrong index: " << _idx << endl;
    }
  
  return retval;
}

void AliFJWrapper::SetupAlgorithmfromOpt(const char *option)
{
  string opt(option);
  
  if (!opt.compare("kt"))                fAlgor_    = fj::kt_algorithm;
  if (!opt.compare("antikt"))            fAlgor_    = fj::antikt_algorithm;
  if (!opt.compare("cambridge"))         fAlgor_    = fj::cambridge_algorithm;
  if (!opt.compare("genkt"))             fAlgor_    = fj::genkt_algorithm;
  if (!opt.compare("cambridge_passive")) fAlgor_    = fj::cambridge_for_passive_algorithm;
  if (!opt.compare("genkt_passive"))     fAlgor_    = fj::genkt_for_passive_algorithm;
  if (!opt.compare("ee_kt"))             fAlgor_    = fj::ee_kt_algorithm;
  if (!opt.compare("ee_genkt"))          fAlgor_    = fj::ee_genkt_algorithm;
  if (!opt.compare("plugin"))            fAlgor_    = fj::plugin_algorithm;


}
void AliFJWrapper::SetupStrategyfromOpt(const char *option)
{
  string opt(option);

  if (!opt.compare("Best"))            fStrategy_ = fj::Best;
  if (!opt.compare("N2MinHeapTiled"))  fStrategy_ = fj::N2MinHeapTiled;
  if (!opt.compare("N2Tiled"))         fStrategy_ = fj::N2Tiled;
  if (!opt.compare("N2PoorTiled"))     fStrategy_ = fj::N2PoorTiled;
  if (!opt.compare("N2Plain"))         fStrategy_ = fj::N2Plain;
  if (!opt.compare("N3Dumb"))          fStrategy_ = fj::N3Dumb;
  if (!opt.compare("NlnN"))            fStrategy_ = fj::NlnN;
  if (!opt.compare("NlnN3pi"))         fStrategy_ = fj::NlnN3pi;
  if (!opt.compare("NlnN4pi"))         fStrategy_ = fj::NlnN4pi;
  if (!opt.compare("NlnNCam4pi"))      fStrategy_ = fj::NlnNCam4pi;
  if (!opt.compare("NlnNCam2pi2R"))    fStrategy_ = fj::NlnNCam2pi2R;
  if (!opt.compare("NlnNCam"))         fStrategy_ = fj::NlnNCam;
  if (!opt.compare("plugin"))          fStrategy_ = fj::plugin_strategy;

}
void AliFJWrapper::SetupSchemefromOpt(const char *option)
{
  string opt(option);

  if (!opt.compare("BIpt"))   fScheme_   = fj::BIpt_scheme;
  if (!opt.compare("BIpt2"))  fScheme_   = fj::BIpt2_scheme;
  if (!opt.compare("E"))      fScheme_   = fj::E_scheme;
  if (!opt.compare("pt"))     fScheme_   = fj::pt_scheme;
  if (!opt.compare("pt2"))    fScheme_   = fj::pt2_scheme;
  if (!opt.compare("Et"))     fScheme_   = fj::Et_scheme;
  if (!opt.compare("Et2"))    fScheme_   = fj::Et2_scheme;

}
void AliFJWrapper::SetupAreaTypefromOpt(const char *option)
{
  string opt(option);

  if (!opt.compare("active"))                      fAreaType_ = fj::active_area;
  if (!opt.compare("invalid"))                     fAreaType_ = fj::invalid_area;
  if (!opt.compare("active_area_explicit_ghosts")) fAreaType_ = fj::active_area_explicit_ghosts;
  if (!opt.compare("one_ghost_passive"))           fAreaType_ = fj::one_ghost_passive_area;
  if (!opt.compare("passive"))                     fAreaType_ = fj::passive_area;
  if (!opt.compare("voronoi"))                     fAreaType_ = fj::voronoi_area;

}
