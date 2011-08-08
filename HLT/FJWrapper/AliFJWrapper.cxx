#include "AliFJWrapper.h"

#include "TString.h"

#include "AliLog.h"
#include "AliFastJetHeader.h"

#include <cmath>
#include <iostream>
/*
@
@
@
@
@
@
@
 */
using namespace std;
namespace fj = fastjet;

AliFJWrapper::AliFJWrapper(const char *name, const char *title)
  : TNamed(name, title)
    
  , fInputVectors      ( )
  , fInclusiveJets     ( )

  , fSubtractedJetsPt  ( )

  , fAreaDef           (0)
  , fVorAreaSpec       (0)
  , fGhostedAreaSpec   (0)
  , fJetDef            (0)
  , fPlugin            (0)
  , fRange             (0)
  , fClustSeq          (0)

  , fStrategy          (fj::Best)
  , fAlgor             (fj::kt_algorithm)
  , fScheme            (fj::BIpt_scheme)
  , fAreaType          (fj::active_area)

  , fNGhostRepeats     (1)
  , fGhostArea         (0.01)
  , fMaxRap            (1.)
  , fR                 (0.4)

  , fGridScatter       (1.0)
  , fKtScatter         (1.0)
  , fMeanGhostKt       (1e-100)
  , fPluginAlgor       (0)

  , fMedianUsedForBgSubtraction (0)
  , fUseArea4Vector    (false)
{
  ;
}

void AliFJWrapper::CopySettingsFrom (const AliFJWrapper& wrapper)
{
  //
  // copy some settings - avoid double typing...
  // you very often want to keep most of the settings 
  // but change only the algorithm or R - do it after call to this function
  //

  fStrategy       = wrapper.fStrategy;
  fAlgor          = wrapper.fAlgor;
  fScheme         = wrapper.fScheme;
  fAreaType       = wrapper.fAreaType;
  
  fNGhostRepeats  = wrapper.fNGhostRepeats;
  fGhostArea      = wrapper.fGhostArea;
  fMaxRap         = wrapper.fMaxRap;
  fR              = wrapper.fR;
  
  fGridScatter    = wrapper.fGridScatter;
  fKtScatter      = wrapper.fKtScatter;
  fMeanGhostKt    = wrapper.fMeanGhostKt;
  fPluginAlgor    = wrapper.fPluginAlgor;
  
  fUseArea4Vector = wrapper.fUseArea4Vector;
  
}

AliFJWrapper::~AliFJWrapper()
{  
  //
  // destructor
  //
  delete fAreaDef;
  delete fVorAreaSpec;
  delete fGhostedAreaSpec;
  delete fJetDef;
  delete fPlugin;
  delete fRange;
  delete fClustSeq;  
}

void AliFJWrapper::Clear(const Option_t *opt)
{
  //
  // simply clear the input vectors
  // make sure done on every event if the instance is reused
  // reset the median to zero
  //

  fInputVectors.clear();
  fMedianUsedForBgSubtraction = 0;

  // for the moment brute force delete everything
  delete fAreaDef;
  delete fVorAreaSpec;
  delete fGhostedAreaSpec;
  delete fJetDef;
  delete fPlugin;
  delete fRange;
  delete fClustSeq;  

  TNamed::Clear(opt);
}

void AliFJWrapper::SetupFromHeader(const AliFastJetHeader *header)
{
  //
  // feed the wrapper from the parameters of the header
  //
  
  fR             = header->GetRparam();
  fStrategy      = header->GetStrategy();
  fScheme        = header->GetRecombScheme();
  fAlgor         = header->GetAlgorithm(); 
  fAreaType      = header->GetAreaType();
  fMaxRap        = header->GetGhostEtaMax(); 
  fGhostArea     = header->GetGhostArea(); 
  fNGhostRepeats = header->GetActiveAreaRepeats(); 

}

void AliFJWrapper::AddInputVector(double px, double py, double pz, double E, int index)
{
  // make the input pseudojet 
  fastjet::PseudoJet inVec(px, py, pz, E);
  
  if (index > -99999)
    {
      inVec.set_user_index(index);
    }
  else
    {
      inVec.set_user_index(fInputVectors.size());
    }

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

void AliFJWrapper::AddInputVector(const fj::PseudoJet& vec, int index)
{
  // add an input pseudojet 
  fj::PseudoJet inVec = vec;
  
  if (index > -99999)
    {
      inVec.set_user_index(index);
    }
  else
    {
      inVec.set_user_index(fInputVectors.size());
    }

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

void AliFJWrapper::AddInputVectors(const vector<fastjet::PseudoJet>& vecs, int offsetIndex)
{
  //
  // add the input from vector of pseudojets
  //

  for (unsigned int i = 0; i < vecs.size(); i++)
    {
      fj::PseudoJet inVec = vecs[i];
      
      if (offsetIndex > -99999)
	inVec.set_user_index(fInputVectors.size() + offsetIndex);
      // add to the fj container of input vectors
      
      fInputVectors.push_back(inVec);
    }
}

void AliFJWrapper::GetMedianAndSigma(double &median, double &sigma, int remove)
{
  //
  // get the median and sigma from fastjet
  // user can also do it on his own because the cluster sequence is exposed (via a getter)
  //

  if (!fClustSeq)
    {
      AliError("[e] Run the jfinder first.");
    }

  double mean_area = 0;
  try 
    {
      if(0 == remove)
	{
	  fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
	}
      else
	{
	  std::vector<fastjet::PseudoJet> input_jets = sorted_by_pt(fClustSeq->inclusive_jets());
	  input_jets.erase(input_jets.begin(), input_jets.begin() + remove);
	  fClustSeq->get_median_rho_and_sigma(input_jets, *fRange, false, median, sigma, mean_area);
	  input_jets.clear();
	}
    }
  
  catch (fj::Error)
    {
      AliError(" [w] FJ Exception caught.");
      median = -1.;
      sigma = -1;
    }
}

int AliFJWrapper::Run()
{
  //
  // run the actual jet finder
  //

  if (fAreaType == fj::voronoi_area)
    {
      // Rfact - check dependence - default is 1.
      // NOTE: hardcoded variable!
      fVorAreaSpec = new fj::VoronoiAreaSpec(1.); 
      fAreaDef     = new fj::AreaDefinition(*fVorAreaSpec);      
    }
  else
    {
      fGhostedAreaSpec = new fj::GhostedAreaSpec(fMaxRap,
						  fNGhostRepeats, 
						  fGhostArea,
						  fGridScatter,
						  fKtScatter,
						  fMeanGhostKt);

      fAreaDef = new fj::AreaDefinition(*fGhostedAreaSpec, fAreaType);
    }

  // this is acceptable by fastjet:
  fRange = new fj::RangeDefinition(fMaxRap - 0.95 * fR);

  if (fAlgor == fj::plugin_algorithm) 
    {
      if (fPluginAlgor == 0)
	{
	  // SIS CONE ALGOR
	  // NOTE: hardcoded split parameter
	  double overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
	  fPlugin = new fj::SISConePlugin(fR, 
					   overlap_threshold,
					   0, //search of stable cones - zero = until no more
					   1.0); // this should be seed effectively for proto jets
	  //__INFO("Define siscone");
	  fJetDef = new fastjet::JetDefinition(fPlugin);
	  //__INFO("SIS cone initialized");
	}
      else
	{
	  AliError("[e] Unrecognized plugin number!");
	}
    }
  else
    {
      //fJetDef = new fj::JetDefinition(fj::kt_algorithm, 0.4, fj::pt_scheme, fj::Best);
      fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);
    }
  
  try 
    {
      fClustSeq = new fj::ClusterSequenceArea(fInputVectors, *fJetDef, *fAreaDef);
    }
  catch (fj::Error)
    {
      AliError(" [w] FJ Exception caught.");
      return -1;
    }

  // inclusive jets:
  fInclusiveJets.clear();
  fInclusiveJets = fClustSeq->inclusive_jets(0.0); //sorted_by_pt(fClustSeq->inclusive_jets(0.0));

  return 0;
}

void AliFJWrapper::SubtractBackground(const double median_pt)
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
  fSubtractedJetsPt.clear();
      
  // check what was specified (default is -1)
  if (median_pt < 0)
    {
      try 
	{
	  fClustSeq->get_median_rho_and_sigma(*fRange, false, median, sigma, mean_area);
	}
      
      catch (fj::Error)
	{
	  AliError(" [w] FJ Exception caught.");
	  median = -9999.;
	  sigma = -1;
	  fMedianUsedForBgSubtraction = median;
	  return;
	}
      fMedianUsedForBgSubtraction = median;
    }
  else
    {
      // we do not know the sigma in this case
      sigma = -1;
      if (0.0 == median_pt)
	{
	  // AliWarning(Form(" [w] Median specified for bg subtraction is ZERO: %f", median_pt )) << endl;
	  fMedianUsedForBgSubtraction = 0.;
	}
      else
	{
	  fMedianUsedForBgSubtraction = median_pt;
	}
    }

  // subtract:
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) 
    {
      //fj::PseudoJet jet_sub = fClustSeq->subtracted_jet(fInclusiveJets[i], fMedianUsedForBgSubtraction_);
      //cout << fInclusiveJets[i].perp() << " (" << fMedianUsedForBgSubtraction_ << ") -> " << jet_sub.perp() << endl;
      
      if (true == fUseArea4Vector)
	{
	  // sutract the background using the area4vector
	  fj::PseudoJet area4v = fClustSeq->area_4vector(fInclusiveJets[i]);
	  fj::PseudoJet jet_sub = fInclusiveJets[i] - area4v * fMedianUsedForBgSubtraction;
	  fSubtractedJetsPt.push_back(jet_sub.perp()); // here we put only the pt of the jet - note: this can be negative
	}
      else
	{
	  // subtract the background using scalars
	  // fj::PseudoJet jet_sub = fInclusiveJets[i] - area * fMedianUsedForBgSubtraction_;

	  double area = fClustSeq->area(fInclusiveJets[i]);

	  // standard subtraction
	  double pt_sub = fInclusiveJets[i].perp() - fMedianUsedForBgSubtraction * area;
	  fSubtractedJetsPt.push_back(pt_sub); // here we put only the pt of the jet - note: this can be negative
	  //cout << "phi incl : " << fInclusiveJets[i].phi() << " subtracted : " << jet_sub.phi() << endl;
	}
    }
}

double AliFJWrapper::GetJetArea(const unsigned int idx)
{
  //
  // get the jet area
  //
  double retval = -1; // really wrong area..
  if ( idx < fInclusiveJets.size() )
    {
      retval = fClustSeq->area(fInclusiveJets[idx]);
    }
  else
    {
      cerr << "[e] ::GetJetArea wrong index: " << idx << endl;
    }
  return retval;
}

std::vector<double> AliFJWrapper::GetSubtractedJetsPts(const double median_pt, const bool sorted)
{ 
  //
  // get subtracted jets' pTs returns vector
  //
  SubtractBackground(median_pt);
  
  if (true == sorted)
    {
      std::sort(fSubtractedJetsPt.begin(), fSubtractedJetsPt.begin());
    }
  return fSubtractedJetsPt;
}

double AliFJWrapper::GetJetSubtractedPt(const unsigned int idx)
{
  //
  // get subtracted jets' pTs returns double (pT)
  //
  double retval = -99999.; // really wrong pt..
  if ( idx < fSubtractedJetsPt.size() )
    {
      retval = fSubtractedJetsPt[idx];
    }
  return retval;
}

std::vector<fastjet::PseudoJet> 
AliFJWrapper::GetJetConstituents(const unsigned int idx)
{
  //
  // get jets' constituents
  //
  std::vector<fastjet::PseudoJet> retval;
  
  if ( idx < fInclusiveJets.size() )
    {
      retval = fClustSeq->constituents(fInclusiveJets[idx]);
    }  
  else
    {
      cerr << "[e] ::GetJetConstituents wrong index: " << idx << endl;
    }
  
  return retval;
}

void AliFJWrapper::SetupAlgorithmfromOpt(const char *option)
{
  //
  // setup algorithm from char
  //
  string opt(option);
  
  if (!opt.compare("kt"))                fAlgor    = fj::kt_algorithm;
  if (!opt.compare("antikt"))            fAlgor    = fj::antikt_algorithm;
  if (!opt.compare("cambridge"))         fAlgor    = fj::cambridge_algorithm;
  if (!opt.compare("genkt"))             fAlgor    = fj::genkt_algorithm;
  if (!opt.compare("cambridge_passive")) fAlgor    = fj::cambridge_for_passive_algorithm;
  if (!opt.compare("genkt_passive"))     fAlgor    = fj::genkt_for_passive_algorithm;
  if (!opt.compare("ee_kt"))             fAlgor    = fj::ee_kt_algorithm;
  if (!opt.compare("ee_genkt"))          fAlgor    = fj::ee_genkt_algorithm;
  if (!opt.compare("plugin"))            fAlgor    = fj::plugin_algorithm;


}
void AliFJWrapper::SetupStrategyfromOpt(const char *option)
{
  //
  // setup strategy from char
  //
  string opt(option);

  if (!opt.compare("Best"))            fStrategy = fj::Best;
  if (!opt.compare("N2MinHeapTiled"))  fStrategy = fj::N2MinHeapTiled;
  if (!opt.compare("N2Tiled"))         fStrategy = fj::N2Tiled;
  if (!opt.compare("N2PoorTiled"))     fStrategy = fj::N2PoorTiled;
  if (!opt.compare("N2Plain"))         fStrategy = fj::N2Plain;
  if (!opt.compare("N3Dumb"))          fStrategy = fj::N3Dumb;
  if (!opt.compare("NlnN"))            fStrategy = fj::NlnN;
  if (!opt.compare("NlnN3pi"))         fStrategy = fj::NlnN3pi;
  if (!opt.compare("NlnN4pi"))         fStrategy = fj::NlnN4pi;
  if (!opt.compare("NlnNCam4pi"))      fStrategy = fj::NlnNCam4pi;
  if (!opt.compare("NlnNCam2pi2R"))    fStrategy = fj::NlnNCam2pi2R;
  if (!opt.compare("NlnNCam"))         fStrategy = fj::NlnNCam;
  if (!opt.compare("plugin"))          fStrategy = fj::plugin_strategy;

}
void AliFJWrapper::SetupSchemefromOpt(const char *option)
{
  //
  // setup scheme from char
  //
  string opt(option);

  if (!opt.compare("BIpt"))   fScheme   = fj::BIpt_scheme;
  if (!opt.compare("BIpt2"))  fScheme   = fj::BIpt2_scheme;
  if (!opt.compare("E"))      fScheme   = fj::E_scheme;
  if (!opt.compare("pt"))     fScheme   = fj::pt_scheme;
  if (!opt.compare("pt2"))    fScheme   = fj::pt2_scheme;
  if (!opt.compare("Et"))     fScheme   = fj::Et_scheme;
  if (!opt.compare("Et2"))    fScheme   = fj::Et2_scheme;

}
void AliFJWrapper::SetupAreaTypefromOpt(const char *option)
{
  //
  // setup area type from char
  //
  string opt(option);

  if (!opt.compare("active"))                      fAreaType = fj::active_area;
  if (!opt.compare("invalid"))                     fAreaType = fj::invalid_area;
  if (!opt.compare("active_area_explicit_ghosts")) fAreaType = fj::active_area_explicit_ghosts;
  if (!opt.compare("one_ghost_passive"))           fAreaType = fj::one_ghost_passive_area;
  if (!opt.compare("passive"))                     fAreaType = fj::passive_area;
  if (!opt.compare("voronoi"))                     fAreaType = fj::voronoi_area;

}
