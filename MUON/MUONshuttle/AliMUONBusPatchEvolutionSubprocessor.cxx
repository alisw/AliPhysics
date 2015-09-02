#include "AliMUONBusPatchEvolutionSubprocessor.h"

#include "AliCDBMetaData.h"
#include "AliMUONBusPatchEvolution.h"
#include "AliMUONPreprocessor.h"
#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include <cassert>
#include <set>
#include "AliMergeableCollection.h"

//-----------------------------------------------------------------------------
/// \class AliMUONBusPatchEvolutionSubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK bus patch evolution.
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONBusPatchEvolutionSubprocessor)
/// \endcond

namespace {
	const Double_t bytes2MB = 1.0 / 1024.0 / 1024.0;
}

//_____________________________________________________________________________
AliMUONBusPatchEvolutionSubprocessor::AliMUONBusPatchEvolutionSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,"Occupancy","Upload MUON Tracker Bus Patch Evolution to OCDB"),
  fBPEVO(0x0), fProductionMode(0)
{
  /// Default ctor
}

//_____________________________________________________________________________
AliMUONBusPatchEvolutionSubprocessor::~AliMUONBusPatchEvolutionSubprocessor()
{
  /// dtor
  delete fBPEVO;
}

//_____________________________________________________________________________
Bool_t 
AliMUONBusPatchEvolutionSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the occupancy ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "BPEVO";
  
  delete fBPEVO;
  fBPEVO = 0x0;
  
  Master()->Log(Form("Reading buspatch evolution file for Run %d startTime %u endTime %u",
                     run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  Int_t n(0);
  
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    if (ReadFile(fileName.Data()))
    {
      ++n;
    }
  }
  
  delete sources;

  if (!n)
  {
    Master()->Log("Failed to read any bus patch evolution");
    delete fBPEVO;
    fBPEVO = 0;

    // return kFALSE; // for the moment, as this one is experimental, does not require it...
  }

  if ( fBPEVO )
  {
	  UInt_t hsizeBefore = fBPEVO->EstimateSize(kFALSE);

	  AliMUONBusPatchEvolution bpe(*fBPEVO);

	  bpe.ShrinkTimeAxis();

	  UInt_t hsizeAfter = fBPEVO->EstimateSize(kFALSE);

	  Master()->Log(Form("Initial collection size shrinked from %7.3f to %7.3f MB",hsizeBefore*bytes2MB,hsizeAfter*bytes2MB));

	  std::vector<int> timeResolutions;
	  bpe.GetTimeResolutions(timeResolutions);

	  if (!timeResolutions.size())
	  {
		  Master()->Log("Input mergeable collection does not seem to have time resolution histograms... Cannot work like that !");
	  }
	  else
	  {
		  TString msg;
		  for ( std::vector<int>::size_type i = 0; i < timeResolutions.size(); ++i )
		  {
			  msg += Form("%d s ",timeResolutions[i]);
		  }
		  Master()->Log(Form("Time resolutions found : %s",msg.Data()));


		  if ( !fProductionMode )
		  {
			  bpe.Augment();

			  hsizeAfter = fBPEVO->EstimateSize(kFALSE);
		  }

		  Master()->Log(Form("Final collection size is %7.3f MB",hsizeAfter*bytes2MB));
	  }
  }

  return kTRUE;
}

//_____________________________________________________________________________
UInt_t 
AliMUONBusPatchEvolutionSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the occupancy map into the CDB
  
  if (!fBPEVO)
  {
	  // as this one is still experimental, do not ever fail...
    return 0;
  }
  
  if ( fBPEVO->NumberOfObjects() )
  {
    Master()->Log("Storing buspatch evolution");
  
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("MUON TRK");
    TString comment("Computed by AliMUONBusPatchEvolutionSubprocessor $Id$");
    comment.ReplaceAll("$","");
    metaData.SetComment(comment.Data());
    
    Bool_t validToInfinity = kFALSE;
    Bool_t result = Master()->Store("Calib", "BPEVO", fBPEVO, &metaData, 0, validToInfinity);
  
    return ( result != kTRUE ); // return 0 if everything is ok.  
  }
  else
  {
    Master()->Log("No buspatch evolution to store");
    return 0;
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONBusPatchEvolutionSubprocessor::ReadFile(const char* filename)
{
  /// Read the occupancy from an ASCII file.                                  \n
  /// Return kFALSE if reading was not successfull.                           \n
  ///
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
  
  TFile* f = TFile::Open(sFilename.Data());
  
  if ( !f || !f->IsOpen() )
  {
      Master()->Log(Form("Could not open %s",sFilename.Data()));
      return kFALSE;
  }
  
  AliMergeableCollection* hc(0x0);

  hc = dynamic_cast<AliMergeableCollection*>(f->Get("bpevo"));

  if (hc)
  {
	  fBPEVO = static_cast<AliMergeableCollection*>(hc->Clone());
  }
  return ( hc != 0x0);
}


//_____________________________________________________________________________
void
AliMUONBusPatchEvolutionSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  if (fBPEVO) fBPEVO->Print(opt);
}
