//
//  CpMacroWithFilter.C
//
//  Created by Laurent Aphecetche
//

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TFile.h"
#include "TObjString.h"

#endif

Int_t CpMacroWithFilter(TString from, const TString& to, 
			TString filterName,
			const TString& filterMacroFullPath,
			const TString& filterRootLogonFullPath,
			const TString& alirootPath)
{
  ///
  /// Copy one file from url "from" into url "to", while filtering it at the same time.
  ///
  
  std::cout << "CpMacroWithFilter from : " << from.Data() << std::endl << " to : " << to.Data() << std::endl
	    << " with filter " << filterName.Data() << std::endl;

  if (from.IsNull() || to.IsNull()) return -1;
  
    if ( gSystem->AccessPathName(filterMacroFullPath.Data()) )
    {
      std::cerr << "Could not find requested filter macro (looking into " << filterMacroFullPath.Data() << ")" << std::endl;
      return -3;
    }

    // check that the companion macro (to load the relevant libraries 
    // and set the additional include paths, if needed) exists

    if ( gSystem->AccessPathName(filterRootLogonFullPath.Data()) )
    {
      std::cerr << "Could not find requested filter companion macro (looking into " << filterRootLogonFullPath.Data() << ")" << std::endl;
      return -4;
    }
  
  if (from.Contains("alien://")) 
  {
      TGrid::Connect("alien://");

      if (!gGrid)
	{
	  std::cerr << "Cannot get gGrid !" << std::endl;
	  return -5;
	}
  }

  if ( !filterName.IsNull() ) 
  {
    // most probably the filter will require AliRoot libs, so add the dynamic path here
    // as well as the include path and the macro path.
    //
    gSystem->AddDynamicPath(Form("%s/lib/tgt_%s",alirootPath.Data(),gSystem->GetBuildArch()));
    gSystem->SetIncludePath(Form("-I%s/include -I%s/PWG/muon",alirootPath.Data(),alirootPath.Data()));
    gROOT->SetMacroPath(Form("%s:%s/PWG/muon",gROOT->GetMacroPath(),alirootPath.Data()));
        
    // execute the companion macro
    std::cout << Form("Will load companion macro %s(\"%s\",\"%s\")",filterRootLogonFullPath.Data(),from.Data(),to.Data()) << std::endl;

    gROOT->Macro(filterRootLogonFullPath.Data());
        
    std::cout << gSystem->GetIncludePath() << std::endl;
        
    // finally delegate the work to the required filter
      
    std::cout << Form("Will compile filter %s+(\"%s\",\"%s\")",filterMacroFullPath.Data(),from.Data(),to.Data()) << std::endl;

    TString compile(filterMacroFullPath.Data());

    compile += "+";

    Int_t fail = gROOT->LoadMacro(compile.Data());
      
    if ( fail )
    {
      std::cout << Form("Failed to load/compile macro %s+",filterMacroFullPath.Data()) << std::endl;
      return -6;
    }
        
    std::cout << Form("Will execute filter %s(\"%s\",\"%s\")",filterName.Data(),from.Data(),to.Data()) << std::endl;

    return (Int_t)gROOT->ProcessLine(Form("%s(\"%s\",\"%s\")",filterName.Data(),from.Data(),to.Data()));
  }

  return 0;
}
