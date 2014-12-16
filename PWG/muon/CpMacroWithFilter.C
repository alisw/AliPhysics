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

Int_t CpMacroWithFilter(TString from, TString to)
{
  /// Copy one file from url "from" into url "to".
  ///
  /// If from contains one of the filter keywords
  /// (in the form of e.g. AliAOD.FILTER_ZZZZ_WITH_ALIROOT_YYYY.root)
  /// then we call the external filter macro
  ///
  /// otherwise we do a simple TFile::Cp(from,to)
  ///

  const char* swBaseDir = "/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliRoot";
  
  if (from.IsNull() || to.IsNull()) return -1;
  
  std::cout << Form("Entering CpMacroWithFilter(\"%s\",\"%s\")",from.Data(),to.Data()) << std::endl;
  
  TString filterName;
  TString alirootVersion;
  TString alirootPath(swBaseDir);
  TString filterMacroFullPath;
  TString filterRootLogonFullPath;

  Int_t ix = from.Index("FILTER_");
  
  if ( ix > 0)
  {
    TString part = from(ix,from.Length()-ix+1);
    TObjArray* tmp = part.Tokenize(".");
    TObjString* ostr = static_cast<TObjString*>(tmp->First());
    if (!ostr)
    {
      std::cerr << "Could not get filter ??? Filename does not look right !!!" << std::endl;
      return -2;
    }
    filterName = ostr->String();
    delete tmp;
    ix = filterName.Index("_WITH_");
    if ( ix > 0 )
    {
      alirootVersion = filterName(ix+strlen("_WITH_ALIROOT_"),part.Length()-ix-strlen("_WITH_ALIROOT_"));
      filterName = filterName(0,ix);
    }

    alirootPath += "/";
    alirootPath += alirootVersion;

    // check the aliroot version required actually exists on cvmfs...
    if ( gSystem->AccessPathName(alirootPath.Data()) ) 
    {
	std::cerr << "Requested AliRoot version not found (looking into " << alirootPath.Data() << ")" << std::endl;
	return -2;
    }

    // check the filter required actually exists
    filterMacroFullPath.Form("%s/PWG/muon/%s.C",alirootPath.Data(),filterName.Data());
   
    if ( gSystem->AccessPathName(filterMacroFullPath.Data()) )
    {
      std::cerr << "Could not find requested filter macro (looking into " << filterMacroFullPath.Data() << ")" << std::endl;
      return -3;
    }

    // check that the companion macro (to load the relevant libraries 
    // and set the additional include paths, if needed) exists

    filterRootLogonFullPath.Form("%s/PWG/muon/%s_rootlogon.C",alirootPath.Data(),filterName.Data());

    if ( gSystem->AccessPathName(filterRootLogonFullPath.Data()) )
    {
      std::cerr << "Could not find requested filter companion macro (looking into " << filterRootLogonFullPath.Data() << ")" << std::endl;
      return -4;
    }

    from.ReplaceAll(filterName.Data(),"");
    from.ReplaceAll("_WITH_ALIROOT_","");
    from.ReplaceAll(alirootVersion,"");
    from.ReplaceAll("..",".");

    std::cout << "Will filter file=" << from.Data() << std::endl;
  }
  
  if (from.Contains("alien://")) TGrid::Connect("alien://");

  if (!gGrid)
    {
      std::cerr << "Cannot get gGrid !" << std::endl;
      return -5;
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

    Int_t fail = gROOT->LoadMacro(Form("%s.C+",filterName.Data()));
      
    if ( fail )
    {
      std::cout << Form("Failed to load/compile macro %s+",filterMacroFullPath.Data()) << std::endl;
      return -6;
    }
        
    std::cout << Form("Will execute filter %s(\"%s\",\"%s\")",filterName.Data(),from.Data(),to.Data()) << std::endl;

    return (Int_t)gROOT->ProcessLine(Form("%s(\"%s\",\"%s\")",filterName.Data(),from.Data(),to.Data()));
  }
  else
  {
    // "normal" case of a simple copy
    //
    // ! operator since TFile::Cp returns kTRUE(1) in case of success
    std::cout << "Performing a simple TFile::Cp" << std::endl;
    return (!TFile::Cp(from.Data(),to.Data()));
  }
  
}
