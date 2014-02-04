//
//  CpMacroWithFilter.C
//
//  Created by Laurent Aphecetche on 30-sep-2013.
//

Int_t CpMacroWithFilter(TString from, TString to)
{
  /// Copy one file from url "from" into url "to".
  ///
  /// If from contains one of the filter keywords
  /// (in the form of e.g. AliAOD.FILTER_ZZZZ.root)
  /// then we call the external filter macro
  ///
  /// otherwise we do a simple TFile::Cp(from,to)
  ///
  
  if (from.IsNull() || to.IsNull()) return 1;
  
  std::cout << Form("Entering CpMacroWithFilter(\"%s\",\"%s\")",from.Data(),to.Data()) << std::endl;
  
  TObjArray filters;
  filters.SetOwner(kTRUE);
  
  filters.Add(new TObjString("FILTER_AODMUONWITHTRACKLETS"));
  filters.Add(new TObjString("FILTER_RAWMUON"));
  
  TIter next(&filters);
  TObjString* s;
  TString filter;
  
  while ( ( s = static_cast<TObjString*>(next())) )
  {
      if ( from.Contains(s->String()) )
      {
        filter = s->String();
        break;
      }
  }
  
  if (from.Contains("alien:\/\/")) TGrid::Connect("alien:\/\/");

  if ( filter.Length() > 0 )
  {
    // check the required filter is actually available
    
    from.ReplaceAll(filter.Data(),"");
    from.ReplaceAll("..",".");
    
    if ( gSystem->AccessPathName(Form("%s/etc/%s.C",gSystem->Getenv("XRDDMSYS"),filter.Data()) ) )
    {
      std::cout << Form("ERROR: could not find a filter named %s.C",filter.Data()) << std::endl;
      return 2;
    }
    else
    {
      // check also we have a companion macro (to load the relevant libraries and
      // set the additional include paths, if needed)
      
      if ( gSystem->AccessPathName(Form("%s/etc/%s_rootlogon.C",gSystem->Getenv("XRDDMSYS"),filter.Data()) ) )
      {
        std::cout << Form("ERROR: could not find the companion macro %s_rootlogon.C for the filter named %s.C",filter.Data(),filter.Data()) << std::endl;
        return 3;
      }

      else
      {
        // most probably the filter will require AliRoot libs, so add the dynamic path here
        // as well as the include path and the macro path.
        //
        gSystem->AddDynamicPath(Form("%s/aliroot/lib/tgt_%s",gSystem->Getenv("ALICE_PROOF_AAF_DIR"),gSystem->GetBuildArch()));
        gSystem->SetIncludePath(Form("-I%s/etc -I%s/aliroot/include",gSystem->Getenv("XRDDMSYS"),gSystem->Getenv("ALICE_PROOF_AAF_DIR")));
        gROOT->SetMacroPath(Form("%s:%s/etc",gROOT->GetMacroPath(),gSystem->Getenv("XRDDMSYS")));
        
//        gSystem->AddDynamicPath(Form("/pool/PROOF-AAF/aliroot/lib/tgt_%s",gSystem->GetBuildArch()));
//        gSystem->SetIncludePath("-I/pool/PROOF-AAF/xrootd_1.0.50/etc -I/pool/PROOF-AAF/aliroot/include");
//        gROOT->SetMacroPath(Form("%s:%s/etc",gROOT->GetMacroPath(),gSystem->Getenv("XRDDMSYS")));

        // execute the companion macro

        std::cout << Form("Will load companion macro %s_rootlogon.C(\"%s\",\"%s\")",filter.Data(),from.Data(),to.Data()) << std::endl;

        gROOT->Macro(Form("%s_rootlogon.C",filter.Data()));
        
        std::cout << gSystem->GetIncludePath() << std::endl;
        
        // finally delegate the work to the required filter
      
        std::cout << Form("Will compile filter %s.C+(\"%s\",\"%s\")",filter.Data(),from.Data(),to.Data()) << std::endl;

        Int_t fail = gROOT->LoadMacro(Form("%s.C+",filter.Data()));
      
        if ( fail )
        {
          std::cout << Form("Failed to load/compile macro %s.C+",filter.Data()) << std::endl;
          return 4;
        }
        
        std::cout << Form("Will execute filter %s(\"%s\",\"%s\")",filter.Data(),from.Data(),to.Data()) << std::endl;

        return (Int_t)gROOT->ProcessLine(Form("%s(\"%s\",\"%s\")",filter.Data(),from.Data(),to.Data()));
      }
    }
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
