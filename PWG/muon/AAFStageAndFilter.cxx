#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TFile.h"
#include "TObjString.h"
#include "TMethodCall.h"

namespace AAF {
  
Int_t StageAndFilter(TString from, const TString& to,
                     TString filterName)
{
  ///
  /// Copy one file from url "from" into url "to", while filtering it at the same time.
  ///
  
  std::cout << "CpMacroWithFilter from : " << from.Data() << std::endl << " to : " << to.Data() << std::endl
	    << " with filter " << filterName.Data() << std::endl;

  if (from.IsNull() || to.IsNull()) return -1;

  if (from.Contains("alien://"))
  {
      TGrid::Connect("alien://");

      if (!gGrid)
      {
        std::cerr << "Cannot get gGrid !" << std::endl;
        return -2;
      }
  }

  if ( !filterName.IsNull() ) 
  {

    // check we have such a filter in our dictionary
    TMethodCall mc;
    
    mc.InitWithPrototype(Form("AAF::FILTER_%s",filterName.Data()),"char*,char*");
    
    if (!mc.IsValid())
    {
      std::cerr << "I don't know this function : AAF::FILTER_" << filterName.Data()
      << std::endl;
      return -3;
    }
    
    char** files = new char*[2];
    
    files[0] = new char[from.Length()+1];
    files[0][from.Length()]='\0';

    files[1] = new char[to.Length()+1];
    files[1][to.Length()]='\0';
    
    strcpy(files[0],from.Data());
    strcpy(files[1],to.Data());

    std::cout << Form("Will execute filter %s(\"%s\",\"%s\")",mc.GetMethodName(),from.Data(),to.Data()) << std::endl;
    
    mc.SetParamPtrs(files);
    
    mc.Execute();
    
    // some clean up
    
    delete[] files[0];
    delete[] files[1];
    delete[] files;
  }

  // we assume the function is successfull if the output file is present
  // and can be opened as a TFile ...
  
  TFile* f = TFile::Open(to.Data());
  
  if (f && f->IsOpen())
  {
    delete f;
    return 0;
  }
  else
  {
    delete f;
    std::cout << "Destination file " << to << " is not OK. Stage failed." << std::endl;
    return -9;
  }
}

}
