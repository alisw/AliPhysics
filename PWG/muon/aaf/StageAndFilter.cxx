#include "Riostream.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TFile.h"
#include "TObjString.h"
#include "TMethodCall.h"

namespace AAF {
  
Int_t StageAndFilter(const TString& from, const TString& to,
                     const TString& filterName, Int_t verboseLevel)
{
  ///
  /// Copy one file from url "from" into url "to", while filtering it at the same time.
  ///
  
	if ( verboseLevel > 0 )
	{
		std::cout << "StageAndFilter from : " << from.Data() << std::endl << " to : " << to.Data() << std::endl
				<< " with filter " << filterName.Data() << std::endl;
	}

  if (from.IsNull())
  {
	  std::cerr << "StageAndFilter : cannot stage from a null source..." << std::endl;
	  return -1;
  }

  if (to.IsNull())
  {
	  std::cerr << "StageAndFilter : cannot stage to a null destination..." << std::endl;
	  return -2;
  }

  if (from.Contains("alien://"))
  {
      TGrid::Connect("alien://",0,0,"t");

      if (!gGrid)
      {
        std::cerr << "Cannot get gGrid !" << std::endl;
        return -3;
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
      return -4;
    }
    
    char** files = new char*[2];
    
    files[0] = new char[from.Length()+1];
    files[0][from.Length()]='\0';

    files[1] = new char[to.Length()+1];
    files[1][to.Length()]='\0';
    
    strcpy(files[0],from.Data());
    strcpy(files[1],to.Data());

    if  (verboseLevel > 0 )
    {
    	std::cout << Form("Will execute filter %s(\"%s\",\"%s\")",mc.GetMethodName(),from.Data(),to.Data()) << std::endl;
    }
    
    mc.SetParamPtrs(files);
    
    mc.Execute();
    
    // some clean up
    
    delete[] files[0];
    delete[] files[1];
    delete[] files;
  }
  else
  {
     // no filter => plain copy
	  TFile::Cp(from.Data(),to.Data());
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
