#ifndef FMDUTILDA_H
#define FMDUTILDA_H
#include "AliFMDBaseDA.h"
#include "daqDA.h"
#include <stdexcept>
#include <iostream>
#ifdef ALI_AMORE
# include <AmoreDA.h>
# include <TH2.h>
# include <TSystem.h>
#endif

void
PostSummaries(AliFMDBaseDA& da, const char* prefix, Int_t runNumber)
{
#ifndef ALI_AMORE
  std::cerr << "AMORE not found - nothing posted" << std::endl;
  return;
#else
  try { 
    TString daName(gSystem->Getenv("AMORE_DA_NAME"));
    if (daName.IsNull()) { 
      daName = gSystem->Getenv("DATE_ROLE_NAME");
      if (daName.IsNull()) 
	gSystem->Setenv("AMORE_DA_NAME", "FMD");
    }
    amore::da::AmoreDA myAmore(amore::da::AmoreDA::kSender);

    for (UShort_t det = 1; det <= 3; det++) {
      if (!da.HasSeenDetector(det)) continue;
      TObject* runNo = new TObject;
      runNo->SetUniqueID(runNumber);
      myAmore.Send(Form("%sRunNoFMD%d", prefix, det), runNo);
    }
    
    TIter     next(&da.GetSummaries());
    TObject*  obj = 0;
    while ((obj = next())) 
      myAmore.Send(obj->GetName(), obj);
    
  }
  catch (std::exception& e) {
    std::cerr << "Failed to make AMORE instance: " << e.what() << std::endl;
  }
#endif
}

Int_t
UploadFiles(const char** files, const char** ids)
{
  const char** pFile = files; 
  const char** pId   = ids; 

  Int_t ret = 0;
  while (*pFile && *pId) { 
    Int_t lret = daqDA_FES_storeFile(*pFile, *pId);
    if (lret != 0) {
      std::cerr << "Failed to upload " << *pFile << " to FES id " 
		<< *pId << std::endl;
      ret++;
    }
    pFile++;
    pId++;
  }
  return ret;
}



#endif // FMDUTILDA_H
// Local Variables: 
//  mode: C++
// End:

