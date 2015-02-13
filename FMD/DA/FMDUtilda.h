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
PostSummaries(AliFMDBaseDA& da, const char* prefix, Int_t runNumber, bool verb=true)
{
#ifndef ALI_AMORE
  std::cerr << "AMORE not found - nothing posted" << std::endl;
  return;
#else
  if (verb) 
    std::cout << "Posting " << prefix << " data to AMORE" << std::endl;
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
      if (verb)
	std::cout << " Posting " << prefix << "RuNoFMD" << det << " ... "
		  << std::flush;
      myAmore.Send(Form("%sRunNoFMD%d", prefix, det), runNo);
      if (verb) 
	std::cout << "done" << std::endl;
    }
    
    TIter     next(&da.GetSummaries());
    if (verb) {
      std::cout << " Will post summaries: " << std::endl;
      da.GetSummaries().ls();
    }
    TObject*  obj = 0;
    while ((obj = next())) {
      if (verb)
	std::cout << " Posting " << obj->GetName() 
		  << " ... " << std::flush;
      myAmore.Send(obj->GetName(), obj);
      if (verb) 
	std::cout << "done" << std::endl;
    }
    
  }
  catch (std::exception& e) {
    std::cerr << "Failed to make AMORE instance: " << e.what() << std::endl;
  }
  if (verb) 
    std::cout << "All objects posted to AMORE" << std::endl;
#endif
}

Int_t
UploadFiles(const char** files, const char** ids, bool verb=true)
{
  const char** pFile = files; 
  const char** pId   = ids; 

  if (verb) 
    std::cout << "Uploading files to FXS" << std::endl;
  Int_t ret = 0;
  while (*pFile && *pId) {
    if (verb) 
      std::cout << " Upload " << *pFile << " (" << *pId << ") ... " << std::flush;
    Int_t lret = daqDA_FES_storeFile(*pFile, *pId);
    if (lret != 0) {
      if (verb) 
	std::cout << "failed" << std::endl;
      std::cerr << "Failed to upload " << *pFile << " to FES id " 
		<< *pId << std::endl;
      ret++;
    }
    if (verb)
      std::cout << "done" << std::endl;
    pFile++;
    pId++;
  }
  if (verb)
    std::cout << "Done uploading with " << ret 
	      << " errors" << std::endl;
  return ret;
}



#endif // FMDUTILDA_H
// Local Variables: 
//  mode: C++
// End:

