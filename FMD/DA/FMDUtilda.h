#ifndef FMDUTILDA_H
#define FMDUTILDA_H
#include "AliFMDBaseDA.h"
#include "daqDA.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#ifdef ALI_AMORE
# include <AmoreDA.h>
# include <TH2.h>
# include <TSystem.h>
#endif

Bool_t
PostSummaries(AliFMDBaseDA& da, const char* prefix, Int_t runNumber, bool verb=true)
{
#ifndef ALI_AMORE
  std::cerr << "AMORE not found - nothing posted" << std::endl;
  return true;
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
    return false;
  }
  if (verb) 
    std::cout << "All objects posted to AMORE" << std::endl;
  return true;
#endif
}

#include <sys/types.h>
#include <sys/wait.h>
#include <cstring>
#include <cerrno>

int system_alternative(const char* pgm, char* const argv[])
{
  // See http://stackoverflow.com/questions/3055924/problems-with-system-calls-in-linux 
  pid_t pid = vfork();
  if (pid > 0) {
    // We're the parent, so wait for child to finish
    int status;
    waitpid(pid, &status, 0);
    return status;
  }
  else if (pid == 0) {
    // We're the child, so run the specified program.  Our exit status will
    // be that of the child program unless the execv() syscall fails.
    std::cout << "Execute \"" << pgm << "\"";
    for (int i = 0; argv[i]; i++) 
      std::cout << " \"" << argv[i] << "\"";  
    std::cout << std::endl;
    int status = execv(pgm, argv);
    if (status == -1) { 
      std::cout << "Failed to exec " << errno 
		<< " (" << strerror(errno) << ")" << std::endl;
    }
    return status;
  }
  else {
    // Something horrible happened, like system out of memory
    std::cout << "Failed to for fork " << errno 
	      << " (" << strerror(errno) << ")" << std::endl;
    return -1;
  }
}

Int_t myFesStore2(const char* file, const char* id, bool verb=true)
{
  TString path(gSystem->Getenv("DAQDALIB_PATH"));
  if (path.IsNull()) {
    std::cout << "DAQDALIB_PATH undefined" << std::endl;
    return -1;
  }
  TString cmd(path); 
  cmd.Append("/daqFES_store");
  std::cout << "Will Execute \"" << cmd 
	    << "\",\"" << file << "\",\"" << id
	    << "\"" << std::endl;
  char* const args[] = { 
    const_cast<char*>(cmd.Data()),
    const_cast<char*>(file), 
    const_cast<char*>(id), 
    0 };
  int ret = system_alternative(cmd.Data(),args);
  if (WIFEXITED(ret)) {
    ret = WEXITSTATUS(ret);
    if (ret != 0) 
      std::cout << "Failed to execute \"" << cmd << "\" -> " 
		<< ret << std::endl;
  }
  else if (WIFSIGNALED(ret)) {
    int sig = WTERMSIG(ret);
    std::cout << "Process \"" << cmd << "\" received signal "
	      << sig << " (" << strsignal(sig) << ")" << std::endl;
    ret = -1;
  }
  return ret;
}
  
Int_t myFesStore(const char* file, const char* id, bool verb=true)
{
  TString path(gSystem->Getenv("DAQDALIB_PATH"));
  if (path.IsNull()) {
    std::cout << "DAQDALIB_PATH undefined" << std::endl;
    return -1;
  }
  // Check if we can get a shell 
  Int_t ret = system(NULL);
  if (ret == 0) { 
    std::cout << "Cannot spawn a shell" << std::endl;
    return -2;
  }
  // Biuld command 
  TString cmd(path); 
  cmd.Append("/daqDetDB_store ");
  cmd.Append(file);
  cmd.Append(" ");
  cmd.Append(id);
  // Execute the command 
  std::cout << "Will now execute \"" << cmd << "\"" << std::endl;
  ret = system(cmd.Data());
  if (WIFEXITED(ret)) {
    ret = WEXITSTATUS(ret);
    if (ret != 0) 
      std::cout << "Failed to execute \"" << cmd << "\" -> " 
		<< ret << std::endl;
  }
  else if (WIFSIGNALED(ret)) {
    int sig = WTERMSIG(ret);
    std::cout << "Process \"" << cmd << "\" received signal "
	      << sig << " (" << strsignal(sig) << ")" << std::endl;
    ret = -1;
  }
  return ret;
}

Int_t
UploadFiles(const char** files, const char** ids, bool useOwn, bool verb=true)
{
  const char** pFile = files; 
  const char** pId   = ids; 

  if (verb) 
    std::cout << "Uploading files to FXS" << std::endl;
  Int_t ret = 0;
  while (*pFile && *pId) {
    if (verb) 
      std::cout << " Upload " << *pFile << " (" << *pId << ") ... " << std::flush;
    Int_t lret = -1;
    if (useOwn) lret = myFesStore2(*pFile, *pId);
    else        lret = daqDA_FES_storeFile(*pFile, *pId);

    if (lret != 0) {
      if (verb) 
	std::cout << "failed (" << lret << ")" << std::endl;
      std::cerr << "Failed to upload " << *pFile << " to FES id " 
		<< *pId << " - exit " << lret << std::endl;
      ret++;
    }
    else {
     if (verb) std::cout << "done" << std::endl;
     std::ofstream touch(Form("%s.upload", *pFile));
     touch << "Done" << std::endl;
     touch.close();
    }
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

