/*
  Utility to do a parallel sync of remote ROOT files locally. As added value
  it allows to only sync the parts (TKeys) of a ROOT file which is match a
  provided regular expression.
 */

#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <cstdlib>
#include <cassert>
#include <string>
#include <cstdio>
#include <cstddef>
#include <stdint.h>
#include <vector>
#include <list>
#include <algorithm>
#include <sys/stat.h>
#include <limits.h>
#include <libgen.h>

#include "RConfig.h"
#include "TFile.h"
#include "THashList.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TTree.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TClass.h"
#include "TSystem.h"
#include "TGrid.h"
#include "TFile.h"

#include "TFileMerger.h"

#include "alisync.h"

static const char *USAGE =
      "usage: alisync [-h] [options] [@]<input file> [<input file>]\n"
      "Synchronizes remote ROOT files locally, atomically and in parallel\n"
      "If <input file> is prepended with @ it will be considered"
      " a list of other files to be merged.\n"
      "\n"
      "The following options can be specified:\n"
      "   -v [<level>]\n"
      "      Set verbosity level to <level>. If <level> is not provided defaults to maximum verbosity (99).\n"
      "   -o <output directory>\n"
      "      Set the output directory, i.e. where the files should be copied to. Defaults to the current one.\n"
      "   -j <# jobs>\n"
      "      Use <# jobs> parallel jobs to do the fetching (1 by default). \n"
      "   -t <# seconds>\n"
      "      Set per-file timeout while fetching to <# seconds> seconds.\n"
      "   -n <# retries>\n"
      "      Retry fetching for <# retries> retries before considering the file as not available.\n"
      "   -s <success rate>\n"
      "      Do not bother fetching more files after success rate has been reached. \n"
      "      E.g. if you are fine with only 90\% of the files being transferred, specify `-s 0.9'.\n"
      "   -i <regex>\n"
      "      Include TKeys whose name matches <regex>. This option can be used multiple times.\n";

const char * OPT_STRING = ":v:fj:t:n:s:i:o:";

static int gVerbosity = 0;
static bool gForce = false;

// Helper to print out verbosity aware messages.
void log(int level, const char *fmt, ...)
{
  if (level > gVerbosity)
    return;

  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
}

// Helper to print out a message and exit with exit code `exit`.
void dieIf(bool condition, int exitCode, const char *fmt, ...)
{
  if (!condition)
    return;
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  exit(exitCode);
}

// Helper for recursive directory creation
static void recursiveLocalMkdir(const char *dir) {
  char tmp[PATH_MAX];
  char *p = NULL;
  size_t len;

  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
          tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++)
          if(*p == '/') {
                  *p = 0;
                  mkdir(tmp, S_IRWXU);
                  *p = '/';
          }
  mkdir(tmp, S_IRWXU);
}


static std::vector<JobInfo> gJobs;

// Returns a string which contains the md5sum of the file.
std::string calculateMD5(const std::string &f) {
  std::string normalizedFilename = normalizePath(f);
  std::string command;
  std::string md5sum = "";
  std::string prefix("alien://");
  std::string out;

  if (normalizedFilename.compare(0, prefix.size(), prefix) == 0) {
    command = "gbbox md5sum ";
    command += normalizedFilename.substr(prefix.size());
  }else {
    command = "md5sum ";
    command += normalizedFilename;
  }

  command += " 2>/dev/null";
  out.resize(PATH_MAX*2+1);
  FILE *of = popen(command.c_str(), "r");
  if (!of)
    return "";
  size_t t = fread((void *)out.c_str(), PATH_MAX*2, 1, of);

  int ret = pclose(of);
  if(WIFEXITED(ret))
    return "";

  return out.substr(0, std::min(out.find("\t", 0), out.find(" ", 0)));
}

// Downloading of files:
// - Connect to alien if the filename starts with "alien".
// - Create a temporary destination file.
// - TFile::Cp the file to the output directory.
void downloadFile(const JobInfo &info, 
                  const std::string &outdir,
                  const std::vector<TRegexp *> &regexps) {
  gSystem->Load("libTreePlayer");

  static const std::string prefix = "alien://"; 
  // Keep track wether destination or source is actually alien
  bool isAlienDest = outdir.compare(0, prefix.size(), prefix) == 0;
  bool isAlienSrc = info.filename.compare(0, prefix.size(), prefix) == 0;

  if (isAlienSrc || isAlienDest)
    TGrid::Connect("alien://");

  char *relDir = strdup(dirname(strdup(info.filename.c_str() + (isAlienSrc ? prefix.size() : 0))));
  char *tmpFilename = strdup(basename(strdup(info.filename.c_str())));
  std::string tmpdest = outdir + "/" + relDir + "/." + tmpFilename + ".tmp";
  char *tmpoutdestCpy = strdup(tmpdest.c_str());
  char *tmpoutdir = strdup(dirname(tmpoutdestCpy));
  std::string dest = outdir + "/" + relDir + "/" + tmpFilename;

  // If destination is a local file, remove old temporary file.
  // else we really copy directly where we want to be.
  if (isAlienDest == false) {
    recursiveLocalMkdir(tmpoutdir);
    unlink(tmpdest.c_str());
  } else {
    assert(gGrid);
    assert(strncmp("alien://", tmpoutdir, prefix.size()) == 0);
    dieIf(gGrid->Mkdir(tmpoutdir+prefix.size(), "-p") == 0,
          1, "Unable to create %s\n", tmpoutdir);
    tmpdest = dest;
  }
  free(tmpoutdir);
  free(relDir);
  free(tmpFilename);

  // Make sure /./ is not part of the destination filename.
  std::string tmpsrc = normalizePath(info.filename);
  tmpdest = normalizePath(tmpdest);
  // If no regexps specified, simply use TFile::Cp. If the checksum between
  // source and dest exists and is the same we skip the copy.
  // If the source file was not modified after the last modification
  // date of the destination file, we skip copying it.
  // 
  // If regular expressions specified, copy only the tkeys that match those
  // files.
  if (regexps.empty())
  {
    std::string srcMD5 = calculateMD5(tmpsrc);
    log(2, "Source MD5: %s\n",srcMD5.c_str());
    std::string destMD5 = calculateMD5(dest);
    log(2, "Dest MD5: %s\n",destMD5.c_str());

    if (!gForce
        && !srcMD5.empty()
        && !destMD5.empty()
        && (srcMD5 == destMD5))
    {
      log(2, "Source and destination MD5 match. Not copying.\n");
      exit(0);
    }

    TFile *dest = TFile::Open(tmpdest.c_str());
    if (dest && !gForce)
    {
      TFile *src = TFile::Open(tmpsrc.c_str());
      if (!src)
        exit(1);
      bool sameSize = src->GetSize() == dest->GetSize();

      if (sameSize && (src->GetModificationDate() <= dest->GetModificationDate()))
      {
        log(2, "Destination file with more modern modification date. Not copying.\n");
        exit(0);
      }
    }
    if (!TFile::Cp(strdup(tmpsrc.c_str()), strdup(tmpdest.c_str())))
      exit(1);
  }
  else
  {
    TFile *src = TFile::Open(tmpsrc.c_str());
    TFile *dest = TFile::Open(tmpdest.c_str(), "RECREATE");
    TIter next(src->GetListOfKeys());
    TKey *key;
    std::vector<std::string> mergeableKeys;
    while ((key = dynamic_cast<TKey *>(next()))) {
      assert(key);
      const char *keyName = key->GetName();
      TString keyNameS(keyName);
      for (size_t rei = 0; rei < regexps.size(); ++rei)
      {
        int len;
        int status = regexps[rei]->Index(keyNameS, &len, 0);
        if (status >= 0 && len == keyNameS.Length())
        {
          log(1, "Key %s will be merged.\n", keyName);
          mergeableKeys.push_back(keyName);
        }
      }
    }
    for (std::vector<std::string>::iterator ki = mergeableKeys.begin();
         ki != mergeableKeys.end() ; ++ki)
    {
      TKey *oldKey = src->GetKey(ki->c_str());
      TObject *content = oldKey->ReadObj();
      TTree *tree = dynamic_cast<TTree *>(content);
      if (tree)
      {
        tree->SetBranchStatus("*",1);
        TTree *cloned = tree->CloneTree();
        dest->WriteTObject(cloned, oldKey->GetName());
        delete cloned;
      }
      else
        dest->WriteTObject(content, oldKey->GetName());
      delete content;
    }
    src->Close();
    dest->WriteKeys();
    dest->Close();
  }
  // In case the file is local, we need to move it in place.
  if (!isAlienDest)
    if (rename(tmpdest.c_str(), dest.c_str()))
    {
      unlink(tmpdest.c_str());
      exit(1);
    }
  exit(0);
}

void
killChildrenAndDie() {
  // Remember to kill all the children.
  exit(0);
}

static int gActiveJobs = 0;

void sigchld_hdl(int sig)
{
  int status;
  pid_t pid;

  while (true) {
    pid = wait3 (&status, WNOHANG, (struct rusage *)NULL);
    if (pid == 0)
      return;
    else if (pid == -1)
      return;
    else {
      std::vector<JobInfo>::iterator it = std::find_if(gJobs.begin(),
                                                       gJobs.end(),
                                                       SamePidAs(pid));
      it->retries--;
      int returnCode;
      if (WIFEXITED(status))
        returnCode = WEXITSTATUS(status);
      else if(WIFSIGNALED(status))
        returnCode = 143;
      else
        assert(false);
    
      if ((returnCode == 0) || (it->retries == 0))
        it->exitCode = returnCode;
      it->running = false;
      gActiveJobs--;
      log(2, "Job %i terminated with exit code %i. Attempts left: %i.\n",
              pid,
              returnCode,
              it->retries);
    }
  }
}

int
viableJobs() {
  for (size_t i = 0; i < gJobs.size(); ++i) {
    JobInfo &job = gJobs[i];
    // If running, do not bother.
    if (job.running)
      continue;
    // If not enough retries, skip it.
    if (!job.retries)
      continue;
    // If it completed successfully, skip it.
    if (job.exitCode == 0)
      continue;
    return i;
  }
  return -1;
}

// Check the array whether we are done or not.
// returns:
// - 0 completed successfully
// - 1 had an error
// - -1 not completed yet
int exitCode(float successrate) {
  int successfull = 0;
  int failed = 0;
  for (size_t i = 0; i < gJobs.size(); ++i) {
    JobInfo &info = gJobs[i];
    if (info.retries <= 0)
      failed++;
    if (info.exitCode == 0)
      successfull++;
  }
  // If we are fully successfull, we exit with 0.
  if (successfull == gJobs.size())
    return 0;

  // If we pass the requested success rate, we are done.
  if (((float)(successfull) / (float)gJobs.size()) > successrate)
    return 0;

  // If we have done all that we needed to do,
  // and we have errors, we exit with 1
  if (failed + successfull == gJobs.size())
    return 1;
  
  // Otherwise we carry on.
  return -1;
}

void
killLateJobs(int timeout) {
  time_t currentTime = time(0);
  for (size_t i = 0;  i < gJobs.size(); ++i)
  {
    JobInfo &info = gJobs[i];
    if (currentTime - info.startTime < timeout)
      continue;
    if (info.running)
    {
      log(2, "Process %i passed timeout, killing and trying again.", info.pid);
      kill(info.pid, SIGTERM);
    }
  }
}

//___________________________________________________________________________
int main( int argc, char **argv )
{
  int gNJobs = 1;
  int gTimeout = 60;
  int gRetries = 5;
  float gSuccessRate = 1.;
  float gSuccessfull = 0.;
  std::string gOutputDirectory = ".";
  std::vector<TRegexp *> gIncludeRE;
  typedef std::vector<std::string> Filenames;
  Filenames gFilenames;

  char ch;
  int intCand;
  float floatCand;

  while ((ch = getopt(argc, argv, OPT_STRING)) != -1) {
    switch (ch) {
      case 'h':
        std::cout << USAGE << std::endl;
        return 1;
      case 'n':
        intCand = atoi(optarg);
        dieIf(intCand <= 0, 1, "Invalid -n argument \"%s\".", optarg);
        gRetries = intCand;
        break;
      case 'o':
      {
        gOutputDirectory = normalizePath(trimTrailingSlashes(optarg));
        break;
      }
      case 'j':
        intCand = atoi(optarg);
        dieIf(intCand <= 0, 1, "Invalid -j argument \"%s\".", optarg);
        gNJobs = intCand;
        break;
      case 't':
        intCand = atoi(optarg);
        dieIf(intCand <= 0, 1, "Invalid -t argument \"%s\".", optarg);
        gTimeout = intCand;
        break;
      case 's':
        floatCand = atof(optarg);
        dieIf(intCand <= 0, 1, "Invalid -s argument \"%s\".", optarg);
        gSuccessRate = floatCand;
        break;
      case 'v':
        char *end;
        intCand = strtol(optarg, &end, 10);
        if (!*end || (intCand < 0)) {
          gVerbosity = 99;
          break;
        }
        gVerbosity = intCand; 
        break;
      case 'i':
        gIncludeRE.push_back(new TRegexp(optarg));
        break;
      case 'f':
        gForce = true;
        break;
      case ':':
        switch (optopt) {
          case 'v':
            gVerbosity = 99;
            break;
        }
        break;
      case '?':
        log(0, "Unknown option -%c.\n%s", optopt, USAGE);
        exit(1);
        break;
    }
  }

  dieIf(argc - optind < 1, 1, "Please, specify at least one file as input.\n\n%s", USAGE);
  
  // Get all the files to merge. 
  // @<filename> should contain a list of files
  for (int i = 0; i < argc - optind ; ++i)
  {
    const char *cf = argv[optind + i];
    if (*cf != '@')
    {
      std::string tcf(cf);
      gFilenames.push_back(trim(tcf));
      continue;
    }
    std::ifstream infile(cf+1);
    std::string line;
    while (std::getline(infile, line))
    {
      line = trim(line);
      if (line.empty())
        continue;
      gFilenames.push_back(trim(line));
    }
  }

  log(2, "Verbosity level: %i\n", gVerbosity);
  log(2, "Number of parallel jobs: %d\n", gNJobs);
  log(2, "Seconds before timeout: %d\n", gTimeout);
  log(2, "Number of retries: %d\n", gRetries);
  log(2, "Success rate %f\n", gSuccessRate);
  log(2, "Output directory: %s\n", gOutputDirectory.c_str());
  log(2, "Filenames:\n");
  for (size_t i = 0; i < gFilenames.size(); ++i)
    log(2, "- %s\n", gFilenames[i].c_str());

  // Setup signal handler for SIGCHLD.
  struct sigaction act;
  memset (&act, 0, sizeof(act));
  act.sa_handler = sigchld_hdl;

  if (sigaction(SIGCHLD, &act, 0)) {
    perror("sigaction");
    return 1;
  }

  // Create the table with all the jobs:
  for (size_t j = 0; j < gFilenames.size(); ++j) {
    JobInfo info;
    info.filename = gFilenames[j];
    info.retries = gRetries;
    info.pid = -1;
    info.exitCode = -1;
    info.running = false;
    gJobs.push_back(info);
  }

  // This loop take care or queueing files for fetching.
  // - Get the index to a viable job (i.e. one which still needs to be done).
  // - If we do not have enough free jobs sleep and try again.
  size_t nextJob;
  int gExitCode;
  while ((gExitCode = exitCode(gSuccessRate)) == -1)
  {
    // If a job takes more than timeout, kill it.
    killLateJobs(gTimeout);

    // When there is too many jobs simply try again in a while.
    if (gActiveJobs >= gNJobs)
    {
      log(2, "%i active jobs. Waiting for some fetcher to finish\n", gActiveJobs);
      sleep(1);
      continue;
    }

    // If there is no viable job, try again in a while, otherwise
    // process it.
    size_t nextJob = viableJobs();
    if (nextJob == -1)
    {
      sleep(1);
      continue;
    }
    JobInfo &current = gJobs[nextJob];

    pid_t pid = fork();
    if (pid)
    {
      log(2, "Scheduled fetching of %s to process %i\n", 
             current.filename.c_str(), pid);
      current.pid = pid;
      current.startTime = time(0);
      current.running = true;
      gActiveJobs++;
    }
    else
    {
      // Child we will not return from here.
      downloadFile(current, gOutputDirectory, gIncludeRE);
    }
  }

  return gExitCode;
}
