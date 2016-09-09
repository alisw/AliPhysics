/*

  This program will add histograms (see note) and Trees from a list of root files and write them
  to a target root file. The target file is newly created and must not be
  identical to one of the source files.

  Syntax:

       hadd targetfile source1 source2 ...
    or
       hadd -f targetfile source1 source2 ...
         (targetfile is overwritten if it exists)

  When the -f option is specified, one can also specify the compression
  level of the target file. By default the compression level is 1, but
  if "-f0" is specified, the target file will not be compressed.
  if "-f6" is specified, the compression level 6 will be used.

  For example assume 3 files f1, f2, f3 containing histograms hn and Trees Tn
    f1 with h1 h2 h3 T1
    f2 with h1 h4 T1 T2
    f3 with h5
   the result of
     hadd -f x.root f1.root f2.root f3.root
   will be a file x.root with h1 h2 h3 h4 h5 T1 T2
   where h1 will be the sum of the 2 histograms in f1 and f2
         T1 will be the merge of the Trees in f1 and f2

   The files may contain sub-directories.

  if the source files contains histograms and Trees, one can skip
  the Trees with
       hadd -T targetfile source1 source2 ...

  Wildcarding and indirect files are also supported
    hadd result.root  myfil*.root
   will merge all files in myfil*.root
    hadd result.root file1.root @list.txt file2. root myfil*.root
    will merge file1. root, file2. root, all files in myfil*.root
    and all files in the indirect text file list.txt ("@" as the first
    character of the file indicates an indirect file. An indirect file
    is a text file containing a list of other files, including other
    indirect files, one line per file).

  If the sources and and target compression levels are identical (default),
  the program uses the TChain::Merge function with option "fast", ie
  the merge will be done without  unzipping or unstreaming the baskets
  (i.e. direct copy of the raw byte on disk). The "fast" mode is typically
  5 times faster than the mode unzipping and unstreaming the baskets.

  NOTE1: By default histograms are added. However hadd does not support the case where
         histograms have their bit TH1::kIsAverage set.

  NOTE2: hadd returns a status code: 0 if OK, -1 otherwise

  Authors: Rene Brun, Dirk Geppert, Sven A. Schmidt, sven.schmidt@cern.ch
         : rewritten from scratch by Rene Brun (30 November 2005)
            to support files with nested directories.
           Toby Burnett implemented the possibility to use indirect files.
           Refactored by Giulio Eulisse in 2015.
 */

#include <unistd.h>
#include <cstdlib>
#include <cassert>
#include <string>
#include <cstdio>
#include <cstddef>
#include <algorithm>
#include <stdint.h>
#include "RConfig.h"
#include "TFile.h"
#include "THashList.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TClass.h"
#include "TSystem.h"

#include "TFileMerger.h"
#include "helpers.h"

static const char *USAGE = 
      "usage: alihadd [-h] [-f[kf123456]] [-a] [-k] [-T] [-O] [-n <max numer of files>]"
      " [-v [<level>]] <output file> [@]<input file> [<input file>]\n"
      "Merge objects found in <input file>s to <output files>\n"
      "If <input file> is prepended with @ it will be considered"
      " a list of other files to be merged.\n"
      "\n"
      "The following options can be specified:\n"
      "   -a\n" 
      "      Append to the output to file if it already exists.\n"
      "   -k\n"
      "      Do not exit on corrupt or non-existing input file but skip to the next one.\n"
      "   -T\n"
      "      Do not merge trees\n"
      "   -O\n"
      "      Optimize basket size when merging TTree"
      "   -n <max number of file>\n"
      "      Specify how many files can be used at once. Defaults to system maximum if not specified.\n"
      "   -r\n"
      "      Randomize input files order\n"
      "   -s <max size of file>\n"
      "      Specify the approximate size of an output file.\n"
      "   -v [<level>]\n"
      "      Set verbosity level to <level>. If <level> is not provided defaults to maximum verbosity (99).\n"
      "   -i <regex>\n"
      "      Include TKeys whose name matches <regex>. This option can be used multiple times.\n"
      "   -f[kf0123456]\n"
      "      Specify compression level between 0 and 6. Default is 1.\n"
      "      0 means no compression.\n"
      "      \"-ff\" will use the compression level of the first file.\n"
      "      \"-fk\" will keep the same compression level for the basket as the inputs unless \"-O\" is specified.\n"
      "      The meta data will be compressed using the compression level specified in the first\n"
      "      input or the compression setting specified follow fk (206 when using -fk206 for example)\n"
      "      If Target and source files have different compression settings\n"
      "      a slower method is used\n";

const char * OPT_STRING = ":hakTOri:n:s:v:f:";

static int gVerbosity = 0;

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
void die(int exitCode, const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  exit(exitCode);
}
//___________________________________________________________________________
int main( int argc, char **argv )
{
  Bool_t append = kFALSE;
  Bool_t force = kFALSE;
  Bool_t skip_errors = kFALSE;
  Bool_t reoptimize = kFALSE;
  Bool_t noTrees = kFALSE;
  Bool_t keepCompressionAsIs = kFALSE;
  Bool_t useFirstInputCompression = kFALSE;
  std::vector<TRegexp *> gIncludeRE;
  size_t gMaxFilesPerJob = 18446744073709551615UL;
  size_t gCostLimit = 18446744073709551615UL;
  bool gRandomizeInput = false;
  std::vector<std::string> mergeableKeys;

  Int_t newcomp = -1;

  char ch;
  int intCand;


  while ((ch = getopt(argc, argv, OPT_STRING)) != -1) {
    switch (ch) {
      case 'h':
        std::cout << USAGE << std::endl;
        return 1;
      case 'T': noTrees = kTRUE;
        break;
      case 'a': append = kTRUE;
        break;
      case 'k': skip_errors = kTRUE;
        break;
      case 'O': reoptimize = kTRUE;
        break;
      case 'n':
        intCand = atoi(optarg);
        if (intCand <= 0) 
          die(1, "Invalid -n argument \"%s\".", optarg);
        gMaxFilesPerJob = intCand;
        break;
      case 'r':
        gRandomizeInput = true;
        break;
      case 's':
        intCand = atoi(optarg);
        if (intCand <= 0) 
          die(1, "Invalid -s argument \"%s\".", optarg);
        gCostLimit = intCand;
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
      case ':':
        switch (optopt) {
          case 'v':
            gVerbosity = 99;
            break;
          case 'f':
            force = kTRUE;
            break;
        }
        break;
      case '?':
        log(0, "Unknown option -%c.\n%s", optopt, USAGE);
        exit(1);
        break;
      case 'f':
        switch(*optarg)
        {
          case 'k':
            keepCompressionAsIs = kTRUE;
            // FIXME: should add support for metadata compression level as
            //        suggested in hadd documentation.
            break;
          case 'f':
            useFirstInputCompression = kTRUE;
            break;
          default:
            force = kTRUE;
            intCand = strtol(optarg, &end, 10);
            log(99, "%i %i", *end, intCand);
            if (*end != 0 || intCand < 0 || intCand > 6)
              die(1, "Bad compression level %s\n", optarg);
            break;
        }
        break;
    }
  }

  log(2, "Verbosity level: %i\n", gVerbosity);

  gSystem->Load("libTreePlayer");

  if (argc - optind < 2)
    die(1, "Please, specify more than one file as input.\n\n%s", USAGE);

  // Payloads to be processed.
  std::vector<MergeJob> mergeJobs;
  std::vector<MergeInput> mergeInputs;

  // First file is always the output.
  const char *outputTemplate = argv[optind];
  log(1, "hadd Target file: %s\n", outputTemplate);

  // Get all the files to merge.
  // @<filename> should contain a list of files
  for (int i = 1; i < argc - optind ; ++i)
  {
    const char *cf = argv[optind + i];
    log(1, "Found option %s.\n", cf);
    if (*cf != '@')
    {
      log(1, "Adding file %s.\n", cf);
      mergeInputs.push_back(MergeInput(cf, rand()));
      continue;
    }
    std::ifstream infile(cf+1);
    std::string line;
    while (std::getline(infile, line))
    {
      log(1, "Adding file %s.\n", line.c_str());
      mergeInputs.push_back(MergeInput(line, rand()));
    }
  }

  assert(mergeInputs.size() > 0);
  // Get the list of all objects in all the files.
  // Calculate the approximate "cost" of the file
  std::vector<std::string> objects;
  std::vector<size_t> sizeCost;

  for (size_t i = 0; i < mergeInputs.size(); ++i)
  {
    MergeInput &input = mergeInputs[i];
    TFile *f = TFile::Open(input.filename.c_str());
    if (!f || f->IsZombie())
    {
      log(1, "Unable to open %s.\n", input.filename.c_str());
      if (f) 
        f->Close();
      continue;
    }
    std::vector<TIterator *> iterators;
    iterators.push_back(f->GetListOfKeys()->MakeIterator());
    while (!iterators.empty())
    {
      TIterator *keys = iterators.back();
      iterators.pop_back();
      while(TObject *obj = keys->Next())
      {
        TKey *key = dynamic_cast<TKey*>(obj);
        if (key)
        {
          const char *keyName = key->GetName();
          int cost  = key->GetNbytes();
          TString keyNameS(keyName); 

          if (gIncludeRE.empty())
          {
            log(1, "Key %s will be merged with a cost of %i.\n", keyName, cost);
            input.cost += cost;
            continue;
          }
          for (size_t rei = 0; rei < gIncludeRE.size(); ++rei)
          {
            int len;
            int status = gIncludeRE[rei]->Index(keyNameS, &len, 0);
            if (status >= 0 && len == keyNameS.Length())
            {
              log(1, "Key %s will be merged.\n", keyName);
              mergeableKeys.push_back(keyName);
              input.cost += cost;
            }
          }
        }
        else
          continue;
      }
      log(1, "Total cost for file %s is %i\n", input.filename.c_str(), input.cost);
      f->Close();
    }
  }

  if (useFirstInputCompression)
  {
    const char *filename = mergeInputs.front().filename.c_str();
    TFile *f = TFile::Open(filename);
    if (f && !f->IsZombie())
    {
      newcomp = f->GetCompressionSettings();
      log(1, "Compression level for file %s is %i. Using that for output.", filename, newcomp);
    }
    else
    {
      log(1, "Unable to determine compression level for first file %s. Defaulting to 1.", filename);
      newcomp = 1;
    }
    f->Close();
  }

  // Never use compression level -1
  if (newcomp == -1)
    newcomp = 1;

  if (gRandomizeInput) {
    std::sort(mergeInputs.begin(), mergeInputs.end(), MergeInputComparator());
  }

  // First calculate how many files we need per job.
  size_t costAcc = 0;
  std::vector<std::pair<size_t, size_t> > splitting;
  splitting.push_back(std::make_pair<size_t, size_t>(0,0));
  for (size_t mi = 0; mi < mergeInputs.size(); ++mi)
  {
    MergeInput &input = mergeInputs[mi];
    std::pair<size_t, size_t> &range = splitting.back();
    costAcc += input.cost;
    // If we are below threshold, we carry on adding to the same job, otherwise
    // we create a new MergeJob.
    range.second += 1;
    if ((costAcc < gCostLimit)
        || ((range.second - range.first) > gMaxFilesPerJob))
      continue;

    splitting.push_back(std::make_pair(range.second, range.second));
  }

  // Remove the last range in case it's empty.
  if (splitting.back().second - splitting.back().first == 0)
    splitting.pop_back();

  // Costruct the merge jobs using the splitting information.
  for (size_t si = 0; si < splitting.size(); ++si)
  {
    MergeJob newJob(outputTemplate, splitting.size());
    for(size_t ii = splitting[si].first ; ii < splitting[si].second; ++ii)
      newJob.inputs.push_back(ii);
    mergeJobs.push_back(newJob);
  }

  // FIXME: handle the case in which we need to keep old compression
  // FIXME: handle the case in which we need to force the new
  //        compression levels.
  log(1, keepCompressionAsIs && !reoptimize
         ? "hadd compression setting for meta data: %i\n"
         : "hadd compression setting for all ouput: %i\n", newcomp);

  for (size_t i = 0; i < mergeJobs.size(); ++i)
  {
    MergeJob &job = mergeJobs[i];

    TFileMerger merger(kFALSE,kFALSE);
    merger.SetMsgPrefix("hadd");
    merger.SetPrintLevel(gVerbosity - 1);

    // Make sure we merge only the keys we want.
    for (size_t mk = 0; mk < mergeableKeys.size(); ++mk)
      merger.AddObjectNames(mergeableKeys[mk].c_str());

    // Set the output file.
    if (append && !merger.OutputFile(job.output.c_str(),"UPDATE",newcomp))
      die(2, "hadd error opening target file for update : %s.\n", job.output.c_str());
    else if (!merger.OutputFile(job.output.c_str(), force, newcomp))
      die(1, "hadd error opening target file (does %s exists?).\n%s", job.output.c_str(), 
             force ? "" : "Pass \"-f\" argument to force re-creation of output file.\n");

    log(1, "Output file set");
    for (size_t i = 0; i < job.inputs.size() ; ++i)
    {
      MergeInput &input = mergeInputs[job.inputs[i]];
      bool success = merger.AddFile(input.filename.c_str());
      if (!success && skip_errors)
        std::cerr << "hadd skipping file with error: " << input.filename << std::endl;
      else if (!success)
        die(1, "hadd exiting due to error in %s.\n", input.filename.c_str());
    }

    if (reoptimize)
      merger.SetFastMethod(kFALSE);
    else if (!keepCompressionAsIs && merger.HasCompressionChange())
        // Don't warn if the user any request re-optimization.
        log(0, "hadd Sources and Target have different compression levels\n"
               "hadd merging will be slower\n");

    merger.SetNotrees(noTrees);
    Bool_t status;
    int doAppend = append ? TFileMerger::kIncremental : TFileMerger::kRegular;
    int doOnlyListed = gIncludeRE.empty() ? 0 : TFileMerger::kOnlyListed;
    status = merger.PartialMerge(TFileMerger::kAll | doAppend | doOnlyListed);

    if (!status) {
      log(1, "hadd failure during the merge of %i input files in %s\n", 
             merger.GetMergeList()->GetEntries(), job.output.c_str());
      exit(1);
    }
    log(1, "hadd merged %i input files in %s\n",
           merger.GetMergeList()->GetEntries(), 
           job.output.c_str());
  }
  return 0;
}
