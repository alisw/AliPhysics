#include "AliRawReader.h"
#include "Riostream.h"
#include "TDatime.h"
#include "TString.h"
#include "TSystem.h"

struct EventInfo {
  Int_t mEventNumber;
  Int_t mTimeStamp;
};

void GenerateEntryListOneChunk(const char* infile, std::ostream& outfile, const TDatime& fromTime,
                               const TDatime& toTime)
{
  std::vector<EventInfo> eventNumbers;
  int n(0);

  std::unique_ptr<AliRawReader> rawReader(AliRawReader::Create(infile));

  while (rawReader->NextEvent()) {
    Int_t ts = rawReader->GetTimestamp();
    TDatime d(ts);
    if (d <= toTime && d >= fromTime) {
      EventInfo ei;
      ei.mEventNumber = n;
      ei.mTimeStamp = ts;
      eventNumbers.push_back(ei);
    }
    ++n;
  }
  if (!eventNumbers.empty()) {
    std::sort(begin(eventNumbers), end(eventNumbers),
              [](EventInfo& i, EventInfo& j) -> bool { return i.mTimeStamp < j.mTimeStamp; });
    outfile << infile << std::endl << eventNumbers.size() << std::endl;
    for (auto ei : eventNumbers) {
      outfile << ei.mEventNumber << " ";
    }
    outfile << std::endl;
  }
}

Int_t GenerateEntryList(const char* filelist, const char* entrylistfilename, const TDatime& fromTime,
                        const TDatime& toTime)
{

  std::string line;
  std::ifstream in(filelist);
  std::vector<std::string> files;

  while (std::getline(in, line)) {
    files.push_back(line.c_str());
  }

  std::ofstream entryListFile(gSystem->ExpandPathName(entrylistfilename));

  for (auto f : files) {
    GenerateEntryListOneChunk(f.c_str(), entryListFile, fromTime, toTime);
  }

  return 0;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Usage : " << argv[0]
              << " --collection filelist.txt --entrylist filename.txt --from [startTime] --to [endTime]" << std::endl;
    return 1;
  }

  TString filelist{ "" };
  TDatime fromTime{ "1995-1-1 00:00:00" };
  TDatime toTime{ "2037-1-1 00:00:00" };
  TString entryList{ "" };

  for (int i = 0; i < argc; ++i) {
    if (TString(argv[i]) == "--entrylist") {
      entryList = argv[i + 1];
    }
    if (TString(argv[i]) == "--filelist") {
      filelist = argv[i + 1];
    }
    if (TString(argv[i]) == "--from") {
      fromTime = TDatime(argv[i + 1]);
    }
    if (TString(argv[i]) == "--to") {
      toTime = TDatime(argv[i + 1]);
    }
  }

  return GenerateEntryList(filelist.Data(), entryList.Data(), fromTime, toTime);
}
