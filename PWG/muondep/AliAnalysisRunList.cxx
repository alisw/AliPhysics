#include "AliAnalysisRunList.h"

#include <algorithm>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>
#include "TSystem.h"
#include <fstream>

/// \ingroup utils

//______________________________________________________________________________
AliAnalysisRunList::AliAnalysisRunList(const char* runlist)
{
    Set(runlist);
}

//______________________________________________________________________________
void AliAnalysisRunList::Append(int runNumber)
{
    fRunList.push_back(runNumber);
    std::set<int> aset(fRunList.begin(),fRunList.end());
    fRunList.assign(aset.begin(),aset.end());
}

//______________________________________________________________________________
std::vector<int> AliAnalysisRunList::AsVector() const
{
    return fRunList;
}

//______________________________________________________________________________
std::vector<int> AliAnalysisRunList::GetRunsInCommon(const std::vector<int>& runlist1,
        const std::vector<int>& runlist2)
{
    /// Build a runlist with the runs that are both in runlist1 and runlist2
    std::vector<int> common;

    // insure both initial vectors are sorted to be able to use the set_intersection algorithm
    std::set<int> s1(runlist1.begin(),runlist1.end());
    std::set<int> s2(runlist2.begin(),runlist2.end());

    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(common));

    return common;
}

//______________________________________________________________________________
void AliAnalysisRunList::Print() const
{
    PrintIntegers(fRunList);
}

//______________________________________________________________________________
void AliAnalysisRunList::PrintIntegers(const std::vector<int>& integers,
                                              char sep,
                                              std::ostream& out)
{
  /// print a list of integers
  for ( std::vector<int>::size_type i = 0; i < integers.size(); ++i )
  {
    out << integers[i] << sep;
  }
  out << std::endl;
}

//______________________________________________________________________________
void AliAnalysisRunList::ReadIntegers(const char* filename,
                                             std::vector<int>& integers,
                                             Bool_t resetVector)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  
  if ( gSystem->AccessPathName(gSystem->ExpandPathName(filename))==kTRUE )
  {
    return;
  }
  std::ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
  if (!resetVector)
  {
    for ( std::vector<int>::size_type j = 0; j < integers.size(); ++ j )
    {
      runset.insert(integers[j]);
    }
  }
  
  char line[10000];
  
  in.getline(line,10000,'\n');
  
  TString sline(line);
  
  if (sline.Contains(","))
  {
    TObjArray* a = sline.Tokenize(",");
    TIter next(a);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    delete a;
  }
  else
  {
    runset.insert(sline.Atoi());
    
    while ( in >> i )
    {
      runset.insert(i);
    }
  }
  
  integers.clear();
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
  {
    integers.push_back((*it));
  }
  
  std::sort(integers.begin(),integers.end());
}

//______________________________________________________________________________
void AliAnalysisRunList::Set(Int_t runNumber)
{
  // Make the runlist be a single run
  fRunList.clear();
  fRunList.push_back(runNumber);
}

//______________________________________________________________________________
void AliAnalysisRunList::Set(const std::vector<int>& runs)
{
  // Make the runlist be a single run
  fRunList = runs;
}

//______________________________________________________________________________
void AliAnalysisRunList::Set(const std::set<int>& runs)
{
  // Make the runlist be a single run
  fRunList.clear();
  for ( std::set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    fRunList.push_back(*it);
  }
}

//______________________________________________________________________________
void AliAnalysisRunList::Set(const char* runlist)
{
  // Read the runlist from an ASCII file or a comma separated list
  // or a space separated list
  
  fRunList.clear();
  
  if ( TString(runlist).Contains(",") || TString(runlist).Contains(" ") )
  {
    TObjArray* runs = 0x0;
    if ( TString(runlist).Contains(",") )
    {
      runs = TString(runlist).Tokenize(",");
    }
    else
    {
      runs = TString(runlist).Tokenize(" ");
    }
    TIter next(runs);
    TObjString* s;
    std::set<int> runset;
    
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    
    for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
    {
      fRunList.push_back((*it));
    }
    
    std::sort(fRunList.begin(),fRunList.end());
    
    delete runs;
  }
  else
  {
    ReadIntegers(runlist,fRunList);
    if (fRunList.empty())
    {
      if ( TString(runlist).IsDigit() )
      {
        Set(TString(runlist).Atoi());
      }
      else
      {
          std::cerr << "Could not set run list !" << std::endl; 
      }
    }
  }
}

