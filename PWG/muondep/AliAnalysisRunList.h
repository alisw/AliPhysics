#ifndef ALIANALYSISRUNLIST_H
#define ALIANALYSISRUNLIST_H

#include <vector>
#include <set>
#include <iostream>

/**

  @ingroup pwg_muondep_misc

  @class AliAnalysisRunList

  @brief Manipulates run lists

  @author Laurent Aphecetche (Subatech)
  */

class AliAnalysisRunList
{
    public:
        AliAnalysisRunList(const char* runlist);

        void Append(int runNumber);

        std::vector<int> AsVector() const;

        void Set(int runNumber);
        void Set(const std::vector<int>& runs);
        void Set(const std::set<int>& runs);
        void Set(const char* runlist);

        void Clear() { fRunList.clear(); }

        void Print() const;

        static void ReadIntegers(const char* filename, std::vector<int>& integers, bool resetVector=true);

        static void PrintIntegers(const std::vector<int>& integers, char sep = '\n', std::ostream& out = std::cout);

        static std::vector<int> GetRunsInCommon(const std::vector<int>& runlist1, const std::vector<int>& runlist2);

    private:

        std::vector<int> fRunList;
};

#endif

