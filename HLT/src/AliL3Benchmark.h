#ifndef AliL3_Benchmark
#define AliL3_Benchmark

#ifndef no_root
#include <Rtypes.h>
class TStopwatch;
class TString;
#else
#include "AliL3RootTypes.h"
#include "AliL3Stopwatch.h"
#endif

class AliL3Benchmark {

private:

   Int_t      fNbench;          //Number of active benchmarks
   Int_t      fNmax;            //Maximum number of benchmarks initialized
#ifndef no_root
   TString    *fNames;          //Names of benchmarks
   TStopwatch *fTimer;          //Timers
#else
   Char_t **fNames;
   AliL3Stopwatch *fTimer;
#endif
   Float_t    *fSum;
   Float_t    *fMin;
   Float_t    *fMax;
   Int_t      *fCount;
   //TStopwatch *fStopwatch;    //Stopwatch

public:
   AliL3Benchmark();
   virtual ~AliL3Benchmark();
   Int_t      GetBench(const char *name);
   void       Start(const char *name);
   void       Stop(const char *name);
   void       Analyze(const char* name);
   
   static Double_t GetCpuTime();
   
   ClassDef(AliL3Benchmark,0)  //L3 benchmark
};

#endif
