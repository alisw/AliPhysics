// @(#) $Id$

#ifndef AliL3BenchmarkH
#define AliL3BenchmarkH

//_____________________________________________________________
//
// AliL3Benchmark
//
//   Benchmark class for level3 code
//  
//

#ifndef no_root
class TStopwatch;
class TString;
#else
class  AliL3Stopwatch;
#endif

class AliL3Benchmark {

public:
   AliL3Benchmark();
   virtual ~AliL3Benchmark();
   Int_t      GetBench(const char *name);
   void       Start(const char *name);
   void       Stop(const char *name);
   void       Analyze(const char* name);
   
   static Double_t GetCpuTime();

private:

   Int_t      fNbench;          //Number of active benchmarks
   Int_t      fNmax;            //Maximum number of benchmarks initialized
#ifndef no_root
   TString    *fNames;          //Names of benchmarks
   TStopwatch *fTimer;          //Timers
#else
   Char_t **fNames;             //Names of benchmarks
   AliL3Stopwatch *fTimer;      //Timers
#endif
   Float_t    *fSum;  //sum of time
   Float_t    *fMin;  //min of time
   Float_t    *fMax;  //max of time
   Int_t      *fCount;// counter

   ClassDef(AliL3Benchmark,0)  //L3 benchmark
};

#endif
