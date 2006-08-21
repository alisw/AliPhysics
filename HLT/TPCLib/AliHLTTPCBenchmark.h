// @(#) $Id$
// Original: AliL3Benchmark.h,v 1.6 2004/06/26 11:39:40 loizides 

#ifndef AliHLTTPCBenchmarkH
#define AliHLTTPCBenchmarkH

//_____________________________________________________________
//
// AliHLTTPCBenchmark
//
//   Benchmark class for level3 code
//  
//

#ifndef no_root
class TStopwatch;
class TString;
#else
class  AliHLTTPCStopwatch;
#endif

class AliHLTTPCBenchmark {

public:
   AliHLTTPCBenchmark();
   virtual ~AliHLTTPCBenchmark();
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
   AliHLTTPCStopwatch *fTimer;      //Timers
#endif
   Float_t    *fSum;  //sum of time
   Float_t    *fMin;  //min of time
   Float_t    *fMax;  //max of time
   Int_t      *fCount;// counter

   ClassDef(AliHLTTPCBenchmark,0)  //HLTTPC benchmark
};

#endif
