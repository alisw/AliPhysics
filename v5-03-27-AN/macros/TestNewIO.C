Int_t TestNewIO(Int_t n = 5,Char_t SlowOrFast='s')
{
   Int_t rc=0;
   
   // Debug and log level
   AliLog::SetGlobalDebugLevel(10);
   AliLog::SetGlobalLogLevel(10);

   /**********************************************/
   /************ G E N E R A T I O N *************/
   /**********************************************/
      
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun(n);


   /**********************************************/
   /******************* T P C ********************/
   /**********************************************/
   //The following part is just core of AliTPCtest.C
 
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2SDigits.C");
   if (rc=AliTPCHits2SDigits(n)) return rc;
  
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCSDigits2Digits.C");
   if (rc=AliTPCSDigits2Digits(n)) return rc;
   
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters(n)) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks(n)) return rc;
  
   /**********************************************/
   /******************* I T S ********************/
   /**********************************************/
   //The following part is just core of AliITStestV2.C
   if (SlowOrFast=='f') {
      cerr<<"Fast AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2FastRecPoints.C");
      AliITSHits2FastRecPoints();
   } else {
      cerr<<"Slow AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2SDigits.C");
      AliITSHits2SDigits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSSDigits2Digits.C");
      AliITSSDigits2Digits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSDigits2RecPoints.C");
      AliITSDigits2RecPoints();
   }
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
   if (rc=AliITSFindClustersV2(SlowOrFast)) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindTracksV2.C");
   if (rc=AliITSFindTracksV2()) return rc;


   /**********************************************/
   /****************** P H O S *******************/
   /**********************************************/

   if (gAlice)
    {
      delete AliRunLoader::Instance();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }
   AliRunLoader* rl =AliRunLoader::Open("galice.root","EVENTFOLDER");
   AliPHOSReconstructioner r("EVENTFOLDER");
   r.ExecuteTask("deb");
   delete rl;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCComparison.C");
   if (rc=AliTPCComparison()) return rc;
   
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSComparisonV2.C");
   if (rc=AliITSComparisonV2()) return rc;



   ::Info("NewIO test","Everything seems to be OK");   
   ::Info("NewIO test","You can try now display.C and TPC/AliTPCDisplayDigits3Dnew.C");
}
