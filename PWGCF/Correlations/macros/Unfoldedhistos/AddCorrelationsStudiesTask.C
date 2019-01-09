#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TROOT.h"
#include "AliLog.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCorrelationsStudies.h"
#endif

AliAnalysisTaskSE *AddCorrelationsStudiesTask(const char *mincenstr, const char *maxcenstr, const char *configstring, const char *corrconfigstring, const char *corrbinning) {

  AliAnalysisTaskCorrelationsStudies* taskCS;
  AliAnalysisManager    *mgr = AliAnalysisManager::GetAnalysisManager();
  TString outfilename = mgr->GetCommonFileName();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  if (mgr != NULL) {

    char suffix[1024];
    char configbuffer[1024];

    TString mincenstring = mincenstr;
    TString maxcenstring = maxcenstr;

    TObjArray *mincentokens = mincenstring.Tokenize("[],");
    TObjArray *maxcentokens = maxcenstring.Tokenize("[],");

    if (mincentokens->GetEntries() != maxcentokens->GetEntries() + 1) {
      /* the first item in mincentokens is the centrality detector qualifier */
      AliFatalGeneral("AddCorrelationsStudiesTask", "number of centrality cuts edges differ. ABORTING");
    }

    int centqualifier = ((TObjString*) mincentokens->At(0))->String().Atoi();
    int *mincen = new int[maxcentokens->GetEntries()];
    int *maxcen = new int[maxcentokens->GetEntries()];

    for (int icen = 0; icen < maxcentokens->GetEntries(); icen++) {
      mincen[icen] = ((TObjString*) mincentokens->At(icen+1))->String().Atoi();
      maxcen[icen] = ((TObjString*) maxcentokens->At(icen))->String().Atoi();
    }

    for (int icen = 0; icen < maxcentokens->GetEntries(); icen++) {
      int centclass;
      int lowedge;
      int upedge;
      /// Sets the type of centrality to handle
      /// \param ctype the centrality type
      ///    |code| detector for centrality estimate, range modifier |
      ///    |:--:|--------|
      ///    |  0 | no centrality cut |
      ///    |  1 | default detector for the concerned system, cut in the range 0-100% in steps of 10% |
      ///    |  2 | alternative detector for the concerned system, cut in the range 0-100% in steps of 10% |
      ///    |  3 | default detector for the concerned system, cut in the range 0-50% in steps of 5% |
      ///    |  4 | default detector for the concerned system, cut in the range 50-100% in steps of 5% |
      ///    |  5 | default detector for the concerned system, cut in the range 0-10% in steps of 1% |
      ///    |  6 | default detector for the concerned system, cut in the range 10-20% in steps of 1% |
      ///    |  7 | alternative detector for the concerned system, cut in the range 0-50% in steps of 5% |
      ///    |  8 | alternative detector for the concerned system, cut in the range 50-100% in steps of 5% |
      if (mincen[icen] < 500 && ((int(mincen[icen] / 100)*100 != mincen[icen]) || (int(maxcen[icen] / 100)*100 != maxcen[icen]))) {
        /* short centrality range */
        switch(centqualifier){
        case 1:
        case 3:
          centclass = 3;
          lowedge = int(mincen[icen]/50);
          upedge = int(maxcen[icen]/50);
          break;
        case 2:
        case 7:
          centclass = 7;
          lowedge = int(mincen[icen]/50);
          upedge = int(maxcen[icen]/50);
          break;
        case 5:
          centclass = 5;
          lowedge = mincen[icen];
          upedge = maxcen[icen];
          if (upedge == 10) upedge = 0;
          break;
        case 6:
          centclass = 6;
          lowedge = mincen[icen] - 10;
          upedge = maxcen[icen] - 10;
          if (upedge == 10) upedge = 0;
          break;
        default:
          AliFatalGeneral("AddCorrelationsStudiesTask", "centrality qualifier not supported. ABORTING");
        }
      }
      else {
        switch(centqualifier){
        case 1:
        case 3:
          centclass = 1;
          lowedge = int(mincen[icen]/100);
          upedge = int(maxcen[icen]/100);
          break;
        case 2:
        case 7:
          centclass = 2;
          lowedge = int(mincen[icen]/100);
          upedge = int(maxcen[icen]/100);
          break;
        default:
          AliFatalGeneral("AddCorrelationsStudiesTask", "centrality qualifier not supported. ABORTING");
        }
      }
      sprintf(suffix,"%dVo%d",int(mincen[icen]),int(maxcen[icen]));
      sprintf(configbuffer,configstring,centclass,lowedge,upedge);

      TString szContainerPrefix;

      taskCS = new AliAnalysisTaskCorrelationsStudies(Form("Correlation Studies %s", suffix));
      taskCS->Configure(configbuffer);
      if (TString(corrconfigstring).Contains("simulate"))
        taskCS->ConfigureCorrelations(corrconfigstring,Form("density_%s_", suffix),szContainerPrefix);
      else
        taskCS->ConfigureCorrelations(corrconfigstring,Form("correction_%s_", suffix),szContainerPrefix);
      taskCS->ConfigureCorrelationsBinning(corrbinning);
      mgr->AddTask(taskCS);

      // create containers for input/output
      AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("CorrelationStudies_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
      AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("CSDptDptCorrelations_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
      AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("CSDptDptCorrelationsMCRecOptions_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
      AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("CSDptDptTrueCorrelations_%s_%s", szContainerPrefix.Data(), suffix), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);

      // connect input/output
      mgr->ConnectInput(taskCS, 0, cinput);
      mgr->ConnectOutput(taskCS, 1, coutput1);
      mgr->ConnectOutput(taskCS, 2, coutput2);
      mgr->ConnectOutput(taskCS, 3, coutput3);
      mgr->ConnectOutput(taskCS, 4, coutput4);
    }
    /* cleaning up */
    delete [] mincen;
    delete [] maxcen;
    delete mincentokens;
    delete maxcentokens;
    /* warning we return the last one so, it should not be used */
    return taskCS;
  }
  else {
    AliFatalGeneral("AddCorrelationsStudiesTask", "Add task needs a previously started analysis manager");
    return NULL;
  }
}
