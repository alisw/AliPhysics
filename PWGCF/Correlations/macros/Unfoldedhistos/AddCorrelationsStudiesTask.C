#include "Riostream.h"
#include "TROOT.h"
#include "AliLog.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCorrelationsStudies.h"

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

    if (mincentokens->GetEntries() != maxcentokens->GetEntries()) {
      AliFatalGeneral("AddCorrelationsStudiesTask", "number of centrality cuts edges differ. ABORTING");
    }

    int *mincen = new int[mincentokens->GetEntries()];
    int *maxcen = new int[mincentokens->GetEntries()];

    for (int icen = 0; icen < mincentokens->GetEntries(); icen++) {
      mincen[icen] = ((TObjString*) mincentokens->At(icen))->String().Atoi();
      maxcen[icen] = ((TObjString*) maxcentokens->At(icen))->String().Atoi();
    }

    for (int icen = 0; icen < mincentokens->GetEntries(); icen++) {
      int centclass;
      int lowedge;
      int upedge;
      if (mincen[icen] < 500 && ((int(mincen[icen] / 100)*100 != mincen[icen]) || (int(maxcen[icen] / 100)*100 != maxcen[icen]))) {
        /* short centrality range */
        centclass = 3;
        lowedge = int(mincen[icen]/50);
        upedge = int(maxcen[icen]/50);
      }
      else {
        centclass = 1;
        lowedge = int(mincen[icen]/100);
        upedge = int(maxcen[icen]/100);
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
