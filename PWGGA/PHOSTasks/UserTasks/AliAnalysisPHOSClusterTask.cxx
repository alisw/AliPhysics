#include <TChain.h>
#include <TParticle.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>



#include "AliAnalysisPHOSClusterTask.h"


//AliRoot include files
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliPHOSGeometry.h"
#include "AliEMCALGeometry.h"

#include <TRefArray.h>
#include <iostream>

using namespace std;


#include <climits>		///What s this?

ClassImp(AliAnalysisPHOSClusterTask);

//using namespace AliAnalysisCODEX;

AliAnalysisPHOSClusterTask::AliAnalysisPHOSClusterTask()
  :AliAnalysisTaskSE(),
  fOutput(0),
  fTEST(0x0)
{
  cout << "Hello World from Constructor 1" << endl;
}


AliAnalysisPHOSClusterTask::AliAnalysisPHOSClusterTask(const char* name)
  :AliAnalysisTaskSE(name),
  fOutput(0),
  fTEST(0x0)
{
  cout << "Hello World from constructor 2" << endl;
  DefineOutput(1, TList::Class());
}


AliAnalysisPHOSClusterTask::~AliAnalysisPHOSClusterTask() {
//  if (mOutput){
//    delete mOutput;
//    mOutput = 0x0;
//  }
}

//_______________________________________ UserCreateOutputObjects _______________________________________

void AliAnalysisPHOSClusterTask::UserCreateOutputObjects() {
  cout << "Hello World from UserCreateOutputObjects" << endl;
// cout << "AliAnalysisPHOSCluster: Input settings" << endl;
// cout << "Z Vertex Cut: " << fZVertex << endl;


// Creation of the histograms, this is called once
//

  if (!fOutput) fOutput = new TList();
  fOutput->SetOwner();


  fTEST = new TH1F("fTEST", "TEST HISTOGRAMM", 1, 0.5, 1.5);
  fOutput ->Add(fTEST);


  PostData(1,fOutput);
}

//_______________________________________ UserExec _______________________________________

void AliAnalysisPHOSClusterTask::UserExec(Option_t *){
  cout << "Hello World from UserExec" << endl;
  fTEST ->Fill(1);


  PostData(1,fOutput);
}
