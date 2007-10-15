
#include "AliHLTPHOSDigitContainerStruct.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSClusterizer.h"
#include "AliPHOSGetter.h"

#define NLOOPS 10

Int_t clusterTest();

int main()
{
  clusterTest();
  return 0;
}

Int_t clusterTest()
{
  
  AliHLTPHOSDigitContainerStruct *digCon = new AliHLTPHOSDigitContainerStruct();  
  AliHLTPHOSRecPointContainerStruct *recCon = new AliHLTPHOSRecPointContainerStruct();
  AliHLTPHOSClusterizer *clusterizer = new AliHLTPHOSClusterizer();
  AliPHOSGetter *getter = AliPHOSGetter::Instance("/home/odjuvsla/Workspace/Simulations/test3/galice.root");
  int nClusters = 0;
  int a = 0;

  clock_t start = 0; 
  clock_t end = 0;
  clock_t total = 0;
  
  clusterizer->SetOfflineMode(getter);
  clusterizer->SetRecPointContainer(recCon);

  while(a < NLOOPS)
    {
      for(int k = 0; k < clusterizer->GetNEvents(); k++)
	{
	  	  
	  clusterizer->GetEvent(k);
  
	  start = clock();
	  clusterizer->SetDigitContainer(digCon);
	  nClusters = clusterizer->ClusterizeEvent();
	  clusterizer->CalculateCenterOfGravity();
	  end = clock();

	  total += end - start;
	}
      a++;
    }
  printf("Clusterizing done!\n");
  
  printf("Time pr event (us): %f\n", (float)total/((float)clusterizer->GetNEvents()*(float)a));
  
  return 0;
}
  
  
