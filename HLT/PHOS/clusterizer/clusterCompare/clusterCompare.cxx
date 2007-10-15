
#include "AliHLTPHOSDigitContainerStruct.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSClusterizer.h"
#include "AliPHOSGetter.h"
#include "TH1I.h"
#include "TObjArray.h"
#include "TFile.h"

#define NLOOPS 1

Int_t clusterCompare();

int main()
{
  clusterCompare();
  return 0;
}

Int_t clusterCompare()
{
  
  AliHLTPHOSDigitContainerStruct *digCon = new AliHLTPHOSDigitContainerStruct();  
  AliHLTPHOSRecPointContainerStruct *recCon = new AliHLTPHOSRecPointContainerStruct();
  AliHLTPHOSClusterizer *clusterizer = new AliHLTPHOSClusterizer();
  AliPHOSGetter *getter = AliPHOSGetter::Instance("/home/odjuvsla/Workspace/Simulations/singlePhoton/galice.root");

  TObjArray *emcRecPoints = getter->EmcRecPoints();

  TH1I *hist = new TH1I("hist", "Difference in number of clusters found", 10, -5, 5);

  int nClusters = 0;
  int a = 0;

  clock_t start = 0; 
  clock_t end = 0;
  clock_t total = 0;
  int count = 0;

  clusterizer->SetOfflineMode(getter);
  clusterizer->SetRecPointContainer(recCon);

  while(a < NLOOPS)
    {
      for(int k = 0; k < clusterizer->GetNEvents(); k++)
	{
	  //  printf("Getting event...");
	  clusterizer->GetEvent(k);
	  //	  getter->Event(k, "R");
	  // printf("Done!\n");

	  //printf("Starting HLT clusterizing... ");
	  start = clock();
	  clusterizer->SetDigitContainer(digCon);
	  nClusters = clusterizer->ClusterizeEvent();
	  clusterizer->CalculateCenterOfGravity();
	  end = clock();
	  //printf("Done!\n");
	  total += end - start;
	  if(nClusters == 2)
	    {
	      printf("WTF! 2 cluster!  ---  Event: %d\n", k);
	      printf("Positions:\n");
	      printf("x = %f -- z = %f\n", recCon->fRecPointArray[0].fX, recCon->fRecPointArray[0].fZ);
	      printf("x = %f -- z = %f\n", recCon->fRecPointArray[1].fX, recCon->fRecPointArray[1].fZ);
	      printf("Energies:\n");
	      printf("%f\n", recCon->fRecPointArray[0].fAmp);
	      printf("%f\n", recCon->fRecPointArray[1].fAmp);
	      printf("Multiplicities:\n");
	      printf("%d\n",recCon->fRecPointArray[0].fMultiplicity);
	      printf("%d\n",recCon->fRecPointArray[1].fMultiplicity);
	      printf("Digits:\n");
	      for(int c = 0; c < recCon->fRecPointArray[0].fMultiplicity; c++)
		{
		  printf("x = %d -- z = %d -- energy = %f\n", recCon->fRecPointArray[0].fDigitsList[c].fX, recCon->fRecPointArray[0].fDigitsList[c].fZ,   recCon->fRecPointArray[0].fDigitsList[c].fEnergy );
		}
	      printf("\n");
	      for(int b = 0; b < recCon->fRecPointArray[1].fMultiplicity; b++)
		{
		  printf("x = %d -- z = %d -- energy = %f\n", recCon->fRecPointArray[1].fDigitsList[b].fX, recCon->fRecPointArray[1].fDigitsList[b].fZ,   recCon->fRecPointArray[1].fDigitsList[b].fEnergy );
		}
	      printf("\n");
	      
		     
	    }
	  else
	    printf("Multiplicity: %d\n\n", recCon->fRecPointArray[0].fMultiplicity);
	  count += nClusters;
	  //hist->Fill(emcRecPoints->GetEntriesFast() - recCon->fNRecPoints);
	}
      a++;
    }
  
  
  TFile *outfile = new TFile("histogram.root","recreate");
  //  hist->Write();
  outfile->Close();

  printf("Clusterizing done!\n");
  printf("Total number of clusters: %d\n", count);
  
  printf("Time pr event (us): %f\n", (float)total/((float)clusterizer->GetNEvents()*(float)a));
  
  return 0;
}
  
  
