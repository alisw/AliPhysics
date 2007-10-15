
#include "AliHLTPHOSDigitContainerStruct.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGetter.h"
#include "AliHLTPHOSClusterizer.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSGeometry.h"
#include "TH2I.h"
#include "TFile.h"

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

  AliPHOSDigit * offDig = 0;
  int nClusters = 0;
  float energy = 0;
  float time = 0;
  Int_t coord[4];
  int a = 0;

  clock_t start = 0; 
  clock_t end = 0;
  clock_t total = 0;

  AliPHOSGeometry *geometry = AliPHOSGeometry::GetInstance();

  TH2F *hist = new TH2F("hist", "hist", 56, 1, 56, 64, 2, 65);
  TFile *outfile = new TFile("histogram.root","recreate");

  AliPHOSGetter *getter = AliPHOSGetter::Instance("/home/odjuvsla/Workspace/Simulations/test4/galice.root");
  TClonesArray* digs = getter->Digits(); 
  
  clusterizer->SetOfflineMode("/home/odjuvsla/Workspace/Simulations/test4/galice.root");
  clusterizer->SetRecPointContainer(recCon);

  while(a < 5000)
    {
      for(int k = 0; k < getter->MaxEvent(); k++)
	{
	  
	  //getter->Event(k,"D");
	  //digCon->fNDigits = digs->GetEntries();
	  
	  clusterizer->GetEvent(k);
	  
	  /*  for(int i = 0; i < digs->GetEntries(); i++)
	    {
	      offDig = (AliPHOSDigit*)digs->At(i);
	      digCon->fDigitDataStruct[i].fID = offDig->GetId();
	      digCon->fDigitDataStruct[i].fEnergy = offDig->GetEnergy();
	      digCon->fDigitDataStruct[i].fTime = offDig->GetTime();
	      hist->Fill(coord[3], coord[2], offDig->GetEnergy());
	    }
	    
	  
  	  outfile->Open("histogram.root", "recreate");
	  hist->Write();
	  outfile->Close("histogram.root");
	  
	  */
	  
	  start = clock();
	  clusterizer->SetDigitContainer(digCon);
	  nClusters = clusterizer->ClusterizeEvent();
	  clusterizer->CalculateCenterOfGravity();
	  end = clock();

	  total += end - start;
	  /*
	    printf("Event #: %d\n", k);
	    printf("Number of clusters found: ", nClusters);

	    for(int m = 0; m < nClusters; m++)
	    {
	    energy = 0;
	    for (int n = 0; n < recCon->fRecPointArray[m].fMultiplicity; n++)
	    {
	    energy += recCon->fRecPointArray[m].fDigitsList[n].fEnergy;
	    }
	    printf("Clusters:\n");
	    printf("Cluster energy: %f\n", energy);
	    printf("Cluster position: x = %f -- z = %f\n", recCon->fRecPointArray[m].fX, recCon->fRecPointArray[m].fZ);
	    printf("Digit multiplicity: %d\n", recCon->fRecPointArray[m].fMultiplicity);
	    }
	  */     
	}
      a++;
    }
  printf("Clusterizing done!\n");
  
  printf("Time pr event (us): %f\n", (float)total/((float)getter->MaxEvent()*(float)a));
  
  return 0;
}
  
  
