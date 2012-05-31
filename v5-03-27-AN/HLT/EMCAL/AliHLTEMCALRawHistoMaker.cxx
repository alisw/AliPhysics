/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 * INFN Laboratori Nazionali di Frascati                                  *
 * Primary Authors: Federico Ronchetti                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * @file   AliHLTEMCALRawHistoMaker.cxx
 * @author Federico Ronchetti
 * @date 
 * @brief  Online Monitoring Histogram maker for EMCAL  
 */
  

#include "AliHLTEMCALRawHistoMaker.h"
#include "AliHLTEMCALConstants.h"
#include "AliHLTEMCALMapper.h"
//#include "AliHLTCaloChannelDataStruct.h"

#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloSharedMemoryInterfacev2.h"

//#include "AliCaloRawAnalyzer.h"
//#include "AliCaloBunchInfo.h"
//#include "AliCaloFitResults.h"




ClassImp(AliHLTEMCALRawHistoMaker);

AliHLTEMCALRawHistoMaker::AliHLTEMCALRawHistoMaker():
  AliHLTCaloConstantsHandler("EMCAL"),
  fShmPtr(0),
  fMapperPtr(0),
  //  fRawCounterMemoryPtr(0),
  //fAltroRawStreamPtr(0),
  //fRawStreamPtr(0),
  //fSTURawStreamPtr(0),
  //fAnalyzerPtr(0),
  fEMCALConstants(NULL),
  hList(0),
  fAmp(0), fTime(0), fAT(0), fCellVsEne(0),
  fClusterReaderPtr(0)

{
  // See header file for documentation

  fShmPtr = new AliHLTCaloSharedMemoryInterfacev2("EMCAL");

  //fRawCounterMemoryPtr = new AliRawReaderMemory();

  //fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawCounterMemoryPtr);

  //fRawStreamPtr = new AliCaloRawStreamV3(fRawCounterMemoryPtr, "EMCAL");

  //fSTURawStreamPtr = new AliEMCALTriggerSTURawStream(fRawCounterMemoryPtr);

  fEMCALConstants = new AliHLTEMCALConstants();

  fClusterReaderPtr = new AliHLTCaloClusterReader();

  // Booking histograms

  char id[100];
  char title[100];
 
  hList = new TObjArray;

  // channel histograms
  fAmp = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fTime = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fAT = new TH2F *[fCaloConstants->GetNMODULES()];
  
  // cluster histograms
  fCellVsEne = new TH2F  *[fCaloConstants->GetNMODULES()];

  //  fCellVsClus = new TH1F ("fCellVsClus","",100,0,50);
  //hList->Add(fCellVsClus);


  for (int i=0; i<fCaloConstants->GetNMODULES(); i++) {



    //sprintf(title, "Row vs Col : SM %d ", i);
    //sprintf(id, "hAmp%d", i);
    //hAmp[i] = new TProfile2D(id,title, 48, -0.5, 47.5, 24, -0.5, 23.5); 
    //sprintf(id, "hTime%d", i);
    //hTime[i] = new TProfile2D(id,title, 48, -0.5, 47.5, 24, -0.5, 23.5); 


    sprintf(title, "X_Z AMP: SM %d ", i);
    sprintf(id, "fAmp%d", i);
    fAmp[i] = new TProfile2D(id,title, 48, -0.5, 47.5, 24, -0.5, 23.5); 
    
    hList->Add(fAmp[i]);
    
    sprintf(title, "X_Z TIME: SM %d ", i);
    sprintf(id, "fTime%d", i);
    fTime[i] = new TProfile2D(id,title, 48, -0.5, 47.5, 24, -0.5, 23.5); 

    hList->Add(fTime[i]);

    sprintf(title, "AMP_TIME: SM %d ", i);
    sprintf(id, "fAT%d", i);
    fAT[i] = new TH2F(id,title, 1024, -0.5, 1023.5, 200, -0.5, 199.5);
    
    hList->Add(fAT[i]);

    sprintf(title, "Cell_Energy: SM %d ", i);
    sprintf(id, "fCellVsEne%d", i);
    fCellVsEne[i] = new TH2F(id,title,50 ,0, 50, 10, 0, 10);
    
    hList->Add(fCellVsEne[i]);

  }
  
  
}

AliHLTEMCALRawHistoMaker::~AliHLTEMCALRawHistoMaker() 
{
  //See header file for documentation
}

// Pointer to histograms objects
TObjArray* AliHLTEMCALRawHistoMaker::GetHistograms()
{
  return hList;
}


Int_t
AliHLTEMCALRawHistoMaker::MakeHisto(AliHLTCaloChannelDataHeaderStruct* channelDataHeader,  
				    AliHLTCaloClusterHeaderStruct *caloClusterHeaderPtr, 
				    int beverbose)
{
	//int tmpsize =  0;
	Int_t crazyness          = 0;
	Int_t nSamples           = 0;
	Short_t channelCount     = 0;


	// Channel variables
	AliHLTCaloCoordinate coord;
	AliHLTCaloChannelDataStruct* currentchannel = 0;


	// Cluster variables
	// Pointer to Cluster struture
	AliHLTCaloClusterDataStruct* caloClusterStructPtr = 0;
	Int_t nClusters = 0;

	if (!caloClusterHeaderPtr) {
	 	  
	} else {
	  
	  // stuff to handle clusters here

	  fClusterReaderPtr->SetMemory(caloClusterHeaderPtr);
	  
	  while((caloClusterStructPtr = fClusterReaderPtr->NextCluster()) != 0)
	  
	    { 
	       
	      cout << "cluster type: " << caloClusterStructPtr->fClusterType << endl;
	       
	      cout << " COORDINATES FROM HISTOMAKER: " << 
		" fX:" << caloClusterStructPtr->fGlobalPos[0] <<
		" fY:" << caloClusterStructPtr->fGlobalPos[1] <<
		" fZ:" << caloClusterStructPtr->fGlobalPos[2] <<
		" fModule: " << caloClusterStructPtr->fModule << 
		" fCell: " << caloClusterStructPtr->fNCells   <<
		" fEnergy " << caloClusterStructPtr->fEnergy << endl;
 	     
	      fCellVsEne[caloClusterStructPtr->fModule]->Fill(caloClusterStructPtr->fNCells, caloClusterStructPtr->fEnergy);
	  	  

	      UShort_t *idArrayPtr = new UShort_t[caloClusterStructPtr->fNCells];
	      Double32_t *ampFracArrayPtr = new Double32_t[caloClusterStructPtr->fNCells];
      
	      for(UInt_t index = 0; index < caloClusterStructPtr->fNCells; index++)
		{
		  fClusterReaderPtr->GetCell(caloClusterStructPtr, idArrayPtr[index], ampFracArrayPtr[index], index);
		  printf("EM: cellId: %d\n", idArrayPtr[index]);;
		}

	      delete [] idArrayPtr;
	      delete [] ampFracArrayPtr;

	      nClusters++;

	    }
	  
	}



	// begin scan channel data and fill histograms


	if (channelDataHeader) {
	 
	  fShmPtr->SetMemory(channelDataHeader);

	  currentchannel = fShmPtr->NextChannel();

	  while(currentchannel != 0) {
	  
	    fMapperPtr->ChannelId2Coordinate(currentchannel->fChannelID, coord);
	    
	    cout << " from histo maker ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << endl;
	    cout << " fX: " << coord.fX << " fZ: " << coord.fZ << endl;
	    cout << " channel ID: " << currentchannel->fChannelID << endl;
	    cout << " channel AMPLITUDE (called energy): " << currentchannel->fEnergy << endl;
	    cout << " channel time: " << currentchannel->fTime << endl;
	    
	    
	    fTime[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fTime);
	    fAmp[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fEnergy);
	    fAT[coord.fModuleId]->Fill(currentchannel->fEnergy, currentchannel->fTime);
	    
	    currentchannel = fShmPtr->NextChannel(); // Get the next channel
	    
	  }
	}
	

return (0); 
}
