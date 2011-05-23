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
  fChannelEMap(0), fChannelTMap(0), fChannelETMap(0), h2DTRU(0), h2DSTU(0),
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

  fChannelEMap = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fChannelTMap = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fChannelETMap = new TH2F *[fCaloConstants->GetNMODULES()];
  
  h2DTRU = new TH2I ("h2DTRU","",24,0,24,4,0,4);
  hList->Add(h2DTRU);

  h2DSTU = new TH2I ("h2DSTU","",24,0,24,4,0,4);
  hList->Add(h2DSTU);

  
  for (int i=0; i<fCaloConstants->GetNMODULES(); i++) {
    sprintf(title, "E(X vs Z): SM %d ", i);
    sprintf(id, "fChannelEMap%d", i);
    fChannelEMap[i] = new TProfile2D(id,title,48 ,0, 47, 24, 0, 23);
    
    hList->Add(fChannelEMap[i]);
    
    sprintf(title, "T(X vs Z): SM %d ", i);
    sprintf(id, "fChannelTMap%d", i);
    fChannelTMap[i] = new TProfile2D(id,title,48 ,0, 47, 24, 0, 23);

    hList->Add(fChannelTMap[i]);

    sprintf(title, "(E vs T): SM %d ", i);
    sprintf(id, "fChannelETMap%d", i);
    fChannelETMap[i] = new TH2F(id,title,100 ,0, 50, 100, 0, 500);
    
    hList->Add(fChannelETMap[i]);
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
	  
	  // NULL pointer
	  cout << "FROM HISTOMAKER NO CLUSTER POINTER: " << caloClusterHeaderPtr << endl;
	  
	} else {
	  
	  // stuff to handle clusters here
	  fClusterReaderPtr->SetMemory(caloClusterHeaderPtr);

	  cout << "FROM HISTOMAKER WE HAVE THE CLUSTER POINTER: " << caloClusterHeaderPtr << endl;


	  while((caloClusterStructPtr = fClusterReaderPtr->NextCluster()) != 0)
	    { 
	       
	      cout << "cluster type: " << caloClusterStructPtr->fClusterType << endl;
	       
	      cout << " COORDINATES FROM HISTOMAKER: " << 
		" fX:" << caloClusterStructPtr->fGlobalPos[0] <<
		" fY:" << caloClusterStructPtr->fGlobalPos[1] <<
		" fZ:" << caloClusterStructPtr->fGlobalPos[2] << endl;
 	     

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

	fShmPtr->SetMemory(channelDataHeader);

	currentchannel = fShmPtr->NextChannel();

	while(currentchannel != 0) {
	  
	  fMapperPtr->ChannelId2Coordinate(currentchannel->fChannelID, coord);

	  cout << " from histo maker ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << currentchannel->fEnergy << endl;
	  cout << " fX: " << coord.fX << " fZ: " << coord.fZ << endl;

	  fChannelTMap[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fTime);
	  fChannelEMap[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fEnergy);
	  fChannelETMap[coord.fModuleId]->Fill(currentchannel->fEnergy, currentchannel->fTime);
	  
	  currentchannel = fShmPtr->NextChannel(); // Get the next channel
	  
	}
	
	

return (0); 
}
