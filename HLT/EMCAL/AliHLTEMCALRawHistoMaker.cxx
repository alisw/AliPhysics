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
 * @brief  Histogram maker/pusher for EMCAL HLT  
 */
  

#include "AliHLTEMCALRawHistoMaker.h"
#include "AliHLTEMCALConstants.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloSharedMemoryInterfacev2.h"
//#include "AliHLTCaloRawAnalyzer.h"
#include "AliCaloRawAnalyzer.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"

ClassImp(AliHLTEMCALRawHistoMaker);

AliHLTEMCALRawHistoMaker::AliHLTEMCALRawHistoMaker():
  AliHLTCaloConstantsHandler("EMCAL"),
  fShmPtr(0),
  fMapperPtr(0),
  fRawCounterMemoryPtr(0),
  fAltroRawStreamPtr(0),
  fRawStreamPtr(0),
  fSTURawStreamPtr(0),
  fAnalyzerPtr(0),
  fEMCALConstants(NULL),
  hList(0),
  fChannelEMap(0), fChannelTMap(0), fChannelETMap(0)
{
  // See header file for documentation

  fShmPtr = new AliHLTCaloSharedMemoryInterfacev2("EMCAL");

  fRawCounterMemoryPtr = new AliRawReaderMemory();

  fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawCounterMemoryPtr);

  fRawStreamPtr = new AliCaloRawStreamV3(fRawCounterMemoryPtr, "EMCAL");

  fSTURawStreamPtr = new AliEMCALTriggerSTURawStream(fRawCounterMemoryPtr);

  fEMCALConstants = new AliHLTEMCALConstants();


  // Booking sample histograms

  char id[100];
  char title[100];
 
  fChannelEMap = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fChannelTMap = new TProfile2D *[fCaloConstants->GetNMODULES()];
  fChannelETMap = new TH2F *[fCaloConstants->GetNMODULES()];
  
  hList = new TObjArray;
  
  for (int i=0; i<fCaloConstants->GetNMODULES(); i++) {
    sprintf(title, "E(X vs Z): SM %d ", i);
    sprintf(id, "fChannelEMap%d", i);
    fChannelEMap[i] = new TProfile2D(id,title,48 ,0, 47, 24, 0, 23);
    
    hList->Add(fChannelEMap[i]);
    
    sprintf(title, "T(X vs Z): SM %d ", i);
    sprintf(id, "fChannelTMap%d", i);
    fChannelTMap[i] = new TProfile2D(id,title,48 ,0, 47, 24, 0, 23);

    hList->Add(fChannelTMap[i]);

    printf(title, "(E vs T): SM %d ", i);
    sprintf(id, "fChannelETMap%d", i);
    fChannelETMap[i] = new TH2F(id,title,100 ,0, 50, 100, 0, 500);
    
    hList->Add(fChannelETMap[i]);
  }
  
  
}

AliHLTEMCALRawHistoMaker::~AliHLTEMCALRawHistoMaker() 
{
  //See header file for documentation
}

// Pointer for histograms objects
TObjArray* AliHLTEMCALRawHistoMaker::GetHistograms()
{
  return hList;
}


Int_t
AliHLTEMCALRawHistoMaker::MakeHisto(AliHLTCaloChannelDataHeaderStruct* channelDataHeader,
		const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr,
		const AliHLTUInt32_t size)
{
	int tmpsize=  0;
	Int_t crazyness          = 0;
	Int_t nSamples           = 0;
	Short_t channelCount     = 0;

	AliHLTCaloCoordinate coord;

	AliHLTCaloChannelDataStruct* currentchannel = 0;

	fShmPtr->SetMemory(channelDataHeader);

	AliHLTCaloChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*>(outputPtr);
	AliHLTCaloChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*>(outputPtr+sizeof(AliHLTCaloChannelDataHeaderStruct));

	fRawCounterMemoryPtr->SetMemory(         reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
	fRawCounterMemoryPtr->SetEquipmentID(    fMapperPtr->GetDDLFromSpec(  iter->fSpecification) + fCaloConstants->GetDDLOFFSET() );
	fRawCounterMemoryPtr->Reset();
	fRawCounterMemoryPtr->NextEvent();

//------
	//fRawDataWriter->NewEvent( );

	if(fAltroRawStreamPtr->NextDDL())
	  {
	    int cnt = 0;
	    int fOffset = 0;
	    while( fAltroRawStreamPtr->NextChannel()  )
	      {
		//	 cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T3"  << endl;
		if(  fAltroRawStreamPtr->GetHWAddress() < 128 || ( fAltroRawStreamPtr->GetHWAddress() ^ 0x800) < 128 )
		  {
		    continue;
		  }
		else
		  {
		    ++ cnt;
		    UShort_t* firstBunchPtr = 0;
		    int chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress());
		    //	    HLTError("Channel HW address: %d", fAltroRawStreamPtr->GetHWAddress());

		    //fRawDataWriter->SetChannelId( chId );

		   cout << fRawStreamPtr->GetCaloFlag() << endl;

		   fRawCounterMemoryPtr->Select("EMCAL",44);
		   fRawCounterMemoryPtr->Reset();

		   int nJ, nG;


		   fSTURawStreamPtr->ReadPayLoad();
		   nJ = fSTURawStreamPtr->GetNL1JetPatch();
		   nG = fSTURawStreamPtr->GetNL1GammaPatch();
		   if(nJ != 0 || nG != 0) {

			   Printf(" %d L1Jet, %d L1Gamma", nJ, nG);
			   Printf("Jet thr : %d,  Gamma thr : %d", fSTURawStreamPtr->GetL1JetThreshold(), fSTURawStreamPtr->GetL1GammaThreshold());
			   Printf("####");
		   }



		   fSTURawStreamPtr->DumpPayLoad("L0");

		    //	    return 1;
		    vector <AliCaloBunchInfo> bvctr;
		    while( fAltroRawStreamPtr->NextBunch() == true )
		      {
			//bvctr.push_back( AliCaloBunchInfo( fAltroRawStreamPtr->GetStartTimeBin(), fAltroRawStreamPtr->GetBunchLength(), fAltroRawStreamPtr->GetSignals() ) );

			nSamples = fAltroRawStreamPtr->GetBunchLength();


			   // fRawDataWriter->WriteBunchData( fAltroRawStreamPtr->GetSignals(), nSamples,  fAltroRawStreamPtr->GetEndTimeBin()  );
			firstBunchPtr = const_cast< UShort_t* >(  fAltroRawStreamPtr->GetSignals()  );
		      }

		    //return 1;


		    //    fAnalyzerPtr->SetData( firstBunchPtr, nSamples);
		  AliCaloFitResults res = fAnalyzerPtr->Evaluate( bvctr,  fAltroRawStreamPtr->GetAltroCFG1(), fAltroRawStreamPtr->GetAltroCFG2() );

		    HLTDebug("Channel energy: %f, max sig: %d, gain = %d, x = %d, z = %d", res.GetAmp(), res.GetMaxSig(), (chId >> 12)&0x1, chId&0x3f, (chId >> 6)&0x3f);

		    //	      if(fAnalyzerPtr->GetTiming() > fMinPeakPosition && fAnalyzerPtr->GetTiming() < fMaxPeakPosition)
		    {
		      channelDataPtr->fChannelID =  chId;
		      channelDataPtr->fEnergy = static_cast<Float_t>( res.GetAmp()  ) - fOffset;

		      channelDataPtr->fTime = static_cast<Float_t>(  res.GetTof() );
		      channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
		      channelCount++;
		      channelDataPtr++; // Updating position of the free output.
		    }
		  }


		  }

		   //fRawDataWriter->NewChannel();

	  }


//-------
	currentchannel = fShmPtr->NextChannel();
 
  while(currentchannel != 0) {

    AliHLTCaloChannelRawDataStruct rawdata = fShmPtr->GetRawData(); 

      fMapperPtr->ChannelId2Coordinate(currentchannel->fChannelID, coord);
      
     	fChannelTMap[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fTime);
     	fChannelEMap[coord.fModuleId]->Fill( coord.fZ,  coord.fX , currentchannel->fEnergy);
     	fChannelETMap[coord.fModuleId]->Fill(currentchannel->fEnergy, currentchannel->fTime);
      
	currentchannel = fShmPtr->NextChannel(); // Get the next channel
     
    }

  
return (0); 
}
