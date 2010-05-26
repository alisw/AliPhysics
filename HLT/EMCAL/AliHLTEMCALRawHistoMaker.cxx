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
  fChannelEMap(0), fChannelTMap(0), fChannelETMap(0), h2DTRU(0), h2DSTU(0)

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

	fRawCounterMemoryPtr->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
	fRawCounterMemoryPtr->SetEquipmentID(    fMapperPtr->GetDDLFromSpec(  iter->fSpecification) + fCaloConstants->GetDDLOFFSET() );

	fRawCounterMemoryPtr->Reset();

	//fRawCounterMemoryPtr->NextEvent();

	//--- from STU macro
   	//while ( rawReader->NextEvent() )

	while (fRawCounterMemoryPtr->NextEvent())
	{
		fRawCounterMemoryPtr->Reset();
		fRawCounterMemoryPtr->Select("EMCAL",44);

		const UInt_t* evtId = fRawCounterMemoryPtr->GetEventId();
		int evno_raw = (int)evtId[0];

		UInt_t eventType = fRawCounterMemoryPtr->GetType();

		//cout << "event type: " << eventType << endl;

//		if (eventType != AliRawEventHeaderBase::kPhysicsEvent) continue;

//		fRawCounterMemoryPtr->DumpData();

//		fRawCounterMemoryPtr->Reset();

		fSTURawStreamPtr->ReadPayLoad();
//		fSTURawStreamPtr->DumpPayLoad("ALL");

		UInt_t adc[96];

		for (Int_t i=0;i<96;i++)
			adc[i] = 0;

		fSTURawStreamPtr->GetADC(29, adc);

		for (Int_t i=0;i<96;i++)
		{
			Int_t x = i / 4;
			Int_t y = i % 4;

			h2DSTU->Fill( x , y , adc[i] );
		}

		fRawCounterMemoryPtr->Reset();
		fRawCounterMemoryPtr->Select("EMCAL",0,43);

		while (fRawStreamPtr->NextDDL())
		{
			Int_t iRCU = fRawStreamPtr->GetDDLNumber() % 2;

			while (fRawStreamPtr->NextChannel())
			{
				Int_t iSM = fRawStreamPtr->GetModule();

				Int_t iTRUId = 3 * iSM + iRCU * fRawStreamPtr->GetBranch() + iRCU;

				if (iSM==0)
				{
					Int_t nsamples = 0;
					vector<AliCaloBunchInfo> bunchlist;
					while (fRawStreamPtr->NextBunch())
					{
						nsamples += fRawStreamPtr->GetBunchLength();
						bunchlist.push_back( AliCaloBunchInfo(fRawStreamPtr->GetStartTimeBin(), fRawStreamPtr->GetBunchLength(), fRawStreamPtr->GetSignals() ) );
					}

					Int_t max = 0;

					for (std::vector<AliCaloBunchInfo>::iterator itVectorData = bunchlist.begin(); itVectorData != bunchlist.end(); itVectorData++)
					{
                        AliCaloBunchInfo bunch = *(itVectorData);

                        const UShort_t* sig = bunch.GetData();
                        Int_t startBin = bunch.GetStartBin();

                        for (Int_t iS = 0; iS < bunch.GetLength(); iS++)
                        {
							if ( sig[iS] > max ) max = sig[iS];
						}
					}

					if (nsamples) // this check is needed for when we have zero-supp. on, but not sparse readout
					{
						if (fRawStreamPtr->IsTRUData() && fRawStreamPtr->GetColumn() < 96)
						{
							if (iTRUId == 2)
							{
								Int_t x = fRawStreamPtr->GetColumn() / 4;
								Int_t y = fRawStreamPtr->GetColumn() % 4;

							h2DTRU->Fill( x , y , max );
							}
						}
					}
				}
			}
		}

//		c_0->cd(1);
//		h2DTRU->Draw("TEXT");
//
//		c_0->cd(2);
//		h2DSTU->Draw("TEXT");
//
//		c_0->Update();
//
//		h2DTRU->Reset();
//		h2DSTU->Reset();
	}

	//--- end of while


//------ FIXME
//------ old stuff to be removed or re-utilized
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
		    //HLTError("Channel HW address: %d", fAltroRawStreamPtr->GetHWAddress());


		   //cout << fRawStreamPtr->GetCaloFlag() << endl;
		    //	    return 1;

		    vector <AliCaloBunchInfo> bvctr;
		    while( fAltroRawStreamPtr->NextBunch() == true )
		      {
			bvctr.push_back( AliCaloBunchInfo( fAltroRawStreamPtr->GetStartTimeBin(), fAltroRawStreamPtr->GetBunchLength(), fAltroRawStreamPtr->GetSignals() ) );

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
