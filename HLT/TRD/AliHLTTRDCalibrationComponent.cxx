// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors:                                                               *
 *          for The ALICE HLT Project.                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//  @file   AliHLTTRDCalibrationComponent.cxx
//  @author 
//  @date
//  @brief  A TRDCalibration processing component for the HLT. 
// 

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH2I.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "AliHLTReadoutList.h"

#include "AliHLTTRDCalibrationComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliRawReaderMemory.h"

#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDtrackV1.h"

#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"

#include <cstdlib>
#include <cerrno>
#include <string>

ClassImp(AliHLTTRDCalibrationComponent);

AliHLTTRDCalibrationComponent::AliHLTTRDCalibrationComponent()
  : AliHLTCalibrationProcessor(),
    fTRDCalibraFillHisto(NULL),
    fOutputSize(500000),
    fTracksArray(NULL),
    fOutArray(NULL),
    fAfterRunArray(NULL),
    fDisplayArray(NULL),
    fSavedTimeBins(kFALSE),
    fTrgStrings(NULL),
    fAccRejTrg(0),
    fMinClusters(0),
    fMinTracklets(0),
    fTakeAllEvents(kFALSE)
{
  // Default constructor
}

AliHLTTRDCalibrationComponent::~AliHLTTRDCalibrationComponent()
{
  // Destructor
}

const char* AliHLTTRDCalibrationComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDCalibration"; // The ID of this component
}

void AliHLTTRDCalibrationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(AliHLTTRDDefinitions::fgkTracksDataType);
}

AliHLTComponentDataType AliHLTTRDCalibrationComponent::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
  //  return AliHLTTRDDefinitions::fgkCalibrationDataType;
 
}

int AliHLTTRDCalibrationComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data type
  tgtList.clear();
  tgtList.push_back(AliHLTTRDDefinitions::fgkCalibrationDataType);
  tgtList.push_back(AliHLTTRDDefinitions::fgkEORCalibrationDataType);
  return tgtList.size();
}

void AliHLTTRDCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDCalibrationComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDCalibrationComponent;
};

Int_t AliHLTTRDCalibrationComponent::ScanArgument( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  int i = 0;
  char* cpErr;
  if(!fTrgStrings)
    fTrgStrings = new TObjArray();

  while ( i < argc )
    {
      HLTDebug("argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_size" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing output_size parameter");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
          fOutputSize = strtoul( argv[i+1], &cpErr, 0 );
          if ( *cpErr )
            {
              HLTError("Cannot convert output_size parameter '%s'", argv[i+1] );
              return EINVAL;
            }
          HLTInfo("Output size set to %lu %%", fOutputSize );
          i += 2;
          continue;
        }
      if ( !strcmp( argv[i], "-minClusters" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing minClusters parameter");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
	  fMinClusters = strtoul( argv[i+1], &cpErr, 0 );
          i += 2;
          continue;
        }
      if ( !strcmp( argv[i], "-minTracklets" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing minTracklets parameter");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
	  fMinTracklets = strtoul( argv[i+1], &cpErr, 0 );
          i += 2;
          continue;
        }
      if ( !strcmp( argv[i], "-TrgStr" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing parameter for mbTriggerString");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
          fTrgStrings->Add(new TObjString(argv[i+1]));
          i += 2;
          continue;
        }

      if ( !strcmp( argv[i], "-acceptTrgStr" ) )
        {
          fAccRejTrg=1;
          i += 1;
          continue;
        }
      if ( !strcmp( argv[i], "-rejectTrgStr" ) )
        {
          fAccRejTrg=-1;
          i += 1;
          continue;
        }
      if ( !strcmp( argv[i], "-takeAllEvents" ) )
        {
	  fTakeAllEvents = kTRUE;
	  fAccRejTrg = 0;
          i += 1;
          continue;
        }

      else {
        HLTError("Unknown option '%s'", argv[i] );
        return EINVAL;
      }
    }
  return i;
}

Int_t AliHLTTRDCalibrationComponent::InitCalibration()
{

  if(!fTrgStrings)
    fTrgStrings = new TObjArray();

  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTError("DefaultStorage is not set in CDBManager");
    return -EINVAL;
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTError("Run Number is not set in CDBManager");
    return -EINVAL;
  }
  HLTInfo("CDB default storage: %s; RunNo: %i", (AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder()).Data(), AliCDBManager::Instance()->GetRun());

  if(fTrgStrings->GetEntriesFast()>0 && !fAccRejTrg){
    HLTError("Trigger string(s) given, but acceptTrgStr or rejectTrgStr not selected");
    return -EINVAL;
  }

  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetIsHLT(kTRUE);
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->SetIsHLT(); // per detector
  //fTRDCalibraFillHisto->SetDebugLevel(1);// debug
  fTRDCalibraFillHisto->SetFillWithZero(kTRUE);
  fTRDCalibraFillHisto->SetLinearFitterOn(kTRUE);
  fTRDCalibraFillHisto->SetNumberBinCharge(100);
  
  fTracksArray = new TClonesArray("AliTRDtrackV1");
  fOutArray = new TObjArray(4);
  fAfterRunArray=new TObjArray(5);
  fDisplayArray=new TObjArray(4);

  HLTDebug("run SetupCTPData");
  SetupCTPData();

  return 0;
}

Int_t AliHLTTRDCalibrationComponent::DeinitCalibration()
{
  
  // Deinitialization of the component
  
  HLTDebug("DeinitCalibration");
  delete fTracksArray; fTracksArray=0;
  //fTRDCalibraFillHisto->Destroy();
  //fOutArray->Delete();
  delete fOutArray; fOutArray=0;
  fAfterRunArray->Delete();
  delete fAfterRunArray; fAfterRunArray=0;
  fDisplayArray->Delete();
  delete fDisplayArray; fDisplayArray=0;
  fTrgStrings->Delete();
  delete fTrgStrings; fTrgStrings=0;
  return 0;
}

Int_t AliHLTTRDCalibrationComponent::ProcessCalibration(const AliHLTComponent_EventData& /*evtData*/,
                                                        const AliHLTComponent_BlockData* /*blocks*/,
                                                        AliHLTComponent_TriggerData& /*trigData*/,
                                                        AliHLTUInt8_t* /*outputPtr*/,
                                                        AliHLTUInt32_t& /*size*/,
                                                        vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
  // Process an event
	
  TClonesArray* TCAarray[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t usedEntries = 0;
  Int_t blockOrObject = 0;
  Int_t nTimeBins = -1;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTRDDefinitions::fgkTracksDataType); pBlock; pBlock=GetNextInputBlock()) 
    {
      TCAarray[0] = fTracksArray;
      AliHLTTRDUtils::ReadTracks(TCAarray[0], pBlock->fPtr, pBlock->fSize, &nTimeBins);
      usedEntries = 1;
      blockOrObject = -1;
    }  

  for(const TObject *iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkHiLvlTracksDataType); iter; iter = GetNextInputObject()) 
    {
      if(blockOrObject<0){
	HLTError("You may not mix high level and low level!");
	return -1;
      }

      TCAarray[usedEntries] = dynamic_cast<TClonesArray*>(const_cast<TObject*>(iter));
      if(!TCAarray[usedEntries])continue;
      TObjString* strg = dynamic_cast<TObjString*>(const_cast<TObject*>(GetNextInputObject()));
      if(!strg)continue;

      nTimeBins = strg->String().Atoi();
      usedEntries++;
      blockOrObject = 1;
    }

  if(!blockOrObject)
    return 0;

  if(!fSavedTimeBins){
    if(nTimeBins<0){
      HLTFatal("Number of timebins is negative!");
      return -1;
    }
    HLTDebug("Saving number of time bins which was read from input block. Value is: %d", nTimeBins);
    fTRDCalibraFillHisto->Init2Dhistos(nTimeBins); // initialise the histos
    fTRDCalibraFillHisto->SetNumberClusters(fMinClusters); // At least fMinClusters clusters
    fTRDCalibraFillHisto->SetNumberClustersf(nTimeBins); // Not more than %d  clusters
    fSavedTimeBins=kTRUE;
  }

  Bool_t TriggerPassed = fTakeAllEvents;

  if(fAccRejTrg){
    if(fAccRejTrg>0){
      TriggerPassed=kFALSE;
      for(int i = 0; i < fTrgStrings->GetEntriesFast(); i++){
	const TObjString *const obString=(TObjString*)fTrgStrings->At(i);
	const TString tString=obString->GetString();
	if(CheckCTPTrigger(tString.Data())>0){TriggerPassed=kTRUE; break;}
      }
    }
    else{
      TriggerPassed=kTRUE;
      for(int i = 0; i < fTrgStrings->GetEntriesFast(); i++){
	const TObjString *const obString=(TObjString*)fTrgStrings->At(i);
	const TString tString=obString->GetString();
	if(CheckCTPTrigger(tString.Data())>0){TriggerPassed=kFALSE; break;}
      }
    }
  }
  
  fTRDCalibraFillHisto->SetCH2dOn(TriggerPassed);
  fTRDCalibraFillHisto->SetPH2dOn(TriggerPassed);
  for(int i=0; i<usedEntries; i++){
    const TClonesArray* inArr = TCAarray[i];
    Int_t nbEntries = inArr->GetEntries();
    HLTDebug(" %i TRDtracks in tracksArray", nbEntries);
    AliTRDtrackV1* trdTrack = 0x0;
    for (Int_t ii = 0; ii < nbEntries; ii++){
      HLTDebug("%i/%i: ", ii+1, nbEntries);
      trdTrack = (AliTRDtrackV1*)inArr->At(ii);
      if(trdTrack->GetNumberOfTracklets()<fMinTracklets)continue;
      fTRDCalibraFillHisto->UpdateHistogramsV1(trdTrack);
      // for(int i3=0; i3<7; i3++)
      //   if(trdTrack->GetTracklet(i3))trdTrack->GetTracklet(i3)->Bootstrap(fReconstructor);
    }
  }

  if(!fOutArray->At(0))FormOutput(0);
  if(!fDisplayArray->At(0))FormOutput(1);
  PushBack(fDisplayArray, AliHLTTRDDefinitions::fgkCalibrationDataType);

  if(blockOrObject<0){
    TCAarray[0]->Delete();
  }

  return 0;

}

/**
 * Form output array of histrograms
 */
//============================================================================
void AliHLTTRDCalibrationComponent::FormOutput(Int_t param)
{
  // gain histo
  TH2I *hCH2d = fTRDCalibraFillHisto->GetCH2d();
  if(!param)fOutArray->Add(hCH2d);
  else fDisplayArray->Add(hCH2d);

  // drift velocity histo
  TProfile2D *hPH2d = fTRDCalibraFillHisto->GetPH2d();
  if(!param)fOutArray->Add(hPH2d);
  else fDisplayArray->Add(hPH2d);

  // PRF histo
  TProfile2D *hPRF2d = fTRDCalibraFillHisto->GetPRF2d();
  if(!param)fOutArray->Add(hPRF2d);
  else fDisplayArray->Add(hPRF2d);

  // Vdrift Linear Fit
  if(!param){
    AliTRDCalibraVdriftLinearFit *hVdriftLinearFitOne=(AliTRDCalibraVdriftLinearFit *)fTRDCalibraFillHisto->GetVdriftLinearFit();
    fOutArray->Add(hVdriftLinearFitOne);
  }
  else{
    TH2S *hVdriftLinearFitOne = (TH2S *)(((AliTRDCalibraVdriftLinearFit *)fTRDCalibraFillHisto->GetVdriftLinearFit())->GetLinearFitterHisto(10,kTRUE));	
    fDisplayArray->Add(hVdriftLinearFitOne);
  }

  HLTDebug("GetCH2d = 0x%x; NEntries = %i; size = %i", hCH2d, hCH2d->GetEntries(), sizeof(*hCH2d));
  hCH2d->Print();
  HLTDebug("GetPH2d = 0x%x; NEntries = %i; size = %i", hPH2d, hPH2d->GetEntries(), sizeof(*hPH2d));
  hPH2d->Print();
  HLTDebug("GetPRF2d = 0x%x; NEntries = %i; size = %i", hPRF2d, hPRF2d->GetEntries(), sizeof(*hPRF2d));
  hPRF2d->Print();
  //HLTDebug("GetVdriftLinearFit = 0x%x; size = %i", hVdriftLinearFitOne, sizeof(hVdriftLinearFitOne)); 
  
  HLTDebug("output Array: pointer = 0x%x; NEntries = %i; size = %i", fOutArray, fOutArray->GetEntries(), sizeof(fOutArray));
   
}

Int_t AliHLTTRDCalibrationComponent::ShipDataToFXS(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  //fTRDCalibraFillHisto->DestroyDebugStreamer();

  AliHLTReadoutList rdList(AliHLTReadoutList::kTRD);

  EORCalibration();
  
  fOutArray->Remove(fOutArray->FindObject("AliTRDCalibraVdriftLinearFit"));
  //fOutArray->Remove(fOutArray->FindObject("PRF2d"));
  //fOutArray->Remove(fOutArray->FindObject("PH2d"));
  //fOutArray->Remove(fOutArray->FindObject("CH2d"));

  if(!(fOutArray->FindObject("CH2d"))) {
    TH2I * ch2d = new TH2I("CH2d","Nz0Nrphi0",100,0.0,300.0,540,0,540);
    fOutArray->Add(ch2d);
  }

  if(!(fOutArray->FindObject("PH2d"))) {
    TProfile2D * ph2d = new TProfile2D("PH2d","Nz0Nrphi0",30,-0.05,2.95,540,0,540);
    fOutArray->Add(ph2d);
  }

  if(!(fOutArray->FindObject("PRF2d"))) {
    TProfile2D * prf2d = new TProfile2D("PRF2d","Nz0Nrphi0Ngp3",60,-9.0,9.0,540,0,540);
    fOutArray->Add(prf2d);
  }

  HLTDebug("Size of the fOutArray is %d\n",fOutArray->GetEntriesFast());

  PushToFXS((TObject*)fOutArray, "TRD", "GAINDRIFTPRF", &rdList );
  //PushToFXS((TObject*)fOutArray->FindObject("CH2d"), "TRD", "GAINDRIFTPRF", rdList.Buffer() );

  return 0;
}

Int_t AliHLTTRDCalibrationComponent::EORCalibration()
{
  //Also Fill histograms for the online display
  TH2I *hCH2d=(TH2I*)fOutArray->FindObject("CH2d");
  TProfile2D *hPH2d=(TProfile2D*)fOutArray->FindObject("PH2d");
  TProfile2D *hPRF2d= (TProfile2D*)fOutArray->FindObject("PRF2d");
  AliTRDCalibraVdriftLinearFit* hVdriftLinearFit = (AliTRDCalibraVdriftLinearFit*)fOutArray->FindObject("AliTRDCalibraVdriftLinearFit");
 

  if(!hCH2d || !hPH2d || !hPRF2d || !hVdriftLinearFit) return 0; 

  //Fit
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();

  //Gain
  calibra->SetMinEntries(100);
  calibra->AnalyseCH(hCH2d);
  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();
  TH1F *coefgain = 0x0;
  // enough statistics
  //if ((nbtg >                  0) && 
  //   (nbfit        >= 0.2*nbE)) {
  // create the cal objects
  //calibra->PutMeanValueOtherVectorFit(1,kTRUE);
  TObjArray object           = calibra->GetVectorFit();
  AliTRDCalDet *objgaindet   = calibra->CreateDetObjectGain(&object,kFALSE);
  coefgain                   = objgaindet->MakeHisto1DAsFunctionOfDet();
  //}
  calibra->ResetVectorFit();

  // vdrift second method
  calibra->SetMinEntries(100); // If there is less than 100
  hVdriftLinearFit->FillPEArray();
  calibra->AnalyseLinearFitters(hVdriftLinearFit);
  nbtg = 540;
  nbfit = calibra->GetNumberFit();
  nbE   = calibra->GetNumberEnt();
  TH1F *coefdriftsecond = 0x0;
  // enough statistics
  //if ((nbtg >                  0) && 
  // (nbfit        >= 0.1*nbE)) {
  // create the cal objects
  //calibra->PutMeanValueOtherVectorFit(1,kTRUE);
  object  = calibra->GetVectorFit();
  AliTRDCalDet *objdriftvelocitydetsecond = calibra->CreateDetObjectVdrift(&object,kTRUE);
  objdriftvelocitydetsecond->SetTitle("secondmethodvdrift");
  coefdriftsecond  = objdriftvelocitydetsecond->MakeHisto1DAsFunctionOfDet();
  //}
  calibra->ResetVectorFit();
  
  // vdrift first method
  calibra->SetMinEntries(100*20); // If there is less than 20000
  calibra->AnalysePH(hPH2d);
  nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
  nbfit        = calibra->GetNumberFit();
  nbE          = calibra->GetNumberEnt();
  TH1F *coefdrift = 0x0;
  TH1F *coeft0 = 0x0;
  // enough statistics
  //if ((nbtg >                  0) && 
  // (nbfit        >= 0.2*nbE)) {
  // create the cal objects
  //calibra->PutMeanValueOtherVectorFit(1,kTRUE);
  //calibra->PutMeanValueOtherVectorFit2(1,kTRUE);
  object  = calibra->GetVectorFit();
  AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
  coefdrift        = objdriftvelocitydet->MakeHisto1DAsFunctionOfDet();
  object              = calibra->GetVectorFit2();
  AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
  coeft0        = objtime0det->MakeHisto1DAsFunctionOfDet();
  //}
  calibra->ResetVectorFit();
           

  //PRF
  calibra->SetMinEntries(200); 
  calibra->AnalysePRFMarianFit(hPRF2d);
  nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
  nbfit        = calibra->GetNumberFit();
  nbE          = calibra->GetNumberEnt();
  TH1F *coefprf = 0x0;
  // enough statistics
  //if ((nbtg >                  0) && 
  //  (nbfit        >= 0.95*nbE)) {
  // create cal pad objects 
  object            = calibra->GetVectorFit();
  TObject *objPRFpad          = calibra->CreatePadObjectPRF(&object);
  coefprf                     = ((AliTRDCalPad *) objPRFpad)->MakeHisto1D();
  //}
  calibra->ResetVectorFit();


  coefgain->SetName("coefgain");
  coefprf->SetName("coefprf");
  coefdrift->SetName("coefdrift");
  coefdriftsecond->SetName("coefdriftsecond");
  coeft0->SetName("coeft0");
  if(coefgain) fAfterRunArray->Add(coefgain);
  if(coefprf) fAfterRunArray->Add(coefprf);
  if(coefdrift) fAfterRunArray->Add(coefdrift);
  if(coefdriftsecond) fAfterRunArray->Add(coefdriftsecond);
  if(coeft0) fAfterRunArray->Add(coeft0);
  

  if(coefgain||coefprf||coefdrift||coeft0||coefdriftsecond) {
    PushBack(fAfterRunArray, AliHLTTRDDefinitions::fgkEORCalibrationDataType);
  }

  // TString fileName="/tmp/CalibHistoDump_run";
  // fileName+=AliCDBManager::Instance()->GetRun();
  // fileName+=".root";
  // HLTInfo("Dumping Histogram file to %s",fileName.Data());
  // TFile* file = TFile::Open(fileName, "RECREATE");
  // fAfterRunArray->Write();
  // fOutArray->Write();
  // file->Close();
  // HLTInfo("Histogram file dumped");

  return 0;
}	

