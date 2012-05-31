// $Id: AliHLTTRDCalibFitComponent.cxx 40282 2010-04-09 13:29:10Z richterm $

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

//  @file   AliHLTTRDCalibFitComponent.cxx
//  @author Theodor Rascanu
//  @date   25.04.2010
//  @brief  A TRDCalibration fitting component for the HLT. 
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

#include "AliHLTTRDCalibFitComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"

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

ClassImp(AliHLTTRDCalibFitComponent);

AliHLTTRDCalibFitComponent::AliHLTTRDCalibFitComponent()
  : AliHLTCalibrationProcessor(),
    fOutputSize(500000),
    fOutArray(NULL),
    fAfterRunArray(NULL),
    fNoOfSM(0),
    fNoOfIncSM(0)
{
  // Default constructor

  for(int i=0; i<18; i++)
    fIncSM[i]=kFALSE;

}

AliHLTTRDCalibFitComponent::~AliHLTTRDCalibFitComponent()
{
  // Destructor
}

const char* AliHLTTRDCalibFitComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDCalibFit"; // The ID of this component
}

void AliHLTTRDCalibFitComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(AliHLTTRDDefinitions::fgkCalibrationDataType);
}

AliHLTComponentDataType AliHLTTRDCalibFitComponent::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
  //  return AliHLTTRDDefinitions::fgkCalibrationDataType;
 
}

int AliHLTTRDCalibFitComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data type
  tgtList.clear();
  tgtList.push_back(AliHLTTRDDefinitions::fgkCalibrationDataType);
  tgtList.push_back(AliHLTTRDDefinitions::fgkEORCalibrationDataType);
  return tgtList.size();
}

void AliHLTTRDCalibFitComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDCalibFitComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDCalibFitComponent;
};

Int_t AliHLTTRDCalibFitComponent::InitCalibration()
{
  for(int i=0; i<18; i++)
    fIncSM[i]=kFALSE;

  fOutArray = new TObjArray(4);
  fAfterRunArray=new TObjArray(5);

  return 0;
}

Int_t AliHLTTRDCalibFitComponent::DeinitCalibration()
{
  
  // Deinitialization of the component
  
  HLTDebug("DeinitCalibration");
  //fOutArray->Delete();
  delete fOutArray; fOutArray=0;
  fAfterRunArray->Delete();
  delete fAfterRunArray; fAfterRunArray=0;
  return 0;
}

Int_t AliHLTTRDCalibFitComponent::ProcessCalibration(const AliHLTComponent_EventData& /*evtData*/,
                                                        const AliHLTComponent_BlockData* /*blocks*/,
                                                        AliHLTComponent_TriggerData& /*trigData*/,
                                                        AliHLTUInt8_t* /*outputPtr*/,
                                                        AliHLTUInt32_t& /*size*/,
                                                        vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
  // Process an event

  if(!IsDataEvent())return 0;

  int lastSM = -1;

  for(const TObject* iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkCalibrationDataType);
      iter != NULL; iter = GetNextInputObject() ) {

    if(!dynamic_cast<const TObjArray*>(iter))
      continue;

    AliHLTUInt32_t spec = GetSpecification(iter);
    int SM = AliHLTTRDUtils::GetSM(spec);

    HLTInfo("Got Data from SM %i", SM);

    if(SM!=lastSM){
      if(fIncSM[SM]){
	if(fNoOfIncSM<fNoOfSM)
	  return 0;
	fNoOfSM=fNoOfIncSM;
	PushBack(fOutArray, AliHLTTRDDefinitions::fgkCalibrationDataType);
	fOutArray->Delete();
	delete fOutArray;
	fOutArray = NULL;
	for(int i=0; i<18; i++)
    	  fIncSM[i]=kFALSE;
	fNoOfIncSM=0;
      }
      lastSM = SM;
      fIncSM[SM]=kTRUE;
      fNoOfIncSM++;
    }

    if(!fOutArray) fOutArray = (TObjArray*)iter->Clone();
    else{
      TObjArray* inArr = (TObjArray*)iter;
      for(int i = inArr->GetEntriesFast(); i--;){
	const TH1* histo = dynamic_cast<const TH1*>(inArr->At(i));
	if(histo){
	  if(fOutArray->At(i)){
	    ((TH1*)fOutArray->At(i))->Add(histo);
	  }else{
	    fOutArray->AddAt(histo->Clone(), i);
	  }
	  continue;
	}
	AliTRDCalibraVdriftLinearFit* obj = dynamic_cast<AliTRDCalibraVdriftLinearFit*>(inArr->At(i));
	if(obj){
	  if(fOutArray->At(i)){
	    ((AliTRDCalibraVdriftLinearFit*)fOutArray->At(i))->Add(obj);
	  }else{
	    fOutArray->AddAt(new AliTRDCalibraVdriftLinearFit(*obj), i);
	  }
	}
      }
    }

  }

  return 0;
}

Int_t AliHLTTRDCalibFitComponent::ShipDataToFXS(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
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

Int_t AliHLTTRDCalibFitComponent::EORCalibration()
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
  //Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
  //  + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
  //Int_t nbfit       = calibra->GetNumberFit();
  //Int_t nbE         = calibra->GetNumberEnt();
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
  //nbtg = 540;
  //nbfit = calibra->GetNumberFit();
  //nbE   = calibra->GetNumberEnt();
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
  //nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
  //  + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
  //nbfit        = calibra->GetNumberFit();
  //nbE          = calibra->GetNumberEnt();
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
  //nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
  //  + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
  //nbfit        = calibra->GetNumberFit();
  //nbE          = calibra->GetNumberEnt();
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
  fAfterRunArray->Add(coefgain);
  fAfterRunArray->Add(coefprf);
  fAfterRunArray->Add(coefdrift);
  fAfterRunArray->Add(coefdriftsecond);
  fAfterRunArray->Add(coeft0);
  
  PushBack(fAfterRunArray, AliHLTTRDDefinitions::fgkEORCalibrationDataType);

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

