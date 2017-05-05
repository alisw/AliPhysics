// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTT0RecoComponent.cxx
    @author  Alla Maevskaya <Alla.Maevskaya@cern.ch>
    @brief   T0 reconstruction component
*/

#include "TTree.h"
#include "TMap.h"
#include "TObjString.h"
#include "TDatime.h"

#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliRawReaderMemory.h"
#include "AliGeomManager.h"

#include "AliT0RecoParam.h"
#include "AliT0Reconstructor.h"
#include "AliT0RawReader.h"

#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTT0RecoComponent.h"
#include "AliESDVertex.h"
#include "AliT0CalibWalk.h"
#include "AliT0CalibTimeEq.h"


using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTT0RecoComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTT0RecoComponent::AliHLTT0RecoComponent() :
  AliHLTProcessor(),
  fRunInfo(NULL),  
  fT0RecoParam(NULL),
  fT0Reconstructor(NULL),
  fRawReader(NULL),
  fESDTZERO(NULL),
  fNevent(0),
  fVertexSPDz(0.),
  fWalk(NULL),
  fT0CalibHisto(NULL)
{
	for (int i = 0;i < 24;i++) fhTimeDiff[i] = NULL;
	for (int i = 0;i < 24;i++) fhCFD[i] = NULL;
	for (int i = 0;i < 3;i++) fhT0[i] = NULL;
	for (int i = 0;i < 24;i++) fMeanCFD[i] = 0.f;
    for (int i = 0;i < 24;i++) fDiffCFD[i] = 0.f;
    for (int i = 0;i < 3;i++) fT0shift[i] = 0.f;

}

// #################################################################################
AliHLTT0RecoComponent::~AliHLTT0RecoComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTT0RecoComponent::GetComponentID() { 
  // see header file for class documentation
  return "T0Reconstruction";
}

// #################################################################################
void AliHLTT0RecoComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0);
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITSSPD);
}

// #################################################################################
AliHLTComponentDataType AliHLTT0RecoComponent::GetOutputDataType() 
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTT0RecoComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{ 
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeESDContent|kAliHLTDataOriginT0);
  return tgtList.size();
}

// #################################################################################
void AliHLTT0RecoComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 3000;
  inputMultiplier = 0.7;
}

// #################################################################################
void AliHLTT0RecoComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  //targetMap->Add(new TObjString("HLT/ConfigT0/T0Reconstruction"),
  //		 new TObjString("configuration object"));

  targetMap->Add(new TObjString("GRP/GRP/Data"),
		 new TObjString("GRP object - run information"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/CTP/TimeAlign"),
		 new TObjString("GRP object - CTP information"));
  targetMap->Add(new TObjString("GRP/Calib/LHCClockPhase"),
		 new TObjString("GRP object - time calibration"));

  targetMap->Add(new TObjString("T0/Calib/TimeDelay"),
		 new TObjString("T0 calibration object"));
  targetMap->Add(new TObjString("T0/Calib/TimeAdjust"),
		 new TObjString("T0 calibration object"));
   targetMap->Add(new TObjString("T0/Calib/Slewing_Walk"),
		 new TObjString("T0 calibration object"));

  return;
}

// #################################################################################
AliHLTComponent* AliHLTT0RecoComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTT0RecoComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTT0RecoComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation
  //cout<<"\n\n\nVZero Reconstruction Init\n\n\n"<<endl;
 
  Int_t iResult=0;

  // -- Load GeomManager
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  
  // -- Read configuration object : HLT/ConfigT0/T0Reconstruction
  TString cdbPath="HLT/ConfigT0/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath, NULL, true);

  // -- Read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  // -- Get AliRunInfo variables
  // -----------------------------
  TObject* pOCDBEntry=LoadAndExtractOCDBObject("GRP/GRP/Data");
  AliGRPObject* pGRP=pOCDBEntry?dynamic_cast<AliGRPObject*>(pOCDBEntry):NULL;
  
  TString beamType = "";
  TString lhcState = "";
  TString runType = "";
  Float_t beamEnergy = 0.;
  UInt_t activeDetectors = 0;
  
  if (pGRP) {
    lhcState        = pGRP->GetLHCState(); 	  	   
    beamType        = pGRP->GetBeamType(); 
    runType         = pGRP->GetRunType(); 
    beamEnergy      = pGRP->GetBeamEnergy();
    activeDetectors = pGRP->GetDetectorMask();
  }
  //Get T0 slewing graphs
  TObject* pOCDBEntryT0slew=LoadAndExtractOCDBObject("T0/Calib/Slewing_Walk");
  AliT0CalibWalk* pT0slew = pOCDBEntryT0slew?dynamic_cast<AliT0CalibWalk*>(pOCDBEntryT0slew):NULL;
  if (pT0slew) {
    fWalk = new TObjArray(24);
   for (Int_t i=0; i<24; i++){
    TGraph* fu = pT0slew->GetWalk(i);
      fWalk->AddAtAndExpand(fu,i);
     }
  }
 //Get T0 time delay
     TObject* pOCDBEntryT0delay=LoadAndExtractOCDBObject("T0/Calib/TimeDelay");
    AliT0CalibTimeEq* pTimeDelay = pOCDBEntryT0delay?dynamic_cast<AliT0CalibTimeEq*>(pOCDBEntryT0delay):NULL;
    for (Int_t i=0; i<24; i++) {
      fMeanCFD[i] = pTimeDelay->GetCFDvalue(i,0);
      fDiffCFD[i] = pTimeDelay->GetTimeEq(i);
    }

  // -- Initialize members
  // -----------------------
  do {
    if (iResult<0) break;

    fRawReader = new AliRawReaderMemory;
    if (!fRawReader) {
      iResult=-ENOMEM;
      break;
    }

    // AliGRPManager grpMan;
    // Bool_t status       = grpMan.ReadGRPEntry(); // Read the corresponding OCDB entry
    // status              = grpMan.SetMagField();  // Set global field instanton
    // AliRunInfo *runInfo = grpMan.GetRunInfo();   // Get instance of run info

    fRunInfo = new AliRunInfo(lhcState.Data(), beamType.Data(),
			      beamEnergy, runType.Data(), activeDetectors);
    if (!fRunInfo) {
      iResult=-ENOMEM;
      break;
    }

    fT0RecoParam = new AliT0RecoParam;
    if (!fT0RecoParam) {
      iResult=-ENOMEM;
      break;
    }  

    fT0Reconstructor = new AliT0Reconstructor;
    if (!fT0Reconstructor) {
      iResult=-ENOMEM;
      break;
    }

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fRawReader) 
      delete fRawReader;
    fRawReader = NULL;

    if (fT0RecoParam)
      delete fT0RecoParam;
    fT0RecoParam = NULL;

    if (fT0Reconstructor)
      delete fT0Reconstructor;
    fT0Reconstructor = NULL;

    if (fRunInfo)
      delete fRunInfo;
    fRunInfo = NULL;
	
	delete fWalk;
	fWalk = NULL;
  }

  if (iResult>=0) {
    fT0Reconstructor->SetRunInfo(fRunInfo);
    fT0Reconstructor->Init();

    //alla    fT0Reconstructor->SetRecoParam(fT0RecoParam);
   fT0CalibHisto  = new TObjArray(100);
  for (Int_t i=0; i<24; i++) {
    fhTimeDiff[i]   = new TH1F (Form("CFD1minCFD%d",i+1),"fTimeDiff",100, -300, 300);
    fhCFD[i]        = new TH1F(Form("CFD%d",i+1),"CFD",100, 2000, 3000);
  }
  
  fhT0[0] = new TH1F("fTzeroORAplusORC","T0A+T0C /2",200,-2000,2000);   //or A plus or C 
  fhT0[1] = new TH1F("fTzeroORA","T0A",200,-2000,2000);// or A spectrum
  fhT0[2] = new TH1F("fTzeroORC","T0C",200,-2000,2000);// or C spectrum
  for (Int_t i=0; i<24; i++) {
    fT0CalibHisto ->AddAtAndExpand( fhTimeDiff[i],i+24); //24 - 48
    fT0CalibHisto ->AddAtAndExpand(fhCFD[i], i); //24 - 48
  }  
  for (Int_t i=0; i<3; i++)
    fT0CalibHisto ->AddAtAndExpand(fhT0[i],i+48); // 48 -52

  fT0CalibHisto->Print();
  fNevent=0;
  }
  
  return iResult;
}

// #################################################################################
Int_t AliHLTT0RecoComponent::ScanConfigurationArgument(Int_t /*argc*/, const Char_t** argv) {
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  Int_t ii =0;
  TString argument=argv[ii];

  if (argument.IsNull()) return 0;

  return 0;
}

// #################################################################################
Int_t AliHLTT0RecoComponent::DoDeinit() {
  // see header file for class documentation

  if (fRawReader) 
    delete fRawReader;
  fRawReader = NULL;
  
  if (fT0RecoParam)
    delete fT0RecoParam;
  fT0RecoParam = NULL;
  
  if (fT0Reconstructor)
    delete fT0Reconstructor;
  fT0Reconstructor = NULL;
  
  if (fRunInfo)
    delete fRunInfo;
  fRunInfo = NULL;
  
  delete fWalk;
  fWalk = NULL;
  
  return 0;
}

// #################################################################################
Int_t AliHLTT0RecoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;
 
  //cout<<"\n\n\nVZero Reconstruction Do Event\n\n\n"<<endl;

  // -- Get T0 raw dat a input block and set up the rawreader
  const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0);
  if (!pBlock) {
    //cout<<"No T0 input block at event"<<endl;
    HLTInfo("No T0 input block at event %d", GetEventCount());
    return 0;
  }
  
  //cout<<"T0 input block found"<<endl;
 
  // -- Add input block to raw reader
  if (!fRawReader->SetMemory((UChar_t*) pBlock->fPtr, pBlock->fSize )){
    cout<<"Could not add buffer of data block to rawreader"<<endl;
    HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	     DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification);
    iResult = -1;
  }
  

    //read ESD vertex   
  AliESDVertex *vtx=0;
  {
    const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
    if (iter == NULL) iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITSSPD);
    if(iter) vtx = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) );   
  }
  //printf("@@@@@@@@ vertex  \n");
  if(vtx) { 
    vtx->Print();
    fVertexSPDz = vtx->GetZ();
  }
  else HLTInfo("No vertex from SPD in this event...");

  if (iResult >= 0) {
    //cout<<"ok 1"<<endl;
    // -- Set T0 EquipmentID
    fRawReader->SetEquipmentID(3328);
     RecT0Raw (fRawReader);
 
    /*
    // -- 1. step T0 reconstruction
    fT0Reconstructor->ConvertDigits(fRawReader, digitsTree);

    // -- 2. step T0 reconstruction -- fill AliESDT0 object
    fT0Reconstructor->FillESD(digitsTree, NULL, NULL);

    AliESDT0 *esdT0 = fT0Reconstructor->GetESDT0();
      
    // Send info every 10 s
    const TDatime time;    
    static UInt_t lastTime=0;
    if (time.Get()-lastTime>10) {
      lastTime=time.Get();
      HLTInfo("T0 Multiplicity A %f - C %f", esdT0->GetMTotV0A(), esdT0->GetMTotV0A() );
    }
    */
    // -- Send AliESDT0 & friend object
     PushBack(fESDTZERO, kAliHLTDataTypeESDContent|kAliHLTDataOriginT0,0);
  }
  
  // -- Clean up

 fRawReader->ClearBuffers();   

 return iResult;
}

// #################################################################################
Int_t AliHLTT0RecoComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigT0/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTT0RecoComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  return 0;
}
void AliHLTT0RecoComponent::RecT0Raw(AliRawReader *rawReader)
{
  //printf("@@@@@@@@@@@@@@2 AliHLTT0RecoComponent::RecT0Raw \n");
  
  Int_t timeCFD[24], chargeQT0[24], chargeQT1[24];
  Float_t  meanA, meanC, meanAC, sigmacfd[24];
  Double32_t qt[24], time[24];
  Int_t fRefPMTA=12, fRefPMTC=0, nent=0;
  Int_t allData[250][5];
  for (Int_t i0=0; i0<250; i0++)
    for (Int_t j0=0; j0<5; j0++)  allData[i0][j0]=0; 
  
  Double32_t besttimeA=9999999; 
  Double32_t besttimeC=9999999; 
  Float_t channelWidth=24.4; 
  Float_t c = 0.0299792458; // cm/ps
  Float_t shift=0;
  Float_t sigmaT0shift[3];
  for(int ipmt=0; ipmt<24; ipmt++)  
    sigmacfd[ipmt]=time[ipmt]=qt[ipmt]=0;
  for(int i=0; i<3; i++) sigmaT0shift[i]=0;
  
  
  AliT0RawReader myrawreader(rawReader);
  UInt_t type =rawReader->GetType();
  if (!myrawreader.Next())
    HLTInfo(" no raw data found!!");
  else
    {  
      for (Int_t i=0; i<24; i++) { 
	timeCFD[i]=0; chargeQT0[i]=0; chargeQT1[i]=0;
      }
      for (Int_t i=0; i<226; i++) {
	for (Int_t iHit=0; iHit<5; iHit++) 
	  {
	    allData[i][iHit] = myrawreader.GetData(i,iHit);
	  }
      }
      //trigger selection      if(allData[50][0]==0) return; 
      //   if (allData[50][0] > 0)
      {
	fNevent++;
	//printf("@@@@@@event %i \n", fNevent);
	for (Int_t in=0; in<12; in++)  
	  {
	    chargeQT0[in]=allData[2*in+25][0];
	    chargeQT1[in]=allData[2*in+26][0];
	    if ((allData[in+1][0] >(fMeanCFD[in]-100)) && 
		(allData[in+1][0]<(fMeanCFD[in]+100)) ) { 
	      timeCFD[in] = allData[in+1][0] ;
	      if(chargeQT1[in] >(fMeanCFD[in]+2000) &&
		 chargeQT1[in] <(fMeanCFD[in]+3000) ) {
		qt[in] = Float_t(chargeQT0[in]  - chargeQT1[in]);
		//printf(" readed Raw  pmt %i QT0 %i QT1 %i qt %f\n", in, chargeQT0[in],chargeQT1[in], qt[in]);
	      }
	    }
	  }
	for (Int_t in=12; in<24; in++)  
	  {
	    timeCFD[in] = allData[in+45][0] ;
	    chargeQT0[in]=allData[2*in+57][0];
	    
	    chargeQT1[in]=allData[2*in+58][0];
	    if ((allData[in+45][0] >(fMeanCFD[in]-100)) && 
		(allData[in+45][0]<(fMeanCFD[in]+100)) ) {
	      timeCFD[in] = allData[in+45][0] ;
	      if (chargeQT1[in] >(fMeanCFD[in]+2000) &&
		  chargeQT1[in] <(fMeanCFD[in]+3000) ) {
		qt[in] = Float_t(chargeQT0[in]  - chargeQT1[in]);
		//printf(" readed Raw  pmt %i QT0 %i QT1 %i qt %f\n", in, chargeQT0[in],chargeQT1[in], qt[in]);
	      }
	    }
	  }
	Double_t walk=0;
	Float_t sigmadiff[24];
	for (Int_t ipmt=0; ipmt<24; ipmt++)
	  {
	    time[ipmt] =0;
	    if(timeCFD[ipmt] >0 && qt[ipmt]>0 )
	      {
		TGraph *fu1=(TGraph*) fWalk->At(ipmt);
		if(fu1 && fu1->GetN()>0) {
		  walk = Int_t(fu1->Eval(qt[ipmt]));
		}
		timeCFD[ipmt] = timeCFD[ipmt] - Int_t (walk) -Int_t (fMeanCFD[ipmt]);
		if(  timeCFD[fRefPMTC] > 0 && ipmt<12) 
		  fhTimeDiff[ipmt] ->Fill(timeCFD[ipmt]-timeCFD[fRefPMTC]); 
		if(  timeCFD[fRefPMTA] > 0 && ipmt>11) 
		  fhTimeDiff[ipmt] ->Fill(timeCFD[ipmt]-timeCFD[fRefPMTA]); 
		fhCFD[ipmt]->Fill( timeCFD[ipmt]);
		nent=fhCFD[ipmt]->GetEntries();
		//printf("ipmt %i  timeCFD %i  walk %f fMeanCFD %f entries %i \n", ipmt, timeCFD[ipmt],walk, fMeanCFD[ipmt], nent);
		
		if (fNevent==1000) {
		  nent=fhCFD[ipmt]->GetEntries();
		  if (nent>200) 
		    GetMeanAndSigma( fhCFD[ipmt], fMeanCFD[ipmt], sigmacfd[ipmt]);
		  //	      GetMeanAndSigma( fhTimeDiff[ipmt], fDiffCFD[ipmt], sigmadiff[ipmt]);
		  //printf("ipmt %i  meancfd %f  sigmacfd %f entries %i \n", ipmt, fMeanCFD[ipmt], sigmacfd[ipmt], nent );
		}
	  
		
		if (fNevent>1000) {
		  if( (timeCFD[ipmt])<besttimeC && ipmt<12) besttimeC = timeCFD[ipmt]; //timeC
		  if( (timeCFD[ipmt])<besttimeA && ipmt>11) besttimeA = timeCFD[ipmt]; //timeA
		  time[ipmt]=timeCFD[ipmt]+fMeanCFD[ipmt];
		}
	      }	    
	  }
  
 	
	if(fNevent<1000)
	  for(int iii=0; iii<3; iii++) fT0shift[iii]=0;
	if (fNevent>1000) 
	  {	  
	    fESDTZERO->SetT0time(time);
	    fESDTZERO->SetT0amplitude(qt);   
	    if(besttimeA < 999999 && besttimeA!=0 &&  besttimeC < 999999 && besttimeC!=0 ) {
	      meanAC=channelWidth* (besttimeA+besttimeC)/2.-fT0shift[0];
	      fhT0[0]->Fill(meanAC);
	      //printf(" BEST A %f C %f AC %f T0shift %f\n", besttimeA, besttimeC, meanAC, fT0shift[0]) ;

	    }
	    shift = fVertexSPDz/c;
	    if(besttimeA < 999999 && besttimeA!=0 ) {
	      meanA=channelWidth*besttimeA +shift-fT0shift[1];
	      //printf(" Best Time A %f SPDshift %f T0shift %f\n",meanA,  shift,fT0shift[1] ); 
	      fhT0[1]->Fill(meanA);
	      if(fNevent>300)  fESDTZERO->SetT0TOF(1,meanA);
	    } 	    
	    
	    if( besttimeC < 999999 && besttimeC!=0) 
	      {
		meanC=channelWidth*besttimeA -shift-fT0shift[2];
		//printf(" Best Time C %f SPDshift %f T0shift %f \n",meanC, shift,fT0shift[2]  );
		fhT0[2]->Fill(meanC);
		if(fNevent>300)  fESDTZERO->SetT0TOF(2,meanC); 	 
	      }
	  }
	
	if(fNevent==2000) {
	  //printf(" @@@event 200 \n");
	  for (int ii=0; ii<3; ii++) 
	    GetMeanAndSigma(fhT0[ii], fT0shift[ii], sigmaT0shift[ii]);
	  //printf(" @@@@@ T0AC %f T0A %f  T0C %f SPDshift %f \n",fT0shift[0] , fT0shift[1],fT0shift[2], shift  );
	}	
	
      }
    }
 
}


//________________________________________________________________________
void  AliHLTT0RecoComponent::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {
  
  const double window = 5.;  //fit window 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  // sigmaEstimate = 10;
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","R","");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}
