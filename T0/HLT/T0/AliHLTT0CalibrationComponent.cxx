// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTSampleCalibrationComponent.cxx
//  @author Matthias Richter
//  @date   2010-04-26
//  @brief  A sample calibration component for the HLT.
//  @ingroup alihlt_tutorial

#include "AliHLTT0CalibrationComponent.h"
#include "AliHLTReadoutList.h"
#include "AliLog.h"
#include "TMap.h"
#include "TObjString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "AliGRPObject.h"
#include "AliT0CalibWalk.h"
#include "AliRawReader.h"
#include "AliT0RawReader.h"
#include "AliRawReaderMemory.h"
#include "AliRunInfo.h"
#include "AliESDVertex.h"
#include <Riostream.h>
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
//#include "AliHLTT0CalibObject.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTT0CalibrationComponent)

AliHLTT0CalibrationComponent::AliHLTT0CalibrationComponent()
  : AliHLTCalibrationProcessor()
{
  // an example component which implements the ALICE HLT calibration
  // processor interface
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
  fNevent = 0;
  fVertexSPDz =0;
  // fT0CalibObject = new AliHLTT0CalibObject();
  
  cout<<"\n\n\nConstructor called!!\n\n\n"<<endl;
}

AliHLTT0CalibrationComponent::~AliHLTT0CalibrationComponent()
{
  // destructor
  //
  // NOTE: implement proper cleanup in DoDeinit()
}

const char* AliHLTT0CalibrationComponent::GetComponentID()
{ 
  // component property: id
  return "T0Calibration";
}

void AliHLTT0CalibrationComponent::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  // component property: list of input data types
  list.push_back(kAliHLTAnyDataType);
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0);
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
}

AliHLTComponentDataType AliHLTT0CalibrationComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTDataTypeFXSCalib|kAliHLTDataOriginT0;
}

void AliHLTT0CalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 10000000;
  inputMultiplier = 0.5;
  //inputMultiplier = 10.;
}

void AliHLTT0CalibrationComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  if (!targetMap) return;
  targetMap->Add(new TObjString("GRP/GRP/Data"),
		 new TObjString("GRP object - run information"));
  targetMap->Add(new TObjString("T0/Calib/Slewing_Walk"),
		 new TObjString("T0 calibration object"));

 }

AliHLTComponent* AliHLTT0CalibrationComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTT0CalibrationComponent;
}

int AliHLTT0CalibrationComponent::DoInit( int argc, const char** argv )
{
  cout<<"\n\n\n  DoInit() called!!\n\n\n"<<endl;
  // see header file for class documentation
  int iResult=0;
  
  // init stage 1: default values for all data members
  
  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  // TString cdbPath="HLT/ConfigT0/";
  //cdbPath+=GetComponentID();
  //iResult=ConfigureFromCDBTObjString(cdbPath);
  
  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }
  TString beamType = "";
  TString lhcState = "";
  TString runType = "";
  Float_t beamEnergy = 0.;
  UInt_t activeDetectors = 0;
  
  
  if (iResult>=0) {
    // implement the component initialization
    // -- Get AliRunInfo information
    // -----------------------------
    TObject* pOCDBEntryGRPdata=LoadAndExtractOCDBObject("GRP/GRP/Data");
    AliGRPObject* pGRP=pOCDBEntryGRPdata?dynamic_cast<AliGRPObject*>(pOCDBEntryGRPdata):NULL;
    
    if (pGRP) {
      lhcState        = pGRP->GetLHCState(); 	  	   
      beamType        = pGRP->GetBeamType(); 
      runType         = pGRP->GetRunType(); 
      beamEnergy      = pGRP->GetBeamEnergy();
      activeDetectors = pGRP->GetDetectorMask();
    }
    //Get T0 slewing graphs
    fWalk    = new TObjArray(0);  
    TObject* pOCDBEntryT0slew=LoadAndExtractOCDBObject("T0/Calib/Slewing_Walk");
    AliT0CalibWalk* pT0slew = pOCDBEntryT0slew?dynamic_cast<AliT0CalibWalk*>(pOCDBEntryT0slew):NULL;
    if (pT0slew) {
      for (Int_t i=0; i<24; i++){
	TGraph* fu = pT0slew->GetWalk(i);
	fWalk->AddAtAndExpand(fu,i);
      } 
    } 
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
    fRunInfo = new AliRunInfo(lhcState.Data(), beamType.Data(),
			      beamEnergy, runType.Data(), activeDetectors);
    if (!fRunInfo) {
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
    if (fRunInfo)
      delete fRunInfo;
    fRunInfo = NULL;    
  }
   fT0CalibHisto  = new TObjArray(100);
  for (Int_t i=0; i<24; i++) {
    fhTimeDiff[i]   = new TH1F (Form("CFD1minCFD%d",i+1),"fTimeDiff",100, -300, 300);
    fhCFD[i]        = new TH1F(Form("CFD%d",i+1),"CFD",100, 2000, 3000);
  }
  
  fhT0[0] =   new TH1F("fTzeroORAplusORC","ORA+ORC /2",200,-2000,2000);   //or A plus or C 
  fhT0[3]      = new TH1F("fResolution","fResolution",200,-1000,1000);// or A minus or C spectrum
  fhT0[1]        = new TH1F("fTzeroORA","fTzeroORA",200,-2000,2000);// or A spectrum
  fhT0[2]        = new TH1F("fTzeroORC","fTzeroORC",200,-2000,2000);// or C spectrum
  for (Int_t i=0; i<24; i++) {
    fT0CalibHisto ->AddAtAndExpand( fhTimeDiff[i],i+24); //24 - 48
    fT0CalibHisto ->AddAtAndExpand(fhCFD[i], i); //24 - 48
  }  
  for (Int_t i=0; i<4; i++)
    fT0CalibHisto ->AddAtAndExpand(fhT0[i],i+48); // 48 -52
  fhSPDvertex = new TH1F("hSPDvertex","SPDvertex", 100, -20,20);
  fhT0vertex = new TH1F("hT0vertex","T0vertex", 100, -20,20);
  fT0CalibHisto ->AddAtAndExpand(fhT0vertex,52); // 48 -52
  fT0CalibHisto ->AddAtAndExpand(fhSPDvertex,53); // 48 -52
  fhT0SPDvertex = new TH2F("hT0SPDvertex","hT0SPDvertex", 100, -20, 20, 100, -20, 20);
  fT0CalibHisto ->AddAtAndExpand(fhT0SPDvertex,54); // 48 -52

  fT0CalibHisto->Print();
  fNevent=0;
  
  printf(" @@@@@@@@@@@@DoInit finished\n");
  
  return iResult;
}

int AliHLTT0CalibrationComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.
  
  int i=0;
  TString argument=argv[i];
  
  if (argument.IsNull()) return 0;
  
  // -mandatory1 arg
  if (argument.CompareTo("-mandatory1")==0) {
    if (++i>=argc) return -EPROTO;
    HLTInfo("got \'-mandatory1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }
  
  // -optional1 arg
  if (argument.CompareTo("-optional1")==0) {
    if (++i>=argc) return -EPROTO;
    HLTInfo("got \'-optional1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -optional2
  if (argument.CompareTo("-optional2")==0) {
    HLTInfo("got \'-optional2\' argument");
    return 1; // only keyword
  }

  return 0;
}

int AliHLTT0CalibrationComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fRawReader) 
    delete fRawReader;
  fRawReader = NULL;
  if (fRunInfo)
    delete fRunInfo;
  fRunInfo = NULL;
  
  
  return 0;
}

int AliHLTT0CalibrationComponent::ProcessCalibration(const AliHLTComponentEventData& /*evtData*/,AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  
  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;
  // see header file for class documentation
  printf("@@@@@@@@  AliHLTT0CalibrationComponent::ProcessCalibration \n");
  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;
  
  AliESDVertex *vtx=0;
  {
    const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
    if(iter) vtx = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) );   
  }
  printf("@@@@@@@@ vertex  \n");
  if(vtx) { 
    vtx->Print();
    fVertexSPDz = vtx->GetZ();
    fhSPDvertex->Fill(fVertexSPDz);
    Printf("@@@ Vertex SPD %f ",fVertexSPDz);
  }
  else printf("\nWhat the fuck? No vertex from SPD...\n\n");
  // -- Get T0 raw dat a input block and set up the rawreader
  const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginT0);
  if (!pBlock) {
    ALIHLTERRORGUARD(1, "No T0 input block at event %d", GetEventCount());
    return 0;
  }
  
  // -- Add input block to raw reader
  if (!fRawReader->SetMemory((UChar_t*) pBlock->fPtr, pBlock->fSize )){
    HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	     DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification);
    iResult = -1;
  }
  
  if (iResult >= 0) {
    
    // -- Set T0 EquipmentID
    fRawReader->SetEquipmentID(3328); //??????????
    RecT0Raw (fRawReader);
    //	cout<<fRawReader<<" "<<myrawreader<<endl;
    
    // write the histogram out
    // this should not be done for every event, however the call can be implemented
    // like that and the publishing steered by the component argument 
    // '-pushback-period=...'
    //    if (PushBack(fhT0[0], kAliHLTDataTypeHistogram)==-ENOSPC) {
    //    printf("@@@@@@@!!!!!!!!!! (PushBack(fT0[0] \n");
    //}
    // increase the output size estimator
    // we add the size of the last object, there might be other blocks to
    // be written in addition to the actual object
    // -- Send AliESDT0
    //PushBack(esdT0, kAliHLTDataTypeESDContent|kAliHLTDataOriginT0,0);
  }
  
// -- Clean up
  fRawReader->ClearBuffers();   
  
  return iResult;
}


int AliHLTT0CalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, 
						     AliHLTComponentTriggerData& /*trigData*/)
{
  // prepare final result and ship to FXS
 
 	      //
  printf(" @@@ ShipDataToFXS fit start\n");  
  Float_t sigmaT0shift=0;

  for (int ii=0; ii<4; ii++) {
    Int_t nent=fhT0[ii]->GetEntries();
    GetMeanAndSigma(fhT0[ii], fT0shift[ii], sigmaT0shift);
    Printf("ii %i nent %i  shift %f sigma %f",ii, nent,  fT0shift[ii],  sigmaT0shift); 
    // fT0CalibObject->SetT0calibParams(ii+48, fT0shift[ii]);
  }
  /*
  for (int i=0; i<24; i++) {
    fT0CalibObject->SetT0calibParams(i, fMeanCFD[i]);
    fT0CalibObject->SetT0calibParams(i, fDiffCFD[i+24]);
  }
*/
  AliHLTReadoutList rdList(AliHLTReadoutList::kHLT);
  Printf("Pushing to FXS...");
  //  PushToFXS(fT0CalibObject, "HLT", "TestHisto", &rdList);
  Printf("... done.");
  TFile* fout = new TFile("foutT0.root", "RECREATE");
  Printf("@@ file open ");
  //  for (int i=0; i<24; i++)     fhCFD[i] ->Write();
  //  for (int ii=0; ii<3; ii++) fhT0[ii]->Write();
   fT0CalibHisto->Write();
  fout->Close();

  return 0;
}

int AliHLTT0CalibrationComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  HLTInfo("reconfigure '%s' from entry %s", chainId, cdbEntry);

  return 0;
}

int AliHLTT0CalibrationComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
// #################################################################################
void AliHLTT0CalibrationComponent::RecT0Raw(AliRawReader *rawReader)
{
  printf("@@@@@@@@@@@@@@2 AliHLTT0CalibrationComponent::RecT0Raw \n");
  
  Int_t timeCFD[24], chargeQT0[24], chargeQT1[24];
  Float_t meanA, meanC, meanAC, sigmacfd[24];
  Double32_t qt[24];
  Int_t fRefPMTA=12, fRefPMTC=0, nent=0;
  Int_t allData[110][5];
  for (Int_t i0=0; i0<110; i0++)
    for (Int_t j0=0; j0<5; j0++)  allData[i0][j0]=0; 
  
  Float_t besttimeA=9999999; 
  Float_t besttimeC=9999999; 
  Float_t channelWidth=24.4; 
  Float_t c = 0.0299792458; // cm/ps
  Float_t shift=0;
  Float_t sigmaT0shift[3];
  for(int ipmt=0; ipmt<24; ipmt++)  
    sigmacfd[ipmt]=qt[ipmt]=0;
  meanA= meanC = meanAC =0;
  AliT0RawReader myrawreader(rawReader);
  UInt_t type =rawReader->GetType();
  if (!myrawreader.Next())
    printf(" no raw data found!!");
  else
    {  
      for (Int_t i=0; i<24; i++) { 
	timeCFD[i]=0; chargeQT0[i]=0; chargeQT1[i]=0;
	if(fNevent<200) fMeanCFD[i]=0;
      }
      for (Int_t i=0; i<107; i++) {
	for (Int_t iHit=0; iHit<5; iHit++) 
	  {
	    allData[i][iHit] = myrawreader.GetData(i,iHit);
	    if(iHit==0 && i==50) printf(" @@@@@@@@@@ tvdc  %i \n",allData[i][iHit]);
	  }
      }
      if (allData[50][0] > 0) {
	  fNevent++;
	  printf("@@@@@@event %i \n", fNevent);
	  for (Int_t in=0; in<12; in++)  
	    {
	      timeCFD[in] = allData[in+1][0] ; 
	      if (timeCFD[in]>0) {
		printf(" readed i %i cfdC %i  \n ", in, timeCFD[in])	; 
		chargeQT0[in]=allData[2*in+25][0];
		chargeQT1[in]=allData[2*in+26][0];
	        qt[in] = Float_t(chargeQT0[in]  - chargeQT1[in]);
  	        printf(" readed Raw pmt %i QT0 %i QT1 %i qt %f\n",
		       in, chargeQT0[in],chargeQT1[in], qt[in]);
	      }
	    }
	  
	  for (Int_t in=12; in<24; in++)  
	    {
	      timeCFD[in] = allData[in+45][0] ;
	      if( timeCFD[in]>0) {
		printf(" readed i %i cfdA %i \n", in, timeCFD[in]);		    
		chargeQT0[in]=allData[2*in+57][0];
		chargeQT1[in]=allData[2*in+58][0];
		qt[in] = Float_t(chargeQT0[in]  - chargeQT1[in]);
		printf(" readed Raw  pmt %i QT0 %i QT1 %i qt %f\n",
		       in, chargeQT0[in],chargeQT1[in], qt[in]);
	      }
	    }
	  Double_t walk=0;
	  Float_t  sigmadiff[24] ;
	  for (Int_t ipmt=0; ipmt<24; ipmt++) {
	    {
	      if(fNevent>300 && 
		 !( (timeCFD[ipmt] >(fMeanCFD[ipmt]-100)) && 
		    (timeCFD[ipmt]<(fMeanCFD[ipmt]+100))) ) continue;
	      
	      if(timeCFD[ipmt] >0 && qt[ipmt]>0)
		{
		  TGraph *fu1=(TGraph*) fWalk->At(ipmt);
		  if(fu1 && fu1->GetN()>0) {
		    walk = Int_t(fu1->Eval(qt[ipmt]));
		  }
		  timeCFD[ipmt] = timeCFD[ipmt] - Int_t (walk) -Int_t (fMeanCFD[ipmt]);
		  if(fNevent<300) {
		    if(  timeCFD[fRefPMTC] > 0 && ipmt<12) 
		      fhTimeDiff[ipmt] ->Fill(timeCFD[ipmt]-timeCFD[fRefPMTC]); 
		    if(  timeCFD[fRefPMTA] > 0 && ipmt>11) 
		      fhTimeDiff[ipmt] ->Fill(timeCFD[ipmt]-timeCFD[fRefPMTA]); 
		    fhCFD[ipmt]->Fill( timeCFD[ipmt]);
		    nent=fhCFD[ipmt]->GetEntries();
		    printf("ipmt %i  timeCFD %i  walk %f fMeanCFD %f entries %i \n",
			   ipmt, timeCFD[ipmt],walk, fMeanCFD[ipmt], nent);
		  }
		  if (fNevent>300) {
		    if( (timeCFD[ipmt])<besttimeC && ipmt<12) besttimeC = timeCFD[ipmt]; //timeC
		    if( (timeCFD[ipmt])<besttimeA && ipmt>11) besttimeA = timeCFD[ipmt]; //timeA
		  }
		}
	    }
	  }
	  if (fNevent==300) {
	    for (Int_t ipmt=0; ipmt<24; ipmt++) {
	      nent=fhCFD[ipmt]->GetEntries();
	      GetMeanAndSigma( fhCFD[ipmt], fMeanCFD[ipmt], sigmacfd[ipmt]);
	      GetMeanAndSigma( fhTimeDiff[ipmt], fDiffCFD[ipmt], sigmadiff[ipmt]);
	      printf("ipmt %i  meancfd %f  sigmacfd %f entries %i \n",
		     ipmt, fMeanCFD[ipmt], sigmacfd[ipmt], nent );
	    }
	  }
	  //	  if(fNevent<400)
	  //    for(int iii=0; iii<3; iii++) fT0shift[iii]=0;
	  if (fNevent>300 ) 
	    {	  
	      //      printf("!!!!fNevent %i \n", fNevent);
	      shift =   fVertexSPDz/c;
	      if(besttimeA < 999999 && besttimeA!=0 ) {
	        meanA=channelWidth*Float_t(besttimeA) + shift;
		printf(" Best Time A %f SPDshift %f T0shift %f\n",meanA,  shift ); 
		fhT0[1]->Fill(meanA);
	      } 	    
	      
	      if( besttimeC < 999999 && besttimeC!=0) 
		{
		   meanC=channelWidth*besttimeC -shift;
		   printf(" Best Time C %f SPDshift %f T0shift %f \n",meanC, shift  );
		  fhT0[2]->Fill(meanC);
		}
	      if(besttimeA < 999999 && besttimeA!=0 &&  besttimeC < 999999 && besttimeC!=0 ) {
		meanAC=channelWidth* (besttimeA+besttimeC)/2.;
		fhT0[0]->Fill(meanAC);
		nent=fhT0[0]->GetEntries();
		printf(" BEST A %f C %f AC %f T0shift %f entries %i\n", besttimeA, besttimeC, meanAC,  nent) ;
		fhT0[3]->Fill( (meanA-meanC)/2.);
		Float_t t0vertex = c* ( besttimeA - besttimeC)* channelWidth/2. ;
		//	Float_t t0vertex = c* ( besttimeC - besttimeA)* channelWidth/2. ;
		fhT0vertex->Fill(t0vertex);
		fhT0SPDvertex->Fill(t0vertex, fVertexSPDz);

	      }
	    }
	  
	  
      }
    }
  
}
//________________________________________________________________________
void  AliHLTT0CalibrationComponent::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {
  
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
