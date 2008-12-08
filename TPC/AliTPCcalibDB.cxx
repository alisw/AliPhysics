/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
// Request an instance with AliTPCcalibDB::Instance()                        //
// If a new event is processed set the event number with SetRun              //
// Then request the calibration data                                         ////
//
//
// Calibration data:
// 0.)  Altro mapping
//          Simulation      - not yet 
//          Reconstruction  - AliTPCclustererMI::Digits2Clusters(AliRawReader* rawReader)
//
// 1.)  pad by pad calibration -  AliTPCCalPad
//      
//      a.) fPadGainFactor
//          Simulation: AliTPCDigitizer::ExecFast - Multiply by gain
//          Reconstruction : AliTPCclustererMI::Digits2Clusters - Divide by gain  
//
//      b.) fPadNoise -
//          Simulation:        AliTPCDigitizer::ExecFast
//          Reconstruction:    AliTPCclustererMI::FindClusters(AliTPCCalROC * noiseROC)
//                             Noise depending cut on clusters charge (n sigma)
//      c.) fPedestal:
//          Simulation:     Not used yet - To be impleneted - Rounding to the nearest integer
//          Reconstruction: Used in AliTPCclustererMI::Digits2Clusters(AliRawReader* rawReader) 
//                          if data taken without zero suppression  
//                          Currently switch in  fRecoParam->GetCalcPedestal();
//      
//      d.) fPadTime0
//          Simulation:      applied in the AliTPC::MakeSector - adding offset
//          Reconstruction:  AliTPCTransform::Transform() - remove offset
//                           AliTPCTransform::Transform() - to be called
//                           in AliTPCtracker::Transform()      
//
// 
// 2.)  Space points transformation:
//
//      a.) General coordinate tranformation - AliTPCtransform (see $ALICE_ROOT/TPC/AliTPCtransform.cxx)
//          Created on fly - use the other calibration components
//                 Unisochronity  - (substract time0 - pad by pad)
//                 Drift velocity - Currently common drift velocity - functionality of AliTPCParam
//                 ExB effect    
//          Simulation     - Not used directly (the effects are applied one by one (see AliTPC::MakeSector)
//          Reconstruction - 
//                           AliTPCclustererMI::AddCluster
//                           AliTPCtrackerMI::Transform
//      b.) ExB effect calibration - 
//             classes (base class AliTPCExB, implementation- AliTPCExBExact.h  AliTPCExBFirst.h)
//             a.a) Simulation:   applied in the AliTPC::MakeSector - 
//                                calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
//             a.b) Reconstruction -  
//                  
//                  in AliTPCtransform::Correct() - called calib->GetExB()->Correct(dxyz0,dxyz1)
//
//  3.)   cluster error, shape and Q parameterization
//
//
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>


#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliLog.h>
#include <AliMagF.h>
#include <AliMagWrapCheb.h>

#include "AliTPCcalibDB.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCExB.h"

#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCSensorTempArray.h"
#include "AliGRPObject.h"
#include "AliTPCTransform.h"

class AliCDBStorage;
class AliTPCCalDet;
//
//

#include "TFile.h"
#include "TKey.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibCE.h"
#include "AliTPCExBFirst.h"
#include "AliTPCTempMap.h"
#include "AliTPCCalibVdrift.h"




ClassImp(AliTPCcalibDB)

AliTPCcalibDB* AliTPCcalibDB::fgInstance = 0;
Bool_t AliTPCcalibDB::fgTerminated = kFALSE;
TObjArray    AliTPCcalibDB::fgExBArray;    // array of ExB corrections


//_ singleton implementation __________________________________________________
AliTPCcalibDB* AliTPCcalibDB::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTPCcalibDB();
  
  return fgInstance;
}

void AliTPCcalibDB::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTPCcalibDB::AliTPCcalibDB():
  TObject(),
  fRun(-1),
  fTransform(0),
  fExB(0),
  fPadGainFactor(0),
  fDedxGainFactor(0),
  fPadTime0(0),
  fPadNoise(0),
  fPedestals(0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),  
  fGRPArray(100000),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(100000),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(100000),         //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(100000),
  fTemperatureArray(100000),    //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(100000),                 //! array of v drift interfaces
  fRunList(100000)              //! run list - indicates try to get the run param 

{
  //
  // constructor
  //  
  //
  Update();    // temporary
}

AliTPCcalibDB::AliTPCcalibDB(const AliTPCcalibDB& ):
  TObject(),
  fRun(-1),
  fTransform(0),
  fExB(0),
  fPadGainFactor(0),
  fDedxGainFactor(0),
  fPadTime0(0),
  fPadNoise(0),
  fPedestals(0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),
  fGRPArray(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(0),        //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(0),
  fTemperatureArray(0),   //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(0),         //! array of v drift interfaces
  fRunList(0)              //! run list - indicates try to get the run param 
{
  //
  // Copy constructor invalid -- singleton implementation
  //
   Error("copy constructor","invalid -- singleton implementation");
}

AliTPCcalibDB& AliTPCcalibDB::operator= (const AliTPCcalibDB& )
{
//
// Singleton implementation - no assignment operator
//
  Error("operator =", "assignment operator not implemented");
  return *this;
}



//_____________________________________________________________________________
AliTPCcalibDB::~AliTPCcalibDB() 
{
  //
  // destructor
  //
  
  // don't delete anything, CDB cache is active!
  //if (fPadGainFactor) delete fPadGainFactor;
  //if (fPadTime0) delete fPadTime0;
  //if (fPadNoise) delete fPadNoise;
}


//_____________________________________________________________________________
AliCDBEntry* AliTPCcalibDB::GetCDBEntry(const char* cdbPath)
{
  // 
  // Retrieves an entry with path <cdbPath> from the CDB.
  //
  char chinfo[1000];
    
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(cdbPath, fRun); 
  if (!entry) 
  { 
    sprintf(chinfo,"AliTPCcalibDB: Failed to get entry:\t%s ", cdbPath);
    AliError(chinfo); 
    return 0; 
  }
  return entry;
}


//_____________________________________________________________________________
void AliTPCcalibDB::SetRun(Long64_t run)
{
  //
  // Sets current run number. Calibration data is read from the corresponding file. 
  //  
  if (fRun == run)
    return;  
  fRun = run;
  Update();
}
  


void AliTPCcalibDB::Update(){
  //
  AliCDBEntry * entry=0;
  
  Bool_t cdbCache = AliCDBManager::Instance()->GetCacheFlag(); // save cache status
  AliCDBManager::Instance()->SetCacheFlag(kTRUE); // activate CDB cache
  
  //
  entry          = GetCDBEntry("TPC/Calib/PadGainFactor");
  if (entry){
    //if (fPadGainFactor) delete fPadGainFactor;
    entry->SetOwner(kTRUE);
    fPadGainFactor = (AliTPCCalPad*)entry->GetObject();
  }
  //
  entry          = GetCDBEntry("TPC/Calib/GainFactorDedx");
  if (entry){
    entry->SetOwner(kTRUE);
    fDedxGainFactor = (AliTPCCalPad*)entry->GetObject();
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadTime0");
  if (entry){
    //if (fPadTime0) delete fPadTime0;
    entry->SetOwner(kTRUE);
    fPadTime0 = (AliTPCCalPad*)entry->GetObject();
  }
  //
  //
  entry          = GetCDBEntry("TPC/Calib/PadNoise");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fPadNoise = (AliTPCCalPad*)entry->GetObject();
  }

  entry          = GetCDBEntry("TPC/Calib/Pedestals");
  if (entry){
    //if (fPedestals) delete fPedestals;
    entry->SetOwner(kTRUE);
    fPedestals = (AliTPCCalPad*)entry->GetObject();
  }

  entry          = GetCDBEntry("TPC/Calib/Temperature");
  if (entry){
    //if (fTemperature) delete fTemperature;
    entry->SetOwner(kTRUE);
    fTemperature = (AliTPCSensorTempArray*)entry->GetObject();
  }

  entry          = GetCDBEntry("TPC/Calib/Parameters");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fParam = (AliTPCParam*)(entry->GetObject()->Clone());
  }

  entry          = GetCDBEntry("TPC/Calib/ClusterParam");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fClusterParam = (AliTPCClusterParam*)(entry->GetObject()->Clone());
  }

  entry          = GetCDBEntry("TPC/Calib/Mapping");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    TObjArray * array = dynamic_cast<TObjArray*>(entry->GetObject());
    if (array && array->GetEntriesFast()==6){
      fMapping = new AliTPCAltroMapping*[6];
      for (Int_t i=0; i<6; i++){
	fMapping[i] =  dynamic_cast<AliTPCAltroMapping*>(array->At(i));
      }
    }
  }



  //entry          = GetCDBEntry("TPC/Calib/ExB");
  //if (entry) {
  //  entry->SetOwner(kTRUE);
  //  fExB=dynamic_cast<AliTPCExB*>(entry->GetObject()->Clone());
  //}
  //
  // ExB  - calculate during initialization 
  //      - 
  fExB =  GetExB(-5,kTRUE);
     //
  if (!fTransform) {
    fTransform=new AliTPCTransform(); 
  }

  //
  AliCDBManager::Instance()->SetCacheFlag(cdbCache); // reset original CDB cache
  
}



void AliTPCcalibDB::CreateObjectList(const Char_t *filename, TObjArray *calibObjects)
{
//
// Create calibration objects and read contents from OCDB
//
   if ( calibObjects == 0x0 ) return;
   ifstream in;
   in.open(filename);
   if ( !in.is_open() ){
      fprintf(stderr,"Error: cannot open list file '%s'", filename);
      return;
   }
   
   AliTPCCalPad *calPad=0x0;
   
   TString sFile;
   sFile.ReadFile(in);
   in.close();
   
   TObjArray *arrFileLine = sFile.Tokenize("\n");
   
   TIter nextLine(arrFileLine);
   
   TObjString *sObjLine=0x0;
   while ( (sObjLine = (TObjString*)nextLine()) ){
      TString sLine(sObjLine->GetString());
      
      TObjArray *arrNextCol = sLine.Tokenize("\t");
      
      TObjString *sObjType     = (TObjString*)(arrNextCol->At(0));
      TObjString *sObjFileName = (TObjString*)(arrNextCol->At(1));
      
      if ( !sObjType || ! sObjFileName ) continue;
      TString sType(sObjType->GetString());
      TString sFileName(sObjFileName->GetString());
      printf("%s\t%s\n",sType.Data(),sFileName.Data());
      
      TFile *fIn = TFile::Open(sFileName);
      if ( !fIn ){
         fprintf(stderr,"File not found: '%s'", sFileName.Data());
         continue;
      }
      
      if ( sType == "CE" ){
         AliTPCCalibCE *ce = (AliTPCCalibCE*)fIn->Get("AliTPCCalibCE");
         
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadT0());         
         calPad->SetNameTitle("CETmean","CETmean");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadQ());         
         calPad->SetNameTitle("CEQmean","CEQmean");
         calibObjects->Add(calPad);        
         
         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadRMS());
         calPad->SetNameTitle("CETrms","CETrms");
         calibObjects->Add(calPad);         
                  
      } else if ( sType == "Pulser") {
         AliTPCCalibPulser *sig = (AliTPCCalibPulser*)fIn->Get("AliTPCCalibPulser");
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadT0());         
         calPad->SetNameTitle("PulserTmean","PulserTmean");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadQ());         
         calPad->SetNameTitle("PulserQmean","PulserQmean");
         calibObjects->Add(calPad);        
         
         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadRMS());
         calPad->SetNameTitle("PulserTrms","PulserTrms");
         calibObjects->Add(calPad);         
      
      } else if ( sType == "Pedestals") {
         AliTPCCalibPedestal *ped = (AliTPCCalibPedestal*)fIn->Get("AliTPCCalibPedestal");
         
         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadPedestal());         
         calPad->SetNameTitle("Pedestals","Pedestals");
         calibObjects->Add(calPad);
         
         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadRMS());         
         calPad->SetNameTitle("Noise","Noise");
         calibObjects->Add(calPad);        
     
      } else {
         fprintf(stderr,"Undefined Type: '%s'",sType.Data());
         
      }
      delete fIn;
   }
}



void AliTPCcalibDB::MakeTree(const char * fileName, TObjArray * array, const char * mapFileName, AliTPCCalPad* outlierPad, Float_t ltmFraction) {
  //
  // Write a tree with all available information
  // if mapFileName is specified, the Map information are also written to the tree
  // pads specified in outlierPad are not used for calculating statistics
  //  - the same function as AliTPCCalPad::MakeTree - 
  //
   AliTPCROC* tpcROCinstance = AliTPCROC::Instance();

   TObjArray* mapIROCs = 0;
   TObjArray* mapOROCs = 0;
   TVectorF *mapIROCArray = 0;
   TVectorF *mapOROCArray = 0;
   Int_t mapEntries = 0;
   TString* mapNames = 0;
   
   if (mapFileName) {
      TFile mapFile(mapFileName, "read");
      
      TList* listOfROCs = mapFile.GetListOfKeys();
      mapEntries = listOfROCs->GetEntries()/2;
      mapIROCs = new TObjArray(mapEntries*2);
      mapOROCs = new TObjArray(mapEntries*2);
      mapIROCArray = new TVectorF[mapEntries];
      mapOROCArray = new TVectorF[mapEntries];
      
      mapNames = new TString[mapEntries];
      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
	TString nameROC(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
         nameROC.Remove(nameROC.Length()-4, 4);
         mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((nameROC + "IROC").Data()), ivalue);
         mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((nameROC + "OROC").Data()), ivalue);
         mapNames[ivalue].Append(nameROC);
      }
      
      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
         mapIROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(0));
         mapOROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(36));
      
         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(0); ichannel++)
            (mapIROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapIROCs->At(ivalue)))->GetValue(ichannel);
         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(36); ichannel++)
            (mapOROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapOROCs->At(ivalue)))->GetValue(ichannel);
      }

   } //  if (mapFileName)
  
   TTreeSRedirector cstream(fileName);
   Int_t arrayEntries = array->GetEntries();
   
   TString* names = new TString[arrayEntries];
   for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
      names[ivalue].Append(((AliTPCCalPad*)array->At(ivalue))->GetName());

   for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++) {
      //
      // get statistic for given sector
      //
      TVectorF median(arrayEntries);
      TVectorF mean(arrayEntries);
      TVectorF rms(arrayEntries);
      TVectorF ltm(arrayEntries);
      TVectorF ltmrms(arrayEntries);
      TVectorF medianWithOut(arrayEntries);
      TVectorF meanWithOut(arrayEntries);
      TVectorF rmsWithOut(arrayEntries);
      TVectorF ltmWithOut(arrayEntries);
      TVectorF ltmrmsWithOut(arrayEntries);
      
      TVectorF *vectorArray = new TVectorF[arrayEntries];
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
         vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));
      
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
         AliTPCCalROC* calROC = calPad->GetCalROC(isector);
         AliTPCCalROC* outlierROC = 0;
         if (outlierPad) outlierROC = outlierPad->GetCalROC(isector);
         if (calROC) {
            median[ivalue] = calROC->GetMedian();
            mean[ivalue] = calROC->GetMean();
            rms[ivalue] = calROC->GetRMS();
            Double_t ltmrmsValue = 0;
            ltm[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction);
            ltmrms[ivalue] = ltmrmsValue;
            if (outlierROC) {
               medianWithOut[ivalue] = calROC->GetMedian(outlierROC);
               meanWithOut[ivalue] = calROC->GetMean(outlierROC);
               rmsWithOut[ivalue] = calROC->GetRMS(outlierROC);
               ltmrmsValue = 0;
               ltmWithOut[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction, outlierROC);
               ltmrmsWithOut[ivalue] = ltmrmsValue;
            }
         }
         else {
            median[ivalue] = 0.;
            mean[ivalue] = 0.;
            rms[ivalue] = 0.;
            ltm[ivalue] = 0.;
            ltmrms[ivalue] = 0.;
            medianWithOut[ivalue] = 0.;
            meanWithOut[ivalue] = 0.;
            rmsWithOut[ivalue] = 0.;
            ltmWithOut[ivalue] = 0.;
            ltmrmsWithOut[ivalue] = 0.;
         }
      }
      
      //
      // fill vectors of variable per pad
      //
      TVectorF *posArray = new TVectorF[8];
      for (Int_t ivalue = 0; ivalue < 8; ivalue++)
         posArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));

      Float_t posG[3] = {0};
      Float_t posL[3] = {0};
      Int_t ichannel = 0;
      for (UInt_t irow = 0; irow < tpcROCinstance->GetNRows(isector); irow++) {
         for (UInt_t ipad = 0; ipad < tpcROCinstance->GetNPads(isector, irow); ipad++) {
            tpcROCinstance->GetPositionLocal(isector, irow, ipad, posL);
            tpcROCinstance->GetPositionGlobal(isector, irow, ipad, posG);
            posArray[0][ichannel] = irow;
            posArray[1][ichannel] = ipad;
            posArray[2][ichannel] = posL[0];
            posArray[3][ichannel] = posL[1];
            posArray[4][ichannel] = posG[0];
            posArray[5][ichannel] = posG[1];
            posArray[6][ichannel] = (Int_t)(ipad - (Double_t)(tpcROCinstance->GetNPads(isector, irow))/2);
            posArray[7][ichannel] = ichannel;
            
            // loop over array containing AliTPCCalPads
            for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
               AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
               AliTPCCalROC* calROC = calPad->GetCalROC(isector);
               if (calROC)
                  (vectorArray[ivalue])[ichannel] = calROC->GetValue(irow, ipad);
               else
                  (vectorArray[ivalue])[ichannel] = 0;
            }
            ichannel++;
         }
      }
      
      cstream << "calPads" <<
         "sector=" << isector;
      
      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + "_Median=").Data()) << median[ivalue] <<
            (Char_t*)((names[ivalue] + "_Mean=").Data()) << mean[ivalue] <<
            (Char_t*)((names[ivalue] + "_RMS=").Data()) << rms[ivalue] <<
            (Char_t*)((names[ivalue] + "_LTM=").Data()) << ltm[ivalue] <<
            (Char_t*)((names[ivalue] + "_RMS_LTM=").Data()) << ltmrms[ivalue];
         if (outlierPad) {
            cstream << "calPads" <<
               (Char_t*)((names[ivalue] + "_Median_OutlierCutted=").Data()) << medianWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_Mean_OutlierCutted=").Data()) << meanWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_RMS_OutlierCutted=").Data()) << rmsWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_LTM_OutlierCutted=").Data()) << ltmWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_RMS_LTM_OutlierCutted=").Data()) << ltmrmsWithOut[ivalue];
         }
      }

      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + ".=").Data()) << &vectorArray[ivalue];
      }

      if (mapFileName) {
         for  (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
            if (isector < 36)
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapIROCArray[ivalue];
            else
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapOROCArray[ivalue];
         }
      }

      cstream << "calPads" <<
         "row.=" << &posArray[0] <<
         "pad.=" << &posArray[1] <<
         "lx.=" << &posArray[2] <<
         "ly.=" << &posArray[3] <<
         "gx.=" << &posArray[4] <<
         "gy.=" << &posArray[5] <<
         "rpad.=" << &posArray[6] <<
         "channel.=" << &posArray[7];
         
      cstream << "calPads" <<
         "\n";

      delete[] posArray;
      delete[] vectorArray;
   }
   

   delete[] names;
   if (mapFileName) {
      delete mapIROCs;
      delete mapOROCs;
      delete[] mapIROCArray;
      delete[] mapOROCArray;
      delete[] mapNames;
   }
}



void AliTPCcalibDB::RegisterExB(Int_t index, Float_t bz, Bool_t bdelete){
  //
  // Register static ExB correction map
  // index - registration index - used for visualization
  // bz    - bz field in kGaus

  Float_t factor =  bz/(-5.);  // default b filed in Cheb with minus sign
  
  AliMagF*   bmap = new AliMagWrapCheb("MapsExB","MapsExB", 2, factor, 10., AliMagWrapCheb::k5kG,kTRUE,"$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  
  AliTPCExBFirst *exb  = new  AliTPCExBFirst(bmap,0.88*2.6400e+04,50,50,50);
  AliTPCExB::SetInstance(exb);
  
  if (bdelete){
    delete bmap;
  }else{
    AliTPCExB::RegisterField(index,bmap);
  }
  if (index>=fgExBArray.GetEntries()) fgExBArray.Expand((index+1)*2+11);
  fgExBArray.AddAt(exb,index);
}


AliTPCExB*    AliTPCcalibDB::GetExB(Float_t bz, Bool_t deleteB) {
  //
  // bz filed in KGaus not in tesla
  // Get ExB correction map
  // if doesn't exist - create it
  //
  Int_t index = TMath::Nint(5+bz);
  if (index>fgExBArray.GetEntries()) fgExBArray.Expand((index+1)*2+11);
  if (!fgExBArray.At(index)) AliTPCcalibDB::RegisterExB(index,bz,deleteB);
  return (AliTPCExB*)fgExBArray.At(index);
}


void  AliTPCcalibDB::SetExBField(Float_t bz){
  //
  // Set magnetic filed for ExB correction
  //
  fExB = GetExB(bz,kFALSE);
}



void AliTPCcalibDB::GetRunInformations( Int_t run){
  //
  // - > Don't use it for reconstruction - Only for Calibration studies
  //
  AliCDBEntry * entry = 0;
  if (run>= fRunList.GetSize()){
    fRunList.Set(run*2+1);
    fGRPArray.Expand(run*2+1);
    fGRPMaps.Expand(run*2+1);
    fGoofieArray.Expand(run*2+1);
    fVoltageArray.Expand(run*2+1); 
    fTemperatureArray.Expand(run*2+1);
    fVdriftArray.Expand(run*2+1);
  }
  if (fRunList[run]>0) return;
  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data",run);
  if (entry)  {
    AliGRPObject * grpRun = dynamic_cast<AliGRPObject*>(entry->GetObject());
    if (!grpRun){
      TMap* map = dynamic_cast<TMap*>(entry->GetObject());
      if (map){
	//grpRun = new AliGRPObject; 
	//grpRun->ReadValuesFromMap(map);
	grpRun =  MakeGRPObjectFromMap(map);

	fGRPMaps.AddAt(map,run);
      }
    }
    fGRPArray.AddAt(grpRun,run);
  }
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Goofie",run);
  if (entry)  fGoofieArray.AddAt(entry->GetObject(),run);
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/HighVoltage",run);
  if (entry)  fVoltageArray.AddAt(entry->GetObject(),run);
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Temperature",run);
  if (entry)  fTemperatureArray.AddAt(entry->GetObject(),run);
  fRunList[run]=1;  // sign as used

  AliDCSSensor * press = GetPressureSensor(run);
  AliTPCSensorTempArray * temp = GetTemperatureSensor(run);
  if (press && temp){
    AliTPCCalibVdrift * vdrift = new AliTPCCalibVdrift(temp, press,0);
    fVdriftArray.AddAt(vdrift,run);
  }
}


Float_t AliTPCcalibDB::GetGain(Int_t sector, Int_t row, Int_t pad){
  //
  //
  AliTPCCalPad *calPad = Instance()->fDedxGainFactor;;
  if (!calPad) return 0;
  return calPad->GetCalROC(sector)->GetValue(row,pad);
}


AliGRPObject *AliTPCcalibDB::GetGRP(Int_t run){
  //
  // Get GRP object for given run 
  //
  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>((Instance()->fGRPArray).At(run));
  if (!grpRun) {
    Instance()->GetRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(Instance()->fGRPArray.At(run));
    if (!grpRun) return 0; 
  }
  return grpRun;
}

TMap *  AliTPCcalibDB::GetGRPMap(Int_t run){
  //
  //
  //
  TMap * grpRun = dynamic_cast<TMap *>((Instance()->fGRPMaps).At(run));
  if (!grpRun) {
    Instance()->GetRunInformations(run);
    grpRun = dynamic_cast<TMap *>(Instance()->fGRPMaps.At(run));
    if (!grpRun) return 0; 
  }
  return grpRun;
}


AliDCSSensor * AliTPCcalibDB::GetPressureSensor(Int_t run, Int_t type){
  //
  // Get Pressure sensor
  //
  //
  // First try to get if trom map - if existing  (Old format of data storing)
  //
  TMap *map = GetGRPMap(run);  
  if (map){
    AliDCSSensor * sensor = 0;
    TObject *osensor=0;
    if (type==0) osensor = ((*map)("fCavernPressure"));
    if (type==1) osensor = ((*map)("fP2Pressure"));
    sensor =dynamic_cast<AliDCSSensor *>(osensor); 
    if (sensor) return sensor;
  }
  //
  // If not map try to get it from the GRPObject
  //
  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.At(run)); 
  if (!grpRun) {
    GetRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.At(run));
    if (!grpRun) return 0; 
  }
  AliDCSSensor * sensor = grpRun->GetCavernAtmosPressure();
  if (type==1) sensor = grpRun->GetSurfaceAtmosPressure();
  return sensor; 
}

AliTPCSensorTempArray * AliTPCcalibDB::GetTemperatureSensor(Int_t run){
  //
  // Get temperature sensor array
  //
  AliTPCSensorTempArray * tempArray = (AliTPCSensorTempArray *)fTemperatureArray.At(run);
  if (!tempArray) {
    GetRunInformations(run);
    tempArray = (AliTPCSensorTempArray *)fTemperatureArray.At(run);
  }
  return tempArray;
}

AliDCSSensorArray * AliTPCcalibDB::GetGoofieSensors(Int_t run){
  //
  // Get temperature sensor array
  //
  AliDCSSensorArray * goofieArray = (AliDCSSensorArray *)fGoofieArray.At(run);
  if (!goofieArray) {
    GetRunInformations(run);
    goofieArray = (AliDCSSensorArray *)fGoofieArray.At(run);
  }
  return goofieArray;
}

AliDCSSensorArray * AliTPCcalibDB::GetVoltageSensors(Int_t run){
  //
  // Get temperature sensor array
  //
  AliDCSSensorArray * voltageArray = (AliDCSSensorArray *)fVoltageArray.At(run);
  if (!voltageArray) {
    GetRunInformations(run);
    voltageArray = (AliDCSSensorArray *)fVoltageArray.At(run);
  }
  return voltageArray;
}

AliTPCCalibVdrift *     AliTPCcalibDB::GetVdrift(Int_t run){
  //
  // Get the interface to the the vdrift 
  //
  AliTPCCalibVdrift  * vdrift = (AliTPCCalibVdrift*)fVdriftArray.At(run);
  if (!vdrift) {
    GetRunInformations(run);
    vdrift= (AliTPCCalibVdrift*)fVdriftArray.At(run);
  }
  return vdrift;
}


Float_t AliTPCcalibDB::GetChamberHighVoltage(Int_t timeStamp, Int_t run, Int_t sector) {
  //
  // return the chamber HV for given run and time: 0-35 IROC, 36-72 OROC
  //
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray) return 0;
  AliDCSSensor *sensor = voltageArray->GetSensor((sector+1)*3);
  if (!sensor) return 0;
  return sensor->GetValue(stamp);
}

Float_t AliTPCcalibDB::GetPressure(Int_t timeStamp, Int_t run, Int_t type){
  //
  // GetPressure for given time stamp and runt
  //
  TTimeStamp stamp(timeStamp);
  AliDCSSensor * sensor = Instance()->GetPressureSensor(run,type);
  if (!sensor) return 0;
  return sensor->GetValue(stamp);
}

Float_t AliTPCcalibDB::GetValueGoofie(Int_t timeStamp, Int_t run, Int_t type){
  //
  // GetPressure for given time stamp and runt
  //
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(run);
  if (!goofieArray) return 0;
  AliDCSSensor *sensor = goofieArray->GetSensor(type);
  return sensor->GetValue(stamp);
}






Bool_t  AliTPCcalibDB::GetTemperatureFit(Int_t timeStamp, Int_t run, Int_t side,TVectorD& fit){
  //
  //
  //
  TTimeStamp tstamp(timeStamp);
  AliTPCSensorTempArray* tempArray  = Instance()->GetTemperatureSensor(run);
  if (! tempArray) return kFALSE;
  AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
  TLinearFitter * fitter = tempMap->GetLinearFitter(3,side,tstamp);
  if (fitter){
    fitter->Eval(); 
    fitter->GetParameters(fit);
  }
  delete fitter;
  delete tempMap;
  if (!fitter) return kFALSE;
  return kTRUE;
}

Float_t AliTPCcalibDB::GetTemperature(Int_t timeStamp, Int_t run, Int_t side){
  //
  //
  //
  TVectorD vec(5);
  if (side==0) {
    GetTemperatureFit(timeStamp,run,0,vec);
    return vec[0];
  }
  if (side==1){
    GetTemperatureFit(timeStamp,run,0,vec);
    return vec[0];
  }
  return 0;
}


Double_t AliTPCcalibDB::GetPTRelative(UInt_t timeSec, Int_t run, Int_t side){
  //
  // Get relative P/T 
  // time - absolute time
  // run  - run number
  // side - 0 - A side   1-C side
  AliTPCCalibVdrift * vdrift =  Instance()->GetVdrift(run);
  if (!vdrift) return 0;
  return vdrift->GetPTRelative(timeSec,side);
}


void AliTPCcalibDB::ProcessEnv(const char * runList){
  //
  // Example test function  - how to use the environment variables
  // runList  -  ascii file with run numbers
  // output   -  dcsTime.root file with tree 

  ifstream in;
  in.open(runList);
  Int_t irun=0;
  TTreeSRedirector *pcstream = new TTreeSRedirector("dcsTime.root");
  while(in.good()) {
    in >> irun;
    if (irun==0) continue;
    printf("Processing run %d\n",irun);
    AliDCSSensor * sensorPressure = AliTPCcalibDB::Instance()->GetPressureSensor(irun);
    if (!sensorPressure) continue;
    AliTPCSensorTempArray * tempArray = AliTPCcalibDB::Instance()->GetTemperatureSensor(irun);
    AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
    AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(irun);
    //
    Int_t startTime = sensorPressure->GetStartTime();
    Int_t endTime = sensorPressure->GetEndTime();
    Int_t dtime = TMath::Max((endTime-startTime)/20,10*60);
    for (Int_t itime=startTime; itime<endTime; itime+=dtime){
      //
      TTimeStamp tstamp(itime);
      Float_t valuePressure = sensorPressure->GetValue(tstamp);

      TLinearFitter * fitter = 0;
      TVectorD vecTemp[10];
      if (itime<tempArray->GetStartTime().GetSec() || itime>tempArray->GetEndTime().GetSec()){	
      }else{
	for (Int_t itype=0; itype<5; itype++)
	  for (Int_t iside=0; iside<2; iside++){
	    fitter= tempMap->GetLinearFitter(itype,iside,tstamp);
	    if (!fitter) continue;
	    fitter->Eval(); fitter->GetParameters(vecTemp[itype+iside*5]);
	    delete fitter;
	  }
      }
      
      TVectorD vecGoofie, vecEntries, vecMean, vecMedian,vecRMS;
      if (goofieArray){	
	vecGoofie.ResizeTo(goofieArray->NumSensors());
	ProcessGoofie(goofieArray, vecEntries ,vecMedian, vecMean, vecRMS);
  //
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor){
	    vecGoofie[isensor] = gsensor->GetValue(tstamp);
	  }
	}
      }


      //tempMap->GetLinearFitter(0,0,itime);
      (*pcstream)<<"dcs"<<
	"run="<<irun<<
	"time="<<itime<<
	"goofie.="<<&vecGoofie<<
	"goofieE.="<<&vecEntries<<
	"goofieMean.="<<&vecMean<<
	"goofieMedian.="<<&vecMedian<<
	"goofieRMS.="<<&vecRMS<<
	"press="<<valuePressure<<
	"temp00.="<<&vecTemp[0]<<
	"temp10.="<<&vecTemp[1]<<
	"temp20.="<<&vecTemp[2]<<
	"temp30.="<<&vecTemp[3]<<
	"temp40.="<<&vecTemp[4]<<
	"temp01.="<<&vecTemp[5]<<
	"temp11.="<<&vecTemp[6]<<
	"temp21.="<<&vecTemp[7]<<
	"temp31.="<<&vecTemp[8]<<
	"temp41.="<<&vecTemp[9]<<
	"\n";
    }
  }
  delete pcstream;
}


void AliTPCcalibDB::ProcessGoofie( AliDCSSensorArray* goofieArray, TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS){
  /*
    
  1       TPC_ANODE_I_A00_STAT
  2       TPC_DVM_CO2
  3       TPC_DVM_DriftVelocity
  4       TPC_DVM_FCageHV
  5       TPC_DVM_GainFar
  6       TPC_DVM_GainNear
  7       TPC_DVM_N2
  8       TPC_DVM_NumberOfSparks
  9       TPC_DVM_PeakAreaFar
  10      TPC_DVM_PeakAreaNear
  11      TPC_DVM_PeakPosFar
  12      TPC_DVM_PeakPosNear
  13      TPC_DVM_PickupHV
  14      TPC_DVM_Pressure
  15      TPC_DVM_T1_Over_P
  16      TPC_DVM_T2_Over_P
  17      TPC_DVM_T_Over_P
  18      TPC_DVM_TemperatureS1
   */
  //
  //
  //  TVectorD  vecMedian; TVectorD  vecEntries; TVectorD  vecMean; TVectorD  vecRMS;
  Double_t kEpsilon=0.0000000001;
  Double_t kBig=100000000000.;
  Int_t nsensors = goofieArray->NumSensors();
  vecEntries.ResizeTo(nsensors);
  vecMedian.ResizeTo(nsensors);
  vecMean.ResizeTo(nsensors);
  vecRMS.ResizeTo(nsensors);
  TVectorF values;
  for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
    AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
    if (gsensor &&  gsensor->GetGraph()){
      Int_t npoints = gsensor->GetGraph()->GetN();
      // filter zeroes
      values.ResizeTo(npoints);
      Int_t nused =0;
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
	if (TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])>kEpsilon && 
	   TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])<kBig ){
	  values[nused]=gsensor->GetGraph()->GetY()[ipoint];
	  nused++;
	}
      }
      //
      vecEntries[isensor]= nused;      
      if (nused>1){
	vecMedian[isensor] = TMath::Median(nused,values.GetMatrixArray());
	vecMean[isensor]   = TMath::Mean(nused,values.GetMatrixArray());
	vecRMS[isensor]    = TMath::RMS(nused,values.GetMatrixArray());
      }
    }
  }
}



AliGRPObject * AliTPCcalibDB::MakeGRPObjectFromMap(TMap *map){
  //
  // Function to covert old GRP run information from TMap to GRPObject
  //
  //  TMap * map = AliTPCcalibDB::GetGRPMap(52406);
  if (!map) return 0;
  AliDCSSensor * sensor = 0;
  TObject *osensor=0;
  osensor = ((*map)("fP2Pressure"));
  sensor  =dynamic_cast<AliDCSSensor *>(osensor); 
  //
  if (!sensor) return 0;
  //
  AliDCSSensor * sensor2 = new AliDCSSensor(*sensor);
  osensor = ((*map)("fCavernPressure"));
  TGraph * gr = new TGraph(2);
  gr->GetX()[0]= -100000.;
  gr->GetX()[1]= 1000000.;
  gr->GetY()[0]= atof(osensor->GetName());
  gr->GetY()[1]= atof(osensor->GetName());
  sensor2->SetGraph(gr);
  sensor2->SetFit(0);
  

  AliGRPObject *grpRun = new AliGRPObject; 
  grpRun->ReadValuesFromMap(map);
  grpRun->SetCavernAtmosPressure(sensor2);
  grpRun->SetSurfaceAtmosPressure(sensor);
  return grpRun;
}



