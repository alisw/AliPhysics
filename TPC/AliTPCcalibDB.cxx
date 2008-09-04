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
  fClusterParam(0)
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
  fClusterParam(0)

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
  
  AliMagF*   bmap = new AliMagWrapCheb("Maps","Maps", 2, factor, 10., AliMagWrapCheb::k5kG,kTRUE,"$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  
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
