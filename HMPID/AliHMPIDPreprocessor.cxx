#include "AliHMPIDPreprocessor.h" //header no includes
#include "AliHMPIDDigit.h"        //ProcPed()
#include "AliHMPIDRawStream.h"    //ProcPed()
#include <Riostream.h>            //ProcPed()  
#include <AliLog.h>               //all
#include <AliCDBMetaData.h>       //ProcPed(), ProcDcs()
#include <AliDCSValue.h>          //ProcDcs()
#include <TObjString.h>           //ProcDcs(), ProcPed()
#include <TTimeStamp.h>           //Initialize()
#include <TF1.h>                  //Process()
#include <TF2.h>                  //Process()
#include <TString.h>
#include <TGraph.h>               //Process()
#include <TMatrix.h>              //ProcPed()
#include <TList.h>                //ProcPed()
#include <TSystem.h>              //ProcPed()
//.
// HMPID Preprocessor base class
//.
//.
//.
ClassImp(AliHMPIDPreprocessor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDPreprocessor::Initialize(Int_t run, UInt_t startTime,UInt_t endTime)
{
// Initialize the parameter coming from AliPreprocessor
//  run -> run number
// startTime -> starting time 
// endTime   -> ending time
  AliPreprocessor::Initialize(run, startTime, endTime);
  
  AliInfo(Form("HMPID started for Run %d \n\tStartTime %s \n\t  EndTime %s", run,TTimeStamp(startTime).AsString(),TTimeStamp(endTime).AsString()));

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliHMPIDPreprocessor::Process(TMap* pMap)
{
// Process all information from DCS and DAQ
// Arguments: pMap- map of DCS aliases
// Returns: 0 on success or 1 on error (opposite to Store!)

  TString runType = GetRunType();
  Log(Form(" AliHMPIDPreprocessor: RunType is %s",runType.Data()));
  Bool_t statusNoise=kFALSE, statusDcs=kFALSE;
// start to check event type and procedures
  
  Log("HMPID - Process in Preprocessor started");
  if(! pMap) {
    Log("HMPID - ERROR - Not map of DCS aliases for HMPID - ");             return kTRUE;   // error in the DCS mapped aliases
  }   
  if (runType == "CALIBRATION"){
    if (!ProcPed()){
    	Log("HMPID - ERROR - Pedestal processing failed!!");                return kTRUE;   // error in pedestal processing
    } else {
    	Log("HMPID - Pedestal processing successful!!");                    return kFALSE;  // ok for pedestals
    }
  }//CALIBRATION
  else if ( runType=="STANDALONE" || runType=="PHYSICS"){   
   statusDcs=ProcDcs(pMap);
   statusNoise=ProcNoiseMap();
   if(!statusDcs || !statusNoise) { Log(Form("HMPID - ERROR - Noise Map(%d) and/or DCS(%d) processing failed!! (0=OK, 1=FAILED)",statusNoise,statusDcs)); return kTRUE; }  // error in Noise Map or DCS processing
   else                           { Log("HMPID - Noise Map and DCS processing successful!!");  return kFALSE;}  // ok
  }//STANDALONE or PHYSICS run
   else {
    Log("HMPID - Nothing to do with preprocessor for HMPID, bye!");         return kFALSE;  // ok - nothing done
  }
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDPreprocessor::ProcNoiseMap()
{
  //
  // Goal: Process the Noise Map created by the HMP Physics DA to mask 
  // noisy channels from reconstruction and follow changes in accepatnce
  // eg. DDL turn on/off after PEDESTAL run and between PHYSICS runs.
  // Returns kFALSE on success
 
  Bool_t stProcNoise=kFALSE;
  TFile  *fNoiseFile;
  TH2F   *hNoiseMap = 0x0;
  
  TList *pNoiseSource=GetFileSources(kDAQ,"HmpPhysicsDaNoiseMap.root"); //get list of DAQ source names containing id "HmpPhysicsDaNoiseMap" --> defined in HMPIDphysda.cxx
  if(!pNoiseSource) {Log(Form("ERROR: Retrieval of sources for noise map: HmpPhysicsDaNoiseMap.root is failed!")); return stProcNoise;}
  if(!(TObjString*)pNoiseSource->At(0)) {Log(Form("ERROR: empty list received from DAQ Source!")); return stProcNoise;}
    
  TString noiseFile = GetFile(kDAQ,Form("HmpPhysicsDaNoiseMap.root"),((TObjString*)pNoiseSource->At(0))->GetName());
  if(noiseFile.Length()==0) {Log(Form("ERROR retrieving noise map file: HmpPhysicsDaNoiseMap.root")); return stProcNoise;}
  
  fNoiseFile = TFile::Open(noiseFile.Data(),"read");
  if(!fNoiseFile) {Log(Form("ERROR cannot open NoiseFile: %s!",noiseFile.Data())); return stProcNoise;}
  hNoiseMap = (TH2F*) fNoiseFile->Get("hHmpNoiseMaps");
  
  AliCDBMetaData metaDataHisto;
  metaDataHisto.SetBeamPeriod(0);
  metaDataHisto.SetResponsible("AliHMPIDPreprocessor"); 
  metaDataHisto.SetComment("AliHMPIDPreprocessor stores the Noise Map object as Reference Data.");
  AliInfo("Storing Reference Data");
  stProcNoise = Store("Calib","NoiseMap",hNoiseMap,&metaDataHisto,0,kTRUE);
  if(!stProcNoise) {
    Log("HMPID - failure to store Noise Map data results in OCDB");   
  }
  return stProcNoise;
}//ProcNoiseMap
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDPreprocessor::ProcDcs(TMap* pMap)
{
// Process: 1 (old). inlet and outlet C6F14 temperature, stores TObjArray of 21 TF1, where TF1 is Nmean=f(t), one per radiator
// Process: 1. inlet and outlet C6F14 temperature, stores TObjArray of 42 TF1, where TF1 are Tin and Tout per radiator
//             + one function for the mean energy photon (in total 43).
//          2. CH4 pressure and HV                 stores TObjArray of 7 TF1 where TF1 is thr=f(t), one per chamber
// Arguments: pDcsMap - map of structure "alias name" - TObjArray of AliDCSValue
// Assume that: HV is the same during the run for a given chamber, different chambers might have different HV
//              P=f(t), different for different chambers
// Returns: kTRUE on success  

  Bool_t stDcsStore=kFALSE;

// Qthr=f(HV,P) [V,mBar]  logA0=k*HV+b is taken from p. 64 TDR plot 2.59 for PC32 
//                           A0=f(P) is taken from DiMauro mail
// Qthr is estimated as 3*A0

  TF2 thr("RthrCH4"  ,"3*10^(3.01e-3*x-4.72)+170745848*exp(-y*0.0162012)"             ,2000,3000,900,1200); 
  
  TObjArray arNmean(43);       arNmean.SetOwner(kTRUE);     //42 Tin and Tout one per radiator + 1 for ePhotMean
  TObjArray arQthre(42);       arQthre.SetOwner(kTRUE);     //42 Qthre=f(time) one per sector
  
  AliDCSValue *pVal; Int_t cnt=0;
  
  Double_t xP,yP;

//  TF1 **pTin  = new TF1*[21];
//  TF1 **pTout = new TF1*[21];
  TF1  *pTin[21];
  TF1 *pTout[21];

// evaluate Environment Pressure
  
  TObjArray *pPenv=(TObjArray*)pMap->GetValue("HMP_DET/HMP_ENV/HMP_ENV_PENV.actual.value");
  if(!pPenv) {
    AliWarning(" No Data Points from HMP_ENV_PENV.actual.value!");
    return kFALSE;
  } else {
    Log(Form(" Environment Pressure data              ---> %3i entries",pPenv->GetEntries()));
    if(pPenv->GetEntries()) {
      TIter nextPenv(pPenv);
      TGraph *pGrPenv=new TGraph; cnt=0;
      while((pVal=(AliDCSValue*)nextPenv())) pGrPenv->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());        //P env
      if( cnt==1) {
        pGrPenv->GetPoint(0,xP,yP);
        new TF1("Penv",Form("%f",yP),fStartTime,fEndTime);
      } else {
        pGrPenv->Fit(new TF1("Penv","1000+x*[0]",fStartTime,fEndTime),"Q");
      }
      delete pGrPenv;
    } else {AliWarning(" No Data Points from HMP_ENV_PENV.actual.value!");return kFALSE;}
  }  
// evaluate Pressure
  
  for(Int_t iCh=0;iCh<7;iCh++){                   
    TObjArray *pP =(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWPC.actual.value",iCh,iCh,iCh));
    if(!pP) {
      AliWarning(Form(" No Data Points from HMP_MP%1i_GAS_PMWPC.actual.value!",iCh));
      return kFALSE;
    } else {
        Log(Form(" Pressure for module %i data             ---> %3i entries",iCh,pP->GetEntries()));
      if(pP->GetEntries()) {
        TIter nextP(pP);    
        TGraph *pGrP=new TGraph; cnt=0; 
        while((pVal=(AliDCSValue*)nextP())) pGrP->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());            //P
        if( cnt==1) {
          pGrP->GetPoint(0,xP,yP);
          new TF1(Form("P%i",iCh),Form("%f",yP),fStartTime,fEndTime);
        } else {
          pGrP->Fit(new TF1(Form("P%i",iCh),"[0] + x*[1]",fStartTime,fEndTime),"Q");
        }
        delete pGrP;
      } else {AliWarning(" No Data Points from HMP_MP0-6_GAS_PMWPC.actual.value!");return kFALSE;}
    }
    
// evaluate High Voltage
    
   for(Int_t iSec=0;iSec<6;iSec++){
      TObjArray *pHV=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC%i/HMP_MP%i_SEC%i_HV.actual.vMon",iCh,iCh,iCh,iSec,iCh,iSec));
      if(!pHV) {
        AliWarning(Form(" No Data Points from HMP_MP%1i_SEC%1i_HV.actual.vMon!",iCh,iSec));
        return kFALSE;
      } else {
        Log(Form(" HV for module %i and secto %i data       ---> %3i entries",iCh,iSec,pHV->GetEntries()));
        if(pHV->GetEntries()) {
          TIter nextHV(pHV);
          TGraph *pGrHV=new TGraph; cnt=0;
          while((pVal=(AliDCSValue*)nextHV())) pGrHV->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());            //HV
          if( cnt==1) {
            pGrHV->GetPoint(0,xP,yP);
            new TF1(Form("HV%i_%i",iCh,iSec),Form("%f",yP),fStartTime,fEndTime);               
          } else {
            pGrHV->Fit(new TF1(Form("HV%i_%i",iCh,iSec),"[0]+x*[1]",fStartTime,fEndTime),"Q");               
          }
          delete pGrHV;
        } else {AliWarning(" No Data Points from HMP_MP0-6_SEC0-5_HV.actual.vMon!");return kFALSE;}
      }   
// evaluate Qthre
     
      arQthre.AddAt(new TF1(Form("HMP_QthreC%iS%i",iCh,iSec),
          Form("3*10^(3.01e-3*HV%i_%i - 4.72)+170745848*exp(-(P%i+Penv)*0.0162012)",iCh,iSec,iCh),fStartTime,fEndTime),6*iCh+iSec);

     //arQthre.AddAt(new TF1(Form("HMP_QthreC%iS%i",iCh,iSec),"100",fStartTime,fEndTime),6*iCh+iSec);  
    }
// evaluate Temperatures: in and out of the radiators    
    // T in
    for(Int_t iRad=0;iRad<3;iRad++){
      
      pTin[3*iCh+iRad]  = new TF1(Form("Tin%i%i" ,iCh,iRad),"[0]+[1]*x",fStartTime,fEndTime);
      pTout[3*iCh+iRad] = new TF1(Form("Tout%i%i",iCh,iRad),"[0]+[1]*x",fStartTime,fEndTime);          

      //pTin[3*iCh+iRad]  = new TF1(Form("Tin%i%i" ,iCh,iRad),"21",fStartTime,fEndTime);
      //pTout[3*iCh+iRad] = new TF1(Form("Tout%i%i",iCh,iRad),"22",fStartTime,fEndTime);          
                  
    TObjArray *pT1=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));
      if(!pT1) {
        AliWarning(Form(" No Data Points from HMP_MP%1i_LIQ_LOOP.actual.sensors.Rad%1iIn_Temp!",iCh,iRad));
        return kFALSE;
      } else {
        Log(Form(" Temperatures for module %i inside data  ---> %3i entries",iCh,pT1->GetEntries()));
        if(pT1->GetEntries()) {
          TIter nextT1(pT1);//Tin
          TGraph *pGrT1=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT1())) pGrT1->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T inlet
          if(cnt==1) { 
            pGrT1->GetPoint(0,xP,yP);
            pTin[3*iCh+iRad]->SetParameter(0,yP);
            pTin[3*iCh+iRad]->SetParameter(1,0);
          } else {
            pGrT1->Fit(pTin[3*iCh+iRad],"Q");
          }
          delete pGrT1;
        } else {AliWarning(" No Data Points from HMP_MP0-6_LIQ_LOOP.actual.sensors.Rad0-2In_Temp!");return kFALSE;}
      }
    // T out
      TObjArray *pT2=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)); 
      if(!pT2) {
        AliWarning(Form(" No Data Points from HMP_MP%1i_LIQ_LOOP.actual.sensors.Rad%1iOut_Temp!",iCh,iRad));
        return kFALSE;
      } else {
        Log(Form(" Temperatures for module %i outside data ---> %3i entries",iCh,pT2->GetEntries()));
        if(pT2->GetEntries()) {
          TIter nextT2(pT2);//Tout      
          TGraph *pGrT2=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT2())) pGrT2->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T outlet 
          if(cnt==1) { 
            pGrT2->GetPoint(0,xP,yP);
            pTout[3*iCh+iRad]->SetParameter(0,yP);
            pTout[3*iCh+iRad]->SetParameter(1,0);
          } else {
            pGrT2->Fit(pTout[3*iCh+iRad],"Q");
          }
          delete pGrT2;
        } else {AliWarning(" No Data Points from HMP_MP0-6_LIQ_LOOP.actual.sensors.Rad0-2Out_Temp!");return kFALSE;}
      }
	
// evaluate Mean Refractive Index
      
      arNmean.AddAt(pTin[3*iCh+iRad] ,6*iCh+2*iRad  ); //Tin =f(t)
      arNmean.AddAt(pTout[3*iCh+iRad],6*iCh+2*iRad+1); //Tout=f(t)
      
    }//radiators loop
  }//chambers loop
  
  Double_t eMean = ProcTrans(pMap);
  arNmean.AddAt(new TF1("HMP_PhotEmean",Form("%f",eMean),fStartTime,fEndTime),42); //Photon energy mean
    
  AliCDBMetaData metaData; 
  metaData.SetBeamPeriod(0); 
  metaData.SetResponsible("AliHMPIDPreprocessor"); 
  metaData.SetComment("HMPID preprocessor fills TObjArrays.");

  stDcsStore =   Store("Calib","Qthre",&arQthre,&metaData,0,kTRUE) &&    // from DCS  0,kTRUE generates the file from Run 0 to Run 99999999
                 Store("Calib","Nmean",&arNmean,&metaData,0,kTRUE);      // from DCS
//  stDcsStore =   Store("Calib","Qthre",&arQthre,&metaData) &&    // from DCS 
//                 Store("Calib","Nmean",&arNmean,&metaData);      // from DCS
  if(!stDcsStore) {
    Log("HMPID - failure to store DCS data results in OCDB");    
  }
  
//  arNmean.Delete();
//  arQthre.Delete();

//  for(Int_t i=0;i<21;i++) delete pTin[i]; delete []pTin;
//  for(Int_t i=0;i<21;i++) delete pTout[i]; delete []pTout;
  
  return stDcsStore;
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDPreprocessor::ProcPed()
{
// Process pedestal files and create 7 M(padx,pady)=sigma, one for each chamber
// Arguments:
// Returns: kTRUE on success
  
  Bool_t stPedStore=kFALSE;
  Bool_t stDeadMaskedStore=kFALSE;
  AliHMPIDDigit dig;
  AliHMPIDRawStream rs;
  Int_t nSigCut,r,d,a,hard;  Float_t mean,sigma;
  Int_t  runNumber,ldcId,timeStamp,nEv,nDdlEv,nBadEv;  Char_t tName[10]; 
  Float_t nBadEvPer;
  
  TObjArray aDaqSig(7);     aDaqSig.SetOwner(kTRUE);     for(Int_t i=0;i<7;i++) aDaqSig.AddAt(new TMatrix(160,144),i);         //TObjArray of 7 TMatrixF, m(padx,pady)=sigma
  TObjArray aDeadMasked(7); aDeadMasked.SetOwner(kTRUE); for(Int_t i=0;i<7;i++) aDeadMasked.AddAt(new TMatrix(160,144),i);     //TObjArray of 7 TMatrixF, m(padx,pady)=pedestal
  
  for(Int_t iddl=0;iddl<AliHMPIDRawStream::kNDDL;iddl++)            //retrieve the files from LDCs independently the DDL<->LDC connection
  {
    TList *pLdc=GetFileSources(kDAQ,Form("HmpidPedDdl%02i.txt",iddl)); //get list of LDC names containing id "pedestals"
    if(!pLdc) {Log(Form("ERROR: Retrieval of sources for pedestals: HmpidPedDdl%02i.txt failed!",iddl));continue;}
    
    Log(Form("HMPID - Pedestal files to be read --> %i LDCs for HMPID",pLdc->GetEntries()));
    for(Int_t i=0;i<pLdc->GetEntries();i++) {//lists of LDCs -- but in general we have 1 LDC for 1 ped file
    TString fileName = GetFile(kDAQ,Form("HmpidPedDdl%02i.txt",iddl),((TObjString*)pLdc->At(i))->GetName());
    if(fileName.Length()==0) {Log(Form("ERROR retrieving pedestal file: HmpidPedDdl%02i.txt!",iddl));continue;}
  
    //reading pedestal file
    ifstream infile(fileName.Data()); 
    
    if(!infile.is_open()) {Log("No pedestal file found for HMPID,bye!");continue;}
    TMatrix *pM=(TMatrixF*)aDaqSig.At(iddl/2);
    TMatrix *pDM=(TMatrixF*)aDeadMasked.At(iddl/2);
    infile>>tName>>runNumber;Printf("Xcheck: reading run %i",runNumber);
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv; 
    infile>>tName>>nDdlEv;
    infile>>tName>>nBadEv;
    infile>>tName>>nBadEvPer;
    infile>>tName>>nSigCut; pM->SetUniqueID(nSigCut); //n. of pedestal distribution sigmas used to create zero suppresion table
    while(!infile.eof()){
      infile>>dec>>r>>d>>a>>mean>>sigma>>hex>>hard;
      if(rs.GetPad(iddl,r,d,a)>=0){			//the GetPad returns meaningful abs pad number								
      dig.SetPad(rs.GetPad(iddl,r,d,a));
      dig.SetQ((Int_t)mean);
      (*pM)(dig.PadChX(),dig.PadChY()) = sigma;
      if( (mean == AliHMPIDParam::kPadMeanZeroCharge && sigma == AliHMPIDParam::kPadSigmaZeroCharge) || 
          (mean == AliHMPIDParam::kPadMeanMasked     && sigma == AliHMPIDParam::kPadSigmaMasked)     ) 
        {(*pDM)(dig.PadChX(),dig.PadChY()) = mean;} 
      }
    }
    infile.close();
    Log(Form("Pedestal file for DDL %i read successfully",iddl));
  
 
    
  }//LDCs reading entries

 }//DDL 

  AliCDBMetaData metaData; 
  metaData.SetBeamPeriod(0); 
  metaData.SetResponsible("AliHMPIDPreprocessor"); 
  metaData.SetComment("HMPID processor fills TObjArrays.");  
  stPedStore = Store("Calib","DaqSig",&aDaqSig,&metaData,0,kTRUE);
  if(!stPedStore) {     Log("HMPID - failure to store PEDESTAL data results in OCDB");      }
  stDeadMaskedStore = Store("Calib","Masked",&aDeadMasked,&metaData,0,kTRUE);
  if(!stDeadMaskedStore) {     Log("HMPID - failure to store DEAD & MASKED channel map  in OCDB");      }
  Bool_t pedRes=stPedStore*stDeadMaskedStore;
  return pedRes;
  
}//ProcPed()  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDPreprocessor::ProcTrans(TMap* pMap)
{
  //  Process transparency monitoring data and calculates Emean  
  
  Double_t sEnergProb=0, sProb=0;

  Double_t tRefCR5 = 19. ;                                      // mean temperature of CR5 where the system is in place

  Double_t eMean = 0;
      
  AliDCSValue *pVal;

  for(Int_t i=0; i<30; i++){

    // evaluate wavelenght 
    TObjArray *pWaveLenght = (TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.waveLenght",i));
    if(!pWaveLenght){ 
        AliWarning(Form("No Data Point values for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.waveLenght  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      } 
        
    pVal=(AliDCSValue*)pWaveLenght->At(0);
    Double_t lambda = pVal->GetFloat();

    if(lambda<150. || lambda>230.){ 
        AliWarning(Form("Wrong value for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.waveLenght  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      } 

    Double_t photEn = 1239.842609/lambda;     // 1239.842609 from nm to eV
    
    if(photEn<AliHMPIDParam::EPhotMin() || photEn>AliHMPIDParam::EPhotMax()) continue;
    
    // evaluate phototube current for argon reference
    TObjArray *pArgonRef  = (TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.argonReference",i));
    if(!pArgonRef){ 
        AliWarning(Form("No Data Point values for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.argonReference  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      } 

    pVal=(AliDCSValue*)pArgonRef->At(0);    
    Double_t aRefArgon = pVal->GetFloat();

    // evaluate phototube current for argon cell
    TObjArray *pArgonCell = (TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.argonCell",i));
    if(!pArgonCell){ 
        AliWarning(Form("No Data Point values for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.argonCell  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      } 
      
    pVal=(AliDCSValue*)pArgonRef->At(0);
    Double_t aCellArgon = pVal->GetFloat();

    //evaluate phototube current for freon reference
    TObjArray *pFreonRef  = (TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.c6f14Reference",i));
    if(!pFreonRef){ 
        AliWarning(Form("No Data Point values for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.c6f14Reference  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      } 
        
    pVal=(AliDCSValue*)pFreonRef->At(0);
    Double_t aRefFreon = pVal->GetFloat();

    //evaluate phototube current for freon cell
    TObjArray *pFreonCell = (TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.c6f14Cell",i));
    if(!pFreonCell){
        AliWarning(Form("No Data Point values for HMP_DET/HMP_INFR/HMP_INFR_TRANPLANT/HMP_INFR_TRANPLANT_MEASURE.mesure%i.c6f14Cell  -----> Default E mean used!!!!!",i));
        return DefaultEMean(); // to be checked
      }

    pVal=(AliDCSValue*)pFreonCell->At(0);
    Double_t aCellFreon = pVal->GetFloat();
 
   //evaluate correction factor to calculate trasparency (Ref. NIMA 486 (2002) 590-609)
    
    Double_t aN1 = AliHMPIDParam::NIdxRad(photEn,tRefCR5);
    Double_t aN2 = AliHMPIDParam::NMgF2Idx(photEn);
    Double_t aN3 = 1;                              // Argon Idx

    Double_t aR1               = ((aN1 - aN2)*(aN1 - aN2))/((aN1 + aN2)*(aN1 + aN2));
    Double_t aR2               = ((aN2 - aN3)*(aN2 - aN3))/((aN2 + aN3)*(aN2 + aN3));
    Double_t aT1               = (1 - aR1);
    Double_t aT2               = (1 - aR2);
    Double_t aCorrFactor       = (aT1*aT1)/(aT2*aT2);

    // evaluate 15 mm of thickness C6F14 Trans
    Double_t aTransRad;
    
    Double_t aConvFactor = 1.0 - 0.3/1.8;         
        
    if(aRefFreon*aRefArgon>0) {
      aTransRad  = TMath::Power((aCellFreon/aRefFreon)/(aCellArgon/aRefArgon)*aCorrFactor,aConvFactor);
    } else {
      return DefaultEMean();
    }

    // evaluate 0.5 mm of thickness SiO2 Trans
    Double_t aTransSiO2 = TMath::Exp(-0.5/AliHMPIDParam::LAbsWin(photEn));

    // evaluate 80 cm of thickness Gap (low density CH4) transparency
    Double_t aTransGap  = TMath::Exp(-80./AliHMPIDParam::LAbsGap(photEn));

    // evaluate CsI quantum efficiency
    Double_t aCsIQE            = AliHMPIDParam::QEffCSI(photEn);

    // evaluate total convolution of all material optical properties
    Double_t aTotConvolution   = aTransRad*aTransSiO2*aTransGap*aCsIQE;

    sEnergProb+=aTotConvolution*photEn;  
  
    sProb+=aTotConvolution;  
}
  if(sProb>0) {
    eMean = sEnergProb/sProb;
  } else {
    return DefaultEMean();
  }
  Log(Form(" Mean energy photon calculated ---> %f eV ",eMean));

  if(eMean<AliHMPIDParam::EPhotMin() || eMean>AliHMPIDParam::EPhotMax()) return DefaultEMean();
  
  return eMean;
}   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDPreprocessor::DefaultEMean()
{
    Double_t eMean = 6.675;      //just set a refractive index for C6F14 at ephot=6.675 eV @ T=25 C
    AliWarning(Form("Mean energy for photons out of range [%f,%f] in Preprocessor. Default value Eph=%f eV taken.",AliHMPIDParam::EPhotMin(),
                                                                                                                   AliHMPIDParam::EPhotMax(),
                                                                                                                   eMean));
    return eMean;
}
