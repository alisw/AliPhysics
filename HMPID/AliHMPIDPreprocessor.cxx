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
  } else if ( runType=="STANDALONE" || runType=="PHYSICS"){
    if (!ProcDcs(pMap)){
    	Log("HMPID - ERROR - DCS processing failed!!");                     return kTRUE;   // error in DCS processing
    } else {
    	Log("HMPID - DCS processing successful!!");                         return kFALSE;  // ok for DCS
    }
  } else {
    Log("HMPID - Nothing to do with preprocessor for HMPID, bye!");         return kFALSE;  // ok - nothing done
  }
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDPreprocessor::ProcDcs(TMap* pMap)
{
// Process: 1. inlet and outlet C6F14 temperature, stores TObjArray of 21 TF1, where TF1 is Nmean=f(t), one per radiator
//          2. CH4 pressure and HV                 stores TObjArray of 7 TF1 where TF1 is thr=f(t), one per chamber
// Arguments: pDcsMap - map of structure "alias name" - TObjArray of AliDCSValue
// Assume that: HV is the same during the run for a given chamber, different chambers might have different HV
//              P=f(t), different for different chambers
// Returns: kTRUE on success  

  Bool_t stDcsStore=kFALSE;

  TF2 idx("RidxC4F14","sqrt(1+0.554*(1239.84/x)^2/((1239.84/x)^2-5769)-0.0005*(y-20))",5.5 ,8.5 ,0  ,50);  //N=f(Ephot,T) [eV,grad C] DiMauro mail
  
// Qthr=f(HV,P) [V,mBar]  logA0=k*HV+b is taken from p. 64 TDR plot 2.59 for PC32 
//                           A0=f(P) is taken from DiMauro mail
// Qthr is estimated as 3*A0
  TF2 thr("RthrCH4"  ,"3*10^(3.01e-3*x-4.72)+170745848*exp(-y*0.0162012)"             ,2000,3000,900,1200); 
  
  TObjArray arTmean(21);       arTmean.SetOwner(kTRUE);     //21 Tmean=f(time) one per radiator
  TObjArray arPress(7);        arPress.SetOwner(kTRUE);     //7  Press=f(time) one per chamber
  TObjArray arNmean(21);       arNmean.SetOwner(kTRUE);     //21 Nmean=f(time) one per radiator
  TObjArray arQthre(42);       arQthre.SetOwner(kTRUE);     //42 Qthre=f(time) one per sector
  
  AliDCSValue *pVal; Int_t cnt=0;
  
  Double_t xP,yP;

// evaluate Environment Pressure
  
  TObjArray *pPenv=(TObjArray*)pMap->GetValue("HMP_DET/HMP_ENV/HMP_ENV_PENV.actual.value");
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
    
// evaluate Pressure
  
  for(Int_t iCh=0;iCh<7;iCh++){                   
    TObjArray *pP =(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWPC.actual.value",iCh,iCh,iCh));
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
    
// evaluate High Voltage
    
    for(Int_t iSec=0;iSec<6;iSec++){
      TObjArray *pHV=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC%i/HMP_MP%i_SEC%i_HV.actual.vMon",iCh,iCh,iCh,iSec,iCh,iSec));
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
     
// evaluate Qthre
     
      arQthre.AddAt(new TF1(Form("HMP_QthreC%iS%i",iCh,iSec),
          Form("3*10^(3.01e-3*HV%i_%i - 4.72)+170745848*exp(-(P%i+Penv)*0.0162012)",iCh,iSec,iCh),fStartTime,fEndTime),6*iCh+iSec);
    }
    
// evaluate Temperatures: in and out of the radiators    
    
    // T in
    for(Int_t iRad=0;iRad<3;iRad++){
      TObjArray *pT1=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));  
      Log(Form(" Temperatures for module %i inside data  ---> %3i entries",iCh,pT1->GetEntries()));
      if(pT1->GetEntries()) {
        TIter nextT1(pT1);//Tin
        TGraph *pGrT1=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT1())) pGrT1->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T inlet
        if(cnt==1) { 
          pGrT1->GetPoint(0,xP,yP);
          new TF1(Form("Tin%i%i",iCh,iRad),Form("%f",yP),fStartTime,fEndTime);
        } else {
          pGrT1->Fit(new TF1(Form("Tin%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",fStartTime,fEndTime),"Q");
        }
        delete pGrT1;
      } else {AliWarning(" No Data Points from HMP_MP0-6_LIQ_LOOP.actual.sensors.Rad0-2In_Temp!");return kFALSE;}
    // T out
      TObjArray *pT2=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)); 
      Log(Form(" Temperatures for module %i outside data ---> %3i entries",iCh,pT2->GetEntries()));
      if(pT2->GetEntries()) {
        TIter nextT2(pT2);//Tout      
        TGraph *pGrT2=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT2())) pGrT2->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T outlet 
        if(cnt==1) { 
          pGrT2->GetPoint(0,xP,yP);
          new TF1(Form("Tou%i%i",iCh,iRad),Form("%f",yP),fStartTime,fEndTime);
        } else {
          pGrT2->Fit(new TF1(Form("Tou%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",fStartTime,fEndTime),"Q");
        }
        delete pGrT2;
      } else {AliWarning(" No Data Points from HMP_MP0-6_LIQ_LOOP.actual.sensors.Rad0-2Out_Temp!");return kFALSE;}
	
// evaluate Mean Refractive Index
      
      arNmean.AddAt(new TF1(Form("HMP_Nmean%i-%i",iCh,iRad),"1.292",fStartTime,fEndTime),3*iCh+iRad); //Nmean=f(t)
      
    }//radiators loop
  }//chambers loop
  
  AliCDBMetaData metaData; 
  metaData.SetBeamPeriod(0); 
  metaData.SetResponsible("AliHMPIDPreprocessor"); 
  metaData.SetComment("HMPID preprocessor fills TObjArrays.");

  stDcsStore =   Store("Calib","Qthre",&arQthre,&metaData) &&    // from DCS 
                 Store("Calib","Nmean",&arNmean,&metaData);      // from DCS
  if(!stDcsStore) {
    Log("HMPID - failure to store DCS data results in OCDB");    
  }
  return stDcsStore;
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDPreprocessor::ProcPed()
{
// Process pedestal files and create 7 M(padx,pady)=sigma, one for each chamber
// Arguments:
// Returns: kTRUE on success
  
  Bool_t stPedStore=kFALSE;
  AliHMPIDDigit dig;
  AliHMPIDRawStream rs;
  Int_t nSigCut,r,d,a,hard;  Float_t mean,sigma;
  Int_t  runNumber,ldcId,timeStamp,nEv,nDdlEv,nBadEv;  Char_t tName[10]; 
  Float_t nBadEvPer;
  
  TObjArray aDaqSig(7); aDaqSig.SetOwner(kTRUE); for(Int_t i=0;i<7;i++) aDaqSig.AddAt(new TMatrix(160,144),i); //TObjArray of 7 TMatrixF, m(padx,pady)=sigma
  
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
  if(!stPedStore) {
    Log("HMPID - failure to store PEDESTAL data results in OCDB");    
  }
  return stPedStore;
  
}//ProcPed()  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t ProcTrans()
{
// Process transparency monitoring data and calculates Emean  
  Double_t eMean=6.67786;     //mean energy of photon defined  by transperancy window
  return eMean;
}   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
