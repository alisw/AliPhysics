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
  AliPreprocessor::Initialize(run, startTime, endTime);
  
  AliInfo(Form("HMPID started for Run %d \n\tStartTime %s \n\t  EndTime %s", run,TTimeStamp(startTime).AsString(),TTimeStamp(endTime).AsString()));

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliHMPIDPreprocessor::Process(TMap* pMap)
{
// Process all information from DCS and DAQ
// Arguments: pMap- map of DCS aliases
// Returns: 0 on success or 1 on error (opposite to Store!)
  
  Log("HMPID - Process in Preprocessor started");
  if(! pMap) {Log(" - Not map of DCS aliases for HMPID - ");return 1;}   
  
  TString runType = GetRunType();
  Log(Form(" AliHMPIDPreprocessor: RunType is %s",runType.Data()));
  Bool_t result1,result2;
  if (runType == "PEDESTAL_RUN"){
    result1 = ProcPed(); return !result1;
  } else if ( runType == "PHYSICS" ){
    result1 = ProcPed(); 
    result2 = ProcDcs(pMap); return !(result1&&result2);
  } else {
    Log("HMPID - Nothing to do with preprocessor for HMPID, bye!");
  return kFALSE;
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
  TObjArray arPress(7);        arPress.SetOwner(kTRUE);     //7  Press=f(time) one pre chamber
  TObjArray arNmean(21);       arNmean.SetOwner(kTRUE);     //21 Nmean=f(time) one per radiator
  TObjArray arQthre(7);        arQthre.SetOwner(kTRUE);     //7  Qthre=f(time) one pre chamber
  TObjArray arUserCut(7);    arUserCut.SetOwner(kTRUE);     //7  user cut in number of sigmas
  
  AliDCSValue *pVal; Int_t cnt=0;
    
  for(Int_t iCh=0;iCh<7;iCh++){                   
// evaluate High Voltage
    TObjArray *pHV=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC0/HMP_MP%i_SEC0_HV.actual.vMon",iCh,iCh,iCh,iCh)); TIter nextHV(pHV);
    TGraph *pGrHV=new TGraph; cnt=0;
    while((pVal=(AliDCSValue*)nextHV())) pGrHV->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());            //P
    if( cnt!=0) pGrHV->Fit(new TF1(Form("HV%i",iCh),"2000+x*[0]",fStartTime,fEndTime),"Q");                       //clm: if no DCS map entry don't fit
    delete pGrHV;

// evaluate Pressure
    TObjArray *pP =(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWC.actual.value"           ,iCh,iCh,iCh));    TIter nextP(pP);
    TGraph *pGrP=new TGraph; cnt=0; 
    while((pVal=(AliDCSValue*)nextP())) pGrP->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());            //P
    if( cnt!=0) pGrP->Fit(new TF1(Form("P%i",iCh),"1005+x*[0]",fStartTime,fEndTime),"Q");                       //clm: if no DCS map entry don't fit
    delete pGrP;
    
// evaluate Qthre
    arQthre.AddAt(new TF1(Form("HMP_Qthre%i",iCh),Form("3*10^(3.01e-3*HV%i - 4.72)+170745848*exp(-P%i*0.0162012)",iCh),fStartTime,fEndTime),iCh);

    
// evaluate UserCut
    Int_t nSigmaUserCut = 3;
    TObject *pUserCut = new TObject();pUserCut->SetUniqueID(nSigmaUserCut);
    arUserCut.AddAt(pUserCut,iCh);    
    
// evaluate Temperatures    
    for(Int_t iRad=0;iRad<3;iRad++){
      TObjArray *pT1=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));  TIter nextT1(pT1);//Tin
      TObjArray *pT2=(TObjArray*)pMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)); TIter nextT2(pT2);//Tout      
      
      TGraph *pGrT1=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT1())) pGrT1->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T inlet
      if(cnt!=0) pGrT1->Fit(new TF1(Form("Tin%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",fStartTime,fEndTime),"Q");       //fit Tin graph -- clm: if DCS entry
      
      TGraph *pGrT2=new TGraph; cnt=0;  while((pVal=(AliDCSValue*)nextT2())) pGrT2->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat()); //T outlet 
      if(cnt!=0) pGrT2->Fit(new TF1(Form("Tou%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",fStartTime,fEndTime),"Q");       //fit Tout graph -- clm: if DCS entry
            
      delete pGrT1;  delete pGrT2;
	
	    
//      arTmean.Add(pRadTempF);  
// evaluate Mean Refractive Index
      arNmean.AddAt(new TF1(Form("HMP_Nmean%i-%i",iCh,iRad),"1.292",fStartTime,fEndTime),3*iCh+iRad); //Nmean=f(t)
    }//radiators loop
  }//chambers loop
  
  AliCDBMetaData metaData; 
  metaData.SetBeamPeriod(0); 
  metaData.SetResponsible("AliHMPIDPreprocessor"); 
  metaData.SetComment("SIMULATED");

  stDcsStore =   Store("Calib","Qthre",&arQthre,&metaData,0,kTRUE) &&    // from DCS 
                 Store("Calib","Nmean",&arNmean,&metaData,0,kTRUE) &&    // from DCS
                 Store("Calib","UserCut",&arUserCut,&metaData,0,kTRUE);  //really not from DCS...a method ProcManual maybe needed
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
//   Returns: kTRUE on success
  
  Bool_t stPedStore=kFALSE;

  TObjArray aDaqSig(7); aDaqSig.SetOwner(kTRUE); for(Int_t i=0;i<7;i++) aDaqSig.AddAt(new TMatrix(160,144),i); //TObjArray of 7 TMatrixF, m(padx,pady)=sigma
  
  TList *pLdc=GetFileSources(kDAQ,"pedestals"); //get list of LDC names containing id "pedestals"
  Log(Form("HMPID - Pedestal files to be read --> %i LDCs for HMPID",pLdc->GetEntries()));
  
  for(Int_t i=0;i<pLdc->GetEntries();i++)//lists of LDCs
    gSystem->Exec(Form("tar xf %s",GetFile(kDAQ,"pedestals",((TObjString*)pLdc->At(i))->GetName()))); //untar pedestal files from current LDC
  
  AliHMPIDDigit dig;
  AliHMPIDRawStream rs;
  Int_t nSigCut,r,d,a,hard;  Float_t mean,sigma;
  
  for(Int_t ddl=0;ddl<14;ddl++){  
    ifstream infile(Form("HmpidPedDdl%02i.txt",ddl));
    if(!infile.is_open()) {Log("No pedestal file found for HMPID,bye!");return kFALSE;}
    TMatrix *pM=(TMatrixF*)aDaqSig.At(ddl/2);
    infile>>nSigCut; pM->SetUniqueID(nSigCut); //n. of pedestal distribution sigmas used to create zero suppresion table
    while(!infile.eof()){
      infile>>dec>>r>>d>>a>>mean>>sigma>>hex>>hard;      
      //dig.Raw(ddl,r,d,a);
      rs.Raw(ddl,r,d,a);
      AliHMPIDDigit dig(rs.GetPad(ddl,r,d,a),rs.GetCharge(ddl,r,d,a));
      (*pM)(dig.PadChX(),dig.PadChY()) = sigma;
    }
    infile.close();
    Log(Form("Pedestal file for DDL %i read successfully",ddl));
  }
  
  AliCDBMetaData metaData; 
  metaData.SetBeamPeriod(0); 
  metaData.SetResponsible("AliHMPIDPreprocessor"); 
  metaData.SetComment("SIMULATED");
  
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
 
