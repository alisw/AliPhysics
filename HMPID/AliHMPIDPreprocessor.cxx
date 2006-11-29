#include "AliHMPIDPreprocessor.h" //header

#include <AliCDBMetaData.h>
#include <AliDCSValue.h>      
#include <TObjArray.h>        //Test()
#include <TObjString.h>       //Test()
#include <AliCDBManager.h>    //Test()
#include <AliCDBEntry.h>      //Test()
//#include <AliTestShuttle.h>   //Test()
#include <TRandom.h>          //Test()
#include <TF1.h>              //Process()
#include <TF2.h>              //Process()
#include <TGraph.h>           //Process()

ClassImp(AliHMPIDPreprocessor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDPreprocessor::Initialize(Int_t run, UInt_t startTime,UInt_t endTime)
{
  AliPreprocessor::Initialize(run, startTime, endTime);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliHMPIDPreprocessor::Process(TMap* pDcsMap)
{
// 
// Argumets: pDcsMap - map of structure "alias name" - TObjArray of AliDCSValue
   
//  TList* list = GetFileSources(kDAQ, "MAP"); //first analyse a set of pedestal files
//  if (list){
//    Log("The following sources produced files with the id MAP");
//    list->Print();
//    delete list;
//  }

  if(!pDcsMap)  return 0;                    //no DCS map provided 

  TF2 idxC6F14("RidxC4F14","sqrt(1+0.554*(1239.84/x)^2/((1239.84/x)^2-5796)-0.0005*(y-20))",5.5,8.5,0,50); //DiMauro mail temp 0-50 degrees C
  Double_t eMean=6.67786;                                                                                  //mean energy of photon defined  by transperancy window
  
  TObjArray radTemp; radTemp.SetOwner(kTRUE);      //store temp versus time as TF1 array for all radiators (21)
  TObjArray meanIdx; meanIdx.SetOwner(kTRUE);      //store ref idx versus time as TF1 array for all radiators (21)
    
  
  for(Int_t iCh=0;iCh<7;iCh++){                             //aliases loop
    for(Int_t iRad=0;iRad<3;iRad++){
      TObjArray *pValLst=(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));//get data points for this alias
      if(!pValLst) continue;                                                                                                                //no data points 
      TF1    *pRadTempF=new TF1(Form("RadTemp%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",0,10); pRadTempF->SetLineColor(iRad+2);             //temp=f(time) 
      TF1    *pMeanIdxF=new TF1(Form("MeanIdx%i%i",iCh,iRad),"[0]+[1]*x+[2]*sin([3]*x)",0,10); pMeanIdxF->SetLineColor(iRad+2);             //idx=f(time) 
      TGraph *pRadTempG=new TGraph;                                                                    //tmp graph of rad temp versus time 
      TGraph *pMeanIdxG=new TGraph;                                                                    //tmp graph of mean ref idx versus time 
      TIter next(pValLst);  AliDCSValue *pDcsVal; Int_t i=0;
      while((pDcsVal=(AliDCSValue*)next())){                                                           //loop over data points for this sensor 
        pRadTempG->SetPoint(i,pDcsVal->GetTimeStamp(),                    pDcsVal->GetFloat());        //and fill the temp graph
        pMeanIdxG->SetPoint(i,pDcsVal->GetTimeStamp(),idxC6F14.Eval(eMean,pDcsVal->GetFloat()));       //and fill the maen ref idx graph
        i++;
      }
      pRadTempG->Fit(pRadTempF,"Q");                                                                   //now fit the temp graph 
      pMeanIdxG->Fit(pMeanIdxF,"Q");                                                                   //now fit the mean idx  graph 
      delete pRadTempG; 
      delete pMeanIdxG;
      radTemp.Add(pRadTempF);
      meanIdx.Add(pMeanIdxF);
    }//radiators loop
  }//chambers loop
  
  AliCDBMetaData metaData; metaData.SetBeamPeriod(0); metaData.SetResponsible("AliHMPIDPreprocessor"); metaData.SetComment("SIMULATED");

  Store("DCS", "RadTemp" , &radTemp , &metaData); //use AliPreprocessor::Store(), not allowed to use AliCDBManager directly
  Store("DCS", "MeanIdx" , &meanIdx , &metaData); 
  
  return 1;

}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
