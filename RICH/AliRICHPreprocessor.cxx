#include "AliRICHPreprocessor.h" //header

#include <AliCDBMetaData.h>
#include <AliDCSValue.h>      
#include <TObjArray.h>        //Test()
#include <TObjString.h>       //Test()
#include <AliCDBManager.h>    //Test()
#include <AliCDBEntry.h>      //Test()
//#include <AliTestShuttle.h>   //Test()
#include <TRandom.h>          //Test()
#include <TF1.h>              //Process()
#include <TGraph.h>           //Process()

ClassImp(AliRICHPreprocessor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHPreprocessor::Initialize(Int_t run, UInt_t startTime,UInt_t endTime)
{
  AliPreprocessor::Initialize(run, startTime, endTime);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliRICHPreprocessor::Process(TMap* pDcsMap)
{
// 
// Argumets: pDcsMap - map of structure "alias name" - TObjArray of AliDCSValue
   
  TList* list = GetFileSources(kDAQ, "MAP"); //first analyse a set of pedestal files
  if (list){
    Log("The following sources produced files with the id MAP");
    list->Print();
    delete list;
  }

  if(!pDcsMap)  return 0;                    //no DCS map provided 

  TObjArray result; result.SetOwner(kTRUE);  //result is a array of TF1
    
  
  for(Int_t iCh=0;iCh<7;iCh++){                             //aliases loop
    for(Int_t iRad=0;iRad<3;iRad++){
      TObjArray *pValLst=(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));//get data serias
      if(!pValLst) continue;                                                                                               //no data serias this alias
      TF1    *pF1=new TF1("t11","[0]+[1]*x+[2]*sin([3]*x)",0,10); pF1->SetLineColor(iRad+2);     //temp=f(time) for fitting data seria 
      TGraph *pGr=new TGraph;                                                  //tmp graph of sensor data versus time 
      TIter next(pValLst);  AliDCSValue *pDcsVal; Int_t i=0;
      while((pDcsVal=(AliDCSValue*)next()))                                    //loop over data points for this sensor 
        pGr->SetPoint(i++,pDcsVal->GetTimeStamp(),pDcsVal->GetFloat());        //and fill the graph
      pGr->Fit(pF1);                                                           //now fit the graph 
      delete pGr;
      result.Add(pF1);
    }//radiators loop
  }//chambers loop
  
  AliCDBMetaData metaData; metaData.SetBeamPeriod(0); metaData.SetResponsible("AliRICHPreprocessor"); metaData.SetComment("SIMULATED");

  return Store("DCS", "RefIdx", &result, &metaData); //use AliPreprocessor::Store(), not allowed to use AliCDBManager directly

}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
