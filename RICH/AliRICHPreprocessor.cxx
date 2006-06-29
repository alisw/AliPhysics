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



const char *AliRICHPreprocessor::fgkAliasName[AliRICHPreprocessor::fgkNalias]={"HMP_DET/HMP_MP0/HMP_MP0_LIQ_LOOP.actual.sensors.Rad1In_Temp",
                                                                               "HMP_DET/HMP_MP0/HMP_MP0_LIQ_LOOP.actual.sensors.Rad1Out_Temp"};




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHPreprocessor::Initialize(Int_t run, UInt_t startTime,UInt_t endTime)
{
  AliPreprocessor::Initialize(run, startTime, endTime);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliRICHPreprocessor::Process(TMap* pDcsMap)
{
//  
  TList* list = GetFileSources(kDAQ, "MAP"); //first analyse a set of pedestal files
  if (list){
    Log("The following sources produced files with the id MAP");
    list->Print();
    delete list;
  }

  if(!pDcsMap)  return 0;                    //no DCS map provided 

  TObjArray result; result.SetOwner(kTRUE);  //result is a array of TF1
    
  
  for(Int_t iAlias=0;iAlias<fgkNalias;iAlias++){//aliases loop
    TObjArray *pOA=(TObjArray*)pDcsMap->GetValue(fgkAliasName[iAlias]);
    if(!pOA) continue;                                                    //no data points for this alias
    TF1 *pF1=new TF1("t11","[0]+[1]*x+[2]*sin([3]*x)",0,10);
    
    TGraph *pGr=new TGraph; pGr->GetXaxis()->SetTimeDisplay(kTRUE);       //tmp graph of sensor data versus time 

    TIter next(pOA);  AliDCSValue *pDcsVal; Int_t i=0;
    while((pDcsVal=(AliDCSValue*)next()))                           //loop over data points for this sensor and fill the graph
      pGr->SetPoint(i++,pDcsVal->GetTimeStamp(),pDcsVal->GetSimpleValue().GetFloat());
    
    
    pGr->Fit(pF1);                                                        //do fit
    delete pGr;
    result.Add(pF1);
  }
  
  AliCDBMetaData metaData; metaData.SetBeamPeriod(0); metaData.SetResponsible("AliRICHPreprocessor"); metaData.SetComment("SIMULATED");

  return Store(&result, &metaData); //use AliPreprocessor::Store(), not allowed to use AliCDBManager directly

}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHPreprocessor::Test()
{
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME/TestCDB"); // initialize location of CDB

//   AliTestShuttle* pShuttle = new AliTestShuttle();   
//   pShuttle->SetDCSInput(SimulateDcsMap());                                           //DCS map format alias->TObjArray of AliDCSValue    
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC1", "map1.root");  //????? real gain map
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC2", "map2.root");  //how to crrespond LDC id and staff from AliRICHDigit ????
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC3", "map3.root");
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC4", "map4.root");
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC5", "map5.root");
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC6", "map6.root");
//   pShuttle->AddInputFile(AliTestShuttle::kDAQ, "RICH", "MAP", "LDC7", "map7.root");
//   
//   AliPreprocessor* pp = new AliRICHPreprocessor(pShuttle);                           //start test, actual invocation of Process will be done from shuttle
//   pShuttle->Process();                                    
//   delete pp;
  
  
  
//read array of TF1 stored in CDB  
  AliCDBEntry *pEntry=AliCDBManager::Instance()->Get("RICH/SHUTTLE/Data",0);
  if(!pEntry) Printf("ERROR file is not retrieved!!!");
}//Test()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TMap* AliRICHPreprocessor::SimulateDcsMap()
{
  TMap*      pDcsMap = new TMap;       pDcsMap->SetOwner(1);
  
  for(Int_t iAlias=0;iAlias<fgkNalias;iAlias++){//loop on aliases
    TObjArray* pOA  = new TObjArray;  pOA->SetOwner(1); //values are supposed to be arranged in TObjArray
    for (Int_t timeStamp=0;timeStamp<1000;timeStamp+=10) {
      AliSimpleValue* pSimVal = new AliSimpleValue(Float_t(20*gRandom->Gaus()));                    //T sensor provides floats
      AliDCSValue*    pDcsVal = new AliDCSValue(*pSimVal, timeStamp);      
      pOA->Add(pDcsVal);                                                                            //add new data point to array
    }
    pDcsMap->Add(new TObjString(fgkAliasName[iAlias]),pOA);                                         //add new array of data points to the map
  }//aliases loop
  return pDcsMap;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
