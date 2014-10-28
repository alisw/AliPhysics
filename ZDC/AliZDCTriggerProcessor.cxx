#include "AliZDCTriggerProcessor.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "AliZDCTriggerParameters.h"

/////////////////////////////////////////////////////////////////////
//								   //
//        Class implementing ZDC  trigger processor.		   //
//								   //
/////////////////////////////////////////////////////////////////////

ClassImp(AliZDCTriggerProcessor)

//______________________________________________________________________________________________
AliZDCTriggerProcessor::AliZDCTriggerProcessor() :
  fSignal(0x0),
  fTriggerParam(0x0)
{
  // default constructor
}

//______________________________________________________________________________________________
AliZDCTriggerProcessor::AliZDCTriggerProcessor(Float_t* signal) :
  fSignal(signal),
  fTriggerParam(GetTriggerParamFromOCDB())
{
  // standard constructor I 
  // gets the trigger parameters directly from OCDB
}

//______________________________________________________________________________________________
AliZDCTriggerProcessor::AliZDCTriggerProcessor(Float_t* signal, AliZDCTriggerParameters* ocdbParam) :
  fSignal(signal),
  fTriggerParam(ocdbParam)
{
  // standard constructor II
}

//_____________________________________________________________________________
AliZDCTriggerProcessor &AliZDCTriggerProcessor::operator =(const AliZDCTriggerProcessor &trig)
{
 // Equal operator.
 this->~AliZDCTriggerProcessor();
 new(this) AliZDCTriggerProcessor(trig);
 return *this;  

}

//______________________________________________________________________________________________
AliZDCTriggerProcessor::AliZDCTriggerProcessor(const AliZDCTriggerProcessor& trigg) :
  TObject(),
  fSignal(trigg.fSignal),
  fTriggerParam(trigg.fTriggerParam)
{
  // copy constructor
}

//______________________________________________________________________________________________
AliZDCTriggerProcessor::~AliZDCTriggerProcessor()
{
  // destructor
}

//______________________________________________________________________________________________
UInt_t AliZDCTriggerProcessor::ProcessEvent()
{
  // process ZDC signals in order to determine the trigger output
  UInt_t ctpInput = 0;
  //
  Bool_t mbTriggered = MBTrigger();
  if(mbTriggered == kTRUE) ctpInput = 0x1;
  Bool_t cenTriggered = CentralTrigger();
  if(cenTriggered == kTRUE) ctpInput = 0x1 << 1;
  Bool_t semicenTriggered = SemicentralTrigger();
  if(semicenTriggered == kTRUE) ctpInput = 0x1 << 2;
  Bool_t emdTriggered = EMDTrigger();
  if(emdTriggered == kTRUE) ctpInput = 0x1 << 3;
  
  if((mbTriggered == kTRUE) || (cenTriggered == kTRUE) || 
     (semicenTriggered == kTRUE) || (emdTriggered == kTRUE)){
      return ctpInput;
  }
  else{
    return 0;
  }
}

//______________________________________________________________________________________________
Bool_t AliZDCTriggerProcessor::MBTrigger()
{
  // is the processed event a MB A-A event?
  Float_t mbTrheshold = 0.;
  if(fTriggerParam) mbTrheshold = fTriggerParam->GetADCMBThreshold();
  // check whether ZDC signal > mbTrheshold
  if((fSignal[0]+fSignal[1]+fSignal[3]+fSignal[4]) > mbTrheshold)
       return kTRUE;
  else return kFALSE;
  
}

//______________________________________________________________________________________________
Bool_t AliZDCTriggerProcessor::CentralTrigger()
{
  // is the processed event a central A-A event?
  
  Float_t zemThr = 0;
  Float_t centralWin[2] = {0,0};
  if(fTriggerParam){
    zemThr = fTriggerParam->GetADCZDCCentralityThr();
    const Float_t* cWin = fTriggerParam->GetADCCentralWindow();
    centralWin[0] = cWin[0];
    centralWin[1] = cWin[1];
  }
  // 
  if((fSignal[0]+fSignal[1]+fSignal[3]+fSignal[4]) > centralWin[0]
     && 
     (fSignal[0]+fSignal[1]+fSignal[3]+fSignal[4]) < centralWin[1]
     &&
     fSignal[2] > zemThr)
  	return kTRUE;
  else return kFALSE;
}

//______________________________________________________________________________________________
Bool_t AliZDCTriggerProcessor::SemicentralTrigger()
{
  // is the processed event a semicentral A-A event?
  
  Float_t zemThr = 0;  
  Float_t semicentralWin[2] = {0,0};
  if(fTriggerParam){
    zemThr  =  fTriggerParam->GetADCZDCCentralityThr();
    const Float_t* cWin =   fTriggerParam->GetADCSemicentralWindow();
    semicentralWin[0] = cWin[0];
    semicentralWin[1] = cWin[1];
  }
  //
  if((fSignal[0]+fSignal[1]+fSignal[3]+fSignal[4]) > semicentralWin[0]
     && 
     (fSignal[0]+fSignal[1]+fSignal[3]+fSignal[4]) < semicentralWin[1]
     &&
     fSignal[2] > zemThr)
  	return kTRUE;
  else return kFALSE;
}

//______________________________________________________________________________________________
Bool_t AliZDCTriggerProcessor::EMDTrigger()
{
  // is the processed an EMD event?
  
  Float_t emdWin[4] = {0,0,0,0};
  if(fTriggerParam){
    const Float_t* eWin =   fTriggerParam->GetADCEMDWindow();
    emdWin[0] = eWin[0];
    emdWin[1] = eWin[1];
    emdWin[2] = eWin[2];
    emdWin[3] = eWin[3];
  }
  // check whether ZNA AND ZNC signals fall into the 
  // 2 distinct windows defined for EMD trigger
  if(fSignal[0] > emdWin[0] && fSignal[0] < emdWin[1]
     &&
     fSignal[3] > emdWin[2] && fSignal[3] < emdWin[3])
  	return kTRUE;
  else return kFALSE;
}

//______________________________________________________________________________________________
AliZDCTriggerParameters* AliZDCTriggerProcessor::GetTriggerParamFromOCDB() const
{
  // retrieving trigger parameter configuration form OCDB
  AliZDCTriggerParameters *trigParam = 0x0;
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Trigger/");
  if(!entry) AliFatal("No calibration data loaded!");  
  else{
    trigParam = dynamic_cast<AliZDCTriggerParameters*>  (entry->GetObject());
    if(!trigParam)  AliFatal("Wrong calibration object in calibration  file!");
  }

  return trigParam;
}
