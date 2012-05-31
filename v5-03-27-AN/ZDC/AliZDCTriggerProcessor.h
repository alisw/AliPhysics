#ifndef ALI_ZDC_TRIGGERPROCESSOR_H
#define ALI_ZDC_TRIGGERPREPROCESSOR_H

//////////////////////////////////////////////////////////////////////
//                                                                  //
// 			ZDC trigger processor                       //
//      It collects data stored by the dedicated DA and writes	    //
// 	an AliZDCTriggerParameters object into the OCDB             //
//                                                                  //
//        Author: Chiara.Oppedisano@to.infn.it                      //
//////////////////////////////////////////////////////////////////////
#include <TObject.h>

class AliZDCTriggerParameters;

class AliZDCTriggerProcessor : public TObject 
{
  public:
    AliZDCTriggerProcessor();
    AliZDCTriggerProcessor(Float_t* signal);
    AliZDCTriggerProcessor(Float_t* signal, AliZDCTriggerParameters* ocdbParam);
    AliZDCTriggerProcessor(const AliZDCTriggerProcessor& trigg);  
    AliZDCTriggerProcessor& operator= (const AliZDCTriggerProcessor &trig);
    virtual ~AliZDCTriggerProcessor();
    
    AliZDCTriggerParameters *GetTriggerParamFromOCDB() const;
    virtual void SetTriggerParam(AliZDCTriggerParameters* ocdbParam) 
    	{fTriggerParam = ocdbParam;}
    
    Float_t* GetSignal() const {return fSignal;}
    Float_t  GetSignal(Int_t idet) const {return fSignal[idet];}
    void  SetSignal(Float_t* signal) 
    	    {for(Int_t i=0; i<6; i++) fSignal[i] = signal[i];}
    void  SetSignal(Int_t idet, Float_t signal) {fSignal[idet] = signal;}

  protected:
    virtual UInt_t ProcessEvent();
    virtual Bool_t MBTrigger();
    virtual Bool_t CentralTrigger();
    virtual Bool_t SemicentralTrigger();
    virtual Bool_t EMDTrigger();

  private:
    Float_t* fSignal;
    AliZDCTriggerParameters *fTriggerParam; // parameters in OCDB

    ClassDef(AliZDCTriggerProcessor, 1);
};

    

#endif
