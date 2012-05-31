#ifndef AliZDCTriggerParameter_H
#define AliZDCTriggerParameter_H

/////////////////////////////////////////////////////////////////
//							       //
//    Class containing the parameters that are configured      //
//    to trigger events with the ZDC (in A-A collisions)       //
//    Use: store the set of parameters needed to calculate     //
//    the trigger function sent to the CTP                     //
//							       //
//    Author: Chiara.Oppedisano@to.infn.it                     //
//							       //
/////////////////////////////////////////////////////////////////                                                             

#include <TObject.h>

class AliZDCTriggerParameters : public TObject{

public:
  AliZDCTriggerParameters();
  AliZDCTriggerParameters(Float_t *adcParam, Float_t *discParam);
  AliZDCTriggerParameters(const AliZDCTriggerParameters &parameters);
  AliZDCTriggerParameters& operator= (const AliZDCTriggerParameters &param);
  virtual ~AliZDCTriggerParameters() {;}
  
  Float_t   GetADCZDCCentralityThr()  const {return fADCZEMCentralityThr;}
  Float_t   GetADCMBThreshold()       const {return fADCMBThreshold;}
  const Float_t * GetADCCentralWindow()     const {return fADCCentralWindow;}
  const Float_t * GetADCSemicentralWindow() const {return fADCSemicentralWindow;}
  const Float_t * GetADCEMDWindow() 	    const {return fADCEMDWindow;}
  //
  Float_t   GetDiscZDCCentralityThr()  const {return fDiscZEMCentralityThr;}
  Float_t   GetDiscMBThreshold()       const {return fDiscMBThreshold;}
  const Float_t * GetDiscCentralWindow()     const {return fDiscCentralWindow;}
  const Float_t * GetDiscEMDWindow() 	     const {return fDiscEMDWindow;}
  const Float_t * GetDiscSemicentralWindow() const {return fDiscSemicentralWindow;}
  
  void SetADCZEMCentralityThr(Float_t zemThrVal) {fADCZEMCentralityThr = zemThrVal;} 
  void SetADCMBThreshold(Float_t mbThrVal) {fADCMBThreshold = mbThrVal;} 
  void SetADCCentralWindow(Float_t* cenThrWin) 
  	{for(int i=0; i<2; i++) fADCCentralWindow[i] = cenThrWin[i];}
  void SetADCSemicentralWindow(Float_t* semicenThrWin)
  	{for(int i=0; i<2; i++) fADCSemicentralWindow[i] = semicenThrWin[i];}
  void SetADCEMDWindow(Float_t* emdWin)
  	{for(int i=0; i<2; i++) fADCEMDWindow[i] = emdWin[i];}
  //
  void SetDiscZEMCentralityThr(Float_t zemThrVal) {fDiscZEMCentralityThr = zemThrVal;} 
  void SetDiscMBThreshold(Float_t mbThrVal) {fDiscMBThreshold = mbThrVal;} 
  void SetDiscCentralWindow(Float_t* cenThrWin) 
  	{for(int i=0; i<2; i++) fDiscCentralWindow[i] = cenThrWin[i];}
  void SetDiscSemicentralWindow(Float_t* semicenThrWin)
  	{for(int i=0; i<2; i++) fDiscSemicentralWindow[i] = semicenThrWin[i];}
  void SetDiscEMDWindow(Float_t* emdWin)
  	{for(int i=0; i<2; i++) fDiscEMDWindow[i] = emdWin[i];}
  
protected:
  // --- Configurable parameters 
  // -> [1] values in ADC channels
  Float_t fADCZEMCentralityThr;     //ZEM ADC value for centrality selection
  Float_t fADCMBThreshold;	    //ZDC ADC value to trigger MB A-A events
  Float_t fADCCentralWindow[2];	    //ZDC ADC value to trigger central A-A events
  Float_t fADCSemicentralWindow[2]; //ZDC ADC value to trigger semicentral A-A events
  Float_t fADCEMDWindow[4];	    //ZDC ADC value to trigger EMD events
  //
  // -> [2] values in discriminator thresholds
  Float_t fDiscZEMCentralityThr;    //ZEM threshold for centrality selection
  Float_t fDiscMBThreshold;	    //ZDC threshold to trigger MB A-A events    
  Float_t fDiscCentralWindow[2];    //ZDC threshold to trigger central A-A events
  Float_t fDiscSemicentralWindow[2];//ZDC threshold to trigger semicentral A-A events
  Float_t fDiscEMDWindow[4];	    //ZDC threshold to trigger EMD events

  ClassDef(AliZDCTriggerParameters, 1)

};

#endif
