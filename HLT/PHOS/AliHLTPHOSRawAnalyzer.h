#ifndef ALIHLTPHOSRAWANALYZER_H
#define ALIHLTPHOSRAWANALYZER_H

#include "AliHLTProcessor.h"


class AliHLTPHOSRawAnalyzer: public AliHLTProcessor
{
 public:
  AliHLTPHOSRawAnalyzer();
  ~AliHLTPHOSRawAnalyzer();
  AliHLTPHOSRawAnalyzer(const AliHLTPHOSRawAnalyzer & );
  AliHLTPHOSRawAnalyzer & operator = (const AliHLTPHOSRawAnalyzer)
    {
      return *this;
    };
  
  virtual int Deinit(){return 0;};
  virtual int DoDeinit(){return 0;}
  

  virtual const char* GetComponentID(){return 0;};
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&){};


  //  virtual AliHLTComponentDataType GetOutputDataType(){return 0;};
 virtual AliHLTComponentDataType GetOutputDataType(){};

 virtual void GetOutputDataSize(long unsigned int&, double&){};
  virtual void GetOutputDataSize(long  int&, double&){};


  virtual AliHLTComponent* Spawn(){return 0;};
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&){return 0;}


  void BaselineCorrection(double *dataPtr, int N);
  void BaselineCorrection(double *dataPtr, double baselineValue);  
  float GetTiming();
  float GetEnergy();
  void SetData(double *data);

  //  virtual void Analyze() const = 0;

 protected:
  double    *fFloatDataPtr;    /**<Float representation of data that should be fitted */
  double     fSampleFrequency; /**<The ADC sample frequency in MHz used under data taking */
  double     fTau;	       /**<The risetime in micro seconds*/		 
  double     fDTof;            /**<Time of flight in entities of sample intervals */
  double     fDAmpl;           /**<Amplitude in entities of ADC levels*/
  int        n;

};


#endif
