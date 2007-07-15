#ifndef ALIHLTPHOSPROCESSOR_H
#define ALIHLTPHOSPROCESSOR_H

#include "AliHLTProcessor.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSCommonDefs.h"
#include "TString.h"
#include "AliHLTPHOSDefinitions.h"
#include "Rtypes.h"
#include <iostream>

using namespace PhosHLTConst;

class AliHLTPHOSProcessor:public AliHLTProcessor
{
 public:
  AliHLTPHOSProcessor();
  virtual ~AliHLTPHOSProcessor();
  virtual int DoInit(int argc, const char** argv) = 0;
  virtual int Deinit() = 0;
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list) =0;
  virtual AliHLTComponentDataType GetOutputDataType() =0;
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) =0;
  virtual AliHLTComponent* Spawn() = 0; 

  template<typename T> 
    void  DumpData(T *array, int N, int nPerLine)
    {
      cout <<   "DumpData N=  " << N <<endl;
      for(int i= 0; i< N; i++)
	{
	  if((i%nPerLine == 0)  &&  (i != 0))
	    {
	      printf("\n");
	    }

	  cout << array[i]<< "\t";

	}
    }

  template<typename T> 
    void  Reset(T *array, int N)
    {
      for(int i= 0; i< N; i++)
	{
	  array[i] = 0;
	}
    }


 protected:
  int fPhosEventCount;                  /**<Global event counter for this component*/
  AliHLTUInt8_t  fModuleID;             /**<ID of the module this component read data from (0-4)*/
  Bool_t fPrintInfo;                    /**<wether or not to print debugg info to std out*/
  int fPrintInfoFrequncy;               /**<Defines the update frequency for information printet to std out*/
  static const AliHLTComponentDataType fgkInputDataTypes[]; /**<List of  datatypes that can be given to this component*/
 private:
    AliHLTPHOSProcessor(const AliHLTPHOSProcessor & );
    AliHLTPHOSProcessor & operator = (const AliHLTPHOSProcessor &);

};


#endif
