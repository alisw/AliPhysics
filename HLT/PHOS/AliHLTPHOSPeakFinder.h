#ifndef ALIHLTPHOSPEAKFINDER_H
#define ALIHLTPHOSPEAKFINDER_H
#include <Rtypes.h>
#include "TObject.h"
#include "AliHLTPHOSRawAnalyzer.h"




/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

class AliHLTPHOSPeakFinder : public  TObject, public AliHLTPHOSRawAnalyzer
//class AliHLTPHOSPeakFinder : public  TObject
{
 public:
  AliHLTPHOSPeakFinder();
  AliHLTPHOSPeakFinder(double *dataPtr, double fs);
  AliHLTPHOSPeakFinder(const AliHLTPHOSPeakFinder & );
  AliHLTPHOSPeakFinder & operator = (const AliHLTPHOSPeakFinder)
    {
      return *this; 
    }
  
  virtual ~AliHLTPHOSPeakFinder();
  void FitPeakFinder(int start = 0, int lenght = 100, double *tVector = 0, double *aVector = 0);
  int FindStartIndex(double treshold);
  //  virtual void Analyze(int start = 0, int lenght = 100, double *tVector = 0, double *aVector = 0) const;
  virtual void Analyze() const;
  virtual const char* GetComponentID(){return 0;};
  //   AliHLTPHOSPeakFinder::GetOutputDataSize(long int&, double&)'

  //void AliHLTComponent::GetOutputDataSize(long unsigned int&, double&)

  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ){};
  virtual void GetOutputDataSize(long  int&, double&) {};

  virtual AliHLTComponent* Spawn(){return 0;};

  virtual AliHLTComponentDataType GetOutputDataType(){};
  virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& ) {};
  //(std::vector<AliHLTComponentDataType, std>&){};


  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&){return 0;};


 private:
  void MakeInitialGuess();
  void MakeInitialGuess(int treshold);
  double     fDTofGuess;       /**<Initial guess for t0*/
  double     fDAmplGuess;      /**<Initial guess for amplitude*/
  double   **kfMCovarPtrPtr;   /**<Covariance matrix of the measurements*/
  double   **fPCovarPtrPtr;    /**<Covariance matrix of the estimated parameters*/
  ClassDef(AliHLTPHOSPeakFinder, 2) 
  
    };

#endif
