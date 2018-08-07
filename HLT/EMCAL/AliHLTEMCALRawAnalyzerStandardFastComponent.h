#ifndef AliHLTEMCALRawAnalyzerStandardFastComponent_H
#define AliHLTEMCALRawAnalyzerStandardFastComponent_H

// Component structure copied from AliHLTEMCALRawAnalyzerStandardComponent

#include  "AliHLTEMCALRawAnalyzerComponent.h"

//AliHLTCALORawAnalyzerCrudeComponent

class  AliHLTEMCALRawAnalyzerStandardFastComponent : public AliHLTEMCALRawAnalyzerComponent
//class  AliHLTEMCALRawAnalyzerStandardFastComponent : public AliHLTCALORawAnalyzerComponent
{
public:
  AliHLTEMCALRawAnalyzerStandardFastComponent();
  virtual ~AliHLTEMCALRawAnalyzerStandardFastComponent();
  virtual int DoDeinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn(); 
  //  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
 private:
  AliHLTEMCALRawAnalyzerStandardFastComponent( const AliHLTEMCALRawAnalyzerStandardFastComponent  & );
  AliHLTEMCALRawAnalyzerStandardFastComponent & operator = (const AliHLTEMCALRawAnalyzerStandardFastComponent  &);
  // bool TestBoolConst() { return false; };
  //  bool TestBool()  {return  false; };

    

};

#endif
