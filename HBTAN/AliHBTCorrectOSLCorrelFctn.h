#ifndef ALIHBTCORRECTOSLCORRELFCTN_H
#define ALIHBTCORRECTOSLCORRELFCTN_H
//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectQ3DCorrelFctn                        //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTCorrectOSLCorrelFctn: public AliHBTOnePairFctn3D
{
  public:
   AliHBTCorrectOSLCorrelFctn(const char* name = "qinvcorrectedCF", 
                               const char* title= "Corrected Q_{inv} Correlation Fonction");
   AliHBTCorrectOSLCorrelFctn(const AliHBTCorrectOSLCorrelFctn& in);
   virtual ~AliHBTCorrectOSLCorrelFctn();
   
  protected:
    TH3D* fMeasCorrelFctn; //!Measured correlation function
    
    TH3D* fSmearedNumer; //! Numerator of smeard q
    TH3D* fSmearedDenom; //! Denominator of smeard q
    TH3D* fMeasNumer;  //! Numerator of ideal q calculated on basis of model equation
    TH3D* fMeasDenom;  //! Denominator of ideal q calculated on basis of model equation

    
  private:
  
    ClassDef(AliHBTCorrectOSLCorrelFctn,1)
};

#endif
