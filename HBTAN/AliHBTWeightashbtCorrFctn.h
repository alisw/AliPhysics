#ifndef ALIHBTWEIGHTASHBTCORRFCTN_H
#define ALIHBTWEIGHTASHBTCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTashbtCorrFctn.h                           //
//                                                   //
// Class for calculating 3D ashbt correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTWeightashbtCorrFctn: public AliHBTOnePairFctn1D
{
  public:
   AliHBTWeightashbtCorrFctn(const char* name = "asejdzbitiCF", 
                         const char* title= "asHBT Correlation Function"){Info("AliHBTWeightashbtCorrFctn","%s %s",name,title);}

   virtual ~AliHBTWeightashbtCorrFctn(){}
  
   ClassDef(AliHBTWeightashbtCorrFctn,2)
};

#endif
