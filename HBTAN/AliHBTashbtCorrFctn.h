#ifndef ALIHBTASHBTCORRFCTN_H
#define ALIHBTASHBTCORRFCTN_H

///////////////////////////////////////////////////////
//                                                   //
// AliHBTashbtCorrFctn.h                           //
//                                                   //
// Class for calculating 3D ashbt correlation       //
// functions                                         //
//                                                   //
///////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTashbtCorrFctn: public AliHBTOnePairFctn1D
{
  public:
   AliHBTashbtCorrFctn(const char* name = "asejdzbitiCF", 
                         const char* title= "asHBT Correlation Function"){Info("AliHBTashbtCorrFctn","%s %s",name,title);}

   virtual ~AliHBTashbtCorrFctn(){}

   ClassDef(AliHBTashbtCorrFctn,2)
};

#endif
