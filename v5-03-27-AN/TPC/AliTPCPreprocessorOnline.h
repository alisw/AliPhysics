#ifndef ALITPCPREPROCESSORONLINE_H
#define ALITPCPREPROCESSORONLINE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCPreprocessorOnline.h,v */

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Preprocessor class for HLT and DAQ                                        //
//  Possible usage: preprocess TPC calibration data                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


class AliTPCCalPad;
class AliTPCCalROC;
class AliTPCCalibViewer;
class TMap;

class AliTPCPreprocessorOnline : public TObject {
public:
   AliTPCPreprocessorOnline();
   AliTPCPreprocessorOnline(const AliTPCPreprocessorOnline &c);
   AliTPCPreprocessorOnline(TMap *map);
   AliTPCPreprocessorOnline &operator = (const AliTPCPreprocessorOnline & param);
   virtual ~AliTPCPreprocessorOnline();
   
   void AddComponent(TObject *obj);
   void DumpToFile(const char* fileName);
   
   TMap *GetMap() {return fMap; } // for debugging
   
         
protected:
   TMap *fMap;       // Map of the AliTPCCalPads
   
   ClassDef(AliTPCPreprocessorOnline,1)    //  TPC preprocessor class
};

#endif



