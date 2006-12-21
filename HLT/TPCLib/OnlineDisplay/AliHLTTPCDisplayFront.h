// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYFRONT_H
#define ALIHLTTPCDISPLAYFRONT_H
/** \class AliHLTTPCDisplayFront
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayFront
//
// Display class for the HLT TPC-Front view 
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <AliHLTTPCDisplayMain.h>

class AliHLTTPCDisplayFront : public AliHLTLogging  {  
    
 public:
    AliHLTTPCDisplayFront(AliHLTTPCDisplayMain* display) ;
    virtual ~AliHLTTPCDisplayFront();
    
    void Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen);
    void Draw();
    void Reset();
    void Save();
    void ExecEvent(Int_t event, Int_t x, Int_t y, TObject *selected);

 private:
    AliHLTTPCDisplayMain* fDisplay;
    TCanvas * fCanvas;

    TH2F *fHistfront;              // histogram for front view of one slice
  TH1F *fHistfrontcl;              // histogram for cluster in front
    Int_t fNTimes;

    Int_t fBinX[2];                // Minimum / Maximum - Bin on X Axis
    Int_t fBinY[2];                // Minimum / Maximum - Bin on Y Axis
    Int_t fTmpEvent;               // Tmp Event for get user range on Axis

    ClassDef(AliHLTTPCDisplayFront,0) 
};

#endif //  ALIHLTTPCDISPLAYFRONT_H
