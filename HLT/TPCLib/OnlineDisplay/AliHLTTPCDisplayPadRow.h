// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYPADROW_H
#define ALIHLTTPCDISPLAYPADROW_H
/** \class AliHLTTPCDisplayPadRow
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayPadRow
//
// Display class for the HLT TPC-PadRow events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TCanvas.h>
#include <TH1.h>
#include <TH2F.h>
#include <AliHLTTPCDisplayMain.h>

class AliHLTTPCDisplayPadRow : public AliHLTLogging  {  
    
 public:
    AliHLTTPCDisplayPadRow(AliHLTTPCDisplayMain* display) ;
    virtual ~AliHLTTPCDisplayPadRow();  

    void Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen);
    void Draw();
    void Reset();
    void Save();
    void Draw3D();
    void ExecEvent(Int_t event, Int_t x, Int_t y, TObject *selected);

 private:
    AliHLTTPCDisplayMain* fDisplay;

    TH1F *fHistrawcl;              // histogram for cluster in padrow
    TH2F *fHistraw;                // histogram for signals in padrow

    Int_t fcolorbin[20];           // number of entries per colorbin
    Int_t fbinct[20];              // index of colorbin
    Float_t *fpmarr[20];           // contains point data
    Int_t fNTimes;                 // number of timebins
    Int_t fBinX[2];                // Minimum / Maximum - Bin on X Axis
    Int_t fBinY[2];                // Minimum / Maximum - Bin on Y Axis
    Int_t fTmpEvent;               // Tmp Event for get user range on Axis

    ClassDef(AliHLTTPCDisplayPadRow,0) 
};

#endif //  ALIHLTTPCDISPLAYPADROW_H
