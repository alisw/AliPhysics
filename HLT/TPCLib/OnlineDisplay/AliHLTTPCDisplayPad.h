// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYPAD_H
#define ALIHLTTPCDISPLAYPAD_H
/** \class AliHLTTPCDisplayPad
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayPad
//
// Display class for the HLT TPC-Pad events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TCanvas.h>
#include <TH1.h>
#include <TH2F.h>
#include <AliHLTTPCDisplayMain.h>
#include <AliHLTTPCDisplayPadRow.h>
#include <AliHLTTPCDisplayFront.h>

class AliHLTTPCDisplayPad : public AliHLTLogging {  
  friend void AliHLTTPCDisplayPadRow::Draw();
  friend void AliHLTTPCDisplayFront::Draw();

 public:
    AliHLTTPCDisplayPad(AliHLTTPCDisplayMain* display) ;
    virtual ~AliHLTTPCDisplayPad();
    
    void Fill(Int_t patch, ULong_t dataBlock, ULong_t dataLen);
    void Draw();
    void Reset();
    void Save();
    void ExecEvent(Int_t event, Int_t x, Int_t y, TObject *selected);

   
 private:
    AliHLTTPCDisplayMain* fDisplay;

    TH1F *fHistpad1;               // histogram for pad in padrow
    TH1F *fHistpad2;               // histogram for pad in padrow
    TH1F *fHistpad3;               // histogram for pad in padrow

    Int_t fNTimes;
    Int_t fBinX[2];                // Minimum / Maximum - Bin on X Axis
    Int_t fTmpEvent;               // Tmp Event for get user range on Axis

    ClassDef(AliHLTTPCDisplayPad,0) 
};

#endif //  ALIHLTTPCDISPLAYPAD_H
