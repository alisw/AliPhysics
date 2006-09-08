// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYCHARGE_H
#define ALIHLTTPCDISPLAYCHARGE_H
/** \class AliHLTTPCDisplayCharge
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayCharge
//
// Display class for the HLT TPC-Charge events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TCanvas.h>
#include <TH1.h>
#include <AliHLTTPCDisplayMain.h>

class AliHLTTPCDisplayCharge : public AliHLTLogging {  
    
 public:
    AliHLTTPCDisplayCharge(AliHLTTPCDisplayMain* display) ;
    virtual ~AliHLTTPCDisplayCharge();
    
    void Fill();
    void Draw();
    void Reset();
    void Save();
    void ExecEvent(Int_t event, Int_t x, Int_t y, TObject *selected);

 private:
    AliHLTTPCDisplayMain* fDisplay;

    TH1F *fHistcharge;             // histogram for clustercharge

    Int_t fMaxCharge;              // Maximum of Charge
    Int_t fBinX[2];                // Minimum / Maximum - Bin on X Axis
    Int_t fTmpEvent;               // Tmp Event for get user range on Axis

    ClassDef(AliHLTTPCDisplayCharge,0) 
};

#endif //  ALIHLTTPCDISPLAYCHARGE_H
