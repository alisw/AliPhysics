// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAYRESIDUALS_H
#define ALIHLTTPCDISPLAYRESIDUALS_H
/** \class AliHLTTPCDisplayResiduals
<pre>
//_____________________________________________________________
// AliHLTTPCDisplayResiduals
//
// Display class for the HLT TPC-Residuals events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TH1.h>
#include <TGraph.h>
#include <AliHLTTPCDisplayMain.h>

class AliHLTTPCDisplayResiduals {  
    
 public:
    AliHLTTPCDisplayResiduals(AliHLTTPCDisplayMain* display) ;
    virtual ~AliHLTTPCDisplayResiduals();
    
    void Draw();
    void Reset();
    void Fill();
    void Save();
 private:
    AliHLTTPCDisplayMain* fDisplay;

    TH1F *fHistallresidualsY;          // histogram for all Y residuals
    TH1F *fHistallresidualsZ;          // histogram for all Z residuals

    TH1F * fHistHits_S;                // histogram for Hits per track length
    TH1F * fHistQ_Track;               // histogram for Charge per track
    TH1F * fHistQ_S;                   // histogram for Charge per track length
    
    TGraph *fGraphresidualsY;          // graph of the Y residuals for one track
    TGraph *fGraphresidualsZ;          // graph of the Z residuals for one track
    TGraph *fGraphresidualsYLength;    // graph of the border of Y residuals for one track

    ClassDef(AliHLTTPCDisplayResiduals,0) 
};

#endif //  ALIHLTTPCDISPLAYRESIDUALS_H
