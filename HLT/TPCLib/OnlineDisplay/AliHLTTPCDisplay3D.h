// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCDISPLAY3D_H
#define ALIHLTTPCDISPLAY3D_H
/** \class AliHLTTPCDisplay3D
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay3D
//
// Display class for the HLT TPC-3D events.
</pre>
*/
// Author: Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>
//*-- Copyright &copy ALICE HLT Group 

#include <TCanvas.h>
#include <TH1.h>
#include <TH2F.h>
#include <AliHLTTPCDisplayMain.h>

class AliHLTTPCDisplay3D : public AliHLTLogging {  
    
 public:
    AliHLTTPCDisplay3D(AliHLTTPCDisplayMain* display, Char_t *gfile ) ;
    virtual ~AliHLTTPCDisplay3D();
    
    void Draw();
    void Save();

    struct AliHLTTPCTrackParameter{
	Int_t nHits;
	Int_t charge;
	Double_t kappa;
	Double_t radius;
	Double_t xyzF[3];
	Double_t xyzL[3];
	Int_t slice;
	Double_t phi0;
	Double_t psi;
	Double_t lambda;
	Double_t pt;
	Int_t id;
	Double_t bfield;
	Double_t s;
    };
    
    AliHLTTPCTrackParameter fTrackParam;

// ---------------------------------------------------
 private:
    void LoadGeometrie(Char_t *gfile);
    void DrawGeomSector(Int_t sector);

    AliHLTTPCDisplayMain* fDisplay;

    TGeometry *fGeom;              // Geometry

    ClassDef(AliHLTTPCDisplay3D,0) 
};

#endif //  ALIHLTTPCDISPLAY3D_H
