#ifndef ALIEMCALSPACEFRAME_H
#define ALIEMCALSPACEFRAME_H
/* Copyright(c) 2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////
//          class for EMCAL Space Frame        //
//              Author: Ryan M. Ward           //                                                   
//     California Polytechnic State Univeristy //                                                   
//               rmward@calpoly.edu            // 
/////////////////////////////////////////////////

#include <TNamed.h> 

class AliEMCALSpaceFrame : public TNamed
{
  
public:

  AliEMCALSpaceFrame();
  virtual      ~AliEMCALSpaceFrame() {}

  AliEMCALSpaceFrame(const AliEMCALSpaceFrame& calFrame);
  AliEMCALSpaceFrame & operator = (const AliEMCALSpaceFrame & /*rvalue*/) {
    Fatal("operator =", "not implemented");
    return *this;
  }
  
  //This method assembles the Geometries and places them into the
  //Alice master volume
  void CreateGeometry();

protected:
  
private:

  //space frame parameters from "SINGLE FRAME ASSEMBLY 27D624H.pdf"
  //provided by Lawrence Berkeley Labs, USA
  Int_t fNumCross;                  // Total number of cross beams including end pieces
  Int_t fNumSubSets;                // Total Number of Cross Beam sections in each Half Section
  Double_t fTotalHalfWidth;         // Half the width of each Half Section 
  Double_t fBeginPhi;               // Begining Phi of Cal Frame
  Double_t fEndPhi;                 // Ending Phi of Cal Frame
  Double_t fTotalPhi;               // Total Phi range of Cal Frame
  Double_t fBeginRadius;            // Begining Radius of Cal Frame
  Double_t fHalfFrameTrans;         // z-direction translation for each half frame
  //flange and rib dimensions
  Double_t fFlangeHeight;           // Ending Radius of Flange (Flange is a TGeoTubeSeg)
  Double_t fFlangeWidth;            // Thickness of Flange in z-direction
  Double_t fRibHeight;              // Ending Radius of Rib
  Double_t fRibWidth;               // Thickness of Rib in z-direction 
  //cross beam sections - Cross beams come in two sections- top and bottom
  Double_t fCrossBottomWidth;       // Width along z direction
  Double_t fCrossTopWidth;          // Width along z direction
  Double_t fCrossBottomHeight;      // Tangental thickness relative to CalFrame arc
  Double_t fCrossBottomRadThick;    // Radial thickness relative to center of geometry
  Double_t fCrossBeamArcLength;     // For calulating placement of
  Double_t fCrossBottomStartRadius; // Radial position relative to center of geometry
  Double_t fCrossTopHeight;         // Tangental thickness relative to the center of geometry
  Double_t fCrossTopRadThick;       // Radial thickness relative to CalFrame arc
  Double_t fCrossTopStart;          // Radial position relative to center of geometry
  Double_t fEndRadius;              // Ending Radius of Mother Volume
  Double_t fEndBeamRadThick;	    // Radial Thickness of the End Beams
  Double_t fEndBeamBeginRadius;	    // Starting Radius for End Beams

  ClassDef(AliEMCALSpaceFrame,1)  //Class for EMCAL Space Frame
};

#endif  //ALIEMCALSPACEFRAME_H
