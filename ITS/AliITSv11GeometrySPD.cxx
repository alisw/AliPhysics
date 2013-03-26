/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// This class Defines the Geometry for the ITS services and support cones
// outside of the central volume (except for the Central support
// cylinders). Other classes define the rest of the ITS, specifically the
// SSD support cone, the SSD Support central cylinder, the SDD support cone,
// the SDD support central cylinder, the SPD Thermal Shield, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC.
//
//     Here is the calling sequence associated with this file
//   SPDSector(TGeoVolume *moth,TGeoManager *mgr)
//   -----CarbonFiberSector(TGeoVolume *moth,Double_t &xAAtubeCenter0,
//                          Double_t &yAAtubeCenter0,TGeoManager *mgr)
//        -----2* SPDsectorShape(Int_t n,const Double_t *xc,const Double_t *yc,
//        |                      const Double_t *r,const Double_t *ths,
//        |                      const Double_t *the,Int_t npr,Int_t &m,
//        |                      Double_t **xp,Double_t **yp)
//        -----StavesInSector(TGeoVolume *moth,TGeoManager *mgr)
//             -----3* CreaeStave(Int_t layer,TArrayD &sizes,Bool_t addClips,
//             |                  TGeoManager *mgr)
//             |    -----2* CreateHalfStave(Boot_t isRight,Int_t layer,
//             |                            Int_t idxCentral,Int_t idxSide,
//             |                            TArrayD &sizes,Bool_t addClips,
//             |                            TGeoManager *mgr)
//             |         -----CreateGrondingFoil(Bool_t isRight,TArrayD &sizes,
//             |         |                       TGeoManager *mgr)
//             |         |    -----4* CreateGroundingFoilSingle(Int_t type,
//             |         |                                     TArrayD &sizes,
//             |         |                                     TGeoManger *mgr)
//             |         |----CreateLadder(Int_t layer, TArrayD &sizes,
//             |         |                 TGeoManager *mgr)
//             |         |----CreateMCM(Bool_t isRight,TArrayD &sizes,
//             |         |              TGeoManger *mgr)
//             |         |----CreatePixelBus(Bool_t isRight,TArrayD &sizes,
//             |         |                   TGeoManager *mgr)
//             |         -----CreateClip(TArrayD &sizes,TGeoManager *mgr)
//             |----GetSectorMountingPoints(Int_t index,Double_t &x0,
//             |                            Double_t &y0,Double_t &x1,
//             |                            Double_t y1)
//             -----3* ParallelPosition(Double_t dist1,Double_t dist2,
//                                      Double_t phi,Double_t &x,Double_t &y)
//
//     Obsoleate or presently unused routines are: setAddStave(Bool_t *mask),
// CreatePixelBusAndExtensions(...) which calles CreateExtender(...).

/* $Id$ */


// General Root includes
#include <Riostream.h>
#include <TMath.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPolyLine.h>
#include <TPolyMarker.h>

// Root Geometry includes
#include <TGeoCompositeShape.h>
#include <TGeoEltu.h>
#include <TGeoGlobalMagField.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoTube.h> // contains TGeoTubeSeg
#include <TGeoVolume.h>
#include <TGeoXtru.h>
#include <TGeoPcon.h>
#include <TGeoPgon.h>
#include <TGeoArb8.h>

// AliRoot includes
#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

// Declaration file
#include "AliITSv11GeometrySPD.h"
#include "AliITSv11GeomCableRound.h"

// Constant definistions
const Double_t AliITSv11GeometrySPD::fgkGapLadder    =
                      AliITSv11Geometry::fgkmicron*75.; //  75 microns
const Double_t AliITSv11GeometrySPD::fgkGapHalfStave =
                     AliITSv11Geometry::fgkmicron*120.; // 120 microns

using std::endl;
using std::cout;
using std::ios;
ClassImp(AliITSv11GeometrySPD)
//______________________________________________________________________
AliITSv11GeometrySPD::AliITSv11GeometrySPD(/*Double_t gap*/):
AliITSv11Geometry(),// Default constructor of base class
fAddStave(),        // [DEBUG] must be TRUE for all staves which will be
                    // mounted in the sector (used to check overlaps)
fSPDsectorX0(0),    // X of first edge of sector plane for stave
fSPDsectorY0(0),    // Y of first edge of sector plane for stave
fSPDsectorX1(0),    // X of second edge of sector plane for stave
fSPDsectorY1(0),    // Y of second edge of sector plane for stave
fTubeEndSector()    // coordinate of cooling tube ends
{
    //
    // Default constructor.
    // This does not initialize anything and is provided just for
    // completeness. It is recommended to use the other one.
    // The alignment gap is specified as argument (default = 0.0075 cm).
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default constructed AliITSv11GeometrySPD class.
    //
    Int_t i = 0,j=0,k=0;

    for (i = 0; i < 6; i++) fAddStave[i] = kTRUE;
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++){
        this->fTubeEndSector[k][0][i][j] = 0.0;
        this->fTubeEndSector[k][1][i][j] = 0.0;
    } // end for i,j
}
//______________________________________________________________________
AliITSv11GeometrySPD::AliITSv11GeometrySPD(Int_t debug/*, Double_t gap*/):
AliITSv11Geometry(debug),// Default constructor of base class
fAddStave(),        // [DEBUG] must be TRUE for all staves which will be
                    // mounted in the sector (used to check overlaps)
fSPDsectorX0(0),    // X of first edge of sector plane for stave
fSPDsectorY0(0),    // Y of first edge of sector plane for stave
fSPDsectorX1(0),    // X of second edge of sector plane for stave
fSPDsectorY1(0),    // Y of second edge of sector plane for stave
fTubeEndSector()    // coordinate of cooling tube ends
{
    //
    // Constructor with debug setting argument
    // This is the constructor which is recommended to be used.
    // It sets a debug level, and initializes the name of the object.
    // The alignment gap is specified as argument (default = 0.0075 cm).
    // Inputs:
    //    Int_t    debug               Debug level, 0= no debug output.
    // Outputs:
    //    none.
    // Return:
    //    A default constructed AliITSv11GeometrySPD class.
    //
    Int_t i = 0,j=0,k=0;

    for (i = 0; i < 6; i++) fAddStave[i] = kTRUE;
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++){
        this->fTubeEndSector[k][0][i][j] = 0.0;
        this->fTubeEndSector[k][1][i][j] = 0.0;
    } // end for i,j
}
//______________________________________________________________________
AliITSv11GeometrySPD::AliITSv11GeometrySPD(const AliITSv11GeometrySPD &s):
AliITSv11Geometry(s),// Base Class Copy constructor
fAddStave(),        // [DEBUG] must be TRUE for all staves which will be
                    // mounted in the sector (used to check overlaps)
fSPDsectorX0(s.fSPDsectorX0),    // X of first edge of sector plane for stave
fSPDsectorY0(s.fSPDsectorY0),    // Y of first edge of sector plane for stave
fSPDsectorX1(s.fSPDsectorX1),    // X of second edge of sector plane for stave
fSPDsectorY1(s.fSPDsectorY1)     // Y of second edge of sector plane for stave
{
    //
    // Copy Constructor
    // Inputs:
    //    AliITSv11GeometrySPD &s      source class
    // Outputs:
    //    none.
    // Return:
    //    A copy of a AliITSv11GeometrySPD class.
    //
    Int_t i=0,j=0,k=0;

    for (i = 0; i < 6; i++) this->fAddStave[i] = s.fAddStave[i];
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++){
        this->fTubeEndSector[k][0][i][j] = s.fTubeEndSector[k][0][i][j];
        this->fTubeEndSector[k][1][i][j] = s.fTubeEndSector[k][1][i][j];
    } // end for i,j
}
//______________________________________________________________________
AliITSv11GeometrySPD& AliITSv11GeometrySPD::operator=(const
                                               AliITSv11GeometrySPD &s)
{
    //
    // = operator
    // Inputs:
    //    AliITSv11GeometrySPD &s      source class
    // Outputs:
    //    none.
    // Return:
    //    A copy of a AliITSv11GeometrySPD class.
    //
    Int_t i=0,j=0,k=0;

    if(this==&s) return *this;
    for (i = 0; i < 6; i++) this->fAddStave[i] = s.fAddStave[i];
    this->fSPDsectorX0=s.fSPDsectorX0;
    this->fSPDsectorY0=s.fSPDsectorY0;
    this->fSPDsectorX1=s.fSPDsectorX1;
    this->fSPDsectorY1=s.fSPDsectorY1;
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++){
        this->fTubeEndSector[k][0][i][j] = s.fTubeEndSector[k][0][i][j];
        this->fTubeEndSector[k][1][i][j] = s.fTubeEndSector[k][1][i][j];
    } // end for i,j
    return *this;
}
//______________________________________________________________________
TGeoMedium* AliITSv11GeometrySPD::GetMedium(const char* mediumName,
                                            const TGeoManager *mgr) const
{
    //
    // This function is used to recovery any medium
    // used to build the geometry volumes.
    // If the required medium does not exists,
    // a NULL pointer is returned, and an error message is written.
    //
     Char_t itsMediumName[30];

     snprintf(itsMediumName, 30, "ITS_%s", mediumName);
     TGeoMedium* medium = mgr->GetMedium(itsMediumName);
     if (!medium) AliError(Form("Medium <%s> not found", mediumName));

     return medium;
}

//______________________________________________________________________
void AliITSv11GeometrySPD::SPDSector(TGeoVolume *moth, TGeoManager *mgr)
{
    //
    // Creates a single SPD carbon fiber sector and places it
    // in a container volume passed as first argument ('moth').
    // Second argument points to the TGeoManager which coordinates
    // the overall volume creation.
    // The position of the sector is based on distance of
    // closest point of SPD stave to beam pipe
    // (figures all-sections-modules.ps) of 7.22mm at section A-A.
    //

    // Begin_Html
    /*
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/assembly.ps"
     title="SPD     Sector    drawing   with all  cross     sections  defined">
     <p>The    SPD  Sector    definition.    In
     <a   href="http://alice.pd.infn.it/latestdr/Geometric-Revision/assembly.hpgl">HPGL</a>    format.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/assembly-10-modules.ps"
     titile="SPD    All  Sectors   end  view with thermal   sheald">
     <p>The    SPD  all  sector    end  view with thermal   sheald.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/assembly.ps"
     title="SPD     side view cross     section">
     <p>SPD    side view cross     section   with condes    and  thermal   shealds.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-A_A.jpg"
     title="Cross   section   A-A"><p>Cross  section   A-A.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-B_B.jpg"
     title="Cross  updated section   A-A"><p>Cross updated section   A-A.
     <img src="http://physics.mps.ohio-state.edu/~nilsen/ITSfigures/Sezione_layerAA.pdf"
     title="Cross   section   B-B"><p>Cross  section   B-B.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-C_C.jpg"
     title-"Cross   section   C-C"><p>Cross  section   C-C.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-D_D.jpg"
     title="Cross   section   D-D"><p>Cross  section   D-D.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-E_E.jpg"
     title="Cross   section   E-E"><p>Cross  section   E-E.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-F_F.jpg"
     title="Cross   section   F-F"><p>Cross  section   F-F.
     <img src="http://alice.pd.infn.it/latestdr/Geometric-Revision/SECTION-G_G.jpg"
     title="Cross   section   G-G"><p>Cross  section   G-G.
    */
    // End_Html

    // Inputs:
    //    TGeoVolume *moth  Pointer to mother volume where this object
    //                      is to be placed in
    //    TGeoManager *mgr  Pointer to the TGeoManager used, defaule is
    //                      gGeoManager.
    // Outputs:
    //    none.
    // Return:
    //    none.
    // Updated values for kSPDclossesStaveAA, kBeamPipeRadius, and
    // staveThicknessAA are taken from
    // http://physics.mps.ohio-state.edu/~nilsen/ITSfigures/Sezione_layerAA.pdf
    //
    const Double_t kSPDclossesStaveAA   =   7.25* fgkmm;
    const Double_t kSectorStartingAngle = -72.0 * fgkDegree;
    const Int_t    kNSectorsTotal       =  10;
    const Double_t kSectorRelativeAngle =  36.0 * fgkDegree;    // = 360.0 / 10
    const Double_t kBeamPipeRadius      =   0.5 * 59.6 * fgkmm; // diam. = 59.6 mm
  //const Double_t staveThicknessAA     =   0.9 *fgkmm;         // nominal thickness
    const Double_t staveThicknessAA     =   1.02 * fgkmm;       // get from stave geometry.

    Int_t i, j, k;
    Double_t angle, radiusSector, xAAtubeCenter0, yAAtubeCenter0;
    TGeoCombiTrans *secRot = new TGeoCombiTrans(), *comrot;
    TGeoVolume *vCarbonFiberSector[10];
    TGeoMedium *medSPDcf;

    // Define an assembly and fill it with the support of
    // a single carbon fiber sector and staves in it
    medSPDcf = GetMedium("SPD C (M55J)$", mgr);
    for(Int_t is=0; is<10; is++)
    {
	    vCarbonFiberSector[is] = new TGeoVolumeAssembly("ITSSPDCarbonFiberSectorV");
	    vCarbonFiberSector[is]->SetMedium(medSPDcf);
	    CarbonFiberSector(vCarbonFiberSector[is], is, xAAtubeCenter0, yAAtubeCenter0, mgr);
    }

    // Compute the radial shift out of the sectors
    radiusSector = kBeamPipeRadius + kSPDclossesStaveAA + staveThicknessAA;
    radiusSector  = GetSPDSectorTranslation(fSPDsectorX0.At(1), fSPDsectorY0.At(1),
                                            fSPDsectorX1.At(1), fSPDsectorY1.At(1), radiusSector);
  //radiusSector *= radiusSector; // squaring;
  //radiusSector -= xAAtubeCenter0 * xAAtubeCenter0;
  //radiusSector  = -yAAtubeCenter0 + TMath::Sqrt(radiusSector);

    AliDebug(1, Form("SPDSector : radiusSector=%f\n",radiusSector));
    i = 1;
    AliDebug(1, Form("i= %d x0=%f y0=%f x1=%f y1=%f\n", i,
                     fSPDsectorX0.At(i), fSPDsectorY0.At(i),
                     fSPDsectorX1.At(i),fSPDsectorY1.At(i)));

    // add 10 single sectors, by replicating the virtual sector defined above
    // and placing at different angles
    Double_t shiftX, shiftY, tub[2][6][3];
    for(i=0;i<2;i++)for(j=0;j<6;j++)for(k=0;k<3;k++) tub[i][j][k] = fTubeEndSector[0][i][j][k];
    angle = kSectorStartingAngle;
    secRot->RotateZ(angle);
    TGeoVolumeAssembly *vcenteral = new TGeoVolumeAssembly("ITSSPD");
    moth->AddNode(vcenteral, 1, 0);
    for(i = 0; i < kNSectorsTotal; i++) {
        shiftX = -radiusSector * TMath::Sin(angle/fgkRadian);
        shiftY =  radiusSector * TMath::Cos(angle/fgkRadian);
        //cout << "ANGLE = " << angle << endl;
        shiftX += 0.1094 * TMath::Cos((angle + 196.)/fgkRadian);
        shiftY += 0.1094 * TMath::Sin((angle + 196.)/fgkRadian);
        //shiftX -= 0.105;
        //shiftY -= 0.031;
        //shiftX -= 0.11 * TMath::Cos(angle/fgkRadian); // add by Alberto
        //shiftY -= 0.11 * TMath::Sin(angle/fgkRadian); // don't ask me where that 0.11 comes from!
        secRot->SetDx(shiftX);
        secRot->SetDy(shiftY);
        comrot  = new TGeoCombiTrans(*secRot);
        vcenteral->AddNode(vCarbonFiberSector[i],i+1,comrot);
        for(j=0;j<2;j++)for(k=0;k<6;k++) // Transform Tube ends for each sector
            comrot->LocalToMaster(tub[j][k],fTubeEndSector[i][j][k]);
        if(GetDebug(5)) {
            AliInfo(Form("i=%d angle=%g angle[rad]=%g radiusSector=%g "
                         "x=%g y=%g \n",i, angle, angle/fgkRadian,
                         radiusSector, shiftX, shiftY));
        } // end if GetDebug(5)
        angle += kSectorRelativeAngle;
        secRot->RotateZ(kSectorRelativeAngle);
    } // end for i
    if(GetDebug(3)) moth->PrintNodes();
    delete secRot;

    CreateCones(moth);
    CreateServices(moth);
}
//______________________________________________________________________
void AliITSv11GeometrySPD::CarbonFiberSector(TGeoVolume *moth, Int_t sect,
     Double_t &xAAtubeCenter0, Double_t &yAAtubeCenter0, TGeoManager *mgr)
{
    // The method has been modified in order to build a support sector
    // whose shape is dependent on the sector number; the aim is to get
    // as close as possible to the shape inferred from alignment
    // and avoid as much as possible overlaps generated by alignment.
    //
    // Define the detail SPD Carbon fiber support Sector geometry.
    // Based on the drawings:
    /*
      http:///QA-construzione-profilo-modulo.ps
     */
    // - ALICE-Pixel "Costruzione Profilo Modulo" (march 25 2004)
    // - ALICE-SUPPORTO "Costruzione Profilo Modulo"
    // ---
    // Define outside radii as negative, where "outside" means that the
    // center of the arc is outside of the object (feb 16 2004).
    // ---
    // Arguments [the one passed by ref contain output values]:
    // Inputs:
    //   TGeoVolume *moth             the voulme which will contain this object
    //   TGeoManager *mgr             TGeo builder defauls is gGeoManager
    // Outputs:
    //   Double_t   &xAAtubeCenter0  (by ref) x location of the outer surface
    //                               of the cooling tube center for tube 0.
    //   Double_t   &yAAtubeCenter0  (by ref) y location of the outer surface
    //                                of the cooling tube center for tube 0.
    // Return:
    //   none.
    // ---
    // Int the two variables passed by reference values will be stored
    // which will then be used to correctly locate this sector.
    // The information used for this is the distance between the
    // center of the #0 detector and the beam pipe.
    // Measurements are taken at cross section A-A.
    //

    //TGeoMedium *medSPDfs      = 0;//SPD support cone inserto stesalite 4411w
    //TGeoMedium *medSPDfo      = 0;//SPD support cone foam, Rohacell 50A.
    //TGeoMedium *medSPDal      = 0;//SPD support cone SDD mounting bracket Al
    TGeoMedium *medSPDcf     = GetMedium("SPD C (M55J)$", mgr);
    TGeoMedium *medSPDss     = GetMedium("INOX$", mgr);
    TGeoMedium *medSPDcoolfl = GetMedium("Freon$", mgr); //ITSspdCoolingFluid
    //
    const Double_t ksecDz           =  0.5 * 500.0 * fgkmm;
    //const Double_t ksecLen        = 30.0 * fgkmm;
    const Double_t ksecCthick       =  0.2 * fgkmm;
    const Double_t ksecDipLength =  3.2 * fgkmm;
    const Double_t ksecDipRadii  =  0.4 * fgkmm;
    //const Double_t ksecCoolingTubeExtraDepth = 0.86 * fgkmm;
    //
    // The following positions ('ksecX#' and 'ksecY#') and radii ('ksecR#')
    // are the centers and radii of curvature of all the rounded corners
    // between the straight borders of the SPD sector shape.
    // To draw this SPD sector, the following steps are followed:
    // 1) the (ksecX, ksecY) points are plotted
    //    and circles of the specified radii are drawn around them.
    // 2) each pair of consecutive circles is connected by a line
    //    tangent to them, in accordance with the radii being "internal"
    //    or "external" with respect to the closed shape which describes
    //    the sector itself.
    // The resulting connected shape is the section
    // of the SPD sector surface in the transverse plane (XY).
    //
    const Double_t ksecX0   = -10.725 * fgkmm;
    const Double_t ksecY0   = -14.853 * fgkmm;
    const Double_t ksecR0   =  -0.8   * fgkmm; // external

    const Double_t ksecR1   =  +0.6   * fgkmm;
    const Double_t ksecR2   =  +0.6   * fgkmm;
    const Double_t ksecR3   =  -0.6   * fgkmm;
    const Double_t ksecR4   =  +0.8   * fgkmm;
    const Double_t ksecR5   =  +0.8   * fgkmm;
    const Double_t ksecR6   =  +0.6   * fgkmm;
    const Double_t ksecR7   =  -0.6   * fgkmm;
    const Double_t ksecR8   =  +0.6   * fgkmm;
    const Double_t ksecR9   =  -0.6   * fgkmm;
    const Double_t ksecR10   =  +0.6   * fgkmm;
    const Double_t ksecR11   =  -0.6   * fgkmm;
    const Double_t ksecR12   =  +0.85   * fgkmm;

//    // IDEAL GEOMETRY
//     const Double_t ksecX1[10] ={-1.3187,-1.3187,-1.3187,-1.3187,-1.3187,-1.3187,-1.3187,-1.3187,-1.3187,-1.3187};
//     const Double_t ksecY1[10] ={-1.9964,-1.9964,-1.9964,-1.9964,-1.9964,-1.9964,-1.9964,-1.9964,-1.9964,-1.9964};
//     const Double_t ksecX2[10] ={-0.3833,-0.3833,-0.3833,-0.3833,-0.3833,-0.3833,-0.3833,-0.3833,-0.3833,-0.3833};
//     const Double_t ksecY2[10] ={-1.7805,-1.7805,-1.7805,-1.7805,-1.7805,-1.7805,-1.7805,-1.7805,-1.7805,-1.7805};
//     const Double_t ksecX3[10] ={-0.3123,-0.3123,-0.3123,-0.3123,-0.3123,-0.3123,-0.3123,-0.3123,-0.3123,-0.3123};
//     const Double_t ksecY3[10] ={-1.4618,-1.4618,-1.4618,-1.4618,-1.4618,-1.4618,-1.4618,-1.4618,-1.4618,-1.4618};
//     const Double_t ksecX4[10] ={+1.1280,+1.1280,+1.1280,+1.1280,+1.1280,+1.1280,+1.1280,+1.1280,+1.1280,+1.1280};
//     const Double_t ksecY4[10] ={-1.4473,-1.4473,-1.4473,-1.4473,-1.4473,-1.4473,-1.4473,-1.4473,-1.4473,-1.4473};
//     const Double_t ksecX5[10] ={+1.9544,+1.9544,+1.9544,+1.9544,+1.9544,+1.9544,+1.9544,+1.9544,+1.9544,+1.9544};
//     const Double_t ksecY5[10] ={+1.0961,+1.0961,+1.0961,+1.0961,+1.0961,+1.0961,+1.0961,+1.0961,+1.0961,+1.0961};
//     const Double_t ksecX6[10] ={+1.0830,+1.0830,+1.0830,+1.0830,+1.0830,+1.0830,+1.0830,+1.0830,+1.0830,+1.0830};
//     const Double_t ksecY6[10] ={+1.6868,+1.6868,+1.6868,+1.6868,+1.6868,+1.6868,+1.6868,+1.6868,+1.6868,+1.6868};
//     const Double_t ksecX7[10] ={+1.1581,+1.1581,+1.1581,+1.1581,+1.1581,+1.1581,+1.1581,+1.1581,+1.1581,+1.1581};
//     const Double_t ksecY7[10] ={+1.3317,+1.3317,+1.3317,+1.3317,+1.3317,+1.3317,+1.3317,+1.3317,+1.3317,+1.3317};
//     const Double_t ksecX8[10] ={-0.0733,-0.0733,-0.0733,-0.0733,-0.0733,-0.0733,-0.0733,-0.0733,-0.0733,-0.0733};
//     const Double_t ksecY8[10] ={+1.7486,+1.7486,+1.7486,+1.7486,+1.7486,+1.7486,+1.7486,+1.7486,+1.7486,+1.7486};
//     const Double_t ksecX9[10] ={+0.0562,+0.0562,+0.0562,+0.0562,+0.0562,+0.0562,+0.0562,+0.0562,+0.0562,+0.0562};
//     const Double_t ksecY9[10] ={+1.4107,+1.4107,+1.4107,+1.4107,+1.4107,+1.4107,+1.4107,+1.4107,+1.4107,+1.4107};
//     const Double_t ksecX10[10]={-1.2252,-1.2252,-1.2252,-1.2252,-1.2252,-1.2252,-1.2252,-1.2252,-1.2252,-1.2252};
//     const Double_t ksecY10[10]={+1.6298,+1.6298,+1.6298,+1.6298,+1.6298,+1.6298,+1.6298,+1.6298,+1.6298,+1.6298};
//     const Double_t ksecX11[10]={-1.0445,-1.0445,-1.0445,-1.0445,-1.0445,-1.0445,-1.0445,-1.0445,-1.0445,-1.0445};
//     const Double_t ksecY11[10]={+1.3162,+1.3162,+1.3162,+1.3162,+1.3162,+1.3162,+1.3162,+1.3162,+1.3162,+1.3162};
//     const Double_t ksecX12[10]={-2.2276,-2.2276,-2.2276,-2.2276,-2.2276,-2.2276,-2.2276,-2.2276,-2.2276,-2.2276};
//     const Double_t ksecY12[10]={+1.2948,+1.2948,+1.2948,+1.2948,+1.2948,+1.2948,+1.2948,+1.2948,+1.2948,+1.2948};
  

//    MODIFIED GEOMETRY according with partial alignment of Staves relative to Sectors
//    last numbers: 2010/06/11 (ML)

    const Double_t ksecX1[10]={-1.305917, -1.322242, -1.300649, -1.298700, -1.290830, -1.274307, -1.276433, -1.286468, -1.274381, -1.314864};
    const Double_t ksecY1[10]={-1.997857, -2.018611, -2.005854, -2.004897, -1.995517, -2.002552, -1.995860, -2.021062, -2.012931, -2.043967};
    const Double_t ksecX2[10]={-0.366115, -0.385562, -0.372689, -0.365682, -0.348432, -0.348442, -0.342468, -0.354071, -0.346900, -0.381275};
    const Double_t ksecY2[10]={-1.801679, -1.808306, -1.759315, -1.778851, -1.811655, -1.747888, -1.773811, -1.792427, -1.764514, -1.820324};
//     const Double_t ksecX1[10]={-1.305917, -1.322242, -1.300649, -1.298700, -1.290830, -1.274307, -1.276433, -1.286468, -1.274381, -1.325864};
//     const Double_t ksecY1[10]={-1.997857, -2.018611, -2.005854, -2.004897, -1.995517, -2.002552, -1.995860, -2.021062, -2.012931, -2.032967};
//     const Double_t ksecX2[10]={-0.366115, -0.385562, -0.372689, -0.365682, -0.348432, -0.348442, -0.342468, -0.354071, -0.346900, -0.392275};
//     const Double_t ksecY2[10]={-1.801679, -1.808306, -1.759315, -1.778851, -1.811655, -1.747888, -1.773811, -1.792427, -1.764514, -1.809324};
    const Double_t ksecX3[10]={-0.314030, -0.315531, -0.347521, -0.337675, -0.300420, -0.378487, -0.330729, -0.330850, -0.362360, -0.321097};
    const Double_t ksecY3[10]={-1.452488, -1.460418, -1.447060, -1.443146, -1.472410, -1.430019, -1.469073, -1.472048, -1.462010, -1.444355};
    const Double_t ksecX4[10]={1.124299, 1.124162, 1.089523, 1.095520, 1.136171, 1.058616, 1.105626, 1.106433, 1.077455, 1.117946};
    const Double_t ksecY4[10]={-1.458714, -1.452649, -1.465297, -1.492717, -1.494665, -1.447732, -1.493369, -1.488126, -1.452925, -1.443447};
    const Double_t ksecX5[10]={1.951621, 1.939284, 1.931830, 1.935235, 1.952206, 1.939082, 1.924822, 1.940114, 1.918160, 1.960017};
    const Double_t ksecY5[10]={1.092731, 1.118870, 1.129765, 1.129422, 1.081511, 1.127387, 1.103960, 1.101784, 1.121428, 1.150110};
    const Double_t ksecX6[10]={1.070070, 1.048297, 1.035920, 1.049049, 1.083621, 1.045882, 1.050399, 1.067823, 1.037967, 1.070850};
    const Double_t ksecY6[10]={1.667590, 1.678571, 1.681383, 1.696892, 1.676520, 1.683470, 1.689988, 1.691111, 1.698432, 1.712770};
    const Double_t ksecX7[10]={1.139398, 1.150471, 1.150074, 1.132807, 1.150192, 1.124064, 1.124335, 1.137723, 1.143056, 1.130568};
    const Double_t ksecY7[10]={1.345588, 1.356062, 1.342468, 1.320467, 1.335807, 1.334477, 1.328622, 1.347184, 1.319861, 1.308420};
    const Double_t ksecX8[10]={-0.096963, -0.098603, -0.095286, -0.099990, -0.075132, -0.121593, -0.108673, -0.104237, -0.092082, -0.104044};
    const Double_t ksecY8[10]={1.751207, 1.731467, 1.726908, 1.734219, 1.766159, 1.718203, 1.741891, 1.739743, 1.728288, 1.718046};
    const Double_t ksecX9[10]={0.047615, 0.087875, 0.034917, 0.071603, 0.026468, 0.091619, 0.051994, 0.059947, 0.079785, 0.043443};
    const Double_t ksecY9[10]={1.414699, 1.403187, 1.399061, 1.403430, 1.435056, 1.384557, 1.397692, 1.420269, 1.391372, 1.398954};
    const Double_t ksecX10[10]={-1.233255, -1.186874, -1.246702, -1.213368, -1.259425, -1.190067, -1.225655, -1.224171, -1.197833, -1.237182};
    const Double_t ksecY10[10]={1.635767, 1.646249, 1.617336, 1.608928, 1.636944, 1.602583, 1.630504, 1.629065, 1.624295, 1.620934};
    const Double_t ksecX11[10]={-1.018270, -1.031317, -0.960524, -1.001155, -1.045437, -0.986867, -1.002685, -1.017369, -1.005614, -0.985385};
    const Double_t ksecY11[10]={1.318108, 1.330683, 1.301572, 1.314410, 1.326680, 1.295226, 1.306372, 1.309414, 1.306542, 1.307086};
    const Double_t ksecX12[10]={-2.199004, -2.214964, -2.139247, -2.180547, -2.224505, -2.165324, -2.175883, -2.193485, -2.183227, -2.161570};
    const Double_t ksecY12[10]={1.317677, 1.303982, 1.317057, 1.324766, 1.339537, 1.312715, 1.359642, 1.343638, 1.330234, 1.340836};


    const Double_t ksecR13  =  -0.8   * fgkmm; // external
    const Double_t ksecAngleSide13 = 36.0 * fgkDegree;
    //
    const Int_t ksecNRadii = 20;
    const Int_t ksecNPointsPerRadii = 4;
    const Int_t ksecNCoolingTubeDips = 6;
    //
    // Since the rounded parts are approximated by a regular polygon
    // and a cooling tube of the propper diameter must fit, a scaling factor
    // increases the size of the polygon for the tube to fit.
    //const Double_t ksecRCoolScale = 1./TMath::Cos(TMath::Pi()/
    //                                      (Double_t)ksecNPointsPerRadii);
    const Double_t ksecZEndLen   = 30.000 * fgkmm;
    //const Double_t ksecZFlangLen = 45.000 * fgkmm;
    const Double_t ksecTl        =  0.860 * fgkmm;
    const Double_t ksecCthick2   =  0.600 * fgkmm;
    //const Double_t ksecCthick3  =  1.80  * fgkmm;
    //const Double_t ksecSidelen  = 22.0   * fgkmm;
    //const Double_t ksecSideD5   =  3.679 * fgkmm;
    //const Double_t ksecSideD12  =  7.066 * fgkmm;
    const Double_t ksecRCoolOut  = 2.400 * fgkmm;
    const Double_t ksecRCoolIn   = 2.000 * fgkmm;
    const Double_t ksecDl1       = 5.900 * fgkmm;
    const Double_t ksecDl2       = 8.035 * fgkmm;
    const Double_t ksecDl3       = 4.553 * fgkmm;
    const Double_t ksecDl4       = 6.978 * fgkmm;
    const Double_t ksecDl5       = 6.978 * fgkmm;
    const Double_t ksecDl6       = 6.978 * fgkmm;
    const Double_t ksecCoolTubeThick  = 0.04  * fgkmm;
    const Double_t ksecCoolTubeROuter = 2.6   * fgkmm;
    const Double_t ksecCoolTubeFlatX  = 3.696 * fgkmm;
    const Double_t ksecCoolTubeFlatY  = 0.68  * fgkmm;
    //const Double_t ksecBeamX0 = 0.0 * fgkmm; // guess
    //const Double_t ksecBeamY0 = (15.223 + 40.) * fgkmm; // guess
    //
    // redefine some of the points already defined above
    // in the format of arrays (???)
    const Int_t ksecNPoints = (ksecNPointsPerRadii + 1) * ksecNRadii + 8;
    Double_t secX[ksecNRadii] = {
        ksecX0,  ksecX1[sect],  -1000.0,
        ksecX2[sect],  ksecX3[sect],  -1000.0,
        ksecX4[sect],  ksecX5[sect],  -1000.0,
        ksecX6[sect],  ksecX7[sect],  -1000.0,
        ksecX8[sect],  ksecX9[sect],  -1000.0,
        ksecX10[sect], ksecX11[sect], -1000.0,
        ksecX12[sect], -1000.0
    };
    Double_t secY[ksecNRadii] = {
        ksecY0,  ksecY1[sect],  -1000.0,
        ksecY2[sect],  ksecY3[sect],  -1000.0,
        ksecY4[sect],  ksecY5[sect],  -1000.0,
        ksecY6[sect],  ksecY7[sect],  -1000.0,
        ksecY8[sect],  ksecY9[sect],  -1000.0,
        ksecY10[sect], ksecY11[sect], -1000.0,
        ksecY12[sect], -1000.0
    };
    Double_t secR[ksecNRadii] = {
        ksecR0,  ksecR1,  -.5 * ksecDipLength - ksecDipRadii,
        ksecR2,  ksecR3,  -.5 * ksecDipLength - ksecDipRadii,
        ksecR4,  ksecR5,  -.5 * ksecDipLength - ksecDipRadii,
        ksecR6,  ksecR7,  -.5 * ksecDipLength - ksecDipRadii,
        ksecR8,  ksecR9,  -.5 * ksecDipLength - ksecDipRadii,
        ksecR10, ksecR11, -.5 * ksecDipLength - ksecDipRadii,
        ksecR12, ksecR13
    };

    Double_t secX2[ksecNRadii];
    Double_t secY2[ksecNRadii];
    Double_t secR2[ksecNRadii] = {
        ksecR0,  ksecR1,  ksecRCoolOut,
        ksecR2,  ksecR3,  ksecRCoolOut,
        ksecR4,  ksecR5,  ksecRCoolOut,
        ksecR6,  ksecR7,  ksecRCoolOut,
        ksecR8,  ksecR9,  ksecRCoolOut,
        ksecR10, ksecR11, ksecRCoolOut,
        ksecR12, ksecR13
    };
    Double_t secDip2[ksecNCoolingTubeDips] = {
        ksecDl1, ksecDl2, ksecDl3,
        ksecDl4, ksecDl5, ksecDl6
    };
    Double_t secX3[ksecNRadii];
    Double_t secY3[ksecNRadii];
    const Int_t ksecDipIndex[ksecNCoolingTubeDips] = {2, 5, 8, 11, 14, 17};
    Double_t secAngleStart[ksecNRadii];
    Double_t secAngleEnd[ksecNRadii];
    for(Int_t i = 0; i < ksecNRadii; i++)secAngleEnd[i] = 0.;
    Double_t secAngleStart2[ksecNRadii];
    Double_t secAngleEnd2[ksecNRadii];
    Double_t secAngleTurbo[ksecNCoolingTubeDips] = {0., 0., 0., 0., 0., 0.0};
    //Double_t secAngleStart3[ksecNRadii];
    //Double_t secAngleEnd3[ksecNRadii];
    Double_t  xpp[ksecNPoints],  ypp[ksecNPoints];
    Double_t  xpp2[ksecNPoints], ypp2[ksecNPoints];
    Double_t *xp[ksecNRadii],   *xp2[ksecNRadii];
    Double_t *yp[ksecNRadii],   *yp2[ksecNRadii];
    TGeoXtru *sA0,  *sA1, *sB0, *sB1;
    TGeoCompositeShape *sA2, *sB2;
    TGeoBBox *sB3;
    TGeoEltu *sTA0, *sTA1;
    TGeoTube *sTB0, *sTB1; //,*sM0;
    TGeoRotation    *rot;
    TGeoTranslation *trans;
    TGeoCombiTrans  *rotrans;
    Double_t t, t0, t1, a, b, x0, y0,z0, x1, y1;
    Int_t i, j, k, m;
    Bool_t tst;

    if(!moth) {
        AliError("Container volume (argument) is NULL");
        return;
    } // end if(!moth)
    for(i = 0; i < ksecNRadii; i++) {
        xp[i]  = &(xpp[i*(ksecNPointsPerRadii+1)]);
        yp[i]  = &(ypp[i*(ksecNPointsPerRadii+1)]);
        xp2[i] = &(xpp2[i*(ksecNPointsPerRadii+1)]);
        yp2[i] = &(ypp2[i*(ksecNPointsPerRadii+1)]);
        secX2[i] = secX[i];
        secY2[i] = secY[i];
        secX3[i] = secX[i];
        secY3[i] = secY[i];
    } // end for i
    //
    // find starting and ending angles for all but cooling tube sections
    secAngleStart[0] = 0.5 * ksecAngleSide13;
    for(i = 0; i < ksecNRadii - 2; i++) {
        tst = kFALSE;
        for(j=0;j<ksecNCoolingTubeDips;j++) tst = (tst||i==ksecDipIndex[j]);
        if (tst) continue;
        tst = kFALSE;
        for(j=0;j<ksecNCoolingTubeDips;j++) tst =(tst||(i+1)==ksecDipIndex[j]);
        if (tst) j = i+2; else j = i+1;
        AnglesForRoundedCorners(secX[i],secY[i],secR[i],secX[j],secY[j],
                                secR[j],t0,t1);
        secAngleEnd[i]   = t0;
        secAngleStart[j] = t1;
        if(secR[i] > 0.0 && secR[j] > 0.0) {
            if(secAngleStart[i] > secAngleEnd[i]) secAngleEnd[i] += 360.0;
        } // end if(secR[i]>0.0 && secR[j]>0.0)
        secAngleStart2[i] = secAngleStart[i];
        secAngleEnd2[i]   = secAngleEnd[i];
    } // end for i
    secAngleEnd[ksecNRadii-2] = secAngleStart[ksecNRadii-2] +
                   (secAngleEnd[ksecNRadii-5] - secAngleStart[ksecNRadii-5]);
    if (secAngleEnd[ksecNRadii-2] < 0.0) secAngleEnd[ksecNRadii-2] += 360.0;
    secAngleStart[ksecNRadii-1]  = secAngleEnd[ksecNRadii-2] - 180.0;
    secAngleEnd[ksecNRadii-1]    = secAngleStart[0];
    secAngleStart2[ksecNRadii-2] = secAngleStart[ksecNRadii-2];
    secAngleEnd2[ksecNRadii-2]   = secAngleEnd[ksecNRadii-2];
    secAngleStart2[ksecNRadii-1] = secAngleStart[ksecNRadii-1];
    secAngleEnd2[ksecNRadii-1]   = secAngleEnd[ksecNRadii-1];
    //
    // find location of circle last rounded corner.
    i = 0;
    j = ksecNRadii - 2;
    t0 = TanD(secAngleStart[i]-90.);
    t1 = TanD(secAngleEnd[j]-90.);
    t  = secY[i] - secY[j];
    // NOTE: secR[i=0] < 0; secR[j=18] > 0; and secR[j+1=19] < 0
    t += (-secR[i]+secR[j+1]) * SinD(secAngleStart[i]);
    t -= (secR[j]-secR[j+1]) * SinD(secAngleEnd[j]);
    t += t1 * secX[j] - t0*secX[i];
    t += t1 * (secR[j] - secR[j+1]) * CosD(secAngleEnd[j]);
    t -= t0 * (-secR[i]+secR[j+1]) * CosD(secAngleStart[i]);
    secX[ksecNRadii-1] = t / (t1-t0);
    secY[ksecNRadii-1] = TanD(90.0+0.5*ksecAngleSide13)*
        (secX[ksecNRadii-1]-secX[0])+secY[0];
    secX2[ksecNRadii-1] = secX[ksecNRadii-1];
    secY2[ksecNRadii-1] = secY[ksecNRadii-1];
    secX3[ksecNRadii-1] = secX[ksecNRadii-1];
    secY3[ksecNRadii-1] = secY[ksecNRadii-1];

    // find location of cooling tube centers
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        j = ksecDipIndex[i];
        x0 = secX[j-1] + TMath::Abs(secR[j-1]) * CosD(secAngleEnd[j-1]);
        y0 = secY[j-1] + TMath::Abs(secR[j-1]) * SinD(secAngleEnd[j-1]);
        x1 = secX[j+1] + TMath::Abs(secR[j+1]) * CosD(secAngleStart[j+1]);
        y1 = secY[j+1] + TMath::Abs(secR[j+1]) * SinD(secAngleStart[j+1]);
        t0 = TMath::Sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
        t  = secDip2[i] / t0;
        a  = x0+(x1-x0) * t;
        b  = y0+(y1-y0) * t;
        if(i == 0) {
            // get location of tube center->Surface for locating
            // this sector around the beam pipe.
            // This needs to be double checked, but I need my notes for that.
            // (Bjorn Nilsen)
            xAAtubeCenter0 = x0 + (x1 - x0) * t * 0.5;
            yAAtubeCenter0 = y0 + (y1 - y0) * t * 0.5;
        }// end if i==0
        if(a + b*(a - x0) / (b - y0) > 0.0) {
            secX[j]  = a + TMath::Abs(y1-y0) * 2.0 * ksecDipRadii/t0;
            secY[j]  = b - TMath::Sign(2.0*ksecDipRadii,y1-y0) * (x1-x0)/t0;
            secX2[j] = a + TMath::Abs(y1-y0) * ksecTl/t0;
            secY2[j] = b - TMath::Sign(ksecTl,y1-y0) * (x1-x0) / t0;
            secX3[j] = a + TMath::Abs(y1-y0) *
                       (2.0*ksecDipRadii-0.5*ksecCoolTubeFlatY)/t0;
            secY3[j] = b - TMath::Sign(2.0*ksecDipRadii-0.5*ksecCoolTubeFlatY,
                                       y1-y0)*(x1-x0)/t0;
        } else {
            secX[j] = a - TMath::Abs(y1-y0)*2.0*ksecDipRadii/t0;
            secY[j] = b + TMath::Sign(2.0*ksecDipRadii,y1-y0)*(x1-x0)/t0;
            secX2[j] = a - TMath::Abs(y1-y0)*ksecTl/t0;
            secY2[j] = b + TMath::Sign(ksecTl,y1-y0)*(x1-x0)/t0;
            secX3[j] = a - TMath::Abs(y1-y0)*(2.0*ksecDipRadii-0.5*
                                                  ksecCoolTubeFlatY)/t0;
            secY3[j] = b + TMath::Sign(2.0*ksecDipRadii-0.5*ksecCoolTubeFlatY,
                                       y1-y0)*(x1-x0)/t0;
        } // end if(a+b*(a-x0)/(b-y0)>0.0)

          // Set up Start and End angles to correspond to start/end of dips.
        t1 = (secDip2[i]-TMath::Abs(secR[j])) / t0;
        secAngleStart[j] =TMath::RadToDeg()*TMath::ATan2(y0+(y1-y0)*t1-secY[j],
                                                        x0+(x1-x0)*t1-secX[j]);
        if (secAngleStart[j]<0.0) secAngleStart[j] += 360.0;
        secAngleStart2[j] = secAngleStart[j];
        t1 = (secDip2[i]+TMath::Abs(secR[j]))/t0;
        secAngleEnd[j] = TMath::RadToDeg()*TMath::ATan2(y0+(y1-y0)*t1-secY[j],
                                                        x0+(x1-x0)*t1-secX[j]);
        if (secAngleEnd[j]<0.0) secAngleEnd[j] += 360.0;
        secAngleEnd2[j] = secAngleEnd[j];
        if (secAngleEnd[j]>secAngleStart[j]) secAngleEnd[j] -= 360.0;
        secR[j] = TMath::Sqrt(secR[j]*secR[j]+4.0*ksecDipRadii*ksecDipRadii);
    } // end for i

    // Special cases
    secAngleStart2[8] -= 360.;
    secAngleStart2[11] -= 360.;

    SPDsectorShape(ksecNRadii, secX, secY, secR, secAngleStart, secAngleEnd,
                   ksecNPointsPerRadii, m, xp, yp);

    //  Fix up dips to be square.
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        j = ksecDipIndex[i];
        t = 0.5*ksecDipLength+ksecDipRadii;
        t0 = TMath::RadToDeg()*TMath::ATan(2.0*ksecDipRadii/t);
        t1 = secAngleEnd[j] + t0;
        t0 = secAngleStart[j] - t0;
        x0 = xp[j][1] = secX[j] + t*CosD(t0);
        y0 = yp[j][1] = secY[j] + t*SinD(t0);
        x1 = xp[j][ksecNPointsPerRadii-1] = secX[j] + t*CosD(t1);
        y1 = yp[j][ksecNPointsPerRadii-1] = secY[j] + t*SinD(t1);
        t0 = 1./((Double_t)(ksecNPointsPerRadii-2));
        for(k = 2; k < ksecNPointsPerRadii - 1; k++) {
            // extra points spread them out.
            t = ((Double_t)(k-1)) * t0;
            xp[j][k] = x0+(x1-x0) * t;
            yp[j][k] = y0+(y1-y0) * t;
        } // end for k
        secAngleTurbo[i] = -TMath::RadToDeg() * TMath::ATan2(y1-y0, x1-x0);
        if(GetDebug(3)) {
            AliInfo(
                Form("i=%d -- angle=%f -- x0,y0=(%f, %f) -- x1,y1=(%f, %f)",
                     i, secAngleTurbo[i], x0, y0, x1, y1));
        } // end if GetDebug(3)
    } // end for i
    sA0 = new TGeoXtru(2);
    sA0->SetName("SectorA0");
    sA0->DefinePolygon(m, xpp, ypp);
    sA0->DefineSection(0, -ksecDz);
    sA0->DefineSection(1,  ksecDz);

    // store the edges of each XY segment which defines
    // one of the plane zones where staves will have to be placed
    fSPDsectorX0.Set(ksecNCoolingTubeDips);
    fSPDsectorY0.Set(ksecNCoolingTubeDips);
    fSPDsectorX1.Set(ksecNCoolingTubeDips);
    fSPDsectorY1.Set(ksecNCoolingTubeDips);
    Int_t ixy0, ixy1;
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        // Find index in xpp[] and ypp[] corresponding to where the
        // SPD ladders are to be attached. Order them according to
        // the ALICE numbering schema. Using array of indexes (+-1 for
        // cooling tubes. For any "bend/dip/edge, there are
        // ksecNPointsPerRadii+1 points involved.
        if(i == 0) j = 1;
        else if (i == 1) j = 0;
        else j = i;
        ixy0 = (ksecDipIndex[j]-1)*(ksecNPointsPerRadii+1)+
            (ksecNPointsPerRadii);
        ixy1 = (ksecDipIndex[j]+1) * (ksecNPointsPerRadii+1);
        fSPDsectorX0[i] = sA0->GetX(ixy0);
        fSPDsectorY0[i] = sA0->GetY(ixy0);
        fSPDsectorX1[i] = sA0->GetX(ixy1);
        fSPDsectorY1[i] = sA0->GetY(ixy1);
    } // end for i

    //printf("SectorA#%d ",0);
    InsidePoint(xpp[m-1],ypp[m-1],xpp[0],ypp[0],xpp[1],ypp[1],ksecCthick,
                xpp2[0],ypp2[0]);
    for(i = 1; i < m - 1; i++) {
        j = i / (ksecNPointsPerRadii+1);
        //printf("SectorA#%d ",i);
        InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],
                    ksecCthick,xpp2[i],ypp2[i]);
    } // end for i
    //printf("SectorA#%d ",m);
    InsidePoint(xpp[m-2],ypp[m-2],xpp[m-1],ypp[m-1],xpp[0],ypp[0],
                ksecCthick,xpp2[m-1],ypp2[m-1]);
    // Fix center value of cooling tube dip and
    // find location of cooling tube centers
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        j = ksecDipIndex[i];
        x0 = xp2[j][1];
        y0 = yp2[j][1];
        x1 = xp2[j][ksecNPointsPerRadii-1];
        y1 = yp2[j][ksecNPointsPerRadii-1];
        t0 = TMath::Sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
        t  = secDip2[i]/t0;
        for(k = 2; k < ksecNPointsPerRadii - 1; k++) {
            // extra points spread them out.
            t = ((Double_t)(k-1)) * t0;
            xp2[j][k] = x0+(x1-x0) * t;
            yp2[j][k] = y0+(y1-y0) * t;
        } // end for k
    } // end for i
    sA1 = new TGeoXtru(2);
    sA1->SetName("SectorA1");
    sA1->DefinePolygon(m, xpp2, ypp2);
    sA1->DefineSection(0, -ksecDz-ksecCthick2);
    sA1->DefineSection(1,  ksecDz+ksecCthick2);

    sA2 = new TGeoCompositeShape("ITS SPD Carbon fiber support Sector A0",
				 "SectorA0-SectorA1");
    //
    // Error in TGeoEltu. Semi-axis X must be < Semi-axis Y (?).
    sTA0 = new TGeoEltu("ITS SPD Cooling Tube TA0", 0.5 * ksecCoolTubeFlatY,
                        0.5 * ksecCoolTubeFlatX, ksecDz);
    sTA1 = new TGeoEltu("ITS SPD Cooling Tube coolant TA1",
                        sTA0->GetA() - ksecCoolTubeThick,
                        sTA0->GetB()-ksecCoolTubeThick,ksecDz);
    SPDsectorShape(ksecNRadii,secX2,secY2,secR2,secAngleStart2,secAngleEnd2,
                   ksecNPointsPerRadii, m, xp, yp);
    sB0 = new TGeoXtru(2);
    sB0->SetName("EndB0");
    sB0->DefinePolygon(m, xpp, ypp);
    sB0->DefineSection(0, ksecDz);
    sB0->DefineSection(1, ksecDz + ksecZEndLen);

    //printf("SectorB#%d ",0);
  // Points around the most sharpened tips have to be avoided - M.S. 24 feb 09
    const Int_t nSpecialPoints = 5;
    const Int_t kSpecialPoints[nSpecialPoints] = {7, 17, 47, 62, 77};
    Int_t i2 = 0;
    InsidePoint(xpp[m-1],ypp[m-1],xpp[0],ypp[0],xpp[1],ypp[1],
                ksecCthick2,xpp2[i2],ypp2[i2]);
    for(i = 1; i < m - 1; i++) {
        t = ksecCthick2;
        for(k = 0; k < ksecNCoolingTubeDips; k++)
            if((i/(ksecNPointsPerRadii+1))==ksecDipIndex[k])
                if(!(ksecDipIndex[k]*(ksecNPointsPerRadii+1) == i ||
                     ksecDipIndex[k]*(ksecNPointsPerRadii+1) +
                     ksecNPointsPerRadii == i))
                    t = ksecRCoolOut-ksecRCoolIn;
        //printf("SectorB#%d ",i);
	Bool_t useThisPoint = kTRUE;
	for(Int_t ii = 0; ii < nSpecialPoints; ii++)
	  if ( (i == kSpecialPoints[ii] - 1) ||
	       (i == kSpecialPoints[ii] + 1)   ) useThisPoint = kFALSE;
	if (useThisPoint) {
	  i2++;
	  InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],t,
		      xpp2[i2],ypp2[i2]);
	}
    }// end for i
    //printf("SectorB#%d ",m);
    i2++;
    InsidePoint(xpp[m-2],ypp[m-2],xpp[m-1],ypp[m-1],xpp[0],ypp[0],
                ksecCthick2,xpp2[i2],ypp2[i2]);
    sB1 = new TGeoXtru(2);
    sB1->SetName("EndB1");
    sB1->DefinePolygon(i2+1, xpp2, ypp2);
    sB1->DefineSection(0,sB0->GetZ(0)-ksecCthick2);
    sB1->DefineSection(1,sB0->GetZ(1)+ksecCthick2);

    sB2 = new TGeoCompositeShape("ITS SPD Carbon fiber support Sector End B0",
				 "EndB0-EndB1");
    // SPD sector mount blocks
    const Double_t kMountBlock[3] = {0.5*(1.8-0.2)*fgkmm,0.5*22.0*fgkmm,
                                     0.5*45.0*fgkmm};
    sB3 = new TGeoBBox((Double_t*)kMountBlock);
    // SPD sector mount block screws and nuts (M.S. - 27 oct 2012)
    const Double_t kMountBlockM3ScrewR = 0.5*3.0*fgkmm; // Metric screw
    const Double_t kMountBlockHead1R   = 0.5*8.0*fgkmm;
    const Double_t kMountBlockHead1H   = 1.0*fgkmm;
    const Double_t kMountBlockHead2R   = 0.5*6.0*fgkmm;
    const Double_t kMountBlockHead2H   = 2.7*fgkmm;
    const Double_t kMountBlockM3NutR   = 1.8*kMountBlockM3ScrewR; // Metric nut
    const Double_t kMountBlockM3NutH   = kMountBlockM3NutR; // Metric nut
    TGeoTube *sM3 = new TGeoTube(0, kMountBlockM3ScrewR, sB3->GetDX());
    TGeoTube *sD1 = new TGeoTube(0, kMountBlockHead1R,kMountBlockHead1H/2);
    TGeoTube *sD2 = new TGeoTube(0, kMountBlockHead2R,kMountBlockHead2H/2);
    TGeoPgon *sN3 = new TGeoPgon(0, 360, 6, 2);
    sN3->DefineSection(0,-kMountBlockM3NutH/2, 0, kMountBlockM3NutR);
    sN3->DefineSection(1, kMountBlockM3NutH/2, 0, kMountBlockM3NutR);
    // SPD sector cooling tubes
    sTB0 = new TGeoTube("ITS SPD Cooling Tube End TB0", 0.0,
                   0.5*ksecCoolTubeROuter,0.5*(sB0->GetZ(1)-sB0->GetZ(0)));
    sTB1 = new TGeoTube("ITS SPD Cooling Tube End coolant TB0", 0.0,
                        sTB0->GetRmax() - ksecCoolTubeThick,sTB0->GetDz());
    //
    if(GetDebug(3)) {
        if(medSPDcf) medSPDcf->Dump(); else AliInfo("medSPDcf = 0");
        if(medSPDss) medSPDss->Dump(); else AliInfo("medSPDss = 0");
        if(medSPDcoolfl) medSPDcoolfl->Dump();else AliInfo("medSPDcoolfl = 0");
        sA0->InspectShape();
        sA1->InspectShape();
        sB0->InspectShape();
        sB1->InspectShape();
        sB2->InspectShape();
    } // end if(GetDebug(3))

    // create the assembly of the support and place staves on it
    TGeoVolumeAssembly *vM0 = new TGeoVolumeAssembly(
                                         "ITSSPDSensitiveVirtualvolumeM0");
    StavesInSector(vM0);
    // create other volumes with some graphical settings
    TGeoVolume *vA0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorA0",
                                     sA2, medSPDcf);
    vA0->SetVisibility(kTRUE);
    vA0->SetLineColor(4); // Blue
    vA0->SetLineWidth(1);
    vA0->SetFillColor(vA0->GetLineColor());
    vA0->SetFillStyle(4010); // 10% transparent
    TGeoVolume *vTA0 = new TGeoVolume("ITSSPDCoolingTubeTA0", sTA0, medSPDss);
    vTA0->SetVisibility(kTRUE);
    vTA0->SetLineColor(15); // gray
    vTA0->SetLineWidth(1);
    vTA0->SetFillColor(vTA0->GetLineColor());
    vTA0->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vTA1 = new TGeoVolume("ITSSPDCoolingTubeFluidTA1",
                                      sTA1, medSPDcoolfl);
    vTA1->SetVisibility(kTRUE);
    vTA1->SetLineColor(6); // Purple
    vTA1->SetLineWidth(1);
    vTA1->SetFillColor(vTA1->GetLineColor());
    vTA1->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vB0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndB0",
                                     sB2, medSPDcf);
    vB0->SetVisibility(kTRUE);
    vB0->SetLineColor(1); // Black
    vB0->SetLineWidth(1);
    vB0->SetFillColor(vB0->GetLineColor());
    vB0->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vB3 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockB3",sB3, medSPDcf);
    vB3->SetVisibility(kTRUE);
    vB3->SetLineColor(26); // Brown shade
    vB3->SetLineWidth(1);
    vB3->SetFillColor(vB3->GetLineColor());
    vB3->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vM3 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockScrewM3",sM3, medSPDss);
    vM3->SetVisibility(kTRUE);
    vM3->SetLineColor(kGray); // Gray
    vM3->SetLineWidth(1);
    vM3->SetFillColor(vM3->GetLineColor());
    vM3->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vD1 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockScrewHead1",sD1, medSPDss);
    vD1->SetVisibility(kTRUE);
    vD1->SetLineColor(kGray); // Gray
    vD1->SetLineWidth(1);
    vD1->SetFillColor(vD1->GetLineColor());
    vD1->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vD2 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockScrewHead2",sD2, medSPDss);
    vD2->SetVisibility(kTRUE);
    vD2->SetLineColor(kGray); // Gray
    vD2->SetLineWidth(1);
    vD2->SetFillColor(vD2->GetLineColor());
    vD2->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vN3 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockScrewNut",sN3, medSPDss);
    vN3->SetVisibility(kTRUE);
    vN3->SetLineColor(kGray); // Gray
    vN3->SetLineWidth(1);
    vN3->SetFillColor(vN3->GetLineColor());
    vN3->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vTB0 = new TGeoVolume("ITSSPDCoolingTubeEndTB0",sTB0,medSPDss);
    vTB0->SetVisibility(kTRUE);
    vTB0->SetLineColor(15); // gray
    vTB0->SetLineWidth(1);
    vTB0->SetFillColor(vTB0->GetLineColor());
    vTB0->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vTB1 = new TGeoVolume("ITSSPDCoolingTubeEndFluidTB1",sTB1,
                                      medSPDcoolfl);
    vTB1->SetVisibility(kTRUE);
    vTB1->SetLineColor(7); // light blue
    vTB1->SetLineWidth(1);
    vTB1->SetFillColor(vTB1->GetLineColor());
    vTB1->SetFillStyle(4050); // 0% transparent

    // add volumes to mother container passed as argument of this method
    moth->AddNode(vM0,1,0); // Add virtual volume to mother
    vTA0->AddNode(vTA1,1,0); // Put cooling liquid indide tube middel.
    vTB0->AddNode(vTB1,1,0); // Put cooling liquid inside tube end.
    Double_t tubeEndLocal[3]={0.0,0.0,sTA0->GetDz()};
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        x0 = secX3[ksecDipIndex[i]];
        y0 = secY3[ksecDipIndex[i]];
        t = 90.0 - secAngleTurbo[i];
	z0 = 0.5*(sB1->GetZ(0)+sB1->GetZ(1));
        trans = new TGeoTranslation("",x0,y0,z0);
        vM0->AddNode(vTB0, i+1, trans);
        // Find location of tube ends for later use.
        trans->LocalToMaster(tubeEndLocal,fTubeEndSector[0][0][i]);
        trans = new TGeoTranslation("",x0,y0,-z0);
        vM0->AddNode(vTB0, i+1+ksecNCoolingTubeDips, trans);
        rot = new TGeoRotation("", 0.0, 0.0, t);
        rotrans = new TGeoCombiTrans("", x0, y0, 0.0, rot);
        vM0->AddNode(vTA0, i+1, rotrans);
    } // end for i
    vM0->AddNode(vA0, 1, 0);
    vM0->AddNode(vB0, 1, 0);
    // Reflection.
    rot = new TGeoRotation("", 90., 0., 90., 90., 180., 0.);
    vM0->AddNode(vB0,2,rot);
    // Find location of tube ends for later use.
    for(i=0;i<ksecNCoolingTubeDips;i++) rot->LocalToMaster(
                            fTubeEndSector[0][0][i],fTubeEndSector[0][1][i]);
    // Put screws inside the mounting block
    const Double_t kMountingBlockScrew1ZPos =  0.7 *fgkcm;
    const Double_t kMountingBlockScrew2ZPos =  2.01*fgkcm;
    const Double_t kMountingBlockScrew34Pos =  0.51*fgkcm;
    vB3->AddNode(vM3, 1, new TGeoCombiTrans(0, 0,
				 (sB3->GetDZ()-kMountingBlockScrew1ZPos),
					    new TGeoRotation("",90,90,90)));
    vB3->AddNode(vM3, 2, new TGeoCombiTrans(0, 0,
				 (sB3->GetDZ()-kMountingBlockScrew2ZPos),
					    new TGeoRotation("",90,90,90)));
    vB3->AddNode(vM3, 3, new TGeoCombiTrans(0,-kMountingBlockScrew34Pos,
				-(sB3->GetDZ()-kMountingBlockScrew34Pos),
					    new TGeoRotation("",90,90,90)));
    vB3->AddNode(vM3, 4, new TGeoCombiTrans(0, kMountingBlockScrew34Pos,
				-(sB3->GetDZ()-kMountingBlockScrew34Pos),
					    new TGeoRotation("",90,90,90)));
    // left side
    t = -TMath::RadToDeg()*TMath::ATan2(
                                   sB0->GetX(0)-sB0->GetX(sB0->GetNvert()-1),
                                   sB0->GetY(0)-sB0->GetY(sB0->GetNvert()-1));
    rot = new TGeoRotation("",t,0.0,0.0);// z axis rotation
    x0 = 0.5*(sB0->GetX(0)+sB0->GetX(sB0->GetNvert()-1))+
        sB3->GetDX()*TMath::Cos(t*TMath::DegToRad());
    y0 = 0.5*(sB0->GetY(0)+sB0->GetY(sB0->GetNvert()-1))+
        sB3->GetDX()*TMath::Sin(t*TMath::DegToRad());
    z0 = sB0->GetZ(0)+sB3->GetDZ();
    rotrans = new TGeoCombiTrans("",x0,y0,z0,rot);
    vM0->AddNode(vB3,1,rotrans); // Put Mounting bracket on sector
    // the screw heads and nuts
    Double_t h = sM3->GetDz() + sD1->GetDz();
    Double_t zt = sB3->GetDZ()-kMountingBlockScrew1ZPos;
    vM0->AddNode(vD1, 1, new TGeoCombiTrans(x0+h*CosD(180+t), y0+h*SinD(180+t),
					    z0+zt,
					    new TGeoRotation("",90+t,90,90)));
    h = sM3->GetDz() + sD2->GetDz() + ksecCthick2 + 0.06;
    zt = sB3->GetDZ()-kMountingBlockScrew2ZPos;
    vM0->AddNode(vD2, 1, new TGeoCombiTrans(x0+h*CosD(180+t), y0+h*SinD(180+t),
					    z0+zt,
					    new TGeoRotation("",90+t,90,90)));
    Double_t loc[3],mas[3];
    loc[0]=0;
    loc[1]=-kMountingBlockScrew34Pos;
    loc[2]=-(sB3->GetDZ()-kMountingBlockScrew34Pos);
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vD2, 2, new TGeoCombiTrans(mas[0]+h*CosD(180+t),
					    mas[1]+h*SinD(180+t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vD2, 3, new TGeoCombiTrans(mas[0]+h*CosD(180+t),
					    mas[1]+h*SinD(180+t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));

    rot = new TGeoRotation("",t,180.0,0.0);// z & x axis rotation
    rotrans = new TGeoCombiTrans("",x0,y0,-z0,rot);
    vM0->AddNode(vB3,2,rotrans); // Put Mounting bracket on sector
    h = sM3->GetDz() + sN3->GetZ(1);
    zt = sB3->GetDZ()-kMountingBlockScrew1ZPos;
    vM0->AddNode(vN3, 1, new TGeoCombiTrans(x0+h*CosD(180+t), y0+h*SinD(180+t),
					   -z0-zt,
					    new TGeoRotation("",90+t,90,90)));
    h += ksecCthick2 + 0.06;
    zt = sB3->GetDZ()-kMountingBlockScrew2ZPos;
    vM0->AddNode(vN3, 2, new TGeoCombiTrans(x0+h*CosD(180+t), y0+h*SinD(180+t),
					   -z0-zt,
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=-kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vN3, 3, new TGeoCombiTrans(mas[0]+h*CosD(180+t),
					    mas[1]+h*SinD(180+t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vN3, 4, new TGeoCombiTrans(mas[0]+h*CosD(180+t),
					    mas[1]+h*SinD(180+t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));

    t *= -1.0;
    rot = new TGeoRotation("",t,0.0,0.0); // z axis rotation
    x0 = -0.5*(sB0->GetX(0)+sB0->GetX(sB0->GetNvert()-1))-3.5*
        sB3->GetDX()*TMath::Cos(t*TMath::DegToRad());
    y0 = 0.5*(sB0->GetY(0)+sB0->GetY(sB0->GetNvert()-1))-3.5*
        sB3->GetDX()*TMath::Sin(t*TMath::DegToRad());
    rotrans = new TGeoCombiTrans("",1.01*x0,y0,z0,rot);
    vM0->AddNode(vB3,3,rotrans); // Put Mounting bracket on sector
    h = sM3->GetDz() + sN3->GetZ(1);
    zt = sB3->GetDZ()-kMountingBlockScrew1ZPos;
    vM0->AddNode(vN3, 5, new TGeoCombiTrans(x0-h*CosD(180-t), y0+h*SinD(180-t),
					    z0+zt,
					    new TGeoRotation("",90+t,90,90)));
    h += ksecCthick2 + 0.02;
    zt = sB3->GetDZ()-kMountingBlockScrew2ZPos;
    vM0->AddNode(vN3, 6, new TGeoCombiTrans(x0-h*CosD(180-t), y0+h*SinD(180-t),
					    z0+zt,
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=-kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vN3, 7, new TGeoCombiTrans(mas[0]-h*CosD(180-t),
					    mas[1]+h*SinD(180-t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vN3, 8, new TGeoCombiTrans(mas[0]-h*CosD(180-t),
					    mas[1]+h*SinD(180-t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));

    rot = new TGeoRotation("",t,180.0,0.0); // z & x axis rotation
    rotrans = new TGeoCombiTrans("",1.01*x0,y0,-z0,rot);
    vM0->AddNode(vB3,4,rotrans); // Put Mounting bracket on sector
    h = sM3->GetDz() + sD1->GetDz();
    zt = sB3->GetDZ()-kMountingBlockScrew1ZPos;
    vM0->AddNode(vD1, 2, new TGeoCombiTrans(x0-h*CosD(180-t), y0+h*SinD(180-t),
					   -z0-zt,
					    new TGeoRotation("",90+t,90,90)));
    h = sM3->GetDz() + sD2->GetDz() + ksecCthick2 + 0.02;
    zt = sB3->GetDZ()-kMountingBlockScrew2ZPos;
    vM0->AddNode(vD2, 4, new TGeoCombiTrans(x0-h*CosD(180-t), y0+h*SinD(180-t),
					   -z0-zt,
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=-kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vD2, 5, new TGeoCombiTrans(mas[0]-h*CosD(180-t),
					    mas[1]+h*SinD(180-t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));
    loc[1]=kMountingBlockScrew34Pos;
    rotrans->LocalToMaster(loc,mas);
    vM0->AddNode(vD2, 6, new TGeoCombiTrans(mas[0]-h*CosD(180-t),
					    mas[1]+h*SinD(180-t),
					    mas[2],
					    new TGeoRotation("",90+t,90,90)));

    if(GetDebug(3)){
        vM0->PrintNodes();
        vA0->PrintNodes();
        vB0->PrintNodes();
        vB3->PrintNodes();
        vTA0->PrintNodes();
        vTA1->PrintNodes();
        vTB0->PrintNodes();
        vTB1->PrintNodes();
    } // end if(GetDebug(3))
}
//______________________________________________________________________
Bool_t AliITSv11GeometrySPD::CFHolePoints(Double_t s,Double_t r1,
                   Double_t r2,Double_t l,Double_t &x,Double_t &y) const
{
    //
    // Step along arck a distancs ds and compute boundry of
    // two holes (radius r1 and r2) a distance l apart (along
    // x-axis).
    // Inputs:
    //   Double_t s   fractional Distance along arcs [0-1]
    //                where 0-> alpha=beta=0, 1-> alpha=90 degrees.
    //   Double_t r1  radius at center circle
    //   Double_t r2  radius of displaced circle
    //   Double_t l   Distance displaced circle is displaces (x-axis)
    // Output:
    //   Double_t x   x coordinate along double circle.
    //   Double_t y   y coordinate along double circle.
    // Return:
    //   logical, kFALSE if an error
    //
    Double_t alpha,beta;
    Double_t ac,bc,scb,sca,t,alphac,betac; // at intersection of two circles

    x=y=0.0;
    ac = r1*r1-l*l-r2*r2;
    bc = 2.*l*r2;
    if(bc==0.0) {printf("bc=0 l=%e r2=%e\n",l,r2);return kFALSE;}
    betac = TMath::ACos(ac/bc);
    alphac = TMath::Sqrt((bc-ac)*(bc+ac))/(2.*l*r1);
    scb = r2*betac;
    sca = r1*alphac;
    t = r1*0.5*TMath::Pi() - sca + scb;
    if(s<= scb/t){
        beta = s*t/r2;
        x = r2*TMath::Cos(beta) + l;
        y = r2*TMath::Sin(beta);
        //printf("betac=%e scb=%e t=%e s=%e beta=%e x=%e y=%e\n",
        //       betac,scb,t,s,beta,x,y);
        return kTRUE;
    }else{
        beta = (s*t-scb+sca)/(r1*0.5*TMath::Pi());
        alpha = beta*0.5*TMath::Pi();
        x = r1*TMath::Cos(alpha);
        y = r1*TMath::Sin(alpha);
        //printf("alphac=%e sca=%e t=%e s=%e beta=%e alpha=%e x=%e y=%e\n",
        //       alphac,sca,t,s,beta,alpha,x,y);
        return kTRUE;
    } // end if
    return kFALSE;
}
//______________________________________________________________________
Bool_t AliITSv11GeometrySPD::GetSectorMountingPoints(Int_t index,Double_t &x0,
                              Double_t &y0, Double_t &x1, Double_t &y1) const
{
    //
    // Returns the edges of the straight borders in the SPD sector shape,
    // which are used to mount staves on them.
    // Coordinate system is that of the carbon fiber sector volume.
    // ---
    // Index numbering is as follows:
    //                         /5
    //                        /\/4
    //                      1\   \/3
    //                      0|___\/2
    // ---
    // Arguments [the ones passed by reference contain output values]:
    //    Int_t    index   --> location index according to above scheme [0-5]
    //    Double_t &x0     --> (by ref) x0 location or the ladder sector [cm]
    //    Double_t &y0     --> (by ref) y0 location of the ladder sector [cm]
    //    Double_t &x1     --> (by ref) x1 location or the ladder sector [cm]
    //    Double_t &y1     --> (by ref) y1 location of the ladder sector [cm]
    //    TGeoManager *mgr --> The TGeo builder
    // ---
    // The location is described by a line going from (x0, y0) to (x1, y1)
    // ---
    // Returns kTRUE if no problems encountered.
    // Returns kFALSE if a problem was encountered (e.g.: shape not found).
    //
    Int_t isize = fSPDsectorX0.GetSize();

    x0 = x1 = y0 = y1 = 0.0;
    if(index < 0 || index > isize) {
      AliError(Form("index = %d: allowed 0 --> %d", index, isize));
      return kFALSE;
    } // end if(index<0||index>isize)
    x0 = fSPDsectorX0[index];
    x1 = fSPDsectorX1[index];
    y0 = fSPDsectorY0[index];
    y1 = fSPDsectorY1[index];
    return kTRUE;
}
//______________________________________________________________________
void AliITSv11GeometrySPD::SPDsectorShape(Int_t n,const Double_t *xc,
                              const Double_t *yc,  const Double_t *r,
                              const Double_t *ths, const Double_t *the,
                      Int_t npr, Int_t &m, Double_t **xp, Double_t **yp) const
{
    //
    // Code to compute the points that make up the shape of the SPD
    // Carbon fiber support sections
    // Inputs:
    //   Int_t n        size of arrays xc,yc, and r.
    //   Double_t *xc   array of x values for radii centers.
    //   Double_t *yc   array of y values for radii centers.
    //   Double_t *r    array of signed radii values.
    //   Double_t *ths  array of starting angles [degrees].
    //   Double_t *the  array of ending angles [degrees].
    //   Int_t     npr  the number of lines segments to aproximate the arc.
    // Outputs (arguments passed by reference):
    //   Int_t       m    the number of enetries in the arrays *xp[npr+1]
    //                    and *yp[npr+1].
    //   Double_t **xp    array of x coordinate values of the line segments
    //                    which make up the SPD support sector shape.
    //   Double_t **yp    array of y coordinate values of the line segments
    //                    which make up the SPD support sector shape.
    //
    Int_t    i, k;
    Double_t t, t0, t1;

    m = n*(npr + 1);
    if(GetDebug(2)) {
        cout <<"  X    \t  Y  \t  R  \t  S  \t  E" << m << endl;
        for(i = 0; i < n; i++) {
            cout << "{"    << xc[i] << ", ";
            cout << yc[i]  << ", ";
            cout << r[i]   << ", ";
            cout << ths[i] << ", ";
            cout << the[i] << "}, " << endl;
        } // end for i
    } // end if(GetDebug(2))
    if (GetDebug(3)) cout << "Double_t sA0 = [" << n*(npr+1)+1<<"][";
    if (GetDebug(4)) cout << "3] {";
    else if(GetDebug(3)) cout <<"2] {";
    t0 = (Double_t)npr;
    for(i = 0; i < n; i++) {
        t1 = (the[i] - ths[i]) / t0;
        if(GetDebug(5)) cout << "t1 = " << t1 << endl;
        for(k = 0; k <= npr; k++) {
            t = ths[i] + ((Double_t)k) * t1;
            xp[i][k] = TMath::Abs(r[i]) * CosD(t) + xc[i];
            yp[i][k] = TMath::Abs(r[i]) * SinD(t) + yc[i];
            if(GetDebug(3)) {
                cout << "{" << xp[i][k] << "," << yp[i][k];
                if (GetDebug(4)) cout << "," << t;
                cout << "},";
            } // end if GetDebug
        } // end for k
        if(GetDebug(3)) cout << endl;
    } // end of i
    if(GetDebug(3)) cout << "{"  << xp[0][0] << ", " << yp[0][0];
    if(GetDebug(4)) cout << ","  << ths[0];
    if(GetDebug(3)) cout << "}}" << endl;
}

//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateLadder(Int_t layer,TArrayD &sizes,
                                               TGeoManager *mgr) const
{
    //
    // Creates the "ladder" = silicon sensor + 5 chips.
    // Returns a TGeoVolume containing the following components:
    //  - the sensor (TGeoBBox), whose name depends on the layer
    //  - 5 identical chips (TGeoBBox)
    //  - a guard ring around the sensor (subtraction of TGeoBBoxes),
    //    which is separated from the rest of sensor because it is not
    //    a sensitive part
    //  - bump bondings (TGeoBBox stripes for the whole width of the
    //    sensor, one per column).
    // ---
    // Arguments:
    //  1 - the owner layer (MUST be 1 or 2 or a fatal error is raised)
    //  2 - a TArrayD passed by reference, which will contain relevant
    //      dimensions related to this object:
    //      size[0] = 'thickness' (the smallest dimension)
    //      size[1] = 'length' (the direction along the ALICE Z axis)
    //      size[2] = 'width' (extension in the direction perp. to the
    //                         above ones)
    //  3 - the used TGeoManager

    // ** CRITICAL CHECK **
    // layer number can be ONLY 1 or 2
    if (layer != 1 && layer != 2) AliFatal("Layer number MUST be 1 or 2");

    // ** MEDIA **
    TGeoMedium *medAir       = GetMedium("AIR$",mgr);
    TGeoMedium *medSPDSiChip = GetMedium("SPD SI CHIP$",mgr); // SPD SI CHIP
    TGeoMedium *medSi        = GetMedium("SI$",mgr);
    TGeoMedium *medBumpBond  = GetMedium("COPPER$",mgr);  // ??? BumpBond

    // ** SIZES **
    Double_t chipThickness  = fgkmm *  0.150;
    Double_t chipWidth      = fgkmm * 15.950;
    Double_t chipLength     = fgkmm * 13.600;
    Double_t chipSpacing    = fgkmm *  0.400; // separation of chips along Z
    Double_t sensThickness  = fgkmm *  0.200;
    Double_t sensLength     = fgkmm * 69.600;
    Double_t sensWidth      = fgkmm * 12.800;
    Double_t guardRingWidth = fgkmm *  0.560; // a border of this thickness
                                              // all around the sensor
    Double_t bbLength       = fgkmm * 0.042;
    Double_t bbWidth        = sensWidth;
    Double_t bbThickness    = fgkmm * 0.012;
    Double_t bbPos          = 0.080;  // Z position w.r. to left pixel edge
    // compute the size of the container volume which
    // will also be returned in the referenced TArrayD;
    // for readability, they are linked by reference to a more meaningful name
    sizes.Set(3);
    Double_t &thickness = sizes[0];
    Double_t &length = sizes[1];
    Double_t &width = sizes[2];
    // the container is a box which exactly enclose all the stuff;
    width = chipWidth;
    length = sensLength + 2.0*guardRingWidth;
    thickness = sensThickness + chipThickness + bbThickness;

    // ** VOLUMES **
    // While creating this volume, since it is a sensitive volume,
    // we must respect some standard criteria for its local reference frame.
    // Local X must correspond to x coordinate of the sensitive volume:
    // this means that we are going to create the container with a local
    // reference system that is **not** in the middle of the box.
    // This is accomplished by calling the shape constructor with an
    // additional option ('originShift'):
    Double_t xSens = 0.5 * (width - sensWidth - 2.0*guardRingWidth);
    Double_t originShift[3] = {-xSens, 0., 0.};
    TGeoBBox *shapeContainer = new TGeoBBox(0.5*width,0.5*thickness,
                                            0.5*length,originShift);
    // then the volume is made of air, and using this shape
    TGeoVolume *container = new TGeoVolume(Form("ITSSPDlay%d-Ladder",layer),
                                           shapeContainer, medAir);
    // the chip is a common box
    TGeoVolume *volChip = mgr->MakeBox("ITSSPDchip",medSPDSiChip,
                              0.5*chipWidth,0.5*chipThickness,0.5*chipLength);
    // the sensor as well
    TGeoVolume *volSens = mgr->MakeBox(GetSenstiveVolumeName(layer),medSi,
                             0.5*sensWidth,0.5*sensThickness,0.5*sensLength);
    // the guard ring shape is the subtraction of two boxes with the
    // same center.
    TGeoBBox  *shIn = new TGeoBBox(0.5*sensWidth,sensThickness,0.5*sensLength);
    TGeoBBox  *shOut = new TGeoBBox(0.5*sensWidth+guardRingWidth,
                              0.5*sensThickness,0.5*sensLength+guardRingWidth);
    shIn->SetName("ITSSPDinnerBox");
    shOut->SetName("ITSSPDouterBox");
    TGeoCompositeShape *shBorder = new TGeoCompositeShape(
      "ITSSPDgaurdRingBorder",Form("%s-%s",shOut->GetName(),shIn->GetName()));
    TGeoVolume *volBorder = new TGeoVolume("ITSSPDgaurdRing",shBorder,medSi);
    // bump bonds for one whole column
    TGeoVolume *volBB = mgr->MakeBox("ITSSPDbb",medBumpBond,0.5*bbWidth,
                                     0.5*bbThickness,0.5*bbLength);
    // set colors of all objects for visualization
    volSens->SetLineColor(kYellow + 1);
    volChip->SetLineColor(kGreen);
    volBorder->SetLineColor(kYellow + 3);
    volBB->SetLineColor(kGray);

    // ** MOVEMENTS **
    // sensor is translated along thickness (X) and width (Y)
    Double_t ySens = 0.5 * (thickness - sensThickness);
    Double_t zSens = 0.0;
    // we want that the x of the ladder is the same as the one of
    // its sensitive volume
    TGeoTranslation *trSens = new TGeoTranslation(0.0, ySens, zSens);
    // bump bonds are translated along all axes:
    // keep same Y used for sensors, but change the Z
    TGeoTranslation *trBB[160];
    Double_t x =  0.0;
    Double_t y =  0.5 * (thickness - bbThickness) - sensThickness;
    Double_t z = -0.5 * sensLength + guardRingWidth + fgkmm*0.425 - bbPos;
    Int_t i;
    for (i = 0; i < 160; i++) {
        trBB[i] = new TGeoTranslation(x, y, z);
        switch(i) {
        case  31:case  63:case  95:case 127:
            z += fgkmm * 0.625 + fgkmm * 0.2;
            break;
        default:
            z += fgkmm * 0.425;
        } // end switch
    } // end for i
    // the chips are translated along the length (Z) and thickness (X)
    TGeoTranslation *trChip[5] = {0, 0, 0, 0, 0};
    x = -xSens;
    y = 0.5 * (chipThickness - thickness);
    z = 0.0;
    for (i = 0; i < 5; i++) {
        z = -0.5*length + guardRingWidth
            + (Double_t)i*chipSpacing + ((Double_t)(i) + 0.5)*chipLength;
        trChip[i] = new TGeoTranslation(x, y, z);
    } // end ofr i

    // add nodes to container
    container->AddNode(volSens, 1, trSens);
    container->AddNode(volBorder, 1, trSens);
    for (i = 0; i < 160; i++) container->AddNode(volBB,i+1,trBB[i]);
    for (i = 0; i < 5; i++) container->AddNode(volChip,i+3,trChip[i]);
    // return the container
    return container;
}

//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateClip(TArrayD &sizes,Bool_t isDummy,
                                             TGeoManager *mgr) const
{
    //
    // Creates the carbon fiber clips which are added to the central ladders.
    // They have a complicated shape which is approximated by a TGeoXtru
    // Implementation of a single clip over an half-stave.
    // It has a complicated shape which is approximated to a section like this:
    //
    //     6
    //     /\   .
    //  7 //\\  5
    //    / 1\\___________________4
    //   0    \___________________
    //        2                   3
    // with a finite thickness for all the shape
    // Its local reference frame is such that point A corresponds to origin.
    //

  // MODIFIED geometry
    Double_t sposty = fgkmm * -0.5; // lower internal side to avoid overlaps with modified geometry

    Double_t fullLength      = fgkmm * 12.6;    // = x4 - x0
    Double_t flatLength      = fgkmm *  5.4;    // = x4 - x3
    Double_t inclLongLength  = fgkmm *  5.0;    // = 5-6
    Double_t inclShortLength = fgkmm *  2.0;    // = 6-7
    Double_t fullHeight      = fgkmm *  2.8;    // = y6 - y3
    Double_t thickness       = fgkmm *  0.18;    // thickness
    Double_t totalLength     = fgkmm * 52.0;    // total length in Z
    Double_t holeSize        = fgkmm *  5.0;    // dimension of cubic
                                                // hole inserted for pt1000
    Double_t angle1          = 27.0;            // supplementary of angle DCB
    Double_t angle2;                            // angle DCB
    Double_t angle3;                            // angle of GH with vertical

    angle2 = 0.5 * (180.0 - angle1);
    angle3 = 90.0 - TMath::ACos(fullLength - flatLength -
                                inclLongLength*TMath::Cos(angle1)) *
                                TMath::RadToDeg();
    angle1 *= TMath::DegToRad();
    angle2 *= TMath::DegToRad();
    angle3 *= TMath::DegToRad();

    Double_t x[8], y[8];

    x[0] =  0.0;
    x[1] = x[0] + fullLength - flatLength - inclLongLength*TMath::Cos(angle1);
    x[2] = x[0] + fullLength - flatLength;
    x[3] = x[0] + fullLength;
    x[4] = x[3];
    x[5] = x[4] - flatLength + thickness * TMath::Cos(angle2);
    x[6] = x[1];
    x[7] = x[0];

    y[0] = 0.0;
    y[1] = y[0] + inclShortLength * TMath::Cos(angle3);
    y[2] = y[1] - inclLongLength * TMath::Sin(angle1);
    y[3] = y[2];
    y[4] = y[3] + thickness;
    y[5] = y[4];
    y[6] = y[1] + thickness;
    y[7] = y[0] + thickness;

    y[0] += sposty;
    y[7] += sposty;

    sizes.Set(7);
    sizes[0] = totalLength;
    sizes[1] = fullHeight;
    sizes[2] = y[2];
    sizes[3] = y[6];
    sizes[4] = x[0];
    sizes[5] = x[3];
    sizes[6] = x[2];

    if(isDummy){// use this argument when on ewant just the
                // positions without create any volume
        return NULL;
    } // end if isDummy

    TGeoXtru *shClip = new TGeoXtru(2);
    shClip->SetName("ITSSPDshclip");
    shClip->DefinePolygon(8, x, y);
    shClip->DefineSection(0, -0.5*totalLength, 0., 0., 1.0);
    shClip->DefineSection(1,  0.5*totalLength, 0., 0., 1.0);

    TGeoBBox *shHole = new TGeoBBox("ITSSPDSHClipHole",0.5*holeSize,
                                    0.5*holeSize,0.5*holeSize);
    TGeoTranslation *tr1 = new TGeoTranslation("ITSSPDTRClipHole1",x[2],0.0,
                                               fgkmm*14.);
    TGeoTranslation *tr2 = new TGeoTranslation("ITSSPDTRClipHole2",x[2],0.0,
                                               0.0);
    TGeoTranslation *tr3 = new TGeoTranslation("ITSSPDTRClipHole3",x[2],0.0,
                                               -fgkmm*14.);
    tr1->RegisterYourself();
    tr2->RegisterYourself();
    tr3->RegisterYourself();

    //TString strExpr("ITSSPDshclip-(");
    TString strExpr(shClip->GetName());
    strExpr.Append("-(");
    strExpr.Append(Form("%s:%s+", shHole->GetName(), tr1->GetName()));
    strExpr.Append(Form("%s:%s+", shHole->GetName(), tr2->GetName()));
    strExpr.Append(Form("%s:%s)", shHole->GetName(), tr3->GetName()));
    TGeoCompositeShape *shClipHole = new TGeoCompositeShape(
        "ITSSPDSHClipHoles",strExpr.Data());

    TGeoMedium *mat = GetMedium("SPD C (M55J)$", mgr);
    TGeoVolume *vClip = new TGeoVolume("ITSSPDclip", shClipHole, mat);
    vClip->SetLineColor(kGray + 2);
    return vClip;
}

//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreatePatchPanel(TArrayD &sizes,
						   TGeoManager *mgr) const
{
    //
    // Creates the patch panel approximated with a "L"-shaped TGeoXtru
    // with a finite thickness for all the shape
    // Its local reference frame is such that point A corresponds to origin.
    //
    Double_t hLength         = fgkmm *  50.0;    // horizontal length
    Double_t vLength         = fgkmm *  50.0;    // vertical length
    Double_t angle           = 88.3;             // angle between hor and vert
    Double_t thickness       = fgkmm *   4.0;    // thickness
    Double_t width           = fgkmm * 100.0;    // width looking from cone

    Double_t x[7], y[7];

    y[0] =  0.0;
    y[1] = y[0] + hLength;
    y[2] = y[1];
    y[3] = y[0] + thickness;
    y[4] = y[3] + vLength * TMath::Cos(angle*TMath::DegToRad());
    y[5] = y[4] - thickness / TMath::Sin(angle*TMath::DegToRad());
    y[6] = y[0];

    x[0] = 0.0;
    x[1] = x[0];
    x[2] = x[1] + thickness;
    x[3] = x[2];
    x[4] = x[3] + vLength * TMath::Sin(angle*TMath::DegToRad());
    x[5] = x[4];
    x[6] = x[0] + thickness;

    sizes.Set(3);
    sizes[0] = hLength;
    sizes[1] = vLength;
    sizes[2] = thickness;

    TGeoXtru *shPatch = new TGeoXtru(2);
    shPatch->SetName("ITSSPDpatchShape1");
    shPatch->DefinePolygon(7, x, y);
    shPatch->DefineSection(0, -0.5*width, 0., 0., 1.0);
    shPatch->DefineSection(1,  0.5*width, 0., 0., 1.0);
    
    /*
    Double_t subThickness = 10.0 * fgkmm;
    Double_t subWidth     = 55.0 * fgkmm;
    new TGeoBBox("ITSSPDpatchShape2", 0.5*subThickness, 60.0 * fgkmm, 0.5*subWidth);
    TGeoRotation *rotSub = new TGeoRotation(*gGeoIdentity);
    rotSub->SetName("shPatchSubRot");
    rotSub->RotateZ(50.0);
    rotSub->RegisterYourself();
    TGeoCombiTrans *trSub = new TGeoCombiTrans(0.26*hLength, 0.26*vLength, 0.0, rotSub);
    trSub->SetName("shPatchSubTr");
    trSub->RegisterYourself();
    
    TGeoCompositeShape *shPatchFinal = new TGeoCompositeShape("ITSSPDpatchShape1-(ITSSPDpatchShape2:shPatchSubTr)");
    */

    TGeoMedium *mat = GetMedium("AL$", mgr);
    //TGeoVolume *vPatch = new TGeoVolume("ITSSPDpatchPanel", shPatchFinal, mat);
    TGeoVolume *vPatch = new TGeoVolume("ITSSPDpatchPanel", shPatch, mat);
    vPatch->SetLineColor(kAzure);
    
    return vPatch;
}

//___________________________________________________________________
TGeoCompositeShape* AliITSv11GeometrySPD::CreateGroundingFoilShape
                       (Int_t itype,Double_t &length,Double_t &width,
                        Double_t thickness,TArrayD &sizes)
{
    //
    // Creates the typical composite shape of the grounding foil:
    //
    //  +---------------------------------------------------------+
    //  |                         5           6      9            |
    //  |                         +-----------+      +------------+ 10
    //  |             O           |           |      |
    //  |                 3 /-----+ 4         +------+
    //  |     1            /                 7        8
    //  |      /----------/
    //  +-----/                2                                  +
    //       0
    //       Z                                                    + 11
    //
    // This shape is used 4 times: two layers of glue, one in kapton
    // and one in aluminum, taking into account that the aliminum
    // layer has small differences in the size of some parts.
    // ---
    // In order to overcome problems apparently due to a large number
    // of points, the shape creation is done according the following
    // steps:
    //    1) a TGeoBBox is created with a size right enough to contain
    //       the whole shape (0-1-X-13)
    //    2) holes are defined as other TGeoBBox which are subtracted
    //       from the main shape
    //    3) a TGeoXtru is defined connecting the points (0-->11-->0)
    //       and is also subtracted from the main shape
    // ---
    // The argument ("type") is used to choose between all these
    // possibilities:
    //   - type = 0 --> kapton layer
    //   - type = 1 --> aluminum layer
    //   - type = 2 --> glue layer between support and GF
    //   - type = 3 --> glue layer between GF and ladders
    // Returns: a TGeoCompositeShape which will then be used to shape
    // several volumes. Since TGeoXtru is used, the local reference
    // frame of this object has X horizontal and Y vertical w.r to
    // the shape drawn above, and Z axis going perpendicularly to the screen.
    // This is not the correct reference for the half stave, for which
    // the "long" dimension is Z and the "short" is X, while Y goes in
    // the direction of thickness. This will imply some rotations when
    // using the volumes created with this shape.

    // suffix to differentiate names
    Char_t type[10];

    // size of the virtual box containing exactly this volume
    length = fgkmm * 243.18;
    width  = fgkmm *  15.95;
    if (itype == 1) {
        length -= fgkmm * 0.4;
        width  -= fgkmm * 0.4;
    } // end if itype==1
    switch (itype) {
    case 0:
        snprintf(type,10,"Kap");
        break;
    case 1:
        snprintf(type,10, "Alu");
        break;
    case 2:
        snprintf(type,10,"Glue1");
        break;
    case 3:
        snprintf(type,10,"Glue2");
        break;
    }
    // we divide the shape in several slices along the horizontal
    // direction (local X) here we define define the length of all
    // sectors (from leftmost to rightmost)
    Int_t i;
    Double_t sliceLength[] = { 140.71,  2.48,  26.78,   4.00,
                                10.00, 24.40,  10.00,  24.81 };
    for (i = 0; i < 8; i++) sliceLength[i] *= fgkmm;
    if (itype == 1) {
        sliceLength[0] -= fgkmm * 0.2;
        sliceLength[4] -= fgkmm * 0.2;
        sliceLength[5] += fgkmm * 0.4;
        sliceLength[6] -= fgkmm * 0.4;
    } // end if itype ==1

    // as shown in the drawing, we have four different widths
    // (along local Y) in this shape:
    Double_t widthMax  = fgkmm * 15.95;
    Double_t widthMed1 = fgkmm * 15.00;
    Double_t widthMed2 = fgkmm * 11.00;
    Double_t widthMin  = fgkmm *  4.40;
    if (itype == 1) {
        widthMax  -= fgkmm * 0.4;
        widthMed1 -= fgkmm * 0.4;
        widthMed2 -= fgkmm * 0.4;
        widthMin  -= fgkmm * 0.4;
    } // end if itype==1

    // create the main shape
    TGeoBBox *shGroundFull = 0;
    shGroundFull = new TGeoBBox(Form("ITSSPDSHgFoil%sFull", type),
                                0.5*length,0.5*width, 0.5*thickness);

    if(GetDebug(5)) shGroundFull->Print(); // Avoid Coverity warning

    // create the polygonal shape to be subtracted to give the correct
    // shape to the borders its vertices are defined in sugh a way that
    // this polygonal will be placed in the correct place considered
    // that the origin of the local reference frame is in the center
    // of the main box: we fix the starting point at the lower-left
    // edge of the shape (point 12), and add all points in order,
    // following a clockwise rotation

    Double_t x[13], y[13];
    x[ 0] = -0.5 * length + sliceLength[0];
    y[ 0] = -0.5 * widthMax;

    x[ 1] = x[0] + sliceLength[1];
    y[ 1] = y[0] + (widthMax - widthMed1);

    x[ 2] = x[1] + sliceLength[2];
    y[ 2] = y[1];

    x[ 3] = x[2] + sliceLength[3];
    y[ 3] = y[2] + (widthMed1 - widthMed2);

    x[ 4] = x[3] + sliceLength[4];
    y[ 4] = y[3];

    x[ 5] = x[4];
    y[ 5] = y[4] + (widthMed2 - widthMin);

    x[ 6] = x[5] + sliceLength[5];
    y[ 6] = y[5];

    x[ 7] = x[6];
    y[ 7] = y[4];

    x[ 8] = x[7] + sliceLength[6];
    y[ 8] = y[7];

    x[ 9] = x[8];
    y[ 9] = y[6];

    x[10] = x[9] + sliceLength[7] + 0.5;
    y[10] = y[9];

    x[11] = x[10];
    y[11] = y[0] - 0.5;

    x[12] = x[0];
    y[12] = y[11];

    // create the shape
    TGeoXtru *shGroundXtru = new TGeoXtru(2);
    shGroundXtru->SetName(Form("ITSSPDSHgFoil%sXtru", type));
    shGroundXtru->DefinePolygon(13, x, y);
    shGroundXtru->DefineSection(0, -thickness, 0., 0., 1.0);
    shGroundXtru->DefineSection(1,  thickness, 0., 0., 1.0);

    // define a string which will express the algebric operations among volumes
    // and add the subtraction of this shape from the main one
    TString strComposite(Form("ITSSPDSHgFoil%sFull-(%s+", type,
                              shGroundXtru->GetName()));

    // define the holes according to size information coming from drawings:
    Double_t holeLength = fgkmm * 10.00;
    Double_t holeWidth  = fgkmm *  7.50;
    Double_t holeSepX0  = fgkmm *  7.05;  // separation between center
                                          // of first hole and left border
    Double_t holeSepXC  = fgkmm * 14.00;  // separation between the centers
                                          // of two consecutive holes
    Double_t holeSepX1  = fgkmm * 15.42;  // separation between centers of
                                          // 5th and 6th hole
    Double_t holeSepX2  = fgkmm * 22.00;  // separation between centers of
                                          // 10th and 11th hole
    if (itype == 1) {
        holeSepX0  -= fgkmm * 0.2;
        holeLength += fgkmm * 0.4;
        holeWidth  += fgkmm * 0.4;
    } // end if itype==1
    sizes.Set(7);
    sizes[0] = holeLength;
    sizes[1] = holeWidth;
    sizes[2] = holeSepX0;
    sizes[3] = holeSepXC;
    sizes[4] = holeSepX1;
    sizes[5] = holeSepX2;
    sizes[6] = fgkmm * 4.40;

    // X position of hole center (will change for each hole)
    Double_t holeX = -0.5*length;
    // Y position of center of all holes (= 4.4 mm from upper border)
    Double_t holeY = 0.5*(width - holeWidth) - widthMin;

    // create a shape for the holes (common)
    new TGeoBBox(Form("ITSSPD%sGfoilHole", type),0.5*holeLength,
                       0.5*holeWidth, thickness);

    // insert the holes in the XTRU shape:
    // starting from the first value of X, they are simply
    // shifted along this axis
    char name[200];
    TGeoTranslation *transHole[11];
    for (i = 0; i < 11; i++) {
        // set the position of the hole, depending on index
        if (i == 0) {
            holeX += holeSepX0;
        }else if (i < 5) {
            holeX += holeSepXC;
        }else if (i == 5) {
            holeX += holeSepX1;
        }else if (i < 10) {
            holeX += holeSepXC;
        }else {
            holeX += holeSepX2;
        } // end if else if's
        //cout << i << " --> X = " << holeX << endl;
        snprintf(name,200,"ITSSPDTRgFoil%sHole%d", type, i);
        transHole[i] = new TGeoTranslation(name, holeX, holeY, 0.0);
        transHole[i]->RegisterYourself();
        strComposite.Append(Form("ITSSPD%sGfoilHole:%s", type, name));
        if (i < 10) strComposite.Append("+"); else strComposite.Append(")");
    } // end for i

    // create composite shape
    TGeoCompositeShape *shGround = new TGeoCompositeShape(
        Form("ITSSPDSHgFoil%s", type), strComposite.Data());

    return shGround;
}
//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreateGroundingFoil(Bool_t isRight,
                                   TArrayD &sizes, TGeoManager *mgr)
{
    //
    // Create a volume containing all parts of the grounding foil a
    // for a half-stave.
    // It consists of 4 layers with the same shape but different thickness:
    // 1) a layer of glue
    // 2) the aluminum layer
    // 3) the kapton layer
    // 4) another layer of glue
    // ---
    // Arguments:
    //  1: a boolean value to know if it is the grounding foir for
    //     the right or left side
    //  2: a TArrayD which will contain the dimension of the container box:
    //       - size[0] = length along Z (the beam line direction)
    //       - size[1] = the 'width' of the stave, which defines, together
    //                   with Z, the plane of the carbon fiber support
    //       - size[2] = 'thickness' (= the direction along which all
    //                    stave components are superimposed)
    //  3: the TGeoManager
    // ---
    // The return value is a TGeoBBox volume containing all grounding
    // foil components.
    // to avoid strange behaviour of the geometry manager,
    // create a suffix to be used in the names of all shapes
    //
    char suf[5];
    if (isRight) strncpy(suf, "R", 5); else strncpy(suf, "L", 5);
    // this volume will be created in order to ease its placement in
    // the half-stave; then, it is added here the small distance of
    // the "central" edge of each volume from the Z=0 plane in the stave
    // reference (which coincides with ALICE one)
    Double_t dist = fgkmm * 0.71;

    // define materials
    TGeoMedium *medKap  = GetMedium("SPD KAPTON(POLYCH2)$", mgr);
    TGeoMedium *medAlu  = GetMedium("AL$", mgr);
    TGeoMedium *medGlue = GetMedium("EPOXY$", mgr); //??? GLUE_GF_SUPPORT

    // compute the volume shapes (thicknesses change from one to the other)
    Double_t kpLength, kpWidth, alLength, alWidth;
    TArrayD  kpSize, alSize, glSize;
    Double_t kpThickness = fgkmm * 0.04;
    Double_t alThickness = fgkmm * 0.01;
//cout << "AL THICKNESS" << alThickness << endl;
    //Double_t g0Thickness = fgkmm * 0.1175 - fgkGapHalfStave;
    //Double_t g1Thickness = fgkmm * 0.1175 - fgkGapLadder;
    Double_t g0Thickness = fgkmm * 0.1275 - fgkGapHalfStave;
    Double_t g1Thickness = fgkmm * 0.1275 - fgkGapLadder;
    TGeoCompositeShape *kpShape = CreateGroundingFoilShape(0,kpLength,kpWidth,
                                                          kpThickness, kpSize);
    TGeoCompositeShape *alShape = CreateGroundingFoilShape(1,alLength,alWidth,
                                                          alThickness, alSize);
    TGeoCompositeShape *g0Shape = CreateGroundingFoilShape(2,kpLength,kpWidth,
                                                          g0Thickness, glSize);
    TGeoCompositeShape *g1Shape = CreateGroundingFoilShape(3,kpLength,kpWidth,
                                                          g1Thickness, glSize);
    // create the component volumes and register their sizes in the
    // passed arrays for readability reasons, some reference variables
    // explicit the meaning of the array slots
    TGeoVolume *kpVol = new TGeoVolume(Form("ITSSPDgFoilKap%s",suf),
                                       kpShape, medKap);
    TGeoVolume *alVol = new TGeoVolume(Form("ITSSPDgFoilAlu%s",suf),
                                       alShape, medAlu);
    TGeoVolume *g0Vol = new TGeoVolume(Form("ITSSPDgFoilGlue%s",suf),
                                       g0Shape, medGlue);
    TGeoVolume *g1Vol = new TGeoVolume(Form("ITSSPDgFoilGlue%s",suf),
                                       g1Shape, medGlue);
    // set colors for the volumes
    kpVol->SetLineColor(kRed);
    alVol->SetLineColor(kGray);
    g0Vol->SetLineColor(kYellow);
    g1Vol->SetLineColor(kYellow);
    // create references for the final size object
    if (sizes.GetSize() != 3) sizes.Set(3);
    Double_t &fullThickness = sizes[0];
    Double_t &fullLength = sizes[1];
    Double_t &fullWidth = sizes[2];
    // kapton leads the larger dimensions of the foil
    // (including the cited small distance from Z=0 stave reference plane)
    // the thickness is the sum of the ones of all components
    fullLength    = kpLength + dist;
    fullWidth     = kpWidth;
    fullThickness = kpThickness + alThickness + g0Thickness + g1Thickness;
    // create the container
//    TGeoMedium *air = GetMedium("AIR$", mgr);
    TGeoVolumeAssembly *container = new TGeoVolumeAssembly(Form("ITSSPDgFOIL-%s",suf));
//    TGeoVolume *container = mgr->MakeBox(Form("ITSSPDgFOIL-%s",suf),
//                 air, 0.5*fullThickness, 0.5*fullWidth, 0.5*fullLength);
    // create the common correction rotation (which depends of what side
    // we are building)
    TGeoRotation *rotCorr = new TGeoRotation(*gGeoIdentity);
    if (isRight) rotCorr->RotateY(90.0);
    else rotCorr->RotateY(-90.0);
    // compute the translations, which are in the length and
    // thickness directions
    Double_t x, y, z, shift = 0.0;
    if (isRight) shift = dist;
    // glue (bottom)
    x = -0.5*(fullThickness - g0Thickness);
    z =  0.5*(fullLength - kpLength) - shift;
    TGeoCombiTrans *glTrans0 = new TGeoCombiTrans(x, 0.0, z, rotCorr);
    // kapton
    x += 0.5*(g0Thickness + kpThickness);
    TGeoCombiTrans *kpTrans  = new TGeoCombiTrans(x, 0.0, z, rotCorr);
    // aluminum
    x += 0.5*(kpThickness + alThickness);
    z  = 0.5*(fullLength - alLength) - shift - 0.5*(kpLength - alLength);
    TGeoCombiTrans *alTrans  = new TGeoCombiTrans(x, 0.0, z, rotCorr);
    // glue (top)
    x += 0.5*(alThickness + g1Thickness);
    z  = 0.5*(fullLength - kpLength) - shift;
    TGeoCombiTrans *glTrans1 = new TGeoCombiTrans(x, 0.0, z, rotCorr);

    //cout << fgkGapHalfStave << endl;
    //cout << g0Thickness << endl;
    //cout << kpThickness << endl;
    //cout << alThickness << endl;
    //cout << g1Thickness << endl;

    // add to container
    container->SetLineColor(kMagenta-10);
    container->AddNode(kpVol, 1, kpTrans);
    container->AddNode(alVol, 1, alTrans);
    container->AddNode(g0Vol, 1, glTrans0);
    container->AddNode(g1Vol, 2, glTrans1);
    // to add the grease we remember the sizes of the holes, stored as
    // additional parameters in the kapton layer size:
    //   - sizes[3] = hole length
    //   - sizes[4] = hole width
    //   - sizes[5] = position of first hole center
    //   - sizes[6] = standard separation between holes
    //   - sizes[7] = separation between 5th and 6th hole
    //   - sizes[8] = separation between 10th and 11th hole
    //   - sizes[9] = separation between the upper hole border and
    //                the foil border
    Double_t holeLength      = kpSize[0];
    Double_t holeWidth       = kpSize[1];
    Double_t holeFirstZ      = kpSize[2];
    Double_t holeSepZ        = kpSize[3];
    Double_t holeSep5th6th   = kpSize[4];
    Double_t holeSep10th11th = kpSize[5];
    Double_t holeSepY        = kpSize[6];
    // volume (common)
    // Grease has not been defined to date. Need much more information
    // no this material!
    TGeoMedium *grease = GetMedium("SPD KAPTON(POLYCH2)$", mgr); // ??? GREASE
    TGeoVolume *hVol   = mgr->MakeBox("ITSSPDGrease", grease,
                           0.5*fullThickness, 0.5*holeWidth, 0.5*holeLength);
    hVol->SetLineColor(kBlue);
    // displacement of volumes in the container
    Int_t    idx = 1;  // copy numbers start from 1.
    x = 0.0;
    y = 0.5*(fullWidth - holeWidth) - holeSepY;
    if (isRight) z = holeFirstZ - 0.5*fullLength + dist;
    else z = 0.5*fullLength - holeFirstZ - dist;
    for (Int_t i = 0; i < 11; i++) {
        TGeoTranslation *t = 0;
        t = new TGeoTranslation(x, y, -z);
        container->AddNode(hVol, idx++, t);
        if (i < 4) shift = holeSepZ;
        else if (i == 4) shift = holeSep5th6th;
        else if (i < 9) shift = holeSepZ;
        else shift = holeSep10th11th;
        if (isRight) z += shift;
        else z -= shift;
    } // end for i
    return container;
}
//___________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreateMCM(Bool_t isRight,
                                   TArrayD &sizes, TGeoManager *mgr) const
{
    //
    // Create a TGeoAssembly containing all the components of the MCM.
    // The TGeoVolume container is rejected due to the possibility of overlaps
    // when placing this object on the carbon fiber sector.
    // The assembly contains:
    //  - the thin part of the MCM (integrated circuit)
    //  - the MCM chips (specifications from EDMS)
    //  - the cap which covers the zone where chips are bound to MCM
    // ---
    // The local reference frame of this assembly is defined in such a way
    // that all volumes are contained in a virtual box whose center
    // is placed exactly in the middle of the occupied space w.r to all
    // directions. This will ease the positioning of this object in the
    // half-stave. The sizes of this virtual box are stored in
    // the array passed by reference.
    // ---
    // Arguments:
    //  - a boolean flag to know if this is the "left" or "right" MCM, when
    //    looking at the stave from above (i.e. the direction from which
    //    one sees bus over ladders over grounding foil) and keeping the
    //    continuous border in the upper part, one sees the thicker part
    //    on the left or right.
    //  - an array passed by reference which will contain the size of
    //    the virtual container.
    //  - a pointer to the used TGeoManager.
    //

    // to distinguish the "left" and "right" objects, a suffix is created
    char suf[5];
    if (isRight) strncpy(suf, "R", 5); else strncpy(suf, "L", 5);

    // ** MEDIA **
    TGeoMedium *medBase = GetMedium("SPD KAPTON(POLYCH2)$",mgr);// ??? MCM BASE
    TGeoMedium *medChip = GetMedium("SPD SI CHIP$",mgr);
    TGeoMedium *medCap  = GetMedium("AL$",mgr);

    // The shape of the MCM is divided into 3 sectors with different
    // widths (Y) and lengths (X), like in this sketch:
    //
    //   0                      1                                   2
    //    +---------------------+-----------------------------------+
    //    |                                    4       sect 2       |
    //    |                    6      sect 1    /-------------------+
    //    |      sect 0         /--------------/                    3
    //    +--------------------/               5
    //   8                     7
    //
    // the inclination of all oblique borders (6-7, 4-5) is always 45 degrees.
    // From drawings we can parametrize the dimensions of all these sectors,
    // then the shape of this part of the MCM is implemented as a
    // TGeoXtru centerd in the virtual XY space.
    // The first step is definig the relevant sizes of this shape:
    Int_t i, j;
    Double_t mcmThickness  = fgkmm * 0.35;
    Double_t sizeXtot      = fgkmm * 105.6;   // total distance (0-2)
    // resp. 7-8, 5-6 and 3-4
    Double_t sizeXsector[3] = {fgkmm * 28.4, fgkmm * 41.4, fgkmm * 28.8};
    // resp. 0-8, 1-6 and 2-3
    Double_t sizeYsector[3] = {fgkmm * 15.0, fgkmm * 11.0, fgkmm *  8.0};
    Double_t sizeSep01 = fgkmm * 4.0;      // x(6)-x(7)
    Double_t sizeSep12 = fgkmm * 3.0;      // x(4)-x(5)

    // define sizes of chips (last is the thickest)
    Double_t chipLength[5]     = { 4.00, 6.15, 3.85, 5.60, 18.00 };
    Double_t chipWidth[5]      = { 3.00, 4.10, 3.85, 5.60,  5.45 };
    Double_t chipThickness[5]  = { 0.60, 0.30, 0.30, 1.00,  1.20 };
    TString  name[5];
    name[0] = "ITSSPDanalog";
    name[1] = "ITSSPDpilot";
    name[2] = "ITSSPDgol";
    name[3] = "ITSSPDrx40";
    name[4] = "ITSSPDoptical";
    Color_t color[5] = { kCyan, kGreen, kYellow, kBlue, kOrange };

    // define the sizes of the cover
    Double_t capThickness = fgkmm * 0.3;
    Double_t capHeight = fgkmm * 1.7;

    // compute the total size of the virtual container box
    sizes.Set(3);
    Double_t &thickness = sizes[0];
    Double_t &length = sizes[1];
    Double_t &width = sizes[2];
    length = sizeXtot;
    width = sizeYsector[0];
    thickness = mcmThickness + capHeight;

    // define all the relevant vertices of the polygon
    // which defines the transverse shape of the MCM.
    // These values are used to several purposes, and
    // for each one, some points must be excluded
    Double_t xRef[9], yRef[9];
    xRef[0] = -0.5*sizeXtot;
    yRef[0] =  0.5*sizeYsector[0];
    xRef[1] =  xRef[0] + sizeXsector[0] + sizeSep01;
    yRef[1] =  yRef[0];
    xRef[2] = -xRef[0];
    yRef[2] =  yRef[0];
    xRef[3] =  xRef[2];
    yRef[3] =  yRef[2] - sizeYsector[2];
    xRef[4] =  xRef[3] - sizeXsector[2];
    yRef[4] =  yRef[3];
    xRef[5] =  xRef[4] - sizeSep12;
    yRef[5] =  yRef[4] - sizeSep12;
    xRef[6] =  xRef[5] - sizeXsector[1];
    yRef[6] =  yRef[5];
    xRef[7] =  xRef[6] - sizeSep01;
    yRef[7] =  yRef[6] - sizeSep01;
    xRef[8] =  xRef[0];
    yRef[8] = -yRef[0];

    // the above points are defined for the "right" MCM (if ve view the
    // stave from above) in order to change to the "left" one, we must
    // change the sign to all X values:
    if (isRight) for (i = 0; i < 9; i++) xRef[i] = -xRef[i];

    // the shape of the MCM and glue layer are done excluding point 1,
    // which is not necessary and cause the geometry builder to get confused
    j = 0;
    Double_t xBase[8], yBase[8];
    for (i = 0; i < 9; i++) {
        if (i == 1) continue;
        xBase[j] = xRef[i];
        yBase[j] = yRef[i];
        j++;
    } // end for i

    // the MCM cover is superimposed over the zones 1 and 2 only
    Double_t xCap[6], yCap[6];
    j = 0;
    for (i = 1; i <= 6; i++) {
        xCap[j] = xRef[i];
        yCap[j] = yRef[i];
        j++;
    } // end for i

    // define positions of chips,
    // which must be added to the bottom-left corner of MCM
    // and divided by 1E4;
    Double_t chipX[5], chipY[5];
    if (isRight) {
        chipX[0] = 666320.;
        chipX[1] = 508320.;
        chipX[2] = 381320.;
        chipX[3] = 295320.;
        chipX[4] = 150320.;
        chipY[0] =  23750.;
        chipY[1] =  27750.;
        chipY[2] =  20750.;
        chipY[3] =  42750.;
        chipY[4] =  39750.;
    } else {
        chipX[0] = 389730.;
        chipX[1] = 548630.;
        chipX[2] = 674930.;
        chipX[3] = 761430.;
        chipX[4] = 905430.;
        chipY[0] =  96250.;
        chipY[1] =  91950.;
        chipY[2] =  99250.;
        chipY[3] = 107250.;
        chipY[4] = 109750.;
    } // end if isRight
    for (i = 0; i < 5; i++) {
        chipX[i] *= 0.00001;
        chipY[i] *= 0.00001;
        if (isRight) {
            chipX[i] += xRef[3];
            chipY[i] += yRef[3];
        } else {
            chipX[i] += xRef[8];
            chipY[i] += yRef[8];
        } // end for isRight
        chipLength[i] *= fgkmm;
        chipWidth[i] *= fgkmm;
        chipThickness[i] *= fgkmm;
    } // end for i

    // create shapes for MCM
    Double_t z1, z2;
    TGeoXtru *shBase = new TGeoXtru(2);
    z1 = -0.5*thickness;
    z2 = z1 + mcmThickness;
    shBase->DefinePolygon(8, xBase, yBase);
    shBase->DefineSection(0, z1, 0., 0., 1.0);
    shBase->DefineSection(1, z2, 0., 0., 1.0);

    // create volumes of MCM
    TGeoVolume *volBase = new TGeoVolume("ITSSPDbase", shBase, medBase);
    volBase->SetLineColor(kRed);

    // to create the border of the MCM cover, it is required the
    // subtraction of two shapes the outer is created using the
    // reference points defined here
    TGeoXtru *shCapOut = new TGeoXtru(2);
    shCapOut->SetName(Form("ITSSPDshCAPOUT%s", suf));
    z1 = z2;
    z2 = z1 + capHeight - capThickness;
    shCapOut->DefinePolygon(6, xCap, yCap);
    shCapOut->DefineSection(0, z1, 0., 0., 1.0);
    shCapOut->DefineSection(1, z2, 0., 0., 1.0);
    // the inner is built similarly but subtracting the thickness
    Double_t angle, cs;
    Double_t xin[6], yin[6];
    if (!isRight) {
        angle = 45.0;
        cs = TMath::Cos( 0.5*(TMath::Pi() - angle*TMath::DegToRad()) );
        xin[0] = xCap[0] + capThickness;
        yin[0] = yCap[0] - capThickness;
        xin[1] = xCap[1] - capThickness;
        yin[1] = yin[0];
        xin[2] = xin[1];
        yin[2] = yCap[2] + capThickness;
        xin[3] = xCap[3] - capThickness*cs;
        yin[3] = yin[2];
        xin[4] = xin[3] - sizeSep12;
        yin[4] = yCap[4] + capThickness;
        xin[5] = xin[0];
        yin[5] = yin[4];
    } else {
        angle = 45.0;
        cs = TMath::Cos( 0.5*(TMath::Pi() - angle*TMath::DegToRad()) );
        xin[0] = xCap[0] - capThickness;
        yin[0] = yCap[0] - capThickness;
        xin[1] = xCap[1] + capThickness;
        yin[1] = yin[0];
        xin[2] = xin[1];
        yin[2] = yCap[2] + capThickness;
        xin[3] = xCap[3] - capThickness*cs;
        yin[3] = yin[2];
        xin[4] = xin[3] + sizeSep12;
        yin[4] = yCap[4] + capThickness;
        xin[5] = xin[0];
        yin[5] = yin[4];
    } // end if !isRight
    TGeoXtru *shCapIn = new TGeoXtru(2);
    shCapIn->SetName(Form("ITSSPDshCAPIN%s", suf));
    shCapIn->DefinePolygon(6, xin, yin);
    shCapIn->DefineSection(0, z1 - 0.01, 0., 0., 1.0);
    shCapIn->DefineSection(1, z2 + 0.01, 0., 0., 1.0);
    // compose shapes
    TGeoCompositeShape *shCapBorder = new TGeoCompositeShape(
                            Form("ITSSPDshBORDER%s", suf),
                            Form("%s-%s", shCapOut->GetName(),
                                 shCapIn->GetName()));
    // create volume
    TGeoVolume *volCapBorder = new TGeoVolume("ITSSPDcapBoarder",
                                              shCapBorder,medCap);
    volCapBorder->SetLineColor(kGreen);
    // finally, we create the top of the cover, which has the same
    // shape of outer border and a thickness equal of the one othe
    // cover border one
    TGeoXtru *shCapTop = new TGeoXtru(2);
    z1 = z2;
    z2 = z1 + capThickness;
    shCapTop->DefinePolygon(6, xCap, yCap);
    shCapTop->DefineSection(0, z1, 0., 0., 1.0);
    shCapTop->DefineSection(1, z2, 0., 0., 1.0);
    TGeoVolume *volCapTop = new TGeoVolume("ITSSPDcapTop", shCapTop, medCap);
    volCapTop->SetLineColor(kBlue);

    // create container assembly with right suffix
    TGeoVolumeAssembly *mcmAssembly = new TGeoVolumeAssembly(
        Form("ITSSPDmcm%s", suf));

    // add mcm layer
    mcmAssembly->AddNode(volBase, 1, gGeoIdentity);
    // add chips
    for (i = 0; i < 5; i++) {
        TGeoVolume *box = gGeoManager->MakeBox(name[i],medChip,
               0.5*chipLength[i], 0.5*chipWidth[i], 0.5*chipThickness[i]);
        TGeoTranslation *tr = new TGeoTranslation(chipX[i],chipY[i],
                      0.5*(-thickness + chipThickness[i]) + mcmThickness);
        box->SetLineColor(color[i]);
        mcmAssembly->AddNode(box, 1, tr);
    } // end for i
    // add cap border
    mcmAssembly->AddNode(volCapBorder, 1, gGeoIdentity);
    // add cap top
    mcmAssembly->AddNode(volCapTop, 1, gGeoIdentity);

    return mcmAssembly;
}

//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreatePixelBus
(Bool_t isRight, Int_t ilayer, TArrayD &sizes, TGeoManager *mgr) const
{
    //
    // The pixel bus is implemented as a TGeoBBox with some objects on it,
    // which could affect the particle energy loss.
    // ---
    // In order to avoid confusion, the bus is directly displaced
    // according to the axis orientations which are used in the final stave:
    // X --> thickness direction
    // Y --> width direction
    // Z --> length direction
    //

    // ** CRITICAL CHECK ******************************************************
    // layer number can be ONLY 1 or 2
    if (ilayer != 1 && ilayer != 2) AliFatal("Layer number MUST be 1 or 2");

    // ** MEDIA **
    //PIXEL BUS
    TGeoMedium *medBus     = GetMedium("SPDBUS(AL+KPT+EPOX)$",mgr);
    TGeoMedium *medPt1000  = GetMedium("CERAMICS$",mgr); // ??? PT1000
    // Capacity
    TGeoMedium *medCap     = GetMedium("SDD X7R capacitors$",mgr);
    // ??? Resistance
    //TGeoMedium *medRes     = GetMedium("SDD X7R capacitors$",mgr);
    TGeoMedium *medRes     = GetMedium("ALUMINUM$",mgr);
    //TGeoMedium *medExt     = GetMedium("SDDKAPTON (POLYCH2)$", mgr);
    TGeoMedium *medExt     = GetMedium("SPD-MIX CU KAPTON$", mgr);
    // ** SIZES & POSITIONS **
    Double_t busLength          = 170.501 * fgkmm; // length of plane part
    Double_t busWidth           =  13.800 * fgkmm; // width
    Double_t busThickness       =   0.280 * fgkmm; // thickness
    Double_t pt1000Length       = fgkmm * 1.50;
    Double_t pt1000Width        = fgkmm * 3.10;
    Double_t pt1000Thickness    = fgkmm * 0.60;
    Double_t pt1000Y, pt1000Z[10];// position of the pt1000's along the bus
    Double_t capLength          = fgkmm * 2.55;
    Double_t capWidth           = fgkmm * 1.50;
    Double_t capThickness       = fgkmm * 1.35;
    Double_t capY[2], capZ[2];

    Double_t resLength          = fgkmm * 2.20;
    Double_t resWidth           = fgkmm * 0.80;
    Double_t resThickness       = fgkmm * 0.35;
    Double_t resY[2], resZ[2];

    Double_t extThickness       = fgkmm * 0.25;
    Double_t ext1Length         = fgkmm * (26.7 - 10.0);
    Double_t ext2Length         = fgkmm * 284.0 - ext1Length + extThickness;
    Double_t ext2LengthL2       = fgkmm * 130.0;
    Double_t ext4Length         = fgkmm * 40.0;
    Double_t ext4Twist          =  66.54; //deg
    Double_t extWidth           = fgkmm * 11.0;
    Double_t extHeight          = fgkmm * 2.5;

    // position of pt1000, resistors and capacitors depends on the
    // bus if it's left or right one
    if (!isRight) {
        pt1000Y    =   64400.;
        pt1000Z[0] =   66160.;
        pt1000Z[1] =  206200.;
        pt1000Z[2] =  346200.;
        pt1000Z[3] =  486200.;
        pt1000Z[4] =  626200.;
        pt1000Z[5] =  776200.;
        pt1000Z[6] =  916200.;
        pt1000Z[7] = 1056200.;
        pt1000Z[8] = 1196200.;
        pt1000Z[9] = 1336200.;
        resZ[0]    = 1397500.;
        resY[0]    =   26900.;
        resZ[1]    =  682500.;
        resY[1]    =   27800.;
        capZ[0]    = 1395700.;
        capY[0]    =   45700.;
        capZ[1]    =  692600.;
        capY[1]    =   45400.;
    } else {
        pt1000Y    =   66100.;
        pt1000Z[0] =  319700.;
        pt1000Z[1] =  459700.;
        pt1000Z[2] =  599700.;
        pt1000Z[3] =  739700.;
        pt1000Z[4] =  879700.;
        pt1000Z[5] = 1029700.;
        pt1000Z[6] = 1169700.;
        pt1000Z[7] = 1309700.;
        pt1000Z[8] = 1449700.;
        pt1000Z[9] = 1589700.;
        capY[0]    =   44500.;
        capZ[0]    =  266700.;
        capY[1]    =   44300.;
        capZ[1]    =  974700.;
        resZ[0]    =  266500.;
        resY[0]    =   29200.;
        resZ[1]    =  974600.;
        resY[1]    =   29900.;
    } // end if isRight
    Int_t i;
    pt1000Y *= 1E-4 * fgkmm;
    for (i = 0; i < 10; i++) {
        pt1000Z[i] *= 1E-4 * fgkmm;
        if (i < 2) {
            capZ[i] *= 1E-4 * fgkmm;
            capY[i] *= 1E-4 * fgkmm;
            resZ[i] *= 1E-4 * fgkmm;
            resY[i] *= 1E-4 * fgkmm;
        }  // end if iM2
    } // end for i

    Double_t &fullLength = sizes[1];
    Double_t &fullWidth = sizes[2];
    Double_t &fullThickness = sizes[0];
    fullLength = busLength;
    fullWidth = busWidth;
    // add the thickness of the thickest component on bus (capacity)
    fullThickness = busThickness + capThickness;

    // ** VOLUMES **
    TGeoVolumeAssembly *container = new TGeoVolumeAssembly("ITSSPDpixelBus");
    TGeoVolume *bus = mgr->MakeBox("ITSSPDbus", medBus, 0.5*busThickness,
                                   0.5*busWidth, 0.5*busLength);
    TGeoVolume *pt1000 = mgr->MakeBox("ITSSPDpt1000",medPt1000,
                        0.5*pt1000Thickness,0.5*pt1000Width, 0.5*pt1000Length);
    TGeoVolume *res = mgr->MakeBox("ITSSPDresistor", medRes, 0.5*resThickness,
                                   0.5*resWidth, 0.5*resLength);
    TGeoVolume *cap = mgr->MakeBox("ITSSPDcapacitor", medCap, 0.5*capThickness,
                                   0.5*capWidth, 0.5*capLength);

    char extname[12];
    snprintf(extname,12,"Extender1l%d",ilayer);
    TGeoVolume *ext1 = mgr->MakeBox(extname, medExt, 0.5*extThickness, 0.5*extWidth, 0.5*ext1Length);
    snprintf(extname,12,"Extender2l%d",ilayer);
    TGeoVolume *ext2 = mgr->MakeBox(extname, medExt, 0.5*extHeight - 2.*extThickness, 0.5*extWidth, 0.5*extThickness);
    TGeoVolume *ext3=0;
    snprintf(extname,12,"Extender3l%d",ilayer);
    TGeoVolume *ext4=0;
    snprintf(extname,12,"Extender3l%d",ilayer);
    if (ilayer==1) {
      Double_t halflen=(0.5*ext2Length + extThickness);
      Double_t xprof[6],yprof[6];
      Double_t alpha=24;
      xprof[0] = -halflen;
      yprof[0] = -0.5*extThickness;
      xprof[1] = halflen/2;
      yprof[1] = yprof[0];
      xprof[2] = xprof[1] + 0.5*halflen*CosD(alpha);
      yprof[2] = yprof[1] + 0.5*halflen*SinD(alpha);
      xprof[3] = xprof[2] - extThickness*SinD(alpha);
      yprof[3] = yprof[2] + extThickness*CosD(alpha);
      InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
		  extThickness, xprof[4], yprof[4]);
      xprof[5] = xprof[0];
      yprof[5] = 0.5*extThickness;
      TGeoXtru *ext3sh = new TGeoXtru(2);
      ext3sh->DefinePolygon(6, xprof, yprof);
      ext3sh->DefineSection(0, -0.5*(extWidth-0.8*fgkmm));
      ext3sh->DefineSection(1,  0.5*(extWidth-0.8*fgkmm));
      ext3 = new TGeoVolume(extname, ext3sh, medExt);
    } else {
      ext3 = mgr->MakeBox(extname, medExt, 0.5*extThickness, 0.5*(extWidth-0.8*fgkmm), 0.5*ext2LengthL2 + extThickness); // Hardcode fix of a small overlap
      ext4= mgr->MakeGtra("Extender4l2", medExt, 0.5*ext4Length, 0, 0, ext4Twist, 0.5*(extWidth-0.8*fgkmm), 0.5*extThickness, 0.5*extThickness, 0, 0.5*(extWidth-0.8*fgkmm), 0.5*extThickness, 0.5*extThickness, 0);
      ext4->SetLineColor(kGray);
    }
    bus->SetLineColor(kYellow + 2);
    pt1000->SetLineColor(kGreen + 3);
    res->SetLineColor(kRed + 1);
    cap->SetLineColor(kBlue - 7);
    ext1->SetLineColor(kGray);
    ext2->SetLineColor(kGray);
    ext3->SetLineColor(kGray);

    // ** MOVEMENTS AND POSITIONEMENT **
    // bus
    TGeoTranslation *trBus = new TGeoTranslation(0.5 * (busThickness -
                                                   fullThickness), 0.0, 0.0);
    container->AddNode(bus, 1, trBus);
    Double_t zRef, yRef, x, y, z;
    if (isRight) {
        zRef = -0.5*fullLength;
        yRef = -0.5*fullWidth;
    } else {
        zRef = -0.5*fullLength;
        yRef = -0.5*fullWidth;
    } // end if isRight
    // pt1000
    x = 0.5*(pt1000Thickness - fullThickness) + busThickness;
    for (i = 0; i < 10; i++) {
        y = yRef + pt1000Y;
        z = zRef + pt1000Z[i];
        TGeoTranslation *tr = new TGeoTranslation(x, y, z);
        container->AddNode(pt1000, i+1, tr);
    } // end for i
    // capacitors
    x = 0.5*(capThickness - fullThickness) + busThickness;
    for (i = 0; i < 2; i++) {
        y = yRef + capY[i];
        z = zRef + capZ[i];
        TGeoTranslation *tr = new TGeoTranslation(x, y, z);
        container->AddNode(cap, i+1, tr);
    } // end for i
    // resistors
    x = 0.5*(resThickness - fullThickness) + busThickness;
    for (i = 0; i < 2; i++) {
        y = yRef + resY[i];
        z = zRef + resZ[i];
        TGeoTranslation *tr = new TGeoTranslation(x, y, z);
        container->AddNode(res, i+1, tr);
    } // end for i

    // extender
        if (ilayer == 2) {
       if (isRight) {
          y = 0.5 * (fullWidth - extWidth) - 0.1;
          z = 0.5 * (-fullLength + fgkmm * 10.0);
       }
       else {
          y = 0.5 * (fullWidth - extWidth) - 0.1;
          z = 0.5 * ( fullLength - fgkmm * 10.0);
       }
        }
        else {
            if (isRight) {
                y = -0.5 * (fullWidth - extWidth);
                z = 0.5 * (-fullLength + fgkmm * 10.0);
            }
            else {
                y = -0.5 * (fullWidth - extWidth);
                z = 0.5 * ( fullLength - fgkmm * 10.0);
            }
        }
    x = 0.5 * (extThickness - fullThickness) + busThickness;
    //y = 0.5 * (fullWidth - extWidth);
    TGeoTranslation *trExt1 = new TGeoTranslation(x, y, z);
    if (isRight) {
        z -= 0.5 * (ext1Length - extThickness);
    }
    else {
        z += 0.5 * (ext1Length - extThickness);
    }
    x += 0.5*(extHeight - 3.*extThickness);
    TGeoTranslation *trExt2 = new TGeoTranslation(x, y, z);
    if (isRight) {
      if (ilayer==1)
        z -= 0.5 * (ext2Length - extThickness) + 2.5*extThickness;
      else
        z -= 0.5 * (ext2LengthL2 - extThickness) + 2.5*extThickness;
    }
    else {
      if (ilayer==1)
        z += 0.5 * (ext2Length - extThickness) + 2.5*extThickness;
      else
        z += 0.5 * (ext2LengthL2 - extThickness) + 2.5*extThickness;
    }
    x += 0.5*(extHeight - extThickness) - 2.*extThickness;
    TGeoCombiTrans *trExt3=0;
    if (ilayer==1) {
      if (isRight)
	trExt3 = new TGeoCombiTrans(x, y, z, new TGeoRotation("",0.,-90.,90.));
      else
	trExt3 = new TGeoCombiTrans(x, y, z, new TGeoRotation("",0., 90.,90.));
    } else
      trExt3 = new TGeoCombiTrans(x, y, z, 0);
    container->AddNode(ext1, 0, trExt1);
    container->AddNode(ext2, 0, trExt2);
    container->AddNode(ext3, 0, trExt3);
    if (ilayer==2) {
      TGeoCombiTrans *trExt4=0;
      if (isRight) {
	z -= ( ((TGeoBBox*)ext3->GetShape())->GetDZ() + ((TGeoGtra*)ext4->GetShape())->GetDZ() );
	trExt4 = new TGeoCombiTrans(x, y, z, new TGeoRotation("", ext4Twist/2,0,0));
      } else {
	z += ( ((TGeoBBox*)ext3->GetShape())->GetDZ() + ((TGeoGtra*)ext4->GetShape())->GetDZ() );
	trExt4 = new TGeoCombiTrans(x, y, z, new TGeoRotation("",-ext4Twist/2,0,0));
      }
      container->AddNode(ext4, 0, trExt4);
    }
    sizes[3] = yRef + pt1000Y;
    sizes[4] = zRef + pt1000Z[2];
    sizes[5] = zRef + pt1000Z[7];

    return container;
}

//______________________________________________________________________
TList* AliITSv11GeometrySPD::CreateConeModule(Bool_t sideC, const Double_t angrot,
					      TGeoManager *mgr) const
{
    //
    // Creates all services modules and places them in a TList
    // angrot is the rotation angle (passed as an argument to avoid
    // defining the same quantity in two different places)
    //
    // Created:      ?? ??? 2008  A. Pulvirenti
    // Updated:      03 May 2010  M. Sitta
    // Updated:      20 Jun 2010  A. Pulvirenti  Optical patch panels
    // Updated:      22 Jun 2010  M. Sitta  Fiber cables
    // Updated:      04 Jul 2010  M. Sitta  Water cooling
    // Updated:      08 Jul 2010  A. Pulvirenti  Air cooling on Side C
    //

    TGeoMedium *medInox  = GetMedium("INOX$",mgr);
    //TGeoMedium *medExt   = GetMedium("SDDKAPTON (POLYCH2)$", mgr);
    TGeoMedium *medExtB  = GetMedium("SPD-BUS CU KAPTON$", mgr);
    TGeoMedium *medExtM  = GetMedium("SPD-MCM CU KAPTON$", mgr);
    TGeoMedium *medPlate = GetMedium("SPD C (M55J)$", mgr);
    TGeoMedium *medFreon = GetMedium("Freon$", mgr);
    TGeoMedium *medGas   = GetMedium("GASEOUS FREON$", mgr);
    TGeoMedium *medFibs  = GetMedium("SDD OPTICFIB$",mgr);
    TGeoMedium *medCopper= GetMedium("COPPER$",mgr);
    TGeoMedium *medPVC   = GetMedium("PVC$",mgr);

    Double_t extThickness = fgkmm * 0.25;
    Double_t ext1Length   = fgkmm * (26.7 - 10.0);
//    Double_t ext2Length   = fgkmm * (285.0 - ext1Length + extThickness);
    Double_t ext2Length   = fgkmm * 285.0 - ext1Length + extThickness;

    const Double_t kCableThickness  =   1.5  *fgkmm;
    Double_t cableL0 =  10.0 * fgkmm;
    Double_t cableL1 = 340.0 * fgkmm - extThickness - ext1Length - ext2Length;
    Double_t cableL2 = 300.0 * fgkmm;
    //Double_t cableL3 = 570.0 * fgkmm;
    Double_t cableL3 = 57.0 * fgkmm;
    Double_t cableW1 =  11.0 * fgkmm;
    Double_t cableW2 =  30.0 * fgkmm;
    Double_t cableW3 =  50.0 * fgkmm;

    const Double_t kMCMLength       =   cableL0 + cableL1 + cableL2 + cableL3;
    const Double_t kMCMWidth        =   cableW1;
    const Double_t kMCMThickness    =   1.2  *fgkmm;

    const Double_t kPlateLength     = 200.0  *fgkmm;
    const Double_t kPlateWidth      =  50.0  *fgkmm;
    const Double_t kPlateThickness  =   5.0  *fgkmm;

    const Double_t kConeTubeRmin    =   2.0  *fgkmm;
    const Double_t kConeTubeRmax    =   3.0  *fgkmm;

    const Double_t kHorizTubeLen    = 150.0  *fgkmm;
    const Double_t kYtoHalfStave    =   9.5  *fgkmm;

    const Double_t kWaterCoolRMax   =   2.6  *fgkmm;
    const Double_t kWaterCoolThick  =   0.04 *fgkmm;
    const Double_t kWaterCoolLen    = 250.0  *fgkmm;
    const Double_t kWCPlateThick    =   0.5  *fgkmm;
    const Double_t kWCPlateWide     =  33.0  *fgkmm;
    const Double_t kWCPlateLen      = 230.0  *fgkmm;
    const Double_t kWCFittingRext1  =   2.4  *fgkmm;
    const Double_t kWCFittingRext2  =   3.7  *fgkmm;
    const Double_t kWCFittingRint1  =   1.9  *fgkmm;
    const Double_t kWCFittingRint2  = kWaterCoolRMax;
    const Double_t kWCFittingLen1   =   7.0  *fgkmm;
    const Double_t kWCFittingLen2   =   8.0  *fgkmm;
    
    const Double_t kCollWidth       =  40.0  *fgkmm;
    const Double_t kCollLength      =  60.0  *fgkmm;
    const Double_t kCollThickness   =  10.0  *fgkmm;
    const Double_t kCollTubeThick   =   1.0  *fgkmm;
    const Double_t kCollTubeRadius  =   7.0  *fgkmm;
    const Double_t kCollTubeLength  = 205.0  *fgkmm;

    const Double_t kOptFibDiamet    =   4.5  *fgkmm;

    Double_t x[12], y[12];
    Double_t xloc, yloc, zloc;

    Int_t kPurple = 6; // Purple (Root does not define it)

    TGeoVolumeAssembly* container[5];
    if (sideC)
    container[0] = new TGeoVolumeAssembly("ITSSPDConeModuleC");
    else
    container[0] = new TGeoVolumeAssembly("ITSSPDConeModuleA");
    container[1] = new TGeoVolumeAssembly("ITSSPDCoolingModuleSideA");
    container[2] = new TGeoVolumeAssembly("ITSSPDCoolingModuleSideC");
    container[3] = new TGeoVolumeAssembly("ITSSPDPatchPanelModule");
    container[4] = new TGeoVolumeAssembly("ITSSPDWaterCooling");

    // The extender on the cone as a Xtru
    x[0] = -cableL0;
    y[0] = 0.0 + 0.5 * cableW1;

    x[1] = x[0] + cableL0 + cableL1 - 0.5*(cableW2 - cableW1);
    y[1] = y[0];

    x[2] = x[0] + cableL0 + cableL1;
    y[2] = y[1] + 0.5*(cableW2 - cableW1);

    x[3] = x[2] + cableL2;
    y[3] = y[2];

    x[4] = x[3] + 0.5*(cableW3 - cableW2);
    y[4] = y[3] + 0.5*(cableW3 - cableW2);

    x[5] = x[4] + cableL3 - 0.5*(cableW3 - cableW2);
    y[5] = y[4];

    for (Int_t i = 6; i < 12; i++) {
        x[i] =  x[11 - i];
        y[i] = -y[11 - i];
    }

    TGeoXtru *shCable = new TGeoXtru(2);
    shCable->DefinePolygon(12, x, y);
    shCable->DefineSection(0, 0.0);
    shCable->DefineSection(1, kCableThickness);

    TGeoVolume *volCable = new TGeoVolume("ITSSPDExtender", shCable, medExtB);
    volCable->SetLineColor(kGreen);

    // The MCM extender on the cone as a Xtru
    TGeoBBox *shMCMExt = new TGeoBBox(0.5*kMCMLength,
				      0.5*kMCMWidth,
				      0.5*kMCMThickness);

    TGeoVolume *volMCMExt = new TGeoVolume("ITSSPDExtenderMCM",
					   shMCMExt, medExtM);
    volMCMExt->SetLineColor(kGreen+3);

    // The support plate on the cone as a composite shape
    Double_t thickness = kCableThickness + kMCMThickness;
    TGeoBBox *shOut = new TGeoBBox("ITSSPD_shape_plateout",
				   0.5*kPlateLength,
				   0.5*kPlateWidth,
				   0.5*kPlateThickness);
    TGeoBBox *shIn  = new TGeoBBox("ITSSPD_shape_platein" ,
				   0.5*kPlateLength,
				   0.5*cableW2,
				   0.5*thickness);
    Char_t string[255];
    snprintf(string, 255, "%s-%s", shOut->GetName(), shIn->GetName());
    TGeoCompositeShape *shPlate = new TGeoCompositeShape("ITSSPDPlate_shape",
				 string);

    TGeoVolume *volPlate = new TGeoVolume("ITSSPDPlate",
					  shPlate, medPlate);
    volPlate->SetLineColor(kRed);
    
    // The air cooling tubes
    TGeoBBox   *shCollBox   = new TGeoBBox("ITSSPD_shape_collector_box", 0.5*kCollLength, 0.5*kCollWidth, 0.5*kCollThickness);
    TGeoTube   *shCollTube  = new TGeoTube("ITSSPD_shape_collector_tube",kCollTubeRadius - kCollTubeThick, kCollTubeRadius, 0.5*kCollTubeLength);
    TGeoVolume *volCollBox  = new TGeoVolume("ITSSPDCollectorBox", shCollBox, medPVC);
    TGeoVolume *volCollTube = new TGeoVolume("ITSSPDCollectorTube", shCollTube, medPVC);
    volCollBox->SetLineColor(kAzure);
    volCollTube->SetLineColor(kAzure);

    // The cooling tube on the cone as a Ctub
    Double_t tubeLength = shCable->GetX(5) - shCable->GetX(0) + kYtoHalfStave -0.85;
    TGeoCtub *shTube = new TGeoCtub(0, kConeTubeRmax, 0.5*tubeLength, 0, 360,
				    0, SinD(angrot/2), -CosD(angrot/2),
				    0,              0,              1);

    TGeoVolume *volTubeA = new TGeoVolume("ITSSPDCoolingTubeOnConeA",
					  shTube, medInox);
    volTubeA->SetLineColor(kGray);

    TGeoVolume *volTubeC = new TGeoVolume("ITSSPDCoolingTubeOnConeC",
					  shTube, medInox);
    volTubeC->SetLineColor(kGray);

    // The freon in the cooling tubes on the cone as a Ctub
    TGeoCtub *shFreon = new TGeoCtub(0, kConeTubeRmin, 0.5*tubeLength, 0, 360,
				     0, SinD(angrot/2), -CosD(angrot/2),
				     0,              0,              1);

    TGeoVolume *volFreon = new TGeoVolume("ITSSPDCoolingFreonOnCone",
					  shFreon, medFreon);
    volFreon->SetLineColor(kPurple);

    TGeoVolume *volGasFr = new TGeoVolume("ITSSPDCoolingFreonGasOnCone",
					  shFreon, medGas);
    volGasFr->SetLineColor(kPurple);

    // The cooling tube inside the cylinder as a Ctub
    TGeoCtub *shCylTub = new TGeoCtub(0, kConeTubeRmax,
				      0.5*kHorizTubeLen, 0, 360,
				      0,            0,           -1,
				      0, SinD(angrot/2), CosD(angrot/2));

    TGeoVolume *volCylTubA = new TGeoVolume("ITSSPDCoolingTubeOnCylA",
					    shCylTub, medInox);
    volCylTubA->SetLineColor(kGray);

    TGeoVolume *volCylTubC = new TGeoVolume("ITSSPDCoolingTubeOnCylC",
					    shCylTub, medInox);
    volCylTubC->SetLineColor(kGray);

    // The freon in the cooling tubes in the cylinder as a Ctub
    TGeoCtub *shCylFr = new TGeoCtub(0, kConeTubeRmin,
				     0.5*kHorizTubeLen, 0, 360,
				     0,            0,           -1,
				     0, SinD(angrot/2), CosD(angrot/2));

    TGeoVolume *volCylFr = new TGeoVolume("ITSSPDCoolingFreonOnCyl",
					  shCylFr, medFreon);
    volCylFr->SetLineColor(kPurple);

    TGeoVolume *volCylGasFr = new TGeoVolume("ITSSPDCoolingFreonGasOnCyl",
					     shCylFr, medGas);
    volCylGasFr->SetLineColor(kPurple);

    // The optical fibers bundle on the cone as a Tube
    Double_t optLength = shCable->GetX(5) - shCable->GetX(0) + kYtoHalfStave -0.85;
    TGeoTube *shOptFibs = new TGeoTube(0., 0.5*kOptFibDiamet, 0.5*optLength);

    TGeoVolume *volOptFibs = new TGeoVolume("ITSSPDOpticalFibersOnCone",
					    shOptFibs, medFibs);
    volOptFibs->SetLineColor(kOrange);

    // The optical patch panels
    TArrayD psizes;
    TGeoVolume *volPatch = CreatePatchPanel(psizes, mgr);

    // The water cooling tube as a Tube
    TGeoTube *shWatCool = new TGeoTube(kWaterCoolRMax-kWaterCoolThick,
				       kWaterCoolRMax, kWaterCoolLen/2);

    TGeoVolume *volWatCool = new TGeoVolume("ITSSPDWaterCoolingOnCone",
					    shWatCool, medInox);
    volWatCool->SetLineColor(kGray);

    // The support plate for the water tubes: a Tubs and a BBox
    TGeoTubeSeg *shWCPltT = new TGeoTubeSeg(kWaterCoolRMax,
					    kWaterCoolRMax+kWCPlateThick,
					    kWCPlateLen/2, 180., 360.);

    Double_t plateBoxWide = (kWCPlateWide - 2*kWaterCoolRMax)/2;
    TGeoBBox *shWCPltB = new TGeoBBox(plateBoxWide/2,
				      kWCPlateThick/2,
				      kWCPlateLen/2);

    TGeoVolume *volWCPltT = new TGeoVolume("ITSSPDWaterCoolingTubsPlate",
					  shWCPltT, medPlate);
    volWCPltT->SetLineColor(kRed);

    TGeoVolume *volWCPltB = new TGeoVolume("ITSSPDWaterCoolingBoxPlate",
					  shWCPltB, medPlate);
    volWCPltB->SetLineColor(kRed);

    // The fitting for the water cooling tube: a Pcon
    TGeoPcon *shFitt = new TGeoPcon(0., 360., 4);
    shFitt->Z(0)    = -kWCFittingLen1;
    shFitt->Rmin(0) =  kWCFittingRint1;
    shFitt->Rmax(0) =  kWCFittingRext1;

    shFitt->Z(1)    =  0;
    shFitt->Rmin(1) =  kWCFittingRint1;
    shFitt->Rmax(1) =  kWCFittingRext1;

    shFitt->Z(2)    =  0;
    shFitt->Rmin(2) =  kWCFittingRint2;
    shFitt->Rmax(2) =  kWCFittingRext2;

    shFitt->Z(3)    =  kWCFittingLen2;
    shFitt->Rmin(3) =  kWCFittingRint2;
    shFitt->Rmax(3) =  kWCFittingRext2;

    TGeoVolume *volFitt = new TGeoVolume("ITSSPDWaterCoolingFitting",
					 shFitt, medCopper);
    volFitt->SetLineColor(kOrange);

    // Now place everything in the containers
    volTubeA->AddNode(volGasFr, 1, 0);
    volTubeC->AddNode(volFreon, 1, 0);

    volCylTubA->AddNode(volCylGasFr, 1, 0);
    volCylTubC->AddNode(volCylFr   , 1, 0);

    container[0]->AddNode(volCable, 1, 0);

    xloc = shMCMExt->GetDX() - cableL0;
    zloc = shMCMExt->GetDZ();
    container[0]->AddNode(volMCMExt, 1,
			  new TGeoTranslation( xloc, 0.,-zloc));

    xloc = shMCMExt->GetDX();
    zloc = shCable->GetZ(1)/2 - shMCMExt->GetDZ();
    container[0]->AddNode(volPlate, 1,
			  new TGeoTranslation( xloc, 0., zloc));

    TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
    rot2->SetName("rotPatch");
    rot2->RotateX(90.0);
    rot2->RotateY(163.0);
    //rot2->RotateZ(132.5);
    
    // add collectors only on side C
    if (sideC)
    {
      TGeoTranslation *trCollBox   = new TGeoTranslation(xloc - 0.5*kPlateLength + 0.5*kCollLength, 0.0, +0.5*(kPlateThickness+1.1*kCollThickness));
      TGeoRotation    *rotCollTube = new TGeoRotation(*gGeoIdentity);
      rotCollTube->RotateY(90.0);
      TGeoCombiTrans  *trCollTube  = new TGeoCombiTrans(xloc + 0.5*kCollTubeLength - (0.5*kPlateLength - kCollLength), 0.0, +0.5*(kPlateThickness+2.0*kCollTubeRadius+kCollTubeThick), rotCollTube);
      container[0]->AddNode(volCollBox, 1, trCollBox);
      container[0]->AddNode(volCollTube, 1, trCollTube);
    }
        
    Double_t dxPatch = 2.75;
    Double_t dzPatch = 2.8;
    TGeoCombiTrans *tr2 = new TGeoCombiTrans(1.7*ext2Length - dxPatch, 0.0, dzPatch, rot2);
    container[3]->AddNode(volPatch, 0, tr2);

    xloc = shTube->GetRmax();
    yloc = shTube->GetRmax();
    zloc = shTube->GetDz() - shTube->GetRmax() - kYtoHalfStave;
    container[1]->AddNode(volTubeA, 1,
			  new TGeoTranslation(-xloc, -yloc, zloc));
    container[2]->AddNode(volTubeC, 1,
			  new TGeoTranslation(-xloc, -yloc, zloc));

    xloc = shTube->GetRmax();
    yloc = (shCylTub->GetDz())*SinD(angrot) - shTube->GetRmax();
    zloc = (shCylTub->GetDz())*CosD(angrot) + shTube->GetRmax() +kYtoHalfStave;
    container[1]->AddNode(volCylTubA, 1,
			  new TGeoCombiTrans(-xloc, yloc,-zloc,
				     new TGeoRotation("",0.,angrot,0.)));
    container[2]->AddNode(volCylTubC, 1,
			  new TGeoCombiTrans(-xloc, yloc,-zloc,
				     new TGeoRotation("",0.,angrot,0.)));

    xloc = shOptFibs->GetRmax() + 2*shTube->GetRmax();
    yloc = 1.6*shOptFibs->GetRmax();
    zloc = shOptFibs->GetDZ() - shTube->GetRmax() - kYtoHalfStave;
    container[1]->AddNode(volOptFibs, 1,
			  new TGeoTranslation(-xloc, -yloc, zloc));
    container[2]->AddNode(volOptFibs, 1,
			  new TGeoTranslation(-xloc, -yloc, zloc));

    yloc = shWatCool->GetRmax();
    zloc = (2*shTube->GetDz() - shTube->GetRmax() - kYtoHalfStave)/2;
    container[4]->AddNode(volWatCool, 1,
			  new TGeoTranslation(0, -yloc, zloc));

    container[4]->AddNode(volWCPltT, 1,
			  new TGeoTranslation(0, -yloc, zloc));

    yloc -= shWCPltB->GetDY();
    xloc = shWatCool->GetRmax() + shWCPltB->GetDX();
    container[4]->AddNode(volWCPltB, 1,
			  new TGeoTranslation( xloc, -yloc, zloc));
    container[4]->AddNode(volWCPltB, 2,
			  new TGeoTranslation(-xloc, -yloc, zloc));

    yloc = shWatCool->GetRmax();
    zloc -= shWatCool->GetDz();
    container[4]->AddNode(volFitt, 1,
			  new TGeoTranslation(0, -yloc, zloc));

    // Finally create the list of assemblies and return it to the caller
    TList* conemodulelist = new TList();
    conemodulelist->Add(container[0]);
    conemodulelist->Add(container[1]);
    conemodulelist->Add(container[2]);
    conemodulelist->Add(container[3]);
    conemodulelist->Add(container[4]);

    return conemodulelist;
}

//______________________________________________________________________
void AliITSv11GeometrySPD::CreateCones(TGeoVolume *moth) const
{
    //
    // Places all services modules in the mother reference system
    //
    // Created:      ?? ??? 2008  Alberto Pulvirenti
    // Updated:      03 May 2010  Mario Sitta
    // Updated:      04 Jul 2010  Mario Sitta  Water cooling
    //

    const Int_t kNumberOfModules    =  10;

    const Double_t kInnerRadius     =  80.775*fgkmm;
    const Double_t kZTrans          = 451.800*fgkmm;
    const Double_t kAlphaRot        =  46.500*fgkDegree;
    const Double_t kAlphaSpaceCool  =   9.200*fgkDegree;

    TList*  modulelistA = CreateConeModule(kFALSE, 90-kAlphaRot);
    TList*  modulelistC = CreateConeModule(kTRUE , 90-kAlphaRot);
    TList* &modulelist  = modulelistC;
    TGeoVolumeAssembly* module, *moduleA, *moduleC;

    Double_t xloc, yloc, zloc;

    //Double_t angle[10] = {18., 54., 90., 126., 162., -18., -54., -90., -126., -162.};
    // anglem for cone modules (cables and cooling tubes)
    // anglep for pathc panels
    Double_t anglem[10] = {18., 54., 90., 126., 162., 198., 234., 270., 306., 342.};
    Double_t anglep[10] = {18., 62., 90., 115., 162., 198., 242., 270., 295., 342.};
//    Double_t angle1m[10] = {23., 53., 90., 127., 157., 203.0, 233.0, 270.0, 307.0, 337.0};
//    Double_t angle2m[10] = {18., 53., 90., 126., 162., 198.0, 233.0, 270.0, 309.0, 342.0};
//    Double_t angle1c[10] = {23., 53., 90., 124., 157., 203.0, 233.0, 270.0, 304.0, 337.0};
//    Double_t angle2c[10] = {18., 44., 90., 126., 162., 198.0, 223.0, 270.0, 309.0, 342.0};

    // First add the cables
    moduleA = (TGeoVolumeAssembly*)modulelistA->At(0);
    moduleC = (TGeoVolumeAssembly*)modulelistC->At(0);
    for (Int_t i = 0; i < kNumberOfModules; i++) {
        TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
	rot1->RotateY(-kAlphaRot);
	rot1->RotateZ(anglem[i]);
        xloc = kInnerRadius*CosD(anglem[i]);
        yloc = kInnerRadius*SinD(anglem[i]);
	zloc = kZTrans;
        moth->AddNode(moduleA, 2*i+2,
		      new TGeoCombiTrans( xloc, yloc, zloc, rot1));

        TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
	rot2->RotateY(180.-kAlphaRot);
	rot2->RotateZ(anglem[i]);
        xloc = kInnerRadius*CosD(anglem[i]);
        yloc = kInnerRadius*SinD(anglem[i]);
	zloc = kZTrans;
        moth->AddNode(moduleC, 2*i+1,
		      new TGeoCombiTrans(-xloc,-yloc,-zloc, rot2));
    }

    // Then the cooling tubes on Side A
    module = (TGeoVolumeAssembly*)modulelist->At(1);
    Double_t anglec;
    for (Int_t i = 0; i < kNumberOfModules; i++) {
        anglec = anglem[i] + kAlphaSpaceCool;
        TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
        rot1->RotateX(-90.0+kAlphaRot-0.04); // 0.04 fixes small overlap
	rot1->RotateZ(-90.0+anglec);
        xloc = kInnerRadius*CosD(anglec);
        yloc = kInnerRadius*SinD(anglec);
	zloc = kZTrans+0.162; // 0.162 fixes small overlap
        moth->AddNode(module, 2*i+2, 
		      new TGeoCombiTrans( xloc, yloc, zloc, rot1));
    }

    // And the cooling tubes on Side C
    module = (TGeoVolumeAssembly*)modulelist->At(2);
    for (Int_t i = 0; i < kNumberOfModules; i++) {
        anglec = anglem[i] - kAlphaSpaceCool;
        TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
        rot2->RotateX(-90.0+kAlphaRot-0.04); // 0.04 fixes small overlap
	rot2->RotateY(180.);
	rot2->RotateZ(90.0+anglec);
        xloc = kInnerRadius*CosD(anglec);
        yloc = kInnerRadius*SinD(anglec);
	zloc = kZTrans+0.162; // 0.162 fixes small overlap
        moth->AddNode(module, 2*i+1,
		      new TGeoCombiTrans(-xloc,-yloc,-zloc, rot2));
    }

    // Then the water cooling tubes
    module = (TGeoVolumeAssembly*)modulelist->At(4);
    for (Int_t i = 1; i < kNumberOfModules; i++) { // i = 1,2,...,9
        if (i != 5) { // There is no tube in this position
	  anglec = (anglem[i-1]+anglem[i])/2;
	    TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
	    rot1->RotateX(-90.0+kAlphaRot);
	    rot1->RotateZ(-90.0+anglec);
	    xloc = kInnerRadius*CosD(anglec);
	    yloc = kInnerRadius*SinD(anglec);
	    zloc = kZTrans;
	    moth->AddNode(module, 2*i+2,
			  new TGeoCombiTrans( xloc, yloc, zloc, rot1));

	    TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
	    rot2->RotateX(-90.0+kAlphaRot);
	    rot2->RotateY(180.);
	    rot2->RotateZ(90.0+anglec);
	    xloc = kInnerRadius*CosD(anglec);
	    yloc = kInnerRadius*SinD(anglec);
	    zloc = kZTrans;
	    moth->AddNode(module, 2*i+1,
			  new TGeoCombiTrans(-xloc,-yloc,-zloc, rot2));
	}
    }

    // Finally the optical patch panels
    module = (TGeoVolumeAssembly*)modulelist->At(3);
    for (Int_t i = 0; i < kNumberOfModules; i++) {
        TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
	rot1->RotateY(-kAlphaRot);
	rot1->RotateZ(anglep[i]);
        xloc = kInnerRadius*CosD(anglep[i]);
        yloc = kInnerRadius*SinD(anglep[i]);
	zloc = kZTrans;
        moth->AddNode(module, 2*i+2,
		      new TGeoCombiTrans( xloc, yloc, zloc, rot1));

        TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
	rot2->RotateY(180.-kAlphaRot);
	rot2->RotateZ(anglep[i]);
        xloc = kInnerRadius*CosD(anglep[i]);
        yloc = kInnerRadius*SinD(anglep[i]);
	zloc = kZTrans;
        moth->AddNode(module, 2*i+1,
		      new TGeoCombiTrans(-xloc,-yloc,-zloc, rot2));
    }

}


//______________________________________________________________________
void AliITSv11GeometrySPD::CreateServices(TGeoVolume *moth) const
{
    //
    // New method to implement SPD services
    //
    // Created:      25 Jul 2012  Mario Sitta
    // Updated:      15 Nov 2012  Mario Sitta
    //
    // Data provided by C.Gargiulo from CAD

    // Cooling manifolds
    const Double_t kCoolManifWidth    = fgkmm * 22.0;
    const Double_t kCoolManifLength   = fgkmm * 50.0;
    const Double_t kCoolManifThick    = fgkmm *  7.0;
    const Double_t kCoolManifFitR1out = fgkmm *  4.0;
    const Double_t kCoolManifFitH1    = fgkmm *  2.5;
    const Double_t kCoolManifFitR2out = fgkmm *  4.0;
    const Double_t kCoolManifFitR2in  = fgkmm *  3.2;
    const Double_t kCoolManifFitH2    = fgkmm *  7.0;
    const Double_t kCoolManifFitZPos  = fgkmm *  2.0; // TO BE CHECKED!
    const Double_t kCoolManifCollR1   = fgkmm *  3.0;
    const Double_t kCoolManifCollH1   = fgkmm *  2.5;
    const Double_t kCoolManifCollR2   = fgkmm *  1.5;
    const Double_t kCoolManifCollH2   = fgkmm *  5.0;
    const Double_t kCoolManifCollXPos = fgkmm *  5.0;
    const Double_t kCoolManifCollDZ   = fgkmm * 13.0;
    const Double_t kCoolManifCollZ0   = fgkmm *  9.0;

    const Double_t kCoolManifRPosCAD  = fgkmm * 76.2;
    const Double_t kCoolManifZPos     = fgkcm * 33.97;// 34.0 - 0.03 toll.
    // Manifold supports
    const Double_t kManifSuppWidth    = fgkmm * 24.0; // TO BE CHECKED!
    const Double_t kManifSuppLen1     = fgkmm * 17.9;
    const Double_t kManifSuppLen2     = fgkmm * 54.2;
    const Double_t kManifSuppLen3     = fgkmm *  7.9;
    const Double_t kManifSuppThick    = fgkmm *  1.5;
    const Double_t kSuppScrewXPos     = fgkmm *  4.0;
    const Double_t kSuppScrewZPos     = fgkmm *  3.0;
    const Double_t kRThermalShield    = fgkcm *  9.9255; // MUST match with GeometrySupport
    // Sector supports
    const Double_t kSectSuppWidth     = fgkmm * 15.0;
    const Double_t kSectSuppLen1      = fgkmm * 16.9; // TO BE CHECKED!
    const Double_t kSectSuppLen2      = fgkmm * 35.1; // TO BE CHECKED!
    const Double_t kSectSuppThick     = fgkmm *  1.5;
    const Double_t kSectSuppDepth     = fgkmm * 17.78; // MUST match with GeometrySupport
    const Double_t kSectScrewZPos     = fgkmm *  5.1; // TO BE CHECKED!

    const Double_t kSectSuppZPos      = fgkcm * 26.5;
    // Sector clips
    const Double_t kSectClipLength    = fgkmm * 30.0;
    const Double_t kSectClipWidth     = fgkmm * 28.53;
    const Double_t kSectClipThick1    = fgkmm *  2.0;
    const Double_t kSectClipThick2    = fgkmm *  0.715;
    const Double_t kSectClipInStave   = fgkmm * 11.0; // Tuned
    const Double_t kSectClipAngle     =         29.0; // Degree. Tuned
    // M3 screws
    const Double_t kScrewM3Diam       = fgkmm *  3.0;
    const Double_t kScrewM3HeadThick  = fgkmm *  2.0;
    const Double_t kScrewM3HeadRmin   = fgkmm *  1.5;
    const Double_t kScrewM3HeadRmax   = fgkmm *  2.5;
    const Double_t kScrewM3OutManifH  = fgkmm *  1.5;
    // Central set pin (in sector support)
    const Double_t kSetPinDiam        = fgkmm *  6.0;
    const Double_t kSetPinHeadDiam    = fgkmm *  8.0;
    const Double_t kSetPinHeadRmin    = fgkmm *  1.5;
    const Double_t kSetPinHeadThick   = fgkmm *  1.5;
    const Double_t kSetPinOutClipH    = fgkmm *  1.0;

    // Local variables
    Double_t xprof[12], yprof[12];
    Double_t radius, theta;
    Double_t xpos, ypos, zpos;
    Double_t tmp;


    // The cooling manifold: an Assembly
    TGeoVolumeAssembly *coolmanifA = new TGeoVolumeAssembly("ITSSPDCoolManifSideA");
    TGeoVolumeAssembly *coolmanifC = new TGeoVolumeAssembly("ITSSPDCoolManifSideC");

    // The various parts of the manifold
    TGeoBBox *manifblksh = new TGeoBBox(kCoolManifWidth/2,
					kCoolManifThick/2,
					kCoolManifLength/2);

    TGeoBBox *manifinscubesh = new TGeoBBox(kCoolManifFitR2out,
					    kCoolManifFitR2out,
					    kCoolManifFitR2out);

    TGeoTube *manifinscyl1sh = new TGeoTube(0, // TO BE CHECKED!
					    kCoolManifFitR1out,
					    kCoolManifFitH1/2);

    TGeoTube *manifinscyl2sh = new TGeoTube(kCoolManifFitR2in,
					    kCoolManifFitR2out,
					    kCoolManifFitH2/2);

    TGeoTube *manifcollcyl1sh = new TGeoTube(0,
					     kCoolManifCollR1,
					     kCoolManifCollH1/2);

    TGeoTube *manifcollcyl2sh = new TGeoTube(0,
					     kCoolManifCollR2,
					     kCoolManifCollH2/2);

    // The cooling manifold supports
    const Double_t kCoolManifRPos = kCoolManifRPosCAD  +
			      (manifinscubesh->GetDY() +
			     2*manifinscyl1sh->GetDz() +
			       manifblksh->GetDY()     );

    const Double_t kManifSuppDepth = kRThermalShield -
				    (kCoolManifRPos + manifblksh->GetDY());

    TGeoXtru *suppmanifsh = new TGeoXtru(2);

    xprof[ 0] = kManifSuppLen2/2 + kManifSuppThick;
    yprof[ 0] = 0;
    xprof[ 1] = xprof[0];
    yprof[ 1] = kManifSuppDepth;
    xprof[ 2] = kManifSuppLen2/2 + kManifSuppLen3;
    yprof[ 2] = yprof[1];
    xprof[ 3] = xprof[2];
    yprof[ 3] = yprof[2] + kManifSuppThick;
    xprof[ 4] = kManifSuppLen2/2;
    yprof[ 4] = yprof[3];
    xprof[ 5] = xprof[4];
    yprof[ 5] = kManifSuppThick;
    xprof[ 6] = -xprof[5];
    yprof[ 6] =  yprof[5];
    xprof[ 7] = -xprof[4];
    yprof[ 7] =  yprof[4];
    xprof[ 8] = -(kManifSuppLen2/2 + kManifSuppLen1);
    yprof[ 8] =  yprof[3];
    xprof[ 9] =  xprof[8];
    yprof[ 9] =  yprof[2];
    xprof[10] = -xprof[1];
    yprof[10] =  yprof[1];
    xprof[11] = -xprof[0];
    yprof[11] =  yprof[0];

    suppmanifsh->DefinePolygon(12,xprof,yprof);
    suppmanifsh->DefineSection(0,-kManifSuppWidth/2);
    suppmanifsh->DefineSection(1, kManifSuppWidth/2);

    // The screw head and body
    TGeoTube *suppscrewbodysh = new TGeoTube(0, kScrewM3Diam/2,
					     kManifSuppThick/2);

    TGeoPcon *suppscrewheadsh = new TGeoPcon(0, 360, 4);
    suppscrewheadsh->DefineSection(0,-kScrewM3HeadThick/2,0, kScrewM3HeadRmax);
    suppscrewheadsh->DefineSection(1, 0,                  0, kScrewM3HeadRmax);
    suppscrewheadsh->DefineSection(2, 0,   kScrewM3HeadRmin, kScrewM3HeadRmax);
    suppscrewheadsh->DefineSection(3, kScrewM3HeadThick/2,
					 kScrewM3HeadRmin, kScrewM3HeadRmax);

    TGeoTube *clipscrewbodysh = new TGeoTube(0, kScrewM3Diam/2,
					     kSectClipThick1/2);

    // The screw segment below the manifold and the sector clip
    TGeoTube *screwoutmanifsh = new TGeoTube(0, kScrewM3Diam/2,
					     kScrewM3OutManifH/2);

    // The sector supports
    TGeoXtru *suppsectsh = new TGeoXtru(2);

    xprof[ 0] = kSectSuppLen2/2 + kSectSuppThick;
    yprof[ 0] = 0;
    xprof[ 1] = xprof[0];
    yprof[ 1] = kSectSuppDepth;
    xprof[ 2] = kSectSuppLen2/2 + kSectSuppLen1;
    yprof[ 2] = yprof[1];
    xprof[ 3] = xprof[2];
    yprof[ 3] = yprof[2] + kSectSuppThick;
    xprof[ 4] = kSectSuppLen2/2;
    yprof[ 4] = yprof[3];
    xprof[ 5] = xprof[4];
    yprof[ 5] = kSectSuppThick;
    xprof[ 6] = -xprof[5];
    yprof[ 6] =  yprof[5];
    xprof[ 7] = -xprof[4];
    yprof[ 7] =  yprof[4];
    xprof[ 8] = -xprof[3];
    yprof[ 8] =  yprof[3];
    xprof[ 9] = -xprof[2];
    yprof[ 9] =  yprof[2];
    xprof[10] = -xprof[1];
    yprof[10] =  yprof[1];
    xprof[11] = -xprof[0];
    yprof[11] =  yprof[0];

    suppsectsh->DefinePolygon(12,xprof,yprof);
    suppsectsh->DefineSection(0,-kSectSuppWidth/2);
    suppsectsh->DefineSection(1, kSectSuppWidth/2);

    // The sector clips
    TGeoXtru *sectclipsh = new TGeoXtru(2);

    xprof[ 0] =  kSectClipWidth/2;
    yprof[ 0] =  0;
    xprof[ 1] = -kSectClipWidth/2;
    yprof[ 1] =  yprof[0];
    xprof[ 2] =  xprof[1];
    yprof[ 2] = -kSectClipThick1;
    xprof[ 3] =  kSectClipWidth/2 - kSectClipThick2;
    yprof[ 3] =  yprof[2];
    xprof[ 4] =  xprof[3] + kSectClipInStave*SinD(kSectClipAngle);
    yprof[ 4] = -kSectClipInStave*CosD(kSectClipAngle);
    xprof[ 5] =  xprof[4] + kSectClipThick2*CosD(kSectClipAngle);
    yprof[ 5] =  yprof[4] + kSectClipThick2*SinD(kSectClipAngle);

    sectclipsh->DefinePolygon(6,xprof,yprof);
    sectclipsh->DefineSection(0,-kSectClipLength/2);
    sectclipsh->DefineSection(1, kSectClipLength/2);

    // The central set pin head and body
    TGeoTube *setpinbodysh = new TGeoTube(0, kSetPinDiam/2,
					  kSectSuppThick/2);

    TGeoTube *setpinheadsh = new TGeoTube(kSetPinHeadRmin, kSetPinHeadDiam/2,
					  kSetPinHeadThick/2);

    TGeoTube *pinclipbodysh = new TGeoTube(0, kSetPinDiam/2,
					   kSectClipThick1/2);

    // The set pin segment below the sector clip
    TGeoTube *setpinoutclipsh = new TGeoTube(0, kSetPinDiam/2,
					     kSetPinOutClipH/2);


    // We have the shapes: now create the real volumes
    TGeoMedium *medInox  = GetMedium("INOX$");
    TGeoMedium *medCu    = GetMedium("COPPER$");
    TGeoMedium *medSPDcf = GetMedium("SPD shield$");

    TGeoVolume *manifblk = new TGeoVolume("ITSSPDBlkManif",
					  manifblksh,medInox);
    manifblk->SetLineColor(kGreen+2);

    TGeoVolume *manifinscube = new TGeoVolume("ITSSPDInsCubeManif",
					      manifinscubesh,medCu);
    manifinscube->SetLineColor(kYellow);

    TGeoVolume *manifinscyl1 = new TGeoVolume("ITSSPDInsCyl1Manif",
					      manifinscyl1sh,medCu);
    manifinscyl1->SetLineColor(kYellow);

    TGeoVolume *manifinscyl2 = new TGeoVolume("ITSSPDInsCyl2Manif",
					      manifinscyl2sh,medCu);
    manifinscyl2->SetLineColor(kYellow);

    TGeoVolume *manifcollcyl1 = new TGeoVolume("ITSSPDCollCyl1Manif",
					       manifcollcyl1sh,medCu);
    manifcollcyl1->SetLineColor(kYellow);

    TGeoVolume *manifcollcyl2 = new TGeoVolume("ITSSPDCollCyl2Manif",
					       manifcollcyl2sh,medCu);
    manifcollcyl2->SetLineColor(kYellow);

    TGeoVolume *suppmanif = new TGeoVolume("ITSSPDCoolManifSupp",
					       suppmanifsh,medSPDcf);
    suppmanif->SetLineColor(7);

    TGeoVolume *suppscrewbody = new TGeoVolume("ITSSPDSuppScrewBody",
					       suppscrewbodysh,medInox);
    suppscrewbody->SetLineColor(kGray);

    xpos = kCoolManifLength/2 - kSuppScrewZPos;
    ypos = suppscrewbodysh->GetDz();
    zpos = kCoolManifWidth/2  - kSuppScrewXPos;
    suppmanif->AddNode(suppscrewbody, 1, new TGeoCombiTrans( xpos, ypos, zpos,
					 new TGeoRotation("",0,90,0)));
    suppmanif->AddNode(suppscrewbody, 2, new TGeoCombiTrans( xpos, ypos,-zpos,
					 new TGeoRotation("",0,90,0)));
    suppmanif->AddNode(suppscrewbody, 3, new TGeoCombiTrans(-xpos, ypos, zpos,
					 new TGeoRotation("",0,90,0)));
    suppmanif->AddNode(suppscrewbody, 4, new TGeoCombiTrans(-xpos, ypos,-zpos,
					 new TGeoRotation("",0,90,0)));

    TGeoVolume *suppscrewhead = new TGeoVolume("ITSSPDSuppScrewHead",
					       suppscrewheadsh,medInox);
    suppscrewhead->SetLineColor(kGray);

    TGeoVolume *screwoutmanif = new TGeoVolume("ITSSPDSuppScrewOutManif",
					       screwoutmanifsh,medInox);
    screwoutmanif->SetLineColor(kGray);

    TGeoVolume *suppsect = new TGeoVolume("ITSSPDCoolSectorSupp",
					  suppsectsh,medSPDcf);
    suppsect->SetLineColor(7);

    xpos = kSectSuppLen2/2 - kSectScrewZPos;
    ypos = suppscrewbodysh->GetDz();
    suppsect->AddNode(suppscrewbody, 1, new TGeoCombiTrans( xpos, ypos, 0,
					new TGeoRotation("",0,90,0)));
    suppsect->AddNode(suppscrewbody, 2, new TGeoCombiTrans(-xpos, ypos, 0,
					new TGeoRotation("",0,90,0)));

    TGeoVolume *setpinbody = new TGeoVolume("ITSSPDSetPinBody",
					    setpinbodysh,medInox);
    setpinbody->SetLineColor(kGray);

    ypos = setpinbodysh->GetDz();
    suppsect->AddNode(setpinbody, 1, new TGeoCombiTrans( 0, ypos, 0,
					new TGeoRotation("",0,90,0)));

    TGeoVolume *setpinhead = new TGeoVolume("ITSSPDSetPinHead",
					    setpinheadsh,medInox);
    setpinhead->SetLineColor(kGray);

    TGeoVolume *sectclip = new TGeoVolume("ITSSPDCoolSectorClip",
					  sectclipsh,medSPDcf);
    sectclip->SetLineColor(7);

    TGeoVolume *clipscrewbody = new TGeoVolume("ITSSPDClipScrewBody",
					       clipscrewbodysh,medInox);
    clipscrewbody->SetLineColor(kGray);

    ypos = -clipscrewbodysh->GetDz();
    zpos = kSectSuppLen2/2 - kSectScrewZPos;
    sectclip->AddNode(clipscrewbody, 1, new TGeoCombiTrans( 0, ypos, zpos,
					new TGeoRotation("",0,90,0)));
    sectclip->AddNode(clipscrewbody, 2, new TGeoCombiTrans( 0, ypos,-zpos,
					new TGeoRotation("",0,90,0)));

    TGeoVolume *pinclipbody = new TGeoVolume("ITSSPDClipPinBody",
					     pinclipbodysh,medInox);
    pinclipbody->SetLineColor(kGray);

    ypos = -pinclipbodysh->GetDz();
    sectclip->AddNode(pinclipbody, 1, new TGeoCombiTrans( 0, ypos, 0,
					new TGeoRotation("",0,90,0)));

    TGeoVolume *setpinoutclip = new TGeoVolume("ITSSPDSetPinOutClip",
					       setpinoutclipsh,medInox);
    setpinoutclip->SetLineColor(kGray);


    // Add all volumes in the assemblies
    coolmanifA->AddNode(manifblk,1,0);
    coolmanifC->AddNode(manifblk,1,0);

    ypos = manifblksh->GetDY() + manifinscyl1sh->GetDz();
    zpos = manifblksh->GetDZ() - manifinscyl1sh->GetRmax() - kCoolManifFitZPos;
    coolmanifA->AddNode(manifinscyl1, 1, new TGeoCombiTrans(0, -ypos, zpos,
					 new TGeoRotation("",0,90,0)));
    coolmanifC->AddNode(manifinscyl1, 1, new TGeoCombiTrans(0, -ypos, zpos,
					 new TGeoRotation("",0,90,0)));

    ypos += (manifinscyl1sh->GetDz() + manifinscubesh->GetDY());
    coolmanifA->AddNode(manifinscube, 1, new TGeoTranslation(0, -ypos, zpos));
    coolmanifC->AddNode(manifinscube, 1, new TGeoTranslation(0, -ypos, zpos));

    zpos += (manifinscubesh->GetDZ() + manifinscyl2sh->GetDz());
    coolmanifA->AddNode(manifinscyl2, 1, new TGeoTranslation(0, -ypos, zpos));
    coolmanifC->AddNode(manifinscyl2, 1, new TGeoTranslation(0, -ypos, zpos));

    ypos = manifblksh->GetDY();
    coolmanifA->AddNode(suppmanif, 1, new TGeoCombiTrans(0, ypos, 0,
					 new TGeoRotation("",-90,90,90)));
    coolmanifC->AddNode(suppmanif, 1, new TGeoCombiTrans(0, ypos, 0,
					 new TGeoRotation("",-90,90,90)));

    ypos += (kManifSuppThick + kScrewM3HeadThick/2);
    xpos = kCoolManifWidth/2   - kSuppScrewXPos;
    zpos = kCoolManifLength/2  - kSuppScrewZPos;
    coolmanifA->AddNode(suppscrewhead, 1, new TGeoCombiTrans( xpos, ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 1, new TGeoCombiTrans( xpos, ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(suppscrewhead, 2, new TGeoCombiTrans( xpos, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 2, new TGeoCombiTrans( xpos, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(suppscrewhead, 3, new TGeoCombiTrans(-xpos, ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 3, new TGeoCombiTrans(-xpos, ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(suppscrewhead, 4, new TGeoCombiTrans(-xpos, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 4, new TGeoCombiTrans(-xpos, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    ypos = manifblksh->GetDY() + screwoutmanifsh->GetDz();
    coolmanifA->AddNode(screwoutmanif, 1, new TGeoCombiTrans( xpos,-ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 1, new TGeoCombiTrans( xpos,-ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(screwoutmanif, 2, new TGeoCombiTrans( xpos,-ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 2, new TGeoCombiTrans( xpos,-ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(screwoutmanif, 3, new TGeoCombiTrans(-xpos,-ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 3, new TGeoCombiTrans(-xpos,-ypos, zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifA->AddNode(screwoutmanif, 4, new TGeoCombiTrans(-xpos,-ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 4, new TGeoCombiTrans(-xpos,-ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    ypos = manifblksh->GetDY() + suppmanifsh->GetY(1) - suppsectsh->GetY(1);
    zpos = manifblksh->GetDZ() + (kCoolManifZPos - kSectSuppZPos);
    coolmanifA->AddNode(suppsect, 1, new TGeoCombiTrans(0, ypos,-zpos,
					 new TGeoRotation("",-90,90,90)));
    coolmanifC->AddNode(suppsect, 1, new TGeoCombiTrans(0, ypos,-zpos,
					 new TGeoRotation("",-90,90,90)));

    tmp = ypos; // Save it to avoid recomputing

    ypos += (kSectSuppThick + kScrewM3HeadThick/2);
    zpos += (kSectSuppLen2/2 - kSectScrewZPos);
    coolmanifA->AddNode(suppscrewhead, 5, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 5, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    zpos -= 2*(kSectSuppLen2/2 - kSectScrewZPos);
    coolmanifA->AddNode(suppscrewhead, 6, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(suppscrewhead, 6, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    ypos = tmp + kSectSuppThick + kSetPinHeadThick/2;
    zpos += (kSectSuppLen2/2 - kSectScrewZPos);
    coolmanifA->AddNode(setpinhead, 1, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(setpinhead, 1, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    ypos = tmp - 8.e-5; // Avoid microscopic overlap
    tmp = ypos;
    coolmanifA->AddNode(sectclip, 1, new TGeoTranslation( 0, ypos,-zpos));
    coolmanifC->AddNode(sectclip, 1, new TGeoCombiTrans ( 0, ypos,-zpos,
					  new TGeoRotation("",-90,180,90)));

    ypos -= (kSectClipThick1 + setpinoutclipsh->GetDz());
    coolmanifA->AddNode(setpinoutclip, 1, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(setpinoutclip, 1, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    ypos = tmp - (kSectClipThick1 + screwoutmanifsh->GetDz());
    zpos += (kSectSuppLen2/2 - kSectScrewZPos);
    coolmanifA->AddNode(screwoutmanif, 5, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 5, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    zpos -= 2*(kSectSuppLen2/2 - kSectScrewZPos);
    coolmanifA->AddNode(screwoutmanif, 6, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));
    coolmanifC->AddNode(screwoutmanif, 6, new TGeoCombiTrans( 0, ypos,-zpos,
					  new TGeoRotation("",0,-90,0)));

    xpos = manifblksh->GetDX() - kCoolManifCollXPos;
    ypos = manifblksh->GetDY() + manifcollcyl1sh->GetDz();
    zpos =-manifblksh->GetDZ() + kCoolManifCollZ0;
    for (Int_t i=0; i<3; i++) {
      coolmanifA->AddNode(manifcollcyl1, 2*i+1,
			  new TGeoCombiTrans( xpos, -ypos, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifA->AddNode(manifcollcyl1, 2*i+2,
			  new TGeoCombiTrans(-xpos, -ypos, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifC->AddNode(manifcollcyl1, 2*i+1,
			  new TGeoCombiTrans( xpos, -ypos, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifC->AddNode(manifcollcyl1, 2*i+2,
			  new TGeoCombiTrans(-xpos, -ypos, zpos,
					     new TGeoRotation("",0,90,0)));
      Double_t y = ypos + manifcollcyl1sh->GetDz() + manifcollcyl2sh->GetDz();
      coolmanifA->AddNode(manifcollcyl2, 2*i+1,
			  new TGeoCombiTrans( xpos, -y, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifA->AddNode(manifcollcyl2, 2*i+2,
			  new TGeoCombiTrans(-xpos, -y, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifC->AddNode(manifcollcyl2, 2*i+1,
			  new TGeoCombiTrans( xpos, -y, zpos,
					     new TGeoRotation("",0,90,0)));
      coolmanifC->AddNode(manifcollcyl2, 2*i+2,
			  new TGeoCombiTrans(-xpos, -y, zpos,
					     new TGeoRotation("",0,90,0)));

      zpos += kCoolManifCollDZ;
    }

    // Now add the cooling tubes to the assembly
    CreateCoolingTubes(coolmanifA, kFALSE);
    CreateCoolingTubes(coolmanifC, kTRUE);


    // Finally put everything in the mother volume
    radius = kCoolManifRPos + 1.e-5; // Avoid microscopic overlap
    zpos = kCoolManifZPos + manifblksh->GetDZ();
    for (Int_t i=0; i<10; i++) {
      theta = 36.*i;
      moth->AddNode(coolmanifA, i+1, new TGeoCombiTrans(radius*SinD(theta),
							radius*CosD(theta),
							zpos,
					  new TGeoRotation("",-theta,0,0)));
      moth->AddNode(coolmanifC, i+1, new TGeoCombiTrans(radius*SinD(theta),
							radius*CosD(theta),
						       -zpos,
					  new TGeoRotation("",90-theta,180,-90)));
    }


}


//______________________________________________________________________
void AliITSv11GeometrySPD::CreateCoolingTubes(TGeoVolume *moth, Bool_t sideC) const
{
    //
    // Private method to implement SPD cooling tubes
    // going from the manifolds to the staves
    // Since their form is quite complicate (especially on Side C
    // where capillaries are located) a separate method is used
    // If sideC is true, the cooling tubes on Side C are created
    // along with the cooling loops (aka "capillaries"), otherwise
    // the (simpler) tubes on Side A get created.
    //
    // In all variables:  L = Left (X > 0)   R = Right (X < 0)
    //
    // Created:      10 Nov 2012  Mario Sitta
    //
    // Data provided by C.Gargiulo from CAD

    // Cooling manifolds - THESE VALUES *MUST* MATCH WITH CALLING METHOD!
    const Double_t kCoolManifWidth    = fgkmm * 22.0;
    const Double_t kCoolManifLength   = fgkmm * 50.0;
    const Double_t kCoolManifThick    = fgkmm *  7.0;
    const Double_t kCoolManifCollH1   = fgkmm *  2.5;
    const Double_t kCoolManifCollH2   = fgkmm *  5.0;
    // Cooling pipes
    const Double_t kCoolPipeSideARin  = fgkmm *  1.5;
    const Double_t kCoolPipeSideARout = fgkmm *  1.8;
    const Double_t kCoolPipeSideCRin  = fgkmm *  0.5;
    const Double_t kCoolPipeSideCRout = fgkmm *  0.85;
    const Double_t kCoolPipeHeight    = fgkmm *  1.923;
    const Double_t kCoolPipeCRadiusL[3] = {11.0, 14.0, 31.34};// TO BE CHECKED!
    const Double_t kCoolPipeCRadiusR[3] = {12.0, 14.0, 35.54};// TO BE CHECKED!
    const Double_t kCoolPipeARadiusL12[2] = {14.0, 30.0};
    const Double_t kCoolPipeARadiusR12[2] = {14.0, 30.0};
    const Double_t kCoolPipeARadiusL34[2] = {22.0, 30.0};
    const Double_t kCoolPipeARadiusR34[2] = {22.0, 30.0};
    const Double_t kCoolPipeARadiusL[3]= {14.0, 14.0, 31.34}; // TO BE CHECKED!
    const Double_t kCoolPipeARadiusR[3]= {14.0, 14.0, 35.54}; // TO BE CHECKED!
    const Double_t kCoolPipeZSPD      = fgkcm *  8.47;
    // Cooling pipes position - THESE VALUES *MUST* MATCH WITH CALLING METHOD!
    const Double_t kCoolManifCollXPos = fgkmm *  5.0;
    const Double_t kCoolManifCollDZ   = fgkmm * 13.0;
    const Double_t kCoolManifCollZ0   = fgkmm *  9.0;

    Int_t kPurple = 6; // Purple (Root does not define it)

    // Local variables
    Double_t xpos, ypos, zpos;
    Char_t pipename[11];

    //
    TGeoMedium *medPhynox  = GetMedium("PHYNOX$");
    TGeoMedium *medFreon   = GetMedium("Freon$");
    TGeoMedium *medGasFr   = GetMedium("GASEOUS FREON$");

    // The cooling tubes are created as CableRound volumes
    // because it's easier to compose them piece by piece
    AliITSv11GeomCableRound *coolpipe[6];

    if (sideC)
      for (Int_t i = 0; i<6; i++) {
	snprintf(pipename,11,"coolPipeC%d",i+1);
	coolpipe[i] = new AliITSv11GeomCableRound(pipename,kCoolPipeSideCRout);
	coolpipe[i]->SetNLayers(2);
	coolpipe[i]->SetLayer(0, kCoolPipeSideCRin, medFreon, kPurple);
	coolpipe[i]->SetLayer(1,(kCoolPipeSideCRout-kCoolPipeSideCRin),
			      medPhynox, kYellow);
      }
    else
      for (Int_t i = 0; i<6; i++) {
	snprintf(pipename,11,"coolPipeA%d",i+1);
	coolpipe[i] = new AliITSv11GeomCableRound(pipename,kCoolPipeSideARout);
	coolpipe[i]->SetNLayers(2);
	coolpipe[i]->SetLayer(0, kCoolPipeSideARin, medGasFr, kPurple);
	coolpipe[i]->SetLayer(1,(kCoolPipeSideARout-kCoolPipeSideARin),
			      medPhynox, kYellow);
      }

     // Now place them in the mother assembly
     xpos = kCoolManifWidth/2  - kCoolManifCollXPos;
     ypos = kCoolManifThick/2  + kCoolManifCollH1 + kCoolManifCollH2;
     zpos =-kCoolManifLength/2 + kCoolManifCollZ0;

     if (sideC) { // On Side C tubes are simpler and can be created in a loop

       for (Int_t i=0; i<3; i++) {

	 Double_t coordL[3] = { xpos,-ypos,zpos};
	 Double_t coordR[3] = {-xpos,-ypos,zpos};
	 Double_t vect[3] = {0, 1, 0};
	 coolpipe[2*i]->AddCheckPoint(moth, 0, coordL, vect);
	 coolpipe[2*i+1]->AddCheckPoint(moth, 0, coordR, vect);
	 coordL[1] -= kCoolPipeHeight;
	 coordR[1] = coordL[1];
	 coolpipe[2*i]->AddCheckPoint(moth, 1, coordL, vect);
	 coolpipe[2*i+1]->AddCheckPoint(moth, 1, coordR, vect);
	 coordL[1] -= kCoolPipeCRadiusL[i]*fgkmm;
	 coordL[2] -= kCoolPipeCRadiusL[i]*fgkmm;
	 coordR[1] -= kCoolPipeCRadiusR[i]*fgkmm;
	 coordR[2] -= kCoolPipeCRadiusR[i]*fgkmm;
	 vect[1] = 0;
	 vect[2] = -1;
	 coolpipe[2*i]->AddCheckPoint(moth, 2, coordL, vect);
	 coolpipe[2*i+1]->AddCheckPoint(moth, 2, coordR, vect);
	 coordL[2] = -kCoolPipeZSPD;
	 coordR[2] = -kCoolPipeZSPD;
	 coolpipe[2*i]->AddCheckPoint(moth, 3, coordL, vect);
	 coolpipe[2*i+1]->AddCheckPoint(moth, 3, coordR, vect);

	 zpos += kCoolManifCollDZ;
       }

       for (Int_t i=0; i<6; i++) {
	 coolpipe[i]->SetInitialNode(moth);
	 
	 coolpipe[i]->CreateAndInsertTubeSegment(1);
	 coolpipe[i]->CreateAndInsertTorusSegment(2,180);
	 coolpipe[i]->CreateAndInsertTubeSegment(3);
       }

     } else { // On Side A tubes are all different so are created one by one

       Double_t coordL[3] = { xpos,-ypos,zpos};
       Double_t coordR[3] = {-xpos,-ypos,zpos};
       Double_t vect[3] = {0, 1, 0};
       coolpipe[0]->AddCheckPoint(moth, 0, coordL, vect);
       coolpipe[1]->AddCheckPoint(moth, 0, coordR, vect);
       coordL[1] -= kCoolPipeHeight;
       coordR[1] = coordL[1];
       coolpipe[0]->AddCheckPoint(moth, 1, coordL, vect);
       coolpipe[1]->AddCheckPoint(moth, 1, coordR, vect);
       coordL[1] -=    SinD(45) *kCoolPipeARadiusL12[0]*fgkmm;
       coordL[2] -= (1+CosD(45))*kCoolPipeARadiusL12[0]*fgkmm;
       coordR[1] -=    SinD(45) *kCoolPipeARadiusR12[0]*fgkmm;
       coordR[2] -= (1+CosD(45))*kCoolPipeARadiusR12[0]*fgkmm;
       vect[1] = TMath::Sqrt(2);
       vect[2] = -vect[1];
       coolpipe[0]->AddCheckPoint(moth, 2, coordL, vect);
       coolpipe[1]->AddCheckPoint(moth, 2, coordR, vect);
       coordL[1] += (1-CosD(45))*kCoolPipeARadiusL12[1]*fgkmm;
       coordL[2] -=    SinD(45) *kCoolPipeARadiusL12[1]*fgkmm;
       coordR[1] += (1-CosD(45))*kCoolPipeARadiusR12[1]*fgkmm;
       coordR[2] -=    SinD(45) *kCoolPipeARadiusR12[1]*fgkmm;
       vect[1] = 0;
       vect[2] = -1;
       coolpipe[0]->AddCheckPoint(moth, 3, coordL, vect);
       coolpipe[1]->AddCheckPoint(moth, 3, coordR, vect);
       coordL[2] = -kCoolPipeZSPD;
       coordR[2] = -kCoolPipeZSPD;
       coolpipe[0]->AddCheckPoint(moth, 4, coordL, vect);
       coolpipe[1]->AddCheckPoint(moth, 4, coordR, vect);

       coolpipe[0]->SetInitialNode(moth); 
       coolpipe[0]->CreateAndInsertTubeSegment(1);
       coolpipe[0]->CreateAndInsertTorusSegment(2,180);
       coolpipe[0]->CreateAndInsertTorusSegment(3,180);
       coolpipe[0]->CreateAndInsertTubeSegment(4);

       coolpipe[1]->SetInitialNode(moth); 
       coolpipe[1]->CreateAndInsertTubeSegment(1);
       coolpipe[1]->CreateAndInsertTorusSegment(2,180);
       coolpipe[1]->CreateAndInsertTorusSegment(3,180);
       coolpipe[1]->CreateAndInsertTubeSegment(4);

       zpos += kCoolManifCollDZ;

       coordL[0] = xpos; coordL[1] = -ypos; coordL[2] = zpos;
       coordR[0] =-xpos; coordR[1] = -ypos; coordR[2] = zpos;
       vect[0] = 0; vect[1] = 1; vect[2] = 0;

       coolpipe[2]->AddCheckPoint(moth, 0, coordL, vect);
       coolpipe[3]->AddCheckPoint(moth, 0, coordR, vect);
       coordL[1] -= kCoolPipeHeight;
       coordR[1] = coordL[1];
       coolpipe[2]->AddCheckPoint(moth, 1, coordL, vect);
       coolpipe[3]->AddCheckPoint(moth, 1, coordR, vect);
       coordL[1] -=    SinD(45) *kCoolPipeARadiusL34[0]*fgkmm;
       coordL[2] -= (1+CosD(45))*kCoolPipeARadiusL34[0]*fgkmm;
       coordR[1] -=    SinD(45) *kCoolPipeARadiusR34[0]*fgkmm;
       coordR[2] -= (1+CosD(45))*kCoolPipeARadiusR34[0]*fgkmm;
       vect[1] = TMath::Sqrt(2);
       vect[2] = -vect[1];
       coolpipe[2]->AddCheckPoint(moth, 2, coordL, vect);
       coolpipe[3]->AddCheckPoint(moth, 2, coordR, vect);
       coordL[1] += (1-CosD(45))*kCoolPipeARadiusL34[1]*fgkmm;
       coordL[2] -=    SinD(45) *kCoolPipeARadiusL34[1]*fgkmm;
       coordR[1] += (1-CosD(45))*kCoolPipeARadiusR34[1]*fgkmm;
       coordR[2] -=    SinD(45) *kCoolPipeARadiusR34[1]*fgkmm;
       vect[1] = 0;
       vect[2] = -1;
       coolpipe[2]->AddCheckPoint(moth, 3, coordL, vect);
       coolpipe[3]->AddCheckPoint(moth, 3, coordR, vect);
       coordL[2] = -kCoolPipeZSPD;
       coordR[2] = -kCoolPipeZSPD;
       coolpipe[2]->AddCheckPoint(moth, 4, coordL, vect);
       coolpipe[3]->AddCheckPoint(moth, 4, coordR, vect);

       coolpipe[2]->SetInitialNode(moth); 
       coolpipe[2]->CreateAndInsertTubeSegment(1);
       coolpipe[2]->CreateAndInsertTorusSegment(2,180);
       coolpipe[2]->CreateAndInsertTorusSegment(3,180);
       coolpipe[2]->CreateAndInsertTubeSegment(4);

       coolpipe[3]->SetInitialNode(moth); 
       coolpipe[3]->CreateAndInsertTubeSegment(1);
       coolpipe[3]->CreateAndInsertTorusSegment(2,180);
       coolpipe[3]->CreateAndInsertTorusSegment(3,180);
       coolpipe[3]->CreateAndInsertTubeSegment(4);

       zpos += kCoolManifCollDZ;

       coordL[0] = xpos; coordL[1] = -ypos; coordL[2] = zpos;
       coordR[0] =-xpos; coordR[1] = -ypos; coordR[2] = zpos;
       vect[0] = 0; vect[1] = 1; vect[2] = 0;

       coolpipe[4]->AddCheckPoint(moth, 0, coordL, vect);
       coolpipe[5]->AddCheckPoint(moth, 0, coordR, vect);
       coordL[1] -= kCoolPipeHeight;
       coordR[1] = coordL[1];
       coolpipe[4]->AddCheckPoint(moth, 1, coordL, vect);
       coolpipe[5]->AddCheckPoint(moth, 1, coordR, vect);
       coordL[1] -= kCoolPipeARadiusL[2]*fgkmm;
       coordL[2] -= kCoolPipeARadiusL[2]*fgkmm;
       coordR[1] -= kCoolPipeARadiusR[2]*fgkmm;
       coordR[2] -= kCoolPipeARadiusR[2]*fgkmm;
       vect[1] = 0;
       vect[2] = -1;
       coolpipe[4]->AddCheckPoint(moth, 2, coordL, vect);
       coolpipe[5]->AddCheckPoint(moth, 2, coordR, vect);
       coordL[2] = -kCoolPipeZSPD;
       coordR[2] = -kCoolPipeZSPD;
       coolpipe[4]->AddCheckPoint(moth, 3, coordL, vect);
       coolpipe[5]->AddCheckPoint(moth, 3, coordR, vect);

       coolpipe[4]->SetInitialNode(moth);
       coolpipe[4]->CreateAndInsertTubeSegment(1);
       coolpipe[4]->CreateAndInsertTorusSegment(2,180);
       coolpipe[4]->CreateAndInsertTubeSegment(3);

       coolpipe[5]->SetInitialNode(moth);
       coolpipe[5]->CreateAndInsertTubeSegment(1);
       coolpipe[5]->CreateAndInsertTorusSegment(2,180);
       coolpipe[5]->CreateAndInsertTubeSegment(3);

     } // if (sideC)

     if(GetDebug(3))
       for (Int_t i=0; i<6; i++)
	 coolpipe[i]->PrintCheckPoints();

}


//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateExtender(
    const Double_t *extenderParams, const TGeoMedium *extenderMedium,
    TArrayD& sizes) const
{
    //
    // ------------------   CREATE AN EXTENDER    ------------------------
    //
    // This function creates the following picture (in plane xOy)
    // Should be useful for the definition of the pixel bus and MCM extenders
    // The origin corresponds to point 0 on the picture, at half-width
    // in Z direction
    //
    //   Y                         7     6                      5
    //   ^                           +---+---------------------+
    //   |                          /                          |
    //   |                         /                           |
    //   0------> X               /      +---------------------+
    //                           /      / 3                     4
    //                          /      /
    //            9          8 /      /
    //            +-----------+      /
    //            |                 /
    //            |                /
    //      --->  +-----------+---+
    //      |     0          1     2
    //      |
    //  origin (0,0,0)
    //
    //
    // Takes 6 parameters in the following order :
    //   |--> par 0 : inner length [0-1] / [9-8]
    //   |--> par 1 : thickness ( = [0-9] / [4-5])
    //   |--> par 2 : angle of the slope
    //   |--> par 3 : total height in local Y direction
    //   |--> par 4 : outer length [3-4] / [6-5]
    //   |--> par 5 : width in local Z direction
    //
    Double_t slopeDeltaX = (extenderParams[3] - extenderParams[1]
                            * TMath::Cos(extenderParams[2])) /
                            TMath::Tan(extenderParams[2]);
    Double_t extenderXtruX[10] = {
        0 ,
        extenderParams[0] ,
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2]) ,
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2])+
                                                              slopeDeltaX ,
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2])+
                                           slopeDeltaX + extenderParams[4],
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2])+
                                           slopeDeltaX + extenderParams[4],
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2])+
                                                              slopeDeltaX ,
        extenderParams[0]+extenderParams[1]*TMath::Sin(extenderParams[2])+
          slopeDeltaX - extenderParams[1] * TMath::Sin(extenderParams[2]) ,
        extenderParams[0] ,
        0
    };
    Double_t extenderXtruY[10] = {
        0 ,
        0 ,
        extenderParams[1] * (1-TMath::Cos(extenderParams[2])) ,
        extenderParams[3] - extenderParams[1] ,
        extenderParams[3] - extenderParams[1] ,
        extenderParams[3] ,
        extenderParams[3] ,
        extenderParams[3]-extenderParams[1]*(1-TMath::Cos(extenderParams[2])) ,
        extenderParams[1] ,
        extenderParams[1]
    };

    if (sizes.GetSize() != 3) sizes.Set(3);
    Double_t &thickness = sizes[0];
    Double_t &length    = sizes[1];
    Double_t &width     = sizes[2];

    thickness = extenderParams[3];
    width     = extenderParams[5];
    length    = extenderParams[0]+extenderParams[1]*
            TMath::Sin(extenderParams[2])+slopeDeltaX+extenderParams[4];

    // creation of the volume
    TGeoXtru   *extenderXtru    = new TGeoXtru(2);
    TGeoVolume *extenderXtruVol = new TGeoVolume("ITSSPDextender",extenderXtru,
                                                 extenderMedium);
    extenderXtru->DefinePolygon(10,extenderXtruX,extenderXtruY);
    extenderXtru->DefineSection(0,-0.5*extenderParams[4]);
    extenderXtru->DefineSection(1, 0.5*extenderParams[4]);
    return extenderXtruVol;
}

//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreateHalfStave(Bool_t isRight,
Int_t layer,Int_t idxCentral,Int_t idxSide,TArrayD &sizes,TGeoManager *mgr)
{
    //
    // Implementation of an half-stave, which depends on the side where
    // we are on the stave. The convention for "left" and "right" is the
    // same as for the MCM. The return value is a TGeoAssembly which is
    // structured in such a way that the origin of its local reference
    // frame coincides with the origin of the whole stave.
    // The TArrayD passed by reference will contain details of the shape:
    //  - sizes[0] = thickness
    //  - sizes[1] = length
    //  - sizes[2] = width
    //  - sizes[3] = common 'x' position for eventual clips
    //  - sizes[4] = common 'y' position for eventual clips
    //  - sizes[5] = 'z' position of first clip
    //  - sizes[6] = 'z' position of second clip
    //

    // ** CHECK **

    // idxCentral and idxSide must be different
    if (idxCentral == idxSide) {
        AliInfo("Ladders must be inserted in half-stave with "
                "different indexes.");
        idxSide = idxCentral + 1;
        AliInfo(Form("Central ladder will be inserted with index %d",
                     idxCentral));
        AliInfo(Form("Side    ladder will be inserted with index %d",idxSide));
    } // end if

    // define the separations along Z direction between the objects
    Double_t sepLadderLadder = fgkmm * 0.2; // sep. btw the 2 ladders
    Double_t sepLadderCenter = fgkmm * 0.4; // sep. btw the "central" ladder
                                            // and the Z=0 plane in stave ref.
    Double_t sepLadderMCM    = fgkmm * 0.3; // sep. btw the "external" ladder
                                            // and MCM
    Double_t sepBusCenter    = fgkmm * 0.3; // sep. btw the bus central edge
                                            // and the Z=0 plane in stave ref.

    // ** VOLUMES **

    // grounding foil
    TArrayD grndSize(3);
    // This one line repalces the 3 bellow, BNS.
    TGeoVolume *grndVol = CreateGroundingFoil(isRight, grndSize, mgr);
    Double_t &grndThickness = grndSize[0];
    Double_t &grndLength = grndSize[1];

    // ladder
    TArrayD ladderSize(3);
    TGeoVolume *ladder = CreateLadder(layer, ladderSize, mgr);
    Double_t ladderThickness = ladderSize[0];
    Double_t ladderLength = ladderSize[1];
    Double_t ladderWidth = ladderSize[2];

    // MCM
    TArrayD mcmSize(3);
    TGeoVolumeAssembly *mcm = CreateMCM(!isRight,mcmSize,mgr);
    Double_t mcmThickness = mcmSize[0];
    Double_t mcmLength = mcmSize[1];
    Double_t mcmWidth = mcmSize[2];

    // bus
    TArrayD busSize(6);
    TGeoVolumeAssembly *bus = CreatePixelBus(isRight, layer, busSize, mgr);
    Double_t busThickness = busSize[0];
    Double_t busLength = busSize[1];
    Double_t busWidth = busSize[2];

    // glue between ladders and pixel bus
    TGeoMedium *medLadGlue = GetMedium("EPOXY$", mgr);
    Double_t ladGlueThickness = fgkmm * 0.1175 - fgkGapLadder;
    TGeoVolume *ladderGlue = mgr->MakeBox("ITSSPDladderGlue",medLadGlue,
                           0.5*ladGlueThickness, 0.5*busWidth, 0.5*busLength);
    ladderGlue->SetLineColor(kYellow + 5);

    // create references for the whole object, as usual
    sizes.Set(7);
    Double_t &fullThickness = sizes[0];
    Double_t &fullLength = sizes[1];
    Double_t &fullWidth = sizes[2];

    // compute the full size of the container
    fullLength    = sepLadderCenter+2.0*ladderLength+sepLadderMCM+
                       sepLadderLadder+mcmLength;
    fullWidth     = ladderWidth;
    fullThickness = grndThickness + fgkGapLadder + mcmThickness + busThickness;
    //cout << "HSTAVE FULL THICKNESS = " << fullThickness << endl;

    // ** MOVEMENTS **

    // grounding foil (shifted only along thickness)
    Double_t xGrnd = -0.5*fullThickness + 0.5*grndThickness;
    Double_t zGrnd = -0.5*grndLength;
    if (!isRight) zGrnd = -zGrnd;
    TGeoTranslation *grndTrans = new TGeoTranslation(xGrnd, 0.0, zGrnd);

    // ladders (translations along thickness and length)
    // layers must be sorted going from the one at largest Z to the
    // one at smallest Z:
    // -|Zmax| ------> |Zmax|
    //      3   2   1   0
    // then, for layer 1 ladders they must be placed exactly this way,
    // and in layer 2 at the opposite. In order to remember the placements,
    // we define as "inner" and "outer" ladder respectively the one close
    // to barrel center, and the one closer to MCM, respectively.
    Double_t xLad, zLadIn, zLadOut;
    xLad    = xGrnd + 0.5*(grndThickness + ladderThickness) +
              0.01175 - fgkGapLadder;
    zLadIn  = -sepLadderCenter - 0.5*ladderLength;
    zLadOut = zLadIn - sepLadderLadder - ladderLength;
    if (!isRight) {
        zLadIn = -zLadIn;
        zLadOut = -zLadOut;
    } // end if !isRight
    TGeoRotation *rotLad = new TGeoRotation(*gGeoIdentity);
    rotLad->RotateZ(90.0);
    rotLad->RotateY(180.0);
    Double_t sensWidth      = fgkmm * 12.800;
    Double_t chipWidth      = fgkmm * 15.950;
    Double_t guardRingWidth = fgkmm *  0.560;
    Double_t ladderShift = 0.5 * (chipWidth - sensWidth - 2.0*guardRingWidth);
    TGeoCombiTrans *trLadIn  = new TGeoCombiTrans(xLad,ladderShift,zLadIn,
                                                  rotLad);
    TGeoCombiTrans *trLadOut = new TGeoCombiTrans(xLad,ladderShift,zLadOut,
                                                  rotLad);

    // MCM (length and thickness direction, placing at same level as the
    // ladder, which implies to recompute the position of center, because
    // ladder and MCM have NOT the same thickness) the two copies of the
    // MCM are placed at the same distance from the center, on both sides
    Double_t xMCM = xGrnd + 0.5*grndThickness + 0.5*mcmThickness +
                    0.01175 - fgkGapLadder;
    Double_t yMCM = 0.5*(fullWidth - mcmWidth);
    Double_t zMCM = zLadOut - 0.5*ladderLength - 0.5*mcmLength - sepLadderMCM;
    if (!isRight) zMCM = zLadOut + 0.5*ladderLength + 0.5*mcmLength +
                         sepLadderMCM;

    // create the correction rotations
    TGeoRotation *rotMCM = new TGeoRotation(*gGeoIdentity);
    rotMCM->RotateY(90.0);
    TGeoCombiTrans *trMCM = new TGeoCombiTrans(xMCM, yMCM, zMCM, rotMCM);

    // glue between ladders and pixel bus
    Double_t xLadGlue = xLad + 0.5*ladderThickness + 0.01175 -
                        fgkGapLadder + 0.5*ladGlueThickness;

    // bus (length and thickness direction)
    Double_t xBus = xLadGlue + 0.5*ladGlueThickness + 0.5*busThickness;
    Double_t yBus  = 0.5*(fullWidth - busWidth) + 0.075; // Hardcode fix of a small overlap
    Double_t zBus = -0.5*busLength - sepBusCenter;
    if (!isRight) zBus = -zBus;
    TGeoTranslation *trBus = new TGeoTranslation(xBus, yBus, zBus);

    TGeoTranslation *trLadGlue = new TGeoTranslation(xLadGlue, 0.0, zBus);

    // create the container
    TGeoVolumeAssembly *container = 0;
    if (idxCentral+idxSide==5) {
        container = new TGeoVolumeAssembly("ITSSPDhalf-Stave1");
    } else {
        container = new TGeoVolumeAssembly("ITSSPDhalf-Stave0");
    } // end if

    // add to container all objects
    container->AddNode(grndVol, 1, grndTrans);
    // ladders are inserted in different order to respect numbering scheme
    // which is inverted when going from outer to inner layer
    container->AddNode(ladder, idxCentral+1, trLadIn);
    container->AddNode(ladder, idxSide+1, trLadOut);
    container->AddNode(ladderGlue, 1, trLadGlue);
    container->AddNode(mcm, 1, trMCM);
    container->AddNode(bus, 1, trBus);

    // since the clips are placed in correspondence of two pt1000s,
    // their position is computed here, but they are not added by default
    // it will be the StavesInSector method which will decide to add them
    // anyway, to recovery some size informations on the clip, it must be
    // created
    TArrayD clipSize;
    //    TGeoVolume *clipDummy = CreateClip(clipSize, kTRUE, mgr);
    CreateClip(clipSize, kTRUE, mgr);
    // define clip movements (width direction)
    sizes[3] = xBus + 0.5*busThickness;
    sizes[4] = 0.5 * (fullWidth - busWidth) - clipSize[6] - fgkmm*0.26;
    sizes[5] = zBus + busSize[4];
    sizes[6] = zBus + busSize[5];

    return container;
}
//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreateStave(Int_t layer,
                                    TArrayD &sizes, TGeoManager *mgr)
{
    //
    // This method uses all other ones which create pieces of the stave
    // and assemblies everything together, in order to return the whole
    // stave implementation, which is returned as a TGeoVolumeAssembly,
    // due to the presence of some parts which could generate fake overlaps
    // when put on the sector.
    // This assembly contains, going from bottom to top in the thickness
    // direction:
    //   - the complete grounding foil, defined by the "CreateGroundingFoil"
    //     method which already joins some glue and real groudning foil
    //     layers for the whole stave (left + right);
    //   - 4 ladders, which are sorted according to the ALICE numbering
    //     scheme, which depends on the layer we are building this stave for;
    //   - 2 MCMs (a left and a right one);
    //   - 2 pixel buses (a left and a right one);
    // ---
    // Arguments:
    //   - the layer number, which determines the displacement and naming
    //     of sensitive volumes
    //   - a TArrayD passed by reference which will contain the size
    //     of virtual box containing the stave
    //   - the TGeoManager
    //

    // create the container
    TGeoVolumeAssembly *container = new TGeoVolumeAssembly(Form(
                                                 "ITSSPDlay%d-Stave",layer));
    // define the indexes of the ladders in order to have the correct order
    // keeping in mind that the staves will be inserted as they are on layer
    // 2, while they are rotated around their local Y axis when inserted
    // on layer 1, so in this case they must be put in the "wrong" order
    // to turn out to be right at the end. The convention is:
    //   -|Zmax| ------> |Zmax|
    //      3   2   1   0
    // with respect to the "native" stave reference frame, "left" is in
    // the positive Z this leads the definition of these indexes:
    Int_t idxCentralL, idxSideL, idxCentralR, idxSideR;

    if (layer == 1) {
        idxSideL = 3;
        idxCentralL = 2;
        idxCentralR = 1;
        idxSideR = 0;
    } else {
        idxSideL = 0;
        idxCentralL = 1;
        idxCentralR = 2;
        idxSideR = 3;
    } // end if layer ==1

     // create the two half-staves
    TArrayD sizeL, sizeR;
    TGeoVolumeAssembly *hstaveL = CreateHalfStave(kFALSE, layer, idxCentralL,
                                             idxSideL, sizeL,mgr);
    TGeoVolumeAssembly *hstaveR = CreateHalfStave(kTRUE, layer, idxCentralR,
                                             idxSideR, sizeR, mgr);
    // copy the size to the stave's one
    sizes.Set(9);
    sizes[0] = sizeL[0];
    sizes[1] = sizeR[1] + sizeL[1];
    sizes[2] = sizeL[2];
    sizes[3] = sizeL[3];
    sizes[4] = sizeL[4];
    sizes[5] = sizeL[5];
    sizes[6] = sizeL[6];
    sizes[7] = sizeR[5];
    sizes[8] = sizeR[6];

    // add to container all objects
    container->AddNode(hstaveL, 1);
    container->AddNode(hstaveR, 1);

    return container;
}
//______________________________________________________________________
void AliITSv11GeometrySPD::SetAddStave(Bool_t *mask)
{
    //
    // Define a mask which states qhich staves must be placed.
    // It is a string which must contain '0' or '1' depending if
    // a stave must be placed or not.
    // Each place is referred to one of the staves, so the first
    // six characters of the string will be checked.
    //
     Int_t i;

     for (i = 0; i < 6; i++) fAddStave[i] = mask[i];
}
//______________________________________________________________________
void AliITSv11GeometrySPD::StavesInSector(TGeoVolume *moth, TGeoManager *mgr)
{
    //
    // Unification of essentially two methods:
    // - the one which creates the sector structure
    // - the one which returns the complete stave
    // ---
    // For compatibility, this method requires the same arguments
    // asked by "CarbonFiberSector" method, which is recalled here.
    // Like this cited method, this one does not return any value,
    // but it inserts in the mother volume (argument 'moth') all the stuff
    // which composes the complete SPD sector.
    // ---
    // In the following, the stave numbering order used for arrays is the
    // same as defined in the GetSectorMountingPoints():
    //                         /5
    //                        /\/4
    //                      1\   \/3
    //                      0|___\/2
    // ---
    // Arguments: see description of "CarbonFiberSector" method.
    //

    Double_t shift[6];  // shift from the innermost position in the
                        // sector placement plane (where the stave
                        // edge is in the point where the rounded
                        // corner begins)

    shift[0] = fgkmm * -0.691;
    shift[1] = fgkmm *  5.041;
    shift[2] = fgkmm *  1.816;
    shift[3] = fgkmm * -0.610;
    shift[4] = fgkmm * -0.610;
    shift[5] = fgkmm * -0.610;

    // corrections after interaction with Andrea and CAD
    Double_t corrX[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Double_t corrY[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    corrX[0] =  0.0046;
    corrX[1] = -0.0041;
    corrX[2] = corrX[3] = corrX[4] = corrX[5] = -0.0016;

    corrY[0] = -0.0007;
    corrY[1] = -0.0009;
    corrY[2] = corrY[3] = corrY[4] = corrY[5] = -0.0003;

    corrX[0] +=  0.00026;
    corrY[0] += -0.00080;

    corrX[1] +=  0.00018;
    corrY[1] += -0.00086;

    corrX[2] +=  0.00020;
    corrY[2] += -0.00062;

    corrX[3] +=  0.00017;
    corrY[3] += -0.00076;

    corrX[4] +=  0.00016;
    corrY[4] += -0.00096;

    corrX[5] +=  0.00018;
    corrY[5] += -0.00107;

    // create stave volumes (different for layer 1 and 2)
    TArrayD staveSizes1(9), staveSizes2(9), clipSize(5);
    Double_t &staveHeight = staveSizes1[2], &staveThickness = staveSizes1[0];
    TGeoVolume *stave1 = CreateStave(1, staveSizes1, mgr);
    TGeoVolume *stave2 = CreateStave(2, staveSizes2, mgr);
    TGeoVolume *clip   = CreateClip(clipSize, kFALSE, mgr);

    Double_t xL, yL;      // leftmost edge of mounting point (XY projection)
    Double_t xR, yR;      // rightmost edge of mounting point (XY projection)
    Double_t xM, yM;      // middle point of the segment L-R
    Double_t dx, dy;      // (xL - xR) and (yL - yR)
    Double_t widthLR;     // width of the segment L-R
    Double_t angle;       // stave rotation angle in degrees
    Double_t diffWidth;   // difference between mounting plane width and
                          // stave width (smaller)
    Double_t xPos, yPos;  // final translation of the stave
    Double_t parMovement; // translation in the LR plane direction

    staveThickness += fgkGapHalfStave;

    // loop on staves
    Int_t i, iclip = 1;
    for (i = 0; i < 6; i++) {
        // in debug mode, if this stave is not required, it is skipped
        if (!fAddStave[i]) continue;
        // retrieve reference points
        GetSectorMountingPoints(i, xL, yL, xR, yR);
        xM = 0.5 * (xL + xR);
        yM = 0.5 * (yL + yR);
        dx = xL - xR;
        dy = yL - yR;
        angle = TMath::ATan2(dy, dx);
        widthLR = TMath::Sqrt(dx*dx + dy*dy);
        diffWidth = 0.5*(widthLR - staveHeight);
        // first, a movement along this plane must be done
        // by an amount equal to the width difference
        // and then the fixed shift must also be added
        parMovement = diffWidth + shift[i];
        // due to stave thickness, another movement must be done
        // in the direction normal to the mounting plane
        // which is computed using an internal method, in a reference
        // frame where the LR segment has its middle point in the origin
        // and axes parallel to the master reference frame
        if (i == 0) {
            ParallelPosition(-0.5*staveThickness, -parMovement, angle,
                                  xPos, yPos);
        } // end if i==0
        if (i == 1) {
            ParallelPosition( 0.5*staveThickness, -parMovement, angle,
                                  xPos, yPos);
        }else {
            ParallelPosition( 0.5*staveThickness,  parMovement, angle,
                                  xPos, yPos);
        } // end if i==1
        // then we go into the true reference frame
        xPos += xM;
        yPos += yM;
        xPos += corrX[i];
        yPos += corrY[i];
        // using the parameters found here, compute the
        // translation and rotation of this stave:
        TGeoRotation *rot = new TGeoRotation(*gGeoIdentity);
        if (i == 0 || i == 1) rot->RotateX(180.0);
        rot->RotateZ(90.0 + angle * TMath::RadToDeg());
        TGeoCombiTrans *trans = new TGeoCombiTrans(xPos, yPos, 0.0, rot);
        if (i == 0 || i == 1) {
            moth->AddNode(stave1, i+1, trans);
        }else {
            moth->AddNode(stave2, i - 1, trans);
            if (i != 2) {
                // except in the case of stave #2,
                // clips must be added, and this is done directly on the sector
                Int_t j;
                //TArrayD clipSize;
                TGeoRotation *rotClip = new TGeoRotation(*gGeoIdentity);
                rotClip->RotateZ(-90.0);
                rotClip->RotateX(180.0);
                Double_t x = staveSizes2[3] + fgkGapHalfStave;
                Double_t y = staveSizes2[4];
                Double_t z[4] = { staveSizes2[5], staveSizes2[6],
                                  staveSizes2[7], staveSizes2[8] };
                for (j = 0; j < 4; j++) {
                    TGeoCombiTrans *trClip = new TGeoCombiTrans(x, y, z[j],
                                                                rotClip);
                    *trClip = *trans * *trClip;
                    moth->AddNode(clip, iclip++, trClip);
                } // end for j
            } // end if i!=2
        } // end if i==0||i==1 else
    } // end for i
    
    
    // Add a box representing the collector for cooling tubes
    // MOVED TO CreateServices() - M.S. 25 jul 12
    
}
//______________________________________________________________________
void AliITSv11GeometrySPD::ParallelPosition(Double_t dist1, Double_t dist2,
                               Double_t phi, Double_t &x, Double_t &y) const
{
    //
    // Performs the following steps:
    // 1 - finds a straight line parallel to the one passing through
    //     the origin and with angle 'phi' with X axis(phi in RADIANS);
    // 2 - finds another line parallel to the previous one, with a
    //     distance 'dist1' from it
    // 3 - takes a reference point in the second line in the intersection
    //     between the normal to both lines  passing through the origin
    // 4 - finds a point whith has distance 'dist2' from this reference,
    //     in the second line (point 2)
    // ----
    // According to the signs given to dist1 and dist2, the point is
    // found in different position w.r. to the origin
    // compute the point
    //
    Double_t cs = TMath::Cos(phi);
    Double_t sn = TMath::Sin(phi);

    x = dist2*cs - dist1*sn;
    y = dist1*cs + dist2*sn;
}
//______________________________________________________________________
Double_t AliITSv11GeometrySPD::GetSPDSectorTranslation(
    Double_t x0,Double_t y0,Double_t x1,Double_t y1,Double_t r) const
{
    //
    // Comutes the radial translation of a sector to give the
    // proper distance between SPD detectors and the beam pipe.
    // Units in are units out.
    //

    //Begin_Html
    /*
      <A HREF="http://www.physics.ohio-state.edu/HIRG/SoftWareDoc/SPD_Sector_Position.png">
      Figure showing the geometry used in the computation below. </A>
     */
    //End_Html

    // Inputs:
    //   Double_t x0  Point x0 on Sector surface for the inner
    //                most detector mounting
    //   Double_t y0  Point y0 on Sector surface for the innor
    //                most detector mounting
    //   Double_t x1  Point x1 on Sector surface for the inner
    //                most detector mounting
    //   Double_t y1  Point y1 on Sector surface for the innor
    //                most detector mounting
    //   Double_t r   The radial distance this mounting surface
    //                should be from the center of the beam pipe.
    // Outputs:
    //   none.
    // Return:
    //   The distance the SPD sector should be displaced radialy.
    //
    Double_t a,b,c;

    a = x0-x1;
    if(a==0.0) return 0.0;
    a = (y0-y1)/a;
    b = TMath::Sqrt(1.0+a*a);
    c = y0-a*x0-r*b;
    return -c;
}

//______________________________________________________________________
void AliITSv11GeometrySPD::PrintAscii(ostream *os) const
{
    //
    // Print out class data values in Ascii Form to output stream
    // Inputs:
    //   ostream *os   Output stream where Ascii data is to be writen
    // Outputs:
    //   none.
    // Return:
    //   none.
    //
    Int_t i,j,k;
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt = cout.flags();
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    *os<< fgkGapLadder <<" "<< fgkGapHalfStave<<" "<< 6 <<" ";
    for(i=0;i<6;i++) *os<< fAddStave[i] <<" "<<fSPDsectorX0.GetSize();
    for(i=0;i<fSPDsectorX0.GetSize();i++) *os<< fSPDsectorX0.GetAt(i) << " ";
    for(i=0;i<fSPDsectorX0.GetSize();i++) *os<< fSPDsectorY0.GetAt(i) << " ";
    for(i=0;i<fSPDsectorX1.GetSize();i++) *os<< fSPDsectorX1.GetAt(i) << " ";
    for(i=0;i<fSPDsectorX1.GetSize();i++) *os<< fSPDsectorY1.GetAt(i) << " ";
    *os<<10<<" "<< 2 <<" " << 6 << " "<< 3 <<" ";
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++)
        *os<<fTubeEndSector[k][0][i][j]<<" ";
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++)
        *os<<fTubeEndSector[k][1][i][j]<<" ";
    os->flags(fmt); // reset back to old Formating.
    return;
}
//
//______________________________________________________________________
void AliITSv11GeometrySPD::ReadAscii(istream* is)
{
    //
    // Read in class data values in Ascii Form to output stream
    // Inputs:
    //   istream *is   Input stream where Ascii data is to be read in from
    // Outputs:
    //   none.
    // Return:
    //   none.
    //
    Int_t i,j,k,n;
    Double_t gapLadder,gapHalfStave;
    const Int_t kLimits = 100;
    *is>>gapLadder>>gapHalfStave>>n;
    if(n!=6){
      AliError(Form("fAddStave Array !=6 n=%d",n));
        return;
    } // end if
    for(i=0;i<n;i++) *is>>fAddStave[i];
    *is>>n;
    if(n<0 || n> kLimits){
      AliError("Anomalous value for parameter n");
      return;
    } 
    fSPDsectorX0.Set(n);
    fSPDsectorY0.Set(n);
    fSPDsectorX1.Set(n);
    fSPDsectorY1.Set(n);
    for(i=0;i<n;i++) *is>>fSPDsectorX0[i];
    for(i=0;i<n;i++) *is>>fSPDsectorY0[i];
    for(i=0;i<n;i++) *is>>fSPDsectorX1[i];
    for(i=0;i<n;i++) *is>>fSPDsectorY1[i];
    *is>> i>>j>>n;
    if(i!=2||j!=6||n!=3){
        Warning("ReadAscii","fTubeEndSector array wrong size [2][6][3],"
                "found [%d][%d][%d]",i,j,n);
        return;
    } // end if
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++)
        *is>>fTubeEndSector[k][0][i][j];
    for(k=0;k<10;k++)for(i=0;i<6;i++)for(j=0;j<3;j++)
        *is>>fTubeEndSector[k][1][i][j];
    return;
}
//
//______________________________________________________________________
ostream &operator<<(ostream &os,const AliITSv11GeometrySPD &s)
{
    //
    // Standard output streaming function
    // Inputs:
    //   ostream            &os  output steam
    //   AliITSvPPRasymmFMD &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer
    //
    s.PrintAscii(&os);
    return os;
}
//
//______________________________________________________________________
istream &operator>>(istream &is,AliITSv11GeometrySPD &s)
{
    //
    // Standard inputput streaming function
    // Inputs:
    //   istream            &is  input steam
    //   AliITSvPPRasymmFMD &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer
    //
    s.ReadAscii(&is);
    return is;
}

