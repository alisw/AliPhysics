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

// AliRoot includes
#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

// Declaration file
#include "AliITSv11GeometrySPD.h"

// Constant definistions
const Double_t AliITSv11GeometrySPD::fgkGapLadder    =
                      AliITSv11Geometry::fgkmicron*75.; //  75 microns
const Double_t AliITSv11GeometrySPD::fgkGapHalfStave =
                     AliITSv11Geometry::fgkmicron*120.; // 120 microns

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
                                            TGeoManager *mgr) const
{
    //
    // This function is used to recovery any medium
    // used to build the geometry volumes.
    // If the required medium does not exists,
    // a NULL pointer is returned, and an error message is written.
    //
     Char_t itsMediumName[30];

     sprintf(itsMediumName, "ITS_%s", mediumName);
     TGeoMedium* medium = mgr->GetMedium(itsMediumName);
     if (!medium) AliError(Form("Medium <%s> not found", mediumName));

     return medium;
}
//______________________________________________________________________
Int_t AliITSv11GeometrySPD::CreateSPDCentralMaterials(Int_t &medOffset,
                                                      Int_t &matOffset) const
{
    //
    // Define the specific materials used for the ITS SPD central detectors.
    // ---
    // NOTE: These are the same old names.
    //       By the ALICE naming conventions, they start with "ITS SPD ...."
    //       Data taken from ** AliITSvPPRasymmFMD::CreateMaterials() **.
    // ---
    // Arguments [the ones passed by reference contain output values]:
    // - medOffset --> (by ref) starting number of the list of media
    // - matOffset --> (by ref) starting number of the list of Materials
    // ---
    // Inputs:
    //   Int_t &medOffset  Starting number of the list of media
    //   Int_t &matOffset  Starting number of the list of materials
    // Outputs:
    //   Int_t &medOffset  Ending number of the list of media
    //   Int_t &matOffset  Ending number of the list of materials
    // Return:
    //   The last material indexused +1. (= next avaiable material index)
    //
    const Double_t ktmaxfd    = 0.1 * fgkDegree; // Degree
    const Double_t kstemax    = 1.0 * fgkcm; // cm
    const Double_t kdeemax    = 0.1;//Fraction of particle's energy 0<deemax<=1
    const Double_t kepsil     = 1.0E-4; //
    const Double_t kstmin     = 0.0 * fgkcm; // cm "Default value used"
    const Double_t ktmaxfdAir = 0.1 * fgkDegree; // Degree
    const Double_t kstemaxAir = 1.0000E+00 * fgkcm; // cm
    const Double_t kdeemaxAir = 0.1;//Fraction of particle's energy 0<deemax<=1
    const Double_t kepsilAir  = 1.0E-4;//
    const Double_t kstminAir  = 0.0 * fgkcm; // cm "Default value used"
    const Double_t ktmaxfdSi  = 0.1 * fgkDegree; // .10000E+01; // Degree
    const Double_t kstemaxSi  = 0.0075 * fgkcm; //  .10000E+01; // cm
    const Double_t kdeemaxSi  = 0.1;//Fraction of particle's energy 0<deemax<=1
    const Double_t kepsilSi   = 1.0E-4;//
    const Double_t kstminSi   = 0.0 * fgkcm; // cm "Default value used"
    //
    Int_t matindex = matOffset;
    Int_t medindex = medOffset;
    TGeoMaterial *mat;
    TGeoMixture  *mix;
    TGeoMedium   *med;
    //
    Int_t    ifield = (((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ());
    Double_t fieldm = (((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max());
    Double_t params[8] = {8 * 0.0};

    params[1] = (Double_t) ifield;
    params[2] = fieldm;
    params[3] = ktmaxfdSi;
    params[4] = kstemaxSi;
    params[5] = kdeemaxSi;
    params[6] = kepsilSi;
    params[7] = kstminSi;

    // Definition of materials and mediums.
    // Last argument in material definition is its pressure,
    // which is initialized to ZERO.
    // For better readability, it is simply set to zero.
    // Then the writing "0.0 * fgkPascal" is replaced by "0."
    // (Alberto)

    // silicon definition for ITS (overall)
    mat = new TGeoMaterial("ITS_SI", 28.086, 14.0, 2.33 * fgkgcm3,
                           TGeoMaterial::kMatStateSolid, 25.0*fgkCelsius, 0.);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SI", medindex++, mat, params);

    // silicon for ladder chips
    mat = new TGeoMaterial("SPD SI CHIP", 28.086, 14.0, 2.33 * fgkgcm3,
                           TGeoMaterial::kMatStateSolid, 25.0*fgkCelsius, 0.);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SPD SI CHIP", medindex++, mat, params);

    // silicon for pixel bus
    mat = new TGeoMaterial("SPD SI BUS", 28.086, 14.0, 2.33 * fgkgcm3,
                           TGeoMaterial::kMatStateSolid, 25.0*fgkCelsius, 0.);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SPD SI BUS", medindex++, mat, params);

    // carbon fiber material is defined as a mix of C-O-N-H
    // defined in terms of fractional weights according to 'C (M55J)'
    // it is used for the support and clips
    mix = new TGeoMixture("C (M55J)", 4, 1.9866 * fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0, 12.01070, 6.0, 0.908508078);// C by fractional weight
    mix->DefineElement(1, 14.00670, 7.0, 0.010387573);// N by fractional weight
    mix->DefineElement(2, 15.99940, 8.0, 0.055957585);// O by fractional weight
    mix->DefineElement(3,  1.00794, 1.0, 0.025146765);// H by fractional weight
    mix->SetPressure(0.0 * fgkPascal);
    mix->SetTemperature(25.0 * fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateSolid);
    params[3] = ktmaxfd;
    params[4] = kstemax;
    params[5] = kdeemax;
    params[6] = kepsil;
    params[7] = kstmin;
    med = new TGeoMedium("ITSspdCarbonFiber", medindex++, mix, params);

    // air defined as a mixture of C-N-O-Ar:
    // it is used to fill all containers
    mix = new TGeoMixture("Air", 4, 1.20479E-3 * fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0, 12.0107,  6.0, 0.000124); // C by fractional weight
    mix->DefineElement(1, 14.0067,  7.0, 0.755267); // N by fractional weight
    mix->DefineElement(2, 15.9994,  8.0, 0.231781); // O by fractional weight
    mix->DefineElement(3, 39.9480, 18.0, 0.012827); // Ar by fractional weight
    mix->SetPressure(101325.0 * fgkPascal); // = 1 atmosphere
    mix->SetTemperature(25.0 * fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateGas);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdAir", medindex++, mix, params);

    // inox stainless steel, defined as a mixture
    // used for all metallic parts
    mix = new TGeoMixture("INOX", 9, 8.03 * fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0, 12.0107,  6., .0003);  // C  by fractional weight
    mix->DefineElement(1, 54.9380, 25., .02);    // Fe by fractional weight
    mix->DefineElement(2, 28.0855, 14., .01);    // Na by fractional weight
    mix->DefineElement(3, 30.9738, 15., .00045); // P  by fractional weight
    mix->DefineElement(4, 32.066 , 16., .0003);  // S  by fractional weight
    mix->DefineElement(5, 58.6928, 28., .12);    // Ni by fractional weight
    mix->DefineElement(6, 55.9961, 24., .17);    //    by fractional weight
    mix->DefineElement(7, 95.84  , 42., .025);   //    by fractional weight
    mix->DefineElement(8, 55.845 , 26., .654);   //    by fractional weight
    mix->SetPressure(0.0 * fgkPascal);
    mix->SetTemperature(25.0 * fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateSolid);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdStainlessSteel", medindex++, mix, params);

    // freon gas which fills the cooling system (C+F)
    mix = new TGeoMixture("Freon", 2, 1.63 * fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0, 12.0107   , 6.0,  4);  // C by fractional weight
    mix->DefineElement(1, 18.9984032, 9.0, 10); // F by fractional weight
    mix->SetPressure(101325.0 * fgkPascal); // = 1 atmosphere
    mix->SetTemperature(25.0 * fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateLiquid);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdCoolingFluid", medindex++, mix, params);

    // return the next index to be used in case of adding new materials
    medOffset = medindex;
    matOffset = matindex;
    return matOffset;
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
    TGeoVolume *vCarbonFiberSector;
    TGeoMedium *medSPDcf;

    // Define an assembly and fill it with the support of
    // a single carbon fiber sector and staves in it
    medSPDcf = GetMedium("SPD C (M55J)$", mgr);
    vCarbonFiberSector = new TGeoVolumeAssembly("ITSSPDCarbonFiberSectorV");
    vCarbonFiberSector->SetMedium(medSPDcf);
    CarbonFiberSector(vCarbonFiberSector, xAAtubeCenter0, yAAtubeCenter0, mgr);

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
        vcenteral->AddNode(vCarbonFiberSector,i+1,comrot);
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
}
//______________________________________________________________________
void AliITSv11GeometrySPD::CarbonFiberSector(TGeoVolume *moth,
     Double_t &xAAtubeCenter0, Double_t &yAAtubeCenter0, TGeoManager *mgr)
{
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
    TGeoMedium *medSPDair    = GetMedium("AIR$", mgr);
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
    const Double_t ksecX1   = -13.187 * fgkmm;
    const Double_t ksecY1   = -19.964 * fgkmm;
    const Double_t ksecR1   =  +0.6   * fgkmm; // internal  // (modif. by Alberto)
    //const Double_t ksecR1   =  +0.8   * fgkmm; // internal  // (modif. by Alberto)

    // const Double_t ksecDip0 = 5.9 * fgkmm;
    //
    //const Double_t ksecX2   =  -3.883 * fgkmm;
    const Double_t ksecX2   =  -3.833 * fgkmm; // (corr. by Alberto)
    const Double_t ksecY2   = -17.805 * fgkmm;
    const Double_t ksecR2   =  +0.6  * fgkmm; // internal (guess)
    const Double_t ksecX3   =  -3.123 * fgkmm;
    const Double_t ksecY3   = -14.618 * fgkmm;
    const Double_t ksecR3   =  -0.6   * fgkmm; // external
    //const Double_t ksecDip1 = 8.035 * fgkmm;
    //
    const Double_t ksecX4   = +11.280 * fgkmm;
    const Double_t ksecY4   = -14.473 * fgkmm;
    const Double_t ksecR4   =  +0.8   * fgkmm; // internal
    const Double_t ksecX5   = +19.544 * fgkmm;
    const Double_t ksecY5   = +10.961 * fgkmm;
    const Double_t ksecR5   =  +0.8   * fgkmm; // internal
    //const Double_t ksecDip2 = 4.553 * fgkmm;
    //
    const Double_t ksecX6   = +10.830 * fgkmm;
    const Double_t ksecY6   = +16.858 * fgkmm;
    const Double_t ksecR6   =  +0.6   * fgkmm; // internal
    const Double_t ksecX7   = +11.581 * fgkmm;
    const Double_t ksecY7   = +13.317 * fgkmm;
    const Double_t ksecR7   =  -0.6   * fgkmm; // external
    //const Double_t ksecDip3 = 6.978 * fgkmm;
    //
    const Double_t ksecX8   =  -0.733 * fgkmm;
    const Double_t ksecY8   = +17.486 * fgkmm;
    const Double_t ksecR8   =  +0.6   * fgkmm; // internal
    const Double_t ksecX9   =  +0.562 * fgkmm;
    //const Double_t ksecY9 = +14.486 * fgkmm; // correction by
    const Double_t ksecY9   = +14.107 * fgkmm; // Alberto
    const Double_t ksecR9   =  -0.6   * fgkmm; // external
    //const Double_t ksecDip4 = 6.978 * fgkmm;
    //
    const Double_t ksecX10  = -12.252 * fgkmm;
    const Double_t ksecY10  = +16.298 * fgkmm;
    const Double_t ksecR10  =  +0.6   * fgkmm; // internal
    const Double_t ksecX11  = -10.445 * fgkmm;
    const Double_t ksecY11  = +13.162 * fgkmm;
    const Double_t ksecR11  =  -0.6   * fgkmm; // external
    //const Double_t ksecDip5 = 6.978 * fgkmm;
    //
    const Double_t ksecX12  = -22.276 * fgkmm;
    const Double_t ksecY12  = +12.948 * fgkmm;
    const Double_t ksecR12  =  +0.85  * fgkmm; // internal
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
        ksecX0,  ksecX1,  -1000.0,
        ksecX2,  ksecX3,  -1000.0,
        ksecX4,  ksecX5,  -1000.0,
        ksecX6,  ksecX7,  -1000.0,
        ksecX8,  ksecX9,  -1000.0,
        ksecX10, ksecX11, -1000.0,
        ksecX12, -1000.0
    };
    Double_t secY[ksecNRadii] = {
        ksecY0,  ksecY1,  -1000.0,
        ksecY2,  ksecY3,  -1000.0,
        ksecY4,  ksecY5,  -1000.0,
        ksecY6,  ksecY7,  -1000.0,
        ksecY8,  ksecY9,  -1000.0,
        ksecY10, ksecY11, -1000.0,
        ksecY12, -1000.0
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
    /*
      Double_t secDip[ksecNRadii] = {
      0., 0., ksecDip0, 0., 0., ksecDip1,
      0., 0., ksecDip2, 0., 0., ksecDip3,
      0., 0., ksecDip4, 0., 0., ksecDip5,
      0., 0.
      };
    */
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
    Double_t secAngleStart2[ksecNRadii];
    Double_t secAngleEnd2[ksecNRadii];
    Double_t secAngleTurbo[ksecNCoolingTubeDips] = {0., 0., 0., 0., 0., 0.0};
    //Double_t secAngleStart3[ksecNRadii];
    //Double_t secAngleEnd3[ksecNRadii];
    Double_t  xpp[ksecNPoints],  ypp[ksecNPoints];
    Double_t  xpp2[ksecNPoints], ypp2[ksecNPoints];
    Double_t *xp[ksecNRadii],   *xp2[ksecNRadii];
    Double_t *yp[ksecNRadii],   *yp2[ksecNRadii];
    TGeoXtru *sA0,  *sA1, *sB0, *sB1,*sB2;
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
    sA0->SetName("ITS SPD Carbon fiber support Sector A0");
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
    sA1->SetName("ITS SPD Carbon fiber support Sector Air A1");
    sA1->DefinePolygon(m, xpp2, ypp2);
    sA1->DefineSection(0, -ksecDz);
    sA1->DefineSection(1,  ksecDz);
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
    sB0->SetName("ITS SPD Carbon fiber support Sector End B0");
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
    sB1->SetName("ITS SPD Carbon fiber support Sector Air End B1");
    sB1->DefinePolygon(i2+1, xpp2, ypp2);
    sB1->DefineSection(0,sB0->GetZ(0));
    sB1->DefineSection(1,sB0->GetZ(1)-ksecCthick2);
    const Double_t kspdEndHoleRadius1=5.698*fgkmm;
    const Double_t kspdEndHoleRadius2=2.336*fgkmm;
    const Double_t kspdEndHoleDisplacement=6.29*fgkmm;
    k = (m-1)/4;
    for(i=0;i<=k;i++){
        t= ((Double_t)i)/((Double_t)(k));
        if(!CFHolePoints(t,kspdEndHoleRadius1,kspdEndHoleRadius2,
                         kspdEndHoleDisplacement,xpp2[i],ypp2[i])){
            Warning("CarbonFiberSector","CFHolePoints failed "
                    "i=%d m=%d k=%d t=%e",i,m,k,t);
        } // end if
        // simitry in each quadrant.
        xpp2[2*k-i] = -xpp2[i];
        ypp2[2*k-i] =  ypp2[i];
        xpp2[2*k+i] = -xpp2[i];
        ypp2[2*k+i] = -ypp2[i];
        xpp2[4*k-i] =  xpp2[i];
        ypp2[4*k-i] = -ypp2[i];
    }// end for i
    //xpp2[m-1] = xpp2[0]; // begining point in
    //ypp2[m-1] = ypp2[0]; // comment with end point
    sB2 = new TGeoXtru(2);
    sB2->SetName("ITS SPD Hole in Carbon fiber support End plate");
    sB2->DefinePolygon(4*k, xpp2, ypp2);
    sB2->DefineSection(0,sB1->GetZ(1));
    sB2->DefineSection(1,sB0->GetZ(1));
    // SPD sector mount blocks
    const Double_t kMountBlock[3] = {0.5*(1.8-0.2)*fgkmm,0.5*22.0*fgkmm,
                                     0.5*45.0*fgkmm};
    sB3 = new TGeoBBox((Double_t*)kMountBlock);
    // SPD sector cooling tubes
    sTB0 = new TGeoTube("ITS SPD Cooling Tube End TB0", 0.0,
                   0.5*ksecCoolTubeROuter,0.5*(sB1->GetZ(1)-sB1->GetZ(0)));
    sTB1 = new TGeoTube("ITS SPD Cooling Tube End coolant TB0", 0.0,
                        sTB0->GetRmax() - ksecCoolTubeThick,sTB0->GetDz());
    //
    if(GetDebug(3)) {
        if(medSPDcf) medSPDcf->Dump(); else AliInfo("medSPDcf = 0");
        if(medSPDss) medSPDss->Dump(); else AliInfo("medSPDss = 0");
        if(medSPDair) medSPDair->Dump(); else AliInfo("medSPDAir = 0");
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
                                     sA0, medSPDcf);
    vA0->SetVisibility(kTRUE);
    vA0->SetLineColor(4); // Blue
    vA0->SetLineWidth(1);
    vA0->SetFillColor(vA0->GetLineColor());
    vA0->SetFillStyle(4010); // 10% transparent
    TGeoVolume *vA1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorAirA1",
                                     sA1, medSPDair);
    vA1->SetVisibility(kTRUE);
    vA1->SetLineColor(7); // light Blue
    vA1->SetLineWidth(1);
    vA1->SetFillColor(vA1->GetLineColor());
    vA1->SetFillStyle(4090); // 90% transparent
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
                                     sB0, medSPDcf);
    vB0->SetVisibility(kTRUE);
    vB0->SetLineColor(1); // Black
    vB0->SetLineWidth(1);
    vB0->SetFillColor(vB0->GetLineColor());
    vB0->SetFillStyle(4000); // 0% transparent
    TGeoVolume *vB1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndAirB1",
                                     sB1, medSPDair);
    vB1->SetVisibility(kTRUE);
    vB1->SetLineColor(0); // white
    vB1->SetLineWidth(1);
    vB1->SetFillColor(vB1->GetLineColor());
    vB1->SetFillStyle(4100); // 100% transparent
    TGeoVolume *vB2 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndAirB2",
                                     sB2, medSPDair);
    vB2->SetVisibility(kTRUE);
    vB2->SetLineColor(0); // white
    vB2->SetLineWidth(1);
    vB2->SetFillColor(vB2->GetLineColor());
    vB2->SetFillStyle(4100); // 100% transparent
    TGeoVolume *vB3 = new TGeoVolume(
        "ITSSPDCarbonFiberSupportSectorMountBlockB3",sB3, medSPDcf);
    vB3->SetVisibility(kTRUE);
    vB3->SetLineColor(1); // Black
    vB3->SetLineWidth(1);
    vB3->SetFillColor(vB3->GetLineColor());
    vB3->SetFillStyle(4000); // 0% transparent
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
    vA0->AddNode(vA1,1,0); // Put air inside carbon fiber.
    vB0->AddNode(vB1,1,0); // Put air inside carbon fiber ends.
    vB0->AddNode(vB2,1,0); // Put air wholes inside carbon fiber ends
    vTA0->AddNode(vTA1,1,0); // Put cooling liquid indide tube middel.
    vTB0->AddNode(vTB1,1,0); // Put cooling liquid inside tube end.
    Double_t tubeEndLocal[3]={0.0,0.0,sTA0->GetDz()};
    for(i = 0; i < ksecNCoolingTubeDips; i++) {
        x0 = secX3[ksecDipIndex[i]];
        y0 = secY3[ksecDipIndex[i]];
        t = 90.0 - secAngleTurbo[i];
        trans = new TGeoTranslation("",x0,y0,0.5*(sB1->GetZ(0)+sB1->GetZ(1)));
        vB1->AddNode(vTB0, i+1, trans);
        // Find location of tube ends for later use.
        trans->LocalToMaster(tubeEndLocal,fTubeEndSector[0][0][i]);
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
    rotrans = new TGeoCombiTrans("",x0,y0,-z0,rot);
    vM0->AddNode(vB3,2,rotrans); // Put Mounting bracket on sector
    /*
    j = 0; // right side, find point with largest x value
    x1 = sB0->GetX(0);
    for(i=1;i<sB0->GetNvert();i++)if(sB0->GetX(i)>x1) {j=i;x1=sB0->GetX(i);}
    j--; // Too big by 1
    //t = -TMath::RadToDeg()*TMath::ATan2(
    //                               sB0->GetX(j)-sB0->GetX(j-1),
    //                               sB0->GetY(j)-sB0->GetY(j-1));
    */
    t *= -1.0;
    rot = new TGeoRotation("",t,0.0,0.0); // z axis rotation
    /*  // this way gets correct orientation but wrong "height"
    x0 = 0.5*(sB0->GetX(j)+sB0->GetX(j-1))+
        sB3->GetDX()*TMath::Cos(t*TMath::DegToRad());
    y0 = 0.5*(sB0->GetY(j)+sB0->GetY(j-1))+
        sB3->GetDX()*TMath::Sin(t*TMath::DegToRad());
    z0 = sB0->GetZ(0)+sB3->GetDZ();
    */ // I don't understand the need for this factor 3.5.
    // posibly the SPD sector as coded isn't symetric which the
    // plans would suggest.
    x0 = -0.5*(sB0->GetX(0)+sB0->GetX(sB0->GetNvert()-1))-3.5*
        sB3->GetDX()*TMath::Cos(t*TMath::DegToRad());
    y0 = 0.5*(sB0->GetY(0)+sB0->GetY(sB0->GetNvert()-1))-3.5*
        sB3->GetDX()*TMath::Sin(t*TMath::DegToRad());
    rotrans = new TGeoCombiTrans("",1.01*x0,y0,z0,rot);
    vM0->AddNode(vB3,3,rotrans); // Put Mounting bracket on sector
    rotrans = new TGeoCombiTrans("",1.01*x0,y0,-z0,rot);
    vM0->AddNode(vB3,4,rotrans); // Put Mounting bracket on sector
    if(GetDebug(3)){
        vM0->PrintNodes();
        vA0->PrintNodes();
        vA1->PrintNodes();
        vB0->PrintNodes();
        vB1->PrintNodes();
        vB2->PrintNodes();
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
        AliError(Form("index = %d: allowed 0 --> %", index, isize));
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

/*
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateLadder
        (Int_t layer, TArrayD &sizes, TGeoManager *mgr) const
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

    // ** CRITICAL CHECK ******************************************************
    // layer number can be ONLY 1 or 2
    if (layer != 1 && layer != 2) AliFatal("Layer number MUST be 1 or 2");

    // ** MEDIA ***************************************************************

    TGeoMedium *medAir       = GetMedium("AIR$",mgr);
    TGeoMedium *medSPDSiChip = GetMedium("SPD SI CHIP$",mgr); // SPD SI CHIP
    TGeoMedium *medSi        = GetMedium("SI$",mgr);
    TGeoMedium *medBumpBond  = GetMedium("COPPER$",mgr);  // ??? BumpBond

    // ** SIZES ***************************************************************

    Double_t chipThickness  = fgkmm *  0.150;
    Double_t chipWidth      = fgkmm * 15.950;
    Double_t chipLength     = fgkmm * 13.600;
    Double_t chipSpacing    = fgkmm *  0.400; // separation of chips along Z
    Double_t sensThickness  = fgkmm *  0.200;
    Double_t sensLength     = fgkmm * 69.600;
    Double_t sensWidth      = fgkmm * 12.800;
    Double_t guardRingWidth = fgkmm *  0.560; // guard ring around sensor
    Double_t bbLength       = fgkmm * 0.042;
    Double_t bbWidth        = sensWidth;
    Double_t bbThickness    = fgkmm * 0.012;
    Double_t bbPos          = 0.080;          // Z position w.r. to left pixel edge

    // the three dimensions of the box which contains the ladder
    // are returned in the 'sizes' argument, and are used for volumes positionement
    // for readability purpose, they are linked by reference to a more meaningful name
    sizes.Set(3);
    Double_t &thickness = sizes[0];
    Double_t &length = sizes[1];
    Double_t &width = sizes[2];
    // the container is a box which exactly enclose all the stuff
    width = chipWidth;
    length = sensLength + 2.0*guardRingWidth;
    thickness = sensThickness + chipThickness + bbThickness;

    // ** VOLUMES *************************************************************

    // This is a sensitive volume.
    // Local X must correspond to x coordinate of the sensitive volume:
    // to respect this, the origin of the local reference system
    // must be shifted from the middle of the box, using
    // an additional option ('originShift') when creating the container shape:
    Double_t xSens = 0.5 * (width - sensWidth - 2.0*guardRingWidth);
    Double_t originShift[3] = {-xSens, 0., 0.};

    // now the container is a TGeoBBox with this shift,
    // and the volume is made of air (it does not exist in reality)
    TGeoBBox *shLadder = new TGeoBBox(0.5*width, 0.5*thickness, 0.5*length, originShift);
    TGeoVolume *vLadder = new TGeoVolume(Form("ITSSPDlay%d-Ladder", layer), shLadder, medAir);

    // the chip is a common box
    TGeoVolume *vChip = mgr->MakeBox("ITSSPDchip", medSPDSiChip,
                                     0.5*chipWidth, 0.5*chipThickness, 0.5*chipLength);

    // to build the sensor with its guard ring, we create a TGeoBBox with the size
    // of the sensor + guard ring, and we insert the true sensor into it as an
    // internal node: this simplifies the implementation with the same result
    TGeoVolume *vSensGuard = mgr->MakeBox(Form("%s-guardRing", GetSenstiveVolumeName(layer)),
                                          medSi,
                                          0.5*sensWidth + guardRingWidth,
                                          0.5*sensThickness,
                                          0.5*sensLength + guardRingWidth);
    TGeoVolume *vSens = mgr->MakeBox(GetSenstiveVolumeName(layer), medSi,
                                     0.5*sensWidth,0.5*sensThickness,0.5*sensLength);
    vSensGuard->AddNode(vSens, 0);
    vSensGuard->SetTransparency(50);

    // bump bond is a common box for one whole column
    TGeoVolume *vBB = mgr->MakeBox("ITSSPDbb", medBumpBond,
                                   0.5*bbWidth, 0.5*bbThickness, 0.5*bbLength);

    // set colors of all objects for visualization
    vLadder->SetLineColor(kRed);
    vSens->SetLineColor(kYellow + 1);
    vChip->SetLineColor(kGreen);
    vSensGuard->SetLineColor(kYellow + 3);
    vBB->SetLineColor(kGray);

    // ** MOVEMENTS **
    // sensor is translated along thickness (Y) and width (X)
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
    vLadder->AddNode(vSensGuard, 1, trSens);
    //vLadderAddNode(volBorder, 1, trSens);
    for (i = 0; i < 160; i++) vLadder->AddNode(vBB,i+1,trBB[i]);
    for (i = 0; i < 5; i++) vLadder->AddNode(vChip,i+3,trChip[i]);
    // return the container
    return vLadder;
}
*/

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
}//______________________________________________________________________
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
        sprintf(type,"Kap");
        break;
    case 1:
        sprintf(type,"Alu");
        break;
    case 2:
        sprintf(type,"Glue1");
        break;
    case 3:
        sprintf(type,"Glue2");
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
    TGeoBBox *shHole = 0;
    shHole = new TGeoBBox(Form("ITSSPD%sGfoilHole", type),0.5*holeLength,
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
        sprintf(name,"ITSSPDTRgFoil%sHole%d", type, i);
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
    if (isRight) strcpy(suf, "R"); else strcpy(suf, "L");
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
    if (isRight) strcpy(suf, "R"); else strcpy(suf, "L");

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

/*
//__________________________________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreatePixelBus
(Bool_t isRight, TArrayD &sizes, TGeoManager *mgr) const
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


    // ** MEDIA **

    //PIXEL BUS
    TGeoMedium *medBus     = GetMedium("SPDBUS(AL+KPT+EPOX)$",mgr);
    TGeoMedium *medPt1000  = GetMedium("CERAMICS$",mgr); // ??? PT1000
    // Capacity
    TGeoMedium *medCap     = GetMedium("SDD X7R capacitors$",mgr);
    // ??? Resistance
    // TGeoMedium *medRes     = GetMedium("SDD X7R capacitors$",mgr);
    TGeoMedium *medRes     = GetMedium("ALUMINUM$",mgr);
    TGeoMedium *medExt     = GetMedium("SDDKAPTON (POLYCH2)$", mgr);
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
    Double_t ext2Length         = fgkmm * (285.0 - ext1Length + extThickness);
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
    TGeoVolumeAssembly *container = new TGeoVolumeAssembly("PixelBus");
    TGeoVolume *bus = mgr->MakeBox("Bus", medBus, 0.5*busThickness, 0.5*busWidth, 0.5*busLength);
    TGeoVolume *pt1000 = mgr->MakeBox("PT1000", medPt1000, 0.5*pt1000Thickness, 0.5*pt1000Width, 0.5*pt1000Length);
    TGeoVolume *res = mgr->MakeBox("Resistor", medRes, 0.5*resThickness, 0.5*resWidth, 0.5*resLength);
    TGeoVolume *cap = mgr->MakeBox("Capacitor", medCap, 0.5*capThickness, 0.5*capWidth, 0.5*capLength);
    TGeoVolume *ext1 = mgr->MakeBox("Extender1", medExt, 0.5*extThickness, 0.5*extWidth, 0.5*ext1Length);
    TGeoVolume *ext2 = mgr->MakeBox("Extender2", medExt, 0.5*extHeight - extThickness, 0.5*extWidth, 0.5*extThickness);
    TGeoVolume *ext3 = mgr->MakeBox("Extender3", medExt, extThickness, 0.5*extWidth, 0.5*ext2Length);
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
    container->AddNode(bus, 0, trBus);
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
        container->AddNode(pt1000, i, tr);
    } // end for i
    // capacitors
    x = 0.5*(capThickness - fullThickness) + busThickness;
    for (i = 0; i < 2; i++) {
        y = yRef + capY[i];
        z = zRef + capZ[i];
        TGeoTranslation *tr = new TGeoTranslation(x, y, z);
        container->AddNode(cap, i, tr);
    } // end for i
    // resistors
    x = 0.5*(resThickness - fullThickness) + busThickness;
    for (i = 0; i < 2; i++) {
        y = yRef + resY[i];
        z = zRef + resZ[i];
        TGeoTranslation *tr = new TGeoTranslation(x, y, z);
        container->AddNode(res, i, tr);
    } // end for i
    // extender
    if (isRight) {
        y = 0.5 * (-fullWidth + extWidth);
        z = 0.5 * (-fullLength + fgkmm * 10.0);
    }
    else {
        y = 0.5 * (fullWidth - extWidth);
        z = 0.5 * ( fullLength - fgkmm * 10.0);
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
    x += 0.5*(extHeight - extThickness);
    TGeoTranslation *trExt2 = new TGeoTranslation(x, y, z);
    if (isRight) {
        z -= 0.5 * (ext2Length - extThickness);
    }
    else {
        z += 0.5 * (ext2Length - extThickness);
    }
    x += 0.5*(extHeight - extThickness) + extThickness;
    TGeoTranslation *trExt3 = new TGeoTranslation(x, y, z);
    container->AddNode(ext1, 0, trExt1);
    container->AddNode(ext2, 0, trExt2);
    container->AddNode(ext3, 0, trExt3);


    sizes[3] = yRef + pt1000Y;
    sizes[4] = zRef + pt1000Z[2];
    sizes[5] = zRef + pt1000Z[7];

    return container;
}
*/

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
    TGeoMedium *medExt     = GetMedium("SDDKAPTON (POLYCH2)$", mgr);
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

    TGeoVolume *ext1 = mgr->MakeBox("Extender1", medExt, 0.5*extThickness, 0.5*extWidth, 0.5*ext1Length);
    TGeoVolume *ext2 = mgr->MakeBox("Extender2", medExt, 0.5*extHeight - 2.*extThickness, 0.5*extWidth, 0.5*extThickness);
    TGeoVolume *ext3 = mgr->MakeBox("Extender3", medExt, 0.5*extThickness, 0.5*(extWidth-0.8*fgkmm), 0.5*ext2Length + extThickness); // Hardcode fix of a small overlap
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
        z -= 0.5 * (ext2Length - extThickness) + 2.5*extThickness;
    }
    else {
        z += 0.5 * (ext2Length - extThickness) + 2.5*extThickness;
    }
    x += 0.5*(extHeight - extThickness) - 2.*extThickness;
    TGeoTranslation *trExt3 = new TGeoTranslation(x, y, z);
    container->AddNode(ext1, 0, trExt1);
    container->AddNode(ext2, 0, trExt2);
    container->AddNode(ext3, 0, trExt3);

    sizes[3] = yRef + pt1000Y;
    sizes[4] = zRef + pt1000Z[2];
    sizes[5] = zRef + pt1000Z[7];

    return container;
}

//______________________________________________________________________
TList* AliITSv11GeometrySPD::CreateConeModule(TGeoManager *mgr) const
{
    TGeoMedium *medInox  = GetMedium("INOX$",mgr);
    TGeoMedium *medExt   = GetMedium("SDDKAPTON (POLYCH2)$", mgr);
    TGeoMedium *medPlate = GetMedium("SPD C (M55J)$", mgr);

    Double_t extThickness = fgkmm * 0.25;
    Double_t ext1Length   = fgkmm * (26.7 - 10.0);
    Double_t ext2Length   = fgkmm * (285.0 - ext1Length + extThickness);

    Double_t cableThickness = 1.5 * fgkmm;
    Double_t cableL1 = 350.0 * fgkmm - extThickness - ext1Length - ext2Length;
    Double_t cableL2 = 340.0 * fgkmm;
    //Double_t cableL3 = 570.0 * fgkmm;
    Double_t cableL3 = 57.0 * fgkmm;
    Double_t cableW1 =  11.0 * fgkmm;
    Double_t cableW2 =  30.0 * fgkmm;
    Double_t cableW3 =  50.0 * fgkmm;

    Double_t mcmThickness = 1.2 *fgkmm;
    Double_t mcmLength = cableL1 + cableL2 + cableL3;
    Double_t mcmWidth = cableW1;

    Double_t plateLength    = 200.0 * fgkmm;
    Double_t plateWidth     =  50.0 * fgkmm;
    Double_t plateThickness =   5.0 * fgkmm;

    Double_t x[12], y[12];

    x[0] = 7.5;
    y[0] = 0.0 + 0.5 * cableW1;

    x[1] = x[0] + cableL1 - 0.5*(cableW2 - cableW1);
    y[1] = y[0];

    x[2] = x[0] + cableL1;
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

    TGeoVolumeAssembly* container[2];
    container[0] = new TGeoVolumeAssembly("ITSSPDConeModule");
    container[1] = new TGeoVolumeAssembly("ITSSPDCoolingModule");

    TGeoXtru *shCable = new TGeoXtru(2);
    shCable->DefinePolygon(12, x, y);
    shCable->DefineSection(0, 0., 0., 0., 1.0);
    shCable->DefineSection(1, cableThickness, 0., 0., 1.0);

    TGeoVolume *volCable = new TGeoVolume("ITSSPDExtender", shCable, medExt);
    volCable->SetLineColor(kGreen);

    TGeoVolume *volTube = gGeoManager->MakeTube("ITSSPDCoolingTubeCone", medInox, 5.*fgkmm, 6.*fgkmm, 0.5*(x[5] - x[0]));
    volTube->SetLineColor(kGray);

    Double_t thickness = cableThickness + mcmThickness;
    TGeoBBox *shOut = new TGeoBBox("ITSSPD_shape_plateout", 0.5*plateThickness, 0.5*plateLength, 0.5*plateWidth);
    TGeoBBox *shIn = new TGeoBBox("ITSSPD_shape_platein", 0.5*thickness, 0.52*plateLength, 0.5*cableW2);
    Char_t string[255];
    sprintf(string, "%s-%s", shOut->GetName(), shIn->GetName());
    TGeoCompositeShape *shPlate = new TGeoCompositeShape("ITSSPDPlate_shape", string);
    TGeoVolume *volPlate = new TGeoVolume("ITSSPDPlate", shPlate, medPlate);
    volPlate->SetLineColor(kRed);

    TGeoVolume *volMCMExt = gGeoManager->MakeBox("ITSSPDextenderMCM", medExt, 0.5*mcmThickness, 0.5*mcmLength, 0.5*mcmWidth);
    volMCMExt->SetLineColor(kGreen+3);

    TGeoRotation *rot = new TGeoRotation(*gGeoIdentity);
    rot->RotateX(90.0);
    rot->RotateZ(90.0);
    container[0]->AddNode(volCable, 0, rot);

    TGeoTranslation *combi = new TGeoTranslation(cableThickness + 0.5*mcmThickness, x[0] + 0.5*mcmLength, 0.0);
    container[0]->AddNode(volMCMExt, 0, combi);

    TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
    rot1->RotateX(87.5);
    TGeoCombiTrans *tr = new TGeoCombiTrans(1.15, x[0] + 0.5*(x[5] - x[0]), -2.95, rot1);
    container[1]->AddNode(volTube, 0, tr);

    TGeoTranslation *tr1 = new TGeoTranslation(0.5*plateThickness - 0.5*(plateThickness-thickness), x[3] - x[0] - 0.52*plateLength, 0.0);
    container[0]->AddNode(volPlate, 0, tr1);

    TList* conemodulelist = new TList();

    conemodulelist->Add(container[0]);
    conemodulelist->Add(container[1]);

    return conemodulelist;
}

//______________________________________________________________________
void AliITSv11GeometrySPD::CreateCones(TGeoVolume *moth) const
{

    TList* modulelist = CreateConeModule(gGeoManager);
    TGeoVolumeAssembly* module;

    //Double_t angle[10] = {18., 54., 90., 126., 162., -18., -54., -90., -126., -162.};
    // angleNm for cone modules (cables), angleNc for cooling tubes
    Double_t angle1m[10] = {23., 53., 90., 127., 157., 203.0, 233.0, 270.0, 307.0, 337.0};
    Double_t angle2m[10] = {18., 53., 90., 126., 162., 198.0, 233.0, 270.0, 309.0, 342.0};
    Double_t angle1c[10] = {23., 53., 90., 124., 157., 203.0, 233.0, 270.0, 304.0, 337.0};
    Double_t angle2c[10] = {18., 44., 90., 126., 162., 198.0, 223.0, 270.0, 309.0, 342.0};

    // First add the cables
    module = (TGeoVolumeAssembly*)modulelist->At(0);
    for (Int_t i = 0; i < 10; i++) {
        TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
        rot1->RotateY(-90.0);
        rot1->RotateX(45.0);
	angle1m[i] -= 1.5;
        rot1->RotateZ(90.0 - angle1m[i]);
        TGeoCombiTrans *tr1 = new TGeoCombiTrans(0.0, 0.0, 38.0, rot1);
        moth->AddNode(module, 2*i, tr1);
        TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
        rot2->RotateY(90.0);
        rot2->RotateX(-45.0);
	angle2m[i] -= 1.5;
        rot2->RotateZ(90.0 - angle2m[i]);
        TGeoCombiTrans *tr2 = new TGeoCombiTrans(0.0, 0.0, -37.9, rot2);
        moth->AddNode(module, 2*i+1, tr2);
    }

    // Then the cooling tubes
    module = (TGeoVolumeAssembly*)modulelist->At(1);
    for (Int_t i = 0; i < 10; i++) {
        TGeoRotation *rot1 = new TGeoRotation(*gGeoIdentity);
        rot1->RotateY(-90.0);
        rot1->RotateX(45.0);
	angle1c[i] -= 1.5;
        rot1->RotateZ(90.0 - angle1c[i]);
        TGeoCombiTrans *tr1 = new TGeoCombiTrans(0.0, 0.0, 38.0, rot1);
        moth->AddNode(module, 2*i, tr1);
        TGeoRotation *rot2 = new TGeoRotation(*gGeoIdentity);
        rot2->RotateY(90.0);
        rot2->RotateX(-45.0);
	angle2c[i] -= 1.5;
        rot2->RotateZ(90.0 - angle2c[i]);
        TGeoCombiTrans *tr2 = new TGeoCombiTrans(0.0, 0.0, -37.9, rot2);
        moth->AddNode(module, 2*i+1, tr2);
    }
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
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreatePixelBusAndExtensions
(Bool_t /*zpos*/, TGeoManager *mgr) const
{
    //
    // Creates an assembly which contains the pixel bus and its extension
    // and the extension of the MCM.
    // By: Renaud Vernet
    // NOTE: to be defined its material and its extension in the outside
    // direction
    //
    // ====   constants   =====
    //get the media
    // PIXEL BUS
    //TGeoMedium   *medPixelBus    = GetMedium("SPDBUS(AL+KPT+EPOX)$",mgr);
    // IXEL BUS EXTENDER
    TGeoMedium *medPBExtender  = GetMedium("SDDKAPTON (POLYCH2)$",mgr);
    //MCM EXTENDER
    TGeoMedium *medMCMExtender = GetMedium("SDDKAPTON (POLYCH2)$",mgr);
    //   //geometrical constants
    const Double_t kPbextenderThickness     =   0.07 * fgkmm;
    //design=?? 70 deg. seems OK
    const Double_t kPbExtenderSlopeAngle    =  70.0  * TMath::Pi()/180.;
    // = 2.6 - (0.28+0.05+0.35) cf design
    const Double_t kPbExtenderHeight        =   1.92 * fgkmm;
    const Double_t kPbExtenderWidthY        =  11.0  * fgkmm;
    //design=?? 70 deg. seems OK
    const Double_t kMcmExtenderSlopeAngle   =  70.0  * TMath::Pi()/180.;
    const Double_t kMcmExtenderThickness    =   0.10 * fgkmm;
    const Double_t kMcmExtenderHeight       =   1.8  * fgkmm;
    const Double_t kMcmExtenderWidthY       =   kPbExtenderWidthY;
    //   const Double_t groundingThickness    =   0.07  * fgkmm;
    //   const Double_t grounding2pixelBusDz  =   0.625 * fgkmm;
    //   const Double_t pixelBusThickness     =   0.28  * fgkmm;
    //   const Double_t groundingWidthX       = 170.501 * fgkmm;
    //   const Double_t pixelBusContactDx     =   1.099 * fgkmm;
    //   const Double_t pixelBusWidthY        =  13.8   * fgkmm;
    //design=20 deg.
    //   const Double_t pixelBusContactPhi    =  20.0   * TMath::Pi()/180.
    //   const Double_t pbExtenderTopZ        =   2.72  * fgkmm;
    //   const Double_t mcmThickness          =   0.35  * fgkmm;
    //   const Double_t halfStaveTotalLength  = 247.64  * fgkmm;
    //   const Double_t deltaYOrigin          =  15.95/2.* fgkmm;
    //   const Double_t deltaXOrigin          =   1.1    * fgkmm;
    //   const Double_t deltaZOrigin          = halfStaveTotalLength / 2.;
    //   const Double_t grounding2pixelBusDz2 = grounding2pixelBusDz+
    //                           groundingThickness/2. + pixelBusThickness/2.;
    //   const Double_t pixelBusWidthX        = groundingWidthX;
    //   const Double_t pixelBusRaiseLength   = (pixelBusContactDx-
    //                  pixelBusThickness*TMath::Sin(pixelBusContactPhi))/
    //                                       TMath::Cos(pixelBusContactPhi);
    //   const Double_t pbExtenderBaseZ       = grounding2pixelBusDz2 +
    //        pixelBusRaiseLength*TMath::Sin(pixelBusContactPhi) +
    //        2*pixelBusThickness*TMath::Sin(pixelBusContactPhi)*
    //        TMath::Tan(pixelBusContactPhi);
    //   const Double_t pbExtenderDeltaZ      = pbExtenderTopZ-pbExtenderBaseZ;
    //   const Double_t pbExtenderEndPointX   = 2*deltaZOrigin -
    //    groundingWidthX - 2*pixelBusThickness*TMath::Sin(pixelBusContactPhi);
    //   const Double_t pbExtenderXtru3L   = 1.5 * fgkmm; //arbitrary ?
    //   const Double_t pbExtenderXtru4L   = (pbExtenderDeltaZ +
    //             pixelBusThickness*(TMath::Cos(extenderSlope)-2))/
    //                                      TMath::Sin(extenderSlope);
    //   const Double_t kMcmExtenderEndPointX  = deltaZOrigin - 48.2 * fgkmm;
    //   const Double_t kMcmExtenderXtru3L     = 1.5  * fgkmm;
    //   //=====  end constants  =====
    const Double_t kPbExtenderInnerLength    = 10. * fgkmm;
    const Double_t kPbExtenderOuterLength    = 15. * fgkmm;
    const Double_t kMcmExtenderInnerLength   = 10. * fgkmm;
    const Double_t kMcmExtenderOuterLength   = 15. * fgkmm;
    Double_t pbExtenderParams[6]  = {kPbExtenderInnerLength,  //0
                                     kPbextenderThickness,    //1
                                     kPbExtenderSlopeAngle,   //2
                                     kPbExtenderHeight,       //3
                                     kPbExtenderOuterLength,  //4
                                     kPbExtenderWidthY};      //5

    Double_t mcmExtenderParams[6] = {kMcmExtenderInnerLength, //0
                                     kMcmExtenderThickness,   //1
                                     kMcmExtenderSlopeAngle,  //2
                                     kMcmExtenderHeight,      //3
                                     kMcmExtenderOuterLength, //4
                                     kMcmExtenderWidthY};     //5

    TArrayD sizes(3);
    TGeoVolume* pbExtender  = CreateExtender(pbExtenderParams,medPBExtender,
                                             sizes);
    if(GetDebug(1))printf("CREATED AN EXTENDER : THICKNESS = %5.5f cm\t"
              "LENGTH=%5.5f cm\tWIDTH=%5.5f cm\n",sizes[0],sizes[1],sizes[2]);
    TGeoVolume* mcmExtender = CreateExtender(mcmExtenderParams,medMCMExtender,
                                             sizes);
    if(GetDebug(1))printf("CREATED AN EXTENDER : THICKNESS = %5.5f cm\t"
             "LENGTH=%5.5f cm\tWIDTH=%5.5f cm\n",sizes[0],sizes[1],sizes[2]);
    //   Double_t pixelBusValues[5]    = {pixelBusWidthX,        //0
    //                     pixelBusThickness,     //1
    //                     pixelBusContactPhi,    //2
    //                     pixelBusRaiseLength,   //3
    //                     pixelBusWidthY};      //4

    //   Double_t pbExtenderValues[8]  = {pixelBusRaiseLength,   //0
    //                     pixelBusContactPhi,     //1
    //                     pbExtenderXtru3L,       //2
    //                     pixelBusThickness,      //3
    //                     extenderSlope,     //4
    //                     pbExtenderXtru4L,      //5
    //                     pbExtenderEndPointX,   //6
    //                     kPbExtenderWidthY};    //7

    //   Double_t mcmExtenderValues[6] = {mcmExtenderXtru3L,     //0
    //                     mcmExtenderThickness,  //1
    //                     extenderSlope,     //2
    //                     deltaMcmMcmExtender,    //3
    //                     mcmExtenderEndPointX,  //4
    //                     mcmExtenderWidthY};    //5
    //   TGeoVolumeAssembly *pixelBus=new TGeoVolumeAssembly("ITSSPDpixelBus");
    //   CreatePixelBus(pixelBus,pixelBusValues,medPixelBus);
    //   TGeoVolumeAssembly *pbExtender = new TGeoVolumeAssembly(
    //                                              "ITSSPDpixelBusExtender");
    //   CreatePixelBusExtender(pbExtender,pbExtenderValues,medPBExtender);
    //   TGeoVolumeAssembly *mcmExtender = new TGeoVolumeAssembly(
    //                                                 "ITSSPDmcmExtender");
    //   CreateMCMExtender(mcmExtender,mcmExtenderValues,medMCMExtender);
    //--------------   DEFINITION OF GEOMETRICAL TRANSFORMATIONS --------
    //   TGeoRotation    * commonRot  = new TGeoRotation("commonRot",0,90,0);
    //   commonRot->MultiplyBy(new TGeoRotation("rot",-90,0,0));
    //   TGeoTranslation * pixelBusTrans   = new TGeoTranslation(
    //                      pixelBusThickness/2. - deltaXOrigin + 0.52*fgkmm ,
    //                                   -pixelBusWidthY/2.   + deltaYOrigin ,
    //                                   -groundingWidthX/2.  + deltaZOrigin);
    //   TGeoRotation    *pixelBusRot     = new TGeoRotation(*commonRot);
    //   TGeoTranslation *pbExtenderTrans =new TGeoTranslation(*pixelBusTrans);
    //   TGeoRotation    *pbExtenderRot   = new TGeoRotation(*pixelBusRot);
    //   pbExtenderTrans->SetDz(*(pbExtenderTrans->GetTranslation()+2) -
    //                          pixelBusWidthX/2. - 2*pixelBusThickness*
    //                                    TMath::Sin(pixelBusContactPhi));
    //   if (!zpos) {
    //     pbExtenderTrans->SetDy(*(pbExtenderTrans->GetTranslation()+1) -
    //                               (pixelBusWidthY - kPbExtenderWidthY)/2.);
    //   } else {
    //     pbExtenderTrans->SetDy(*(pbExtenderTrans->GetTranslation()+1) +
    //                            (pixelBusWidthY - kPbExtenderWidthY)/2.);
    //   }
    //   pbExtenderTrans->SetDx(*(pbExtenderTrans->GetTranslation()) +
    //                      pixelBusThickness/2 + 2*pixelBusThickness*
    //                      TMath::Sin(pixelBusContactPhi)*
    //                      TMath::Tan(pixelBusContactPhi));
    //   TGeoTranslation * mcmExtenderTrans = new TGeoTranslation(0.12*fgkmm +
    //                                    mcmThickness - deltaXOrigin,
    //                                    pbExtenderTrans->GetTranslation()[1],
    //                                    -4.82);
    //   TGeoRotation    * mcmExtenderRot   = new TGeoRotation(*pbExtenderRot);
    //   // add pt1000 components
    //   Double_t pt1000Z = fgkmm * 64400. * 1E-4;
    //   //Double_t pt1000X[10] = {319700.,  459700.,  599700.,  739700.,
    //                             879700., 1029700., 1169700., 1309700.,
    //                            1449700., 1589700.};
    //   Double_t pt1000X[10] ={66160., 206200.,  346200.,  486200.,  626200.,
    //                         776200., 916200., 1056200., 1196200., 1336200.};
    //   Double_t pt1000size[3] = {fgkmm*1.5, fgkmm*0.6, fgkmm*3.1};
    //   Int_t i;
    //   for (i = 0; i < 10; i++) {
    //     pt1000X[i] *= fgkmm * 1E-4;
    //   }
    //   TGeoVolume *pt1000 = mgr->MakeBox("ITSSPDpt1000",0,0.5*pt1000size[0],
    //                              0.5*pt1000size[1], 0.5*pt1000size[2]);
    //   pt1000->SetLineColor(kGray);
    //   Double_t refThickness = - pixelBusThickness;
    //   for (i = 0; i < 10; i++) {
    //     TGeoTranslation *tr = new TGeoTranslation(pt1000X[i]-
    //          0.5*pixelBusWidthX, 0.002+0.5*(-3.*refThickness+pt1000size[3]),
    //                                            pt1000Z -0.5*pixelBusWidthY);
    //     pixelBus->AddNode(pt1000, i+1, tr);
    //   }

    //CREATE FINAL VOLUME ASSEMBLY AND ROTATE IT
    TGeoVolumeAssembly *assembly = new TGeoVolumeAssembly("ITSSPDextenders");
    //   assembly->AddNode((TGeoVolume*)pixelBus,1,
    //          new TGeoCombiTrans(*pixelBusTrans,*pixelBusRot));
    //   assembly->AddNode((TGeoVolume*)pbExtender,1,
    //           new TGeoCombiTrans(*pbExtenderTrans,*pbExtenderRot));
    //   assembly->AddNode((TGeoVolume*)mcmExtender,1,
    //         new TGeoCombiTrans(*mcmExtenderTrans,*mcmExtenderRot));
    //   assembly->AddNode(mcmExtender,1,new TGeoIdentity());
    assembly->AddNode(pbExtender,1);
    assembly->AddNode(mcmExtender,1);
    //   assembly->SetTransparency(50);

    return assembly;
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
    sizes[4] = 0.5 * (fullWidth - busWidth) - clipSize[6] - fgkmm*0.48;
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
void AliITSv11GeometrySPD::CreateFigure0(const Char_t *filepath,
                                         const Char_t *type,
                                         TGeoManager *mgr) const
{
    //
    // Creates Figure 0 for the documentation of this class. In this
    // specific case, it creates the X,Y cross section of the SPD suport
    // section, center and ends. The output is written to a standard
    // file name to the path specificed.
    // Inputs:
    //   const Char_t *filepath  Path where the figure is to be drawn
    //   const Char_t *type      The type of file, default is gif.
    //   TGeoManager  *mgr       The TGeoManager default gGeoManager
    // Output:
    //   none.
    // Return:
    //   none.
    //
    TGeoXtru *sA0,*sA1,*sB0,*sB1;
    //TPolyMarker *pmA,*pmB;
    TPolyLine plA0,plA1,plB0,plB1;
    TCanvas *canvas;
    TLatex txt;
    Double_t x=0.0,y=0.0;
    Int_t i,kNRadii=6;

    if(strcmp(filepath,"")){
        Error("CreateFigure0","filepath=%s type=%s",filepath,type);
    } // end if
    //
    sA0 = (TGeoXtru*) mgr->GetVolume("ITSSPDCarbonFiberSupportSectorA0_1")->
              GetShape();
    sA1 = (TGeoXtru*) mgr->GetVolume("ITSSPDCarbonFiberSupportSectorAirA1_1")->
              GetShape();
    sB0 = (TGeoXtru*) mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndB0_1")->
             GetShape();
    sB1 = (TGeoXtru*) mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndAirB1_1"
           )->GetShape();
    //pmA = new TPolyMarker();
    //pmA.SetMarkerStyle(2); // +
    //pmA.SetMarkerColor(7); // light blue
    //pmB = new TPolyMarker();
    //pmB.SetMarkerStyle(5); // X
    //pmB.SetMarkerColor(6); // purple
    plA0.SetPolyLine(sA0->GetNvert());
    plA0.SetLineColor(1); // black
    plA0.SetLineStyle(1);
    plA1.SetPolyLine(sA1->GetNvert());
    plA1.SetLineColor(2); // red
    plA1.SetLineStyle(1);
    plB0.SetPolyLine(sB0->GetNvert());
    plB0.SetLineColor(3); // Green
    plB0.SetLineStyle(2);
    plB1.SetPolyLine(sB1->GetNvert());
    plB1.SetLineColor(4); // Blue
    plB1.SetLineStyle(2);
    //for(i=0;i<kNRadii;i++) pmA.SetPoint(i,xyB1p[i][0],xyB1p[i][1]);
    //for(i=0;i<kNRadii;i++) pmB.SetPoint(i,xyB1p[i][0],xyB1p[i][1]);
    for(i=0;i<sA0->GetNvert();i++) plA0.SetPoint(i,sA0->GetX(i),sA0->GetY(i));
    for(i=0;i<sA1->GetNvert();i++) plA1.SetPoint(i,sA1->GetX(i),sA1->GetY(i));
    for(i=0;i<sB0->GetNvert();i++) plB0.SetPoint(i,sB0->GetX(i),sB0->GetY(i));
    for(i=0;i<sB1->GetNvert();i++) plB1.SetPoint(i,sB1->GetX(i),sB1->GetY(i));
    canvas = new TCanvas("AliITSv11GeometrySPDFig0","",1000,1000);
    canvas->Range(-3.,-3.,3.,3.);
    txt.SetTextSize(0.05);
    txt.SetTextAlign(33);
    txt.SetTextColor(1);
    txt.DrawLatex(2.9,2.9,"Section A-A outer Carbon Fiber surface");
    txt.SetTextColor(2);
    txt.DrawLatex(2.9,2.5,"Section A-A Inner Carbon Fiber surface");
    txt.SetTextColor(3);
    txt.DrawLatex(2.9,2.1,"Section E-E outer Carbon Fiber surface");
    txt.SetTextColor(4);
    txt.DrawLatex(2.9,1.7,"Section E-E Inner Carbon Fiber surface");
    plA0.Draw();
    plA1.Draw();
    plB0.Draw();
    plB1.Draw();
    //pmA.Draw();
    //pmB.Draw();
    //
    x = 1.0;
    y = -2.5;
    Char_t chr[3];
    for(i=0;i<kNRadii;i++){
        sprintf(chr,"%2d",i);txt.DrawLatex(x-0.1,y,chr);
        sprintf(chr,"%8.4f",5.000);txt.DrawLatex(x,y,chr);
        sprintf(chr,"%8.4f",5.000);txt.DrawLatex(x+0.5,y,chr);
        sprintf(chr,"%8.4f",5.000);txt.DrawLatex(x+1.0,y,chr);
        sprintf(chr,"%8.4f",5.000);txt.DrawLatex(x+1.5,y,chr);
        sprintf(chr,"%8.4f",5.000);txt.DrawLatex(x+2.0,y,chr);
        if(kTRUE) txt.DrawLatex(x+2.5,y,"A-A/E-E");
        else txt.DrawLatex(x+2.5,y,"E-E");
    } // end for i
    txt.DrawLatex(x,y,"x_{c} mm");
    txt.DrawLatex(x+0.5,y,"y_{c} mm");
    txt.DrawLatex(x+1.0,y,"R mm");
    txt.DrawLatex(x+1.5,y,"#theta_{start}^{#circle}");
    txt.DrawLatex(x+2.0,y,"#theta_{end}^{#circle}");
    txt.DrawLatex(x+2.5,y,"Section");
    //
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
    Double_t gapLadder,GapHalfStave;

    *is>>gapLadder>>GapHalfStave>>n;
    if(n!=6){
        Warning("ReadAscii","fAddStave Array !=6 n=%d",n);
        return;
    } // end if
    for(i=0;i<n;i++) *is>>fAddStave[i];
    *is>>n;
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
//
//______________________________________________________________________
Bool_t AliITSv11GeometrySPD::Make2DCrossSections(TPolyLine &a0,TPolyLine &a1,
                             TPolyLine &b0,TPolyLine &b1,TPolyMarker &p)const
{
    //
    // Fill the objects with the points representing
    // a0 the outer carbon fiber SPD sector shape Cross Section A
    // a1 the inner carbon fiber SPD sector shape Cross Section A
    // b0 the outer carbon fiber SPD sector shape Cross Section B
    // b1 the inner carbon fiber SPD sector shape Cross Section B
    //
    // Inputs:
    //   TPolyLine &a0   The outer carbon fiber SPD sector shape
    //   TPolyLine &a1   The Inner carbon fiber SPD sector shape
    //   TPolyLine &b0   The outer carbon fiber SPD sector shape
    //   TPolyLine &b1   The Inner carbon fiber SPD sector shape
    //   TPolyMarker &p  The points where the ladders are to be placed
    // Outputs:
    //   TPolyLine &a0   The shape filled with the points
    //   TPolyLine &a1   The shape filled with the points
    //   TPolyLine &b0   The shape filled with the points
    //   TPolyLine &b1   The shape filled with the points
    //   TPolyMarker &p  The filled array of points
    // Return:
    //     An error flag.
    //
    Int_t n0,n1,i;
    Double_t x,y;
    TGeoVolume *a0V,*a1V,*b0V,*b1V;
    TGeoXtru *a0S,*a1S,*b0S,*b1S;
    TGeoManager *mgr = gGeoManager;

    a0V = mgr->GetVolume("ITS SPD Carbon fiber support Sector A0");
    a0S = dynamic_cast<TGeoXtru*>(a0V->GetShape());
    n0 = a0S->GetNvert();
    a0.SetPolyLine(n0+1);
    //for(i=0;i<fSPDsectorPoints0.GetSize();i++)
    //  printf("%d %d %d\n",i,fSPDsectorPoints0[i],fSPDsectorPoints1[i]);
    for(i=0;i<n0;i++){
        x = a0S->GetX(i);
          y = a0S->GetY(i);
          //printf("%d %g %g\n",i,x,y);
        a0.SetPoint(i,x,y);
          if(i==0) a0.SetPoint(n0,x,y);
    } // end for i
    a1V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorAirA1");
    a1S = dynamic_cast<TGeoXtru*>(a1V->GetShape());
    n1 = a1S->GetNvert();
    a1.SetPolyLine(n1+1);
    for(i=0;i<n1;i++){
        x = a1S->GetX(i);
          y = a1S->GetY(i);
        a1.SetPoint(i,x,y);
          if(i==0) a1.SetPoint(n1,x,y);
    } // end for i
    // Cross Section B
    b0V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndB0");
    b0S = dynamic_cast<TGeoXtru*>(b0V->GetShape());
    n0 = b0S->GetNvert();
    b0.SetPolyLine(n0+1);
    for(i=0;i<n0;i++){
        x = b0S->GetX(i);
          y = b0S->GetY(i);
        b0.SetPoint(i,x,y);
          if(i==0) b0.SetPoint(n0,x,y);
    } // end for i
    b1V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndAirB1");
    b1S = dynamic_cast<TGeoXtru*>(b1V->GetShape());
    n1 = b1S->GetNvert();
    b1.SetPolyLine(n1+1);
    for(i=0;i<n1;i++){
        x = b1S->GetX(i);
          y = b1S->GetY(i);
        b1.SetPoint(i,x,y);
          if(i==0) b1.SetPoint(n1,x,y);
    } // end for i
    //
    Double_t x0,y0,x1,y1;
    p.SetPolyMarker(2*fSPDsectorX0.GetSize());
    for(i=0;i<fSPDsectorX0.GetSize();i++){
          GetSectorMountingPoints(i,x0,y0,x1,y1);
          p.SetPoint(2*i,x0,y0);
          p.SetPoint(2*i+1,x1,y1);
    } // end for i
    return kTRUE;
}
