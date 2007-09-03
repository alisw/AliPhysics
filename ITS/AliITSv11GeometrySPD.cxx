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
// outside of the ceneteral volume (except for the Ceneteral support 
// cylinders. Other classes define the rest of the ITS. Specificaly the ITS
// The SSD support cone, SSD Support centeral cylinder, SDD support cone,
// The SDD cupport centeral cylinder, the SPD Thermal Sheald, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC.
//

/* $Id$ */

// General Root includes
#include <Riostream.h>
#include <TMath.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPolyLine.h>
#include <TPolyMarker.h>
// Root Geometry includes
#include <TGeoVolume.h>
#include <TGeoTube.h> // contains TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoEltu.h>
#include <TGeoXtru.h>
#include <TGeoMatrix.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoCompositeShape.h>
// AliRoot includes
#include "AliMagF.h"
#include "AliRun.h"
// Declaration file
#include "AliITSv11GeometrySPD.h"

ClassImp(AliITSv11GeometrySPD)

#define SQ(A) (A)*(A)

//______________________________________________________________________
Int_t AliITSv11GeometrySPD::CreateSPDCentralMaterials(Int_t &medOffset, Int_t &matOffset) const
{
    // Define the specific materials used for the ITS SPD central
    // detectors. Note, These are the same old names. By the ALICE
    // naming convension, these should start out at ITS SPD ....
    // This data has been taken from AliITSvPPRasymmFMD::CreateMaterials().
    // Intputs:
    //    Int_t  &medOffset   The starting number of the list of media
    //    Int_t  &matOffset   The starting number of the list of Materials
    // Outputs:
    //    Int_t  &medOffset   The ending number of the list of media
    //    Int_t  &matOffset   The ending number of the list of Materials
    // Return:
    //    the last material number used +1 (the next avaiable material number).
    //Begin_Html
    /*
      <img src="http://alice.pd.infn.it/latestdr/all-sections-module.ps"
       title="SPD Sector drawing with all cross sections defined">
       <p>The SPD Sector definition.
      <img src="http://alice.pd.infn.it/latestdr/assembly-10-modules.ps"
      titile="SPD All Sectors end view with thermal sheald">
      <p>The SPD all sector end view with thermal sheald.
      <img src="http://alice.pd.infn.it/latestdr/assembly.ps"
      title="SPD side view cross section">
      <p>SPD side view cross section with condes and thermal shealds.
      <img src="http://alice.pd.infn.it/latestdr/SECTION-A_A.jpg"
      title="Cross setion A-A"><p>Cross section A-A
      <img src="http://alice.pd.infn.it/latestdr/SECTION-B_B.jpg"
      title="Cross section B-B"><p>Cross section B-B
      <img src="http://alice.pd.infn.it/latestdr/SECTION-C_C.jpg"
      title-"Cross section C-C"><p>Cross section C-C
      <img src="http://alice.pd.infn.it/latestdr/SECTION-D_D.jpg"
      title="Cross section D-D"><p>Cross section D-D
      <img src="http://alice.pd.infn.it/latestdr/SECTION-F_F.jpg"
      title="Cross section F-F"><p>Cross section F-F
      <img src="http://alice.pd.infn.it/latestdr/SECTION-G_G.jpg"
      title="Cross section G-G"><p>Cross section G-G
     */
    //End_Html
    const Double_t ktmaxfd = 0.1*fgkDegree; // Degree
    const Double_t kstemax = 1.0*fgkcm; // cm
    const Double_t kdeemax = 0.1; // Fraction of particle's energy 0<deemax<=1
    const Double_t kepsil  = 1.0E-4; //
    const Double_t kstmin  = 0.0*fgkcm; // cm "Default value used"
    const Double_t ktmaxfdAir = 0.1*fgkDegree; // Degree
    const Double_t kstemaxAir = 1.0000E+00*fgkcm; // cm
    const Double_t kdeemaxAir = 0.1; //Fraction of particle's energy 0<deemax<=1
    const Double_t kepsilAir  = 1.0E-4;//
    const Double_t kstminAir  = 0.0*fgkcm; // cm "Default value used"
    const Double_t ktmaxfdSi = 0.1*fgkDegree; // .10000E+01; // Degree
    const Double_t kstemaxSi = 0.0075*fgkcm; //  .10000E+01; // cm
    const Double_t kdeemaxSi = 0.1; // Fraction of particle's energy 0<deemax<=1
    const Double_t kepsilSi  = 1.0E-4;//
    const Double_t kstminSi  = 0.0*fgkcm; // cm "Default value used"
    //
    Int_t matindex=matOffset;
    Int_t medindex=medOffset;
    Double_t params[8]={8*0.0};
    TGeoMaterial *mat;
    TGeoMixture  *mix;
    TGeoMedium   *med;
    //
    Int_t    ifield = (gAlice->Field()->Integ());
    Double_t fieldm = (gAlice->Field()->Max());
    params[1] = (Double_t) ifield;
    params[2] = fieldm;
    params[3] = ktmaxfdSi;
    params[4] = kstemaxSi;
    params[5] = kdeemaxSi;
    params[6] = kepsilSi;
    params[7] = kstminSi;

    mat = new TGeoMaterial("SI",28.086,14.0,2.33*fgkgcm3,
			   TGeoMaterial::kMatStateSolid,25.0*fgkCelsius,
			   0.0*fgkPascal);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SI",medindex++,mat,params);
    //med = new TGeoMedium("SI",medindex++,matindex++,0,ifield,
    //		 fieldm,ktmaxfdSi,kstemaxSi,kdeemaxSi,kepsilSi,kstminSi);
    //
    mat = new TGeoMaterial("SPD SI CHIP",28.086,14.0,2.33*fgkgcm3,
			   TGeoMaterial::kMatStateSolid,25.0*fgkCelsius,
			   0.0*fgkPascal);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SPD SI CHIP",medindex++,mat,params);
    //med = new TGeoMedium("SPD SI CHIP",medindex++,matindex++,0,ifield,
    //		 fieldm,ktmaxfdSi,kstemaxSi,kdeemaxSi,kepsilSi,kstminSi);
    //
    mat = new TGeoMaterial("SPD SI BUS",28.086,14.0,2.33*fgkgcm3,
			   TGeoMaterial::kMatStateSolid,25.0*fgkCelsius,
			   0.0*fgkPascal);
    mat->SetIndex(matindex);
    med = new TGeoMedium("SPD SI BUS",medindex++,mat,params);
    //med = new TGeoMedium("SPD SI BUS",medindex++,matindex++,0,ifield,
    //		 fieldm,ktmaxfdSi,kstemaxSi,kdeemaxSi,kepsilSi,kstminSi);
    //
    // Carbon fiber by fractional weight "C (M55J)"
    mix = new TGeoMixture("C (M55J)",4,1.9866*fgkgcm3);
    mix->SetIndex(matindex);
     // Carbon by fractional weight
    mix->DefineElement(0,12.0107,6.0,0.908508078);
    // Nitrogen by fractional weight
    mix->DefineElement(1,14.0067,7.0,0.010387573); 
    // Oxigen by fractional weight
    mix->DefineElement(2,15.9994,8.0,0.055957585); 
    // Hydrogen by fractional weight
    mix->DefineElement(3,1.00794,1.0,0.025146765); 
    mix->SetPressure(0.0*fgkPascal);
    mix->SetTemperature(25.0*fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateSolid);
    params[3] = ktmaxfd;
    params[4] = kstemax;
    params[5] = kdeemax;
    params[6] = kepsil;
    params[7] = kstmin;
    med = new TGeoMedium("ITSspdCarbonFiber",medindex++,mix,params);
    //med = new TGeoMedium("ITSspdCarbonFiber",medindex++,matindex++,0,ifield,
    //		 fieldm,ktmaxfd,kstemax,kdeemax,kepsil,kstmin);
    //
    // Carbon fiber by fractional weight
    mix = new TGeoMixture("Air",4,1.20479E-3*fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0,12.0107,6.0,0.000124); // Carbon by fractional weight
    mix->DefineElement(1,14.0067,7.0,0.755267); // Nitrogen by fractional weight
    mix->DefineElement(2,15.9994,8.0,0.231781); // Oxigen by fractional weight
    mix->DefineElement(3,39.948,18.0,0.012827); // Argon by fractional weight
    mix->SetPressure(101325.0*fgkPascal); // 1 atmosphere
    mix->SetTemperature(25.0*fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateGas);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdAir",medindex++,mix,params);
    //med = new TGeoMedium("ITSspdAir",medindex++,matindex++,0,ifield,
    //	       fieldm,ktmaxfdAir,kstemaxAir,kdeemaxAir,kepsilAir,kstminAir);
    //
    // Carbon fiber by fractional weight
    mix = new TGeoMixture("INOX",9,8.03*fgkgcm3);
    mix->SetIndex(matindex);
    mix->DefineElement(0,12.0107, 6.0,0.0003); // Carbon by fractional weight
    mix->DefineElement(1,54.9380,25.0,0.02); // Iron by fractional weight
    mix->DefineElement(2,28.0855,14.0,0.01); // Sodium by fractional weight
    mix->DefineElement(3,30.9738,15.0,0.00045); //  by fractional weight
    mix->DefineElement(4,32.066 ,16.0,0.0003); // by fractional weight
    mix->DefineElement(5,58.6928,28.0,0.12); // Nickel by fractional weight
    mix->DefineElement(6,55.9961,24.0,0.17); // by fractional weight
    mix->DefineElement(7,95.84  ,42.0,0.025); // by fractional weight
    mix->DefineElement(8,55.845 ,26.0,0.654); // by fractional weight
    mix->SetPressure(0.0*fgkPascal); //
    mix->SetTemperature(25.0*fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateSolid);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdStainlessSteel",medindex++,mix,params);
    //med =new TGeoMedium("ITSspdStainlessSteel",medindex++,matindex++,0,ifield,
    //	       fieldm,ktmaxfdAir,kstemaxAir,kdeemaxAir,kepsilAir,kstminAir);
    //
    // Carbon fiber by fractional weight
    mix->SetIndex(matindex);
    mix = new TGeoMixture("Freon",2,1.63*fgkgcm3);
    mix->DefineElement(0,12.0107,6.0,4); // Carbon by fractional weight
    mix->DefineElement(1,18.9984032,9.0,10); // Florine by fractional weight
    mix->SetPressure(101325.0*fgkPascal); // 1 atmosphere
    mix->SetTemperature(25.0*fgkCelsius);
    mix->SetState(TGeoMaterial::kMatStateLiquid);
    params[3] = ktmaxfdAir;
    params[4] = kstemaxAir;
    params[5] = kdeemaxAir;
    params[6] = kepsilAir;
    params[7] = kstminAir;
    med = new TGeoMedium("ITSspdCoolingFluid",medindex++,mix,params);
    //med = new TGeoMedium("ITSspdCoolingFluid",medindex++,matindex++,0,ifield,
    //	       fieldm,ktmaxfdAir,kstemaxAir,kdeemaxAir,kepsilAir,kstminAir);
    //
    medOffset = medindex;
    matOffset = matindex;
    return matOffset;
}
//______________________________________________________________________
void AliITSv11GeometrySPD::InitSPDCentral(Int_t offset,TVirtualMC *vmc) const {
    // Do any SPD Central detector related initilizations, setting
    // transport cuts for example.
    // Some GEANT3 Physics switches
    // "MULTS"
    // Multiple scattering. The variable IMULS controls this process. For 
    // more information see [PHYS320 or 325 or 328].
    // 0 - No multiple scattering.
    // 1 - Multiple scattering according to Moli�re theory. Default setting.
    // 2 - Same as 1. Kept for backward compatibility.
    // 3 - Pure Gaussian scattering according to the Rossi formula.
    // "DRAY"
    // delta ray production. The variable IDRAY controls this process. 
    // See [PHYS430]
    // 0 - No delta rays production.
    // 1 - delta rays production with generation of . Default setting.
    // 2 - delta rays production without generation of .
    // "LOSS"
    // Continuous energy loss. The variable ILOSS controls this process.
    // 0 - No continuous energy loss, IDRAY is set to 0.
    // 1 - Continuous energy loss with generation of delta rays above 
    //     DCUTE (common/GCUTS/) and restricted Landau fluctuations below DCUTE.
    // 2 - Continuous energy loss without generation of delta rays and full 
    //     Landau-Vavilov-Gauss fluctuations. In this case the variable IDRAY 
    //     is forced to 0 to avoid double counting of fluctuations. Default 
    //     setting.
    // 3 - Same as 1, kept for backward compatibility.
    // 4 - Energy loss without fluctuation. The value obtained from the 
    //     tables is used directly.
    // Intputs:
    //    Int_t       offset The material/medium index offset.
    //    TVirturalMC *vmc The pointer to the virtual Monte Carlo default gMC.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,n=4;

    for(i=0;i<n;i++){
      vmc->Gstpar(i+offset,"CUTGAM",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"CUTELE",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"CUTNEU",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"CUTHAD",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"CUTMUO",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"BCUTE",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"BCUTM",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"DCUTE",30.0*fgkKeV);
      vmc->Gstpar(i+offset,"DCUTM",30.0*fgkKeV);
      //vmc->Gstpar(i+offset,"PPCUTM",);
      //vmc->Gstpar(i+offset,"PAIR",);
      //vmc->Gstpar(i+offset,"COMPT",);
      //vmc->Gstpar(i+offset,"PHOT",);
      //vmc->Gstpar(i+offset,"PFIS",);
      vmc->Gstpar(i+offset,"DRAY",1);
      //vmc->Gstpar(i+offset,"ANNI",);
      //vmc->Gstpar(i+offset,"BREM",);
      //vmc->Gstpar(i+offset,"HADR",);
      //vmc->Gstpar(i+offset,"MUNU",);
      //vmc->Gstpar(i+offset,"DCAY",);
      vmc->Gstpar(i+offset,"LOSS",1);
      //vmc->Gstpar(i+offset,"MULS",);
      //vmc->Gstpar(i+offset,"GHCOR1",);
      //vmc->Gstpar(i+offset,"BIRK1",);
      //vmc->Gstpar(i+offset,"BRIK2",);
      //vmc->Gstpar(i+offset,"BRIK3",);
      //vmc->Gstpar(i+offset,"LABS",);
      //vmc->Gstpar(i+offset,"SYNC",);
      //vmc->Gstpar(i+offset,"STRA",);
    } // end for i
}
//______________________________________________________________________
void AliITSv11GeometrySPD::SPDSector(TGeoVolume *moth,TGeoManager *mgr){
    // Position of the Carbon Fiber Assembly based on distance
    // of closest point of SPD stave to beam pipe figures
    // all-sections-modules.ps of 7.22mm at section A-A.
    // Inputs:
    //   TGeoVolume *moth   the mother volume which this
    //                      object/volume is to be placed in.
    // Outputs:
    //   none.
    // Return:
    //   none.
    const Double_t kSPDclossesStaveAA    = 7.22*fgkmm;
    const Double_t kSectorStartingAngle  = -72.0*fgkDegree;
    const Double_t kNSectorsTotal        = 10.; // number
    const Double_t kSectorRelativeAngle  = 360./kNSectorsTotal*fgkDegree;
    const Double_t kBeamPipeRadius       = 0.5*60.0*fgkmm;
    //
    Int_t i;
    Double_t angle,radiusSector,xAAtubeCenter0,yAAtubeCenter0;
    Double_t staveThicknessAA=1.03*fgkmm; // get from stave geometry.
    TGeoCombiTrans *secRot=new TGeoCombiTrans();
    TGeoVolume *vCarbonFiberSector;
    TGeoMedium *medSPDcf;

    medSPDcf = mgr->GetMedium("ITSspdCarbonFiber");
    vCarbonFiberSector = new TGeoVolumeAssembly("ITSSPDCarbonFiberSectorV");
    vCarbonFiberSector->SetMedium(medSPDcf);
    CarbonFiberSector(vCarbonFiberSector,xAAtubeCenter0,yAAtubeCenter0);
    //SectorPlusStaves(vCarbonFiberSector,xAAtubeCenter0,yAAtubeCenter0);
    vCarbonFiberSector->SetVisibility(kTRUE); // logical volume
    // Compute the radial shift out of the sectors.
    radiusSector = kBeamPipeRadius+kSPDclossesStaveAA+staveThicknessAA;
    radiusSector *= radiusSector; // squaring;
    radiusSector -= xAAtubeCenter0*xAAtubeCenter0;
    radiusSector = -yAAtubeCenter0+TMath::Sqrt(radiusSector);
    angle = kSectorStartingAngle;
    secRot->RotateZ(angle);
    for(i=0;i<(Int_t)kNSectorsTotal;i++){
        secRot->SetDx(-radiusSector*TMath::Sin(angle/fgkRadian));
        secRot->SetDy(radiusSector*TMath::Cos(angle/fgkRadian));
        //secRot->RegisterYourself();
        moth->AddNode(vCarbonFiberSector,i+1,new TGeoCombiTrans(*secRot));
        if(GetDebug(5)){
            printf("i=%d angle=%g angle[rad]=%g radiusSector=%g x=%g y=%g \n",
                   i,angle,angle/fgkRadian,radiusSector,
                   -radiusSector*TMath::Sin(angle/fgkRadian),
                   radiusSector*TMath::Cos(angle/fgkRadian));
        } // end if GetDebug(5)
        angle += kSectorRelativeAngle;
        secRot->RotateZ(kSectorRelativeAngle);
    } // end for i
    if(GetDebug(3)){
        moth->PrintNodes();
    } // end if GetDebug().
    delete secRot;
}
//______________________________________________________________________
void AliITSv11GeometrySPD::CarbonFiberSector(TGeoVolume *moth,
					     Double_t &xAAtubeCenter0,
					     Double_t &yAAtubeCenter0,
					     TGeoManager *mgr){
    // Define the detail SPD Carbon fiber support Sector geometry.
    // Based on the drawings ALICE-Pixel "Construzione Profilo Modulo"
    // March 25 2004 and ALICE-SUPPORTO "construzione Profilo Modulo"
    // Define Outside radii as negitive, Outside in the sence that the
    // center of the arc is outside of the object.
    // February 16 2004.
    // Inputs:
    //   TGeoVolume *moth  The mother volume to put this object
    // Outputs:
    //  Double_t &xAAtubeCenter0  The x location of the outer surface
    //                            of the cooling tube center for tube 0.
    //                            This location helps determine where 
    //                            this sector is to be located (information
    //                            used for this is the distance the
    //                            center of the #0 detector is from the
    //                            beam pipe. Measurements taken at 
    //                            cross section A-A.
    //  Double_t &yAAtubeCenter0  The y location of the outer surface
    //                            of the cooling tube center for tube 0
    //                            This location helps determine where
    //                            this sector is to be located (information
    //                            used for this is the distance the 
    //                            center of the #0 detector is from the
    //                            beam pipe. Measurements taken at 
    //                            cross section A-A.
    //   TGeoManager *mgr         The TGeoManager as needed, default is
    //                            gGeoManager.
    // Return:
    //  none.
    TGeoMedium *medSPDcf  = 0; // SPD support cone Carbon Fiber materal number.
    //TGeoMedium *medSPDfs  = 0; // SPD support cone inserto stesalite 4411w.
    //TGeoMedium *medSPDfo  = 0; // SPD support cone foam, Rohacell 50A.
    TGeoMedium *medSPDss  = 0; // SPD support cone screw material,Stainless
    TGeoMedium *medSPDair = 0; // SPD support cone Air
    //TGeoMedium *medSPDal  = 0; // SPD support cone SDD mounting bracket Al
    TGeoMedium *medSPDcoolfl  = 0; // SPD cooling fluid, Freeon
    medSPDcf = mgr->GetMedium("ITSspdCarbonFiber");
    //medSPDfs = mgr->GetMedium("ITSspdStaselite4411w");
    //medSPDfo = mgr->GetMedium("ITSspdRohacell50A");
    medSPDss = mgr->GetMedium("ITSspdStainlessSteel");
    medSPDair= mgr->GetMedium("ITSspdAir");
    medSPDcoolfl= mgr->GetMedium("ITSspdCoolingFluid");
    //
    const Double_t ksecDz        = 0.5*500.0*fgkmm;
    const Double_t ksecLen       = 30.0*fgkmm;
    const Double_t ksecCthick    = 0.20*fgkmm;
    const Double_t ksecDipLength = 3.2*fgkmm;
    const Double_t ksecDipRadii  = 0.4*fgkmm;
    //const Double_t ksecCoolingTubeExtraDepth = 0.86*fgkmm;
    // These positions, ksecX*,ksecY* are the center of curvatures
    // for the different point around the SPD sector. The radii,
    // inner and outer, are the radous of curvature about the centers
    // ksecX* and ksecY*. To draw this SPD sector, first plot all of
    // the ksecX and ksecY points and draw circles of the specified
    // radius about these points. Connect the circles, such that the
    // lines are tangent to the circles, in accordance with the
    // radii being "Inside" or "Outside". These lines and the 
    // corresponding arc's are the surface of this SPD sector.
    const Double_t ksecX0   = -10.725*fgkmm;
    const Double_t ksecY0   = -14.853*fgkmm;
    const Double_t ksecR0   = -0.8*fgkmm; // Outside
    const Double_t ksecX1   = -13.187*fgkmm;
    const Double_t ksecY1   = -19.964*fgkmm;
    const Double_t ksecR1   = +0.6*fgkmm; // Inside
    //const Double_t ksecDip0 = 5.9*fgkmm;
    //
    const Double_t ksecX2   = -3.883*fgkmm;
    const Double_t ksecY2   = -17.805*fgkmm;
    const Double_t ksecR2   = +0.80*fgkmm; // Inside Guess. 
    const Double_t ksecX3   = -3.123*fgkmm;
    const Double_t ksecY3   = -14.618*fgkmm;
    const Double_t ksecR3   = -0.6*fgkmm; // Outside
    //const Double_t ksecDip1 = 8.035*fgkmm;
    //
    const Double_t ksecX4   = +11.280*fgkmm;
    const Double_t ksecY4   = -14.473*fgkmm;
    const Double_t ksecR4   = +0.8*fgkmm; // Inside
    const Double_t ksecX5   = +19.544*fgkmm;
    const Double_t ksecY5   = +10.961*fgkmm;
    const Double_t ksecR5   = +0.8*fgkmm; // Inside
    //const Double_t ksecDip2 = 4.553*fgkmm;
    //
    const Double_t ksecX6   = +10.830*fgkmm;
    const Double_t ksecY6   = +16.858*fgkmm;
    const Double_t ksecR6   = +0.6*fgkmm; // Inside
    const Double_t ksecX7   = +11.581*fgkmm;
    const Double_t ksecY7   = +13.317*fgkmm;
    const Double_t ksecR7   = -0.6*fgkmm; // Outside
    //const Double_t ksecDip3 = 6.978*fgkmm;
    //
    const Double_t ksecX8   = -0.733*fgkmm;
    const Double_t ksecY8   = +17.486*fgkmm;
    const Double_t ksecR8   = +0.6*fgkmm; // Inside
    const Double_t ksecX9   = +0.562*fgkmm;
    //const Double_t ksecY9   = +14.486*fgkmm;  // correction by
    const Double_t ksecY9   = +14.107*fgkmm;    // Alberto
    const Double_t ksecR9   = -0.6*fgkmm; // Outside
    //const Double_t ksecDip4 = 6.978*fgkmm;
    //
    const Double_t ksecX10  = -12.252*fgkmm;
    const Double_t ksecY10  = +16.298*fgkmm;
    const Double_t ksecR10  = +0.6*fgkmm; // Inside
    const Double_t ksecX11  = -10.445*fgkmm;
    const Double_t ksecY11  = +13.162*fgkmm;
    const Double_t ksecR11  = -0.6*fgkmm; // Outside
    //const Double_t ksecDip5 = 6.978*fgkmm;
    //
    const Double_t ksecX12  = -22.276*fgkmm;
    const Double_t ksecY12  = +12.948*fgkmm;
    const Double_t ksecR12  = +0.85*fgkmm; // Inside
    //const Double_t ksecX13 = *fgkmm;
    //const Double_t ksecY13 = *fgkmm;
    const Double_t ksecR13  = -0.8*fgkmm; // Outside
    const Double_t ksecAngleSide13 = 36.0*fgkDegree;
    //
    const Int_t ksecNRadii = 20;
    const Int_t ksecNPointsPerRadii = 4;
    const Int_t ksecNCoolingTubeDips = 6;
    // Since the Rounded parts are aproximated by a regular polygon and
    // a cooling tube of the propper diameter must fit, a scaling factor
    // increases the size of the polygon for the tube to fit.
    //const Double_t ksecRCoolScale = 1./TMath::Cos(TMath::Pi()/
    //                                          (Double_t)ksecNPointsPerRadii);
    const Double_t ksecZEndLen  = 30.00*fgkmm;
    //const Double_t ksecZFlangLen= 45.00*fgkmm;
    const Double_t ksecTl       = 0.860*fgkmm;
    const Double_t ksecCthick2  = 0.600*fgkmm;
    //const Double_t ksecCthick3  = 1.800*fgkmm;
    //const Double_t ksecSidelen  = 22.00*fgkmm;
    //const Double_t ksecSideD5   = 3.679*fgkmm;
    //const Double_t ksecSideD12  = 7.066*fgkmm;
    const Double_t ksecRCoolOut = 2.400*fgkmm;
    const Double_t ksecRCoolIn  = 2.000*fgkmm;
    const Double_t ksecDl1      = 5.900*fgkmm;
    const Double_t ksecDl2      = 8.035*fgkmm;
    const Double_t ksecDl3      = 4.553*fgkmm;
    const Double_t ksecDl4      = 6.978*fgkmm;
    const Double_t ksecDl5      = 6.978*fgkmm;
    const Double_t ksecDl6      = 6.978*fgkmm;
    const Double_t ksecCoolTubeThick  = 0.04*fgkmm;
    const Double_t ksecCoolTubeROuter = 2.6*fgkmm;
    const Double_t ksecCoolTubeFlatX  = 3.696*fgkmm;
    const Double_t ksecCoolTubeFlatY  = 0.68*fgkmm;
    //const Double_t ksecBeamX0   = 0.0*fgkmm; // guess
    //const Double_t ksecBeamY0   = (15.223+40.)*fgkmm; // guess
    //
    const Int_t ksecNPoints = (ksecNPointsPerRadii+1)*ksecNRadii + 8;
    Double_t secX[ksecNRadii] = {ksecX0,ksecX1,-1000.0,ksecX2 ,ksecX3 ,-1000.0,
				 ksecX4,ksecX5,-1000.0,ksecX6 ,ksecX7 ,-1000.0,
				 ksecX8,ksecX9,-1000.0,ksecX10,ksecX11,-1000.0,
				 ksecX12,-1000.0};
    Double_t secY[ksecNRadii] = {ksecY0,ksecY1,-1000.0,ksecY2 ,ksecY3 ,-1000.0,
				 ksecY4,ksecY5,-1000.0,ksecY6 ,ksecY7 ,-1000.0,
				 ksecY8,ksecY9,-1000.0,ksecY10,ksecY11,-1000.0,
				 ksecY12,-1000.0};
    Double_t secR[ksecNRadii] ={ksecR0 ,ksecR1 ,-.5*ksecDipLength-ksecDipRadii,
				ksecR2 ,ksecR3 ,-.5*ksecDipLength-ksecDipRadii,
				ksecR4 ,ksecR5 ,-.5*ksecDipLength-ksecDipRadii,
				ksecR6 ,ksecR7 ,-.5*ksecDipLength-ksecDipRadii,
				ksecR8 ,ksecR9 ,-.5*ksecDipLength-ksecDipRadii,
				ksecR10,ksecR11,-.5*ksecDipLength-ksecDipRadii,
				ksecR12,ksecR13};/*
    Double_t secDip[ksecNRadii]={0.0,0.0,ksecDip0,0.0,0.0,ksecDip1,
				 0.0,0.0,ksecDip2,0.0,0.0,ksecDip3,
				 0.0,0.0,ksecDip4,0.0,0.0,ksecDip5,
				 0.0,0.0};*/
    Double_t secX2[ksecNRadii];
    Double_t secY2[ksecNRadii];
    Double_t secR2[ksecNRadii] = {
	ksecR0,ksecR1,ksecRCoolOut,ksecR2,ksecR3,ksecRCoolOut,ksecR4,ksecR5,
	ksecRCoolOut,ksecR6,ksecR7,ksecRCoolOut,ksecR8,ksecR9,ksecRCoolOut,
	ksecR10,ksecR11,ksecRCoolOut,ksecR12,ksecR13};
    Double_t secDip2[ksecNCoolingTubeDips]={ksecDl1,ksecDl2,ksecDl3,
					    ksecDl4,ksecDl5,ksecDl6};
    Double_t secX3[ksecNRadii];
    Double_t secY3[ksecNRadii];
    const Int_t ksecDipIndex[ksecNCoolingTubeDips] = {2,5,8,11,14,17};
    Double_t secAngleStart[ksecNRadii];
    Double_t secAngleEnd[ksecNRadii];
    Double_t secAngleStart2[ksecNRadii];
    Double_t secAngleEnd2[ksecNRadii];
    Double_t secAngleTurbo[ksecNCoolingTubeDips] = {0.0,0.0,0.0,0.0,0.0,0.0};
    //Double_t secAngleStart3[ksecNRadii];
    //Double_t secAngleEnd3[ksecNRadii];
    Double_t xpp[ksecNPoints],ypp[ksecNPoints];
    Double_t xpp2[ksecNPoints],ypp2[ksecNPoints];
    Double_t *xp[ksecNRadii],*xp2[ksecNRadii];
    Double_t *yp[ksecNRadii],*yp2[ksecNRadii];
    TGeoXtru *sA0,*sA1,*sB0,*sB1;
    TGeoEltu *sTA0,*sTA1;
    TGeoTube *sTB0,*sTB1; //,*sM0;
    TGeoRotation    *rot;
    TGeoTranslation *trans;
    TGeoCombiTrans  *rotrans;
    Double_t t,t0,t1,a,b,x0,y0,x1,y1;
    Int_t i,j,k,m;
    Bool_t tst;

    if(moth==0){
	Error("CarbonFiberSector","moth=%p",moth);
	return;
    } // end if moth==0
    //SetDebug(3);
    for(i=0;i<ksecNRadii;i++){
	xp[i]  = &(xpp[i*(ksecNPointsPerRadii+1)]);
	yp[i]  = &(ypp[i*(ksecNPointsPerRadii+1)]);
	xp2[i] = &(xpp2[i*(ksecNPointsPerRadii+1)]);
	yp2[i] = &(ypp2[i*(ksecNPointsPerRadii+1)]);
	secX2[i] = secX[i];
	secY2[i] = secY[i];
	secX3[i] = secX[i];
	secY3[i] = secY[i];
    } // end for i

    // Find starting and ending angles for all but cooling tube sections
    secAngleStart[0] = 0.5*ksecAngleSide13;
    for(i=0;i<ksecNRadii-2;i++){
	tst = kFALSE;
	for(j=0;j<ksecNCoolingTubeDips;j++) tst = tst||i==ksecDipIndex[j];
	if(tst) continue;
	tst = kFALSE;
	for(j=0;j<ksecNCoolingTubeDips;j++) tst = tst||(i+1)==ksecDipIndex[j];
	if(tst) j = i+2;
	else j = i+1;
	AnglesForRoundedCorners(secX[i],secY[i],secR[i],
				secX[j],secY[j],secR[j],t0,t1);
	secAngleEnd[i]   = t0;
	secAngleStart[j] = t1;
	if(secR[i]>0.0&&secR[j]>0.0)if(secAngleStart[i]>secAngleEnd[i])
	    secAngleEnd[i] += 360.0;
	secAngleStart2[i] = secAngleStart[i];
	secAngleEnd2[i]   = secAngleEnd[i];
    } // end for i
    secAngleEnd[ksecNRadii-2]   = secAngleStart[ksecNRadii-2] + 
				     (secAngleEnd[ksecNRadii-5]-
				      secAngleStart[ksecNRadii-5]);
    if(secAngleEnd[ksecNRadii-2]<0.0) secAngleEnd[ksecNRadii-2] += 360.0;
    secAngleStart[ksecNRadii-1] = secAngleEnd[ksecNRadii-2] - 180.0;
    secAngleEnd[ksecNRadii-1]   = secAngleStart[0];
    secAngleStart2[ksecNRadii-2] = secAngleStart[ksecNRadii-2];
    secAngleEnd2[ksecNRadii-2]   = secAngleEnd[ksecNRadii-2];
    secAngleStart2[ksecNRadii-1] = secAngleStart[ksecNRadii-1];
    secAngleEnd2[ksecNRadii-1]   = secAngleEnd[ksecNRadii-1];
    // Find location of circle last rounded corner.
    i = 0;
    j = ksecNRadii-2;
    t0 = TanD(secAngleStart[i]-90.);
    t1 = TanD(secAngleEnd[j]-90.);
    t  = secY[i] - secY[j];
    // Note, secR[i=0] <0; secR[j=18]>0; and secR[j+1=19] <0
    t += (-secR[i]+secR[j+1])*SinD(secAngleStart[i]);
    t -= (secR[j]-secR[j+1])*SinD(secAngleEnd[j]);
    t += t1*secX[j] - t0*secX[i];
    t += t1*(secR[j]-secR[j+1])*CosD(secAngleEnd[j]);
    t -= t0*(-secR[i]+secR[j+1])*CosD(secAngleStart[i]);
    secX[ksecNRadii-1] = t/(t1-t0);
    secY[ksecNRadii-1] = TanD(90.+0.5*ksecAngleSide13)*
			  (secX[ksecNRadii-1]-secX[0]) + secY[0];
    secX2[ksecNRadii-1] = secX[ksecNRadii-1];
    secY2[ksecNRadii-1] = secY[ksecNRadii-1];
    secX3[ksecNRadii-1] = secX[ksecNRadii-1];
    secY3[ksecNRadii-1] = secY[ksecNRadii-1];
    // find location of cooling tube centers
    for(i=0;i<ksecNCoolingTubeDips;i++){
	j = ksecDipIndex[i];
	x0 = secX[j-1] + TMath::Abs(secR[j-1])*CosD(secAngleEnd[j-1]);
	y0 = secY[j-1] + TMath::Abs(secR[j-1])*SinD(secAngleEnd[j-1]);
	x1 = secX[j+1] + TMath::Abs(secR[j+1])*CosD(secAngleStart[j+1]);
	y1 = secY[j+1] + TMath::Abs(secR[j+1])*SinD(secAngleStart[j+1]);
	t0 = TMath::Sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	t  = secDip2[i]/t0;
	a  = x0+(x1-x0)*t;
	b  = y0+(y1-y0)*t;
	if(i==0){ // get location of tube center->Surface for locating
		  // this sector around the beam pipe. This needs to be
	          // double checked, but I need my notes for that, Bjorn Nilsen
	    xAAtubeCenter0 = x0+(x1-x0)*t*0.5;
	    yAAtubeCenter0 = y0+(y1-y0)*t*0.5;
	} // end if i==0
	if(a+b*(a-x0)/(b-y0)>0.0){
	    secX[j] = a + TMath::Abs(y1-y0)*2.0*ksecDipRadii/t0;
	    secY[j] = b - TMath::Sign(2.0*ksecDipRadii,y1-y0)*(x1-x0)/t0;
	    secX2[j] = a + TMath::Abs(y1-y0)*ksecTl/t0;
	    secY2[j] = b - TMath::Sign(ksecTl,y1-y0)*(x1-x0)/t0;
	    secX3[j] = a + TMath::Abs(y1-y0)*(2.0*ksecDipRadii-
					  0.5*ksecCoolTubeFlatY)/t0;
	    secY3[j] = b - TMath::Sign(2.0*ksecDipRadii-0.5*ksecCoolTubeFlatY,
				   y1-y0)*(x1-x0)/t0;
	}else{
	    secX[j] = a - TMath::Abs(y1-y0)*2.0*ksecDipRadii/t0;
	    secY[j] = b + TMath::Sign(2.0*ksecDipRadii,y1-y0)*(x1-x0)/t0;
	    secX2[j] = a - TMath::Abs(y1-y0)*ksecTl/t0;
	    secY2[j] = b + TMath::Sign(ksecTl,y1-y0)*(x1-x0)/t0;
	    secX3[j] = a - TMath::Abs(y1-y0)*(2.0*ksecDipRadii-
					  0.5*ksecCoolTubeFlatY)/t0;
	    secY3[j] = b + TMath::Sign(2.0*ksecDipRadii-0.5*ksecCoolTubeFlatY,
				      y1-y0)*(x1-x0)/t0;
        } // end if
        // Set up Start and End angles to correspond to start/end of dips.
        t1 = (secDip2[i]-TMath::Abs(secR[j]))/t0;
        secAngleStart[j] = TMath::RadToDeg()*TMath::ATan2(
                               y0+(y1-y0)*t1-secY[j],x0+(x1-x0)*t1-secX[j]);
        if(secAngleStart[j]<0.0) secAngleStart[j] += 360.0;
        secAngleStart2[j] = secAngleStart[j];
        t1 = (secDip2[i]+TMath::Abs(secR[j]))/t0;
        secAngleEnd[j] = TMath::RadToDeg()*TMath::ATan2(
                               y0+(y1-y0)*t1-secY[j],x0+(x1-x0)*t1-secX[j]);
        if(secAngleEnd[j]<0.0) secAngleEnd[j] += 360.0;
        secAngleEnd2[j]   = secAngleEnd[j];
        if(secAngleEnd[j]>secAngleStart[j]) secAngleEnd[j] -= 360.0;
        secR[j] = TMath::Sqrt(secR[j]*secR[j]+4.0*ksecDipRadii*ksecDipRadii);
    } // end for i
    // Special cases
    secAngleStart2[8] -= 360.;
    secAngleStart2[11] -= 360.;
    //
    fSPDsectorPoints0.Set(ksecNCoolingTubeDips);
    fSPDsectorPoints1.Set(ksecNCoolingTubeDips);
    //
    for(i=0;i<ksecNCoolingTubeDips;i++){
        // Find index in xpp[] and ypp[] corresponding to where the
        // SPD ladders are to be attached. Order them according to
        // the ALICE numbering schema. Using array of indexes (+-1 for
        // cooling tubes. For any "bend/dip/edge, there are 
        // ksecNPointsPerRadii+1 points involved.
        if(i==0) j=1;
        else if(i==1) j=0;
	else j=i;
	fSPDsectorPoints0[i] = (ksecDipIndex[j]-1)*(ksecNPointsPerRadii+1)+
                               (ksecNPointsPerRadii);
	fSPDsectorPoints1[i] = (ksecDipIndex[j]+1)*(ksecNPointsPerRadii+1);
    } // end for i
    SPDsectorShape(ksecNRadii,secX,secY,secR,secAngleStart,secAngleEnd,
                   ksecNPointsPerRadii,m,xp,yp);
    //  Fix up dips to be square.
    for(i=0;i<ksecNCoolingTubeDips;i++){
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
        for(k=2;k<ksecNPointsPerRadii-1;k++){// extra points spread them out.
            t = ((Double_t)(k-1))*t0;
            xp[j][k] = x0+(x1-x0)*t;
            yp[j][k] = y0+(y1-y0)*t;
        } // end for k
        secAngleTurbo[i] = -TMath::RadToDeg()*TMath::ATan2(y1-y0,x1-x0);
        if(GetDebug(3)){ 
           cout <<"i="<<i<<" angle="<<secAngleTurbo[i]<<" x0,y0{"
                <<x0<<","<<y0<<"} x1y1={"<<x1<<","<<y1<<"}"<<endl;
        } // end if
    } // end for i
    sA0 = new TGeoXtru(2);
    // This shape needs to be access later to mount the SPD sector to.
    //fSPDsectorShapeName = "ITS SPD Carbon fiber support Sector A0";
    //sA0->SetName(fSPDsectorShapeName.Data());
    sA0->SetName("ITS SPD Carbon fiber support Sector A0");
    sA0->DefinePolygon(m,xpp,ypp);
    sA0->DefineSection(0,-ksecDz);
    sA0->DefineSection(1,ksecDz);
    //
    //printf("SectorA#%d ",0);
    InsidePoint(xpp[m-1],ypp[m-1],xpp[0],ypp[0],xpp[1],ypp[1],
                ksecCthick,xpp2[0],ypp2[0]);
    for(i=1;i<m-1;i++){
        j = i/(ksecNPointsPerRadii+1);
        //printf("SectorA#%d ",i);
        InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],
                    ksecCthick,xpp2[i],ypp2[i]);
    } // end for i
    //printf("SectorA#%d ",m);
    InsidePoint(xpp[m-2],ypp[m-2],xpp[m-1],ypp[m-1],xpp[0],ypp[0],
                ksecCthick,xpp2[m-1],ypp2[m-1]);
    // Fix center value of cooling tube dip.
    // find location of cooling tube centers
    for(i=0;i<ksecNCoolingTubeDips;i++){
        j = ksecDipIndex[i];
        x0 = xp2[j][1];
        y0 = yp2[j][1];
        x1 = xp2[j][ksecNPointsPerRadii-1];
        y1 = yp2[j][ksecNPointsPerRadii-1];
        t0 = TMath::Sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
        t  = secDip2[i]/t0;
        for(k=2;k<ksecNPointsPerRadii-1;k++){// extra points spread them out.
            t = ((Double_t)(k-1))*t0;
            xp2[j][k] = x0+(x1-x0)*t;
            yp2[j][k] = y0+(y1-y0)*t;
        } // end for k
    } // end for i
    sA1 = new TGeoXtru(2);
    sA1->SetName("ITS SPD Carbon fiber support Sector Air A1");
    sA1->DefinePolygon(m,xpp2,ypp2);
    sA1->DefineSection(0,-ksecDz);
    sA1->DefineSection(1,ksecDz);
    //
    // Error in TGeoEltu. Semi-axis X must be < Semi-axis Y (?).
    sTA0 = new TGeoEltu("ITS SPD Cooling Tube TA0",
                      0.5* ksecCoolTubeFlatY, 0.5* ksecCoolTubeFlatX,ksecDz);
    sTA1 = new TGeoEltu("ITS SPD Cooling Tube coolant TA1",
                        sTA0->GetA()-ksecCoolTubeThick,
                        sTA0->GetB()-ksecCoolTubeThick,ksecDz);
    //
    SPDsectorShape(ksecNRadii,secX2,secY2,secR2,secAngleStart2,secAngleEnd2,
                   ksecNPointsPerRadii,m,xp,yp);
    //
    sB0 = new TGeoXtru(2);
    sB0->SetName("ITS SPD Carbon fiber support Sector End B0");
    sB0->DefinePolygon(m,xpp,ypp);
    sB0->DefineSection(0,ksecDz);
    sB0->DefineSection(1,ksecDz+ksecZEndLen);
    //
    //printf("SectorB#%d ",0);
    InsidePoint(xpp[m-1],ypp[m-1],xpp[0],ypp[0],xpp[1],ypp[1],
                ksecCthick2,xpp2[0],ypp2[0]);
    for(i=1;i<m-1;i++){
        t = ksecCthick2;
        for(k=0;k<ksecNCoolingTubeDips;k++)
            if((i/(ksecNPointsPerRadii+1))==ksecDipIndex[k]) 
                if(!(ksecDipIndex[k]*(ksecNPointsPerRadii+1)==i || 
                     ksecDipIndex[k]*(ksecNPointsPerRadii+1)+
                     ksecNPointsPerRadii==i   )) 
                    t = ksecRCoolOut-ksecRCoolIn;
        //printf("SectorB#%d ",i);
        InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],
                    t,xpp2[i],ypp2[i]);
    } // end for
    //printf("SectorB#%d ",m);
    InsidePoint(xpp[m-2],ypp[m-2],xpp[m-1],ypp[m-1],xpp[0],ypp[0],
                ksecCthick2,xpp2[m-1],ypp2[m-1]);
    sB1 = new TGeoXtru(2);
    sB1->SetName("ITS SPD Carbon fiber support Sector Air End B1");
    sB1->DefinePolygon(m,xpp2,ypp2);
    sB1->DefineSection(0,ksecDz);
    sB1->DefineSection(1,ksecDz+ksecLen);
    sTB0 = new TGeoTube("ITS SPD Cooling Tube End TB0",0.0,
                       0.5*ksecCoolTubeROuter,0.5*ksecLen);
    sTB1 = new TGeoTube("ITS SPD Cooling Tube End coolant TB0",0.0,
                       sTB0->GetRmax()-ksecCoolTubeThick,0.5*ksecLen);
    //
    //sM0 = new TGeoTube("ITS SPD Sensitive Virutual Volume M0",0.0,8.0,
    //                   sA0->GetZ(1)+sB0->GetZ(1));
    //
    if(GetDebug(3)){
        if(medSPDcf) medSPDcf->Dump();
        else printf("medSPDcf=0\n");
        if(medSPDss) medSPDss->Dump();
        else printf("medSPDss=0\n");
        if(medSPDair) medSPDair->Dump();
        else printf("medSPDAir=0\n");
        if(medSPDcoolfl) medSPDcoolfl->Dump();
        else printf("medSPDcoolfl=0\n");
        //sM0->InspectShape();
        sA0->InspectShape();
        sA1->InspectShape();
        sB0->InspectShape();
        sB1->InspectShape();
    } // end if GetDebug
    //
    TGeoVolume *vA0,*vA1,*vTA0,*vTA1,*vB0,*vB1,*vTB0,*vTB1;
    TGeoVolumeAssembly *vM0;
    vM0 = new TGeoVolumeAssembly("ITSSPDSensitiveVirtualvolumeM0");
    //vM0 = new TGeoVolume("ITSSPDSensitiveVirtualvolumeM0",sM0,medSPDair);
    //vM0->SetVisibility(kTRUE);
    //vM0->SetLineColor(7); // light Blue
    //vM0->SetLineWidth(1);
    //vM0->SetFillColor(vM0->GetLineColor());
    //vM0->SetFillStyle(4090); // 90% transparent
    // ALBERTO
    fSPDsectorShapeName = "ITSSPDCarbonFiberSupportSectorA0";
    vA0 = new TGeoVolume(fSPDsectorShapeName,sA0,medSPDcf);
    //vA0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorA0",sA0,medSPDcf);
    vA0->SetVisibility(kTRUE);
    vA0->SetLineColor(4); // Blue
    vA0->SetLineWidth(1);
    vA0->SetFillColor(vA0->GetLineColor());
    vA0->SetFillStyle(4010); // 10% transparent
    vA1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorAirA1",sA1,medSPDair);
    vA1->SetVisibility(kTRUE);
    vA1->SetLineColor(7); // light Blue
    vA1->SetLineWidth(1);
    vA1->SetFillColor(vA1->GetLineColor());
    vA1->SetFillStyle(4090); // 90% transparent
    vTA0 = new TGeoVolume("ITSSPDCoolingTubeTA0",sTA0,medSPDss);
    vTA0->SetVisibility(kTRUE);
    vTA0->SetLineColor(1); // Black
    vTA0->SetLineWidth(1);
    vTA0->SetFillColor(vTA0->GetLineColor());
    vTA0->SetFillStyle(4000); // 0% transparent
    vTA1 = new TGeoVolume("ITSSPDCoolingTubeFluidTA1",sTA1,medSPDcoolfl);
    vTA1->SetVisibility(kTRUE);
    vTA1->SetLineColor(6); // Purple
    vTA1->SetLineWidth(1);
    vTA1->SetFillColor(vTA1->GetLineColor());
    vTA1->SetFillStyle(4000); // 0% transparent
    vB0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndB0",sB0,medSPDcf);
    vB0->SetVisibility(kTRUE);
    vB0->SetLineColor(4); // Blue
    vB0->SetLineWidth(1);
    vB0->SetFillColor(vB0->GetLineColor());
    vB0->SetFillStyle(4010); // 10% transparent
    vB1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndAirB1",
                         sB1,medSPDair);
    vB1->SetVisibility(kTRUE);
    vB1->SetLineColor(7); // light Blue
    vB1->SetLineWidth(1);
    vB1->SetFillColor(vB1->GetLineColor());
    vB1->SetFillStyle(4090); // 90% transparent
    vTB0 = new TGeoVolume("ITSSPDCoolingTubeEndTB0",sTB0,medSPDss);
    vTB0->SetVisibility(kTRUE);
    vTB0->SetLineColor(1); // Black
    vTB0->SetLineWidth(1);
    vTB0->SetFillColor(vTB0->GetLineColor());
    vTB0->SetFillStyle(4000); // 0% transparent
    vTB1 = new TGeoVolume("ITSSPDCoolingTubeEndFluidTB1",sTB1,medSPDcoolfl);
    vTB1->SetVisibility(kTRUE);
    vTB1->SetLineColor(6); // Purple
    vTB1->SetLineWidth(1);
    vTB1->SetFillColor(vTB1->GetLineColor());
    vTB1->SetFillStyle(4000); // 0% transparent
    //
    StavesInSector(vM0);
    moth->AddNode(vM0,1,0); // Add virtual volume to mother
    vA0->AddNode(vA1,1,0); // Put air inside carbon fiber.
    vB0->AddNode(vB1,1,0); // Put air inside carbon fiber.
    vTA0->AddNode(vTA1,1,0); // Put air inside carbon fiber.
    vTB0->AddNode(vTB1,1,0); // Put air inside carbon fiber.
    for(i=0;i<ksecNCoolingTubeDips;i++){
        x0 = secX3[ksecDipIndex[i]];
        y0 = secY3[ksecDipIndex[i]];
        t = 90.0-secAngleTurbo[i];
        trans = new TGeoTranslation("",x0,y0,0.5*(sB1->GetZ(0)+sB1->GetZ(1)));
        vB1->AddNode(vTB0,i+1,trans);
        rot = new TGeoRotation("",0.0,0.0,t);
        rotrans = new TGeoCombiTrans("",x0,y0,0.0,rot);
        vM0->AddNode(vTA0,i+1,rotrans);
        //delete rot; // rot owned by AliITSv11GeometerySPD::CarbonFiberSector
    } // end for i
    vM0->AddNode(vA0,1,0);
    vM0->AddNode(vB0,1,0);
    // Reflection.
    vM0->AddNode(vB0,2,new TGeoRotation("",90.,0.,90.,90.,180.,0.));
    if(GetDebug(3)){
        vM0->PrintNodes();
        vA0->PrintNodes();
        vA1->PrintNodes();
        vB0->PrintNodes();
        vB1->PrintNodes();
        vTA0->PrintNodes();
        vTA1->PrintNodes();
        vTB0->PrintNodes();
        vTB1->PrintNodes();
    } // end if GetDebug
    //
}
//----------------------------------------------------------------------
Bool_t AliITSv11GeometrySPD::GetSectorMountingPoints(Int_t index,
			     Double_t &x0,Double_t &y0,
			     Double_t &x1,Double_t &y1,TGeoManager *mgr)const{
    // Return's the mounting locations needed to mount the SPD ladders 
    // on the SPD Carbon fiber Sectors (A cross section). Coordinate
    // system is that of the carbon fiber sector TVolume 
    // "ITSSPDCarbonFiberSupportSectorA0". Index numbering is as follows
    //                         /5
    //                        /\/4
    //                      1\   \/3
    //                      0|___\/2
    // Inputs:
    //    Int_t index   the index for which location on the SPD sector [0-5]
    // Outputs:
    //    Double_t &x0     The x0 location or the ladder sector [cm]
    //    Double_t &y0     The y0 location of the ladder sector [cm]
    //    Double_t &x1     The x1 location or the ladder sector [cm]
    //    Double_t &y1     The y1 location of the ladder sector [cm]
    //    TGeoManager *mgr The Geometry manager to use [gGeoManager]
    // Return:
    //     Returns kTRUE if no problems incountered. Returns kFALSE
    //     if a problem was incountered (for example the shape has 
    //     not been found.
    TGeoVolume *spdSectorV=0;
    TGeoXtru *spdSector=0;
    Int_t ixy0,ixy1;

    x0 = x1 = y0 = y1 = 0.0;
    if(index<0 || index>fSPDsectorPoints0.GetSize()){
      Error("GetSectorMountingPoints","index=%d size=%d",index,
	    fSPDsectorPoints0.GetSize());
      return kFALSE;
    }// end if
    spdSectorV = mgr->GetVolume(fSPDsectorShapeName.Data());
    if(spdSectorV==0){
      Error("GetSectorMountingPoints","spdSectorV==0 name=%s",
	    fSPDsectorShapeName.Data());
      return kFALSE;
    } // end if
    spdSector = dynamic_cast<TGeoXtru*>(spdSectorV->GetShape());
    if(spdSector==0){
      Error("GetSectorMountingPoints","spdSector==0");
      return kFALSE;
    } // end if
    ixy0 = fSPDsectorPoints0.At(index);
    ixy1 = fSPDsectorPoints1.At(index);
    x0 = spdSector->GetX(ixy0);
    y0 = spdSector->GetY(ixy0);
    x1 = spdSector->GetX(ixy1);
    y1 = spdSector->GetY(ixy1);
    return kTRUE;
}
//----------------------------------------------------------------------
void AliITSv11GeometrySPD::SPDsectorShape(Int_t n,const Double_t *xc,
const Double_t *yc,const Double_t *r,const Double_t *ths,const Double_t *the,
                               Int_t npr,Int_t &m,Double_t **xp,Double_t **yp){
    // Code to compute the points that make up the shape of the SPD
    // Carbon fiber support sections
    // Inputs:
    //    Int_t    n       Size of arrays xc,yc, and r.
    //    Double_t *xc     Array of x values for radii centers.
    //    Double_t *yc     Array of y values for radii centers.
    //    Double_t *r      Array of signed radii values.
    //    Double_t *ths    Array of starting angles [degrees].
    //    Double_t *the    Array of ending angles [degrees].
    //    Int_t    npr     The number of lines segments to aproximate the arc.
    // Outputs:
    //    Int_t    m       The number of enetries in the arrays *xp[npr+1] 
    //                     and *yp[npr+1].
    //    Double_t **xp    Array of x coordinate values of the line segments
    //                     which make up the SPD support sector shape.
    //    Double_t **yp    Array of y coordinate values of the line segments
    //                     which make up the SPD support sector shape.
    // Return:
    //    none.
    Int_t i,k;
    Double_t t,t0,t1;

    m = n*(npr+1);
    if(GetDebug(2)){
        cout <<"    X    \t  Y  \t  R  \t  S  \t  E"<< m <<endl;
        for(i=0;i<n;i++){
            cout <<"{"<< xc[i] <<",";
            cout << yc[i] <<",";
            cout << r[i] <<",";
            cout << ths[i] <<",";
            cout << the[i] <<"},"<< endl;
        } // end for i
    } // end if GetDebug
    //
    if(GetDebug(3)) cout <<"Double_t sA0 = ["<< n*(npr+1)+1<<"][";
    if(GetDebug(4)) cout <<"3]{";
    else if(GetDebug(3)) cout <<"2]{";
    t0 = (Double_t)npr;
    for(i=0;i<n;i++){
        t1 = (the[i]-ths[i])/t0;
        if(GetDebug(5)) cout<<"t1="<< t1<<endl;
        for(k=0;k<=npr;k++){
            t=ths[i]+((Double_t)k)*t1;
            xp[i][k] = TMath::Abs(r[i])*CosD(t)+xc[i];
            yp[i][k] = TMath::Abs(r[i])*SinD(t)+yc[i];
            if(GetDebug(3)){
                cout << "{"<<xp[i][k]<<","<<yp[i][k];
                if(GetDebug(4)) cout <<","<<t;
                cout <<"},";
            } // end if GetDebug
        } // end for k
        if(GetDebug(3)) cout << endl;
    } // end of i
    if(GetDebug(3)) cout<<"{"<<xp[0][0]<<","<<yp[0][0];
    if(GetDebug(4)) cout<<","<< ths[0];
    if(GetDebug(3)) cout<<"}}"<<endl;
    //
    return;
}
//----------------------------------------------------------------------
void AliITSv11GeometrySPD::CreateFigure0(const Char_t *filepath,
                                         const Char_t *type,
					 TGeoManager *mgr){
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
    sA0 = (TGeoXtru*) mgr->GetVolume(
        "ITSSPDCarbonFiberSupportSectorA0_1")->GetShape();
    sA1 = (TGeoXtru*) mgr->GetVolume(
        "ITSSPDCarbonFiberSupportSectorAirA1_1")->GetShape();
    sB0 = (TGeoXtru*) mgr->GetVolume(
        "ITSSPDCarbonFiberSupportSectorEndB0_1")->GetShape();
    sB1 = (TGeoXtru*) mgr->GetVolume(
        "ITSSPDCarbonFiberSupportSectorEndAirB1_1")->GetShape();
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
TGeoVolume* AliITSv11GeometrySPD::CreateLadder(Int_t layer,Double_t &length,
		     Double_t &width,Double_t &thickness, TGeoManager *mgr){
	// Creates the "ladder" = silicon sensor + 5 chips.
	// All parts are implemented as TGeoBBox and inserted 
	// into a container which is the return value of this method.
	// The sizes of the components come from drawings 
	// of the Technical office of INFN Padova.
	// Due to the requirement to specify the sensitive volume 
        // separately from the rest, the sensor is implemented as the 
        // sum of a central sensitive part + a guard ring.
	// Also the bump-bondings are added in form of small cylinders.
	// ---
	// Arguments:
	//  - the layer which will own this ladder (MUST be 1 or 2)
	//  - the used TGeoManager
	// ---
	// Returns:
	//  - the container TGeoBBox (return value)
	//  - the size of the container box (arguments passed by reference)
	// ---
	// NOTE 1
	// Here and in the other methods which contribute to the stave 
        // definition a convention is used for the naming of the three 
        // dimensions of the volumes:
	//  - 'length'    refers to the size in the Z direction of the 
        //                ALICE reference frame
	//  - 'width'     refers to the "large" dimension orthogonal to 
        //                Z axis in the local reference frame of the 
        //                object being implemented (e.g., 15.95 mm for 
        //                the chips)
	//  - 'thickness' refers to the "small" dimension orthogonal to 
        //                Z axis, which is also the direction along which 
        //                the components are superimposed on each other
	// ---
	// NOTE 2
	// all sizes taken are expressed in mm in drawings, and this is 
        // kept as is, to avoid confusion the conversion is made 
        // multiplying by the conversion factor
	
	// ** CRITICAL CHECK **
	// layer number can be ONLY 1 or 2
	if (layer != 1 && layer != 2) AliFatal("Layer number MUST be 1 or 2");
	
	// instantiate all required media
	TGeoMedium *medAir       = mgr->GetMedium("Air");
	TGeoMedium *medSPDSiChip = mgr->GetMedium("SPD SI CHIP");
	TGeoMedium *medSi        = mgr->GetMedium("Si");
	TGeoMedium *medBumpBond  = mgr->GetMedium("BumpBond");
	
	// ** Define sizes **
	// they are expressed in mm in the drawings so they require conversion
	// 'length'    is in the direction of the detector length (Z axis)
	// 'thickness' is obvious
	// 'width'     is in the direction orthogonal to 'width' and 'thickness'

	// for the chip, also the spacing between them is required
	Double_t chipThickness  = fgkmm *  0.150;
	Double_t chipWidth      = fgkmm * 15.950;
	Double_t chipLength     = fgkmm * 13.600;
	Double_t chipSpacing    = fgkmm *  0.400;
	
	// for the sensor, we define the area of sensitive volume
	// while the guard ring is added as a separate piece
	Double_t sensThickness  = fgkmm *  0.200;
	Double_t sensLength     = fgkmm * 69.600;
	Double_t sensWidth      = fgkmm * 13.920;
	Double_t guardRingWidth = fgkmm *  0.560;
	
	// bump bond is defined as a small stripe of height = 0.012 mm
	// and a suitable width to keep the same volume it has 
	// before being compressed (a line of spheres of 0.025 mm radius)
	Double_t bbLength    = fgkmm * 0.042;
	Double_t bbWidth     = sensWidth;
	Double_t bbThickness = fgkmm * 0.012;
	Double_t bbPos       = 0.080;   // Z position w.r. to left pixel edge
		
	// ** Create volumes **
	// the container is the return value, and is built as a box
	// whose edges exactly enclose the stuff we inserted here, 
	// filled with air. Its name depends on the layer number.
	width = chipWidth;
	length = sensLength + 2.0*guardRingWidth;
	thickness = sensThickness + chipThickness + bbThickness;
	TGeoVolume *container = mgr->MakeBox(Form("LAY%d_LADDER", layer),
                               medAir, 0.5*thickness, 0.5*width, 0.5*length);
	// the chip is a simple box:
	TGeoVolume *volChip = mgr->MakeBox("CHIP", medSPDSiChip,
                           0.5*chipThickness, 0.5*chipWidth, 0.5*chipLength);
	
	// the sensor is the union of a box and a border, to separate 
        // sensitive part from the rest the sensitive volume (inner part) 
        // is named according to the owner layer. To compute the shape 
	// subtraction which is needed for this we create two shapes,
	// which are two boxes with the same center. The smaller one is 
	// then used to define the sensor, while the subtraction of the two
	// is used for the guard ring.
	TGeoBBox  *shSens = new TGeoBBox(0.5*sensThickness, 0.5*sensWidth, 
					 0.5*sensLength);
	TGeoBBox  *shIn   = new TGeoBBox(sensThickness, 0.5*sensWidth, 
					 0.5*sensLength);
	TGeoBBox  *shOut  = new TGeoBBox(0.5*sensThickness, 
					 0.5*sensWidth + guardRingWidth, 
					 0.5*sensLength + guardRingWidth);
	shIn->SetName("innerBox");
	shOut->SetName("outerBox");
	TGeoCompositeShape *shBorder = new TGeoCompositeShape("",
							 "outerBox-innerBox");
	TGeoVolume *volSens = new TGeoVolume(Form("LAY%d_SENSOR", layer),
					     shSens, medSi);
	TGeoVolume *volBorder = new TGeoVolume("GUARD_RING", shBorder, medSi);

	// one line of bumpbonds
	TGeoVolume *volBB = mgr->MakeBox("BB", medBumpBond, 0.5*bbThickness,
					 0.5*bbWidth, 0.5*bbLength);
		
	// set colors of all objects for visualization	
	volSens->SetLineColor(kYellow + 1);
	volChip->SetLineColor(kGreen);
	volBorder->SetLineColor(kYellow + 3);

	// translations for the chip box: direction of length and 
	// thickness (moved down)
	TGeoTranslation *trChip[5] = {0, 0, 0, 0, 0};
	Double_t x = 0.5 * (chipThickness - thickness);
	Double_t y = 0.0;
	Double_t z = 0.0;
	Int_t i;
	for (i = 0; i < 5; i++) {
		z = -0.5*length + guardRingWidth + (Double_t)i*chipSpacing + 
		  ((Double_t)(i) + 0.5)*chipLength;
		trChip[i] = new TGeoTranslation(x, y, z);
	} // end for i
	
	// translation for the sensor parts: direction of width (moved 
	// to edge of container) and thickness (moved up)
	x = 0.5 * (thickness - sensThickness);
	y = 0.5 * (width - sensWidth - 2.0*guardRingWidth);
	z = 0.0;
	TGeoTranslation *trSens = new TGeoTranslation(x, y, z);
	
	// translation for the bump bonds:
	// keep same y used for sensors, but change the Z
	TGeoTranslation *trBB[160];
	//x = 0.5 * (thickness - bbThickness) + 0.5*sensThickness;
	x = 0.5 * (thickness - bbThickness) - sensThickness;
	z = -0.5 * sensLength + guardRingWidth + fgkmm*0.425 - bbPos;
	for (i = 0; i < 160; i++) {
		trBB[i] = new TGeoTranslation(x, y, z);
		switch(i) {
			case  31:
			case  63:
			case  95:
			case 127:
				z += fgkmm * 0.625 + fgkmm * 0.2;
				break;
			default:
				z += fgkmm * 0.425;
		} // end switch
	} // end for i
		
	// add nodes to container
	container->AddNode(volSens, 1, trSens);
	container->AddNode(volBorder, 1, trSens);
	for (i = 0; i < 160; i++) container->AddNode(volBB, i, trBB[i]);
	for (i = 0; i < 5; i++){
	    container->AddNode(volChip, i + 2, trChip[i]);
	} // end for i
	
	// return the container
	return container;
}
/*
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateGroundingFoilSingle(Bool_t kaptonLayer,
    Double_t &length, Double_t &width, Double_t &thickness, TGeoManager *mgr){
	//
	// Creates the grounding foil layer made in Kapton.
	// Both layers of the grounding foil have the same shape, but 
	// with small differences in the size of some parts (holes, 
	// overall size). The Kapton layer is a little bit wider and 
	// has smaller holes.
	// ---
	// The complete object is created as the superimposition of 
	// an XTRU with some holes
	// ---
	// Whenever possible, the size of the parts is parameterized with 
	// variable names, even if their value is fixed according 
	// to the design parameters given by engineers' drawings.
	// ---
	// Returns: a TGeoVolume object which contains all parts of this layer
	//
	
	// The shape of the grounding foil is an irregular polygon, which 
	// can easily be implemented as a TGeoXtru using the corners as 
	// reference points:
	// 
	// 0                                                                             1
	//  +-----------------------------------------------------------------------------+
	//  |                                    7              6      3                |
	//  |                                     +--------------+      +----------------+ 2
	//  |                         O           |              |      |
	//  |                             9 /-----+ 8            +------+ 4
	//  |                              /                    5
	//  |           11 /--------------/ 10
	//  +-------------/ 
	// 13           12
	//
	// in total: 14 points (X is just a referencem but is unused in 
	// the implementation. The whole shape can be subdivided into 
	// sectors delimited by vertical lines passing througth the 
	// points in the lower part of the shape. This convention is 
	// used to names their length which is different for each one 
	// (the widths, instead, are common for some)	
	// instantiate the media:
	// - kapton/aluminum for the pysical volumes
	TGeoMedium *material = kaptonLayer ? mgr->GetMedium("KAPTON") : 
	                                             mgr->GetMedium("AL");
	
	// label
	char type[3];
	if (kaptonLayer) {
		strcpy(type, "KP"); 
		thickness = fgkmm * 0.05;
	}
	else {
		strcpy(type, "AL");
		thickness = fgkmm * 0.02;
	}
	
	// define the length of all sectors (from leftmost to rightmost)
	Int_t i;
	Double_t sectorLength[] = {140.71,2.48,26.78,4.0,10.0,24.4,10.0,24.81};
	if (!kaptonLayer) {
		sectorLength[0] -= 0.2;
		sectorLength[4] -= 0.2;
		sectorLength[5] += 0.4;
		sectorLength[6] -= 0.4;
	}
	length = 0.0;
	for (i = 0; i < 8; i++) {
		sectorLength[i] *= fgkmm;
		length += sectorLength[i];
	}
		
	// as shown in the drawing, we have three different widths in 
	// this shape:
	Double_t widthMax  = fgkmm * 15.95;
	Double_t widthMed1 = fgkmm * 15.00;
	Double_t widthMed2 = fgkmm * 11.00;
	Double_t widthMin  = fgkmm *  4.40;
	if (!kaptonLayer) {
		widthMax  -= fgkmm * 0.4;
		widthMed1 -= fgkmm * 0.4;
		widthMed2 -= fgkmm * 0.4;
		widthMin  -= fgkmm * 0.4;
	}
	width = widthMax;
	
	// the vertices of the polygon are arrays correctly ordered in 
	// the counterclockwise direction: initially we place the point 
	// 0 in the origin, and all others will be defined accordingly
	Double_t x[14], y[14];
	x[ 0] = 0.0;
	y[ 0] = 0.0;
	
	x[ 1] = x[0] + length;
	y[ 1] = 0.0;
	
	x[ 2] = x[1];
	y[ 2] = -widthMin;
	
	x[ 3] = x[2] - sectorLength[7];
	y[ 3] = y[2];
	
	x[ 4] = x[3];
	y[ 4] = -widthMed2;
	
	x[ 5] = x[4] - sectorLength[6];
	y[ 5] = y[4];
	
	x[ 6] = x[5];
	y[ 6] = -widthMin;
	
	x[ 7] = x[6] - sectorLength[5];
	y[ 7] = y[6];
	
	x[ 8] = x[7];
	y[ 8] = -widthMed2;
	
	x[ 9] = x[8] - sectorLength[4];
	y[ 9] = y[8];
	
	x[10] = x[9] - sectorLength[3];
	y[10] = -widthMed1;
	 
	x[11] = x[10] - sectorLength[2];
	y[11] = y[10];
	
	x[12] = x[11] - sectorLength[1];
	y[12] = -widthMax;
	
	x[13] = x[0];
	y[13] = -widthMax;
	
	// then, we shift all points in such a way that the origin will 
	// be at the centers
	for (i = 0; i < 14; i++) {
		x[i] -= 0.5*length;
		y[i] += 0.5*width;
	}
	
	// create the shape
	char shName[200];
	sprintf(shName, "SH_%sGFOIL_FULL", type);
	TGeoXtru *shGroundFull = new TGeoXtru(2);
	shGroundFull->SetName(shName);
	shGroundFull->DefinePolygon(14, x, y);
	shGroundFull->DefineSection(0, -0.5*thickness, 0., 0., 1.0);
	shGroundFull->DefineSection(1,  0.5*thickness, 0., 0., 1.0);
	
	// this volume contains some holes which are here implemented 
	// as simple boxes of fixed size, which are displaced along the 
	// shape itself and then composed using the facilities of the 
	// TGeo package

	Double_t holeLength = fgkmm * 10.00;
	Double_t holeWidth  = fgkmm *  7.50;
	Double_t holeSepX0  = fgkmm *  7.05;//separation between center 
	                                    //of first hole and left border
	Double_t holeSepXC  = fgkmm * 14.00;//separation between the 
                                            // centers of two consecutive holes
	Double_t holeSepX1  = fgkmm * 15.42;//separation between centers 
                                            // of 5th and 6th hole
	Double_t holeSepX2  = fgkmm * 22.00;//separation between centers 
                                            // of 10th and 11th hole
	if (!kaptonLayer) {
		holeSepX0  -= fgkmm * 0.2;
		holeLength += fgkmm * 0.4;
		holeWidth  += fgkmm * 0.4;
	}
	
	// X position of hole center (will change for each hole)
	Double_t holeX = -0.5*length;
	// Y position of center of all holes (= 4.4 mm from upper border)
	Double_t holeY = 0.5*(width - holeWidth) - widthMin;
	//if (!kaptonLayer) holeY += 0.02;
		
	// create a shape for the holes (common)
	char holeName[200];
	sprintf(holeName, "%sHOLE", type);
	TGeoBBox *shHole = 0;
	shHole = new TGeoBBox(holeName,0.5*holeLength,0.5*holeWidth,thickness);
	
	// insert the holes in the XTRU shape:
	// starting from the first value of X, they are simply shifted 
	// along this axis
	char trName[200];
	TGeoTranslation *transHole[11];
	TString strComposite(shName);
	strComposite.Append("-(");
	for (Int_t i = 0; i < 11; i++) {
		// set the position of the hole, depending on index
		if (i == 0) {
			holeX += holeSepX0;
		}
		else if (i < 4) {
			holeX += holeSepXC;
		}
		else if (i == 4) {
			holeX += holeSepX1;
		}
		else if (i < 10) {
			holeX += holeSepXC;
		}
		else {
			holeX += holeSepX2;
		}
		sprintf(trName, "%sTR%d", type, i);
		transHole[i] = new TGeoTranslation(trName, holeX, holeY, 0.0);
		transHole[i]->RegisterYourself();
		strComposite.Append(holeName);
		strComposite.Append(":");
		strComposite.Append(trName);
		if (i < 10) strComposite.Append("+");
		//MM		cout << holeX << endl;
	}
	strComposite.Append(")");
	//MM	cout << strComposite.Data() << endl;
	
	// create composite shape (with holes)
	TGeoCompositeShape *shGround = new TGeoCompositeShape(
                                 Form("SH_%sGFOIL", type), strComposite.Data());
	
	// create the volume
	TGeoVolume *vol = new TGeoVolume(Form("%sGFOIL",type),shGround,
	                                                       material);
	return vol;
}
*/
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateGroundingFoilSingle(
      Bool_t kaptonLayer, Double_t &length, Double_t &width, 
      Double_t &thickness, TGeoManager *mgr){
	// Creates the grounding foil layer made in Kapton.
	// Both layers of the grounding foil have the same shape, but with small
	// differences in the size of some parts (holes, overall size).
	// The Kapton layer is a little bit wider and has smaller holes.
	// ---
	// The complete object is created as the sum of the following parts:
	// 1) the part which is connected to the chips, which is a 
	//    simple BOX with some box-shaped holes at regular intervals
	// 2) a trapezoidal connection where the Y size changes
	// 3) another box with a unique hole of the same shape and size as above
	// 4) another trapezoidal connection where the Y size changes
	// 5) a final part which is built as a sequence of 4 BOX volumes
	//    where the first and the third are equal and the others have 
        //    same size in Y.
	// ---
	// Whenever possible, the size of the parts is parameterized with 
	// variable names, even if their value is fixed according 
	// to the design parameters given by engineers' drawings.
	// ---
	// Returns: a TGeoVolume object which contanis all parts of this layer
	//

	// instantiate the media:
	// - vacuum for the container volume
	// - kapton for the pysical volumes
	TGeoMedium *vacuum   = mgr->GetMedium("VACUUM");
	TGeoMedium *material = mgr->GetMedium("KAPTON");
	
	// === Define size of all elements ===
	Double_t sizeZ      = fgkmm *   0.05;
	
	Double_t part1X     = fgkmm * 140.71;
	Double_t part2X     = fgkmm *   2.48;
	Double_t part3X     = fgkmm *  26.78;
	Double_t part4X     = fgkmm *   4.00;
	Double_t part5X     = fgkmm *  10.00;
	Double_t part6X     = fgkmm *  24.40;
	Double_t part7X     = fgkmm *  10.00;
	Double_t part8X     = fgkmm *  24.81;
	
	Double_t sizeYMax   = fgkmm *  15.95;
	Double_t sizeYMed1  = fgkmm *  15.00;
	Double_t sizeYMed2  = fgkmm *  11.00;
	Double_t sizeYMin   = fgkmm *   4.40;
	
	Double_t holeX      = fgkmm *  10.00;
	Double_t holeY      = fgkmm *   7.50;
	Double_t holeSepX   = fgkmm *  14.00;  // separation between the 
	                                       // centers of two consecutive 
	                                       // holes
	Double_t holeSepX1  = fgkmm *   1.42;  // to be added after 4th hole 
	                                       // in volume 1
	Double_t holeFirstX = fgkmm *   7.05;  // position of center of first 
	                                       // hole
	Double_t holeSepY   = fgkmm *   4.40;  // dist between hole's and 
	                                       // volume's upper border
	Double_t holeAloneX = fgkmm *  13.28;  // position of hole center 
	                                       // in box "part 3"

	// correct data in case we are on Aluminum foil
	if (!kaptonLayer) {
		material = mgr->GetMedium("AL");
		sizeZ       = fgkmm * 0.02;
		part1X     -= fgkmm * 0.2;
		part5X     -= fgkmm * 0.2;
		part6X     += fgkmm * 0.4;
		part7X     -= fgkmm * 0.4;
			
		sizeYMax   -= fgkmm * 0.4;
		sizeYMed1  -= fgkmm * 0.4;
		sizeYMed2  -= fgkmm * 0.4;
		sizeYMin   -= fgkmm * 0.4;
	
		holeX      += fgkmm * 0.4;
		holeY      += fgkmm * 0.4;
		holeFirstX -= fgkmm * 0.2;
		holeSepY   -= fgkmm * 0.4;
	}
	
	// define names for the object
	char type[4];
	if (kaptonLayer) strcpy(type, "KAP"); else strcpy(type, "ALU");
	
	// compute full length and width
	length = part1X + part2X + part3X + part4X + part5X + part6X + 
	         part7X + part8X;
	width = sizeYMax;
	thickness = sizeZ;
		
	// grounding foil world, bounded exactly around the limits 
	// of the structure
	TGeoVolume *container = mgr->MakeBox(Form("GFOIL_%s", type), 
			      vacuum, 0.5*length, 0.5*sizeYMax, 0.5*sizeZ);

	// === PART 1: box with holes ===

	TGeoBBox *shBox1 = 0, *shHole = 0;
	shBox1 = new TGeoBBox(Form("GF%s_BOX1", type), 0.5*part1X, 
			      0.5*sizeYMax, 0.5*sizeZ);
	shHole = new TGeoBBox(Form("GF%s_HOLE", type), 0.5*holeX, 0.5*holeY, 
			      0.5*sizeZ + 0.01);

	// define the position of all holes and compose the expression
	// to define the composite shape (box - holes)
	Double_t firstX = -0.5*part1X + holeFirstX;
	Double_t transY =  0.5*sizeYMax - holeSepY - 0.5*holeY;
	Double_t transX;
	TGeoTranslation *transHole[10];
	TString strComposite(Form("%s - (", shBox1->GetName()));
	for (Int_t i = 0; i < 10; i++) {
		transX = firstX + (Double_t)i * holeSepX;
		if (i > 4) transX += holeSepX1;
		transHole[i] = new TGeoTranslation(Form("TGF%s_HOLE%d",type,i),
						   transX, transY, 0.0);
		transHole[i]->RegisterYourself();
		strComposite.Append(Form("%s:%s", shHole->GetName(),
					 transHole[i]->GetName()));
		if (i < 9) strComposite.Append("+"); 
		else strComposite.Append(")");
	} // end for i
	// create composite shape
	TGeoCompositeShape *shPart1 = new TGeoCompositeShape(
                          Form("GF%s_PART1_SHAPE", type), strComposite.Data());
	// create the volume
	TGeoVolume *volPart1 = new TGeoVolume(Form("GF%s_PART1", type),
					      shPart1, material);

	// === PART 2: first trapezoidal connection
	
	TGeoArb8 *shTrap1 = new TGeoArb8(0.5*sizeZ);
	shTrap1->SetVertex(0, -0.5*part2X,  0.5*sizeYMax);
	shTrap1->SetVertex(1,  0.5*part2X,  0.5*sizeYMax);
	shTrap1->SetVertex(2,  0.5*part2X,  0.5*sizeYMax - sizeYMed1);
	shTrap1->SetVertex(3, -0.5*part2X, -0.5*sizeYMax);
	shTrap1->SetVertex(4, -0.5*part2X,  0.5*sizeYMax);
	shTrap1->SetVertex(5,  0.5*part2X,  0.5*sizeYMax);
	shTrap1->SetVertex(6,  0.5*part2X,  0.5*sizeYMax - sizeYMed1);
	shTrap1->SetVertex(7, -0.5*part2X, -0.5*sizeYMax);
	TGeoVolume *volPart2 = new TGeoVolume(Form("GF%s_PART2", type),
					      shTrap1, material);
	
	// === PART 3: other box with one hole
	
	TGeoBBox *shBox2 = 0;
	shBox2 = new TGeoBBox(Form("GF%s_BOX2", type), 0.5*part3X,
			      0.5*sizeYMed1, 0.5*sizeZ);
		
	// define the position of the hole
	transX = holeAloneX - 0.5*part3X;
	TGeoTranslation *transHoleAlone = new TGeoTranslation(
                        Form("TGF%s_HOLE_ALONE", type), transX, transY, 0.0);
	transHoleAlone->RegisterYourself();
	// create composite shape
	TGeoCompositeShape *shPart3 = new TGeoCompositeShape(
                        Form("GF%sPART3_SHAPE", type),
                        Form("%s - %s:%s", shBox2->GetName(),
			     shHole->GetName(), transHoleAlone->GetName()));
	// create the volume
	TGeoVolume *volPart3 = new TGeoVolume(Form("GF%s_PART3", type),
					      shPart3, material);
		
	// === PART 4: second trapezoidal connection
	
	TGeoArb8 *shTrap2 = new TGeoArb8(0.5*sizeZ);
	shTrap2->SetVertex(0, -0.5*part4X,  0.5*sizeYMed1);
	shTrap2->SetVertex(1,  0.5*part4X,  0.5*sizeYMed1);
	shTrap2->SetVertex(2,  0.5*part4X,  0.5*sizeYMed1 - sizeYMed2);
	shTrap2->SetVertex(3, -0.5*part4X, -0.5*sizeYMed1);
	shTrap2->SetVertex(4, -0.5*part4X,  0.5*sizeYMed1);
	shTrap2->SetVertex(5,  0.5*part4X,  0.5*sizeYMed1);
	shTrap2->SetVertex(6,  0.5*part4X,  0.5*sizeYMed1 - sizeYMed2);
	shTrap2->SetVertex(7, -0.5*part4X, -0.5*sizeYMed1);
	TGeoVolume *volPart4 = new TGeoVolume(Form("GF%s_PART4", type),
					      shTrap2, material);
		
	// === PART 5 --> 8: sequence of boxes ===
	
	TGeoVolume *volPart5 = mgr->MakeBox(Form("GF%s_BOX3", type),
			     material, 0.5*part5X, 0.5*sizeYMed2, 0.5*sizeZ);
	TGeoVolume *volPart6 = mgr->MakeBox(Form("GF%s_BOX4", type),
		             material, 0.5*part6X, 0.5*sizeYMin , 0.5*sizeZ);
	TGeoVolume *volPart7 = mgr->MakeBox(Form("GF%s_BOX5", type),
			     material, 0.5*part7X, 0.5*sizeYMed2, 0.5*sizeZ);
	TGeoVolume *volPart8 = mgr->MakeBox(Form("GF%s_BOX6", type), 
			      material, 0.5*part8X, 0.5*sizeYMin , 0.5*sizeZ);
	
	// === SET COLOR ===
	if (kaptonLayer) {
		volPart1->SetLineColor(kRed + 3);
		volPart2->SetLineColor(kRed + 3);
		volPart3->SetLineColor(kRed + 3);
		volPart4->SetLineColor(kRed + 3);
		volPart5->SetLineColor(kRed + 3);
		volPart6->SetLineColor(kRed + 3);
		volPart7->SetLineColor(kRed + 3);
		volPart8->SetLineColor(kRed + 3);
	}else{
		volPart1->SetLineColor(kGreen);
		volPart2->SetLineColor(kGreen);
		volPart3->SetLineColor(kGreen);
		volPart4->SetLineColor(kGreen);
		volPart5->SetLineColor(kGreen);
		volPart6->SetLineColor(kGreen);
		volPart7->SetLineColor(kGreen);
		volPart8->SetLineColor(kGreen);
	} // end if (kaptonLayer)
		
	// === TRANSLATION OF ALL PARTS ===
	
	transX = 0.5*(part1X - length);
	TGeoTranslation *transPart1 = new TGeoTranslation(transX, 0.0, 0.0);
	transX += 0.5*(part1X + part2X);
	TGeoTranslation *transPart2 = new TGeoTranslation(transX, 0.0, 0.0);
	transX += 0.5*(part2X + part3X);
	transY  = 0.5*(sizeYMax - sizeYMed1);
	TGeoTranslation *transPart3 = new TGeoTranslation(transX, transY, 0.0);
	transX += 0.5*(part3X + part4X);
	TGeoTranslation *transPart4 = new TGeoTranslation(transX, transY, 0.0);
	transX += 0.5*(part4X + part5X);
	transY  = 0.5*(sizeYMax - sizeYMed2);
	TGeoTranslation *transPart5 = new TGeoTranslation(transX, transY, 0.0);
	transX += 0.5*(part5X + part6X);
	transY  = 0.5*(sizeYMax - sizeYMin);
	TGeoTranslation *transPart6 = new TGeoTranslation(transX, transY, 0.0);
	transX += 0.5*(part6X + part7X);
	transY  = 0.5*(sizeYMax - sizeYMed2);
	TGeoTranslation *transPart7 = new TGeoTranslation(transX, transY, 0.0);
	transX += 0.5*(part7X + part8X);
	transY  = 0.5*(sizeYMax - sizeYMin);
	TGeoTranslation *transPart8 = new TGeoTranslation(transX, transY, 0.0);
	
	// add the partial volumes to the container
	container->AddNode(volPart1, 1, transPart1);
	container->AddNode(volPart2, 2, transPart2);
	container->AddNode(volPart3, 3, transPart3);
	container->AddNode(volPart4, 4, transPart4);
	container->AddNode(volPart5, 5, transPart5);
	container->AddNode(volPart6, 6, transPart6);
	container->AddNode(volPart7, 7, transPart7);
	container->AddNode(volPart8, 8, transPart8);
			
	return container;
}

//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateGroundingFoil(Double_t &thickness,
                                                      TGeoManager *mgr){
	// Joins two Kapton and two Aluminum layers of the grounding foil
	// in order to create the complete grounding foil for a whole stave.
	// into a unique container volume, which is returned as output.
	// The use of the TGeoXtru shape requires that in the separate 
        // foils, the Z axis lies perpendicularly to the polygonal basis 
        // of this shape;  this caused the components to have their Z 
        // axis corresponding to the X axis of the ALICE reference frame 
        // and vieceversa; to correct this, a rotation is necessary 
        // around their middle axis, to exchange X and Z axes and displace 
        // the object correctly in the ALICE frame.
	// ---
	// Arguments:
	//  - the sizes of the container box (passed by reference and 
        //    filled here)
	//  - the TGeoManager
	// ---
	// Returns: 
	//  - the container TGeoBBox (return value)
	//  - the size of the container (reference variables)
	//

	// sizes of the added volumes, which are filled by passing them 
	// to the volume creation methods
	Double_t kpLength, kpWidth, kpThick;
	Double_t alLength, alWidth, alThick;
	// separation between left and right volumes
	Double_t separation = fgkmm * 1.42;
	
	// create the two component volumes (each one will be replicated 
	// twice) this gives also the size of their virtual container 
	// boxes (just a reference, not a volume)
	TGeoVolume *kVol = CreateGroundingFoilSingle(kTRUE, kpLength, 
						     kpWidth, kpThick, mgr);
	TGeoVolume *aVol = CreateGroundingFoilSingle(kFALSE, alLength, 
						     alWidth, alThick, mgr);
	kVol->SetLineColor(kRed);
	aVol->SetLineColor(kGray);
	
	// kapton leads the total size of the foil (including spagcing 
	// of 1.42 mm between them in the center)
	Double_t length, width;
	length    = 2.0 * kpLength + separation;
	width     = kpWidth;
	thickness = kpThick + alThick;
	
	// create the container
	TGeoMedium *vacuum = mgr->GetMedium("VACUUM");
	TGeoVolume *container = mgr->MakeBox("GFOIL", vacuum, 0.5*thickness, 
					     0.5*width, 0.5*length);
	
	// create the common correction rotations
	TGeoRotation *rotCorr1 = new TGeoRotation(*gGeoIdentity);
	TGeoRotation *rotCorr2 = new TGeoRotation(*gGeoIdentity);
	rotCorr1->RotateY(-90.0);
	rotCorr2->RotateY( 90.0);
		
	// compute the translations to place the objects at the edges of 
	// the volume the kapton foils are also shifted down, and the 
	// aluminum foils are shifted up with respect to the thickness 
	// direction
	TGeoTranslation *kTrans1 = new TGeoTranslation(0.5*(-thickness+kpThick),
						       0.0,
						       0.5*(length-kpLength));
	TGeoTranslation *kTrans2 = new TGeoTranslation(0.5*(-thickness+kpThick),
						       0.0,
						       0.5*(-length+kpLength));
	TGeoTranslation *aTrans1 = new TGeoTranslation(0.5*(thickness-alThick),
						       0.0,
						    0.5*(length-alLength)-0.02);
	TGeoTranslation *aTrans2 = new TGeoTranslation(0.5*(thickness-alThick),
						       0.0,
						   0.5*(-length+alLength)+0.02);
	
	// combine translations and rotations
	TGeoCombiTrans *kCombi1 = new TGeoCombiTrans(*kTrans1, *rotCorr1);
	TGeoCombiTrans *kCombi2 = new TGeoCombiTrans(*kTrans2, *rotCorr2);
	TGeoCombiTrans *aCombi1 = new TGeoCombiTrans(*aTrans1, *rotCorr1);
	TGeoCombiTrans *aCombi2 = new TGeoCombiTrans(*aTrans2, *rotCorr2);
		
	// add to container
	container->AddNode(kVol, 0, kCombi1);
	container->AddNode(kVol, 1, kCombi2);
	container->AddNode(aVol, 0, aCombi1);
	container->AddNode(aVol, 1, aCombi2);
	
	return container;
}

//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateMCMBase(TGeoManager *geom) const{
	// Creates the MCM basis volume.
	// It is a little bit more complicated because this is a plain base
	// with a poly shape similar to the one of grounding foil but 
        // there are also some chips glued to its base and covered with 
        // a cave cap.
	// ---
	// The complete MCM object is created as the sum of the following parts:
	// 1) a planar basis shaped according to the MCM typical shape
	// 2) some boxes which represent the chips and devices mounted on 
        //    this base
	// 3) a cave cap which covers the portion of MCM containing these chips
	// ---
	// Due to the different widths of MCM, it is implemented in a 
        // more complicated way:
	// - cap and chips will define a sub-volume of this structure, 
        //   which can be bounded by a complete box
	// - base of MCM will be a separate volume
	// - these two objects will need to be glued together into an 
        //  upper-level volume
	// ---
	// This metod creates only the thin base (point 1 in the list)
	//
	
	// medium
	TGeoMedium *medBase = geom->GetMedium("MCM BASE");
	
	// parameterize the interesting sizes of MCM
	// it is divided into 3 sectors which have different size in 
	// X and Y and are connected by trapezoidal-based shapes, 
	// where the oblique angle makes a 45 degrees angle with the 
	// vertical, so that the X size and Y size of these "intermezzo"'s 
	// is the same
	// +--------------------------------+
	// |                   sect 2       |
	// | sect 1     --------------------+
	// +-----------/
	Double_t sizeZ = fgkmm * 0.35;
	Double_t sizeXtot = fgkmm * 105.6;
	Double_t sizeXsector[3] = {fgkmm * 28.4, fgkmm * 41.4, fgkmm * 28.8};
	Double_t sizeYsector[3] = {fgkmm * 15.0, fgkmm * 11.0, fgkmm *  8.0};
	Double_t sizeSep01 = fgkmm * 4.0, sizeSep12 = fgkmm * 3.0;
	Double_t sizeHole = fgkmm * 1.0;
	Double_t posHoleX = fgkmm * -0.5*sizeXtot + 26.7 + 0.5*sizeHole;
	Double_t posHoleY = fgkmm * -0.5*sizeYsector[0] + 0.5*sizeHole;
	
	// define the shape of base volume as an XTRU with two identical faces 
	// distantiated by the width of the  itself
	Double_t x[8], y[8];
	x[0] = -0.5*sizeXtot;
	y[0] =  0.5*sizeYsector[0];
	x[1] = -x[0];
	y[1] =  y[0];
	x[2] =  x[1];
	y[2] =  y[1] - sizeYsector[2];
	x[3] =  x[2] - sizeXsector[2];
	y[3] =  y[2];
	x[4] =  x[3] - sizeSep12;
	y[4] =  y[3] - sizeSep12;
	x[5] =  x[4] - sizeXsector[1];
	y[5] =  y[4];
	x[6] =  x[5] - sizeSep01;
	y[6] =  y[5] - sizeSep01;
	x[7] =  x[0];
	y[7] = -y[0];
	
	// create shape
	TGeoXtru *shPoly = new TGeoXtru(2);
	shPoly->SetName("SH_MCMBASE_POLY");
	shPoly->DefinePolygon(8, x, y);
	shPoly->DefineSection(0, -0.5*sizeZ, 0., 0., 1.0);
	shPoly->DefineSection(1,  0.5*sizeZ, 0., 0., 1.0);
	
	// create small hole
	TGeoBBox *shHole = 0;
	shHole = new TGeoBBox("SH_MCMBASE_HOLE", 0.5*sizeHole, 0.5*sizeHole, 
			      0.5*sizeZ+0.01);
	TGeoTranslation *transHole = new TGeoTranslation("TR_MCMBASE_HOLE", 
						    posHoleX, posHoleY, 0.0);
	transHole->RegisterYourself(); 
	
	// create shape intersection
	TGeoCompositeShape *shBase = new TGeoCompositeShape("SH_MCMBASE",
                          "SH_MCMBASE_POLY - SH_MCMBASE_HOLE:TR_MCMBASE_HOLE");
	
	// create volume
	TGeoVolume *volBase = new TGeoVolume("VOL_MCMBASE", shBase, medBase);
	volBase->SetLineColor(kRed);
	
	return volBase;
}
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateMCMCoverBorder(TGeoManager *geom){
	// Creates the MCM basis volume.
	// It is a little bit more complicated because this is a plain base
	// with a poly shape similar to the one of grounding foil but there 
        // are also some chips glued to its base and covered with a cave cap.
	// ---
	// The complete MCM object is created as the sum of the following parts:
	// 1) a planar basis shaped according to the MCM typical shape
	// 2) some boxes which represent the chips and devices mounted on 
        //  this base
	// 3) a cave cap which covers the portion of MCM containing these chips
	// ---
	// Due to the different widths of MCM, it is implemented in a more 
        // complicated way:
	// - cap and chips will define a sub-volume of this structure, 
        //   which can be bounded by a complete box
	// - base of MCM will be a separate volume
	// - these two objects will need to be glued together into an 
        //   upper-level volume
	// ---
	// This metod creates the thicker cap and its contents (points 2-3 
        // in the list). Since it covers only two of the three sectors of 
        // the MCM base with different width
	// the computations and variables related to the largest sector 
        // are removed, while
	// the other are the same as the other part of the MCM.
	//
	
	// media
	TGeoMedium *medCap  = geom->GetMedium("MCM COVER");
	
	// parameterize the interesting sizes of MCM
	// it is divided into 3 sectors which have different size in 
	// X and Y and  are connected by trapezoidal-based shapes, 
	// where the oblique angle makes a 45 degrees angle with the 
	// vertical, so that the X size and Y size of these "intermezzo"'s 
	// is the same
	// +--------------------------------+
	// |                   sect 2       |
	// | sect 1     --------------------+
	// +-----------/
	Double_t sizeZ = fgkmm * 0.3;
	Double_t capHeight = fgkmm * 1.7 - sizeZ;
	Double_t sizeXtot = fgkmm * 73.2;
	Double_t sizeXsector[2] = {fgkmm * 41.4, fgkmm * 28.8};
	Double_t sizeYsector[2] = {fgkmm * 11.0, fgkmm *  8.0};
	Double_t sizeSep = fgkmm * 3.0;
	
	// === PART 1: border ===
	
	// define the shape of base volume as an XTRU with two identical faces 
	// distantiated by the width of the  itself
	Double_t x[6], y[6];
	x[0] = -0.5*sizeXtot;
	y[0] =  0.5*sizeYsector[0];
	x[1] = -x[0];
	y[1] =  y[0];
	x[2] =  x[1];
	y[2] =  y[1] - sizeYsector[1];
	x[3] =  x[2] - sizeXsector[1];
	y[3] =  y[2];
	x[4] =  x[3] - sizeSep;
	y[4] =  y[3] - sizeSep;
	x[5] =  x[0];
	y[5] = -y[0];
	
	// create outer border shape with above coordinates
	TGeoXtru *capOut = new TGeoXtru(2);
	capOut->SetName("SH_MCMCAPOUT");
	capOut->DefinePolygon(6, x, y);
	capOut->DefineSection(0, -0.5*capHeight, 0., 0., 1.0);
	capOut->DefineSection(1,  0.5*capHeight, 0., 0., 1.0);
	
	// the inner border is built similarly but subtracting the thickness
	Double_t angle = 45.0;
	Double_t cs = TMath::Cos( 0.5*(TMath::Pi() - angle*TMath::DegToRad()) );
	Double_t xin[6], yin[6];
	xin[0] = x[0] + sizeZ;
	yin[0] = y[0] - sizeZ;
	xin[1] = x[1] - sizeZ;
	yin[1] = yin[0];
	xin[2] = xin[1];
	yin[2] = y[2] + sizeZ;
	xin[3] = x[3] - sizeZ*cs;
	yin[3] = yin[2];
	xin[4] = xin[3] - sizeSep;
	yin[4] = y[4] + sizeZ;
	xin[5] = xin[0];
	yin[5] = yin[4];
		
	// create inner border shape
	TGeoXtru *capIn = new TGeoXtru(2);
	capIn->SetName("SH_MCMCAPIN");
	capIn->DefinePolygon(6, xin, yin);
	capIn->DefineSection(0, -0.5*capHeight-0.01, 0., 0., 1.0);
	capIn->DefineSection(1,  0.5*capHeight+0.01, 0., 0., 1.0);
	
	// compose shape
	TGeoCompositeShape *shBorder = new TGeoCompositeShape("SH_MCMCAPBORDER",
                                "SH_MCMCAPOUT-SH_MCMCAPIN");
	
	// create volume
	TGeoVolume *volBorder = new TGeoVolume("VOL_MCMCAPBORDER", shBorder,
					       medCap);
	volBorder->SetLineColor(kGreen);
	
	return volBorder;
}
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateMCMCoverTop(TGeoManager *geom){
	// Creates the MCM basis volume.
	// It is a little bit more complicated because this is a plain base
	// with a poly shape similar to the one of grounding foil but 
        // there are also
	// some chips glued to its base and covered with a cave cap.
	// ---
	// The complete MCM object is created as the sum of the following parts:
	// 1) a planar basis shaped according to the MCM typical shape
	// 2) some boxes which represent the chips and devices mounted on 
        //    this base
	// 3) a cave cap which covers the portion of MCM containing these chips
	// ---
	// Due to the different widths of MCM, it is implemented in a 
        // more complicated way:
	// - cap and chips will define a sub-volume of this structure, 
        //   which can be bounded by a complete box
	// - base of MCM will be a separate volume
	// - these two objects will need to be glued together into an 
        //   upper-level volume
	// ---
	// This metod creates the thicker cap and its contents (points 
        // 2-3 in the list). Since it covers only two of the three 
        // sectors of the MCM base with different width
	// the computations and variables related to the largest sector 
        // are removed, while the other are the same as the other part 
        // of the MCM.
	//

	// media
	TGeoMedium *medCap  = geom->GetMedium("MCM COVER");
	
	// parameterize the interesting sizes of MCM
	// it is divided into 3 sectors which have different size in X 
	// and Y and  are connected by trapezoidal-based shapes, where 
	// the oblique angle makes a 45 degrees angle with the vertical, 
	// so that the X size and Y size of these "intermezzo"'s is the same
	// +--------------------------------+
	// |                   sect 2       |
	// | sect 1     --------------------+
	// +-----------/
	Double_t sizeZ = fgkmm * 0.3;
	Double_t sizeXtot = fgkmm * 73.2;
	Double_t sizeXsector[2] = {fgkmm * 41.4, fgkmm * 28.8};
	Double_t sizeYsector[2] = {fgkmm * 11.0, fgkmm *  8.0};
	Double_t sizeSep = fgkmm * 3.0;
	
	// === PART 1: border ===
	
	// define the shape of base volume as an XTRU with two identical faces 
	// distantiated by the width of the  itself
	Double_t x[6], y[6];
	x[0] = -0.5*sizeXtot;
	y[0] =  0.5*sizeYsector[0];
	x[1] = -x[0];
	y[1] =  y[0];
	x[2] =  x[1];
	y[2] =  y[1] - sizeYsector[1];
	x[3] =  x[2] - sizeXsector[1];
	y[3] =  y[2];
	x[4] =  x[3] - sizeSep;
	y[4] =  y[3] - sizeSep;
	x[5] =  x[0];
	y[5] = -y[0];
	
	// create outer border shape with above coordinates
	TGeoXtru *capOut = new TGeoXtru(2);
	capOut->SetName("SH_MCMCAPOUT");
	capOut->DefinePolygon(6, x, y);
	capOut->DefineSection(0, -0.5*sizeZ, 0., 0., 1.0);
	capOut->DefineSection(1,  0.5*sizeZ, 0., 0., 1.0);
	
	// the inner border is built similarly but subtracting the thickness
	Double_t angle = 45.0;
	Double_t cs = TMath::Cos( 0.5*(TMath::Pi() - angle*TMath::DegToRad()) );
	Double_t xin[6], yin[6];
	xin[0] = x[0] + sizeZ;
	yin[0] = y[0] - sizeZ;
	xin[1] = x[1] - sizeZ;
	yin[1] = yin[0];
	xin[2] = xin[1];
	yin[2] = y[2] + sizeZ;
	xin[3] = x[3] - sizeZ*cs;
	yin[3] = yin[2];
	xin[4] = xin[3] - sizeSep;
	yin[4] = y[4] + sizeZ;
	xin[5] = xin[0];
	yin[5] = yin[4];
		
	// coverage of upper part (equal to external border, but full)
	TGeoXtru *shCover = new TGeoXtru(2);
	shCover->SetName("SH_MCMCAPCOVER");
	shCover->DefinePolygon(6, x, y);
	shCover->DefineSection(0, -0.5*sizeZ, 0., 0., 1.0);
	shCover->DefineSection(1,  0.5*sizeZ, 0., 0., 1.0);
	
	// create volume
	TGeoVolume *volCover  = new TGeoVolume("VOL_MCMCAPCOVER", shCover, 
					       medCap);
	volCover->SetLineColor(kBlue);
	
	return volCover;
}
//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreateStave(Int_t layer, 
                             Double_t &fullThickness, TGeoManager *mgr){
	// Creates the complete stave as an assembly which contains all 
        // the stuff defined in the "CreateStaveBase" method (which are 
        // the thin part of the structure) and adds to this the thick 
        // cover of the MCM and the Pixel bus. This is done as an 
        // assembly to avoid the problem of a "ghost" overlap which occurs
	// when putting the stave on the carbon fiber sector, in the case 
        // that we define it as a volume container.
	// ---
	// Arguments:
	//     - the layer where the stave has to be put (hard check on this)
	//     - the geometry manager
	//
	
	// ** CRITICAL CHECK **
	// layer number can be ONLY 1 or 2
	if (layer != 1 && layer != 2) 
	  AliFatal("Required that layer number be 1 or 2");
	
	// sizes regarding the components
	Double_t baseWidth, baseHeight, baseThickness;
	Double_t mcmCapBorderThickness = fgkmm *  0.3;
	Double_t mcmCapThickness       = fgkmm *  1.7 - mcmCapBorderThickness;
	Double_t mcmCapHeight          = fgkmm * 11.0;
	Double_t mcmCapWidth           = fgkmm * 73.2;
	
	// create container
	TGeoVolumeAssembly *container = new TGeoVolumeAssembly(
                                         Form("LAY%d_FULLSTAVE", layer));
	
	// create subvolumes
	TGeoVolume *staveBase = CreateStaveBase(layer, baseWidth, baseHeight, 
						baseThickness, mgr);
	TGeoVolume *mcmCapBorder = CreateMCMCoverBorder(mgr);
	TGeoVolume *mcmCapTop = CreateMCMCoverTop(mgr);
	// bus in z > 0
	TGeoVolumeAssembly *bus0 = CreatePixelBusAndExtensions(kTRUE, mgr);
	// bus in z < 0
	TGeoVolumeAssembly *bus1 = CreatePixelBusAndExtensions(kFALSE, mgr);
	
	// the full width and height of the area which contains all 
	// components corresponds to the one of the stave base built with 
	// the "CreateStaveBase" method while the thickness must be 
	// computed as the sum of this base + the cover
	fullThickness = baseThickness + mcmCapThickness + mcmCapBorderThickness;
	
	// 1 - MCM cover	
		
	// translations (in the X direction, MCM is at the same level as ladder)
	Double_t xBase = -0.5*fullThickness + 0.5*baseThickness;
	TGeoTranslation *trBase = new TGeoTranslation(xBase, 0.0, 0.0);
	Double_t xMCMCapB = xBase + 0.5*baseThickness + 0.5*mcmCapThickness;
	Double_t xMCMCapT = xMCMCapB + 0.5*mcmCapThickness + 
	                                          0.5*mcmCapBorderThickness;
	Double_t yMCMCap  = 0.5*(baseHeight - mcmCapHeight);
	Double_t zMCMCap1 = 0.5*baseWidth - 0.5*mcmCapWidth;
	Double_t zMCMCap0 = -zMCMCap1;
	// correction rotations
	TGeoRotation *rotCorr0 = new TGeoRotation(*gGeoIdentity);
	TGeoRotation *rotCorr1 = new TGeoRotation(*gGeoIdentity);
	rotCorr0->RotateY( 90.0);
	rotCorr1->RotateY(-90.0);
	TGeoCombiTrans  *trMCMCapBorder0 = new TGeoCombiTrans(xMCMCapB, 
                                             yMCMCap, zMCMCap0, rotCorr0);
	TGeoCombiTrans  *trMCMCapBorder1 = new TGeoCombiTrans(xMCMCapB, 
                                             yMCMCap, zMCMCap1, rotCorr1);
	TGeoCombiTrans  *trMCMCapTop0 = new TGeoCombiTrans(xMCMCapT, 
                                             yMCMCap, zMCMCap0, rotCorr0);
	TGeoCombiTrans  *trMCMCapTop1 = new TGeoCombiTrans(xMCMCapT, 
                                             yMCMCap, zMCMCap1, rotCorr1);
	// add to container
	container->AddNode(staveBase, 0, trBase);
	container->AddNode(mcmCapBorder, 0, trMCMCapBorder0);
	container->AddNode(mcmCapBorder, 1, trMCMCapBorder1);
	container->AddNode(mcmCapTop, 0, trMCMCapTop0);
	container->AddNode(mcmCapTop, 1, trMCMCapTop1);
	
	// 2 - Pixel Bus
	
	// translations
	// for the moment, a correction amount of 0.04 is required to 
	// place correctly the object in X and another correction of 
	// 0.015 in Z
	Double_t busHeight  = fgkmm * 13.8;
	Double_t xPixelBus  = xBase + baseThickness + 0.04;
	Double_t yPixelBus1 = 0.5*baseHeight - 0.5*busHeight + 
	                               0.5*(baseHeight - busHeight);
	Double_t zPixelBus0 = -0.25*baseWidth + 0.015 - 0.03;
	//Double_t zPixelBus0 = -0.5*(0.5*baseWidth - 0.04);
	Double_t zPixelBus1 = -zPixelBus0;
	// correction rotations
	TGeoRotation *rotCorrBus1 = new TGeoRotation(*gGeoIdentity);
	rotCorrBus1->RotateX(180.0);
	//TGeoCombiTrans *trBus0 = new TGeoCombiTrans(xPixelBus, 0.0, 
	//                                     zPixelBus0, rotCorrBus);
	TGeoTranslation *trBus0 = new TGeoTranslation(xPixelBus, 0.0, 
                                                           zPixelBus0);
	//TGeoTranslation *trBus1 = new TGeoTranslation(xPixelBus, 0.0, 
	//                                                 zPixelBus1);
	TGeoCombiTrans *trBus1 = new TGeoCombiTrans(xPixelBus, yPixelBus1, 
						    zPixelBus1, rotCorrBus1);

	// add to container
	container->AddNode(bus0, 0, trBus0);
	container->AddNode(bus1, 1, trBus1);
	
	return container;
}
//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySPD::CreatePixelBusAndExtensions(
                                     Bool_t zpos, TGeoManager *mgr){
  // Creates an assembly which contains the pixel bus and its extension
  // and the extension of the MCM.
  // By: Renaud Vernet
  // NOTE: to be defined its material and its extension in the outside direction
  //
  
  // ====   constants   =====

  //get the media
  TGeoMedium   *medPixelBus    = mgr->GetMedium("PIXEL BUS") ;
  TGeoMedium   *medPBExtender  = mgr->GetMedium("PIXEL BUS EXTENDER") ;
  TGeoMedium   *medMCMExtender = mgr->GetMedium("MCM EXTENDER") ;

  //geometrical constants
  const Double_t kGroundingThickness    =   0.07  * fgkmm ;
  const Double_t kGrounding2pixelBusDz  =   0.625 * fgkmm ;
  const Double_t kPixelBusThickness     =   0.28  * fgkmm ;
  const Double_t kGroundingWidthX       = 170.501 * fgkmm ;
  const Double_t kPixelBusContactDx     =   1.099 * fgkmm ;
  const Double_t kPixelBusWidthY        =  13.8   * fgkmm ;
  //design=20 deg.
  const Double_t kPixelBusContactPhi    =  20.0   * TMath::DegToRad();
  //design=?? 70 deg. seems OK
  const Double_t kPbExtenderPsi         =  70.0   * TMath::DegToRad(); 
  const Double_t kPbExtenderWidthY      =  11.0   * fgkmm ;
  const Double_t kPbExtenderTopZ        =   2.72  * fgkmm ;
  const Double_t kMcmThickness          =   0.35  * fgkmm ;
  const Double_t kMcmExtenderThickness  =   0.20  * fgkmm ;
  const Double_t kDeltaMcmMcmextender   =   1.6   * fgkmm ;
  const Double_t kHalfStaveTotalLength  = 247.64  * fgkmm ;
  const Double_t kDeltaYOrigin          =  15.95/2.* fgkmm ;
  const Double_t kDeltaXOrigin          =   1.1    * fgkmm ;
  const Double_t kDeltaZOrigin          = kHalfStaveTotalLength / 2. ;

  const Double_t kGrounding2pixelBusDz2 = kGrounding2pixelBusDz+
                            kGroundingThickness/2. + kPixelBusThickness/2. ;
  const Double_t kPixelBusWidthX        = kGroundingWidthX ;
  const Double_t kPixelBusRaiseLength   = (kPixelBusContactDx-
	      kPixelBusThickness*TMath::Sin(kPixelBusContactPhi))/
                                      TMath::Cos(kPixelBusContactPhi) ;
  const Double_t kPbExtenderBaseZ       = kGrounding2pixelBusDz2 + 
              kPixelBusRaiseLength*TMath::Sin(kPixelBusContactPhi) + 
              2*kPixelBusThickness*TMath::Sin(kPixelBusContactPhi)*
              TMath::Tan(kPixelBusContactPhi) ;
  const Double_t kPbExtenderDeltaZ      = kPbExtenderTopZ-kPbExtenderBaseZ ;
  const Double_t kPbExtenderEndPointX   = 2*kDeltaZOrigin - kGroundingWidthX - 
                        2*kPixelBusThickness*TMath::Sin(kPixelBusContactPhi) ;
  const Double_t kMcmextenderEndPointX  = kDeltaZOrigin - 48.2 * fgkmm ;
  const Double_t kMcmExtenderWidthY     = kPbExtenderWidthY ;

  //=====  end constants  =====

  
  /*
  // -----------------   CREATE THE PIXEL BUS --------------------------
  // At the end of the pixel bus, a small piece is added for the contact 
  // with the pixel bus extender.
  // The whole piece is made with an extrusion, using 7 points
  //
  //                                   4
  //                                  /\
  // 6                            5  /  \ 3
  //  +-----------------------------+    /
  //  |                                 /
  //  +-----------------------------+--+
  // 0                              1   2
  //
  // The length of the pixel bus is defined (170.501mm) by the technical design
  // this length corresponds to distance [0-1] and [6-5]

  */

  TGeoVolumeAssembly *pixelBus = new TGeoVolumeAssembly("PIXEL BUS");

  // definition of the 7 points for the extrusion
  Double_t pixelBusXtruX[7] = {
    -kPixelBusWidthX/2. ,
    kPixelBusWidthX/2. ,
    kPixelBusWidthX/2. + kPixelBusThickness * TMath::Sin(kPixelBusContactPhi) ,
    kPixelBusWidthX/2. + kPixelBusThickness * TMath::Sin(kPixelBusContactPhi) +
                       kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) ,
    kPixelBusWidthX/2. + kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi),
    kPixelBusWidthX/2. ,
    -kPixelBusWidthX/2.
  } ;
  Double_t pixelBusXtruY[7] = {
    -kPixelBusThickness/2. ,
    -kPixelBusThickness/2. ,
    -kPixelBusThickness/2. + kPixelBusThickness *
                                 (1 - TMath::Cos(kPixelBusContactPhi)) ,
    -kPixelBusThickness/2. + kPixelBusThickness *
                                 (1 - TMath::Cos(kPixelBusContactPhi)) +
                  kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) ,
    kPixelBusThickness/2.  + kPixelBusRaiseLength * 
                                  TMath::Sin(kPixelBusContactPhi) ,
    kPixelBusThickness/2. ,
    kPixelBusThickness/2.
  } ;

  // creation of the volume
  TGeoXtru   *pixelBusXtru    = new TGeoXtru(2);
  TGeoVolume* pixelBusXtruVol = new TGeoVolume("pixelBusXtru",
                                                pixelBusXtru,medPixelBus) ;
  pixelBusXtru->DefinePolygon(7,pixelBusXtruX,pixelBusXtruY);
  pixelBusXtru->DefineSection(0,-kPixelBusWidthY/2.);
  pixelBusXtru->DefineSection(1, kPixelBusWidthY/2.);
  // --------------- END PIXEL BUS -------------------------------------


  // ------------------------- CREATE THE PIXEL BUS EXTENDER -----------
  // The geometry of the extender is a bit complicated sinceit is constrained
  // to be in contact with the pixel bus.
  // It consists of an extrusion using 13 points as shows the scheme below :
  //
  //                             8     7                       6
  //                               +---+---------------------+
  //                              /                          |
  //                             /                           |
  //                            /      +---------------------+
  //                           /      / 4                     5
  //                          /      /
  //       11  10          9 /      /
  //        +---+-----------+      /
  //       /                      /
  //      /                      /
  //     /      +-----------+---+
  // 12 +      / 1         2     3
  //     \    /
  //      \  /
  //        +
  //        0

  // ====   constants   =====
  const Double_t kPbExtenderXtru3L   = 1.5 * fgkmm ; //arbitrary ?
  const Double_t kPbExtenderXtru4L   = (kPbExtenderDeltaZ + 
                 kPixelBusThickness*(TMath::Cos(kPbExtenderPsi)-2))/
                     TMath::Sin(kPbExtenderPsi) ;
  //=====  end constants  =====

  TGeoVolumeAssembly *pbExtender = new TGeoVolumeAssembly("PIXEL BUS EXTENDER");

  Double_t pbExtenderXtruX[13] = {
    0, 
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) , 
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) + kPbExtenderXtru3L ,
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) + 
           kPbExtenderXtru3L + kPixelBusThickness * TMath::Sin(kPbExtenderPsi) ,
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) + 
           kPbExtenderXtru3L + kPixelBusThickness * 
           TMath::Sin(kPbExtenderPsi) + kPbExtenderXtru4L * 
           TMath::Cos(kPbExtenderPsi) ,
    kPbExtenderEndPointX ,
    kPbExtenderEndPointX ,
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) + 
           kPbExtenderXtru3L + kPixelBusThickness * TMath::Sin(kPbExtenderPsi)+
           kPbExtenderXtru4L * TMath::Cos(kPbExtenderPsi) ,
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi)  + 
           kPbExtenderXtru3L + kPixelBusThickness * TMath::Sin(kPbExtenderPsi)+
           kPbExtenderXtru4L * TMath::Cos(kPbExtenderPsi) - kPixelBusThickness*
           TMath::Sin(kPbExtenderPsi),
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) + kPbExtenderXtru3L ,
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) , 
    kPixelBusRaiseLength * TMath::Cos(kPixelBusContactPhi) - 
          kPixelBusThickness*TMath::Sin(kPixelBusContactPhi) , 
    -kPixelBusThickness * TMath::Sin(kPixelBusContactPhi)
  } ;
  Double_t pbExtenderXtruY[13] = {
    0, 
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) , 
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness * (1-TMath::Cos(kPbExtenderPsi)) ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness * (1-TMath::Cos(kPbExtenderPsi)) + 
                  kPbExtenderXtru4L * TMath::Sin(kPbExtenderPsi) ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness * (1-TMath::Cos(kPbExtenderPsi)) + 
                  kPbExtenderXtru4L * TMath::Sin(kPbExtenderPsi) ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness * (1-TMath::Cos(kPbExtenderPsi)) + 
                  kPbExtenderXtru4L * TMath::Sin(kPbExtenderPsi) + 
                  kPixelBusThickness ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness * (1-TMath::Cos(kPbExtenderPsi)) + 
                  kPbExtenderXtru4L * TMath::Sin(kPbExtenderPsi) + 
                  kPixelBusThickness ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi) + 
                  kPixelBusThickness + kPbExtenderXtru4L * 
                  TMath::Sin(kPbExtenderPsi),
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi)+kPixelBusThickness ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi)+kPixelBusThickness ,
    kPixelBusRaiseLength * TMath::Sin(kPixelBusContactPhi)+kPixelBusThickness*
                           TMath::Cos(kPixelBusContactPhi) ,
    kPixelBusThickness * TMath::Cos(kPixelBusContactPhi)
  } ;
  
  // creation of the volume
  TGeoXtru   *pbExtenderXtru    = new TGeoXtru(2);
  TGeoVolume *pbExtenderXtruVol = new TGeoVolume("pbExtenderXtru",
                                         pbExtenderXtru,medPBExtender) ;
  pbExtenderXtru->DefinePolygon(13,pbExtenderXtruX,pbExtenderXtruY);
  pbExtenderXtru->DefineSection(0,-kPbExtenderWidthY/2.);
  pbExtenderXtru->DefineSection(1, kPbExtenderWidthY/2.);
  // -------------- END PIXEL BUS EXTENDER -----------------------------


  // ------------------   CREATE THE MCM EXTENDER    -------------------
  // 
  // The MCM extender is located betwen the MCM and the Pixel Bus Extender
  // It consists of an extrusion using 10 points as shows the scheme below :
  //
  //                             7     6                       5
  //                               +---+---------------------+
  //                              /                          |
  //                             /                           |
  //                            /      +---------------------+
  //                           /      / 3                     4
  //                          /      /
  //            9          8 /      /
  //            +-----------+      /
  //            |                 /
  //            |                /
  //            +-----------+---+
  //            0          1     2

  //constants
  const Double_t kMcmExtenderXtru3L  = 1.5  * fgkmm ;
  //end constants

  TGeoVolumeAssembly *mcmExtender   = new TGeoVolumeAssembly("MCM EXTENDER");
  Double_t mcmExtenderXtruX[10] = {
    0 ,
    kMcmExtenderXtru3L ,
    kMcmExtenderXtru3L + kMcmExtenderThickness * TMath::Sin(kPbExtenderPsi) , 
    kMcmExtenderXtru3L + kMcmExtenderThickness * TMath::Sin(kPbExtenderPsi) + 
                           kDeltaMcmMcmextender / TMath::Tan(kPbExtenderPsi) ,
    kMcmextenderEndPointX ,
    kMcmextenderEndPointX ,
    kMcmExtenderXtru3L + kMcmExtenderThickness * TMath::Sin(kPbExtenderPsi) + 
                           kDeltaMcmMcmextender / TMath::Tan(kPbExtenderPsi) ,
    kMcmExtenderXtru3L + kDeltaMcmMcmextender / TMath::Tan(kPbExtenderPsi) ,
    kMcmExtenderXtru3L ,
    0
  } ;

  Double_t mcmExtenderXtruY[10] = {
    0 ,
    0 ,
    kMcmExtenderThickness*(1.-TMath::Cos(kPbExtenderPsi)),
    kMcmExtenderThickness*(1.-TMath::Cos(kPbExtenderPsi))+kDeltaMcmMcmextender,
    kMcmExtenderThickness*(1.-TMath::Cos(kPbExtenderPsi))+kDeltaMcmMcmextender,
    kMcmExtenderThickness*(2.-TMath::Cos(kPbExtenderPsi))+kDeltaMcmMcmextender,
    kMcmExtenderThickness*(2.-TMath::Cos(kPbExtenderPsi))+kDeltaMcmMcmextender,
    kMcmExtenderThickness + kDeltaMcmMcmextender ,
    kMcmExtenderThickness ,
    kMcmExtenderThickness ,
  } ;

  // creation of the volume
  TGeoXtru   *mcmExtenderXtru    = new TGeoXtru(2);
  TGeoVolume *mcmExtenderXtruVol = new TGeoVolume("mcmExtenderXtru",
					     mcmExtenderXtru,medMCMExtender) ;
  mcmExtenderXtru->DefinePolygon(10,mcmExtenderXtruX,mcmExtenderXtruY);
  mcmExtenderXtru->DefineSection(0,-kMcmExtenderWidthY/2.);
  mcmExtenderXtru->DefineSection(1, kMcmExtenderWidthY/2.);


  //--------------   DEFINITION OF GEOMETRICAL TRANSFORMATIONS ---------
  TGeoRotation    * commonRot       = new TGeoRotation("commonRot",0,90,0);
  commonRot->MultiplyBy(new TGeoRotation("rot",-90,0,0)) ;
  TGeoTranslation * pixelBusTrans   = new TGeoTranslation(kPixelBusThickness/2.
						 - kDeltaXOrigin + 0.52*fgkmm ,
					  -kPixelBusWidthY/2.+ kDeltaYOrigin , 
					  -kGroundingWidthX/2.+ kDeltaZOrigin) ;
  TGeoRotation    * pixelBusRot     = new TGeoRotation(*commonRot);
  TGeoTranslation * pbExtenderTrans = new TGeoTranslation(*pixelBusTrans) ;
  TGeoRotation    * pbExtenderRot   = new TGeoRotation(*pixelBusRot) ;
  pbExtenderTrans->SetDz(*(pbExtenderTrans->GetTranslation()+2)-
		   kPixelBusWidthX/2.-2.*kPixelBusThickness*
			 TMath::Sin(kPixelBusContactPhi)) ;  
  if (!zpos) {
    pbExtenderTrans->SetDy(*(pbExtenderTrans->GetTranslation()+1) - 
			   (kPixelBusWidthY - kPbExtenderWidthY)/2.);
  }else {
    pbExtenderTrans->SetDy(*(pbExtenderTrans->GetTranslation()+1) + 
			   (kPixelBusWidthY - kPbExtenderWidthY)/2.);
  } // end if !zpos
  pbExtenderTrans->SetDx(*(pbExtenderTrans->GetTranslation()) +
			 kPixelBusThickness/2 + 2*kPixelBusThickness*
			 TMath::Sin(kPixelBusContactPhi)*
			 TMath::Tan(kPixelBusContactPhi)) ;
  TGeoTranslation * mcmExtenderTrans = new TGeoTranslation(0.12*fgkmm + 
					kMcmThickness - kDeltaXOrigin,
				   pbExtenderTrans->GetTranslation()[1],-4.82);
  TGeoRotation    * mcmExtenderRot   = new TGeoRotation(*pbExtenderRot);

  //ADD NODES TO ASSEMBLIES
  pixelBus    ->AddNode((TGeoVolume*)pixelBusXtruVol,0);
  pbExtender  ->AddNode((TGeoVolume*)pbExtenderXtruVol,0);
  mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtruVol,0);
//   mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtru3Vol,0);
//   mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtru3PrimVol,1);
//   mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtru4Vol,2);
//   mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtru4PrimVol,3);
//   mcmExtender ->AddNode((TGeoVolume*)mcmExtenderXtru5Vol,4);


  //CREATE FINAL VOLUME ASSEMBLY AND ROTATE IT
  TGeoVolumeAssembly *assembly = new TGeoVolumeAssembly("EXTENDERS");
  assembly->AddNode((TGeoVolume*)pixelBus    ,0, 
		    new TGeoCombiTrans(*pixelBusTrans,*pixelBusRot));
  assembly->AddNode((TGeoVolume*)pbExtender  ,0, 
		    new TGeoCombiTrans(*pbExtenderTrans,*pbExtenderRot));
  assembly->AddNode((TGeoVolume*)mcmExtender ,0, 
		    new TGeoCombiTrans(*mcmExtenderTrans,*mcmExtenderRot));
  assembly->SetTransparency(50);
  return assembly ;
}
//______________________________________________________________________
TGeoVolume* AliITSv11GeometrySPD::CreateStaveBase(Int_t layer,
	   Double_t &fullWidth, Double_t &fullHeight, Double_t &fullThickness,
	   TGeoManager *mgr){
	// Creates a box which contains the followin parts of the whole stave:
	// - the two layers of grounding foil
	// - the ladders
	// - the thin base of the MCM (except its thick cover)
	// - the pixel bus
	// ---
	// Since it is required by detector numbering conventions, 
	// it is required as argument the layer which owns this stave.
	// This number will be used to define the name of the ladder volume, 
	// which must be different for layer1 and layer2 objects.
	// ---
	// Arguments:
	//    - layer number (will be checked to be 1 or 2)
	//    - geometry manager
	// ---
	// Returns:
	//    - a TGeoBBox volume containing all this stuff
	//    - the size of the container box are stored in the 
        //      reference-passed variables
	//
	// sizes of all objects to be inserted
	// these values are used to compute the total volume of the container
	// and to compute parametrically the position of each piece, instead
	// of putting hard-coded number (this helps in eventually modifying 
        // everything)
	Double_t mcmThickness    = fgkmm * 0.35;
	Double_t grndThickness   = fgkmm * 0.07; // = 0.05 + 0.02
	Double_t sepThickness    = fgkmm * 0.05;

	Double_t ladderWidth     = fgkmm *  70.72;
	Double_t mcmWidth        = fgkmm * 105.60;
	Double_t sepLaddersWidth = fgkmm *   0.20;
	Double_t sepMCMWidth     = fgkmm *   0.30;
	// separations between central ladders in the two half-staves 
	Double_t sepLaddersCtr   = fgkmm *   0.40; 

	Double_t mcmHeight       = fgkmm *  15.00;
	
	// compute the size of the container
	fullWidth     = 2.0*sepLaddersCtr + 4.0*ladderWidth + 
	            2.0*sepMCMWidth + 2.0*sepLaddersWidth + 2.0*mcmWidth;
	fullHeight    = fgkmm * 15.95;
	fullThickness = grndThickness + sepThickness + mcmThickness;
	
	// create the container
	TGeoVolume *container = mgr->MakeBox(Form("LAY%d_STAVE", layer),
			  mgr->GetMedium("VACUUM"), 0.5*fullThickness,
			  0.5*fullHeight, 0.5*fullWidth);
		
	// fill the container going from bottom to top 
	// with respect to the thickness direction
	
	// 1 - Grounding foil
	// volume
	TGeoVolume *grndVol = CreateGroundingFoil(grndThickness);
	// translation
	Double_t xGrnd = -0.5*fullThickness + 0.5*grndThickness;
	TGeoTranslation *grndTrans = new TGeoTranslation(xGrnd, 0.0, 0.0);
	// add to container
	container->AddNode(grndVol, 1, grndTrans);
	
	// 2 - Ladders
	// volume (will be replicated 4 times)
	Double_t ladderLength, ladderThickness;
	TGeoVolume *ladder = CreateLadder(layer, ladderLength, ladderWidth, 
					  ladderThickness, mgr);
	// translations (in thickness direction, the MCM thickness is used)
	// layers are sorted going from the one at largest Z to the one 
	// at smallest Z:
	// -|Zmax| ------> |Zmax|
	//      0   1   2   3
	// but it is more comfortable to start defining their Z position 
	// from center
	Double_t xLad  = xGrnd + 0.5*grndThickness + 0.5*mcmThickness + 
	                 sepThickness;
	Double_t zLad1 = -0.5*ladderWidth - sepLaddersCtr;
	Double_t zLad0 = zLad1 - ladderWidth - sepLaddersWidth;
	Double_t zLad2 = -zLad1;
	Double_t zLad3 = -zLad0;
	// rotLad->RotateZ(180.0);
	TGeoRotation   *rotLad = new TGeoRotation(*gGeoIdentity);
	TGeoCombiTrans *trLad0 = new TGeoCombiTrans(xLad, 0.0, zLad0, rotLad);
	TGeoCombiTrans *trLad1 = new TGeoCombiTrans(xLad, 0.0, zLad1, rotLad);
	TGeoCombiTrans *trLad2 = new TGeoCombiTrans(xLad, 0.0, zLad2, rotLad);
	TGeoCombiTrans *trLad3 = new TGeoCombiTrans(xLad, 0.0, zLad3, rotLad);
	// add to container
	container->AddNode(ladder, 0, trLad0);
	container->AddNode(ladder, 1, trLad1);
	container->AddNode(ladder, 2, trLad2);
	container->AddNode(ladder, 3, trLad3);
	
	// 3 - MCM (only the base, the cover is added as a separate 
	// volume in a more global 'stave' assembly volume (will be 
	// replicated twice)
	TGeoVolume *mcm = CreateMCMBase(mgr);
	// translations (in the X direction, MCM is at the same 
	// level as ladder) the two copies of the MCM are placed at 
	// the same distance from the center, on both sides and their 
	// sorting is the same as ladders' one (MCM0 is at Z < 0, 
	// MCM1 at Z > 0);
	Double_t xMCM  = xLad;
	Double_t yMCM  = 0.5*(fullHeight - mcmHeight);
	Double_t zMCM1 = zLad3 + 0.5*ladderWidth + 0.5*mcmWidth + sepMCMWidth;
	Double_t zMCM0 = -zMCM1;
	// create the common correction rotations
	TGeoRotation *rotCorr0 = new TGeoRotation(*gGeoIdentity);
	TGeoRotation *rotCorr1 = new TGeoRotation(*gGeoIdentity);
	rotCorr0->RotateY( 90.0);
	rotCorr1->RotateY(-90.0);
	TGeoCombiTrans *trMCM0 = new TGeoCombiTrans(xMCM,yMCM,zMCM0,rotCorr0);
	TGeoCombiTrans *trMCM1 = new TGeoCombiTrans(xMCM,yMCM,zMCM1,rotCorr1);
	// add to container
	container->AddNode(mcm, 0, trMCM0);
	container->AddNode(mcm, 1, trMCM1);
		
	return container;
}
//______________________________________________________________________
void AliITSv11GeometrySPD::StavesInSector(TGeoVolume *moth, TGeoManager *mgr){
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
	// Arguments: see description of "CarbonFiberSector" method.
	//
	
	// This service class is useful to this method only
	// to store in a meaningful way the data about the 
	// rounded corners of the support, and some computations
	// which could turn out to be useful for stave placement
	// 'left' and 'right' (L/R) here are considered looking the support
	// from the positive Z side.
	// The sign of the radius is used to know what kind of tangent
	// must be found for the two circles which describe the rounded angles.
	class clsSupportPlane {
	public:
                // curvature center and radius (with sign) of left corner
		Double_t xL, yL, rL, sL;
                // curvature center and radius (with sign) of right corner
		Double_t xR, yR, rR, sR;
                // shift from the innermost position (where the stave edge is
	        // in the point where the rounded corner begins
		Double_t shift;
		
		// Constructor with arguments which allow to set 
	        // directly everything since the values are given in 
	        // millimiters from drawings, they must be converted to cm
		clsSupportPlane
		(Double_t xLin, Double_t yLin, Double_t rLin, Double_t sLin, 
		 Double_t xRin, Double_t yRin, Double_t rRin, Double_t sRin,
                 Double_t shiftin) :
		 xL(xLin), yL(yLin), rL(rLin), sL(sLin), xR(xRin), yR(yRin), 
                 rR(rRin), sR(sRin), shift(shiftin) {
			xL *= fgkmm;
			yL *= fgkmm;
			rL *= fgkmm;
			xR *= fgkmm;
			yR *= fgkmm;
			rR *= fgkmm;
			shift *= fgkmm;
		} // end group.
		
		// Computation of the line tangent to both circles 
	        // defined here which is taken above or below the center 
	        // according to the radius sign. This method returns:
		//   - the mid-popint of the segment between the two 
	        //     points where the tangent touches the two circles, 
		//   - the inclination of this segment
		//   - the half-length of this segment
		Double_t TangentSegment(Double_t &midX, Double_t &midY, 
					Double_t &phi){
			// compute the straight line which is tangent to 
		        // the two circles and extract its inclination 
		        // 'phi' w.r. to X axis
			Double_t dx = xL - xR;
			Double_t dy = yL - yR;
			Double_t R  = rL*sL + rR*sR;
			Double_t delta = dy*dy + dx*dx - R*R;
			Double_t tan05phi = (-dy+TMath::Sqrt(delta))/(R - dx);
			phi = 2.0 * TMath::ATan(tan05phi);
			// compute the points where this line touchs the 
			// two circles
			Double_t leftX  = xL + sL*rL*TMath::Cos(phi);
			Double_t leftY  = yL + sL*rL*TMath::Sin(phi);
			Double_t rightX = xR + sR*rR*TMath::Cos(phi);
			Double_t rightY = yR + sR*rR*TMath::Sin(phi);
			// compute the mid point
			midX = 0.5 * (leftX + rightX);
			midY = 0.5 * (leftY + rightY);
			// compute angular coefficient for the line joining
			// the two points found using the above method
			dx = rightX - leftX;
			dy = rightY - leftY;
			phi = TMath::ATan2(dy, dx);
			// compute the half-length of this segment
			Double_t len = 0.5*TMath::Sqrt((rightX-leftX)*
                                            (rightX-leftX) + (rightY-leftY)*
                                                             (rightY-leftY));
			//MM			cout << 2.0*len << endl;
			return len;
		} // end function
	}; // end class	
	// instantiate this class for each layer1 and layer2 corners
	clsSupportPlane *plane[6] = {0, 0, 0, 0, 0, 0};
	
	// layer 2
	plane[0] = new clsSupportPlane( 10.830,  16.858, 0.60,  1.,  19.544,  
					10.961, 0.8,  1.,  1.816);
	plane[1] = new clsSupportPlane(- 0.733,  17.486, 0.60,  1.,  11.581,  
				       13.371, 0.6, -1., -0.610);
	plane[2] = new clsSupportPlane(-12.252,  16.298, 0.60,  1.,   0.562,  
				       14.107, 0.6, -1., -0.610);
	plane[3] = new clsSupportPlane(-22.276,  12.948, 0.85,  1., -10.445,  
				       13.162, 0.6, -1., -0.610);
	// layer 1
	plane[4] = new clsSupportPlane(- 3.123, -14.618, 0.50,  1.,  11.280, 
				       -14.473, 0.9, -1., -0.691);
	plane[5] = new clsSupportPlane(-13.187, -19.964, 0.50, -1., - 3.833, 
				       -17.805, 0.6, -1.,  1.300);
	// put the sector in the container
	//CarbonFiberSector(moth, xAAtubeCenter0, yAAtubeCenter0, mgr);
	
	// create stave volume
	Double_t staveHeight = 1.595, staveThickness;
	TGeoVolume *stave1 = CreateStave(1, staveThickness,mgr);
	TGeoVolume *stave2 = CreateStave(2, staveThickness,mgr);
		
	// compute positions and rotation angles
	Double_t xm, ym, halfPlaneHeight, heightDiff, position, phi, xPos, yPos;
	for (Int_t i = 0; i < 6; i++) {
	    //
	    // This functioninserted here for test. Added By Bjorn Nilsen
	    // August 29 2007.
	    Double_t x0,y0,x1,y1; // should be move out of loop
	    Bool_t lreturn;
	    lreturn = GetSectorMountingPoints(i,x0,y0,x1,y1,mgr);
	    //
		// recall the geometry computations defined for the classe
		halfPlaneHeight = plane[i]->TangentSegment(xm, ym, phi);
		// compute the difference between plane and stave heights
		heightDiff = halfPlaneHeight - 0.5*staveHeight;
		// It is necessary to shift the stave by at least 
		// an amount equal to this difference
		// to avoid overlaps.
		// Moreover, some more shift is done for building reasons,
		// and it depends on the single plane (data-member 'shift')
		position = heightDiff + plane[i]->shift;
		// taking into account this shift plus another in the direction
		// normal to the support plane, due to the stave thickness,
		// the final position of the stave is computed in a temporary 
		// reference frame where the mid-point of the support plane 
		// is in the origin
		if (i < 4) {
			ParallelPosition(0.5*staveThickness, position, phi, 
					 xPos, yPos);
		}else if (i == 4) {
			ParallelPosition(-0.5*staveThickness, -position, phi, 
					 xPos, yPos);
		}else {
			ParallelPosition(-0.5*staveThickness, -position, phi, 
					 xPos, yPos);
		}
		// then we go into the true reference frame
		xPos += xm;
		yPos += ym;
		/*
		// TEMP
		TGeoVolume *tubeTemp1 = mgr->MakeTube("tubeTemp1", NULL, 
		                                          0.0, 0.01, 50.0);
		TGeoTranslation *trTemp1 = new TGeoTranslation(xm, ym, 0.0);
		tubeTemp1->SetLineColor(kRed);
		moth->AddNode(tubeTemp1, i + 1, trTemp1);
		TGeoVolume *tubeTemp2 = mgr->MakeTube("tubeTemp2", NULL, 
                                                            0.0, 0.01, 50.0);
		TGeoTranslation *trTemp2 = new TGeoTranslation(xPos, yPos, 0.0);
		tubeTemp2->SetLineColor(kBlue);
		moth->AddNode(tubeTemp2, i + 1, trTemp2);
		// END TEMP
		*/
		// using the parameters found here, compute the 
		// translation and rotation of this stave:
		TGeoRotation *rot = new TGeoRotation(*gGeoIdentity);
		if (i >= 4) rot->RotateY(180.0);
		rot->RotateZ(90.0 + phi * TMath::RadToDeg());
		TGeoCombiTrans *trans = new TGeoCombiTrans(xPos,yPos,0.0,rot);
		if (i < 4) {
			moth->AddNode(stave2, i, trans);
		}else {
			moth->AddNode(stave1, i - 4, trans);
		} // end if i<4
	} // for i
}
//______________________________________________________________________
void AliITSv11GeometrySPD::ParallelPosition(Double_t dist1, Double_t dist2, 
				     Double_t phi, Double_t &x, Double_t &y){
	// Performs the following steps:
	// 1 - finds a straight line parallel to the one passing through 
        //     the origin and with angle 'phi' with X axis (phi in RADIANS);
	// 2 - finds another line parallel to the previous one, with a 
        //     distance 'dist1' from it
	// 3 - takes a reference point in the second line in the 
        //     intersection between the normal to both lines passing 
        //     through the origin
	// 4 - finds a point whith has distance 'dist2' from this 
        //     reference, in the second line (point 2)
	// ----
	// According to the signs given to dist1 and dist2, the point 
        // is found in different position w.r. to the origin

	// compute the point
	Double_t cs = TMath::Cos(phi);
	Double_t sn = TMath::Sin(phi);
	
	x = dist2*cs - dist1*sn;
	y = dist1*cs + dist2*sn;
}
//----------------------------------------------------------------------
Bool_t AliITSv11GeometrySPD::Make2DcrossSections(TPolyLine &a0,TPolyLine &a1,
			   TPolyLine &b0,TPolyLine &b1,TPolyMarker &p)const{
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
    Int_t n0,n1,i;
    Double_t x,y;
    TGeoVolume *a0V,*a1V,*b0V,*b1V;
    TGeoXtru *a0S,*a1S,*b0S,*b1S;
    TGeoManager *mgr = gGeoManager;

    a0V = mgr->GetVolume(fSPDsectorShapeName.Data());
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
    p.SetPolyMarker(2*fSPDsectorPoints0.GetSize());
    for(i=0;i<fSPDsectorPoints0.GetSize();i++){
      GetSectorMountingPoints(i,x0,y0,x1,y1);
      p.SetPoint(2*i,x0,y0);
      p.SetPoint(2*i+1,x1,y1);
    } // end for i
    return kTRUE;
}

