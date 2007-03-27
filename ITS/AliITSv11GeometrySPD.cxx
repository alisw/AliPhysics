/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// This class Defines the Geometry for the ITS services and support cones
// outside of the ceneteral volume (except for the Ceneteral support 
// cylinders. Other classes define the rest of the ITS. Specificaly the ITS
// The SSD support cone,SSD Support centeral cylinder, SDD support cone,
// The SDD cupport centeral cylinder, the SPD Thermal Sheald, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC. 

// General Root includes
#include <Riostream.h>
#include <TMath.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPolyLine.h>
// Root Geometry includes
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoEltu.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include "AliMagF.h"
#include "AliRun.h"
//#include <TGeoRotation.h>
//#include <TGeoCombiTrans.h>
//#include <TGeoTranslation.h>
#include "AliITSv11GeometrySPD.h"

ClassImp(AliITSv11GeometrySPD)

#define SQ(A) (A)*(A)

//______________________________________________________________________
Int_t AliITSv11GeometrySPD::CreateSPDCenteralMaterials(Int_t &medOffset,
						       Int_t &matOffset){
    // Define the specific materials used for the ITS SPD centeral
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
    const Double_t kdeemaxAir = 0.1; // Fraction of particle's energy 0<deemax<=1
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
    mix = new TGeoMixture("C (M55J)",4,1.9866*fgkgcm3);// Carbon fiber by fractional weight "C (M55J)"
    mix->SetIndex(matindex);
    mix->DefineElement(0,12.0107,6.0,0.908508078); // Carbon by fractional weight
    mix->DefineElement(1,14.0067,7.0,0.010387573); // Nitrogen by fractional weight
    mix->DefineElement(2,15.9994,8.0,0.055957585); // Oxigen by fractional weight
    mix->DefineElement(3,1.00794,1.0,0.025146765); // Hydrogen by fractional weight
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
    mix = new TGeoMixture("Air",4,1.20479E-3*fgkgcm3);// Carbon fiber by fractional weight
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
    mix = new TGeoMixture("INOX",9,8.03*fgkgcm3);// Carbon fiber by fractional weight
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
    //med = new TGeoMedium("ITSspdStainlessSteel",medindex++,matindex++,0,ifield,
    //	       fieldm,ktmaxfdAir,kstemaxAir,kdeemaxAir,kepsilAir,kstminAir);
    //
    mix = new TGeoMixture("Freon",2,1.63*fgkgcm3);// Carbon fiber by fractional weight
    mix->SetIndex(matindex);
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
void AliITSv11GeometrySPD::InitSPDCenteral(Int_t offset,TVirtualMC *vmc){
    // Do any SPD Centeral detector related initilizations, setting
    // transport cuts for example.
    // Some GEANT3 Physics switches
    // "MULTS"
    // Multiple scattering. The variable IMULS controls this process. For 
    // more information see [PHYS320 or 325 or 328].
    // 0 - No multiple scattering.
    // 1 - Multiple scattering according to Molière theory. Default setting.
    // 2 - Same as 1. Kept for backward compatibility.
    // 3 - Pure Gaussian scattering according to the Rossi formula.
    // "DRAY"
    // delta ray production. The variable IDRAY controls this process. See [PHYS430]
    // 0 - No delta rays production.
    // 1 - delta rays production with generation of . Default setting.
    // 2 - delta rays production without generation of .
    // "LOSS"
    // Continuous energy loss. The variable ILOSS controls this process.
    // 0 - No continuous energy loss, IDRAY is set to 0.
    // 1 - Continuous energy loss with generation of delta rays above 
    //     DCUTE (common/GCUTS/) and restricted Landau fluctuations below  DCUTE.
    // 2 - Continuous energy loss without generation of delta rays and full 
    //     Landau-Vavilov-Gauss fluctuations. In this case the variable IDRAY 
    //     is forced to 0 to avoid double counting of fluctuations. Default setting.
    // 3 - Same as 1, kept for backward compatibility.
    // 4 - Energy loss without fluctuation. The value obtained from the tables is 
    //     used directly.
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
    vCarbonFiberSector->SetVisibility(kFALSE); // logical volume
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
        printf("i=%d angle=%g angle[rad]=%g radiusSector=%g x=%g y=%g \n",
               i,angle,angle/fgkRadian,radiusSector,
               -radiusSector*TMath::Sin(angle/fgkRadian),
               radiusSector*TMath::Cos(angle/fgkRadian));
        angle += kSectorRelativeAngle;
        secRot->RotateZ(kSectorRelativeAngle);
    } // end for i
    if(GetDebug()){
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
    const Double_t ksecY9   = +14.486*fgkmm;
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
    TGeoTube *sTB0,*sTB1,*sM0;
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
    // Spcial cases
    secAngleStart2[8] -= 360.;
    secAngleStart2[11] -= 360.;
    //
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
    sA0->SetName("ITS SPD Carbon fiber support Sector A0");
    sA0->DefinePolygon(m,xpp,ypp);
    sA0->DefineSection(0,-ksecDz);
    sA0->DefineSection(1,ksecDz);
    //
    InsidePoint(xpp[m-1],ypp[m-1],xpp[0],ypp[0],xpp[1],ypp[1],
                ksecCthick,xpp2[0],ypp2[0]);
    for(i=1;i<m-1;i++){
        j = i/(ksecNPointsPerRadii+1);
        InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],
                    ksecCthick,xpp2[i],ypp2[i]);
    } // end for i
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
        InsidePoint(xpp[i-1],ypp[i-1],xpp[i],ypp[i],xpp[i+1],ypp[i+1],
                    t,xpp2[i],ypp2[i]);
    } // end for i
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
    sM0 = new TGeoTube("ITS SPD Sensitive Virutual Volume M0",0.0,8.0,
                       sA0->GetZ(1)+sB0->GetZ(1));
    //
    if(GetDebug()){
        cout<<"medSPDcf= "<<medSPDcf<<endl;
	//        printf("medSPDcf=%x\n",medSPDcf);
        if(medSPDcf) medSPDcf->Dump();
        cout<<"medSPDss= "<<medSPDss<<endl;
	//       printf("medSPDss=%x\n",medSPDss);
        if(medSPDss) medSPDss->Dump();
	cout<<"medSPDair= "<<medSPDair<<endl;
	//       printf("medSPDair=%x\n",medSPDair);
        if(medSPDair) medSPDair->Dump();
	cout<<"medSPDcoolfl= "<<medSPDcoolfl<<endl;
	//       printf("medSPDcoolfl=%x\n",medSPDcoolfl);
        if(medSPDcoolfl) medSPDcoolfl->Dump();
        sM0->InspectShape();
        sA0->InspectShape();
        sA1->InspectShape();
        sB0->InspectShape();
        sB1->InspectShape();
    } // end if GetDebug
    //
    TGeoVolume *vM0,*vA0,*vA1,*vTA0,*vTA1,*vB0,*vB1,*vTB0,*vTB1;
    vM0 = new TGeoVolume("ITSSPDSensitiveVirtualvolumeM0",sM0,medSPDair);
    vM0->SetVisibility(kTRUE);
    vM0->SetLineColor(7); // light Blue
    vM0->SetLineWidth(1);
    vM0->SetFillColor(vM0->GetLineColor());
    vM0->SetFillStyle(4090); // 90% transparent
    vA0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorA0",sA0,medSPDcf);
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
    if(GetDebug()){
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
//______________________________________________________________________
void AliITSv11GeometrySPD::HalfStave(TGeoVolume *moth,Double_t &thicknessAA,
				     TGeoManager *mgr){
    // Define the detail SPD Half Stave geometry.
    // Inputs:
    //   TGeoVolume  *moth  The mother volume to place this object.
    //   Int_t      &thicknessAA Thickness of stave at section A-A
    //   TGeoManager *mgr   TGeoManager default gGeoManager
    // Outputs:
    //  none.
    // Return:
    //  none.

    thicknessAA = 1.03*fgkmm; // Default value
    if(moth==0){
      Error("HalfStave","moth=%p mgr=%p",moth,mgr);
        return;
    } // end if moth==0
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
