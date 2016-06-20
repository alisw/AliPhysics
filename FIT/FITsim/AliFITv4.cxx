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


/////////////////////////////////////////////////////////////////////
//                                                                 //
// FIT detector full geometry  version 4                             //
//
//Begin Html       
/*
<img src="gif/AliFITv4Class.gif">
*/
//End Html
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>

#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"


#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVirtualMC.h>
#include <TString.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "AliFITHits.h"
#include "AliFITv4.h"

#include "AliMC.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackReference.h"

ClassImp(AliFITv4)


//--------------------------------------------------------------------
AliFITv4::AliFITv4():  AliFIT(),
		     fIdSens1(0),
		     fIdSens2(0),
		     fPMTeff(0x0)

{
  //
  // Standart constructor for T0 Detector version 0
}
//--------------------------------------------------------------------
AliFITv4::AliFITv4(const char *name, const char *title):
  AliFIT(name,title),
  fIdSens1(0),
  fIdSens2(0),
  fPMTeff(0x0)

{
  //
  // Standart constructor for FIT Detector version 0
  //
  fIshunt = 2; 
  SetPMTeff();
}
//_____________________________________________________________________________

AliFITv4::~AliFITv4() 
{
  // desctructor  
}

//-------------------------------------------------------------------------
void AliFITv4::CreateGeometry()
{
  //
  // Create the geometry of FIT Detector version 1 full geometry
  //
  // begin Html
  //

  Int_t *idtmed = fIdtmed->GetArray();
  Float_t zdetC = 85; //center of mother volume
  Float_t zdetA = 333;
  
  Int_t idrotm[999];
  Double_t x,y,z;
  Float_t pstartC[3] = {6., 20 ,5};
  Float_t pstartA[3] = {2.55, 20 ,5};
  Float_t pinstart[3] = {2.95,2.95,4.};
  Float_t  pmcp[3] = {2.949, 2.949, 2.8}; //MCP
  Float_t ptop[3] = {1.324, 1.324, 1.};//cherenkov radiator
  Float_t preg[3] = {1.324, 1.324, 0.05};//photcathode 
 
  Float_t zV0A = 325.;
  Float_t zV0C = 89;
  Float_t pV0Amother[3] = {4.25, 41.25, 0.15};
  Float_t pV0A[3] = {4.3, 41.2, 0.1};

  AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  
  //-------------------------------------------------------------------
  //  T0 volume 
  //-------------------------------------------------------------------
  //C side

  Float_t xc[28] = {-8.95, -3.05, 3.05, 8.95, -14.85,
		    -8.95, -3.05, 3.05, 8.95,  14.85,
		    -14.85, -8.95, 8.95, 14.85, -14.85, 
		    -8.95, 8.95, 14.85, -14.85,-8.95,
		    -3.05, 3.05, 8.95, 14.85, -8.95,
		    -3.05, 3.05, 8.95 };
  
  Float_t yc[28] = {14.85, 14.85,14.85,14.85,8.95,
		    8.95,8.95,8.95,8.95,8.95, 
		    3.05, 3.05,3.05,3.05,-3.05,
		    -3.05,-3.05,-3.05,-8.95,-8.95,
		    -8.95,-8.95,-8.95,-8.95,-14.85, 
		    -14.85,-14.85,-14.85};
  
  /*
  Float_t xa[24] = {-11.8, -5.9, 0, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 0, 5.9, 11.8 };
  
  Float_t ya[24] = { 11.9, 11.9, 11.9, 11.9, 11.9,
		     6.0,   6.0,  6.0, 6.0,  6.0,
		    -0.1, -0.1, 0.1, 0.1,
		    -6.0, -6.0, -6.0, -6.0, -6.0,
		     -11.9, -11.9, -11.9,  -11.9, -11.9 };  
  */  
 Float_t xa[24] = {-11.8, -5.9, 0, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -12.8, -6.9, 6.9, 12.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 0, 5.9, 11.8 };
  
  Float_t ya[24] = { 11.9, 11.9, 12.9, 11.9, 11.9,
		     6.0,   6.0,  7.0, 6.0,  6.0,
		    -0.1, -0.1, 0.1, 0.1,
		    -6.0, -6.0, -7.0, -6.0, -6.0,
		     -11.9, -11.9, -12.9,  -11.9, -11.9 };  
  
 
  
   
  TGeoVolumeAssembly * stlinA = new TGeoVolumeAssembly("0STL");  // A side mother 
  TGeoVolumeAssembly * stlinC = new TGeoVolumeAssembly("0STR");  // C side mother 
 //T0 interior
  TVirtualMC::GetMC()->Gsvolu("0INS","BOX",idtmed[kAir],pinstart,3);
  TGeoVolume *ins = gGeoManager->GetVolume("0INS");
 // 
  TGeoTranslation *tr[52];
  TString nameTr;

  for (Int_t itr=0; itr<24; itr++) {
    nameTr = Form("0TR%i",itr+1);
    z=-pstartA[2]+pinstart[2];
    tr[itr] = new TGeoTranslation(nameTr.Data(),xa[itr],ya[itr], z );
    printf(" itr %i A %f %f %f \n",itr, xa[itr], ya[itr], z+zdetA);
    tr[itr]->RegisterYourself();
    stlinA->AddNode(ins,itr,tr[itr]);
  }
  for (Int_t itr=24; itr<52; itr++) {
    z=-pstartC[2]+pinstart[2];
    tr[itr] = new TGeoTranslation(nameTr.Data(),xc[itr-24],yc[itr-24], z );
    tr[itr]->RegisterYourself();
    stlinC->AddNode(ins,itr,tr[itr]);
    printf(" itr %i C %f %f %f \n",itr, xc[itr-24], yc[itr-24], z+zdetC);
  }
  TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(stlinA,1,new TGeoTranslation(0,0, zdetA ) );
  //  alice->AddNode(stlinC,1,new TGeoTranslation(0,0, zdetC ) );
  TGeoRotation * rotC = new TGeoRotation( "rotC",90., 0., 90., 90., 180., 0.);
  alice->AddNode(stlinC,1, new TGeoCombiTrans(0., 0., -zdetC , rotC) );
  
  x=0;
  y=0;
   
  // Entry window (glass)
  TVirtualMC::GetMC()->Gsvolu("0TOP","BOX",idtmed[kOpGlass],ptop,3); //glass
  TGeoVolume *top = gGeoManager->GetVolume("0TOP");
  TVirtualMC::GetMC()->Gsvolu ("0REG", "BOX", idtmed[kOpGlassCathode], preg, 3); 
  TGeoVolume *cat = gGeoManager->GetVolume("0REG");
  TVirtualMC::GetMC()->Gsvolu("0MCP","BOX",idtmed[kGlass],pmcp,3); //glass
  TGeoVolume *mcp = gGeoManager->GetVolume("0MCP");
 
  Int_t ntops=0;
   Float_t xin=0, yin=0;
   for (Int_t ix=0; ix<2; ix++) {
     xin = - pinstart[0] + 0.3 + (ix+0.5)*2*ptop[0] ;
     for (Int_t iy=0; iy<2 ; iy++) {
       z = - pinstart[2]+ptop[2];
       yin = - pinstart[1] + 0.3 + (iy+0.5)*2*ptop[1];
       ntops++;
       ins->AddNode(top, ntops, new TGeoTranslation(xin,yin,z) );
       // printf(" 0TOP  full x %f y %f z %f \n", xin, yin, z);
       z = -pinstart[2] + 2 * ptop[2] + preg[2];
       ins->AddNode(cat, ntops, new TGeoTranslation(xin,yin,z) );

       // printf(" GEOGEO  %i %i %i %f %f %f %f %f %f \n", ntops, ix, iy,
       //  xin,yin,x1[ntops],y1[ntops],x1[ntops]+xin,y1[ntops]+yin);
     }
   }
// MCP
   z=-pinstart[2] + 2*ptop[2] + 2*preg[2] + pmcp[2];
  ins->AddNode(mcp, 1 , new TGeoTranslation(0,0,z) );

 //V0A && V0C
   TVirtualMC::GetMC()->Gsvolu("0V0AM","TUBE",idtmed[kAir],pV0Amother,3);
   TVirtualMC::GetMC()->Gspos ("0V0AM",1, "ALIC", 0,0,zV0A , 0, "ONLY");
   TVirtualMC::GetMC()->Gspos ("0V0AM",2, "ALIC", 0,0,-zV0C , 0, "ONLY");
   TVirtualMC::GetMC()->Gsvolu("0V0A","TUBE",idtmed[kSensAir],pV0A,3);
   TVirtualMC::GetMC()->Gspos ("0V0A",1, "0V0AM", 0, 0, 0, 0, "ONLY");


 
}    
//------------------------------------------------------------------------
void AliFITv4::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //
  TString volPath;
  TString symName, sn;
  TString vpAalign = "/ALIC_1/0STL_1";
  TString vpCalign = "/ALIC_1/0STR_1";
  for (Int_t imod=0; imod<2; imod++)  {
    if (imod==0) {volPath  = vpCalign; symName="/ALIC_1/0STL"; }
    if (imod==1) {volPath  = vpAalign; symName="/ALIC_1/0STR"; }
    
    AliDebug(2,"--------------------------------------------");
    AliDebug(2,Form("volPath=%s\n",volPath.Data()));
    AliDebug(2,Form("symName=%s\n",symName.Data()));
    AliDebug(2,"--------------------------------------------");
    if(!gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data()))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",
symName.Data(),volPath.Data()));
   }
}   
//------------------------------------------------------------------------
void AliFITv4::CreateMaterials()
{

   Int_t isxfld   = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
   Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
   //   Float_t a,z,d,radl,absl,buf[1];
   // Int_t nbuf;
// AIR
                                                                                
   Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
   Float_t zAir[4]={6.,7.,8.,18.};
   Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
   Float_t dAir = 1.20479E-3;
   Float_t dAir1 = 1.20479E-11;
 // Radiator  glass SiO2
   Float_t aglass[2]={28.0855,15.9994};
   Float_t zglass[2]={14.,8.};
   Float_t wglass[2]={1.,2.};
   Float_t dglass=2.65;
 // MCP glass SiO2
   Float_t dglass_mcp=1.3;
//*** Definition Of avaible T0 materials ***
   AliMixture(1, "Vacuum$", aAir, zAir, dAir1,4,wAir);
   AliMixture(2, "Air$", aAir, zAir, dAir,4,wAir);
   AliMixture( 4, "MCP glass   $",aglass,zglass,dglass_mcp,-2,wglass);
   AliMixture( 24, "Radiator Optical glass$",aglass,zglass,dglass,-2,wglass);
   
   AliMedium(1, "Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(6, "Glass$", 4, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
  
   AliMedium(16, "OpticalGlass$", 24, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
    AliMedium(19, "OpticalGlassCathode$", 24, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(22, "SensAir$", 2, 1, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   
   AliDebugClass(1,": ++++++++++++++Medium set++++++++++");
   
   
}

//-------------------------------------------------------------------
void AliFITv4::DefineOpticalProperties()
{


// Optical properties definition.
   Int_t *idtmed = fIdtmed->GetArray();
// Definition Cherenkov parameters
   int i;
   const Int_t kNbins=31;
   
   Float_t rindexSiO2[kNbins],  efficAll[kNbins], rindexAir[kNbins], absorAir[kNbins],rindexCathodeNext[kNbins], absorbCathodeNext[kNbins];
   Double_t efficMet[kNbins], aReflMet[kNbins];
           
   // quartz 20mm
   Float_t aAbsSiO2[kNbins]={29.0, 28.6, 28.3, 27.7, 27.3, 26.7, 26.4, 
			     25.9, 25.3, 24.9, 24.5, 23.7, 
			     23.2, 22.8, 22.4, 21.8, 21.3,
			     22.8, 22.1, 21.7, 21.2, 20.5, 
			     19.9, 19.3, 18.7, 18.0, 17.1, 
			     16.3, 15.3, 14.3, 14.3   };
          
   Float_t aPckov[kNbins]  ={3.87, 3.94, 4.02, 4.11, 4.19, 4.29, 4.38,
			     4.48, 4.58, 4.69, 4.81, 4.93, 
			     5.05, 5.19, 5.33, 5.48, 5.63,
			     5.8,  5.97, 6.16, 6.36, 6.57, 
			     6.8,  7.04, 7.3,  7.58, 7.89, 
			     8.22, 8.57, 8.97, 9.39 };  
  Double_t dPckov[kNbins]  ={3.87, 3.94, 4.02, 4.11, 4.19, 4.29, 4.38,
                             4.48, 4.58, 4.69, 4.81, 4.93,
                             5.05, 5.19, 5.33, 5.48, 5.63,
                             5.8,  5.97, 6.16, 6.36, 6.57,
                             6.8,  7.04, 7.3,  7.58, 7.89,
                             8.22, 8.57, 8.97, 9.39 };

   /*     
   Float_t effCathode[kNbins]={0.11, 0.13, 0.15, 0.16, 0.18, 0.19, 0.20,
			      0.21, 0.22, 0.23, 0.24, 0.26, 
			      0.27, 0.29, 0.30, 0.29, 0.29, 
			      0.28, 0.28, 0.27, 0.26, 0.25, 
			      0.25, 0.23, 0.20, 0.19, 0.17,
			      0.17, 0.17, 0.2, 0.23};
   */     
   //  Float_t aAbsSiO2[kNbins]; //quartz 30mm
 for(i=0;i<kNbins;i++)

    {
      aPckov[i]=aPckov[i]*1e-9;//Photons energy bins 4 eV - 8.5 eV step 0.1 eV
      dPckov[i]=dPckov[i]*1e-9;//Photons energy bins 4 eV - 8.5 eV step 0.1 eV 
      //      rindexAir[i]=0.0001;
      rindexAir[i] = 1.;
      rindexSiO2[i]=1.458; //refractive index for qwarts
      rindexCathodeNext[i]=0;
      efficAll[i]=1.;
      efficMet[i]=0.;
      aReflMet[i]=1.;
      //      aAbsSiO2[i]=28.5; //quartz 30mm
      absorAir[i]=0.3;      
      absorbCathodeNext[i]=0;

    }
  
  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlass], kNbins, aPckov, aAbsSiO2, efficAll, rindexSiO2 );
   // TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2, effCathode, rindexSiO2 );
   TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2,efficAll , rindexSiO2 );
   //  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAir], kNbins, aPckov,absorAir , efficAll,rindexAir );
   //   TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAirNext], kNbins, aPckov,absorbCathodeNext , efficAll, rindexCathodeNext);

   //Define a boarder for radiator optical properties
   TVirtualMC::GetMC()->DefineOpSurface("surfRd", kUnified /*kGlisur*/,kDielectric_metal,kPolished, 0.);
   TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "EFFICIENCY", kNbins, dPckov, efficMet);
   TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "REFLECTIVITY", kNbins, dPckov, aReflMet);


}

//-------------------------------------------------------------------
void AliFITv4::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
  AliFIT::Init();
  fIdSens1=TVirtualMC::GetMC()->VolId("0REG");
  fIdSens2=TVirtualMC::GetMC()->VolId("0V0A");

   AliDebug(1,Form("%s: *** FIT version 1 initialized ***\n",ClassName()));
}

//-------------------------------------------------------------------

void AliFITv4::StepManager()
{
  //
  // Called for every step in the T0 Detector
  //
  Int_t id,copy,copy1;
  static Float_t hits[6];
  static Int_t vol[3];
  TLorentzVector pos;
  TLorentzVector mom;
  //   TClonesArray &lhits = *fHits;
  
  if(!TVirtualMC::GetMC()->IsTrackAlive()) return; // particle has disappeared
  
  id=TVirtualMC::GetMC()->CurrentVolID(copy);  
  // Check the sensetive volume
  if(id==fIdSens1 ) { 
    if(TVirtualMC::GetMC()->IsTrackEntering()) {
      TVirtualMC::GetMC()->CurrentVolOffID(1,copy1);
      vol[1] = copy1;
      vol[0]=copy;
      TVirtualMC::GetMC()->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      if(pos[2]<0) vol[2] = 0;
      else vol[2] = 1 ;
      //      printf(" volumes pmt %i mcp %i side %i x %f y %f z %f\n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2] );
      
      Float_t etot=TVirtualMC::GetMC()->Etot();
      hits[3]=etot;
      Int_t iPart= TVirtualMC::GetMC()->TrackPid();
      Int_t partID=TVirtualMC::GetMC()->IdFromPDG(iPart);
      hits[4]=partID;
      Float_t ttime=TVirtualMC::GetMC()->TrackTime();
      hits[5]=ttime*1e12;
      if (TVirtualMC::GetMC()->TrackPid() == 50000050)   // If particles is photon then ...
	{
	  if(RegisterPhotoE(hits[3])) {
	    fIshunt = 2;
	    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
	    // Create a track reference at the exit of photocatode
	  }	    
	}
      //charge particle HITS
      if ( TVirtualMC::GetMC()->TrackCharge() ) {
	fIshunt = 0;
	AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
	//	printf("T0 :::volumes pmt %i mcp %i vol %i x %f y %f z %f particle %f all \n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2], hits[4]);
      }
            
      //charge particle TrackReference
      if ( TVirtualMC::GetMC()->TrackCharge() )
	AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kFIT);
      
    }// trck entering		
  } //sensitive
  //V0A
  if(id==fIdSens2 ) { 
    if ( TVirtualMC::GetMC()->TrackCharge()  ) {
      if(TVirtualMC::GetMC()->IsTrackEntering()) {
	TVirtualMC::GetMC()->TrackPosition(pos);
	hits[0] = pos[0];
	hits[1] = pos[1];
	hits[2] = pos[2];
	vol[0]=0;
	vol[1]=0;
	vol[2]=2;
	
	Float_t etot=TVirtualMC::GetMC()->Etot();
	hits[3]=etot;
	Int_t iPart= TVirtualMC::GetMC()->TrackPid();
	Int_t partID=TVirtualMC::GetMC()->IdFromPDG(iPart);
	hits[4]=partID;
	Float_t ttime=TVirtualMC::GetMC()->TrackTime();
	hits[5]=ttime*1e12;
	fIshunt = 0;
	AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
	//	printf(" V0 ::: volumes pmt %i mcp %i vol %i x %f y %f z %f particle %f all \n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2], hits[4]);
      //charge particle TrackReference
      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kFIT);
    }
  }    
}

}


//------------------------------------------------------------------------
Bool_t AliFITv4::RegisterPhotoE(Double_t energy)
{
  
  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Double_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/energy;
  Float_t eff = fPMTeff->Eval(lambda);
  Double_t  p = gRandom->Rndm();
  
  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}

//----------------------------------------------------------------------------

void AliFITv4::SetPMTeff()
{
  Float_t lambda[50];
  Float_t eff[50 ] = {0,        0,       0.23619,  0.202909, 0.177913, 
		    0.175667, 0.17856, 0.190769, 0.206667, 0.230286,
		    0.252276, 0.256267,0.26,     0.27125,  0.281818,
		    0.288118, 0.294057,0.296222, 0.301622, 0.290421, 
		    0.276615, 0.2666,  0.248,    0.23619,  0.227814, 
		    0.219818, 0.206667,0.194087, 0.184681, 0.167917, 
		    0.154367, 0.1364,  0.109412, 0.0834615,0.0725283, 
		    0.0642963,0.05861, 0.0465,   0.0413333,0.032069, 
		    0.0252203,0.02066, 0.016262, 0.012,    0.00590476,
		    0.003875, 0.00190, 0,        0,        0          } ;
  for (Int_t i=0; i<50; i++) lambda[i]=200+10*i; 

  fPMTeff = new TGraph(50,lambda,eff);
 
}





