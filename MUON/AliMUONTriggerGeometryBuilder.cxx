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

// $Id$
//
// Class AliMUONTriggerGeometryBuilder
// -----------------------------------
// MUON Trigger stations geometry construction class.
//
// Author: Philippe Crochette, LPC Clermont-Ferrand

#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliMUONTriggerGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONChamberGeometry.h"

ClassImp(AliMUONTriggerGeometryBuilder)

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(&muon->Chamber(10), &muon->Chamber(11),&muon->Chamber(12),&muon->Chamber(13)),
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder(const AliMUONTriggerGeometryBuilder& rhs)
  : AliMUONVGeometryBuilder(rhs)
{
// Protected copy constructor

  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::~AliMUONTriggerGeometryBuilder() {
//
}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder& 
AliMUONTriggerGeometryBuilder::operator = (const AliMUONTriggerGeometryBuilder& rhs) 
{
// Protected assignement operator

  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::CreateGeometry()
{
// From AliMUONv1::CreateGeometry()

    /* 
       zpos1 and zpos2 are the middle of the first and second
       planes of station 1 (+1m for second station):
       zpos1=(zpos1m+zpos1p)/2=(15999+16071)/2=16035 mm, thick/2=40 mm
       zpos2=(zpos2m+zpos2p)/2=(16169+16241)/2=16205 mm, thick/2=40 mm
       zposxm and zposxp= middles of gaz gaps within a detection plane
       rem: the total thickness accounts for 1 mm of al on both
       side of the RPCs (see zpos1 and zpos2)
    */
    
    Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
    Int_t idAir= idtmed[1100]; // medium 1
    Int_t idAlu1=idtmed[1103]; // medium 4
    Int_t detElementNumber=0;          // Detection Element Number    
    Float_t tpar[3];
    Double_t dpar[3];    
    
// vertical gap between right and left chambers (kDXZERO*2=4cm)
    const Float_t kDXZERO=2.; 
// main distances for chamber definition in first plane/first station
    const Float_t kXMIN=34.;       
    const Float_t kXMED=51.;                                
    const Float_t kXMAX=255.; 
// 090704 kXMAX changed from 272 to 255.
// (see fig.2-4 & 2-5 of Local Trigger Board PRR)
// segmentation updated accordingly
    const Float_t kYMIN=34.;                              
    const Float_t kYMAX=51.;                              
// inner/outer radius of flange between beam shield. and chambers (1/station)
    const Float_t kRMIN[2]={50.,50.};
    const Float_t kRMAX[2]={64.,68.};
// z position of the middle of the gas gap in mother vol 
    const Float_t kZm=-3.6;
    const Float_t kZp=+3.6;     

    AliMUONChamber *iChamber, *iChamber1;
    iChamber1 = GetChamber(10);
    Float_t zpos1=-iChamber1->Z(); 
    
// ratio of zpos1m/zpos1p and inverse for first plane
    Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
    Float_t zpm=1./zmp;
    
    Int_t icount=0;  // chamber counter (0 1 2 3)
    
    for (Int_t istation=0; istation<2; istation++) { // loop on stations
	for (Int_t iplane=0; iplane<2; iplane++) { // loop on detection planes
	    
	    Int_t iVolNum=1; // counter Volume Number
	    icount = Int_t(iplane*TMath::Power(2,0))+
		Int_t(istation*TMath::Power(2,1));
	    
	    iChamber = GetChamber(10+icount);
	    Float_t zpos = - iChamber->Z();	     
	    
// Flange between beam shielding and RPC 
	    tpar[0]= kRMIN[istation];
	    tpar[1]= kRMAX[istation];
	    tpar[2]= 4.0;	    
	    char volFlange[5];
	    sprintf(volFlange,"SF%dA",icount+1);	 
	    gMC->Gsvolu(volFlange,"TUBE",idAlu1,tpar,3);     // Al
            // changed by ivana
	    //gMC->Gspos(volFlange,1,"ALIC",0.,0.,zpos,0,"MANY");
	    iChamber->GetGeometry()->AddEnvelope(volFlange, false, "MANY");
	    
// scaling factor
	    Float_t zRatio = zpos / zpos1;
	    
// envelopes (same size except line 5, all virtual)
	    char volEnv[18][5];
	    tpar[1] = kYMIN * zRatio; 
	    tpar[2] = 0.4;
	    Int_t i=0;    // counter
	    for (Int_t icolumn=0; icolumn<2; icolumn++) {
		for (Int_t iline=1; iline<10; iline++){
		    tpar[0] = (kXMAX/2.) * zRatio;
		    if (iline==5) tpar[0] = ((kXMAX-kXMED)/2.)*zRatio;
		    if (icolumn==0) 
			sprintf(volEnv[i],"S%dR%d",icount,iline);
		    else
			sprintf(volEnv[i],"S%dL%d",icount,iline);
		    gMC->Gsvolu(volEnv[i],"BOX",idAir,tpar,0); 
		    i++;
		}
	    }

// chamber prototype
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    
	    char volAlu[5];     // Alu 
	    char volBak[5];     // Bakelite
	    char volGaz[5];     // Gas streamer	    
	    sprintf(volAlu,"SC%dA",icount+1);
	    sprintf(volBak,"SB%dA",icount+1);
	    sprintf(volGaz,"SG%dA",icount+1);
	    gMC->Gsvolu(volAlu,"BOX",idAlu1,tpar,0);         // Al
	    gMC->Gsvolu(volBak,"BOX",idtmed[1107],tpar,0);   // Bakelite
	    gMC->Gsvolu(volGaz,"BOX",idtmed[1106],tpar,0);   // Gas streamer
	    tpar[0] = -1.;
	    tpar[1] = -1.;
	    tpar[2] = 0.1;    
	    gMC->Gsposp(volGaz,1,volBak,0.,0.,0.,0,"ONLY",tpar,3);
	    tpar[2] = 0.3;
	    gMC->Gsposp(volBak,1,volAlu,0.,0.,0.,0,"ONLY",tpar,3);

// chamber type A
	    Float_t xEnv = (kDXZERO+kXMED+(kXMAX-kXMED)/2.)*zRatio;
	    Float_t yEnvM = 0.;	 // y low position of envelope in chamber
	    Float_t yEnvP = 0.;	 // y up position of envelope in chamber
	    Float_t yEnvPsave = 0.; // tmp data
	    Float_t yEnvMsave = 0.; // tmp data
	    Float_t xpos = 0.; // x position of RPC in envelope	    
	    Float_t ypos = 0.; // y position of RPC in envelope
	    dpar[2] = 0.4;	    
	    dpar[0] = ((kXMAX-kXMED)/2.)*zRatio;
	    dpar[1] = kYMIN * zRatio;

	    detElementNumber = (10+icount+1)*100+4;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[4], true, TGeoTranslation(xEnv,yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+50+4;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[13], true, TGeoTranslation(-xEnv,yEnvP,kZp));

	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[4],iVolNum++,3, dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[13],iVolNum++,3, dpar);	    

// chamber type B (plus envelope chambers B & C)   
	    xEnv = (kDXZERO+kXMAX/2.)*zRatio;
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;
	    dpar[0] = ((kXMAX-kXMIN)/2.) * zRatio;
	    dpar[1] = ((kYMAX-kYMIN)/2.) * zRatio;
	    xpos = kXMIN/2. * zRatio;
	    ypos = (kYMIN - kYMIN/4.) * zRatio;

	    detElementNumber = (10+icount+1)*100+3;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[3], true, TGeoTranslation( xEnv,-yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+5;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[5], true, TGeoTranslation( xEnv, yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+50+3;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[12], true, TGeoTranslation(-xEnv,-yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+50+5;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[14], true, TGeoTranslation(-xEnv, yEnvM,kZm));

	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[3],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[5],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[12],iVolNum++,TGeoTranslation(-xpos, ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[14],iVolNum++,TGeoTranslation(-xpos,-ypos,0.),3,dpar);

// chamber type C (note: same Z than type B)
	    dpar[0] = (kXMAX/2)*zRatio;
	    dpar[1] = (kYMAX/2)*zRatio;
	    xpos = 0.;	    
	    ypos = ((kYMAX - kYMIN)/2.) * zRatio;

	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[3],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[5],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[12],iVolNum++,TGeoTranslation(-xpos,-ypos,0.),3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[14],iVolNum++,TGeoTranslation(-xpos, ypos,0.),3,dpar);
    
// chamber type D, E and F (same size)
// D	    
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;
	    dpar[0] = (kXMAX/2.)*zRatio;
	    dpar[1] =  kYMIN*zRatio;

	    detElementNumber = (10+icount+1)*100+2;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[2], true, TGeoTranslation(xEnv,-yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+6;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[6], true, TGeoTranslation(xEnv, yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+50+2;
            GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[11], true, TGeoTranslation(-xEnv,-yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+50+6;
            GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[15], true, TGeoTranslation(-xEnv, yEnvP,kZp));

	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[2],iVolNum++,3, dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[6],iVolNum++,3, dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[11],iVolNum++,3, dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[15],iVolNum++,3, dpar);

// E
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;

	    detElementNumber = (10+icount+1)*100+1;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[1], true, TGeoTranslation(xEnv,-yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+7;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[7], true, TGeoTranslation(xEnv, yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+50+1;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[10], true, TGeoTranslation(-xEnv,-yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+50+7;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[16], true, TGeoTranslation(-xEnv, yEnvM,kZm));

	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[1],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[7],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[10],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[16],iVolNum++,3,dpar);


// F
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;

	    detElementNumber = (10+icount+1)*100;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[0], true, TGeoTranslation(xEnv,-yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+8;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[8], true, TGeoTranslation(xEnv, yEnvM,kZm));
	    detElementNumber = (10+icount+1)*100+50;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[9], true, TGeoTranslation(-xEnv,-yEnvP,kZp));
	    detElementNumber = (10+icount+1)*100+50+8;
	    GetChamber(10+icount)->GetGeometry()->AddEnvelope(volEnv[17], true, TGeoTranslation(-xEnv, yEnvP,kZp));
	    
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[0],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[8],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[9],iVolNum++,3,dpar);
	    GetChamber(10+icount)->GetGeometry()->AddEnvelopeConstituentParam(volAlu,volEnv[17],iVolNum++,3,dpar);

	} // end loop on detection planes
    } // end loop on stations    
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::SetTransformations() 
{
// Defines the transformations for the trigger chambers.
// ---
    Double_t zpos1, zpos2;    
    AliMUONChamber *iChamber1, *iChamber2;

    iChamber1 = GetChamber(10);
    zpos1= - iChamber1->Z(); 
    iChamber1->GetGeometry()
	->SetTranslation(TGeoTranslation(0., 0., zpos1));
    
    iChamber2 = GetChamber(11);
    zpos2 = - iChamber2->Z(); 
    iChamber2->GetGeometry()
	->SetTranslation(TGeoTranslation(0., 0., zpos2));

    iChamber1 = GetChamber(12);
    zpos1 = - iChamber1->Z(); 
    iChamber1->GetGeometry()
	->SetTranslation(TGeoTranslation(0., 0., zpos1));
    
    iChamber2 = GetChamber(13);
    zpos2 = - iChamber2->Z(); 
    iChamber2->GetGeometry()
	->SetTranslation(TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::SetSensitiveVolumes()
{
// Defines the sensitive volumes for trigger station chambers.
// ---

  GetChamber(10)->GetGeometry()->SetSensitiveVolume("SG1A");
  GetChamber(11)->GetGeometry()->SetSensitiveVolume("SG2A");
  GetChamber(12)->GetGeometry()->SetSensitiveVolume("SG3A");
  GetChamber(13)->GetGeometry()->SetSensitiveVolume("SG4A");
}

