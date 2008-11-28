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

//-----------------------------------------------------------------------------
/// Class AliMUONTriggerGeometryBuilder
// -----------------------------------
// MUON Trigger stations geometry 
// construction class.
// Author: Philippe Crochet (LPCCFd)
// Support for trigger chambers added April 07 by Enrico Scomparin (INFN To)
//-----------------------------------------------------------------------------

#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliLog.h"
#include "AliRun.h"

#include "AliMUONTriggerGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerGeometryBuilder)
/// \endcond

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(16, 4),
   fMUON(muon)
{
/// Standard constructor

}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::~AliMUONTriggerGeometryBuilder() 
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::CreateGeometry()
{
/// From AliMUONv1::CreateGeometry()

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
    Int_t idInox = idtmed[1128];       // medium 29 Stainless Steel (18%Cr,9%Ni,Fe) 
    
    Int_t detElemId=0;          // Detection Element Number    
    Float_t tpar[3];
    Double_t dpar[3];    
    Double_t spar[3];    
    Double_t ppar[3];    
   
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
//    const Float_t kRMIN[2]={50.,50.};
//    const Float_t kRMAX[2]={64.,68.};
// z position of the middle of the gas gap in mother vol 
    const Float_t kZm=-3.6;
    const Float_t kZp=+3.6;
    
// y positions of vertical supports
    const Float_t kYVSup[4]={61.45,122.45,192.95,236.95}; 
// dimensions of vertical supports 
    const Float_t kSizeVSupExt[3]={1.5,1.5,306.+5.}; 
    const Float_t kSizeVSupInt[3]={1.2,1.2,306.+5.};  
// transverse dimensions of angular supports 
    const Float_t kSizeSupport1V[3]={0.,1.5,0.1}; 
    const Float_t kSizeSupport1H[3]={0.,0.1,1.15}; // z should be 1.4 in the installed set-up 
    const Float_t kSizeSupport2V[3]={0.,3.0,0.1}; 
    const Float_t kSizeSupport2H[3]={0.,0.1,1.9}; 
    const Float_t kSizeSupportXV[3]={0.,1.25,0.25}; 
    const Float_t kSizeSupportXH[3]={0.,0.25,1.5}; 
// transverse dimensions of horizontal cable supports
    const Float_t kSizeSupportCable[3]={0.,2.,3.}; 
// dimensions of gas pipes (inner and outer radius)
    const Float_t kSizeGasPipe[3]={0.2,0.4,0.}; 
// Position of gas pipe with respect to angular support
    const Float_t kOffsetGasPipe=0.75; 
// Small cut on some volumes to avoid extrusion from SC1x
    const Float_t kAvoidExtrusion=2.9;    
   
    Float_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    Double_t dstation =  ( ( - AliMUONConstants::DefaultChamberZ(11)) - 
                           ( - AliMUONConstants::DefaultChamberZ(10)) ) /2.1;
    Float_t par[3];
    par[2] = dstation;

// ratio of zpos1m/zpos1p and inverse for first plane
    Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
    Float_t zpm=1./zmp;
    
    Int_t icount=0;  // chamber counter (0 1 2 3)
    
    for (Int_t istation=0; istation<2; istation++) { // loop on stations
	for (Int_t iplane=0; iplane<2; iplane++) { // loop on detection planes
	    
	    Int_t iVolNum=1; // counter Volume Number
	    icount = Int_t(iplane<<0)+Int_t(istation<<1);
	    
	    par[0] = AliMUONConstants::Rmin(5+istation); 
	    par[1] = AliMUONConstants::Rmax(5+istation);
	    Char_t volName[6];
	    sprintf(volName,"%s%d", "SC",11+icount);
 	    gMC->Gsvolu(volName,"TUBE", idAir, par, 3);
 	    //SetVolume(10+icount, volName);
//	    Float_t zpos =  AliMUONConstants::DefaultChamberZ(10+icount);

/* removed 03/18/05
// Flange between beam shielding and RPC 
	    tpar[0]= kRMIN[istation];
	    tpar[1]= kRMAX[istation];
	    tpar[2]= 4.0;	    
	    char volFlange[5];
	    sprintf(volFlange,"SF%dA",icount+1);	 
	    gMC->Gsvolu(volFlange,"TUBE",idAlu1,tpar,3);     // Al
            // changed by ivana
	    //gMC->Gspos(volFlange,1,"ALIC",0.,0.,zpos,0,"MANY");
	    iChamber->GetGeometry()->GetEnvelopeStore()
	      ->AddEnvelope(volFlange, 0, false, "MANY");
*/
	    
// scaling factor
//	    Float_t zRatio = zpos / zpos1;
	    Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);	    	    



    
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
		    // gMC->Gsvolu(volEnv[i],"BOX",idAir,tpar,0); 
		    i++;
		}
	    }

// chamber prototype
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    char volAlu[5];     // Alu 
	    char volBak[5];     // Bakelite
	    char volGaz[5];     // Gas streamer	    
	    sprintf(volAlu,"SC%dA",icount+1);
	    sprintf(volBak,"SB%dA",icount+1);
	    sprintf(volGaz,"S%dG",icount+11);
	    gMC->Gsvolu(volAlu,"BOX",idAlu1,tpar,0);         // Al
	    gMC->Gsvolu(volBak,"BOX",idtmed[1107],tpar,0);   // Bakelite
	    gMC->Gsvolu(volGaz,"BOX",idtmed[1106],tpar,0);   // Gas streamer
	    tpar[0] = -1.;
	    tpar[1] = -1.;
	    tpar[2] = 0.1;    
	    gMC->Gsposp(volGaz,1,volBak,0.,0.,0.,0,"ONLY",tpar,3);
	    tpar[2] = 0.3;
	    gMC->Gsposp(volBak,1,volAlu,0.,0.,0.,0,"ONLY",tpar,3);

// RPC supports (vertical)

            char volAluSupport[5],volAirSupport[5];
	    sprintf(volAluSupport,"SAL%d",icount+1);
	    sprintf(volAirSupport,"SAI%d",icount+1);
	    char volEnvSupport[12][7];
	    for(Int_t ii=0;ii<8;ii++){
	      sprintf(volEnvSupport[ii],"SEA%dV%d",icount+1,ii);
	    }
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    
            gMC->Gsvolu(volAluSupport,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAirSupport,"BOX",idAir,tpar,0);
	    tpar[0]=kSizeVSupInt[0];
	    tpar[1]=kSizeVSupInt[1];
	    tpar[2]=-1.;
	    gMC->Gsposp(volAirSupport,1,volAluSupport,0.,0.,0.,0,"ONLY",tpar,3);
	    
	    TGeoRotation rsupportv;
	    rsupportv.SetAngles(0.,90.,0.);
	    dpar[0]=kSizeVSupExt[0];
	    dpar[1]=kSizeVSupExt[1];
	    dpar[2]=kSizeVSupExt[2]*zRatio;
	    for(Int_t ii=0;ii<4;ii++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupport[ii], 0, true,
	      TGeoTranslation(-kYVSup[ii]*zRatio,0.,0.),rsupportv);
	      GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupport,volEnvSupport[ii],iVolNum++,3, dpar);
	    }	    
	    for(Int_t ii=4;ii<8;ii++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupport[ii], 0, true,
	      TGeoTranslation(kYVSup[ii-4]*zRatio,0.,0.),rsupportv);
	      GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupport,volEnvSupport[ii],iVolNum++,3, dpar);
	    }	
	        
// RPC supports (horizontal)

// supports for cables
            char volAluSupportH[6],volAirSupportH[6];
	    sprintf(volAluSupportH,"SALH%d",icount+1);
	    sprintf(volAirSupportH,"SAIH%d",icount+1);
	    char volEnvSupportHA[6][8],volEnvSupportHBC[12][8],volEnvSupportHD[12][8],volEnvSupportHE[12][8],volEnvSupportHF[12][8];
            for(Int_t jj=0;jj<2;jj++){
	      for(Int_t ii=0;ii<6;ii++){
	        if(ii<3)sprintf(volEnvSupportHA[3*jj+ii],"SA%dHA%d",icount+1,3*jj+ii);
	        sprintf(volEnvSupportHBC[6*jj+ii],"SA%dHB%d",icount+1,6*jj+ii);
	        sprintf(volEnvSupportHD[6*jj+ii],"SA%dHD%d",icount+1,6*jj+ii);
	        sprintf(volEnvSupportHE[6*jj+ii],"SA%dHE%d",icount+1,6*jj+ii);
	        sprintf(volEnvSupportHF[6*jj+ii],"SA%dHF%d",icount+1,6*jj+ii);
	      }
	    }
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    
            gMC->Gsvolu(volAluSupportH,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAirSupportH,"BOX",idAir,tpar,0);
	    tpar[0]=-1.;
	    tpar[1]=1.9;
	    tpar[2]=2.8;
	    gMC->Gsposp(volAirSupportH,1,volAluSupportH,0.,0.,0.,0,"ONLY",tpar,3);
	    
// Angular supports for chambers
	    char volAluAngSupport1V[6],volAluAngSupport1H[6];
	    sprintf(volAluAngSupport1H,"SA1H%d",icount+1);
	    sprintf(volAluAngSupport1V,"SA1V%d",icount+1);
	    char volAluAngSupport2V[6],volAluAngSupport2H[6];
	    sprintf(volAluAngSupport2H,"SA2H%d",icount+1);
	    sprintf(volAluAngSupport2V,"SA2V%d",icount+1);
	    char volAluAngSupport3V[6],volAluAngSupport3H[6];
	    sprintf(volAluAngSupport3H,"SA3H%d",icount+1);
	    sprintf(volAluAngSupport3V,"SA3V%d",icount+1);
	    char volAluAngSupport4V[6],volAluAngSupport4H[6];
	    sprintf(volAluAngSupport4H,"SA4H%d",icount+1);
	    sprintf(volAluAngSupport4V,"SA4V%d",icount+1);
	    char volAluAngSupportXV[6],volAluAngSupportXH[6];
	    sprintf(volAluAngSupportXH,"SAXH%d",icount+1);
	    sprintf(volAluAngSupportXV,"SAXV%d",icount+1);
	    char volEnvSuppAng1HA[2][7],volEnvSuppAng1HBC[4][7],volEnvSuppAng1HD[4][7],volEnvSuppAng1HE[4][7],volEnvSuppAng1HF[4][7];
	    char volEnvSuppAng1VA[2][7],volEnvSuppAng1VBC[4][7],volEnvSuppAng1VD[4][7],volEnvSuppAng1VE[4][7],volEnvSuppAng1VF[4][7];
	    char volEnvSuppAng2HA[2][7],volEnvSuppAng2HBC[4][7],volEnvSuppAng2HD[4][7],volEnvSuppAng2HE[4][7],volEnvSuppAng2HF[4][7];
	    char volEnvSuppAng2VA[2][7],volEnvSuppAng2VBC[4][7],volEnvSuppAng2VD[4][7],volEnvSuppAng2VE[4][7],volEnvSuppAng2VF[4][7];
	    char volEnvSuppAng3HA[2][7],volEnvSuppAng3HBC[4][7],volEnvSuppAng3HD[4][7],volEnvSuppAng3HE[4][7],volEnvSuppAng3HF[4][7];
	    char volEnvSuppAng3VA[2][7],volEnvSuppAng3VBC[4][7],volEnvSuppAng3VD[4][7],volEnvSuppAng3VE[4][7],volEnvSuppAng3VF[4][7];
	    char volEnvSuppAng4HA[2][7],volEnvSuppAng4HBC[4][7],volEnvSuppAng4HD[4][7],volEnvSuppAng4HE[4][7],volEnvSuppAng4HF[4][7];
	    char volEnvSuppAng4VA[2][7],volEnvSuppAng4VBC[4][7],volEnvSuppAng4VD[4][7],volEnvSuppAng4VE[4][7],volEnvSuppAng4VF[4][7];
	    char volEnvSuppAngXHA[2][7],volEnvSuppAngXHBC[4][7],volEnvSuppAngXHD[4][7],volEnvSuppAngXHE[4][7],volEnvSuppAngXHF[4][7];
	    char volEnvSuppAngXVA[2][7],volEnvSuppAngXVBC[4][7],volEnvSuppAngXVD[4][7],volEnvSuppAngXVE[4][7],volEnvSuppAngXVF[4][7];
	    for(Int_t ii=0;ii<4;ii++){
	      if(ii<2)sprintf(volEnvSuppAng1HA[ii],"SH1%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng1HBC[ii],"SH1%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng1HD[ii],"SH1%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng1HE[ii],"SH1%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng1HF[ii],"SH1%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng1VA[ii],"SV1%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng1VBC[ii],"SV1%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng1VD[ii],"SV1%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng1VE[ii],"SV1%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng1VF[ii],"SV1%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng2HA[ii],"SH2%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng2HBC[ii],"SH2%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng2HD[ii],"SH2%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng2HE[ii],"SH2%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng2HF[ii],"SH2%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng2VA[ii],"SV2%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng2VBC[ii],"SV2%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng2VD[ii],"SV2%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng2VE[ii],"SV2%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng2VF[ii],"SV2%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng3HA[ii],"SH3%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng3HBC[ii],"SH3%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng3HD[ii],"SH3%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng3HE[ii],"SH3%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng3HF[ii],"SH3%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng3VA[ii],"SV3%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng3VBC[ii],"SV3%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng3VD[ii],"SV3%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng3VE[ii],"SV3%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng3VF[ii],"SV3%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng4HA[ii],"SH4%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng4HBC[ii],"SH4%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng4HD[ii],"SH4%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng4HE[ii],"SH4%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng4HF[ii],"SH4%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAng4VA[ii],"SV4%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAng4VBC[ii],"SV4%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAng4VD[ii],"SV4%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAng4VE[ii],"SV4%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAng4VF[ii],"SV4%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAngXHA[ii],"SHX%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAngXHBC[ii],"SHX%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAngXHD[ii],"SHX%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAngXHE[ii],"SHX%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAngXHF[ii],"SHX%dF%d",icount+1,ii);
	      if(ii<2)sprintf(volEnvSuppAngXVA[ii],"SVX%dA%d",icount+1,ii);
	      sprintf(volEnvSuppAngXVBC[ii],"SVX%dB%d",icount+1,ii);
	      sprintf(volEnvSuppAngXVD[ii],"SVX%dD%d",icount+1,ii);
	      sprintf(volEnvSuppAngXVE[ii],"SVX%dE%d",icount+1,ii);
	      sprintf(volEnvSuppAngXVF[ii],"SVX%dF%d",icount+1,ii);
	    }
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    
            gMC->Gsvolu(volAluAngSupport1V,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport1H,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport2V,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport2H,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport3V,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport3H,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport4V,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupport4H,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupportXV,"BOX",idAlu1,tpar,0);
            gMC->Gsvolu(volAluAngSupportXH,"BOX",idAlu1,tpar,0);
	    
// gas pipes
	    char volInoxGasPipe[7];
	    sprintf(volInoxGasPipe,"SPINO%d",icount+1);
            char volEnvInoxGasPipe1A[2][7],volEnvInoxGasPipe1BC[4][8],volEnvInoxGasPipe1D[4][7],volEnvInoxGasPipe1E[4][7],volEnvInoxGasPipe1F[4][7];
            char volEnvInoxGasPipe2A[2][7],volEnvInoxGasPipe2BC[4][8],volEnvInoxGasPipe2D[4][7],volEnvInoxGasPipe2E[4][7],volEnvInoxGasPipe2F[4][7];
	    for(Int_t ii=0;ii<4;ii++){
	      if(ii<2)sprintf(volEnvInoxGasPipe1A[ii],"SP1%dA%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe1BC[ii],"SP1%dBC%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe1D[ii],"SP1%dD%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe1E[ii],"SP1%dE%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe1F[ii],"SP1%dF%d",icount+1,ii);
	    }
	    for(Int_t ii=0;ii<4;ii++){
	      if(ii<2)sprintf(volEnvInoxGasPipe2A[ii],"SP2%dA%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe2BC[ii],"SP2%dBC%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe2D[ii],"SP2%dD%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe2E[ii],"SP2%dE%d",icount+1,ii);
	      sprintf(volEnvInoxGasPipe2F[ii],"SP2%dF%d",icount+1,ii);
	    }
	    tpar[0]= 0.;
	    tpar[1]= 0.;
	    tpar[2]= 0.;	    
            gMC->Gsvolu(volInoxGasPipe,"TUBE",idInox,tpar,0);
	    TGeoRotation rsupportpipe;
	    rsupportpipe.SetAngles(90.,90.,0.);
	   

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

	    detElemId = (10+icount+1)*100;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[4], detElemId, true, TGeoTranslation(xEnv,yEnvP,kZp));
	    detElemId = (10+icount+1)*100+9;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[13], detElemId, true, TGeoTranslation(-xEnv,yEnvM,kZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[4],iVolNum++,3, dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[13],iVolNum++,3, dpar);	    

// horizontal cable support chamber type A

	    spar[0]=((kXMAX/2)-kYVSup[0]/2.)*zRatio;
	    spar[1]=kSizeSupportCable[1];
	    spar[2]=kSizeSupportCable[2];
	    Float_t offsetSuppA = ((kXMAX-kXMED)/2.)*zRatio-(((kXMAX/2)-kYVSup[0]/2.)*zRatio);
	    for(Int_t in=0;in<3;in++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHA[in], 0, true,
	      TGeoTranslation(xEnv+offsetSuppA/2.,yEnvP+dpar[1]/2.*(in-1),-(kSizeVSupExt[0]+spar[2])));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHA[in+3], 0, true,
 	      TGeoTranslation(-(xEnv+offsetSuppA/2.),yEnvM+dpar[1]/2.*(in-1),kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    }
	    for(Int_t ii=0;ii<6;ii++) 
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHA[ii],iVolNum++,3, spar);
	    
// angular supports chamber type A
// 1 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    Float_t sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1],kZp+dpar[2]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VA[1],iVolNum++,3, spar);	    

// 1 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HA[1],iVolNum++,3, spar);	    

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1A[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1A[1], 0, true,
	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1A[0],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1A[1],iVolNum++,3, ppar);

// 2 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VA[1],iVolNum++,3, spar);	    

// 2 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2]; 
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HA[1],iVolNum++,3, spar);	    

// 3 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1],kZp+dpar[2]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VA[1],iVolNum++,3, spar);	    

// 3 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HA[1],iVolNum++,3, spar);	
	        
// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2A[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2A[1], 0, true,
	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2A[0],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2A[1],iVolNum++,3, ppar);

// 4 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VA[1],iVolNum++,3, spar);	    

// 4 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HA[0], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HA[1], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HA[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HA[1],iVolNum++,3, spar);	    

// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    Float_t sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHA[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHA[1], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHA[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHA[1],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVA[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVA[1], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVA[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVA[1],iVolNum++,3, spar);	    

// chamber type B (plus envelope chambers B & C)   
	    xEnv = (kDXZERO+kXMAX/2.)*zRatio;
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;
	    dpar[0] = ((kXMAX-kXMIN)/2.) * zRatio;
	    dpar[1] = ((kYMAX-kYMIN)/2.) * zRatio;
	    Float_t dysave = dpar[1];
	    Float_t dxsave = dpar[0];
	    xpos = kXMIN/2. * zRatio;
	    ypos = (kYMIN - kYMIN/4.) * zRatio;
	    Float_t xpossave = xpos;

	    detElemId = (10+icount+1)*100+17;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[3], detElemId, true, TGeoTranslation( xEnv,-yEnvM,kZm));	    
	    detElemId = (10+icount+1)*100+1;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[5], detElemId, true, TGeoTranslation( xEnv, yEnvM,kZm));
	    detElemId = (10+icount+1)*100+10;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[12], detElemId, true, TGeoTranslation(-xEnv,-yEnvP,kZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    detElemId = (10+icount+1)*100+8;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[14], detElemId, true, TGeoTranslation(-xEnv, yEnvP,kZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[3],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[5],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[12],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[14],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);

	    	    
// chamber type C (note: same Z than type B)
	    dpar[0] = (kXMAX/2)*zRatio;
	    dpar[1] = (kYMAX/2)*zRatio;
	    xpos = 0.;	    
	    ypos = ((kYMAX - kYMIN)/2.) * zRatio;

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[3],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[5],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[12],iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[14],iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);

// horizontal cable support chamber type B+C

	    spar[0]=dpar[0]-kYVSup[0]/2.;
	    spar[1]=kSizeSupportCable[1];
	    spar[2]=kSizeSupportCable[2];
	    for(Int_t in=0;in<3;in++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio,-(yEnvM+(kYMAX-kYMIN/2.*zRatio)/2.*(in-1)),kSizeVSupExt[0]+spar[2]));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+3], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio, yEnvM+(kYMAX-kYMIN/2.*zRatio)/2.*(in-1),kSizeVSupExt[0]+spar[2]));
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+6], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio),-(yEnvP+(kYMAX-kYMIN/2.*zRatio)/2.*(in-1)),-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+9], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio), yEnvP+(kYMAX-kYMIN/2.*zRatio)/2.*(in-1),-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    }
	    for(Int_t ii=0;ii<12;ii++) 
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHBC[ii],iVolNum++,3, spar);

// angular supports chamber type B and C
// C	   
// 1 vertical
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];

       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-dysave,kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-dysave,kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC[2],iVolNum++,3, spar);

// 1 horizontal
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
	    
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-dysave-sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-dysave-sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC[2],iVolNum++,3, spar);	    

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-dysave-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-dysave-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC[0],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC[2],iVolNum++,3, ppar);

// 2 vertical	   
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-dysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-dysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC[2],iVolNum++,3, spar);	    

// 2 horizontal	   
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-dysave-sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-dysave-sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC[2],iVolNum++,3, spar);	    

// 3 vertical
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC[0], 0, true,
	    TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+dysave,kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC[2], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+dysave,kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC[2],iVolNum++,3, spar);

// 3 horizontal
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC[0], 0, true,
	    TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+dysave+sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC[2], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+dysave+sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC[0],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC[2],iVolNum++,3, spar);	    
	        
// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dxsave-kAvoidExtrusion;
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC[0], 0, true,
	    TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+dysave+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC[2], 0, true,
	    TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+dysave+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC[0],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC[2],iVolNum++,3, ppar);

// 4 vertical	   
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC[0], 0, true,
	    TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+dysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC[2], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+dysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC[2],iVolNum++,3, spar);	    

// 4 horizontal	   
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC[0], 0, true,
	    TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+dysave+sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC[2], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+dysave+sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC[2],iVolNum++,3, spar);	    

// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvM+dpar[1]+dysave+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvP+dpar[1]+dysave+sparysave+1.0,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC[2],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvM+dpar[1]+dysave+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvP+dpar[1]+dysave+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC[0],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC[2],iVolNum++,3, spar);	    

// B
// 1 vertical
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
            GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC[1], 0, true,
	    TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-dysave,kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC[3], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-dysave,kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC[1],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC[3],iVolNum++,3, spar);


// 1 horizontal
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
	    
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC[1], 0, true,
	    TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-dysave-sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC[3], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-dysave-sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC[1],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC[3],iVolNum++,3, spar);	    

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dxsave-kAvoidExtrusion;
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC[1], 0, true,
	    TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-dysave-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC[3], 0, true,
	    TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-dysave-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC[1],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC[3],iVolNum++,3, ppar);

// 2 vertical	   
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC[1], 0, true,
	    TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-dysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC[3], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-dysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC[3],iVolNum++,3, spar);	    

// 2 horizontal	   
	    spar[0]=dxsave-kAvoidExtrusion;
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC[1], 0, true,
	    TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-dysave-sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC[3], 0, true,
 	    TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-dysave-sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC[3],iVolNum++,3, spar);	    

// 3 vertical
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
            GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+dysave,kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+dysave,kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC[1],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC[3],iVolNum++,3, spar);

// 3 horizontal
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
	    
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+dysave+sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+dysave+sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC[1],iVolNum++,3, spar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC[3],iVolNum++,3, spar);	    

// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+dysave+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+dysave+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC[1],iVolNum++,3, ppar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC[3],iVolNum++,3, ppar);

// 4 vertical	   
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+dysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+dysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC[3],iVolNum++,3, spar);	    

// 4 horizontal	   
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+dysave+sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+dysave+sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC[3],iVolNum++,3, spar);	    

// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvM+dpar[1]+dysave+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvP+dpar[1]+dysave+sparysave+1.0,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC[3],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvM+dpar[1]+dysave+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvP+dpar[1]+dysave+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC[1],iVolNum++,3, spar);	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC[3],iVolNum++,3, spar);	    


// chamber type D, E and F (same size)
// D	    
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;
	    dpar[0] = (kXMAX/2.)*zRatio;
	    dpar[1] =  kYMIN*zRatio;

	    detElemId = (10+icount+1)*100+16;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[2], detElemId, true, TGeoTranslation(xEnv,-yEnvP,kZp));
	    detElemId = (10+icount+1)*100+2;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[6], detElemId, true, TGeoTranslation(xEnv, yEnvP,kZp));
	    detElemId = (10+icount+1)*100+11;
            GetEnvelopes(16+icount)->AddEnvelope(volEnv[11], detElemId, true, TGeoTranslation(-xEnv,-yEnvM,kZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    detElemId = (10+icount+1)*100+7;
            GetEnvelopes(16+icount)->AddEnvelope(volEnv[15], detElemId, true, TGeoTranslation(-xEnv, yEnvM,kZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[2],iVolNum++,3, dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[6],iVolNum++,3, dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[11],iVolNum++,3, dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[15],iVolNum++,3, dpar);

// horizontal cable support chamber type D

	    spar[0]=dpar[0]-(kYVSup[0]/2.)*zRatio;
	    spar[1]=kSizeSupportCable[1];
	    spar[2]=kSizeSupportCable[2];
	    for(Int_t in=0;in<3;in++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio,-(yEnvP+dpar[1]/2.*(in-1)),-(kSizeVSupExt[0]+spar[2])));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+3], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio, yEnvP+dpar[1]/2.*(in-1),-(kSizeVSupExt[0]+spar[2])));
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+6], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio),-(yEnvM+dpar[1]/2.*(in-1)),kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+9], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio),yEnvM+dpar[1]/2.*(in-1),kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    }
	    for(Int_t ii=0;ii<12;ii++) 
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHD[ii],iVolNum++,3, spar);
	    
// angular supports chamber type D
// 1 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1],kZp+dpar[2]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1],kZp+dpar[2]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VD[i],iVolNum++,3, spar);


// 1 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HD[i],iVolNum++,3, spar);

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1D[i],iVolNum++,3, ppar);

// 2 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VD[i],iVolNum++,3, spar);	    

// 2 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HD[i],iVolNum++,3, spar);	    

// 3 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1],kZp+dpar[2]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1],kZp+dpar[2]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VD[i],iVolNum++,3, spar);


// 3 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HD[i],iVolNum++,3, spar);
	        
// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2D[i],iVolNum++,3, ppar);

// 4 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VD[i],iVolNum++,3, spar);	    

// 4 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HD[i],iVolNum++,3, spar);	    

// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHD[i],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVD[i],iVolNum++,3, spar);	    

// E
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;

	    detElemId = (10+icount+1)*100+15;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[1], detElemId, true, TGeoTranslation(xEnv,-yEnvM,kZm));
	    detElemId = (10+icount+1)*100+3;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[7], detElemId, true, TGeoTranslation(xEnv, yEnvM,kZm));
	    detElemId = (10+icount+1)*100+12;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[10], detElemId, true, TGeoTranslation(-xEnv,-yEnvP,kZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    detElemId = (10+icount+1)*100+6;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[16], detElemId, true, TGeoTranslation(-xEnv, yEnvP,kZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[1],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[7],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[10],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[16],iVolNum++,3,dpar);

// horizontal cable support chamber type E

	    spar[0]=dpar[0]-(kYVSup[0]/2.)*zRatio;
	    spar[1]=kSizeSupportCable[1];
	    spar[2]=kSizeSupportCable[2];
	    for(Int_t in=0;in<3;in++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio,-(yEnvM+dpar[1]/2.*(in-1)),kSizeVSupExt[0]+spar[2]));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+3], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio, yEnvM+dpar[1]/2.*(in-1),kSizeVSupExt[0]+spar[2]));
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+6], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio),-(yEnvP+dpar[1]/2.*(in-1)),-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+9], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio), yEnvP+dpar[1]/2.*(in-1),-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    }
	    for(Int_t ii=0;ii<12;ii++) 
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHE[ii],iVolNum++,3, spar);
	    
// angular supports chamber type E
// 1 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1],kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM-dpar[1],kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvP-dpar[1],kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvP-dpar[1],kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VE[i],iVolNum++,3, spar);


// 1 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HE[i],iVolNum++,3, spar);

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1E[i],iVolNum++,3, ppar);

// 2 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VE[i],iVolNum++,3, spar);	    

// 2 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HE[i],iVolNum++,3, spar);	    

// 3 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM+dpar[1],kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1],kZm-dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvP+dpar[1],kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvP+dpar[1],kZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VE[i],iVolNum++,3, spar);


// 3 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HE[i],iVolNum++,3, spar);

// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2E[i],iVolNum++,3, ppar);

// 4 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VE[i],iVolNum++,3, spar);	    

// 4 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HE[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HE[1], 0, true,
	    TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HE[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HE[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HE[i],iVolNum++,3, spar);	    

// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHE[i],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVE[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVE[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVE[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVE[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVE[i],iVolNum++,3, spar);	    

// F
	    yEnvPsave = yEnvP;
	    yEnvMsave = yEnvM;
	    yEnvP = (yEnvMsave + kYMIN * zRatio ) * zpm + kYMIN * zRatio;
	    yEnvM = (yEnvPsave + kYMIN * zRatio ) * zmp + kYMIN * zRatio;

	    detElemId = (10+icount+1)*100+14;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[0], detElemId, true, TGeoTranslation(xEnv,-yEnvP,kZp));
	    detElemId = (10+icount+1)*100+4;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[8], detElemId, true, TGeoTranslation(xEnv, yEnvP,kZp));
	    detElemId = (10+icount+1)*100+13;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[9], detElemId, true, TGeoTranslation(-xEnv,-yEnvM,kZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    detElemId = (10+icount+1)*100+5;
	    GetEnvelopes(16+icount)->AddEnvelope(volEnv[17], detElemId, true, TGeoTranslation(-xEnv, yEnvM,kZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[0],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[8],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[9],iVolNum++,3,dpar);
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv[17],iVolNum++,3,dpar);

// horizontal cable support chamber type F

	    spar[0]=dpar[0]-(kYVSup[0]/2.)*zRatio;
	    spar[1]=kSizeSupportCable[1];
	    spar[2]=kSizeSupportCable[2];
	    for(Int_t in=0;in<3;in++){
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHF[in], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio,-(yEnvP+dpar[1]/2.*(in-1)),-(kSizeVSupExt[0]+spar[2])));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHF[in+3], 0, true, TGeoTranslation(xEnv+kYVSup[0]/2.*zRatio,yEnvP+dpar[1]/2.*(in-1),-(kSizeVSupExt[0]+spar[2])));
       	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHF[in+6], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio),-(yEnvM+dpar[1]/2.*(in-1)),kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	      GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHF[in+9], 0, true, TGeoTranslation(-(xEnv+kYVSup[0]/2.*zRatio), yEnvM+dpar[1]/2.*(in-1),kSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    }
	    for(Int_t ii=0;ii<12;ii++) 
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHF[ii],iVolNum++,3, spar);

// angular supports chamber type F
// 1 vertical	   
 
 	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[0], 0, true,
 	    TGeoTranslation(xEnv,-yEnvP-dpar[1],kZp+dpar[2]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[1], 0, true,
 	    TGeoTranslation(xEnv,yEnvP-dpar[1],kZp+dpar[2]+spar[2]));
 	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[2], 0, true,
  	    TGeoTranslation(-xEnv,-yEnvM-dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[3], 0, true,
  	    TGeoTranslation(-xEnv,yEnvM-dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
            for(i=0;i<4;i++)
 	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VF[i],iVolNum++,3, spar);

// 1 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HF[i],iVolNum++,3, spar);

// gas pipe (low)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-kOffsetGasPipe,kZm),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1F[i],iVolNum++,3, ppar);

// 2 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VF[i],iVolNum++,3, spar);	    

// 2 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HF[i],iVolNum++,3, spar);	    

// 3 vertical	   
 
 	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1V[1];
	    spar[2]=kSizeSupport1V[2];
	    sparysave=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[0], 0, true,
 	    TGeoTranslation(xEnv,-yEnvP+dpar[1],kZp+dpar[2]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[1], 0, true,
 	    TGeoTranslation(xEnv,yEnvP+dpar[1],kZp+dpar[2]+spar[2]));
 	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[2], 0, true,
  	    TGeoTranslation(-xEnv,-yEnvM+dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[3], 0, true,
  	    TGeoTranslation(-xEnv,yEnvM+dpar[1],kZm-dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
            for(i=0;i<4;i++)
 	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VF[i],iVolNum++,3, spar);

// 3 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport1H[1];
	    spar[2]=kSizeSupport1H[2];
      	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kZp+dpar[2]-spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,kZm-(dpar[2]-spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HF[i],iVolNum++,3, spar);
	    
// gas pipe (high)
	    ppar[0]=kSizeGasPipe[0];
	    ppar[1]=kSizeGasPipe[1];
	    ppar[2]=dpar[0];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+kOffsetGasPipe,kZp),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[2], 0, true,
	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[3], 0, true,
	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+kOffsetGasPipe,kZm),rsupportpipe);

            for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2F[i],iVolNum++,3, ppar);

// 4 vertical	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2V[1];
	    spar[2]=kSizeSupport2V[2];
	    sparysave=spar[1]+kSizeSupport2H[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    
	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VF[i],iVolNum++,3, spar);	    

// 4 horizontal	   
 
	    spar[0]=dpar[0];
	    spar[1]=kSizeSupport2H[1];
	    spar[2]=kSizeSupport2H[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[0], 0, true,
	    TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[1], 0, true,
	    TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[2], 0, true,
 	    TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[3], 0, true,
 	    TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HF[i],iVolNum++,3, spar);	    


// X horizontal	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXH[1];
	    spar[2]=kSizeSupportXH[2];
	    sparysavex=spar[1];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHF[i],iVolNum++,3, spar);	    

// X vertical	   
	    spar[0]=(kYVSup[3]-kYVSup[0])/2.*zRatio;
	    spar[1]=kSizeSupportXV[1];
	    spar[2]=kSizeSupportXV[2];
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[0], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]));
       	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[1], 0, true,
	    TGeoTranslation(spar[0]+kYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],kSizeVSupExt[0]+spar[2]));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[2], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
	    GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[3], 0, true,
 	    TGeoTranslation(-(spar[0]+kYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],-(kSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

	    for(i=0;i<4;i++)
	    GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVF[i],iVolNum++,3, spar);
	    	    
	} // end loop on detection planes
    } // end loop on stations    
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::SetVolumes() 
{
/// Defines the volumes for the trigger chambers.

    if (gAlice->GetModule("SHIL")) {
      SetMotherVolume(16, "YOUT2");
      SetMotherVolume(17, "YOUT2");
      SetMotherVolume(18, "YOUT2");
      SetMotherVolume(19, "YOUT2");
    }  

    SetVolume(16, "SC11");
    SetVolume(17, "SC12");
    SetVolume(18, "SC13");
    SetVolume(19, "SC14");
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::SetTransformations() 
{
/// Defines the transformations for the trigger chambers.

    Double_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    SetTranslation(16, TGeoTranslation(0., 0., zpos1));
    
    zpos1= AliMUONConstants::DefaultChamberZ(11); 
    SetTranslation(17, TGeoTranslation(0., 0., zpos1));

    zpos1= AliMUONConstants::DefaultChamberZ(12); 
    SetTranslation(18, TGeoTranslation(0., 0., zpos1));

    zpos1= AliMUONConstants::DefaultChamberZ(13); 
    SetTranslation(19, TGeoTranslation(0., 0., zpos1));
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::SetSensitiveVolumes()
{
/// Defines the sensitive volumes for trigger station chambers.

  GetGeometry(16)->SetSensitiveVolume("S11G");
  GetGeometry(17)->SetSensitiveVolume("S12G");
  GetGeometry(18)->SetSensitiveVolume("S13G");
  GetGeometry(19)->SetSensitiveVolume("S14G");
}

