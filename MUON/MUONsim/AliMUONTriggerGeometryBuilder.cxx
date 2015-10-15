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
#include <TArrayI.h>

#include "AliLog.h"
#include "AliRun.h"

#include "AliMUONTriggerGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include <iostream>

using std::endl;
using std::cout;
/// \cond CLASSIMP
ClassImp(AliMUONTriggerGeometryBuilder)
/// \endcond

// vertical gap between right and left chambers (kDXZERO*2=4cm)
const Float_t AliMUONTriggerGeometryBuilder::fgkDXZERO=2.; 
// main distances for chamber definition in first plane/first station
const Float_t AliMUONTriggerGeometryBuilder::fgkXMIN=34.;       
const Float_t AliMUONTriggerGeometryBuilder::fgkXMED=51.;                                
const Float_t AliMUONTriggerGeometryBuilder::fgkXMAX=255.; 
// 090704 fgkXMAX changed from 272 to 255.
// (see fig.2-4 & 2-5 of Local Trigger Board PRR)
// segmentation updated accordingly
const Float_t AliMUONTriggerGeometryBuilder::fgkYMIN=34.;                              
const Float_t AliMUONTriggerGeometryBuilder::fgkYMAX=51.;                              
// inner/outer radius of flange between beam shield. and chambers (1/station)
//const Float_t AliMUONTriggerGeometryBuilder::fgkRMIN[2]={50.,50.};
//const Float_t AliMUONTriggerGeometryBuilder::fgkRMAX[2]={64.,68.};
// z position of the middle of the gas gap in mother vol 
const Float_t AliMUONTriggerGeometryBuilder::fgkZm=-3.6;
const Float_t AliMUONTriggerGeometryBuilder::fgkZp=+3.6;
    
// y positions of vertical supports
const Float_t AliMUONTriggerGeometryBuilder::fgkYVSup[4]={61.45,122.45,192.95,236.95}; 
// dimensions of vertical supports 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeVSupExt[3]={1.5,1.5,306.+5.}; 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeVSupInt[3]={1.2,1.2,306.+5.};  
// transverse dimensions of angular supports 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupport1V[3]={0.,1.5,0.1}; 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupport1H[3]={0.,0.1,1.15}; // z should be 1.4 in the installed set-up 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupport2V[3]={0.,3.0,0.1}; 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupport2H[3]={0.,0.1,1.9}; 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupportXV[3]={0.,1.25,0.25}; 
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupportXH[3]={0.,0.25,1.5}; 
// transverse dimensions of horizontal cable supports
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeSupportCable[3]={0.,2.,3.}; 
// dimensions of gas pipes (inner and outer radius)
const Float_t AliMUONTriggerGeometryBuilder::fgkSizeGasPipe[3]={0.2,0.4,0.}; 
// Position of gas pipe with respect to angular support
const Float_t AliMUONTriggerGeometryBuilder::fgkOffsetGasPipe=0.75; 
// Small cut on some volumes to avoid extrusion from SC1x
const Float_t AliMUONTriggerGeometryBuilder::fgkAvoidExtrusion=2.9;    

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(16, 4),
   fMUON(muon),
   fIdtmed(0),
   fIdAir(0),  
   fIdAlu1(0),
   fIdInox(0),
   fYEnvPsave(0.), 
   fYEnvMsave(0.), 
   fDYsave(0.),
   fDXsave(0.),
   fRsupportpipe()
{
/// Standard constructor
   fRsupportpipe.SetAngles(90.,90.,0.);
}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::AliMUONTriggerGeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0),
   fIdtmed(0),
   fIdAir(0),  
   fIdAlu1(0),
   fIdInox(0),
   fYEnvPsave(0.), 
   fYEnvMsave(0.), 
   fDYsave(0.),
   fDXsave(0.),
   fRsupportpipe()
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONTriggerGeometryBuilder::~AliMUONTriggerGeometryBuilder() 
{
/// Destructor
}

//______________________________________________________________________________
TString AliMUONTriggerGeometryBuilder::GetVolumeName(const TString& volume, Int_t icount) const
{
// Function to generate a short volume name from its long variant

  if ( volume == "volAlu" ) {
    char volAlu[5];     // Alu 
    snprintf(volAlu,5,"SC%dA",icount+1);
    return volAlu;
  }
  else if ( volume == "volAluSupportH") {
    char volAluSupportH[6];
    snprintf(volAluSupportH,6,"SALH%d",icount+1);
    return  volAluSupportH;
  }         
  else if ( volume == "volAirSupportH") {
    char volAirSupportH[6];
    snprintf(volAirSupportH,6,"SAIH%d",icount+1);
    return  volAirSupportH;
  }
  else if ( volume == "volInoxGasPipe") {
    char volInoxGasPipe[7];
    snprintf(volInoxGasPipe,7,"SPINO%d",icount+1);
    return volInoxGasPipe;
  }
  
  AliErrorStream() << "Volume " << volume << " name is not defined." << endl; 
  return "";
}

//______________________________________________________________________________
TString AliMUONTriggerGeometryBuilder::GetVolEnvName(Int_t icount, Int_t ienv) const
{
/// Compose envelope names as:
/// S0R1, S0R2, ..., S0R9, S0L1, S0L2, ..., S0L9
/// where ienv = 0, .., 17

  TString name = "S";
  name += icount;
  if ( ienv < 9 ) {
    name += "R";
    name += (ienv + 1);
  }  
  else {
    name += "L";
    name += (ienv - 8) ;
  }  
  return name;
}     

//______________________________________________________________________________
TString AliMUONTriggerGeometryBuilder::GetVolAluAngSuppName(
                                          const TString& type1234X, 
                                          const TString& typeHV,
                                          Int_t icount) const
{
/// Utility function to generate volume name 

  TString name = "SA";
  name += type1234X;
  name += typeHV;
  name += icount+1;
  return name;
}  
   
//______________________________________________________________________________
TString AliMUONTriggerGeometryBuilder::GetVolEnvSuppAngName(
                                          const TString& type1234X, 
                                          const TString& typeHV, 
                                          const TString& typeABDEF,
                                          Int_t icount, Int_t ivol) const
{
/// Utility function to generate volume name 

  TString name = "S";
  name += typeHV;
  name += type1234X;
  name += icount+1;
  name += typeABDEF;
  name += ivol;
  return name;
}  

//______________________________________________________________________________
TString AliMUONTriggerGeometryBuilder::GetVolEnvInoxGasPipeName(
                                          const TString& type12, 
                                          const TString& typeABCDEF,
                                          Int_t icount, Int_t ivol) const
{
/// Utility function to generate volume name 

  TString name = "SP";
  name += type12;
  name += icount+1;
  name += typeABCDEF;
  name += ivol;
  return name;
}  

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberPrototype(Int_t icount) const
{   
/// Build chamber prototype

    Float_t tpar[3];
    tpar[0]= 0.;
    tpar[1]= 0.;
    tpar[2]= 0.;	    
    char volBak[5];     // Bakelite
    char volGaz[5];     // Gas streamer	    
    snprintf(volBak,5,"SB%dA",icount+1);
    snprintf(volGaz,5,"S%dG",icount+11);
    TVirtualMC::GetMC()->Gsvolu(GetVolumeName("volAlu", icount),"BOX",fIdAlu1,tpar,0);         // Al
    TVirtualMC::GetMC()->Gsvolu(volBak,"BOX",fIdtmed[1107],tpar,0);   // Bakelite
    TVirtualMC::GetMC()->Gsvolu(volGaz,"BOX",fIdtmed[1106],tpar,0);   // Gas streamer
    tpar[0] = -1.;
    tpar[1] = -1.;
    tpar[2] = 0.1;    
    TVirtualMC::GetMC()->Gsposp(volGaz,1,volBak,0.,0.,0.,0,"ONLY",tpar,3);
    tpar[2] = 0.3;
    TVirtualMC::GetMC()->Gsposp(volBak,1,GetVolumeName("volAlu", icount),0.,0.,0.,0,"ONLY",tpar,3);
}    

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildRPCSupportsVertical(Int_t& iVolNum, Int_t icount) const
{
/// Build RPC vertical supports

   Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);	    	    
   char volAluSupport[5],volAirSupport[5];
   snprintf(volAluSupport,5,"SAL%d",icount+1);
   snprintf(volAirSupport,5,"SAI%d",icount+1);
   char volEnvSupport[12][7];
   for(Int_t ii=0;ii<8;ii++){
     snprintf(volEnvSupport[ii],7,"SEA%dV%d",icount+1,ii);
   }
   Float_t tpar[3];
   tpar[0]= 0.;
   tpar[1]= 0.;
   tpar[2]= 0.;            
   TVirtualMC::GetMC()->Gsvolu(volAluSupport,"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(volAirSupport,"BOX",fIdAir,tpar,0);
   tpar[0]=fgkSizeVSupInt[0];
   tpar[1]=fgkSizeVSupInt[1];
   tpar[2]=-1.;
   TVirtualMC::GetMC()->Gsposp(volAirSupport,1,volAluSupport,0.,0.,0.,0,"ONLY",tpar,3);
   
   TGeoRotation rsupportv;
   rsupportv.SetAngles(0.,90.,0.);
   Double_t dpar[3];
   dpar[0]=fgkSizeVSupExt[0];
   dpar[1]=fgkSizeVSupExt[1];
   dpar[2]=fgkSizeVSupExt[2]*zRatio;
   for(Int_t ii=0;ii<4;ii++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupport[ii], 0, true,
     TGeoTranslation(-fgkYVSup[ii]*zRatio,0.,0.),rsupportv);
     GetEnvelopes(16+icount)
     ->AddEnvelopeConstituentParam(volAluSupport,volEnvSupport[ii],iVolNum++,3, dpar);
   }       
   for(Int_t ii=4;ii<8;ii++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupport[ii], 0, true,
     TGeoTranslation(fgkYVSup[ii-4]*zRatio,0.,0.),rsupportv);
     GetEnvelopes(16+icount)
     ->AddEnvelopeConstituentParam(volAluSupport,volEnvSupport[ii],iVolNum++,3, dpar);
   }   
}	        

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildRPCSupportsHorizontal(Int_t icount) const
{ 
/// Build RPC horizontal supports   

// RPC supports (horizontal)

// supports for cables

   Float_t tpar[3];
   tpar[0]= 0.;
   tpar[1]= 0.;
   tpar[2]= 0.;  
   TString volAluSupportH = GetVolumeName("volAluSupportH", icount);        
   TString volAirSupportH = GetVolumeName("volAirSupportH", icount);        
   TVirtualMC::GetMC()->Gsvolu(volAluSupportH,"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(volAirSupportH,"BOX",fIdAir,tpar,0);
   tpar[0]=-1.;
   tpar[1]=1.9;
   tpar[2]=2.8;
   TVirtualMC::GetMC()->Gsposp(volAirSupportH,1,volAluSupportH,0.,0.,0.,0,"ONLY",tpar,3);
}	    

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildAngularSupportForChambers(Int_t icount) const
{
/// Build angular supports for chambers

   Float_t tpar[3];
   tpar[0]= 0.;
   tpar[1]= 0.;
   tpar[2]= 0.;            
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("1","V",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("1","H",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("2","V",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("2","H",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("3","V",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("3","H",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("4","V",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("4","H",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("X","V",icount),"BOX",fIdAlu1,tpar,0);
   TVirtualMC::GetMC()->Gsvolu(GetVolAluAngSuppName("X","H",icount),"BOX",fIdAlu1,tpar,0);
}	    

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildGasPipes(Int_t icount) const
{        
/// Build gas pipes
   TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   Float_t tpar[3];
   tpar[0]= 0.;
   tpar[1]= 0.;
   tpar[2]= 0.;            
   TVirtualMC::GetMC()->Gsvolu(volInoxGasPipe,"TUBE",fIdInox,tpar,0);
}	   
  
//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberTypeA(Int_t& iVolNum, Int_t icount) 
{
/// Build chamber type A and horizontal support 

   Double_t dpar[3];    
   Double_t spar[3];    
   Double_t ppar[3];    

// chamber type A 
   Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);	    	    
   Float_t xEnv = (fgkDXZERO+fgkXMED+(fgkXMAX-fgkXMED)/2.)*zRatio;
   Float_t yEnvM = 0.;	 // y low position of envelope in chamber
   Float_t yEnvP = 0.;	 // y up position of envelope in chamber
   fYEnvPsave = 0.; // tmp data
   fYEnvMsave = 0.; // tmp data
   //Float_t xpos = 0.; // x position of RPC in envelope	    
   //Float_t ypos = 0.; // y position of RPC in envelope
   dpar[2] = 0.4;	    
   dpar[0] = ((fgkXMAX-fgkXMED)/2.)*zRatio;
   dpar[1] = fgkYMIN * zRatio;

   Int_t detElemId = (10+icount+1)*100;
   TString volEnv4 = GetVolEnvName(icount, 4);
   TString volEnv13 = GetVolEnvName(icount, 13);
   GetEnvelopes(16+icount)->AddEnvelope(volEnv4, detElemId, true,  
   TGeoTranslation(xEnv,yEnvP,fgkZp));
   detElemId = (10+icount+1)*100+9;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv13, detElemId, true,    
   TGeoTranslation(-xEnv,yEnvM,fgkZm),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAlu = GetVolumeName("volAlu", icount);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv4,iVolNum++,3, dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv13,iVolNum++,3, dpar);	    

// horizontal cable support chamber type A
   char volEnvSupportHA[6][8];
   for(Int_t jj=0;jj<2;jj++){
     for(Int_t ii=0;ii<6;ii++){
       if(ii<3)snprintf(volEnvSupportHA[3*jj+ii],8,"SA%dHA%d",icount+1,3*jj+ii);
     }
   }

   spar[0]=((fgkXMAX/2)-fgkYVSup[0]/2.)*zRatio;
   spar[1]=fgkSizeSupportCable[1];
   spar[2]=fgkSizeSupportCable[2];
   Float_t offsetSuppA = ((fgkXMAX-fgkXMED)/2.)*zRatio-(((fgkXMAX/2)-fgkYVSup[0]/2.)*zRatio);
   for(Int_t in=0;in<3;in++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHA[in], 0, true, 
     TGeoTranslation(xEnv+offsetSuppA/2.,yEnvP+dpar[1]/2.*(in-1),-(fgkSizeVSupExt[0]+spar[2])));    
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHA[in+3], 0, true, 
     TGeoTranslation(-(xEnv+offsetSuppA/2.),yEnvM+dpar[1]/2.*(in-1),fgkSizeVSupExt[0]+spar[2]),
     TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   }
   for(Int_t ii=0;ii<6;ii++) {
     TString volAluSupportH = GetVolumeName("volAluSupportH", icount);
     GetEnvelopes(16+icount)
     ->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHA[ii],iVolNum++,3, spar);
   }  

   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   Float_t sparysave=spar[1];
   TString volEnvSuppAng1VA0 = GetVolEnvSuppAngName("1", "V", "A", icount, 0);
   TString volEnvSuppAng1VA1 = GetVolEnvSuppAngName("1", "V", "A", icount, 1);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VA0, 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VA1, 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VA1,iVolNum++,3, spar);	    

// 1 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   TString volEnvSuppAng1HA0 = GetVolEnvSuppAngName("1", "H", "A", icount, 0);
   TString volEnvSuppAng1HA1 = GetVolEnvSuppAngName("1", "H", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HA0, 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HA1, 0, true, 
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HA1,iVolNum++,3, spar);	    

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe1A0 = GetVolEnvInoxGasPipeName("1", "A", icount, 0);
   TString volEnvInoxGasPipe1A1 = GetVolEnvInoxGasPipeName("1", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1A0, 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1A1, 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1A0,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1A1,iVolNum++,3, ppar);

// 2 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng2VA0 = GetVolEnvSuppAngName("2", "V", "A", icount, 0);
   TString volEnvSuppAng2VA1 = GetVolEnvSuppAngName("2", "V", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VA0, 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VA1, 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VA1,iVolNum++,3, spar);	    

// 2 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2]; 
   TString volEnvSuppAng2HA0 = GetVolEnvSuppAngName("2", "H", "A", icount, 0);
   TString volEnvSuppAng2HA1 = GetVolEnvSuppAngName("2", "H", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HA0, 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HA1, 0, true,   
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HA1,iVolNum++,3, spar);	    

// 3 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng3VA0 = GetVolEnvSuppAngName("3", "V", "A", icount, 0);
   TString volEnvSuppAng3VA1 = GetVolEnvSuppAngName("3", "V", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VA0, 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VA1, 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VA1,iVolNum++,3, spar);	    

// 3 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   TString volEnvSuppAng3HA0 = GetVolEnvSuppAngName("3", "H", "A", icount, 0);
   TString volEnvSuppAng3HA1 = GetVolEnvSuppAngName("3", "H", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HA0, 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HA1, 0, true,   
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3H = GetVolAluAngSuppName("3", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HA1,iVolNum++,3, spar);	
       
// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe2A0 = GetVolEnvInoxGasPipeName("2", "A", icount, 0);
   TString volEnvInoxGasPipe2A1 = GetVolEnvInoxGasPipeName("2", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2A0, 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2A1, 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2A0,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2A1,iVolNum++,3, ppar);

// 4 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng4VA0 = GetVolEnvSuppAngName("4", "V", "A", icount, 0);
   TString volEnvSuppAng4VA1 = GetVolEnvSuppAngName("4", "V", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VA0, 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VA1, 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4V = GetVolAluAngSuppName("4", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VA1,iVolNum++,3, spar);   

// 4 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng4HA0 = GetVolEnvSuppAngName("4", "H", "A", icount, 0);
   TString volEnvSuppAng4HA1 = GetVolEnvSuppAngName("4", "H", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HA0, 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HA1, 0, true, 
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4H = GetVolAluAngSuppName("4", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HA0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HA1,iVolNum++,3, spar);   

// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   Float_t sparysavex=spar[1];
   TString volEnvSuppAngXHA0 = GetVolEnvSuppAngName("X", "H", "A", icount, 0);
   TString volEnvSuppAngXHA1 = GetVolEnvSuppAngName("X", "H", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHA0, 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHA1, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXH = GetVolAluAngSuppName("X", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHA0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHA1,iVolNum++,3, spar);   

// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];
   TString volEnvSuppAngXVA0 = GetVolEnvSuppAngName("X", "V", "A", icount, 0);
   TString volEnvSuppAngXVA1 = GetVolEnvSuppAngName("X", "V", "A", icount, 1);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVA0, 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVA1, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXV = GetVolAluAngSuppName("X", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVA0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVA1,iVolNum++,3, spar);	   

   // keep values of yEnvP, yEnvM
   fYEnvPsave = yEnvP;
   fYEnvMsave = yEnvM;   
} 

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberTypeB(Int_t& iVolNum, Int_t icount)
{
// ratio of zpos1m/zpos1p and inverse for first plane
   Float_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
   Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
   Float_t zpm=1./zmp;
    
// chamber type B (plus envelope chambers B & C)   
   Double_t dpar[3];    
   Double_t spar[3];    
   Double_t ppar[3];    

// chamber type B
   Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);      
   Float_t xEnv = (fgkDXZERO+fgkXMAX/2.)*zRatio;
   Float_t yEnvP = 0;
   Float_t yEnvM = 0;
   yEnvP = (fYEnvMsave + fgkYMIN * zRatio ) * zpm + fgkYMIN * zRatio;
   yEnvM = (fYEnvPsave + fgkYMIN * zRatio ) * zmp + fgkYMIN * zRatio;
   dpar[0] = ((fgkXMAX-fgkXMIN)/2.) * zRatio;
   dpar[1] = ((fgkYMAX-fgkYMIN)/2.) * zRatio;
   dpar[2] = 0.4;   
   fDYsave = dpar[1];
   fDXsave = dpar[0];
   Float_t xpos = fgkXMIN/2. * zRatio;
   Float_t ypos = (fgkYMIN - fgkYMIN/4.) * zRatio;
   Float_t xpossave = xpos;

   Int_t detElemId = (10+icount+1)*100+17;            
   TString volEnv3 = GetVolEnvName(icount, 3);
   TString volEnv5 = GetVolEnvName(icount, 5);
   TString volEnv12 = GetVolEnvName(icount, 12);
   TString volEnv14 = GetVolEnvName(icount, 14);
   TString volAlu = GetVolumeName("volAlu", icount);

   GetEnvelopes(16+icount)->AddEnvelope(volEnv3, detElemId, true,
   TGeoTranslation(xEnv,-yEnvM,fgkZm));   
   detElemId = (10+icount+1)*100+1;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv5, detElemId, true, 
   TGeoTranslation( xEnv, yEnvM,fgkZm));
   detElemId = (10+icount+1)*100+10;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv12, detElemId, true,
   TGeoTranslation(-xEnv,-yEnvP,fgkZp),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   detElemId = (10+icount+1)*100+8;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv14, detElemId, true, 
   TGeoTranslation(-xEnv, yEnvP,fgkZp),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv3,iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv5,iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv12,iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv14,iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);

// chamber type C (note: same Z than type B)
   dpar[0] = (fgkXMAX/2)*zRatio;
   dpar[1] = (fgkYMAX/2)*zRatio;
   xpos = 0.;   
   ypos = ((fgkYMAX - fgkYMIN)/2.) * zRatio;

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv3,iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv5,iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv12,iVolNum++,TGeoTranslation(xpos,-ypos,0.),3,dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv14,iVolNum++,TGeoTranslation(xpos, ypos,0.),3,dpar);

// horizontal cable support chamber type B+C

   char volEnvSupportHBC[12][8];
   for(Int_t jj=0;jj<2;jj++){
     for(Int_t ii=0;ii<6;ii++){
       snprintf(volEnvSupportHBC[6*jj+ii],8,"SA%dHB%d",icount+1,6*jj+ii);
     }
   }
      
   spar[0]=dpar[0]-fgkYVSup[0]/2.;
   spar[1]=fgkSizeSupportCable[1];
   spar[2]=fgkSizeSupportCable[2];
   for(Int_t in=0;in<3;in++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in], 0, true, 
     TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio,-(yEnvM+(fgkYMAX-fgkYMIN/2.*zRatio)/2.*(in-1)),
     fgkSizeVSupExt[0]+spar[2]));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+3], 0, true, 
     TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio, yEnvM+(fgkYMAX-fgkYMIN/2.*zRatio)/2.*(in-1),
     fgkSizeVSupExt[0]+spar[2]));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+6], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio),-(yEnvP+(fgkYMAX-fgkYMIN/2.*zRatio)/2.*(in-1)),
     -(fgkSizeVSupExt[0]+spar[2])),
     TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHBC[in+9], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio), yEnvP+(fgkYMAX-fgkYMIN/2.*zRatio)/2.*(in-1),
     -(fgkSizeVSupExt[0]+spar[2])),
     TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   }
   for(Int_t ii=0;ii<12;ii++) { 
     TString volAluSupportH = GetVolumeName("volAluSupportH", icount);
     GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHBC[ii],iVolNum++,3, spar);
   }  

// angular supports chamber type B and C
// C	   
// 1 vertical
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   Float_t sparysave=spar[1];

   TString volEnvSuppAng1VBC0 = GetVolEnvSuppAngName("1", "V", "B", icount, 0); 
   TString volEnvSuppAng1VBC2 = GetVolEnvSuppAngName("1", "V", "B", icount, 2); 
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC0, 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-fDYsave,fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC2, 0, true,
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-fDYsave,fgkZp+dpar[2]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC2,iVolNum++,3, spar);

// 1 horizontal
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   
   TString volEnvSuppAng1HBC0 = GetVolEnvSuppAngName("1", "H", "B", icount, 0);
   TString volEnvSuppAng1HBC2 = GetVolEnvSuppAngName("1", "H", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC0, 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-fDYsave-sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC2, 0, true, 
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-fDYsave-sparysave,fgkZp+dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC2,iVolNum++,3, spar);   

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe1BC0 = GetVolEnvInoxGasPipeName("1", "BC", icount, 0);
   TString volEnvInoxGasPipe1BC2 = GetVolEnvInoxGasPipeName("1", "BC", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC0, 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-fDYsave-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC2, 0, true,
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-fDYsave-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC0,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC2,iVolNum++,3, ppar);

// 2 vertical	   
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng2VBC0 = GetVolEnvSuppAngName("2", "V", "B", icount, 0);
   TString volEnvSuppAng2VBC2 = GetVolEnvSuppAngName("2", "V", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC0, 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-fDYsave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC2, 0, true, 
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-fDYsave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC2,iVolNum++,3, spar);   

// 2 horizontal	   
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng2HBC0 = GetVolEnvSuppAngName("2", "H", "B", icount, 0);
   TString volEnvSuppAng2HBC2 = GetVolEnvSuppAngName("2", "H", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC0, 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-fDYsave-sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC2, 0, true, 
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-fDYsave-sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC2,iVolNum++,3, spar);   

// 3 vertical
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng3VBC0 = GetVolEnvSuppAngName("3", "V", "B", icount, 0);
   TString volEnvSuppAng3VBC2 = GetVolEnvSuppAngName("3", "V", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC0, 0, true,
   TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+fDYsave,fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC2, 0, true, 
   TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+fDYsave,fgkZp+dpar[2]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC2,iVolNum++,3, spar);

// 3 horizontal
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   TString volEnvSuppAng3HBC0 = GetVolEnvSuppAngName("3", "H", "B", icount, 0);
   TString volEnvSuppAng3HBC2 = GetVolEnvSuppAngName("3", "H", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC0, 0, true,
   TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+fDYsave+sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC2, 0, true, 
   TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+fDYsave+sparysave,fgkZp+dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3H = GetVolAluAngSuppName("3", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC0,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC2,iVolNum++,3, spar);   
       
// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=fDXsave-fgkAvoidExtrusion;
   TString volEnvInoxGasPipe2BC0 = GetVolEnvInoxGasPipeName("2", "BC", icount, 0);
   TString volEnvInoxGasPipe2BC2 = GetVolEnvInoxGasPipeName("2", "BC", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC0, 0, true,
   TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+fDYsave+sparysave+fgkOffsetGasPipe,fgkZm),
   fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC2, 0, true,
   TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+fDYsave+sparysave+fgkOffsetGasPipe,fgkZp),
   fRsupportpipe);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC0,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC2,iVolNum++,3, ppar);

// 4 vertical	   
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng4VBC0 = GetVolEnvSuppAngName("4", "V", "B", icount, 0);
   TString volEnvSuppAng4VBC2 = GetVolEnvSuppAngName("4", "V", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC0, 0, true,
   TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+fDYsave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC2, 0, true, 
   TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+fDYsave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport4V = GetVolAluAngSuppName("4", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC2,iVolNum++,3, spar);   

// 4 horizontal	   
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng4HBC0 = GetVolEnvSuppAngName("4", "H", "B", icount, 0);
   TString volEnvSuppAng4HBC2 = GetVolEnvSuppAngName("4", "H", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC0, 0, true,
   TGeoTranslation(xEnv+xpossave,-yEnvM+dpar[1]+fDYsave+sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC2, 0, true, 
   TGeoTranslation(-xEnv-xpossave,-yEnvP+dpar[1]+fDYsave+sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4H = GetVolAluAngSuppName("4", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC2,iVolNum++,3, spar);   

// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   Float_t sparysavex=spar[1];
   TString volEnvSuppAngXHBC0 = GetVolEnvSuppAngName("X", "H", "B", icount, 0);
   TString volEnvSuppAngXHBC2 = GetVolEnvSuppAngName("X", "H", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC0, 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvM+dpar[1]+fDYsave+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC2, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvP+dpar[1]+fDYsave+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXH = GetVolAluAngSuppName("X", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC2,iVolNum++,3, spar);	
// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];
   TString volEnvSuppAngXVBC0 = GetVolEnvSuppAngName("X", "V", "B", icount, 0);
   TString volEnvSuppAngXVBC2 = GetVolEnvSuppAngName("X", "V", "B", icount, 2);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC0, 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,
   -yEnvM+dpar[1]+fDYsave+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC2, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),
   -yEnvP+dpar[1]+fDYsave+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXV = GetVolAluAngSuppName("X", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC0,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC2,iVolNum++,3, spar);   

// B
// 1 vertical
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng1VBC1 = GetVolEnvSuppAngName("1", "V", "B", icount, 1);
   TString volEnvSuppAng1VBC3 = GetVolEnvSuppAngName("1", "V", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC1, 0, true,
   TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-fDYsave,fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VBC3, 0, true, 
   TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-fDYsave,fgkZp+dpar[2]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   //TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC1,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VBC3,iVolNum++,3, spar);


// 1 horizontal
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   
   TString volEnvSuppAng1HBC1 = GetVolEnvSuppAngName("1", "H", "B", icount, 1);
   TString volEnvSuppAng1HBC3 = GetVolEnvSuppAngName("1", "H", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC1, 0, true,
   TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-fDYsave-sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HBC3, 0, true, 
   TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-fDYsave-sparysave,fgkZp+dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   //TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC1,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HBC3,iVolNum++,3, spar);   

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=fDXsave-fgkAvoidExtrusion;
   TString volEnvInoxGasPipe1BC1 = GetVolEnvInoxGasPipeName("1", "BC", icount, 1);
   TString volEnvInoxGasPipe1BC3 = GetVolEnvInoxGasPipeName("1", "BC", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC1, 0, true,
   TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-fDYsave-sparysave-fgkOffsetGasPipe,fgkZm),
   fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1BC3, 0, true,
   TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-fDYsave-sparysave-fgkOffsetGasPipe,fgkZp),
   fRsupportpipe);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC1,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1BC3,iVolNum++,3, ppar);

// 2 vertical	   
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng2VBC1 = GetVolEnvSuppAngName("2", "V", "B", icount, 1);
   TString volEnvSuppAng2VBC3 = GetVolEnvSuppAngName("2", "V", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC1, 0, true,
   TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-fDYsave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VBC3, 0, true, 
   TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-fDYsave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   //TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VBC3,iVolNum++,3, spar);   
// 2 horizontal	   
   spar[0]=fDXsave-fgkAvoidExtrusion;
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng2HBC1 = GetVolEnvSuppAngName("2", "H", "B", icount, 1);
   TString volEnvSuppAng2HBC3 = GetVolEnvSuppAngName("2", "H", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC1, 0, true,
   TGeoTranslation(xEnv+xpossave,yEnvM-dpar[1]-fDYsave-sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HBC3, 0, true, 
   TGeoTranslation(-xEnv-xpossave,yEnvP-dpar[1]-fDYsave-sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   //TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HBC3,iVolNum++,3, spar);   

// 3 vertical
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng3VBC1 = GetVolEnvSuppAngName("3", "V", "B", icount, 1);
   TString volEnvSuppAng3VBC3 = GetVolEnvSuppAngName("3", "V", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC1, 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+fDYsave,fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VBC3, 0, true,
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+fDYsave,fgkZp+dpar[2]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   //TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC1,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VBC3,iVolNum++,3, spar);

// 3 horizontal
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   
   TString volEnvSuppAng3HBC1 = GetVolEnvSuppAngName("3", "H", "B", icount, 1);
   TString volEnvSuppAng3HBC3 = GetVolEnvSuppAngName("3", "H", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC1, 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+fDYsave+sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HBC3, 0, true, 
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+fDYsave+sparysave,fgkZp+dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC1,iVolNum++,3, spar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HBC3,iVolNum++,3, spar);   

// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe2BC1 = GetVolEnvInoxGasPipeName("2", "BC", icount, 1);
   TString volEnvInoxGasPipe2BC3 = GetVolEnvInoxGasPipeName("2", "BC", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC1, 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+fDYsave+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2BC3, 0, true,
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+fDYsave+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC1,iVolNum++,3, ppar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2BC3,iVolNum++,3, ppar);

// 4 vertical	   
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng4VBC1 = GetVolEnvSuppAngName("4", "V", "B", icount, 1);
   TString volEnvSuppAng4VBC3 = GetVolEnvSuppAngName("4", "V", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC1, 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+fDYsave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VBC3, 0, true, 
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+fDYsave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VBC3,iVolNum++,3, spar);   

// 4 horizontal	   
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng4HBC1 = GetVolEnvSuppAngName("4", "H", "B", icount, 1);
   TString volEnvSuppAng4HBC3 = GetVolEnvSuppAngName("4", "H", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC1, 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+fDYsave+sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HBC3, 0, true, 
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+fDYsave+sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HBC3,iVolNum++,3, spar);   

// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   sparysavex=spar[1];
   TString volEnvSuppAngXHBC1 = GetVolEnvSuppAngName("X", "H", "B", icount, 1);
   TString volEnvSuppAngXHBC3 = GetVolEnvSuppAngName("X", "H", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC1, 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvM+dpar[1]+fDYsave+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHBC3, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvP+dpar[1]+fDYsave+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHBC3,iVolNum++,3, spar);   

// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];
   TString volEnvSuppAngXVBC1 = GetVolEnvSuppAngName("X", "V", "B", icount, 1);
   TString volEnvSuppAngXVBC3 = GetVolEnvSuppAngName("X", "V", "B", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC1, 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,
   yEnvM+dpar[1]+fDYsave+sparysave+1.0+sparysavex+spar[1],-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVBC3, 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),
   yEnvP+dpar[1]+fDYsave+sparysave+1.0+sparysavex+spar[1],fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC1,iVolNum++,3, spar);   
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVBC3,iVolNum++,3, spar);   

   // keep values of yEnvP, yEnvM
   fYEnvPsave = yEnvP;
   fYEnvMsave = yEnvM;     
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberTypeD(Int_t& iVolNum, Int_t icount)
{
// ratio of zpos1m/zpos1p and inverse for first plane
    Float_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
    Float_t zpm=1./zmp;
    Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);      
    Float_t xEnv = (fgkDXZERO+fgkXMAX/2.)*zRatio;

   Double_t dpar[3];    
   Double_t spar[3];    
   Double_t ppar[3];    

// D   
   Float_t yEnvP = 0;
   Float_t yEnvM = 0;
   yEnvP = (fYEnvMsave + fgkYMIN * zRatio ) * zpm + fgkYMIN * zRatio;
   yEnvM = (fYEnvPsave + fgkYMIN * zRatio ) * zmp + fgkYMIN * zRatio;
   dpar[0] = (fgkXMAX/2.)*zRatio;
   dpar[1] =  fgkYMIN*zRatio;
   dpar[2] = 0.4;   

   Int_t detElemId = (10+icount+1)*100+16;   
   TString volEnv2 = GetVolEnvName(icount, 2);
   TString volEnv6 = GetVolEnvName(icount, 6);
   TString volEnv11 = GetVolEnvName(icount, 11);
   TString volEnv15 = GetVolEnvName(icount, 15);
   TString volAlu = GetVolumeName("volAlu", icount);

   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv2, detElemId, true, TGeoTranslation(xEnv,-yEnvP,fgkZp));
   detElemId = (10+icount+1)*100+2;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv6, detElemId, true, TGeoTranslation(xEnv, yEnvP,fgkZp));
   detElemId = (10+icount+1)*100+11;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv11, detElemId, true, TGeoTranslation(-xEnv,-yEnvM,fgkZm),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   detElemId = (10+icount+1)*100+7;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv15, detElemId, true, TGeoTranslation(-xEnv, yEnvM,fgkZm),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv2,iVolNum++,3, dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv6,iVolNum++,3, dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv11,iVolNum++,3, dpar);
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAlu,volEnv15,iVolNum++,3, dpar);

// horizontal cable support chamber type D

   char volEnvSupportHD[12][8];
   for(Int_t jj=0;jj<2;jj++){
     for(Int_t ii=0;ii<6;ii++){
       snprintf(volEnvSupportHD[6*jj+ii],8,"SA%dHD%d",icount+1,6*jj+ii);
     }
   }
      
   spar[0]=dpar[0]-(fgkYVSup[0]/2.)*zRatio;
   spar[1]=fgkSizeSupportCable[1];
   spar[2]=fgkSizeSupportCable[2];
   for(Int_t in=0;in<3;in++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in], 0, true,
      TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio,-(yEnvP+dpar[1]/2.*(in-1)),
      -(fgkSizeVSupExt[0]+spar[2])));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+3], 0, true, 
     TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio, yEnvP+dpar[1]/2.*(in-1),
     -(fgkSizeVSupExt[0]+spar[2])));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+6], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio),-(yEnvM+dpar[1]/2.*(in-1)),
     fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHD[in+9], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio),yEnvM+dpar[1]/2.*(in-1),
     fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   }
   for(Int_t ii=0;ii<12;ii++) { 
     TString volAluSupportH = GetVolumeName("volAluSupportH", icount);
     GetEnvelopes(16+icount)
     ->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHD[ii],iVolNum++,3, spar);
   }  
   
// angular supports chamber type D
// 1 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   Double_t sparysave=spar[1];
   
   TString volEnvSuppAng1VD[4];
   volEnvSuppAng1VD[0] =  GetVolEnvSuppAngName("1", "V", "D", icount, 0);
   volEnvSuppAng1VD[1] =  GetVolEnvSuppAngName("1", "V", "D", icount, 1);
   volEnvSuppAng1VD[2] =  GetVolEnvSuppAngName("1", "V", "D", icount, 2);
   volEnvSuppAng1VD[3] =  GetVolEnvSuppAngName("1", "V", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VD[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

       
   TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
   for (Int_t i=0;i<4;i++) 
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VD[i],iVolNum++,3, spar);


// 1 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];

   TString volEnvSuppAng1HD[4];
   volEnvSuppAng1HD[0] =  GetVolEnvSuppAngName("1", "H", "D", icount, 0);
   volEnvSuppAng1HD[1] =  GetVolEnvSuppAngName("1", "H", "D", icount, 1);
   volEnvSuppAng1HD[2] =  GetVolEnvSuppAngName("1", "H", "D", icount, 2);
   volEnvSuppAng1HD[3] =  GetVolEnvSuppAngName("1", "H", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HD[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HD[i],iVolNum++,3, spar);

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe1D[4];
   volEnvInoxGasPipe1D[0] = GetVolEnvInoxGasPipeName("1", "D", icount, 0);
   volEnvInoxGasPipe1D[1] = GetVolEnvInoxGasPipeName("1", "D", icount, 1);
   volEnvInoxGasPipe1D[2] = GetVolEnvInoxGasPipeName("1", "D", icount, 2);
   volEnvInoxGasPipe1D[3] = GetVolEnvInoxGasPipeName("1", "D", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1D[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);

   TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1D[i],iVolNum++,3, ppar);

// 2 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng2VD[4];
   volEnvSuppAng2VD[0] =  GetVolEnvSuppAngName("2", "V", "D", icount, 0);
   volEnvSuppAng2VD[1] =  GetVolEnvSuppAngName("2", "V", "D", icount, 1);
   volEnvSuppAng2VD[2] =  GetVolEnvSuppAngName("2", "V", "D", icount, 2);
   volEnvSuppAng2VD[3] =  GetVolEnvSuppAngName("2", "V", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VD[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VD[i],iVolNum++,3, spar);   

// 2 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng2HD[4];
   volEnvSuppAng2HD[0] =  GetVolEnvSuppAngName("2", "H", "D", icount, 0);
   volEnvSuppAng2HD[1] =  GetVolEnvSuppAngName("2", "H", "D", icount, 1);
   volEnvSuppAng2HD[2] =  GetVolEnvSuppAngName("2", "H", "D", icount, 2);
   volEnvSuppAng2HD[3] =  GetVolEnvSuppAngName("2", "H", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HD[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HD[i],iVolNum++,3, spar);   

// 3 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng3VD[4];
   volEnvSuppAng3VD[0] =  GetVolEnvSuppAngName("3", "V", "D", icount, 0);
   volEnvSuppAng3VD[1] =  GetVolEnvSuppAngName("3", "V", "D", icount, 1);
   volEnvSuppAng3VD[2] =  GetVolEnvSuppAngName("3", "V", "D", icount, 2);
   volEnvSuppAng3VD[3] =  GetVolEnvSuppAngName("3", "V", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VD[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VD[i],iVolNum++,3, spar);


// 3 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
   TString volEnvSuppAng3HD[4];
   volEnvSuppAng3HD[0] =  GetVolEnvSuppAngName("3", "H", "D", icount, 0);
   volEnvSuppAng3HD[1] =  GetVolEnvSuppAngName("3", "H", "D", icount, 1);
   volEnvSuppAng3HD[2] =  GetVolEnvSuppAngName("3", "H", "D", icount, 2);
   volEnvSuppAng3HD[3] =  GetVolEnvSuppAngName("3", "H", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[2], 0, true,   
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HD[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3H = GetVolAluAngSuppName("3", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HD[i],iVolNum++,3, spar);
       
// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe2D[4];
   volEnvInoxGasPipe2D[0] = GetVolEnvInoxGasPipeName("2", "D", icount, 0);
   volEnvInoxGasPipe2D[1] = GetVolEnvInoxGasPipeName("2", "D", icount, 1);
   volEnvInoxGasPipe2D[2] = GetVolEnvInoxGasPipeName("2", "D", icount, 2);
   volEnvInoxGasPipe2D[3] = GetVolEnvInoxGasPipeName("2", "D", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2D[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);

   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2D[i],iVolNum++,3, ppar);

// 4 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng4VD[4];
   volEnvSuppAng4VD[0] =  GetVolEnvSuppAngName("4", "V", "D", icount, 0);
   volEnvSuppAng4VD[1] =  GetVolEnvSuppAngName("4", "V", "D", icount, 1);
   volEnvSuppAng4VD[2] =  GetVolEnvSuppAngName("4", "V", "D", icount, 2);
   volEnvSuppAng4VD[3] =  GetVolEnvSuppAngName("4", "V", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VD[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport4V = GetVolAluAngSuppName("4", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VD[i],iVolNum++,3, spar);   

// 4 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];
   TString volEnvSuppAng4HD[4];
   volEnvSuppAng4HD[0] =  GetVolEnvSuppAngName("4", "H", "D", icount, 0);
   volEnvSuppAng4HD[1] =  GetVolEnvSuppAngName("4", "H", "D", icount, 1);
   volEnvSuppAng4HD[2] =  GetVolEnvSuppAngName("4", "H", "D", icount, 2);
   volEnvSuppAng4HD[3] =  GetVolEnvSuppAngName("4", "H", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HD[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4H = GetVolAluAngSuppName("4", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HD[i],iVolNum++,3, spar);   

// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   Double_t sparysavex=spar[1];
   TString volEnvSuppAngXHD[4];
   volEnvSuppAngXHD[0] =  GetVolEnvSuppAngName("X", "H", "D", icount, 0);
   volEnvSuppAngXHD[1] =  GetVolEnvSuppAngName("X", "H", "D", icount, 1);
   volEnvSuppAngXHD[2] =  GetVolEnvSuppAngName("X", "H", "D", icount, 2);
   volEnvSuppAngXHD[3] =  GetVolEnvSuppAngName("X", "H", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[0], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[1], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[2], 0, true,
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHD[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXH = GetVolAluAngSuppName("X", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHD[i],iVolNum++,3, spar);   

// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];
   TString volEnvSuppAngXVD[4];
   volEnvSuppAngXVD[0] =  GetVolEnvSuppAngName("X", "V", "D", icount, 0);
   volEnvSuppAngXVD[1] =  GetVolEnvSuppAngName("X", "V", "D", icount, 1);
   volEnvSuppAngXVD[2] =  GetVolEnvSuppAngName("X", "V", "D", icount, 2);
   volEnvSuppAngXVD[3] =  GetVolEnvSuppAngName("X", "V", "D", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[0], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[1], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[2], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVD[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXV = GetVolAluAngSuppName("X", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVD[i],iVolNum++,3, spar);	

   // keep values of yEnvP, yEnvM
   fYEnvPsave = yEnvP;
   fYEnvMsave = yEnvM;   
}    

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberTypeE(Int_t& iVolNum, Int_t icount)
{
// ratio of zpos1m/zpos1p and inverse for first plane
    Float_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
    Float_t zpm=1./zmp;
    Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);      
    Float_t xEnv = (fgkDXZERO+fgkXMAX/2.)*zRatio;
   
   Double_t dpar[3];    
   Double_t spar[3];    
   Double_t ppar[3];    

// E
   Float_t yEnvP = 0;
   Float_t yEnvM = 0;
   yEnvP = (fYEnvMsave + fgkYMIN * zRatio ) * zpm + fgkYMIN * zRatio;
   yEnvM = (fYEnvPsave + fgkYMIN * zRatio ) * zmp + fgkYMIN * zRatio;

   Int_t detElemId = (10+icount+1)*100+15;
   TString volEnv1 = GetVolEnvName(icount, 1);
   TString volEnv7 = GetVolEnvName(icount, 7);
   TString volEnv10 = GetVolEnvName(icount, 10);
   TString volEnv16 = GetVolEnvName(icount, 16);
   TString volAlu = GetVolumeName("volAlu", icount);

   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv1, detElemId, true, TGeoTranslation(xEnv,-yEnvM,fgkZm));
   detElemId = (10+icount+1)*100+3;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv7, detElemId, true, TGeoTranslation(xEnv, yEnvM,fgkZm));
   detElemId = (10+icount+1)*100+12;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv10, detElemId, true, TGeoTranslation(-xEnv,-yEnvP,fgkZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   detElemId = (10+icount+1)*100+6;
   GetEnvelopes(16+icount)->AddEnvelope(volEnv16, detElemId, true, TGeoTranslation(-xEnv, yEnvP,fgkZp),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv1,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv7,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv10,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv16,iVolNum++,3,dpar);

// horizontal cable support chamber type E

   char volEnvSupportHE[12][8];
   for(Int_t jj=0;jj<2;jj++){
     for(Int_t ii=0;ii<6;ii++){
       snprintf(volEnvSupportHE[6*jj+ii],8,"SA%dHE%d",icount+1,6*jj+ii);
     }
   }
      
   spar[0]=dpar[0]-(fgkYVSup[0]/2.)*zRatio;
   spar[1]=fgkSizeSupportCable[1];
   spar[2]=fgkSizeSupportCable[2];
   for(Int_t in=0;in<3;in++){
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in], 0, true, TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio,-(yEnvM+dpar[1]/2.*(in-1)),fgkSizeVSupExt[0]+spar[2]));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+3], 0, true, TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio, yEnvM+dpar[1]/2.*(in-1),fgkSizeVSupExt[0]+spar[2]));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+6], 0, true, TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio),-(yEnvP+dpar[1]/2.*(in-1)),-(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
     GetEnvelopes(16+icount)->AddEnvelope(volEnvSupportHE[in+9], 0, true, TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio), yEnvP+dpar[1]/2.*(in-1),-(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   }
   for(Int_t ii=0;ii<12;ii++) {
     TString volAluSupportH = GetVolumeName("volAluSupportH", icount);
     GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHE[ii],iVolNum++,3, spar);
   }  
   
// angular supports chamber type E
// 1 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   Double_t sparysave=spar[1];

   TString volEnvSuppAng1VE[4];
   volEnvSuppAng1VE[0] =  GetVolEnvSuppAngName("1", "V", "E", icount, 0);
   volEnvSuppAng1VE[1] =  GetVolEnvSuppAngName("1", "V", "E", icount, 1);
   volEnvSuppAng1VE[2] =  GetVolEnvSuppAngName("1", "V", "E", icount, 2);
   volEnvSuppAng1VE[3] =  GetVolEnvSuppAngName("1", "V", "E", icount, 3);

   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VE[3], 0, true,
   TGeoTranslation(-xEnv,yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
   TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
    for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VE[i],iVolNum++,3, spar);


// 1 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];

   TString volEnvSuppAng1HE[4];
   volEnvSuppAng1HE[0] =  GetVolEnvSuppAngName("1", "H", "E", icount, 0);
   volEnvSuppAng1HE[1] =  GetVolEnvSuppAngName("1", "H", "E", icount, 1);
   volEnvSuppAng1HE[2] =  GetVolEnvSuppAngName("1", "H", "E", icount, 2);
   volEnvSuppAng1HE[3] =  GetVolEnvSuppAngName("1", "H", "E", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[2], 0, true, TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HE[3], 0, true, TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HE[i],iVolNum++,3, spar);

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe1E[4];
   volEnvInoxGasPipe1E[0] = GetVolEnvInoxGasPipeName("1", "E", icount, 0);
   volEnvInoxGasPipe1E[1] = GetVolEnvInoxGasPipeName("1", "E", icount, 1);
   volEnvInoxGasPipe1E[2] = GetVolEnvInoxGasPipeName("1", "E", icount, 2);
   volEnvInoxGasPipe1E[3] = GetVolEnvInoxGasPipeName("1", "E", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[1], 0, true,
   TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1E[3], 0, true,
   TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);

    TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1E[i],iVolNum++,3, ppar);

// 2 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
   TString volEnvSuppAng2VE[4];
   volEnvSuppAng2VE[0] =  GetVolEnvSuppAngName("2", "V", "E", icount, 0);
   volEnvSuppAng2VE[1] =  GetVolEnvSuppAngName("2", "V", "E", icount, 1);
   volEnvSuppAng2VE[2] =  GetVolEnvSuppAngName("2", "V", "E", icount, 2);
   volEnvSuppAng2VE[3] =  GetVolEnvSuppAngName("2", "V", "E", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[2], 0, true,
    TGeoTranslation(-xEnv,-yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VE[3], 0, true,
    TGeoTranslation(-xEnv,yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VE[i],iVolNum++,3, spar);   

// 2 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];

   TString volEnvSuppAng2HE[4];
   volEnvSuppAng2HE[0] =  GetVolEnvSuppAngName("2", "H", "E", icount, 0);
   volEnvSuppAng2HE[1] =  GetVolEnvSuppAngName("2", "H", "E", icount, 1);
   volEnvSuppAng2HE[2] =  GetVolEnvSuppAngName("2", "H", "E", icount, 2);
   volEnvSuppAng2HE[3] =  GetVolEnvSuppAngName("2", "H", "E", icount, 3);

   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[2], 0, true, TGeoTranslation(-xEnv,-yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HE[3], 0, true, TGeoTranslation(-xEnv,yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HE[i],iVolNum++,3, spar);   

// 3 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];
   TString volEnvSuppAng3VE[4];
   volEnvSuppAng3VE[0] =  GetVolEnvSuppAngName("3", "V", "E", icount, 0);
   volEnvSuppAng3VE[1] =  GetVolEnvSuppAngName("3", "V", "E", icount, 1);
   volEnvSuppAng3VE[2] =  GetVolEnvSuppAngName("3", "V", "E", icount, 2);
   volEnvSuppAng3VE[3] =  GetVolEnvSuppAngName("3", "V", "E", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VE[3], 0, true,
   TGeoTranslation(-xEnv,yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
   TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VE[i],iVolNum++,3, spar);


// 3 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];

   TString volEnvSuppAng3HE[4];
   volEnvSuppAng3HE[0] =  GetVolEnvSuppAngName("3", "H", "E", icount, 0);
   volEnvSuppAng3HE[1] =  GetVolEnvSuppAngName("3", "H", "E", icount, 1);
   volEnvSuppAng3HE[2] =  GetVolEnvSuppAngName("3", "H", "E", icount, 2);
   volEnvSuppAng3HE[3] =  GetVolEnvSuppAngName("3", "H", "E", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[2], 0, true, TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HE[3], 0, true, TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3H = GetVolAluAngSuppName("3", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HE[i],iVolNum++,3, spar);

// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe2E[4];
   volEnvInoxGasPipe2E[0] = GetVolEnvInoxGasPipeName("2", "E", icount, 0);
   volEnvInoxGasPipe2E[1] = GetVolEnvInoxGasPipeName("2", "E", icount, 1);
   volEnvInoxGasPipe2E[2] = GetVolEnvInoxGasPipeName("2", "E", icount, 2);
   volEnvInoxGasPipe2E[3] = GetVolEnvInoxGasPipeName("2", "E", icount, 3);
         GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[1], 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2E[3], 0, true,
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);

   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2E[i],iVolNum++,3, ppar);

// 4 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];

   TString volEnvSuppAng4VE[4];
   volEnvSuppAng4VE[0] =  GetVolEnvSuppAngName("4", "V", "E", icount, 0);
   volEnvSuppAng4VE[1] =  GetVolEnvSuppAngName("4", "V", "E", icount, 1);
   volEnvSuppAng4VE[2] =  GetVolEnvSuppAngName("4", "V", "E", icount, 2);
   volEnvSuppAng4VE[3] =  GetVolEnvSuppAngName("4", "V", "E", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAng4VE[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VE[3], 0, true,
   TGeoTranslation(-xEnv,yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport4V = GetVolAluAngSuppName("4", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VE[i],iVolNum++,3, spar);   

// 4 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];

   TString volEnvSuppAng4HE[4];
   volEnvSuppAng4HE[0] =  GetVolEnvSuppAngName("4", "H", "E", icount, 0);
   volEnvSuppAng4HE[1] =  GetVolEnvSuppAngName("4", "H", "E", icount, 1);
   volEnvSuppAng4HE[2] =  GetVolEnvSuppAngName("4", "H", "E", icount, 2);
   volEnvSuppAng4HE[3] =  GetVolEnvSuppAngName("4", "H", "E", icount, 3);
   
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAng4HE[0], 0, true,
   TGeoTranslation(xEnv,-yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAng4HE[1], 0, true,
   TGeoTranslation(xEnv,yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAng4HE[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HE[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4H = GetVolAluAngSuppName("4", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HE[i],iVolNum++,3, spar);   

// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   Double_t sparysavex=spar[1];

   TString volEnvSuppAngXHE[4];
   volEnvSuppAngXHE[0] =  GetVolEnvSuppAngName("X", "H", "E", icount, 0);
   volEnvSuppAngXHE[1] =  GetVolEnvSuppAngName("X", "H", "E", icount, 1);
   volEnvSuppAngXHE[2] =  GetVolEnvSuppAngName("X", "H", "E", icount, 2);
   volEnvSuppAngXHE[3] =  GetVolEnvSuppAngName("X", "H", "E", icount, 3);
   
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXHE[0], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXHE[1], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[2], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHE[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXH = GetVolAluAngSuppName("X", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHE[i],iVolNum++,3, spar);   

// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];

   TString volEnvSuppAngXVE[4];
   volEnvSuppAngXVE[0] =  GetVolEnvSuppAngName("X", "V", "E", icount, 0);
   volEnvSuppAngXVE[1] =  GetVolEnvSuppAngName("X", "V", "E", icount, 1);
   volEnvSuppAngXVE[2] =  GetVolEnvSuppAngName("X", "V", "E", icount, 2);
   volEnvSuppAngXVE[3] =  GetVolEnvSuppAngName("X", "V", "E", icount, 3);
   
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXVE[0], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXVE[1], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXVE[2], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnvSuppAngXVE[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXV = GetVolAluAngSuppName("X", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVE[i],iVolNum++,3, spar);
   
   // keep values of yEnvP, yEnvM
   fYEnvPsave = yEnvP;
   fYEnvMsave = yEnvM;   
}

//______________________________________________________________________________
void AliMUONTriggerGeometryBuilder::BuildChamberTypeF(Int_t& iVolNum, Int_t icount)
{ 
// ratio of zpos1m/zpos1p and inverse for first plane
    Float_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
    Float_t zpm=1./zmp;
    Float_t zRatio = AliMUONConstants::DefaultRatioTriggerChamber(icount);      
    Float_t xEnv = (fgkDXZERO+fgkXMAX/2.)*zRatio;
    
   Double_t dpar[3];    
   Double_t spar[3];    
   Double_t ppar[3];    

// F
   dpar[0] = (fgkXMAX/2.)*zRatio;
   dpar[1] =  fgkYMIN*zRatio;
   dpar[2] = 0.4;   

   Float_t yEnvP = 0;
   Float_t yEnvM = 0;
   yEnvP = (fYEnvMsave + fgkYMIN * zRatio ) * zpm + fgkYMIN * zRatio;
   yEnvM = (fYEnvPsave + fgkYMIN * zRatio ) * zmp + fgkYMIN * zRatio;

   Int_t detElemId = (10+icount+1)*100+14;
   TString volEnv0 = GetVolEnvName(icount, 0);
   TString volEnv8 = GetVolEnvName(icount, 8);
   TString volEnv9 = GetVolEnvName(icount, 9);
   TString volEnv17 = GetVolEnvName(icount, 17);
   TString volAlu = GetVolumeName("volAlu", icount);
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv0, detElemId, true, TGeoTranslation(xEnv,-yEnvP,fgkZp));
   detElemId = (10+icount+1)*100+4;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv8, detElemId, true, TGeoTranslation(xEnv, yEnvP,fgkZp));
   detElemId = (10+icount+1)*100+13;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv9, detElemId, true, TGeoTranslation(-xEnv,-yEnvM,fgkZm),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   detElemId = (10+icount+1)*100+5;
   GetEnvelopes(16+icount)
   ->AddEnvelope(volEnv17, detElemId, true, TGeoTranslation(-xEnv, yEnvM,fgkZm),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv0,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv8,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv9,iVolNum++,3,dpar);
   GetEnvelopes(16+icount)->AddEnvelopeConstituentParam(volAlu,volEnv17,iVolNum++,3,dpar);

// horizontal cable support chamber type F

   char volEnvSupportHF[12][8];
   for(Int_t jj=0;jj<2;jj++){
     for(Int_t ii=0;ii<6;ii++){
       snprintf(volEnvSupportHF[6*jj+ii],8,"SA%dHF%d",icount+1,6*jj+ii);
     }
   }

   spar[0]=dpar[0]-(fgkYVSup[0]/2.)*zRatio;
   spar[1]=fgkSizeSupportCable[1];
   spar[2]=fgkSizeSupportCable[2];
   for(Int_t in=0;in<3;in++){
     GetEnvelopes(16+icount)
     ->AddEnvelope(volEnvSupportHF[in], 0, true, 
     TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio,-(yEnvP+dpar[1]/2.*(in-1)),
     -(fgkSizeVSupExt[0]+spar[2])));
     GetEnvelopes(16+icount)
     ->AddEnvelope(volEnvSupportHF[in+3], 0, true, 
     TGeoTranslation(xEnv+fgkYVSup[0]/2.*zRatio,yEnvP+dpar[1]/2.*(in-1),
     -(fgkSizeVSupExt[0]+spar[2])));
     GetEnvelopes(16+icount)
     ->AddEnvelope(volEnvSupportHF[in+6], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio),-(yEnvM+dpar[1]/2.*(in-1)),
     fgkSizeVSupExt[0]+spar[2]),
     TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
     GetEnvelopes(16+icount)
     ->AddEnvelope(volEnvSupportHF[in+9], 0, true, 
     TGeoTranslation(-(xEnv+fgkYVSup[0]/2.*zRatio), yEnvM+dpar[1]/2.*(in-1),
     fgkSizeVSupExt[0]+spar[2]),
     TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   }
   for(Int_t ii=0;ii<12;ii++) {
     TString volAluSupportH = GetVolumeName("volAluSupportH", icount);
     GetEnvelopes(16+icount)
     ->AddEnvelopeConstituentParam(volAluSupportH,volEnvSupportHF[ii],iVolNum++,3, spar);
   }  

// angular supports chamber type F
// 1 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   Double_t sparysave=spar[1];
   TString volEnvSuppAng1VF[4];
   volEnvSuppAng1VF[0] =  GetVolEnvSuppAngName("1", "V", "F", icount, 0);
   volEnvSuppAng1VF[1] =  GetVolEnvSuppAngName("1", "V", "F", icount, 1);
   volEnvSuppAng1VF[2] =  GetVolEnvSuppAngName("1", "V", "F", icount, 2);
   volEnvSuppAng1VF[3] =  GetVolEnvSuppAngName("1", "V", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1VF[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
 
   TString volAluAngSupport1V = GetVolAluAngSuppName("1", "V", icount);
   for (Int_t i=0;i<4;i++)
    GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1V,volEnvSuppAng1VF[i],iVolNum++,3, spar);

// 1 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];

   TString volEnvSuppAng1HF[4];
   volEnvSuppAng1HF[0] =  GetVolEnvSuppAngName("1", "H", "F", icount, 0);
   volEnvSuppAng1HF[1] =  GetVolEnvSuppAngName("1", "H", "F", icount, 1);
   volEnvSuppAng1HF[2] =  GetVolEnvSuppAngName("1", "H", "F", icount, 2);
   volEnvSuppAng1HF[3] =  GetVolEnvSuppAngName("1", "H", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng1HF[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport1H = GetVolAluAngSuppName("1", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport1H,volEnvSuppAng1HF[i],iVolNum++,3, spar);

// gas pipe (low)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe1F[4];
   volEnvInoxGasPipe1F[0] = GetVolEnvInoxGasPipeName("1", "F", icount, 0);
   volEnvInoxGasPipe1F[1] = GetVolEnvInoxGasPipeName("1", "F", icount, 1);
   volEnvInoxGasPipe1F[2] = GetVolEnvInoxGasPipeName("1", "F", icount, 2);
   volEnvInoxGasPipe1F[3] = GetVolEnvInoxGasPipeName("1", "F", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe1F[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave-fgkOffsetGasPipe,fgkZm),fRsupportpipe);

   TString volInoxGasPipe = GetVolumeName("volInoxGasPipe", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe1F[i],iVolNum++,3, ppar);

// 2 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];
 
   TString volEnvSuppAng2VF[4];
   volEnvSuppAng2VF[0] =  GetVolEnvSuppAngName("2", "V", "F", icount, 0);
   volEnvSuppAng2VF[1] =  GetVolEnvSuppAngName("2", "V", "F", icount, 1);
   volEnvSuppAng2VF[2] =  GetVolEnvSuppAngName("2", "V", "F", icount, 2);
   volEnvSuppAng2VF[3] =  GetVolEnvSuppAngName("2", "V", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[2], 0, true,
    TGeoTranslation(-xEnv,-yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2VF[3], 0, true,
    TGeoTranslation(-xEnv,yEnvM-dpar[1],-(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport2V = GetVolAluAngSuppName("2", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2V,volEnvSuppAng2VF[i],iVolNum++,3, spar);   

// 2 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];

   TString volEnvSuppAng2HF[4];
   volEnvSuppAng2HF[0] =  GetVolEnvSuppAngName("2", "H", "F", icount, 0);
   volEnvSuppAng2HF[1] =  GetVolEnvSuppAngName("2", "H", "F", icount, 1);
   volEnvSuppAng2HF[2] =  GetVolEnvSuppAngName("2", "H", "F", icount, 2);
   volEnvSuppAng2HF[3] =  GetVolEnvSuppAngName("2", "H", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP-dpar[1]-sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng2HF[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM-dpar[1]-sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport2H = GetVolAluAngSuppName("2", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport2H,volEnvSuppAng2HF[i],iVolNum++,3, spar);   

// 3 vertical	   
 
    spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1V[1];
   spar[2]=fgkSizeSupport1V[2];
   sparysave=spar[1];

   TString volEnvSuppAng3VF[4];
   volEnvSuppAng3VF[0] =  GetVolEnvSuppAngName("3", "V", "F", icount, 0);
   volEnvSuppAng3VF[1] =  GetVolEnvSuppAngName("3", "V", "F", icount, 1);
   volEnvSuppAng3VF[2] =  GetVolEnvSuppAngName("3", "V", "F", icount, 2);
   volEnvSuppAng3VF[3] =  GetVolEnvSuppAngName("3", "V", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkZp+dpar[2]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3VF[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],fgkZm-dpar[2]-spar[2]),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3V = GetVolAluAngSuppName("3", "V", icount);
   for (Int_t i=0;i<4;i++)
    GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3V,volEnvSuppAng3VF[i],iVolNum++,3, spar);

// 3 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport1H[1];
   spar[2]=fgkSizeSupport1H[2];
 
   TString volEnvSuppAng3HF[4];
   volEnvSuppAng3HF[0] =  GetVolEnvSuppAngName("3", "H", "F", icount, 0);
   volEnvSuppAng3HF[1] =  GetVolEnvSuppAngName("3", "H", "F", icount, 1);
   volEnvSuppAng3HF[2] =  GetVolEnvSuppAngName("3", "H", "F", icount, 2);
   volEnvSuppAng3HF[3] =  GetVolEnvSuppAngName("3", "H", "F", icount, 3);
  
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkZp+dpar[2]-spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng3HF[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,fgkZm-(dpar[2]-spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport3H = GetVolAluAngSuppName("3", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport3H,volEnvSuppAng3HF[i],iVolNum++,3, spar);
   
// gas pipe (high)
   ppar[0]=fgkSizeGasPipe[0];
   ppar[1]=fgkSizeGasPipe[1];
   ppar[2]=dpar[0];
   TString volEnvInoxGasPipe2F[4];
   volEnvInoxGasPipe2F[0] = GetVolEnvInoxGasPipeName("2", "F", icount, 0);
   volEnvInoxGasPipe2F[1] = GetVolEnvInoxGasPipeName("2", "F", icount, 1);
   volEnvInoxGasPipe2F[2] = GetVolEnvInoxGasPipeName("2", "F", icount, 2);
   volEnvInoxGasPipe2F[3] = GetVolEnvInoxGasPipeName("2", "F", icount, 3);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZp),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);
   GetEnvelopes(16+icount)->AddEnvelope(volEnvInoxGasPipe2F[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave+fgkOffsetGasPipe,fgkZm),fRsupportpipe);

   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volInoxGasPipe,volEnvInoxGasPipe2F[i],iVolNum++,3, ppar);

// 4 vertical	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2V[1];
   spar[2]=fgkSizeSupport2V[2];
   sparysave=spar[1]+fgkSizeSupport2H[1];

   TString volEnvSuppAng4VF[4];
   volEnvSuppAng4VF[0] =  GetVolEnvSuppAngName("4", "V", "F", icount, 0);
   volEnvSuppAng4VF[1] =  GetVolEnvSuppAngName("4", "V", "F", icount, 1);
   volEnvSuppAng4VF[2] =  GetVolEnvSuppAngName("4", "V", "F", icount, 2);
   volEnvSuppAng4VF[3] =  GetVolEnvSuppAngName("4", "V", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1],fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[2], 0, true,
   TGeoTranslation(-xEnv,-yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4VF[3], 0, true,
   TGeoTranslation(-xEnv,yEnvM+dpar[1],-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   
   TString volAluAngSupport4V = GetVolAluAngSuppName("4", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4V,volEnvSuppAng4VF[i],iVolNum++,3, spar);   

// 4 horizontal	   
 
   spar[0]=dpar[0];
   spar[1]=fgkSizeSupport2H[1];
   spar[2]=fgkSizeSupport2H[2];

   TString volEnvSuppAng4HF[4];
   volEnvSuppAng4HF[0] =  GetVolEnvSuppAngName("4", "H", "F", icount, 0);
   volEnvSuppAng4HF[1] =  GetVolEnvSuppAngName("4", "H", "F", icount, 1);
   volEnvSuppAng4HF[2] =  GetVolEnvSuppAngName("4", "H", "F", icount, 2);
   volEnvSuppAng4HF[3] =  GetVolEnvSuppAngName("4", "H", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[0], 0, true,
   TGeoTranslation(xEnv,-yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[1], 0, true,
   TGeoTranslation(xEnv,yEnvP+dpar[1]+sparysave,fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[2], 0, true, 
   TGeoTranslation(-xEnv,-yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAng4HF[3], 0, true, 
   TGeoTranslation(-xEnv,yEnvM+dpar[1]+sparysave,-(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupport4H = GetVolAluAngSuppName("4", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupport4H,volEnvSuppAng4HF[i],iVolNum++,3, spar);   


// X horizontal	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXH[1];
   spar[2]=fgkSizeSupportXH[2];
   Double_t sparysavex=spar[1];

   TString volEnvSuppAngXHF[4];
   volEnvSuppAngXHF[0] =  GetVolEnvSuppAngName("X", "H", "F", icount, 0);
   volEnvSuppAngXHF[1] =  GetVolEnvSuppAngName("X", "H", "F", icount, 1);
   volEnvSuppAngXHF[2] =  GetVolEnvSuppAngName("X", "H", "F", icount, 2);
   volEnvSuppAngXHF[3] =  GetVolEnvSuppAngName("X", "H", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[0], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[1], 0, true,
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0,
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[2], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXHF[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0,
   -(fgkSizeVSupExt[0]+spar[2])),
   TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXH = GetVolAluAngSuppName("X", "H", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXH,volEnvSuppAngXHF[i],iVolNum++,3, spar);   

// X vertical	   
   spar[0]=(fgkYVSup[3]-fgkYVSup[0])/2.*zRatio;
   spar[1]=fgkSizeSupportXV[1];
   spar[2]=fgkSizeSupportXV[2];

   TString volEnvSuppAngXVF[4];
   volEnvSuppAngXVF[0] =  GetVolEnvSuppAngName("X", "V", "F", icount, 0);
   volEnvSuppAngXVF[1] =  GetVolEnvSuppAngName("X", "V", "F", icount, 1);
   volEnvSuppAngXVF[2] =  GetVolEnvSuppAngName("X", "V", "F", icount, 2);
   volEnvSuppAngXVF[3] =  GetVolEnvSuppAngName("X", "V", "F", icount, 3);
   
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[0], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,-yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[1], 0, true, 
   TGeoTranslation(spar[0]+fgkYVSup[0]*zRatio,yEnvP+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   fgkSizeVSupExt[0]+spar[2]));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[2], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),-yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));
   GetEnvelopes(16+icount)->AddEnvelope(volEnvSuppAngXVF[3], 0, true, 
   TGeoTranslation(-(spar[0]+fgkYVSup[0]*zRatio),yEnvM+dpar[1]+sparysave+1.0+sparysavex+spar[1],
   -(fgkSizeVSupExt[0]+spar[2])),TGeoRotation("rot1",90.,180.,90.,90.,180.,0.));

   TString volAluAngSupportXV = GetVolAluAngSuppName("X", "V", icount);
   for (Int_t i=0;i<4;i++)
   GetEnvelopes(16+icount)
   ->AddEnvelopeConstituentParam(volAluAngSupportXV,volEnvSuppAngXVF[i],iVolNum++,3, spar);
   
      // keep values of yEnvP, yEnvM
   fYEnvPsave = yEnvP;
   fYEnvMsave = yEnvM;   
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
    
    fIdtmed = fMUON->GetIdtmed()->GetArray()-1099;
    fIdAir= fIdtmed[1100]; // medium 1
    fIdAlu1=fIdtmed[1103]; // medium 4
    fIdInox = fIdtmed[1128];       // medium 29 Stainless Steel (18%Cr,9%Ni,Fe) 
    
    Double_t dstation =  ( ( - AliMUONConstants::DefaultChamberZ(11)) - 
         ( - AliMUONConstants::DefaultChamberZ(10)) ) /2.1;
    Float_t par[3];
    par[2] = dstation;
 
    Int_t icount=0;  // chamber counter (0 1 2 3)
    
    for (Int_t istation=0; istation<2; istation++) { // loop on stations
      for (Int_t iplane=0; iplane<2; iplane++) { // loop on detection planes
   
          Int_t iVolNum=1; // counter Volume Number
          icount = Int_t(iplane<<0)+Int_t(istation<<1);
          
          cout << "## In AliMUONTriggerGeometryBuilder " << icount << endl;
    
          par[0] = AliMUONConstants::Rmin(5+istation); 
          par[1] = AliMUONConstants::Rmax(5+istation);
          Char_t volName[6];
          snprintf(volName,6,"%s%d", "SC",11+icount);
          TVirtualMC::GetMC()->Gsvolu(volName,"TUBE", fIdAir, par, 3);
          
// chamber prototype
          BuildChamberPrototype(icount);

// RPC supports (vertical)
          BuildRPCSupportsVertical(iVolNum, icount);
            
// RPC supports (horizontal)

// supports for cables
          BuildRPCSupportsHorizontal(icount);
          
// Angular supports for chambers
          BuildAngularSupportForChambers(icount);

// gas pipes
          BuildGasPipes(icount);

// chamber type A
          BuildChamberTypeA(iVolNum, icount);

// chamber type B (plus envelope chambers B & C)   
          BuildChamberTypeB(iVolNum, icount);

// chamber type D, E and F (same size)
          BuildChamberTypeD(iVolNum, icount);
          BuildChamberTypeE(iVolNum, icount);
          BuildChamberTypeF(iVolNum, icount);
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

    TGeoRotation st345inclination("rotbeam");
    st345inclination.RotateX(-AliMUONConstants::St345Inclination());

    Double_t zpos1= AliMUONConstants::DefaultChamberZ(10); 
    SetTransformation(16, TGeoTranslation(0., 0, zpos1), st345inclination);
    
    zpos1= AliMUONConstants::DefaultChamberZ(11); 
    SetTransformation(17, TGeoTranslation(0., 0, zpos1), st345inclination);

    zpos1= AliMUONConstants::DefaultChamberZ(12); 
    SetTransformation(18, TGeoTranslation(0., 0, zpos1), st345inclination);

    zpos1= AliMUONConstants::DefaultChamberZ(13); 
    SetTransformation(19, TGeoTranslation(0., 0, zpos1), st345inclination);
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

