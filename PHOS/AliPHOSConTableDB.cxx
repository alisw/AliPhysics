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

/* $Id$ */

//_________________________________________________________________________
// Short description  
//
//*-- Author :  (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TArrayS.h"
#include "TH2S.h"

// --- Standard library ---
#include <iostream.h>
// --- AliRoot header files ---
#include "AliPHOSGeometry.h"
#include "AliPHOSConTableDB.h"

ClassImp(AliPHOSConTableDB)


//____________________________________________________________________________ 
  AliPHOSConTableDB::AliPHOSConTableDB():TNamed("AliPHOSConTableDB","Beamtest2002") 
{
  fNcrInProto = 0 ;
  fProtoRaws = 0 ;
  fProtoColumns = 0 ;
  fRawOffset = 0 ;
  fColOffset = 0 ;
  fGeom = 0;
  fAbsIdMap = 0 ;
}

//____________________________________________________________________________ 
  AliPHOSConTableDB::AliPHOSConTableDB(const char * title):TNamed("AliPHOSConTableDB",title) 
{
  fNcrInProto = 0 ;
  fProtoRaws = 0 ;
  fProtoColumns = 0 ;
  fRawOffset = 0 ;
  fColOffset = 0 ;

  fGeom = AliPHOSGeometry::GetInstance("GPS2","") ;

}

//____________________________________________________________________________ 
  AliPHOSConTableDB::~AliPHOSConTableDB()
{
  if(fAbsIdMap)
    delete [] fAbsIdMap ;
}

//____________________________________________________________________________ 
void  AliPHOSConTableDB::BuildDB(void)
{ 
  //Make a map between Protopype cristalls and PHOS crystalls
  //assuming, that prototype is centered in the third module of the PHOS
  fNcrInProto =fProtoRaws*fProtoColumns ;
  if(!fNcrInProto){
    cout << "configuratio of prototype is not known!!!" << endl ;
    cout << "specify number of raws and columns in prototype" << endl ;
    return ;
  }
  fRawOffset = (fGeom->GetNPhi() - fProtoRaws)/2 ;
  fColOffset = (fGeom->GetNZ() - fProtoColumns )/ 2 ;
  fAbsIdMap = new TArrayS(fNcrInProto) ;
  for(Int_t raw =0; raw < fProtoRaws; raw ++){
    for(Int_t col = 0; col < fProtoColumns ; col ++){
      Int_t rawId = col*fProtoRaws + raw ;
      Int_t rel[4] = {3,0,0,0} ; //We assume, that we deal with third module
      rel[2]=raw + fRawOffset ;
      rel[3]=col + fColOffset ;
      Int_t absId ;
      fGeom->RelToAbsNumbering(rel,absId) ;
      fAbsIdMap->AddAt(static_cast<UInt_t>(absId),rawId) ;
    }
  }

}
//____________________________________________________________________________ 
void AliPHOSConTableDB::PlotProtoMap(Option_t * opt)
{
  //Visualyse connection table

  TH2S * hMapProto = new TH2S("hMap","Map of Prototype ids",
			      fGeom->GetNPhi(),0,fGeom->GetNPhi(),
			      fGeom->GetNZ(),0,fGeom->GetNZ()) ;
  TH2S * hMapPHOS = new TH2S("hMapPHOS","Map of PHOS ids",
			     fGeom->GetNPhi(),0,fGeom->GetNPhi(),
			     fGeom->GetNZ(),0,fGeom->GetNZ()) ;
  TH2C * hMapBox = new TH2C("hMapBox","Map of Prototype ids",
			      fGeom->GetNPhi(),0,fGeom->GetNPhi(),
			      fGeom->GetNZ(),0,fGeom->GetNZ()) ; 
  for(Int_t raw =0; raw <fGeom->GetNPhi() ; raw ++)
    for(Int_t col = 0; col <fGeom->GetNZ() ; col ++)
      hMapBox->SetBinContent(raw+1,col+1,1) ;
  
  for(Int_t raw =0; raw < fProtoRaws; raw ++){
    for(Int_t col = 0; col < fProtoColumns ; col ++){
      Int_t rawId = col*fProtoRaws + raw ;
      Int_t rel[4] = {3,0,0,0} ; //We assume, that we deal with third module
      rel[2]=raw + fRawOffset ;
      rel[3]=col + fColOffset ;
      hMapProto->SetBinContent(rel[2]+1,rel[3]+1,rawId);
      Int_t absId ;
      fGeom->RelToAbsNumbering(rel,absId) ;
      hMapPHOS->SetBinContent(rel[2]+1,rel[3]+1,absId) ;
    }
  }


  if(strstr(opt,"Zoom")||strstr(opt,"zoom")){
    static_cast<TAxis *>(hMapBox->GetXaxis())->SetRange(fRawOffset+1,fGeom->GetNPhi()-fRawOffset) ;
    static_cast<TAxis *>(hMapBox->GetYaxis())->SetRange(fColOffset+1,fGeom->GetNZ()-fColOffset) ;    
  }
   hMapBox->Draw("box") ;
   if(strstr(opt,"PHOS"))
     hMapPHOS->Draw("textsame") ;
   else
     hMapProto->Draw("textsame") ;

} 
//____________________________________________________________________________ 
Int_t AliPHOSConTableDB::Raw2AbsId(Int_t rawId){
  //converts numbering of modules in prototipe into
  //numbering in PHOS
  if(rawId >= 0 && rawId<fNcrInProto)
    return fAbsIdMap->At(rawId) ;
  else
    return 0 ;
}
//____________________________________________________________________________ 
void AliPHOSConTableDB::Print(Option_t * option)const {

  cout << GetName() <<  " " << GetTitle() << endl ;
  cout << "PHOS Geometry configured for " ; 
  if(fGeom)
    cout << fGeom->GetName() << " " << fGeom->GetTitle() <<  endl ;
  else
    cout << " null " << endl ;
  cout << "-------Prototype parameters--------" << endl ;
  cout << "    number of columns: " << fProtoColumns << endl ;
  cout << "    number of raws:    " << fProtoRaws << endl ;
  cout << "    centered in third PHOS module with offsets: " <<endl ;
  cout << "    raw: " << fRawOffset << " of " << fGeom->GetNPhi() << endl ;
  cout << "    col: " << fColOffset << " of " << fGeom->GetNZ() << endl ;
  cout << "------------------------------------" << endl ;
}
