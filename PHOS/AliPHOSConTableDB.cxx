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

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
// Class provides correspondence between "raw numbers" i.e. number of crustall 
// in prototype and PHOT AbsId numer, used in reconstruction.
// First it calculates correspondence automatically, assuming, that 
// prototype, having N raws and M columns is situated in the center 
// of middle (third) PHOS module. Then this correspondence can be edited 
// manually. One can convert Raw->AbsId and visa versa AbsId->RawId.
//
//*-- Author :  D.Peressounko ("RRC Kurchatov Institute") 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TArrayS.h"
#include "TH2S.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSConTableDB.h"

ClassImp(AliPHOSConTableDB)


//____________________________________________________________________________ 
  AliPHOSConTableDB::AliPHOSConTableDB():TNamed("AliPHOSConTableDB","Beamtest2002") 
{
//default constructor, nothing created.
  fNcrInProto = 0 ;
  fProtoRaws = 0 ;
  fProtoColumns = 0 ;
  fRawOffset = 0 ;
  fColOffset = 0 ;
  fGeom = 0;
  fAbsIdMap = 0 ;
  fRawIdMap = 0 ;
}

//____________________________________________________________________________ 
  AliPHOSConTableDB::AliPHOSConTableDB(const char * title):TNamed("AliPHOSConTableDB",title) 
{
 //Normally used constructor 
  fNcrInProto = 0 ;
  fProtoRaws = 0 ;
  fProtoColumns = 0 ;
  fRawOffset = 0 ;
  fColOffset = 0 ;
  fAbsIdMap = 0 ;
  fRawIdMap = 0 ;

  fGeom = AliPHOSGeometry::GetInstance("IHEP","") ;

}

//____________________________________________________________________________ 
AliPHOSConTableDB::AliPHOSConTableDB(const AliPHOSConTableDB& cdb):TNamed(cdb.GetName(), cdb.GetTitle()) 
{
  //Copy constructor
  
  fProtoRaws=cdb.fProtoRaws ;        //  Parameters
  fProtoColumns=cdb.fProtoColumns ;     //  used to calculate
  fRawOffset=cdb.fRawOffset ;        //  correspondance
  fColOffset=cdb.fColOffset ;        //  map
  fNcrInProto=cdb.fNcrInProto ;       //Number of channels in prototype
  fMinAbsId=cdb.fMinAbsId ;         //Minimal AbsId, corresponding to some prototype cristall.
  fMaxAbsId=cdb.fMaxAbsId ;         //Maximal AbsId, corresponding to some prototype cristall
  fAbsIdMap=new TArrayS(*(cdb.fAbsIdMap)) ;         //Map of correspondance between Raw and PHOS ID
  fRawIdMap=new TArrayS(*(cdb.fRawIdMap)) ;         //Map of correspondance between AbsId and Raw
}

//____________________________________________________________________________ 
  AliPHOSConTableDB::~AliPHOSConTableDB()
{
  if(fAbsIdMap)
    delete fAbsIdMap ;
  if(fRawIdMap)
    delete fRawIdMap ;
}

//____________________________________________________________________________ 
void  AliPHOSConTableDB::BuildDB(void)
{ 
  //Make a map between Protopype cristalls and PHOS crystalls
  //assuming, that prototype is centered in the third module of the PHOS
  fNcrInProto =fProtoRaws*fProtoColumns ;
  if(!fNcrInProto){
    AliError(Form("configuratio of prototype is not known!!!\n Specify number of raws and columns in prototype"));
    return ;
  }
  fRawOffset = (fGeom->GetNPhi() - fProtoRaws)/2 ;
  fColOffset = (fGeom->GetNZ() - fProtoColumns )/ 2 ;
  fAbsIdMap = new TArrayS(fNcrInProto) ;
  fMinAbsId = fGeom->GetNCristalsInModule()*2 +
    fRawOffset*fGeom->GetNZ()+fColOffset+1 ;
  fMaxAbsId = fGeom->GetNCristalsInModule()*2 +
    (fRawOffset + fProtoRaws)*fGeom->GetNZ()- 
     fColOffset ;
  fRawIdMap = new TArrayS(fMaxAbsId-fMinAbsId+1) ;
  for(Int_t raw =0; raw < fProtoRaws ; raw ++){
    for(Int_t col = 0; col < fProtoColumns ; col ++){
      Int_t rawId = raw*fProtoColumns + col ;
      Int_t rel[4] = {3,0,0,0} ; //We assume, that we deal with third module
      rel[2]=raw + fRawOffset+1 ;
      rel[3]=col + fColOffset+1 ;
      Int_t absId ;
      fGeom->RelToAbsNumbering(rel,absId) ;
      fAbsIdMap->AddAt(static_cast<UInt_t>(absId),rawId) ;
      fRawIdMap->AddAt(static_cast<UInt_t>(rawId),absId-fMinAbsId) ;
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
Int_t AliPHOSConTableDB::AbsId2Raw(Int_t absId)const{
  //converts numbering of modules in PHOS into
  //numbering in prototype
  if(absId >= fMinAbsId && absId<=fMaxAbsId){    
    return fRawIdMap->At(absId-fMinAbsId) ;
  }
  else
    return -1 ;
}
//____________________________________________________________________________ 
Int_t AliPHOSConTableDB::Raw2AbsId(Int_t rawId)const{
  //converts numbering of modules in prototipe into
  //numbering in PHOS
  if(rawId >= 0 && rawId<fNcrInProto)
    return fAbsIdMap->At(rawId) ;
  else
    return 0 ;
}
//____________________________________________________________________________ 
void AliPHOSConTableDB::Print(const Option_t *)const {
//prints configuraion

  TString message ; 
  message  = " %s %s\n" ;
  message += "PHOS Geometry configured for " ; 
  if(fGeom)
    message += "%s %s \n" ;
  else
    message += " null \n"  ;

  AliInfo(Form(message.Data(), GetName(), GetTitle(), fGeom->GetName(), fGeom->GetTitle() )) ; 

  message  = "\n-------Prototype parameters--------\n" ;
  message += "    number of columns: %d \n" ; 
  message += "    number of raws:    %d \n" ;
  message += "    centered in third PHOS module with offsets: \n " ;
  message += "    raw: %d of %d\n" ;
  message += "    col: %d of %d\n" ; 
  message += "------------------------------------ \n" ;

  AliInfo(Form(message.Data(), fProtoColumns, fProtoRaws, fRawOffset, fGeom->GetNPhi(), fColOffset,fGeom->GetNZ() ));   
}
//____________________________________________________________________________
AliPHOSConTableDB& AliPHOSConTableDB::operator=(const AliPHOSConTableDB& cdb){
//Operator for coding convetion
  fGeom=cdb.fGeom ;   //! PHOS geometry class
  fProtoRaws=cdb.fProtoRaws ;        //  Parameters
  fProtoColumns=cdb.fProtoColumns ;     //  used to calculate
  fRawOffset=cdb.fRawOffset ;        //  correspondance
  fColOffset=cdb.fColOffset ;        //  map
  fNcrInProto=cdb.fNcrInProto ;       //Number of channels in prototype
  fMinAbsId=cdb.fMinAbsId ;         //Minimal AbsId, corresponding to some prototype cristall.
  fMaxAbsId=cdb.fMaxAbsId ;         //Maximal AbsId, corresponding to some prototype cristall
  fAbsIdMap=new TArrayS(*(cdb.fAbsIdMap)) ;         //Map of correspondance between Raw and PHOS ID
  fRawIdMap=new TArrayS(*(cdb.fRawIdMap)) ;         //Map of correspondance between AbsId and Raw
  return *this ;
}



