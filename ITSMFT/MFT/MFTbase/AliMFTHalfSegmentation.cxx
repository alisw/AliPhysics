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

//====================================================================================================================================================
//
//      Segmentation class for each half of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TNamed.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTGeometry.h"

ClassImp(AliMFTHalfSegmentation)

//====================================================================================================================================================

AliMFTHalfSegmentation::AliMFTHalfSegmentation(): 
  AliMFTVSegmentation(),
  fMFTHalfDisks(NULL)
{ 

  // default constructor

}

//====================================================================================================================================================

AliMFTHalfSegmentation::AliMFTHalfSegmentation(const AliMFTHalfSegmentation& source):
AliMFTVSegmentation(source),
fMFTHalfDisks(NULL){
  
  if (source.fMFTHalfDisks) fMFTHalfDisks = new TClonesArray(*(source.fMFTHalfDisks));

	
}
//====================================================================================================================================================

AliMFTHalfSegmentation::AliMFTHalfSegmentation(const Char_t *nameGeomFile, const Short_t id):
  AliMFTVSegmentation(),
  fMFTHalfDisks(NULL)
{
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  UInt_t halfUniqueID = mftGeom->GetObjectID(AliMFTGeometry::kHalfMFTType, id);
  SetUniqueID(halfUniqueID);
  SetName(Form("MFT_H_%d",id));
  
  
  fMFTHalfDisks = new TClonesArray("AliMFTHalfDiskSegmentation", AliMFTConstants::kNDisks);
  fMFTHalfDisks -> SetOwner(kTRUE);

  // Create XML engine
  TXMLEngine* geomFile = new TXMLEngine;
  
  // take access to main node
  XMLDocPointer_t  xmldoc   = geomFile->ParseFile(nameGeomFile);
  if (xmldoc==0) {
    delete geomFile;
    AliFatal(Form("Could not parse Geometry XML File named %s ", nameGeomFile));
  }
  XMLNodePointer_t mainnode = geomFile->DocGetRootElement(xmldoc);
  
  // Find  HAlf-MFT node in the XML file
  XMLNodePointer_t halfnode ;
  FindHalf(geomFile, mainnode, halfnode);

  // Create Half Disks belonging to that Half-MFT
  CreateHalfDisks(geomFile, halfnode);

  
  // Release memory
  geomFile->FreeDoc(xmldoc);
  delete geomFile;
  
//
//  Float_t zCenter, planeType, nLadderPerPlane, rMinSupport, rMaxSupport, rMinLadders, rMaxLadders, supportThickness, interspaceThickness;
//  Float_t ladderStiffenerThickness, ladderFlexThickness, chipActiveHeight, chipReadoutHeight, chipActiveSuperposition, sensorInterspace;
//  Float_t pixelSizeX, pixelSizeY;
//
//  Float_t equivalentSilicon, equivalentSiliconBeforeFront, equivalentSiliconBeforeBack, hasPixelRectangularPatternAlongY;
//
//  TNtuple *geomNtuple = (TNtuple*) initFile->Get("AliMFTGeometry");
//
//  geomNtuple -> SetBranchAddress("zCenter",                          &zCenter);
//  geomNtuple -> SetBranchAddress("planeType",                        &planeType);                       
//  geomNtuple -> SetBranchAddress("nLadderPerPlane",             &nLadderPerPlane);
//  geomNtuple -> SetBranchAddress("rMinSupport",		       	     &rMinSupport);
//  geomNtuple -> SetBranchAddress("rMaxSupport",		       	     &rMaxSupport);			   
//  geomNtuple -> SetBranchAddress("rMinLadders", 	       	     &rMinLadders);			   
//  geomNtuple -> SetBranchAddress("rMaxLadders",		       	     &rMaxLadders);			   
//  geomNtuple -> SetBranchAddress("supportThickness",	       	     &supportThickness);		   
//  geomNtuple -> SetBranchAddress("interspaceThickness",	       	     &interspaceThickness);		   
//  geomNtuple -> SetBranchAddress("ladderStiffenerThickness",   	     &ladderStiffenerThickness);	   
//  geomNtuple -> SetBranchAddress("ladderFlexThickness",	       	     &ladderFlexThickness);		   
//  geomNtuple -> SetBranchAddress("chipActiveHeight",	       	     &chipActiveHeight);
//  geomNtuple -> SetBranchAddress("chipReadoutHeight",	       	     &chipReadoutHeight);
//  geomNtuple -> SetBranchAddress("chipActiveSuperposition",    	     &chipActiveSuperposition);
//  geomNtuple -> SetBranchAddress("sensorInterspace",    	     &sensorInterspace);
//  geomNtuple -> SetBranchAddress("pixelSizeX",		       	     &pixelSizeX);
//  geomNtuple -> SetBranchAddress("pixelSizeY", 		       	     &pixelSizeY);			   
//
//  geomNtuple -> SetBranchAddress("equivalentSilicon",                &equivalentSilicon);
//  geomNtuple -> SetBranchAddress("equivalentSiliconBeforeFront",     &equivalentSiliconBeforeFront);
//  geomNtuple -> SetBranchAddress("equivalentSiliconBeforeBack",      &equivalentSiliconBeforeBack);
//
//  if (geomNtuple -> GetBranch("hasPixelRectangularPatternAlongY")) {
//    geomNtuple -> SetBranchAddress("hasPixelRectangularPatternAlongY", &hasPixelRectangularPatternAlongY);
//  }
//  else hasPixelRectangularPatternAlongY = 0.;
//  
//  Int_t nDisks = geomNtuple->GetEntries();
//  AliInfo(Form("%d half-disks will be created \n", nDisks));
//
//  Bool_t isDown = (GetID()==AliMFTSegmentation::kBottom);
//  AliInfo(Form("IsDown = %d \n", isDown));
//
//  for (Int_t iDisk=0; iDisk<nDisks; iDisk++) {
//
//    // Create new half disk
//
//
//    geomNtuple -> GetEntry(iDisk);
//    AliInfo(Form("Setting segmentation for MFT half disk #%02d with %2d ladders \n", iDisk,TMath::Nint(nLadderPerPlane)));
//    // zCenter = TMath::Abs(zCenter); ---- > Why ????
//
//    AliMFTHalfDiskSegmentation *halfDisk = new AliMFTHalfDiskSegmentation(Form("MFT_D_%1d_%1d", GetID(), iDisk), Form("MFT_D_%1d_%1d", GetID(), iDisk));
//
//    UInt_t halfDiskID = (1<<13) // Half-Disk Type
//                      + (GetID()<<12) //Half-MFT ID
//                      + (iDisk<<9);
//    AliInfo(Form("Init segmentation for MFT half disk #%02d with Unique ID = %d \n",iDisk, halfDiskID));
//
//    halfDisk -> Init(halfDiskID,
//		     isDown,
//		     zCenter, 
//		     TMath::Nint(planeType),
//         TMath::Nint(nLadderPerPlane),
//		     rMinSupport,
//		     rMaxSupport,
//		     rMinLadders, 
//		     rMaxLadders,
//		     supportThickness,
//		     interspaceThickness,                       // only for the kHollow type
//		     ladderStiffenerThickness,
//		     ladderFlexThickness,
//		     chipActiveHeight,
//		     chipReadoutHeight,
//		     chipActiveSuperposition,
//         sensorInterspace,
//		     pixelSizeX,
//		     pixelSizeY, 
//		     (hasPixelRectangularPatternAlongY>0.5));
////
////    halfDisk -> SetEquivalentSilicon(equivalentSilicon);
////    halfDisk -> SetEquivalentSiliconBeforeFront(equivalentSiliconBeforeFront);
////    halfDisk -> SetEquivalentSiliconBeforeBack(equivalentSiliconBeforeBack);
////    
////    new ((*fMFTHalfDisks)[fMFTHalfDisks->GetEntries()]) AliMFTHalfDiskSegmentation(*halfDisk);
////    delete halfDisk;
//  
//    AliInfo(Form("Done Setting segmentation for MFT half disk #%02d\n", iDisk));
//
//  }
  
  

}

//====================================================================================================================================================

AliMFTHalfSegmentation::~AliMFTHalfSegmentation() {

  if (fMFTHalfDisks) fMFTHalfDisks->Delete();
  delete fMFTHalfDisks; 
  
}

//====================================================================================================================================================

void AliMFTHalfSegmentation::Clear(const Option_t* /*opt*/) {

  if (fMFTHalfDisks) fMFTHalfDisks->Delete();
  delete fMFTHalfDisks; 
  fMFTHalfDisks = NULL;
  
}

//====================================================================================================================================================

THnSparseC* AliMFTHalfSegmentation::GetDetElem(Int_t detElemID) const {
      
  // Find det elem

  Int_t planeNb = detElemID/fNMaxDetElemPerPlane;
  Int_t detElemNb = detElemID - planeNb*fNMaxDetElemPerPlane;
  
  /// To do fixed !!!!!!!
  //THnSparseC *detElem = GetHalfDisk(planeNb)->GetActiveElement(detElemNb);
  THnSparseC *detElem=0;
  ///
  return detElem;

}

//====================================================================================================================================================

Bool_t AliMFTHalfSegmentation::Hit2PixelID(Double_t xHit, Double_t yHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel) {

  // xPixel and yPixel start from 0

  THnSparseC *detElem = GetDetElem(detElemID);

  if ( xHit<detElem->GetAxis(0)->GetXmin() ||
       xHit>detElem->GetAxis(0)->GetXmax() ||
       yHit<detElem->GetAxis(1)->GetXmin() ||
       yHit>detElem->GetAxis(1)->GetXmax() ) return kFALSE;

  xPixel = detElem->GetAxis(0)->FindBin(xHit) - 1;
  yPixel = detElem->GetAxis(1)->FindBin(yHit) - 1;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTHalfSegmentation::DoesPixelExist(Int_t detElemID, Int_t xPixel, Int_t yPixel) {

  THnSparseC *detElem = GetDetElem(detElemID);

  if (xPixel>=0 && xPixel<detElem->GetAxis(0)->GetNbins() && yPixel>=0 && yPixel<detElem->GetAxis(1)->GetNbins()) return kTRUE;
  else return kFALSE;

}

//====================================================================================
void AliMFTHalfSegmentation::CreateHalfDisks(TXMLEngine* xml, XMLNodePointer_t node)
{
  // this function display all accessible information about xml node and its children
  Int_t idisk;
  Int_t nladder;
  Double_t pos[3]={0., 0., 0.};
  Double_t ang[3]={0., 0., 0.};

  TString nodeName = xml->GetNodeName(node);
  if (!nodeName.CompareTo("disk")) {
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    while (attr!=0) {
      TString attrName = xml->GetAttrName(attr);
      TString attrVal  = xml->GetAttrValue(attr);
      if(!attrName.CompareTo("idisk")){
        idisk = attrVal.Atoi();
        if (idisk>=AliMFTConstants::kNDisks || idisk<0) {
          AliFatal(Form(" Wrong Disk number :  %d ",idisk));
        }
      } else
      if(!attrName.CompareTo("nladder")){
        nladder = attrVal.Atoi();
      } else
      if(!attrName.CompareTo("xpos")){
        pos[0] = attrVal.Atof();
      } else
      if(!attrName.CompareTo("ypos")){
        pos[1] = attrVal.Atof();
      } else
      if(!attrName.CompareTo("zpos")){
        pos[2] = attrVal.Atof();
      } else
        if(!attrName.CompareTo("phi")){
          ang[0] = attrVal.Atof();
        } else
          if(!attrName.CompareTo("theta")){
            ang[1] = attrVal.Atof();
          } else
            if(!attrName.CompareTo("psi")){
              ang[2] = attrVal.Atof();
            } else{
        AliError(Form(" Unknwon Attribute name %s ",xml->GetAttrName(attr)));
      }
      
      attr = xml->GetNextAttr(attr);
    }
    
    AliDebug(1,Form("Creating Half-Disk %d with %d Ladders at the position (%.2f,%.2f,%.2f) with angles  (%.2f,%.2f,%.2f)",idisk,nladder,pos[0],pos[1],pos[2],ang[0],ang[1],ang[2]));
    
    AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
    
    UInt_t diskUniqueID = mftGeom->GetObjectID(AliMFTGeometry::kHalfDiskType,
                                                 mftGeom->GetHalfMFTID(GetUniqueID()),
                                                 idisk );

    
    AliMFTHalfDiskSegmentation *halfDisk = new AliMFTHalfDiskSegmentation(diskUniqueID);
    halfDisk->SetPosition(pos);
    halfDisk->SetRotationAngles(ang);
    halfDisk->SetNLadders(nladder);
    halfDisk->CreateLadders(xml, node);
    if(halfDisk->GetNLaddersBuild() != halfDisk->GetNLadders()){
      AliFatal(Form("Number of ladder build (%d) does not correspond to the number declared (%d) : Check XML file",halfDisk->GetNLaddersBuild(), halfDisk->GetNLadders()));
    }
    new ((*fMFTHalfDisks)[idisk]) AliMFTHalfDiskSegmentation(*halfDisk);
    delete halfDisk;
    //GetHalfDisk(idisk)->Print("ls");

  }

  // display all child nodes
  XMLNodePointer_t child = xml->GetChild(node);
  while (child!=0) {
    CreateHalfDisks(xml, child);
    child = xml->GetNext(child);
  }
}
//====================================================================================

void AliMFTHalfSegmentation::FindHalf(TXMLEngine* xml, XMLNodePointer_t node, XMLNodePointer_t &retnode){
  // Find in the XML Geometry File the node corresponding to the Half-MFT being build
  // Set Position and Orientation of the Half-MFT
  Int_t isTop;
  Int_t ndisk;
  Double_t pos[3] = {0., 0., 0.};
  Double_t ang[3] = {0., 0., 0.};

  TString nodeName = xml->GetNodeName(node);
  if (!nodeName.CompareTo("half")) {
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    while (attr!=0) {
      TString attrName = xml->GetAttrName(attr);
      TString attrVal  = xml->GetAttrValue(attr);
      if(!attrName.CompareTo("top")){
        isTop = attrVal.Atoi();
        if (isTop>1 || isTop<0) {
          AliFatal(Form(" Wrong Half MFT number  :  %d ",isTop));
        }
      } else
        if(!attrName.CompareTo("ndisk")){
          ndisk = attrVal.Atoi();
          if (ndisk>5 || ndisk<0) {
            AliError(Form(" Wrong number of disk :  %d ",ndisk));
          }
          
        } else
          if(!attrName.CompareTo("xpos")){
            pos[0] = attrVal.Atof();
          } else
            if(!attrName.CompareTo("ypos")){
              pos[1] = attrVal.Atof();
            } else
              if(!attrName.CompareTo("zpos")){
                pos[2] = attrVal.Atof();
              } else
                if(!attrName.CompareTo("phi")){
                  ang[0] = attrVal.Atof();
                } else
                  if(!attrName.CompareTo("theta")){
                    ang[1] = attrVal.Atof();
                  } else
                    if(!attrName.CompareTo("psi")){
                      ang[2] = attrVal.Atof();
                    } else{
                AliError(Form(" Unknwon Attribute name %s ",xml->GetAttrName(attr)));
              }
      
      attr = xml->GetNextAttr(attr);
    }
    
    AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
    if(isTop == mftGeom->GetHalfMFTID(GetUniqueID())) {
      AliDebug(1,Form("Setting up %s Half-MFT  %d Disk(s) at the position (%.2f,%.2f,%.2f) with angles (%.2f,%.2f,%.2f)",(isTop?"Top":"Bottom"),ndisk,pos[0],pos[1],pos[2],ang[0],ang[1],ang[2]));
      SetPosition(pos);
      SetRotationAngles(ang);
      retnode = node;
      return;
    }
    
  }
  
  // display all child nodes
  XMLNodePointer_t child = xml->GetChild(node);
  while (child!=0) {
    FindHalf(xml, child, retnode);
    child = xml->GetNext(child);
  }
}

