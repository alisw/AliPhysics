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
/// \class AliMFTHalfSegmentation
///
/// Segmentation class for each half of the ALICE Muon Forward Tracker
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TClonesArray.h"

#include "AliLog.h"

#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTGeometry.h"

/// \cond CLASSIMP
ClassImp(AliMFTHalfSegmentation);
/// \endcond

//====================================================================================================================================================
/// Default constructor
AliMFTHalfSegmentation::AliMFTHalfSegmentation():
  AliMFTVSegmentation(),
  fMFTHalfDisks(NULL)
{ 


}

//====================================================================================================================================================
/// Copy constructor

AliMFTHalfSegmentation::AliMFTHalfSegmentation(const AliMFTHalfSegmentation& source):
AliMFTVSegmentation(source),
fMFTHalfDisks(NULL){
  
  if (source.fMFTHalfDisks) fMFTHalfDisks = new TClonesArray(*(source.fMFTHalfDisks));

	
}
//====================================================================================================================================================
/// Constructor
/// \param nameGeomFile Char_t * : name of the XML geometry file.
/// By default it is : $ALICE_ROOT/ITSMFT/MFT/data/AliMFTGeometry.xml
/// \param id Short_t : ID Of the Half-MFT to build (0=Bottom; 1=Top)
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
  
  // Find  Half-MFT node in the XML file
  XMLNodePointer_t halfnode ;
  FindHalf(geomFile, mainnode, halfnode);

  // Create Half Disks belonging to that Half-MFT
  CreateHalfDisks(geomFile, halfnode);

  
  // Release memory
  geomFile->FreeDoc(xmldoc);
  delete geomFile;
  

}

//====================================================================================================================================================

AliMFTHalfSegmentation::~AliMFTHalfSegmentation() {

  if (fMFTHalfDisks) fMFTHalfDisks->Delete();
  delete fMFTHalfDisks; 
  
}

//====================================================================================================================================================
///Clear the TClonesArray holding the AliMFTHalfDiskSegmentation objects

void AliMFTHalfSegmentation::Clear(const Option_t* /*opt*/) {

  if (fMFTHalfDisks) fMFTHalfDisks->Delete();
  delete fMFTHalfDisks; 
  fMFTHalfDisks = NULL;
  
}

//====================================================================================
///Create the Half-Disks 
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
/// Find Half-Disk in the XML file (private)

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

