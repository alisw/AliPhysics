/*************************************************************************
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
/// \class AliMFTHalfDiskSegmentation
///
/// Class for the description of the structure a Half-Disk
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "AliLog.h"
#include "AliMFTConstants.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTGeometry.h"

/// \cond CLASSIMP
ClassImp(AliMFTHalfDiskSegmentation);
/// \endcond


//====================================================================================================================================================
/// Default constructor
AliMFTHalfDiskSegmentation::AliMFTHalfDiskSegmentation():
  AliMFTVSegmentation(),
  fNLadders(0),
  fLadders(NULL)
{


}

//====================================================================================================================================================
/// Constructor
/// \param [in] uniqueID UInt_t: Unique ID of the Half-Disk to build
AliMFTHalfDiskSegmentation::AliMFTHalfDiskSegmentation(UInt_t uniqueID):
  AliMFTVSegmentation(),
  fNLadders(0),
  fLadders(NULL)
{

  // constructor
  SetUniqueID(uniqueID);

  AliDebug(1,Form("Start Creating Half-Disk UniqueID = %d",GetUniqueID()));
  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  SetName(Form("MFT_D_%d_%d",mftGeom->GetHalfMFTID(GetUniqueID()), mftGeom->GetHalfDiskID(GetUniqueID()) ));
  
  fLadders  = new TClonesArray("AliMFTLadderSegmentation");
  fLadders -> SetOwner(kTRUE);
  
  
  
}

//====================================================================================================================================================
/// Copy Constructor
AliMFTHalfDiskSegmentation::AliMFTHalfDiskSegmentation(const AliMFTHalfDiskSegmentation& input):
  AliMFTVSegmentation(input),
  fNLadders(input.fNLadders)
{
  
  // copy constructor
  if(input.fLadders)  fLadders  = new TClonesArray(*(input.fLadders));
  else   fLadders  = new TClonesArray("AliMFTLadderSegmentation");
  fLadders -> SetOwner(kTRUE);

}

//====================================================================================================================================================

AliMFTHalfDiskSegmentation::~AliMFTHalfDiskSegmentation() {

  Clear("");

}

//====================================================================================================================================================
/// Clear the TClonesArray holding the ladder segmentations

void AliMFTHalfDiskSegmentation::Clear(const Option_t* /*opt*/) {

  if (fLadders) fLadders->Delete();
  delete fLadders; 
  fLadders = NULL;

}


//====================================================================================================================================================
/// Creates the Ladders on this half-Disk based on the information contained in the XML file
void AliMFTHalfDiskSegmentation::CreateLadders(TXMLEngine* xml, XMLNodePointer_t node)
{
  Int_t iladder;
  Int_t nsensor;
  Double_t pos[3];
  Double_t ang[3]={0.,0.,0.};

  TString nodeName = xml->GetNodeName(node);
  if (!nodeName.CompareTo("ladder")) {
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    while (attr!=0) {
      TString attrName = xml->GetAttrName(attr);
      TString attrVal  = xml->GetAttrValue(attr);
      if(!attrName.CompareTo("iladder")){
        iladder = attrVal.Atoi();
        if (iladder>=GetNLadders() || iladder<0) {
          AliFatal(Form(" Wrong Ladder number :  %d ",iladder));
        }

      } else
        if(!attrName.CompareTo("nsensor")){
          nsensor = attrVal.Atoi();
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
    
    Int_t plane = -1;
    Int_t ladderID=iladder;
    if( iladder < GetNLadders()/2) {
      plane = 0;
    } else {
      plane = 1;
      // ladderID -= GetNLadders()/2;
    }
    
//    if ((plane==0 && pos[2]<0.) || (plane==1 && pos[2]>0.))
//      AliFatal(Form(" Wrong Z Position or ladder number ???  :  z= %f ladder id = %d",pos[2],ladderID));

    AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();

    UInt_t ladderUniqueID = mftGeom->GetObjectID(AliMFTGeometry::kLadderType,
                                                 mftGeom->GetHalfMFTID(GetUniqueID()),
                                                 mftGeom->GetHalfDiskID(GetUniqueID()),
                                                 ladderID);

//    UInt_t ladderUniqueID = (AliMFTGeometry::kLadderType<<13) +  (((GetUniqueID()>>9) & 0xF)<<9) + (plane<<8) + (ladderID<<3);

    
    AliMFTLadderSegmentation * ladder = new AliMFTLadderSegmentation(ladderUniqueID);
    ladder->SetNSensors(nsensor);
    ladder->SetPosition(pos);
    ladder->SetRotationAngles(ang);

    /// @todo : In the XML geometry file, the position of the top-left corner of the chip closest to the pipe is given in the Halfdisk coordinate system.
    /// Need to put in the XML file the position of the ladder coordinate center
    // Find the position of the corner of the flex which is the ladder corrdinate system center.
    
    pos[0] = -AliMFTConstants::kSensorSideOffset;
    pos[1] = -AliMFTConstants::kSensorTopOffset - AliMFTConstants::kSensorHeight;
    pos[2] = -AliMFTConstants::kFlexThickness - AliMFTConstants::kSensorThickness;
    Double_t master[3];
    ladder->GetTransformation()->LocalToMaster(pos, master);
    ladder->SetPosition(master);
    AliDebug(2,Form("Creating Ladder %2d with %d Sensors at the position (%.2f,%.2f,%.2f) with angles (%.2f,%.2f,%.2f) and ID = %d",iladder,nsensor,master[0],master[1],master[2],ang[0],ang[1],ang[2], ladderUniqueID ) );

    
    ladder->CreateSensors();

    new ((*fLadders)[iladder]) AliMFTLadderSegmentation(*ladder);
    delete ladder;

    //GetLadder(iladder)->Print();

    
  }
  
  // display all child nodes
  XMLNodePointer_t child = xml->GetChild(node);
  while (child!=0) {
    CreateLadders(xml, child);
    child = xml->GetNext(child);
  }
}



//==================================================================================================================
/// Returns the number of sensors on the Half-Disk
Int_t AliMFTHalfDiskSegmentation::GetNChips() {

  Int_t nChips = 0;

  for (Int_t iLadder=0; iLadder<fLadders->GetEntries(); iLadder++) {

    AliMFTLadderSegmentation *ladder = (AliMFTLadderSegmentation*) fLadders->At(iLadder);
    nChips += ladder -> GetNSensors();

  }

  return nChips;

}

//==================================================================================================================
/// Print out Half-Disk information
/// \param [in] opt "l" or "ladder" -> The ladder information will be printed out as well

void AliMFTHalfDiskSegmentation::Print(Option_t* opt){

  AliInfo(Form("Half-Disk %s (Unique ID = %d)",GetName(),GetUniqueID()));
  GetTransformation()->Print();
  AliInfo(Form("N Ladders = %d",fNLadders));
  if(opt && (strstr(opt,"ladder")||strstr(opt,"l"))){
    for (int i=0; i<GetNLadders(); i++)  GetLadder(i)->Print(opt);
    
  }
}
