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

//-----------------------------------------------------------------
//   Implementation of the alignment object class through the abstract
//  class AliAlignObj. From it two derived concrete representation of
//  alignment object class (AliAlignObjAngles, AliAlignObjMatrix) are
//  derived in separate files.
//-----------------------------------------------------------------
/*****************************************************************************
 * AliAlignObjAngles: derived alignment class storing alignment information  *
 *   for a single volume in form of three doubles for the translation        *
 *   and three doubles for the rotation expressed with the euler angles      *
 *   in the xyz-convention (http://mathworld.wolfram.com/EulerAngles.html),  *
 *   also known as roll, pitch, yaw. PLEASE NOTE THE ANGLES SIGNS ARE        *
 *   INVERSE WITH RESPECT TO THIS REFERENCE!!! In this way the representation*
 *   is fully consistent with the TGeo Rotation methods.                     *
 *****************************************************************************/

#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliAlignObjAngles.h"
 
ClassImp(AliAlignObj)

Int_t AliAlignObj::fgLayerSize[kLastLayer - kFirstLayer] = {
  80, 160,  // ITS SPD
  84, 176,  // ITS SDD
  748, 950, // ITS SSD
  36, 36,   // TPC
  90, 90, 90, 90, 90, 90,  // TRD
  1674,     // TOF
  1, 1,     // PHOS ??
  7,        // RICH ??
  1         // MUON ??
};

const char* AliAlignObj::fgLayerName[kLastLayer - kFirstLayer] = {
  "ITS inner pixels layer", "ITS outer pixels layer",
  "ITS inner drifts layer", "ITS outer drifts layer",
  "ITS inner strips layer", "ITS outer strips layer",
  "TPC inner chambers layer", "TPC outer chambers layer",
  "TRD chambers layer 1", "TRD chambers layer 2", "TRD chambers layer 3",
  "TRD chambers layer 4", "TRD chambers layer 5", "TRD chambers layer 6",
  "TOF layer",
  "?","?",
  "RICH layer",
  "?"
};

TString* AliAlignObj::fgVolPath[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

AliAlignObj** AliAlignObj::fgAlignObjs[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

//_____________________________________________________________________________
AliAlignObj::AliAlignObj():
  fVolUID(0)
{
  // default constructor
  InitVolPaths();
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const AliAlignObj& theAlignObj) :
  TObject(theAlignObj)
{
  //copy constructor
  fVolPath = theAlignObj.GetVolPath();
  fVolUID = theAlignObj.GetVolUID();
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  if(this==&theAlignObj) return *this;
  fVolPath = theAlignObj.GetVolPath();
  fVolUID = theAlignObj.GetVolUID();
  return *this;
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator*=(const AliAlignObj& theAlignObj)
{
  // multiplication operator
  // The operator can be used to 'combine'
  // two alignment objects
  TGeoHMatrix m1;
  GetMatrix(m1);
  TGeoHMatrix m2;
  theAlignObj.GetMatrix(m2);
  m1.MultiplyLeft(&m2);
  SetMatrix(m1);
  return *this;
}

//_____________________________________________________________________________
AliAlignObj::~AliAlignObj()
{
  // dummy destructor
}

//_____________________________________________________________________________
void AliAlignObj::SetVolUID(ELayerID detId, Int_t modId)
{
  // From detector name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for detID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  fVolUID = LayerToVolUID(detId,modId);
}

//_____________________________________________________________________________
void AliAlignObj::GetVolUID(ELayerID &layerId, Int_t &modId) const
{
  // From detector name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for detID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  layerId = VolUIDToLayer(fVolUID,modId);
}

//_____________________________________________________________________________
void AliAlignObj::AnglesToMatrix(const Double_t *angles, Double_t *rot) const
{
  // Calculates the rotation matrix using the 
  // Euler angles in "x y z" notation
  Double_t degrad = TMath::DegToRad();
  Double_t sinpsi = TMath::Sin(degrad*angles[0]);
  Double_t cospsi = TMath::Cos(degrad*angles[0]);
  Double_t sinthe = TMath::Sin(degrad*angles[1]);
  Double_t costhe = TMath::Cos(degrad*angles[1]);
  Double_t sinphi = TMath::Sin(degrad*angles[2]);
  Double_t cosphi = TMath::Cos(degrad*angles[2]);

  rot[0] =  costhe*cosphi;
  rot[1] = -costhe*sinphi;
  rot[2] =  sinthe;
  rot[3] =  sinpsi*sinthe*cosphi + cospsi*sinphi;
  rot[4] = -sinpsi*sinthe*sinphi + cospsi*cosphi;
  rot[5] = -costhe*sinpsi;
  rot[6] = -cospsi*sinthe*cosphi + sinpsi*sinphi;
  rot[7] =  cospsi*sinthe*sinphi + sinpsi*cosphi;
  rot[8] =  costhe*cospsi;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::MatrixToAngles(const Double_t *rot, Double_t *angles) const
{
  // Calculates the Euler angles in "x y z" notation
  // using the rotation matrix
  if(rot[0]<1e-7 || rot[8]<1e-7) return kFALSE;
  Double_t raddeg = TMath::RadToDeg();
  angles[0]=raddeg*TMath::ATan2(-rot[5],rot[8]);
  angles[1]=raddeg*TMath::ASin(rot[2]);
  angles[2]=raddeg*TMath::ATan2(-rot[1],rot[0]);
  return kTRUE;
}

//______________________________________________________________________________
void AliAlignObj::Transform(AliTrackPoint &p) const
{
  // The method transforms the space-point coordinates using the
  // transformation matrix provided by the AliAlignObj
  // The covariance matrix is not affected since we assume
  // that the transformations are sufficiently small

  if (fVolUID != p.GetVolumeID())
    AliWarning(Form("Alignment object ID is not equal to the space-point ID (%d != %d)",fVolUID,p.GetVolumeID())); 

  TGeoHMatrix m;
  GetMatrix(m);
  Double_t *rot = m.GetRotationMatrix();
  Double_t *tr  = m.GetTranslation();

  Float_t xyzin[3],xyzout[3];
  p.GetXYZ(xyzin);
  for (Int_t i = 0; i < 3; i++)
    xyzout[i] = tr[i]+
                xyzin[0]*rot[3*i]+
                xyzin[1]*rot[3*i+1]+
                xyzin[2]*rot[3*i+2];
  p.SetXYZ(xyzout);
  
}

//______________________________________________________________________________
void AliAlignObj::Transform(AliTrackPointArray &array) const
{
  AliTrackPoint p;
  for (Int_t i = 0; i < array.GetNPoints(); i++) {
    array.GetPoint(p,i);
    Transform(p);
    array.AddPoint(i,&p);
  }
}

//_____________________________________________________________________________
void AliAlignObj::Print(Option_t *) const
{
  // Print the contents of the
  // alignment object in angles and
  // matrix representations
  Double_t tr[3];
  GetTranslation(tr);
  Double_t angles[3];
  GetAngles(angles);
  TGeoHMatrix m;
  GetMatrix(m);
  const Double_t *rot = m.GetRotationMatrix();
//   printf("Volume=%s ID=%u\n", GetVolPath(),GetVolUID());
  ELayerID layerId;
  Int_t modId;
  GetVolUID(layerId,modId);
  printf("Volume=%s LayerID=%d ModuleID=%d\n", GetVolPath(),layerId,modId);
  printf("%12.6f%12.6f%12.6f    Tx = %12.6f    Psi   = %12.6f\n", rot[0], rot[1], rot[2], tr[0], angles[0]);
  printf("%12.6f%12.6f%12.6f    Ty = %12.6f    Theta = %12.6f\n", rot[3], rot[4], rot[5], tr[1], angles[1]);
  printf("%12.6f%12.6f%12.6f    Tz = %12.6f    Phi   = %12.6f\n", rot[6], rot[7], rot[8], tr[2], angles[2]);

}

//_____________________________________________________________________________
UShort_t AliAlignObj::LayerToVolUID(ELayerID layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
UShort_t AliAlignObj::LayerToVolUID(Int_t   layerId, Int_t modId)
{
  // From detector (layer) index and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
AliAlignObj::ELayerID AliAlignObj::VolUIDToLayer(UShort_t voluid, Int_t &modId)
{
  // From detector (layer) name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  modId = voluid & 0x7ff;

  return VolUIDToLayer(voluid);
}

//_____________________________________________________________________________
AliAlignObj::ELayerID AliAlignObj::VolUIDToLayer(UShort_t voluid)
{
  // From detector (layer) name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  return ELayerID((voluid >> 11) & 0x1f);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalPars(Double_t x, Double_t y, Double_t z,
				 Double_t psi, Double_t theta, Double_t phi)
{
  // Set the translations and angles by using parameters
  // defined in the local (in TGeo means) coordinate system
  // of the alignable volume. In case that the TGeo was
  // initialized, returns false and the object parameters are
  // not set.
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't set the alignment object parameters! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  const char* volpath = GetVolPath();
  TGeoPhysicalNode* node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(volpath);
  if (!node) {
    AliError(Form("Volume path %s not valid!",volpath));
    return kFALSE;
  }
  if (node->IsAligned())
    AliWarning(Form("Volume %s has been already misaligned!",volpath));

  TGeoHMatrix m;
  Double_t tr[3];
  tr[0]=x; tr[1]=y; tr[2]=z;
  m.SetTranslation(tr);
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  m.SetRotation(rot);

  TGeoHMatrix align,gprime,gprimeinv;
  gprime = *node->GetMatrix();
  gprimeinv = gprime.Inverse();
  m.Multiply(&gprimeinv);
  m.MultiplyLeft(&gprime);

  SetMatrix(m);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::ApplyToGeometry()
{
  // Apply the current alignment object
  // to the TGeo geometry

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't apply the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }
  
  const char* volpath = GetVolPath();
  TGeoPhysicalNode* node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(volpath);
  if (!node) {
    AliError(Form("Volume path %s not valid!",volpath));
    return kFALSE;
  }
  if (node->IsAligned()) {
    AliWarning(Form("Volume %s has been already misaligned!",volpath));
    return kFALSE;
  }

  TGeoHMatrix align,gprime;
  gprime = *node->GetMatrix();
  GetMatrix(align);
  gprime.MultiplyLeft(&align);
  TGeoHMatrix *ginv = new TGeoHMatrix;
  TGeoHMatrix *g = node->GetMatrix(node->GetLevel()-1);
  *ginv = g->Inverse();
  *ginv *= gprime;
  AliAlignObj::ELayerID layerId; // unique identity for volume in the alobj
  Int_t modId; // unique identity for volume in the alobj
  GetVolUID(layerId, modId);
  AliInfo(Form("Aligning volume %s of detector layer %d with local ID %d",volpath,layerId,modId));
  node->Align(ginv);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetFromGeometry(const char *path, AliAlignObj &alobj)
{
  // Get the alignment object which correspond
  // to the TGeo volume defined by the 'path'.
  // The method is extremely slow due to the
  // searching by string. Therefore it should
  // be used with great care!!

  // Reset the alignment object
  alobj.SetPars(0,0,0,0,0,0);
  alobj.SetVolPath(path);

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (!gGeoManager->GetListOfPhysicalNodes()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't contain any aligned nodes!");
    return kFALSE;
  }

  TObjArray* nodesArr = gGeoManager->GetListOfPhysicalNodes();
  TGeoPhysicalNode* node = NULL;
  for (Int_t iNode = 0; iNode < nodesArr->GetEntriesFast(); iNode++) {
    node = (TGeoPhysicalNode*) nodesArr->UncheckedAt(iNode);
    const char *nodePath = node->GetName();
    if (strcmp(path,nodePath) == 0) break;
  }
  if (!node) {
    if (!gGeoManager->cd(path)) {
      AliErrorClass(Form("Volume path %s not found!",path));
      return kFALSE;
    }
    else {
      AliWarningClass(Form("Volume (%s) has not been misaligned!",path));
      return kTRUE;
    }
  }

  TGeoHMatrix align,gprime,g,ginv,l;
  gprime = *node->GetMatrix();
  l = *node->GetOriginalMatrix();
  g = *node->GetMatrix(node->GetLevel()-1);
  g *= l;
  ginv = g.Inverse();
  align = gprime * ginv;
  alobj.SetMatrix(align);

  return kTRUE;
}

void  AliAlignObj::InitAlignObjFromGeometry()
{
  // Loop over all alignable volumes and extract
  // the corresponding alignment objects from
  // the TGeo geometry

  if(fgAlignObjs[0]) return;
  
  InitVolPaths();

  for (Int_t iLayer = 0; iLayer < (AliAlignObj::kLastLayer - AliAlignObj::kFirstLayer); iLayer++) {
    fgAlignObjs[iLayer] = new AliAlignObj*[AliAlignObj::LayerSize(iLayer)];
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer+ AliAlignObj::kFirstLayer,iModule);
      fgAlignObjs[iLayer][iModule] = new AliAlignObjAngles("",volid,0,0,0,0,0,0);
      const char *path = GetVolPath(volid);
      if (!GetFromGeometry(path, *fgAlignObjs[iLayer][iModule]))
	AliErrorClass(Form("Failed to extract the alignment object for the volume (ID=%d and path=%s) !",volid,path));
    }
  }
  
}

//_____________________________________________________________________________
AliAlignObj* AliAlignObj::GetAlignObj(ELayerID layerId, Int_t modId)
{
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  return fgAlignObjs[layerId-kFirstLayer][modId];
}

//_____________________________________________________________________________
const char* AliAlignObj::GetVolPath(ELayerID layerId, Int_t modId)
{
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  return fgVolPath[layerId-kFirstLayer][modId].Data();
}

//_____________________________________________________________________________
void AliAlignObj::InitVolPaths()
{
  // Initialize the LUTs which contain
  // the TGeo volume paths for each
  // alignable volume. The LUTs are
  // static, so they are created during
  // the creation of the first intance
  // of AliAlignObj

  if (fgVolPath[0]) return;

  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++)
    fgVolPath[iLayer] = new TString[fgLayerSize[iLayer]];

  /*********************       SPD layer1  ***********************/
  {
    Int_t modnum = 0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT12_1/I12B_"; //".../I12A_"
    TString str1 = "/I10B_";    //"/I10A_";
    TString str2 = "/I107_";    //"/I103_"
    //    TString str3 = "/I101_1/ITS1_1";
    TString volpath, volpath1, volpath2;

    for(Int_t c1 = 1; c1<=10; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 =1; c2<=2; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	volpath1 += str2;
	for(Int_t c3 =1; c3<=4; c3++){
	  volpath2 = volpath1;
	  volpath2 += c3;
	  //	  volpath2 += str3;
	  fgVolPath[kSPD1-kFirstLayer][modnum] = volpath2.Data();
	  modnum++;
	}
      }
    }
  }
  
  /*********************       SPD layer2  ***********************/
  {
    Int_t modnum = 0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT12_1/I12B_";  //".../I12A_"
    TString str1 = "/I20B_";  //"/I20A"
    TString str2 = "/I1D7_";  //"/I1D3"
    //    TString str3 = "/I1D1_1/ITS2_1";
    TString volpath, volpath1, volpath2;

    for(Int_t c1 = 1; c1<=10; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 =1; c2<=4; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	volpath1 += str2;
	for(Int_t c3 =1; c3<=4; c3++){
	  volpath2 = volpath1;
	  volpath2 += c3;
	  //	  volpath2 += str3;
	  fgVolPath[kSPD2-kFirstLayer][modnum] = volpath2.Data();
	  modnum++;
	}
      }
    }
  }

  /*********************       SDD layer1  ***********************/
  {
    Int_t modnum=0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT34_1/I004_";
    TString str1 = "/I302_";
    //    TString str2 = "/ITS3_1";
    TString volpath, volpath1;

    for(Int_t c1 = 1; c1<=14; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 =1; c2<=6; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	//	volpath1 += str2;
	fgVolPath[kSDD1-kFirstLayer][modnum] = volpath1.Data();
	modnum++;
      }
    }
  }

  /*********************       SDD layer2  ***********************/
  {
    Int_t modnum=0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT34_1/I005_";
    TString str1 = "/I402_";
    //    TString str2 = "/ITS4_1";
    TString volpath, volpath1;

    for(Int_t c1 = 1; c1<=22; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 = 1; c2<=8; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	//	volpath1 += str2;
	fgVolPath[kSDD2-kFirstLayer][modnum] = volpath1.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer1  ***********************/
  {
    Int_t modnum=0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT56_1/I565_";
    TString str1 = "/I562_";
    //    TString str2 = "/ITS5_1";
    TString volpath, volpath1;

    for(Int_t c1 = 1; c1<=34; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 = 1; c2<=22; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	//	volpath1 += str2;
	fgVolPath[kSSD1-kFirstLayer][modnum] = volpath1.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer1  ***********************/
  {
    Int_t modnum=0;
    TString str0 = "ALIC_1/ITSV_1/ITSD_1/IT56_1/I569_";
    TString str1 = "/I566_";
    //    TString str2 = "/ITS6_1";
    TString volpath, volpath1;

    for(Int_t c1 = 1; c1<=38; c1++){
      volpath = str0;
      volpath += c1;
      volpath += str1;
      for(Int_t c2 = 1; c2<=25; c2++){
	volpath1 = volpath;
	volpath1 += c2;
	//	volpath1 += str2;
	fgVolPath[kSSD2-kFirstLayer][modnum] = volpath1.Data();
	modnum++;
      }
    }
  }

  /***************    TPC inner chambers' layer    ****************/
  {
    Int_t modnum = 0;
    TString str1 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_1/TPC_SECT_";
    TString str2 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_2/TPC_SECT_";
    TString str_in = "/TPC_IROC_1";
    TString volpath;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      volpath = str1;
      volpath += cnt;
      volpath += str_in;
      fgVolPath[kTPC1-kFirstLayer][modnum] = volpath.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      volpath = str2;
      volpath += cnt;
      volpath += str_in;
      fgVolPath[kTPC1-kFirstLayer][modnum] = volpath.Data();
      modnum++;
    }
  }

  /***************    TPC outer chambers' layer    ****************/
  {
    Int_t modnum = 0;
    TString str1 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_1/TPC_SECT_";
    TString str2 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_2/TPC_SECT_";
    TString str_out = "/TPC_OROC_1";
    TString volpath;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      volpath = str1;
      volpath += cnt;
      volpath += str_out;
      fgVolPath[kTPC2-kFirstLayer][modnum] = volpath.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      volpath = str2;
      volpath += cnt;
      volpath += str_out;
      fgVolPath[kTPC2-kFirstLayer][modnum] = volpath.Data();
      modnum++;
    }
  }    

  /*********************       TOF layer   ***********************/
  {
    Int_t nstrA=15;
    Int_t nstrB=19;
    Int_t nstrC=20;
    Int_t nStripSec=nstrA+2*nstrB+2*nstrC;

    for (Int_t modnum=0; modnum < 1674; modnum++) {

      Int_t sector = modnum/nStripSec;
      Char_t  string1[100];
      Char_t  string2[100];

      Int_t icopy=-1;

      if(sector<3){
	icopy=sector+1;
	sprintf(string1,"/ALIC_1/B077_1/B075_%i/BTO3_1",icopy);
      }
      else if(sector<11){
	icopy=sector-2;
	sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1",icopy);
      }
      else if(sector==11 || sector==12){
	icopy=sector-10;
	sprintf(string1,"/ALIC_1/B077_1/B074_%i/BTO2_1",icopy);
      }
      else {
	icopy=sector-4;
	sprintf(string1,"/ALIC_1/B077_1/B071_%i/BTO1_1",icopy);
      }

      Int_t strInSec=modnum%nStripSec;

      if( strInSec < nstrC){
	icopy= nstrC - (strInSec+1) + 1;
	sprintf(string2,"FTOC_1/FLTC_0/FSTR_%i",icopy);
      }
      else if(strInSec< nstrC+nstrB){
 
	icopy= nstrB - (strInSec-nstrC+1) + 1;
	sprintf(string2,"FTOB_1/FLTB_0/FSTR_%i",icopy);

      }
      else if(strInSec< nstrC+nstrB+nstrA){   

	icopy= strInSec-(nstrC+nstrB)+1;
   	sprintf(string2,"FTOA_0/FLTA_0/FSTR_%i",icopy); 
      }
      else if(strInSec< nstrC+2*nstrB+nstrA){ 

	icopy= strInSec-(nstrC+nstrB+nstrA)+1;
 	sprintf(string2,"FTOB_2/FLTB_0/FSTR_%i",icopy);

      }
      else  { 

	icopy= strInSec-(nstrC+2*nstrB+nstrA)+1;
	sprintf(string2,"FTOC_2/FLTC_0/FSTR_%i",icopy);

      }
  
      Char_t  path[100];
      sprintf(path,"%s/%s",string1,string2); 
      //      printf("%d  %s\n",modnum,path);
      fgVolPath[kTOF-kFirstLayer][modnum] = path;
    }
  } 

  /*********************      RICH layer   ***********************/
  {
    TString str = "ALIC_1/RICH_";
    TString volpath;

    for (Int_t modnum=0; modnum < 7; modnum++) {
      volpath = str;
      volpath += (modnum+1);
      fgVolPath[kRICH-kFirstLayer][modnum] = volpath.Data();
    }
  }
}
