/* *************************************************************************
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

/*
$Log$
Revision 1.8  2000/07/12 08:56:32  fca
Coding convention correction and warning removal

Revision 1.7  2000/06/28 21:27:45  morsch
Most coding rule violations corrected.
Still to do: Split the file (on file per class) ? Avoid the global variables.
Copy constructors and assignment operators (dummy ?)

Revision 1.6  2000/04/14 11:07:46  morsch
Correct volume to medium assignment in case several media are asigned to the
same material.

Revision 1.5  2000/03/20 15:11:03  fca
Mods to make the code compile on HP

Revision 1.4  2000/01/18 16:12:08  morsch
Bug in calculation of number of volume divisions and number of positionings corrected
Browser for Material and Media properties added

Revision 1.3  1999/11/14 14:31:14  fca
Correct small error and remove compilation warnings on HP

Revision 1.2  1999/11/10 16:53:35  fca
The new geometry viewer from A.Morsch

*/

/* 
 *  Version: 0
 *  Written by Andreas Morsch
 *  
 * 
 *
 * For questions critics and suggestions to this part of the code
 * contact andreas.morsch@cern.ch
 * 
 **************************************************************************/


#include "AliGeant3GeometryGUI.h"
#include "AliDrawVolume.h"
#include "AliGUIMaterial.h"
#include "AliGUIMedium.h"
#include "AliGuiGeomMain.h"

AliDrawVolume  *gCurrentVolume   = new AliDrawVolume("NULL");
AliGUIMaterial *gCurrentMaterial = new AliGUIMaterial();
AliGUIMedium   *gCurrentMedium   = new AliGUIMedium();

ClassImp(AliGeant3GeometryGUI)

    AliGeant3GeometryGUI::AliGeant3GeometryGUI()
{
// Constructor
    fPanel  =   new AliGuiGeomMain(gClient->GetRoot(), 500, 500);
    fNstack = 0;
    fVolumes =   new TClonesArray("AliDrawVolume",1000);
    fMaterials = new TClonesArray("AliGUIMaterial",1000);
    fMedia =     new TClonesArray("AliGUIMedium",1000);
//  Store local copy of zebra bank entries
    TGeant3 *geant3 = (TGeant3*) gMC;
    if (geant3) {
	fZlq=geant3->Lq();
	fZq=geant3->Q();
	fZiq=geant3->Iq();
	fGclink=geant3->Gclink();
	fGcnum=geant3->Gcnum();
//
	ReadGeometryTree();
	ReadMaterials();
    }
}
void AliGeant3GeometryGUI::Streamer(TBuffer &)
{
// Dummy Streamer
;
}


void AliGeant3GeometryGUI::ReadGeometryTree()
{
//
// Copy zebra volume tree into ROOT LisTree
//
    char *vname;
    char /* *namec, */ *tmp;
    char namec[30];
//  Icons for 
//  Closed folder (=volume containing children)
    const TGPicture* kFolder     = gClient->GetPicture("folder_t.xpm");
//  Open folder   (=volume containing children)
    const TGPicture* kOpenFolder = gClient->GetPicture("ofolder_t.xpm");
//  Empty object
    const TGPicture* kDocument   = gClient->GetPicture("doc_t.xpm");

    AliDrawVolume  *volume;

    Int_t nst=0;
    Int_t nlevel=1;
    Int_t newlevel=nlevel;

    volume = new AliDrawVolume("ALIC");
    volume->SetIdVolume(((TGeant3*)gMC)->VolId("ALIC"));
    volume->SetIdCopy(0);
    volume->SetItem(NULL);
    (*fVolumes)[0]=volume;
//
//  Loop over volumes for which information has been collected
    while(nlevel>nst) {
	for (Int_t i=nst; i<nlevel; i++) 
	{
	    TGListTreeItem *itemi, *item2;
// GEANT3 volume number
	    Int_t ivol=TMath::Abs(Volume(i)->GetIdVolume());
// Copy number
// icopy=1 normal positioning
// icopy>1 positioning with parposp
// icopy<0 division
	    Int_t icopy = Volume(i)->GetIdCopy();
// Medium and material number, handle special case of divisions
	    Int_t imat, imed;
	    
	    if (icopy <0) {
		imed=Medium(-ivol);
		imat=Material(-ivol);
	    } else {
		imed=Medium(ivol);
		imat=Material(ivol);
	    }
//
// Number of children
	    Int_t nch = NChildren(ivol);
	    strcpy(namec,((TGeant3*)gMC)->VolName(ivol));
	    if (nch >= 0) {
		printf("\n %s has %d children  \n ", namec,  nch);
	    } else {
		printf("\n %s has  divisions \n ", namec);
	    }
//
// Name to be used in ListTree
	    vname = new char[5];
	    strncpy(vname,namec, 4);
	    vname[4]='\0';

	    if (icopy >1) {
		sprintf(namec,"%s*%3dPos",namec,icopy);
	    } else if (icopy <0) {
		sprintf(namec,"%s*%3dDiv",namec,-icopy);
	    }
	    if (i>0) {
		itemi=Volume(i)->GetItem();
	    } else {
		itemi=NULL;
	    }
	    
//
// Add volume to list tree
	    
	    if (nch!=0) {
		item2 = fPanel->AddItem(new AliDrawVolume(vname), 
				    itemi, namec, kOpenFolder, kFolder);
	    } else {
		item2 = fPanel->AddItem(new AliDrawVolume(vname), 
				    itemi, namec, kDocument, kDocument);
	    }
//
// Add medium information to list tree item
	    ((AliDrawVolume *) item2->GetUserData())->SetIdVolume(ivol);
	    ((AliDrawVolume *) item2->GetUserData())->SetIdMaterial(imat);
	    ((AliDrawVolume *) item2->GetUserData())->SetIdMedium(imed);
//
// Set current volume to first list tree entry
	    if (!i) gCurrentVolume= ((AliDrawVolume *) item2->GetUserData());
	    
//
// Collect children information
//
// nch < 0: Children by division	    
	    if (nch < 0) {
//
// Geant volume number
		Int_t icvol=Child(ivol,1);
// Name
		strcpy(namec,((TGeant3*)gMC)->VolName(-icvol));
		tmp = new char[4];
		strncpy(tmp,(char *) &namec, 4);
		volume = new AliDrawVolume(namec);
		volume->SetIdVolume(-icvol);
// Number of divisions
		Int_t jvo  = fZlq[fGclink->jvolum-ivol];
		Int_t jdiv = fZlq[jvo-1];
		Int_t ndiv = Int_t (fZq[jdiv+3]);
		volume->SetIdCopy(-ndiv);
// Link to mother
		volume->SetItem(item2);
		(*fVolumes)[newlevel]=volume;
		printf("\n volume %s %d %d %d", namec, icvol, nch, ndiv);
		newlevel++;
//
// Children by positioning
	    } else {
		Int_t nnew=0;
// Loop over children 
		for (Int_t j=0; j<nch; j++) 
		{
//
// Find out if this volume was already positioned and count copies 
		    Int_t icvol=Child(ivol,j+1);
		    icvol = TMath::Abs(icvol);
		    Bool_t inList=kFALSE;
		    for (Int_t k=0; k<nnew; k++) {
			if (icvol==
			    Volume(newlevel-k-1)->GetIdVolume()) 
			{
			    Volume(newlevel-k-1)->AddCopy();
			    inList=kTRUE;
			}
		    }
//
// New child
		    if (!inList) {
//
// Name
			strcpy(namec,((TGeant3*)gMC)->VolName(icvol));
			tmp = new char[4];
			strncpy(tmp,(char *) &namec, 4);
			volume = new AliDrawVolume(namec);
			volume->SetIdVolume(icvol);
			volume->SetIdCopy(1);
// Link to mother
			volume->SetItem(item2);
			(*fVolumes)[newlevel]=volume;
			printf("\n volume %s %d %d", namec, icvol, nch);
			newlevel++;
			nnew++;
		    }
		}
	    }
	}
// Move one level deaper
	nst=nlevel;
	nlevel=newlevel;
    }
}

void AliGeant3GeometryGUI::ReadMaterials()
{
//
// Puts media and material names into ComboBox and 
// collects material information
// 
    Float_t a, z, dens, radl, absl;
    a=z=dens=radl=absl=-1;
    Int_t npar=0;
    Int_t imat, isvol, ifield;
    Float_t fieldm, tmaxfd, stemax, deemax, epsil, stmin;
    Float_t par[50];
    Float_t ubuf[50];
    Int_t nwbuf;
    
    char natmed[21], namate[21];
//
// Loop over media
    Int_t nEntries=0;
    
    for(Int_t itm=1;itm<=fGcnum->ntmed;itm++) {
	Int_t jtm  = fZlq[fGclink->jtmed-itm];
	if (jtm > 0) {
	    nEntries++;
// 
// Get medium parameters
	    ((TGeant3*)(gMC))->Gftmed(itm, natmed, imat, isvol, ifield, fieldm, 
				      tmaxfd, stemax, deemax, epsil, stmin, ubuf, &nwbuf);
	    strncpy(natmed,(char*)&fZiq[jtm+1],20);
	    natmed[20]='\0';
//
// Create new medium object 
	    AliGUIMedium * medium = 
		new AliGUIMedium(itm, imat, natmed, isvol, ifield, fieldm, 
				      tmaxfd, stemax, deemax, epsil, stmin);
	    (*fMedia)[nEntries-1]=medium;
        { //Begin local scope for j
          for (Int_t j=1; j<=22; j++) {
            medium->SetPar(j,Cut(itm,j));
          }
        } //End local scope for j
        { //Begin local scope for j
          for (Int_t j=23; j<=26; j++) {
            medium->SetPar(j,Cut(itm,j+3));
          }
        } //End local scope for j
        { //Begin local scope for j
          for (Int_t j=27; j<=29; j++) {
            medium->SetPar(j,Cut(itm,j+4));
          }
        } //End local scope for j
//
// Add to ComboBox
	    fPanel->AddMedium(medium, nEntries);
//
// Associated material
	    imat =  Int_t (fZq[jtm+6]);
	    Int_t jma  =  Int_t (fZlq[fGclink->jmate-imat]);
//
// Get material parameters
	    ((TGeant3*)(gMC))->Gfmate (imat,namate,a,z,dens,radl,absl,par,npar);
	    strncpy(namate,(char *)&fZiq[jma+1],20);
	    namate[20]='\0';
//
// Create new material object
	    AliGUIMaterial * material = 
		new AliGUIMaterial(imat,namate,a,z,dens,radl,absl);
	    (*fMaterials)[nEntries-1]=material;
	    material->Dump();
//
// Add to combo box
	    fPanel->AddMaterial(material, nEntries);
	    gCurrentMaterial=material;
	}
    }
    fPanel->SetComboEntries(fMaterials);
    fPanel->SetMediaComboEntries(fMedia);
    fPanel->Update();
}

Int_t AliGeant3GeometryGUI::NChildren(Int_t idvol)
{
//
// Return number of children for volume idvol
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin = Int_t(fZq[jvo+3]);
    return nin;
}

Int_t AliGeant3GeometryGUI::Child(Int_t idvol, Int_t idc)
{
//
// Return GEANT id of child number idc of volume idvol
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin=idc;
    Int_t jin = fZlq[jvo-nin];
    Int_t numb =  Int_t (fZq[jin +3]);
    if (numb > 1) {
	return -Int_t(fZq[jin+2]);
    } else {
	return Int_t(fZq[jin+2]);
    }
}

Int_t AliGeant3GeometryGUI::Medium(Int_t idvol)
{
//
// Return medium number for volume idvol.
// If idvol is negative the volume results from a division.
    Int_t imed;
    if (idvol > 0) {
	Int_t jvo = fZlq[fGclink->jvolum-idvol];
	imed = Int_t(fZq[jvo+4]);
    } else {
	idvol=-idvol;
	Int_t jdiv = fZlq[fGclink->jvolum-idvol];
	Int_t ivin = Int_t ( fZq[jdiv+2]);
	Int_t jvin = fZlq[fGclink->jvolum-ivin];
	imed = Int_t (fZq[jvin+4]);
    }
    return imed;
}

Int_t AliGeant3GeometryGUI::Material(Int_t idvol)
{
// Return material number for volume idvol.
// If idvol is negative the volume results from a division.

    Int_t imed=Medium(idvol);
    Int_t jtm  = fZlq[fGclink->jtmed-imed];
    return Int_t (fZq[jtm+6]);
}


Float_t AliGeant3GeometryGUI::Cut(Int_t imed, Int_t icut)
{
// Return cuts icut for medium idmed 
// 
    Int_t jtm  = fZlq[fGclink->jtmed-imed];
//
//  Have the defaults been modified ??
    Int_t jtmn = fZlq[jtm];
    if (jtmn >0) {
	jtm=jtmn;
	
    } else {
	jtm=fGclink->jtmed;
    }
    
    return Float_t (fZq[jtm+icut]);
}


