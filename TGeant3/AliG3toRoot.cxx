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
*/


#include "AliG3toRoot.h"
#include "AliG3Volume.h"
#include "AliG3Medium.h"
#include "AliG3Material.h"
#include "AliNode.h"
#include "AliRun.h"

#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>
#include <TFolder.h>
#include <TRotMatrix.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TBRIK.h>

ClassImp(AliG3toRoot)
    AliG3toRoot::AliG3toRoot() 
{
//  Constructor
    TFolder* top = gROOT->GetRootFolder ()->AddFolder ("G3toRoot", "Geometry Tree");
    fTopFolder=top->AddFolder("ALIC", "Top Volume");
    gROOT->GetListOfBrowsables ()->Add (top, "G3toRoot Folders");

//
    fVolumes   = new TClonesArray("AliG3Volume",1000);
    fMaterials = new TClonesArray("AliG3Material",1000);
    fMedia     = new TClonesArray("AliG3Medium",1000);
    fRotations = 0;
    fExpand    = 0;
    
//  Store local copy of zebra bank entries
    TGeant3 *geant3 = (TGeant3*) gMC;
    if (geant3) {
	fZlq=geant3->Lq();
	fZq=geant3->Q();
	fZiq=geant3->Iq();
	fGclink=geant3->Gclink();
	fGcnum=geant3->Gcnum();
    }
// Geometry
    fGeometry = new TGeometry("AliceGeom","Detailed Geometry");
}
void AliG3toRoot::G3toRoot()
{
    ReadGeometryTree();
    ReadMaterials();
    ReadRotations();
    ConvertToRootShapes();

//
    TFolder* topFolder = gROOT->GetRootFolder ()->AddFolder ("G3", "G3");
    TFolder* geoFolder = topFolder->AddFolder("G3toRoot Shapes", "  ");
    geoFolder->Add(fGeometry);
    gROOT->GetListOfBrowsables ()->Add(geoFolder, "G3toRoot Shapes");

    TFile* file = new TFile("g3toRootShapes.root", "recreate");
    geoFolder->Write();
    file->Close();

    file = new TFile("g3toRootFolders.root", "recreate");
    fTopFolder->Write();
    file->Close();

}

void AliG3toRoot::ReadGeometryTree()
{
//
// Copy zebra volume tree into ROOT LisTree
//
    char *vname;
    char /* *namec, */ *tmp;
    char namec[30];

    AliG3Volume  *volume;

    Int_t nst=0;
    Int_t nlevel=1;
    Int_t newlevel=nlevel;
    Float_t mpar[3] = {2000., 2000., 2000.};
    
    
    volume = new AliG3Volume("ALIC");
    volume->SetIdVolume(((TGeant3*)gMC)->VolId("ALIC"));
    volume->SetIdCopy(0);
    volume->SetShape(1);
    volume->SetPosition(0.,0.,0.);
    volume->SetParameters(3,mpar);
    volume->SetItem(fTopFolder);
    (*fVolumes)[0]=volume;
//
//  Loop over volumes for which information has been collected
    while(nlevel>nst) {
	for (Int_t i=nst; i<nlevel; i++) 
	{
	    TFolder *itemi, *item2;
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

	    TArrayF pos;
	    pos = Volume(i)->Position(0);
	    
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
		itemi=(TFolder*) Volume(i)->GetItem();
	    } else {
		itemi=fTopFolder;
	    }
	    Volume(i)->SetName(namec);
	    Volume(i)->SetIdMedium(imed);
	    Volume(i)->SetIdMaterial(imat);
	    
//
// Add volume to list tree
	    
	    if (i>0) {
		item2 = AddVolume(new AliG3Volume(*Volume(i)), itemi, namec);
	    } else {
		item2 = fTopFolder;
		fTopFolder->Add(new AliG3Volume(*Volume(i)));
	    }
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
		volume = new AliG3Volume(namec);
		volume->SetIdVolume(-icvol);

		Int_t jvo  = fZlq[fGclink->jvolum-ivol];
		Int_t jdiv = fZlq[jvo-1];
// Axis
		Int_t iaxe = Int_t (fZq[jdiv+1]);
// System volume number
//		Int_t ivin = Int_t (fZq[jdiv+2]);
// Number of divisions
		Int_t ndiv = Int_t (fZq[jdiv+3]);
// Start value
		Float_t startc =    fZq[jdiv+4];
// Step    
		Float_t step   =    fZq[jdiv+5];

		Int_t ish   = Volume(i)->Shape();
		Int_t npar  = Volume(i)->NParam();
// Copy here shape parameters from mother
// Expand divisions later (when needed)
		TArrayF apar;
		apar.Set(npar);
		Volume(i)->Parameters(0, apar);
		Float_t* par = new Float_t[npar];
		for (Int_t jj=0; jj<npar; jj++) par[jj]=apar.At(jj);
		volume->SetIdCopy(-ndiv);
		volume->SetDivision(ndiv, iaxe, startc, step);
		volume->SetParameters(npar, par);
		volume->SetPosition(0.,0.,0.);
		volume->SetShape(ish);
		volume->SetRotMatrix(0);
// Link to mother
		volume->SetItem(item2);
		(*fVolumes)[newlevel]=volume;
		newlevel++;
//
// Children by positioning
	    } else if (nch > 0) {
		Int_t nnew=0;
// Loop over children 
		for (Int_t j=0; j<nch; j++) 
		{
		    Int_t jvo  = fZlq[fGclink->jvolum-ivol];
		    Int_t jin  = fZlq[jvo-j-1];
		    Int_t ivin = Int_t (fZq[jin + 2 ]);
		    Int_t jvin = fZlq[fGclink->jvolum - ivin];
		    Int_t ish  = Int_t (fZq[jvin + 2]);
		    Int_t irot = Int_t(fZq[jin +4]);
		    Float_t x    =  fZq[jin +5];
		    Float_t y    =  fZq[jin +6];
		    Float_t z    =  fZq[jin +7];
		    Int_t ndata = Int_t(fZiq[jin-1]);
		    Int_t npar;
		    Float_t* par;
//
// Get shape parameters
//
// parp or normal positioning ?		    
		    if (ndata == 8) {
			npar  = Int_t (fZq[jvin+5]);
			par = new Float_t[npar];
			if (npar>0) {
			    for(Int_t jj=0; jj<npar; jj++) 
				par[jj] = fZq[jvin+7+jj];
			}
		    } else {
			npar =  Int_t (fZq[jin+9]);
			par = new Float_t[npar];
			for(Int_t jj=0;jj<npar;jj++) 
			    par[jj]=fZq[jin+10+jj];
		    }
//
// Find out if this volume was already positioned and count copies 
		    Int_t icvol=Child(ivol,j+1);
		    icvol = TMath::Abs(icvol);
		    Bool_t inList=kFALSE;
		    AliG3Volume* copy0=0;
		    
		    for (Int_t k=0; k<nnew; k++) {
			if (icvol==
			    Volume(newlevel-k-1)->GetIdVolume()) 
			{
			    copy0  = Volume(newlevel-k-1);
			    Volume(newlevel-k-1)->AddCopy();
			    inList=kTRUE;
			}
		    }
//
// New child
//
// Name
		    strcpy(namec,((TGeant3*)gMC)->VolName(icvol));
		    tmp = new char[4];
		    strncpy(tmp,(char *) &namec, 4);
		    volume = new AliG3Volume(namec);
		    volume->SetPosition(x,y,z);
		    volume->SetShape(ish);
		    volume->SetParameters(npar, par);
		    volume->SetRotMatrix(irot);
		    if (ndata != 8) volume->SetPosp(kTRUE);
		    volume->SetIdVolume(icvol);
		    volume->SetIdCopy(1);
// Link to mother
		    volume->SetItem(item2);
		    if (!inList) {
			(*fVolumes)[newlevel]=volume;
			newlevel++;
			nnew++;
		    } else {
			copy0->AddCopy(volume);
		    }
		}
	    }
	}
// Move one level deaper
	nst=nlevel;
	nlevel=newlevel;
    }
}

void AliG3toRoot::ReadMaterials()
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
	    ((TGeant3*)(gMC))->
		Gftmed(itm, natmed, imat, isvol, ifield, fieldm, 
		       tmaxfd, stemax, deemax, epsil, stmin, ubuf, &nwbuf);
	    strncpy(natmed,(char*)&fZiq[jtm+1],20);
	    natmed[20]='\0';
//
// Create new medium object 
	    AliG3Medium * medium = 
		new AliG3Medium(itm, imat, natmed, isvol, ifield, fieldm, 
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
	    AliG3Material* material = new AliG3Material(namate, " ", a,z,dens,radl,absl);
	    material->SetId(imat);
	    
	    (*fMaterials)[nEntries-1] = material;
	    material->Dump();
	}
    }
}

void AliG3toRoot::ReadRotations()
{
// Read rotation matrices
    Int_t nrot = 0;
    Int_t irm, jrm;
    char name[9];
    if(fGclink->jrotm) nrot = Int_t(fZiq[fGclink->jrotm-2]);
    if (nrot) {
	 fRotations = new TObjArray(nrot);
	 for(irm=1;irm<=nrot;irm++) {
	     jrm  = fZlq[fGclink->jrotm-irm];
	     sprintf(name, "R%d", irm);
	     (*fRotations)[irm-1] =
		 new TRotMatrix(name, "Rotation", fZq[jrm+11], fZq[jrm+12], fZq[jrm+13], 
				fZq[jrm+14], fZq[jrm+15], fZq[jrm+16]);
	 }
    }
}

Int_t AliG3toRoot::NChildren(Int_t idvol)
{
//
// Return number of children for volume idvol
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin = Int_t(fZq[jvo+3]);
    return nin;
}

Int_t AliG3toRoot::Child(Int_t idvol, Int_t idc)
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

Int_t AliG3toRoot::Medium(Int_t idvol)
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

Int_t AliG3toRoot::Material(Int_t idvol)
{
// Return material number for volume idvol.
// If idvol is negative the volume results from a division.

    Int_t imed=Medium(idvol);
    Int_t jtm  = fZlq[fGclink->jtmed-imed];
    return Int_t (fZq[jtm+6]);
}


Float_t AliG3toRoot::Cut(Int_t imed, Int_t icut)
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

TFolder* AliG3toRoot::AddVolume(TObject * obj, TFolder *parent, const char* name)
{
// Add item to the list tree
    TFolder* newFolder = parent->AddFolder(name, "volume");
    newFolder->Add(obj);
    return newFolder;
}

AliG3Volume* AliG3toRoot::Volume(Int_t id)
	
{
// Volume at id
    return (AliG3Volume *) (fVolumes->UncheckedAt(id));
}


void AliG3toRoot::ConvertToRootShapes(TFolder *item, AliNode** node, Int_t nNodes)
{
// Convert the geometry represented by TFolders into root shapes
//
    AliG3Volume* volume;
    TArrayF pos;
    Int_t irot, imat;
    Int_t npos = 0;
    Bool_t top=kFALSE;
    
    char nameV[20];
    char nameN[20];
    
    if (!item) {
	item = fTopFolder;
	top  = kTRUE;
    }
    if (nNodes == 0) nNodes = 1;
    
    TList* folders = (TList*) item->GetListOfFolders();
    TIter next(folders);
 
    npos = 0;
    volume = ((AliG3Volume *) item->FindObject(item->GetName()));
    Int_t ncopy = volume->NCopies();
    AliNode** newnode = new AliNode*[(ncopy+1)*nNodes];

    for (Int_t ino=0; ino<nNodes; ino++) {
		    
	pos    = volume->Position(0);
	irot   = volume->RotMatrix();
	imat   = volume->Material();
	
	sprintf(nameV,"v%s%d", volume->GetName(), npos);
	sprintf(nameN,"%s%d" , volume->GetName(), npos);
	
	
	imat   = volume->Material();
	Int_t nmat = fMaterials->GetEntriesFast();
	AliG3Material* mat=0;

	for (Int_t mati=0; mati<nmat; mati++) {
	    mat = (AliG3Material*) 
		(fMaterials->UncheckedAt(mati));
	    if ((mat->Id())==imat) break;
	}
	    

	volume->CreateTShape(nameV, mat);

	if (!top) {
	    node[ino]->cd();
	    newnode[npos] = new AliNode(nameN," ",nameV, pos[0], pos[1], pos[2]);
	    newnode[npos]->SetDivision(volume->Ndiv(), volume->Axis(), 
				       volume->Step(), volume->StartC());
	    if (irot > 0) newnode[npos]->SetMatrix((TRotMatrix*) (*fRotations)[irot-1]);
	    npos++;
		
	} else {
	    new TMaterial("void","Vacuum",0,0,0);  //Everything is void
	    TBRIK *brik   = new TBRIK("D_alice","alice volume",
				      "void",2000,2000,3000);
	    brik->SetVisibility(0);
	    newnode[npos] = new AliNode("alice","alice","D_alice");
	    npos++;
	}
	
	if (ncopy) {
	    for (Int_t icop = 0; icop < ncopy; icop++) {
		AliG3Volume* copy = volume->Copy(icop);
		
		sprintf(nameN,"%s%d" , volume->GetName(), icop+1);
		if (copy->Posp()) {
		    sprintf(nameV,"v%s%d", volume->GetName(), icop+1);
		    volume->CreateTShape(nameV, mat);
		}
		node[ino]->cd();
		pos    = copy->Position(0);
		irot   = copy->RotMatrix();
		newnode[npos] = new AliNode(nameN," ",nameV, pos[0], pos[1], pos[2]);
		newnode[npos]->SetDivision(volume->Ndiv(), volume->Axis(), 
					   volume->Step(), volume->StartC());
		if (irot >0) newnode[npos]->SetMatrix((TRotMatrix*) (*fRotations)[irot-1]);
		npos++;
	    } // copy loop
	} // copies exist
    } // node loop
    TObject* obj;
    
    while ((obj = next()))
    {
	if ((AliG3Volume*) obj == volume) continue;
	item = (TFolder*) obj;
	
	ConvertToRootShapes(item, newnode, npos);
    }

//  Expand divisions of demand
    if (fExpand) ExpandDivisions();
}

void AliG3toRoot::ExpandDivisions(AliNode* node)
{
// Expand volume divisions
//

    Bool_t top = kFALSE;
    
    if (!node) {
	node = (AliNode*) fGeometry->GetNode("alice");
	top = kTRUE;
    }
    if (!top && node->Ndiv() > 0 && node->Axis()>0) {
	
	node->ExpandDivisions();
	node->GetParent()->RecursiveRemove(node);
	ExpandDivisions();

    } else {
	TList* sons = node->GetListOfNodes();
	TIter next(sons);
	AliNode* son=0;

	while((son = (AliNode*)next())) {
	    ExpandDivisions(son);
	}  
    }
} 







