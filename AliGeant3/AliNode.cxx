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
Revision 1.1  2001/07/09 11:44:01  morsch
Node class derived from TNode to allow expansion of volume divisions.

*/


#include "AliNode.h"
#include "TShape.h"
#include "TTUBE.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTRAP.h"
#include "TTUBS.h"
#include "TCONE.h"
#include "TCONS.h"
#include "TSPHE.h"
#include "TPARA.h"
#include "TPGON.h"
#include "TPCON.h"
#include "TTUBS.h"
#include "TELTU.h"
#include "THYPE.h"
#include "TGTRA.h"
#include "TCTUB.h"

ClassImp(AliNode)
    AliNode::AliNode(const char* name, const char* title, const char* shapename,
		     Double_t x, Double_t y, Double_t z, const char* matrixname,
		     Option_t* option) : 
	TNode(name, title, shapename, x, y, z, matrixname, option)
{
    fNDivision   = -1;
    fAxis   =  0;
    fStartC = 0.;
    fStep   = 0.;
}

  AliNode::AliNode(const char* name, const char* title, TShape* shape,
		     Double_t x, Double_t y, Double_t z, TRotMatrix* matrix,
		     Option_t* option) : 
	TNode(name, title, shape, x, y, z, matrix, option)
{
    fNDivision   = -1;
    fAxis   =  0;
    fStartC = 0.;
    fStep   = 0.;
}

void   AliNode::SetDivision(Int_t ndiv, Int_t axis, Float_t start, Float_t step)
{
    fNDivision   = ndiv;
    fAxis   = axis;
    fStartC = start;
    fStep   = step;
}

void AliNode::ExpandDivisions()
{
    Int_t i;
    char vName[20];
    char nName[20];
    
    char tmp[4];
    AliNode* node;
    TShape*  parent =  fParent->GetShape();
    TShape*  newsh;
    
    strcpy(tmp, parent->GetTitle());
    Int_t ndiv = fNDivision;

//                                      TUBE
    
    if (strcmp(tmp, "TUBE")==0) {
	TTUBE * shape = (TTUBE*) parent;

	Float_t dZ   = shape->GetDz();
	Float_t rMin = shape->GetRmin();
	Float_t rMax = shape->GetRmax();
	
	if (fAxis == 1) {
//  radial-division
	    Float_t dr = (rMax-rMin)/Float_t(fNDivision);
	    Float_t r1, r2;
	    for (i=0; i<ndiv; i++) {
		r1 = rMin+Float_t(i)*dr;
		r2 = r1+dr;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TTUBE(vName, "TUBE", "void", r1, r2, dZ);
		fParent->cd();
		node = new AliNode(nName,"", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
    
	} else if (fAxis == 2) {
//  phi-division
	    Float_t dPhi = 360./Float_t(fNDivision);
	    Float_t phi1, phi2;
	    for (i=0; i<ndiv; i++) {
		phi1 = Float_t(i)*dPhi;
		phi2 = phi1+dPhi;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TTUBS(vName, "TUBS", "void", rMin, rMax, dZ, phi1, phi2);
		fParent->cd();
		node = new AliNode(nName, "", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else {
//  z-division
	    Float_t delZ = dZ/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TTUBE(vName, "TUBE", "void", rMin, rMax, delZ);
		fParent->cd();
		Float_t zpos = -dZ+delZ*(2.*Float_t(i)+1.);
		node = new AliNode(nName, "",newsh, 0., 0., zpos);
		node->AddSons(fNodes);
		cd();
	    }
	}
//
//                            TUBS
//
    } else if (strcmp(tmp, "TUBS")==0) {
	TTUBS * shape = (TTUBS*) parent;
	Float_t dZ   = shape->GetDz();
	Float_t rMin = shape->GetRmin();
	Float_t rMax = shape->GetRmax();
	Float_t phi1 = shape->GetPhi1();
	Float_t phi2 = shape->GetPhi2();
	
	if (fAxis == 1) {
//  radial-division
	    Float_t dr = (rMax-rMin)/Float_t(fNDivision);
	    Float_t r1, r2;
	    Int_t ndiv = fNDivision;
	    for (i=0; i<ndiv; i++) {
		r1 = rMin+Float_t(i)*dr;
		r2 = r1+dr;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TTUBS(vName, "TUBS", "void", r1, r2, dZ, phi1, phi2);
		fParent->cd();
		node = new AliNode(nName,"", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	    
	} else if (fAxis == 2) {
//  phi-division
	    Float_t dPhi = (phi2-phi1)/Float_t(fNDivision);
	    Float_t nphi1, nphi2;
	    
	    for (i=0; i<fNDivision; i++) {
		nphi1 = phi1+Float_t(i)*dPhi;
		nphi2 = nphi1+dPhi;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);

		newsh = new TTUBS(vName, "TUBS", "void", rMin, rMax, dZ, nphi1, nphi2);
		fParent->cd();
		node = new AliNode(nName, "", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else {
//  z-division
	    Float_t delZ = dZ/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TTUBS(vName, "TUBS", "void", rMin, rMax, delZ, phi1, phi2);
		fParent->cd();
		Float_t zpos = -dZ+delZ*(2.*Float_t(i)+1.);
		node = new AliNode(nName, "",newsh, 0., 0., zpos);
		node->AddSons(fNodes);
		cd();
	    }
	}
	
    } else if (strcmp(tmp, "CONE")==0) {
	TCONE * shape = (TCONE*) parent;

	Float_t dZ   = shape->GetDz();
	Float_t rMin1 = shape->GetRmin();
	Float_t rMax1 = shape->GetRmax();
	Float_t rMin2 = shape->GetRmin2();
	Float_t rMax2 = shape->GetRmax2();
	
	if (fAxis == 1) {
//  radial-division
	    Float_t dr1 = (rMax1-rMin1)/Float_t(fNDivision);
	    Float_t dr2 = (rMax2-rMin2)/Float_t(fNDivision);
	    Float_t r11, r12, r21, r22;
	    for (i=0; i<ndiv; i++) {
		r11 = rMin1+Float_t(i)*dr1;
		r12 = r11+dr1;
		r21 = rMin2+Float_t(i)*dr2;
		r22 = r21+dr2;

		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TCONE(vName, "CONE", "void", dZ, r11, r12, r21, r22);
		fParent->cd();
		node = new AliNode(nName,"", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else if (fAxis == 2) {
//  phi-division
	    Float_t dPhi = 360./Float_t(fNDivision);
	    Float_t phi1, phi2;
	    for (i=0; i<ndiv; i++) {
		phi1 = Float_t(i)*dPhi;
		phi2 = phi1+dPhi;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TCONS(vName, "CONS", "void", dZ, rMin1, rMax1, 
			  rMin2, rMax2, phi1, phi2);
		fParent->cd();
		node = new AliNode(nName, "",newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else {
//  z-division
	    Float_t delZ = dZ/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TCONE(vName, "CONE", "void", delZ, rMin1, rMax1, rMin2, rMax2);
		fParent->cd();
		Float_t zpos = -dZ+delZ*(2.*Float_t(i)+1.);
		node = new AliNode(nName, "",newsh, 0., 0., zpos);
		node->AddSons(fNodes);
		cd();
	    }
	}
	
    } else if (strcmp(tmp, "CONS")==0) {
	TCONS * shape = (TCONS*) parent;
	Float_t dZ   = shape->GetDz();
	Float_t rMin1 = shape->GetRmin();
	Float_t rMax1 = shape->GetRmax();
	Float_t rMin2 = shape->GetRmin2();
	Float_t rMax2 = shape->GetRmax2();
	Float_t phi1 = shape->GetPhi1();
	Float_t phi2 = shape->GetPhi2();
	if (fAxis == 1) {
//  radial-division
	    Float_t dr1 = (rMax1-rMin1)/Float_t(fNDivision);
	    Float_t dr2 = (rMax2-rMin2)/Float_t(fNDivision);
	    Float_t r11, r12, r21, r22;
	    for (i=0; i<ndiv; i++) {
		r11 = rMin1+Float_t(i)*dr1;
		r12 = r11+dr1;
		r21 = rMin2+Float_t(i)*dr2;
		r22 = r21+dr2;

		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TCONS(vName, "CONS", "void", dZ, r11, r12, r21, r22, phi1, phi2);
		fParent->cd();
		node = new AliNode(nName,"", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	    
	} else if (fAxis == 2) {
//  phi-division
	    Float_t dPhi = (phi2-phi1)/Float_t(fNDivision);
	    Float_t nphi1, nphi2;
	    
	    for (i=0; i<fNDivision; i++) {
		nphi1 = phi1+Float_t(i)*dPhi;
		nphi2 = nphi1+dPhi;
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);

		newsh = new TCONS(vName, "CONS", "void", dZ, rMin1, rMax1, rMin2, rMax2, nphi1, nphi2);
		fParent->cd();
		node = new AliNode(nName, "", newsh, 0., 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else {
//  z-division
	    Float_t delZ = dZ/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TCONS(vName, "CONS", "void", delZ, rMin1, rMax1, rMin2, rMax2, phi1, phi2);
		fParent->cd();
		Float_t zpos = -dZ+delZ*(2.*Float_t(i)+1.);
		node = new AliNode(nName,"",newsh, 0., 0., zpos);
		node->AddSons(fNodes);
		cd();
	    }
	}
    } else if (strcmp(tmp, "BRIK")==0) {
//
//                          BRIK
//
	TBRIK * shape = (TBRIK*) parent;
	Float_t dX   = shape->GetDx();
	Float_t dY   = shape->GetDy();
	Float_t dZ   = shape->GetDz();

	if (fAxis == 1) {
// division in x
	    Float_t delX = dX/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TBRIK(vName, "BRIK", "void", delX, dY, dZ);
		fParent->cd();
		Float_t xpos = -dX+delX*(2.*Float_t(i)+1.);
		node = new AliNode(nName,"",newsh, xpos, 0., 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else if (fAxis == 2) {
// division in y
	    Float_t delY = dY/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TBRIK(vName, "BRIK", "void", dX, delY, dZ);
		fParent->cd();
		Float_t ypos = -dY+delY*(2.*Float_t(i)+1.);
		node = new AliNode(nName,"",newsh, 0., ypos, 0.);
		node->AddSons(fNodes);
		cd();
	    }
	} else {
// division in z
	    Float_t delZ = dZ/Float_t(fNDivision);
	    for (i=0; i<ndiv; i++) {
		sprintf(vName, "%sD%d", fShape->GetName(), i);
		sprintf(nName, "%sD%d", GetName(), i);
		newsh = new TBRIK(vName, "BRIK", "void", dX, dY, delZ);
		fParent->cd();
		Float_t zpos = -dZ+delZ*(2.*Float_t(i)+1.);
		node = new AliNode(nName,"",newsh, 0., 0., zpos);
		node->AddSons(fNodes);
		cd();
	    }
	}
    }
}

void AliNode::AddSons(TList* list)
{
    if (!list) return;
    if (!fNodes) fNodes = new TList();
    
    TIter next(list);
    AliNode* node;
    
    while((node = (AliNode*)next())) {
	AliNode* newNode = new AliNode(*node, this);
	fNodes->Add(newNode);
	newNode->SetParent(this);
	newNode->AddSons(node->GetListOfNodes());
    }
}


AliNode::AliNode(const AliNode &node, AliNode* parent)
{
    fNDivision   = node.Ndiv();
    fAxis        = node.Axis();
    fStartC      = node.StartC();
    fStep        = node.Step();
    fName        = node.GetName();
    fTitle       = node.GetTitle();
    fX           = node.GetX();
    fY           = node.GetY();
    fZ           = node.GetZ();
    fMatrix      = node.GetMatrix();
    fNodes       = 0;
    fParent      = parent;
    
    
    if (fNDivision > 0) {
	fShape = new TShape(*node.GetShape());
    } else {
	fShape = (TShape*) node.GetShape()->Clone();
    }
    
    
//    fShape = (TShape*) shape->Clone();
/*
    char tmp[4];        
    strcpy(tmp, shape->ClassName());

    if (strcmp(tmp, "TTUBE")==0) 
	fShape =  (TTUBE*)shape->Clone();
    else if (strcmp(tmp, "TBRIK")==0) 
	fShape =  (TBRIK*)shape->Clone();
    else if (strcmp(tmp, "TCONE")==0) 
	fShape =  (TCONE*)shape->Clone();
    else if (strcmp(tmp, "TCONS")==0) 
	fShape =  (TCONS*)shape->Clone();
    else if (strcmp(tmp, "TTUBS")==0) 
	fShape =  (TTUBS*)shape->Clone();
    else if (strcmp(tmp, "TTRAP")==0) 
	fShape =  (TTRAP*)shape->Clone();
    else if (strcmp(tmp, "TTRD1")==0) 
	fShape =  (TTRD1*)shape->Clone();
    else if (strcmp(tmp, "TTRD2")==0)
	fShape =  (TTRD2*)shape->Clone();
    else if (strcmp(tmp, "TSPHE")==0) 
	fShape =  (TSPHE*)shape->Clone();
    else if (strcmp(tmp, "TPGON")==0) 
	fShape =  (TPGON*)shape->Clone();
    else if (strcmp(tmp, "TPCON")==0) 
	fShape =  (TPCON*)shape->Clone();
*/
}

void AliNode::AddSon(AliNode* node)
{
    fNodes->Add(node);
}















