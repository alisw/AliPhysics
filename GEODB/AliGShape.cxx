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

/*
$Log$
*/

// -*- C++ -*-
// 
// 1998/10/22
// ---------------------------------------------------------------------------
//
// AliGShape Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------

#include <TView.h>
#include <TVirtualPad.h>
#include <iostream.h>
#include "AliGShape.h"
#include "AliGNode.h"
#include "AliGeometry.h"

//---------------------------------------------------------------------------

ClassImp(AliGShape)

AliGShape::AliGShape()
{
    /* Default Constructor */
    fName  = "";
    fTitle = "";
}

//---------------------------------------------------------------------------

AliGShape::AliGShape(Text_t* name, Text_t* title): TNamed(name, title), TAttLine(), TAttFill()
{
   /* Constructor */
}

//---------------------------------------------------------------------------

AliGShape::AliGShape( AliGShape* shape )
{
    /* Copy Constructor */
    fColor = shape->fColor;
}


//---------------------------------------------------------------------------

void AliGShape::Paint(Option_t *)
{
   // This method must be overridden by the real shape implementation.
   // AbstractMethod("Paint");
}

// ---------------------------------------------------------------------------

Int_t AliGShape::DistancetoPrimitive(Int_t, Int_t)
{
    TView *view = gPad->GetView();
    gPad->SetSelected(view);
    return 0;
}

// ---------------------------------------------------------------------------

void AliGShape::PaintShape(X3DBuffer *buff, Bool_t rangeView)
{
//*-*-*-*-*Paint 3-D shape in current pad with its current attributes*-*-*-*-*
//*-*      ==========================================================

    if (!buff) return;

    Float_t *point = &(buff->points[0]);
    for (Int_t j = 0; j < buff->numPoints; j++) {
      //           printf("Before %f %f %f\n",point[3*j],point[3*j+1],point[3*j+2]);
           gAliGeometry->Local2Master(&point[3*j],&point[3*j]);
	   //           printf("After %f %f %f\n",point[3*j],point[3*j+1],point[3*j+2]);
    }
          
    Float_t points[6], x0, y0, z0, x1, y1, z1;
    const Int_t kExpandView = 2;
    int i0;

    x0 = y0 = z0 = x1 = y1 = z1 = buff->points[0];

    TAttLine::Modify();  //Change line attributes only if necessary
    TAttFill::Modify();  //Change fill area attributes only if necessary

    for (Int_t i = 0; i < buff->numSegs; i++) {
        i0 = 3*buff->segs[3*i+1];
        points[0] = buff->points[i0++];
        points[1] = buff->points[i0++];
        points[2] = buff->points[i0];

        i0 = 3*buff->segs[3*i+2];
        points[3] = buff->points[i0++];
        points[4] = buff->points[i0++];
        points[5] = buff->points[i0];

        x0 = points[0] < x0 ? points[0] : x0;
        y0 = points[1] < y0 ? points[1] : y0;
        z0 = points[2] < z0 ? points[2] : z0;
        x1 = points[0] > x1 ? points[0] : x1;
        y1 = points[1] > y1 ? points[1] : y1;
        z1 = points[2] > z1 ? points[2] : z1;

        Float_t *ptpoints_0 = &points[0];
        Float_t *ptpoints_3 = &points[3];

	//printf("Painting from %f to %f\n",*ptpoints_0, *ptpoints_3);
	
        gPad->PaintLine3D(ptpoints_0, ptpoints_3);
//        gPad->PaintLine3D(&points[0], &points[3]);
    }
    
    TView *view = gPad->GetView();
    //cout << "x0,y0,x1,y1=" << x0 << y0 << x1 << y1 << endl;
    if (view->GetAutoRange()) view->SetRange(-500,-500,-500,500,500,500,kExpandView);
    
    // if (view->GetAutoRange()) view->SetRange(-1,-1,-1,1,1,1,kExpandView);
}



// ---------------------------------------------------------------------------

void AliGShape::SetPoints(Float_t *) {
    AbstractMethod("SetPoints(Float_t *buffer)");
}

// ---------------------------------------------------------------------------


