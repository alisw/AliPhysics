///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Interface Class to the Geant3.21 MonteCarlo                              //
//  Author :    Rene Brun   08/12/98                                         //
//                                                                           //
//  A TPaveTree is a TPaveLabel  specialized to draw Geant GDtree            //
//                                                                           //
//  When a TPaveTree object is drawn in a canvas, one can:                   //
//    - double click to draw the spec                                        //
//    - select pop-up item DrawTree to draw its tree                         //
//    - select pop-up item DrawTreeParent to draw its parent tree            //
//                                                                           //
//Begin_Html
/*
<img src="picts/TPaveTreeClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <fstream.h>
#include <iostream.h>

#include <TROOT.h>
#include <TVirtualPad.h>
#include <Buttons.h>
#include "TPaveTree.h"
#include "TGeant3.h"


ClassImp(TPaveTree)

//_____________________________________________________________________________
TPaveTree::TPaveTree(): TPaveLabel()
{
  //
  // Default Constructor
  //
}

//_____________________________________________________________________________
TPaveTree::TPaveTree(Coord_t x1, Coord_t y1,Coord_t x2, Coord_t  y2, 
		     const Text_t *label)
  :TPaveLabel(x1,y1,x2,y2,label,"br")
{
  //
  // TPaveTree normal constructor
  //
  SetName(label);
}

//_____________________________________________________________________________
TPaveTree::~TPaveTree()
{
  //
  // Standard Destructor
  //
}

//_____________________________________________________________________________
TPaveTree::TPaveTree(const TPaveTree &PaveTree)
{
  //
  // Copy Constructor
  //
  ((TPaveTree&)PaveTree).Copy(*this);
}

//_____________________________________________________________________________
void TPaveTree::Copy(TObject &obj)
{
  //
  //  Copy this PaveTree to PaveTree
  //
  TPaveLabel::Copy(obj);
}

//_____________________________________________________________________________
void TPaveTree::DrawSpec()
{
  //
  //  Draw specs of the volume in this TPaveTree
  //
  gMC->DrawOneSpec(GetLabel());
}

//_____________________________________________________________________________
void TPaveTree::DrawTree(Int_t levmax, Int_t isel)
{
  //
  //  Draw tree of the volume in this TPaveTree
  //
  ((TGeant3*)gMC)->Gdtree(GetLabel(),levmax,isel);
}

//_____________________________________________________________________________
void TPaveTree::DrawTreeParent(Int_t levmax, Int_t isel)
{
  //
  //  Draw parent tree of the volume in this TPaveTree
  //
  ((TGeant3*)gMC)->GdtreeParent(GetLabel(),levmax,isel);
}

//_____________________________________________________________________________
void TPaveTree::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //
  // Process mouse events.
  // Invokes TPabeLabel::ExecuteEvent. In case of a double click
  // draw specs of volume corresponding to this TPaveTree
  //
  
  TPaveLabel::ExecuteEvent(event,px,py);
  
  if (event == kButton1Double) {
    //printf("TPaveTree::ExecuteEvent\n");
    gPad->SetCursor(kWatch);
    gMC->DrawOneSpec(GetLabel());
  }
}

//_____________________________________________________________________________
void TPaveTree::SavePrimitive(ofstream &out, Option_t *)
{
  //
  // Save primitive as a C++ statement(s) on output stream out
  //
  char quote = '"';
  out<<"   "<<endl;
  if (gROOT->ClassSaved(TPaveTree::Class())) {
    out<<"   ";
  } else {
    out<<"   TPaveTree *";
  }
  out<<"pvar = new TPaveTree("<<fX1<<","<<fY1<<","<<fX2<<","<<fY2
     <<","<<quote<<fLabel<<quote<<","<<quote<<fOption<<quote<<");"<<endl;
  
  SaveFillAttributes(out,"pvar",0,1001);
  SaveLineAttributes(out,"pvar",1,1,1);
  SaveTextAttributes(out,"pvar",22,0,1,62,0);
  
  out<<"   pvar->Draw();"<<endl;
}
