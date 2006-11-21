// $Header$

#include "Line.h"

using namespace Reve;

//______________________________________________________________________
// Line
//

ClassImp(Line)


Line::Line(Int_t n_points, TreeVarType_e tv_type) :
  PointSet(n_points, tv_type),
  fRnrLine   (kTRUE),
  fRnrPoints (kFALSE)
{
  fMainColorPtr = &fLineColor;
}

Line::Line(const Text_t* name, Int_t n_points, TreeVarType_e tv_type) :
  PointSet(name, n_points, tv_type),
  fRnrLine   (kTRUE),
  fRnrPoints (kFALSE)
{
  fMainColorPtr = &fLineColor;
}

Line::Line(const Text_t* name, TTree* tree, TreeVarType_e tv_type) :
  PointSet(name, tree, tv_type),
  fRnrLine   (kTRUE),
  fRnrPoints (kFALSE)
{
  fMainColorPtr = &fLineColor;
}

Line::~Line()
{}
