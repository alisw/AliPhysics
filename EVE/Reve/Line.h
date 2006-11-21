// $Header$

#ifndef REVE_Line_H
#define REVE_Line_H

#include <Reve/Reve.h>
#include <Reve/PointSet.h>

#include <TAttLine.h>

namespace Reve {

class Line : public PointSet,
	     public TAttLine
{
  friend class LineEditor;
  friend class LineGL;

private:
  Line(const Line&);            // Not implemented
  Line& operator=(const Line&); // Not implemented

protected:
  Bool_t  fRnrLine;
  Bool_t  fRnrPoints;

public:
  Line(Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  Line(const Text_t* name, Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  Line(const Text_t* name, TTree* tree, TreeVarType_e tv_type=TVT_XYZ);
  virtual ~Line();

  virtual void SetMarkerColor(Color_t col)
  { TAttMarker::SetMarkerColor(col); }
  virtual void SetLineColor(Color_t col)
  { SetMainColor(col); }

  Bool_t GetRnrLine() const   { return fRnrLine;   }
  void SetRnrLine(Bool_t r)   { fRnrLine = r;      }
  Bool_t GetRnrPoints() const { return fRnrPoints; }
  void SetRnrPoints(Bool_t r) { fRnrPoints = r;    }

  ClassDef(Line, 1);
}; // endclass Line

}

#endif
