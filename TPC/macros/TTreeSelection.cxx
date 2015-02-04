/// \class TTreeSelection
///
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/macros/TTreeSelection.cxx++
///
/// TFile f("eveTree.root");
/// Tracks;
///
/// TTreeDraw draw("draw",Tracks);
/// draw.AddSelection("TPC clusters","Tr.fTPCncls>0");
/// draw.AddSelection("TRD clusters","Tr.fTRDncls>0");
/// ~~~

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TPolyMarker3D.h"
#include "TVectorD.h"
#include "TObjString.h"


class TTreePoint: public TNamed{
public:
  TTreePoint();
  TTreePoint(const char *alias, const char * px, const char * py, const char *pz, Int_t mColor, Int_t mSize, Int_t mType);
  //  TString GetSelection();
public:
  TString fPx;   ///< point X
  TString fPy;   ///< point Y
  TString fPz;   ///< point Z
  Int_t   fMColor; ///< color
  Int_t   fMSize;  ///< marker size
  Int_t   fMType;  ///< marker type
  Bool_t  fIsOn;   ///< is On
  ClassDef(TTreePoint,1)
};


class TTreeCutAtom: public TNamed{
public:
  enum ExprType { kBool=0, kInt=1, kRange=2};
  TTreeCutAtom();
  TTreeCutAtom(const char *alias, const char *expr, ExprType type, Double_t val0, Double_t val1=0);
  TString GetSelection();
public:
  Double_t   fInt0;    ///< interval  value 0
  Double_t   fInt1;    ///< interval  value 1
  Double_t   fVal0;    ///< selection value 0
  Double_t   fVal1;    ///< selection value 1
  ExprType   fType;    ///< selection type
  ClassDef(TTreeCutAtom,1)
};

class TTreeDraw: public TNamed{
public:
  TTreeDraw(const char *name, TTree * tree);
  ~TTreeDraw(){;}
  TString  MakeSelection();
  void    AddSelectionRange(const char *alias, const char*expr, Float_t min, Float_t max);
  void    AddSelection(const char *alias, const char*expr);
  void    AddDraw(const char *alias, const char * px, const char * py, const char *pz, Int_t mColor, Int_t mSize, Int_t mType);
  //  TGCompositeFrame * MakeFrame();
public:
  TTree     * fTree;          ///< tree
  TObjArray   fCutAtoms;      ///< array of axpressions
  TObjArray   fDraws;         ///< array of draw experssions
private:
  TTreeDraw();
  ClassDef(TTreeDraw,1)
};

ClassImp(TTreePoint)
ClassImp(TTreeCutAtom)
/// \cond CLASSIMP
ClassImp(TTreeDraw)
/// \endcond

TTreePoint::TTreePoint():
  TNamed(),
  fPx(),   // point X
  fPy(),   // point Y
  fPz(),   // point Z
  fMColor(1), //color
  fMSize(1),  //marker size
  fMType(22),  //marker type
  fIsOn(kTRUE)   //is On
{
}

TTreePoint::TTreePoint(const char *alias, const char * px, const char * py, const char *pz, Int_t mColor, Int_t mSize, Int_t mType)
:
  TNamed(alias,alias),
  fPx(px),   // point X
  fPy(py),   // point Y
  fPz(pz),   // point Z
  fMColor(mColor), //color
  fMSize(mSize),  //marker size
  fMType(mType),  //marker type
  fIsOn(kTRUE)   //is On
{
}


TTreeCutAtom::TTreeCutAtom():
    TNamed(),
    fVal0(0),
    fVal1(0),
    fType(kBool)
{
  //
  //
  //
}

TTreeCutAtom::TTreeCutAtom(const char *alias, const char *expr, ExprType type, Double_t val0, Double_t val1):
  TNamed(alias,expr),
  fInt0(val0),
  fInt1(val1),
  fVal0(val0),
  fVal1(val1),
  fType(type)
{
  //
  //
  //
}

TString TTreeCutAtom::GetSelection(){
  //
  //
  //
  TString str;
  char  command[1000];
  if (fType==kBool) str = fTitle;
  if (fType==kInt){
    sprintf(command,"(%s==%d)",GetTitle(), TMath::Nint(fVal0));
    str = command;
  }
  if (fType==kRange){
    sprintf(command,"((%s>%f) &&(%s<%f))",GetTitle(), fVal0,GetTitle(),fVal1);
    str = command;
  }
  return str;
}


TTreeDraw::TTreeDraw():
    TNamed(),
    fTree(0),
    fCutAtoms(0),
    fDraws(0)
{
}

TTreeDraw::TTreeDraw(const char *name, TTree * tree):
  TNamed(name,name),
  fTree(tree),
  fCutAtoms(100),
  fDraws(100)
{
}

void    TTreeDraw::AddSelection(const char *alias, const char*expr){
  //
  // add string selection
  //
  TTreeCutAtom * atom = new TTreeCutAtom(alias,expr,TTreeCutAtom::kBool,0,1);
  fCutAtoms.AddLast(atom);
}

void    TTreeDraw::AddSelectionRange(const char *alias, const char*expr, Float_t min, Float_t max){
  //
  //
  //
  TTreeCutAtom * atom = new TTreeCutAtom(alias,expr,TTreeCutAtom::kRange,min,max);
  fCutAtoms.AddLast(atom);
}

TString  TTreeDraw::MakeSelection(){
  //
  // Make selection string
  //
  TString res;
  for (Int_t i=0; i<fCutAtoms.GetEntries(); i++){
    TTreeCutAtom * atom = (TTreeCutAtom*)fCutAtoms.At(i);
    if (!atom) continue;
    if (res.Length()>0) res+="&&";
    res+=atom->GetSelection();
  }
  return res;
}
