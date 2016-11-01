#ifndef ALITREETRENDING_H
#define ALITREETRENDING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
// class for visualization of trending trees  //
//
////////////////////////////////////////////////

 
class TTree;
class TObjArray;
class TNamed;
class TStyle;

class AliTreeTrending : public TNamed {

public:
  AliTreeTrending();
  AliTreeTrending(const char *name, const char *title);
  ~AliTreeTrending(){;}
  void SetTree(TTree * tree) {fTree=tree;};
  void AddUserDescription(TNamed *description);  //
  Bool_t  InitSummaryTrending(TString statusDescription[3], Float_t descriptionSize=0.015);
  void SetDefaultStyle();  // own style to be used - not yet
  void AppendStatusPad(Float_t padratio, Float_t bottomMargin, Float_t rightMargin);
public:
  TObjArray *  GetTexDescription(TLatex *latex);  // Currently only standard variables
  TTree     *  fTree;              // working tree with friends
  TObjArray *  fUserDescription;   // user defined description
  TObjArray *  fLatexDescription;  // description of process  (time/aliroot/root/user)
  TCanvas   *  fWorkingCanvas;     // default canvas
  TObjArray *  fStatusGraph;       // fStatusGraph
  //TVectorF   fPadMargin;         // pad margin
  TStyle    *  fDrawStyle;         // TreeTrending owner of drawing style - currently gStyle used
  ClassDef(AliTreeTrending,1)       
};


#endif

