#ifndef ALIDNDPTEVENTCUTS_H
#define ALIDNDPTEVENTCUTS_H

//------------------------------------------------------------------------------
// Class to keep event selection cuts for dNdPt analysis. 
// 
// Author: J.Otwinowski 01/11/2008 
//------------------------------------------------------------------------------

#include "AliAnalysisCuts.h"

class AliESDEvent;
class AliESDVertex;
class AliMCEvent;
class AliHeader;
class AliGenEventHeader;

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"

class AlidNdPtEventCuts : public AliAnalysisCuts
{
public:
  AlidNdPtEventCuts(const Char_t* name ="AlidNdPtEventCuts", const Char_t *title ="");
  virtual ~AlidNdPtEventCuts(); 
 
  // setters 
  void SetTriggerRequired(const Bool_t bFlag=kTRUE)  {fTriggerRequired=bFlag;}
  void SetRecVertexRequired(const Bool_t bFlag=kTRUE)  {fRecVertexRequired=bFlag;}
  void SetEventProcessType(AliPWG0Helper::MCProcessType type=AliPWG0Helper::kInvalidProcess)  {fEventProcessType=type;}
  void SetNContributorsRange(const Float_t min=0.,const Float_t max=1e99) {fMinNContributors=min; fMaxNContributors=max;}
  void SetMaxR(const Float_t max=1e99) {fMaxR=max;}
  void SetZvRange(const Float_t min=-1e99, const Float_t max=1e99) {fMinZv=min; fMaxZv=max;}

  void SetMeanXYZv(const Float_t xv=0.0, const Float_t yv=0.0, const Float_t zv=0.0) {
    fMeanXv = xv; fMeanYv = yv; fMeanZv = zv;
  }

  void SetSigmaMeanXYZv(const Float_t sxv=1.0, const Float_t syv=1.0, const Float_t szv=10.0) {
    fSigmaMeanXv = sxv; fSigmaMeanYv = syv; fSigmaMeanZv = szv;
  }


  void SetRedoTPCVertex(const Bool_t redo = kTRUE) {fRedoTPCVertex = redo;}
  void SetUseBeamSpotConstraint(const Bool_t useConstr = kTRUE) {fUseBeamSpotConstraint = useConstr;}
  void SetEventSelectedRequired(const Bool_t evtSel = kTRUE) {fEventSelectedRequired = evtSel;} 


  // getters 
  Bool_t  IsEventSelectedRequired() const {return fEventSelectedRequired;}
  Bool_t  IsTriggerRequired() const {return fTriggerRequired;}
  Bool_t  IsRecVertexRequired() const {return fRecVertexRequired;}
  Int_t   GetEventProcessType() const {return fEventProcessType;}  
  Float_t GetMinNContributors() const {return fMinNContributors;}
  Float_t GetMaxNContributors() const {return fMaxNContributors;}
  Float_t GetMaxR() const {return fMaxR;}
  Float_t GetMinZv() const {return fMinZv;}
  Float_t GetMaxZv() const {return fMaxZv;}

  Float_t GetMeanXv() const {return fMeanXv;}
  Float_t GetMeanYv() const {return fMeanYv;}
  Float_t GetMeanZv() const {return fMeanZv;}

  Float_t GetSigmaMeanXv() const {return fSigmaMeanXv;}
  Float_t GetSigmaMeanYv() const {return fSigmaMeanYv;}
  Float_t GetSigmaMeanZv() const {return fSigmaMeanZv;}
 
  Bool_t IsRedoTPCVertex() const {return fRedoTPCVertex;}
  Bool_t IsUseBeamSpotConstraint() const {return fUseBeamSpotConstraint;}


  // cuts init function
  void Init();

  // check MC tracks
  Bool_t IsSelected(TObject *) {return kTRUE;}
  Bool_t IsSelected(TList *) {return kTRUE;}

  // accept event
  Bool_t AcceptEvent(AliESDEvent *event=0, AliMCEvent *mcEvent=0, const AliESDVertex *vtx=0);
  Bool_t AcceptMCEvent(AliMCEvent *mcEvent=0);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

private:
  Bool_t fTriggerRequired; // trigger required  
  Bool_t fRecVertexRequired; // reconstructed event vertex required  
  Int_t fEventProcessType;   // select MC event process type (ND, SD, DD)
  Float_t fMinNContributors; // min. number of contributing vertex tracks
  Float_t fMaxNContributors; // max. number of contributing vertex tracks
  Float_t fMaxR;             // max. vertex radii (R = sqrt(Xv^2+Yv^2) 
  Float_t fMinZv;            // min. Zv vertex
  Float_t fMaxZv;            // max. Zv vertex

  // interaction spot constraint
  Float_t fMeanXv; // mean Xv position
  Float_t fMeanYv; // mean Yv position
  Float_t fMeanZv; // mean Zv position

  Float_t fSigmaMeanXv; // sigma mean Xv position 
  Float_t fSigmaMeanYv; // sigma mean Yv position
  Float_t fSigmaMeanZv; // sigma mean Zv position
 
  Bool_t fRedoTPCVertex;         // redo vertex
  Bool_t fUseBeamSpotConstraint; // use beam spot contraints  

  Bool_t fEventSelectedRequired; // event with at least one track (pT>0.5 GeV, |eta|<0.8) required

  AlidNdPtEventCuts(const AlidNdPtEventCuts&); // not implemented
  AlidNdPtEventCuts& operator=(const AlidNdPtEventCuts&); // not implemented

  ClassDef(AlidNdPtEventCuts, 2)
};

#endif // ALIDNDPTEVENTCUTS_H
