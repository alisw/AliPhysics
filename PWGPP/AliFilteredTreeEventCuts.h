#ifndef ALIFILTEREDTREEEVENTCUTS_H
#define ALIFILTEREDTREEEVENTCUTS_H

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
class AliStack;

class AliFilteredTreeEventCuts : public AliAnalysisCuts
{
public:
  AliFilteredTreeEventCuts(const Char_t* name ="AliFilteredTreeEventCuts", const Char_t *title ="");
  virtual ~AliFilteredTreeEventCuts(); 

  //TODO: copied from AliPWG0Helper, find a more central place for these
  enum AnalysisMode { kInvalid = -1, kSPD = 0x1, kTPC = 0x2, kTPCITS = 0x4, kFieldOn = 0x8, kSPDOnlyL0 = 0x10, kTPCSPD = 0x20};
  enum MCProcessType { kInvalidProcess = -1, kND = 0x1, kDD = 0x2, kSD = 0x4, kOnePart = 0x8 };
  enum DiffTreatment { kMCFlags = 0, kUA5Cuts = 1, kE710Cuts, kALICEHadronLevel };
  
  // setters 
  void SetTriggerRequired(Bool_t bFlag=kTRUE)  {fTriggerRequired=bFlag;}
  void SetRecVertexRequired(Bool_t bFlag=kTRUE)  {fRecVertexRequired=bFlag;}
  void SetEventProcessType(MCProcessType type=kInvalidProcess)  {fEventProcessType=type;}
  void SetNContributorsRange(Float_t min=0.,Float_t max=1e99) {fMinNContributors=min; fMaxNContributors=max;}
  void SetMaxR(Float_t max=1e99) {fMaxR=max;}
  void SetZvRange(Float_t min=-1e99, Float_t max=1e99) {fMinZv=min; fMaxZv=max;}

  void SetMeanXYZv(Float_t xv=0.0, Float_t yv=0.0, Float_t zv=0.0) {
    fMeanXv = xv; fMeanYv = yv; fMeanZv = zv;
  }

  void SetSigmaMeanXYZv(Float_t sxv=1.0, Float_t syv=1.0, Float_t szv=10.0) {
    fSigmaMeanXv = sxv; fSigmaMeanYv = syv; fSigmaMeanZv = szv;
  }


  void SetRedoTPCVertex(Bool_t redo = kTRUE) {fRedoTPCVertex = redo;}
  void SetUseBeamSpotConstraint(Bool_t useConstr = kTRUE) {fUseBeamSpotConstraint = useConstr;}
  void SetEventSelectedRequired(Bool_t evtSel = kTRUE) {fEventSelectedRequired = evtSel;} 


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

  //statics copied from AliPWG0Helper
  static MCProcessType GetEventProcessType(AliESDEvent* esd, AliHeader* header, AliStack* stack, DiffTreatment diffTreatment);
  static MCProcessType GetEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
  static MCProcessType GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
  static MCProcessType GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
  static Bool_t IsHadronLevelSingleDiffractive(AliStack* stack, Float_t cms, Float_t xiMin, Float_t xiMax);
  static Double_t Rapidity(Double_t pt, Double_t pz, Double_t m);
protected:
  static Int_t fgLastProcessType;    // stores the raw value of the last process type extracted
  //
 
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

  AliFilteredTreeEventCuts(const AliFilteredTreeEventCuts&); // not implemented
  AliFilteredTreeEventCuts& operator=(const AliFilteredTreeEventCuts&); // not implemented

  ClassDef(AliFilteredTreeEventCuts, 1)
};

#endif // ALIFILTEREDTREEEVENTCUTS_H
