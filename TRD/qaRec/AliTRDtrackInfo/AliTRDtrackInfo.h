#ifndef ALITRDTRACKINFO_H
#define ALITRDTRACKINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackInfo.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef Root_TObject
#include "TObject.h"
#endif

class AliTRDseedV1;
class AliTRDtrackV1;
class AliTrackReference;
class AliExternalTrackParam;
class AliTRDtrackInfo : public TObject{
public:
  AliTRDtrackInfo();
  AliTRDtrackInfo(Int_t pdg);
  AliTRDtrackInfo(const AliTRDtrackInfo &);
  ~AliTRDtrackInfo();
  
//  void               Clear(const Option_t *){}
  void               Delete(const Option_t *);
  
  AliTRDtrackInfo&   operator=(const AliTRDtrackInfo &);
  
  void               AddTrackRef(const AliTrackReference *trackRef);
  
  Int_t              GetTrackId() { return fId;}
  Int_t              GetNumberOfClusters() const;
  Int_t              GetNumberOfClustersRefit() const {return fNClusters;}
  Int_t              GetNTracklets() const;
  Int_t              GetNTrackRefs() const {return fNTrackRefs;} 
  Int_t              GetLabel() const { return fLabel; }
  Int_t              GetPDG() const { return fPDG; }
  ULong_t            GetStatus() const {return fStatus;}
  AliTRDseedV1*      GetTracklet(Int_t entry) const;
  AliTRDtrackV1 *	 	 GetTRDtrack() const { return fTRDtrack; }
  AliTrackReference* GetTrackRef(Int_t entry) const;
  AliExternalTrackParam* GetOuterParam() const {return fOP;}

  Bool_t             IsCurved() const {return TestBit(kCurv);}
  Bool_t			 			 IsPrimary() const {return TestBit(kPrim);}
	Bool_t						 HasESDtrack() const{return ((fTRDtrack != 0x0) ||(fOP != 0));}
	Bool_t             HasMCinfo() const { return fNTrackRefs > 0; }

  void               SetCurved(Bool_t curv = kTRUE) {SetBit(kCurv, curv);}
  void               SetLabel(Int_t lab) { fLabel = lab; }
  void               SetNumberOfClustersRefit(Int_t n) {fNClusters = n;}
  void               SetPDG(Int_t pdg) { fPDG = pdg; }
  void				 			 SetPrimary(Bool_t prim = kTRUE) {SetBit(kPrim, prim);}
  void               SetOuterParam(const AliExternalTrackParam *op);
  void               SetStatus(ULong_t stat) {fStatus = stat;}
  void               SetTrackId(Int_t id) {fId = id;}
  void               SetTRDtrack(const AliTRDtrackV1 *track);
  
private:
  	enum{
  		kCurv = 14,
  		kPrim = 15
	};
  Int_t              fPDG;           	// particle code
  ULong_t            fStatus;        	// ESD track status
  Int_t              fId;            	// ESD track id
  Int_t              fLabel;         	// MC label  
  Int_t              fNClusters;     	// Numer of clusters from refit
  Int_t              fNTrackRefs;    	// number of track refs
  AliTrackReference  *fTrackRefs[12];	// no of track refs
  AliTRDtrackV1      *fTRDtrack; 	// no of tracklets
  AliExternalTrackParam *fOP;       	// outer param if no tracklets
 
  ClassDef(AliTRDtrackInfo, 1)          // TRD track info
};
#endif
