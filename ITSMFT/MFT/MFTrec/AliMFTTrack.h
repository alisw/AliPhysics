#ifndef AliMFTTrack_H
#define AliMFTTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTrack
/// \brief Description of an ALICE Standalone MFT track
///
///
/// Description of an ALICE Standalone MFT track
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date April 27th, 2015

#include "TObject.h"
class TObjArray;
class AliMFTCATrack;
class AliMFTTrackParam;

//=============================================================================================

class AliMFTTrack : public TObject {
	
public:
	
	AliMFTTrack();
	AliMFTTrack(AliMFTCATrack *catrack);
	virtual ~AliMFTTrack();
	AliMFTTrack (const AliMFTTrack& track); // copy constructor
	AliMFTTrack& operator=(const AliMFTTrack& track); // assignment operator
	
	/// return the minimum value of the function minimized by the fit
	Double_t GetChi2() const {return fChi2;}
	/// set the minimum value of the function minimized by the fit
	void SetChi2(Double_t chi2) { fChi2 = chi2;}
	Int_t    GetNDF() const;
	Double_t GetNormalizedChi2() const;
	
	/// return pointer to track parameters at vertex (can be 0x0)
	AliMFTTrackParam* GetTrackParamAtVertex() const {return fTrackParamAtVertex;}
	void SetTrackParamAtVertex(  AliMFTTrackParam* trackParam){};
	
	/// return the number of clusters attached to the track
	Int_t GetNClusters() const {return fTrackParamAtCluster ? fTrackParamAtCluster->GetEntriesFast() : 0;}
	
	
	/// return pointer to track found by Track Finder (includes clusters)
	AliMFTCATrack* GetCATrack() const {return fCATrack;}
	void SetCATrack(  AliMFTCATrack* track) {};
	
	void AddTrackParamAtCluster(const AliMFTTrackParam &trackParam);
	TObjArray*    GetTrackParamAtCluster() const;
	
	/// set the corresponding MC track number
	void  SetMCLabel(Int_t label) {fTrackID = label;}
	/// return the corresponding MC track number
	Int_t GetMCLabel() const {return fTrackID;}
	
	virtual void Print(Option_t* opt="") const;
	virtual void Clear(Option_t* opt="");
	
	void SetP(double val) {fP = val;};
	void SetPt(double val) {fPt = val;};
	void SetTheta(double val) {fTheta = val;};
	void SetPhi(double val) {fPhi = val;};
	
	Double_t P() {return fP ;};
	Double_t Pt() {return fPt ;};
	Double_t Theta() {return fTheta ;};
	Double_t Phi() {return fPhi ;};
	
	
	
protected:
	
	Double_t fChi2;                         ///<  Chi2 of the track
	mutable TObjArray* fTrackParamAtCluster; //! ///< Track parameters at clusters
	AliMFTTrackParam* fTrackParamAtVertex;   //! ///<  Track parameters at vertex
	AliMFTCATrack* fCATrack;                //! ///<  Track found by Track Finder (includes clusters)
	Int_t fTrackID;													///< Point to the corresponding MC track
	Double_t fP;                            ///<  Temporary to be removed....
	Double_t fTheta;                            ///<  Temporary to be removed....
	Double_t fPhi;                            ///<  Temporary to be removed....
	Double_t fPt;                            ///<  Temporary to be removed.... From sagitta estimation
	
	/// \cond CLASSIMP
	ClassDef(AliMFTTrack,1);
	/// \endcond
};

//=============================================================================================

#endif
