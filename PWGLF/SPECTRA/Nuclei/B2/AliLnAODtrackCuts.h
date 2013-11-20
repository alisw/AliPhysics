#ifndef ALILNAODTRACKCUTS_H
#define ALILNAODTRACKCUTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// AOD track cuts for B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

class AliAODEvent;
class AliAODMCParticle;
class AliAODTrack;
class TString;

class AliLnAODtrackCuts: public TObject
{
  public:
	
	AliLnAODtrackCuts();
	virtual ~AliLnAODtrackCuts();
	
	Bool_t IsWithinGeoAcceptance(const AliAODMCParticle* prt) const;
	Bool_t IsWithinGeoAcceptance(const AliAODTrack* trk) const;
	Bool_t IsWithinGeoAcceptance(Double_t p[3]) const;
	
	Bool_t IsKinkDaughter(const AliAODTrack* trk) const;
	
	Bool_t AcceptItsTpcNSigma(const AliAODTrack* trk, Double_t b[2], Double_t bCov[3]) const;
	Bool_t AcceptItsTpcDCA(const AliAODTrack* trk, Double_t b[2]) const;
	Bool_t AcceptItsTpcStdCut(const AliAODTrack* trk, Double_t b[2]) const;
	
	Bool_t AcceptTOF(const AliAODTrack* trk) const;
	
	Bool_t AcceptTrack(const AliAODTrack* trk, Double_t b[2], Double_t bCov[3]) const;
	
	Bool_t TOFmatch() const { return fTOFmatch; }
	
	Double_t GetIntegratedLength(const AliAODTrack* trk, Int_t pid=0) const;
	Double_t GetNSigmaToVertex(Double_t b[2], Double_t bCov[3]) const;
	
	Double_t GetNTPCXRowsOverFindable(const AliAODTrack* trk) const;
	
	void SetSelectionCriteria(const TString& trksel);
	
	void SetMaxDCAxy(Double_t max) { fMaxDCAxy=max; }
	void SetMaxDCAz(Double_t max) { fMaxDCAz=max; }
	void SetMaxNSigmaToVtx(Double_t max) { fMaxNSigma=max; }
	void SetMaxEta(Double_t max) { fMaxEta=max; }
	void SetTPCXRowsCut(Bool_t xrows=1) { fTPCXRows=xrows; }
	void SetMinTPCnClsOrXRows(Int_t min) { fMinTPCnClsOrXRows=min; }
	
  private:
	
	TString fTrackSel; // selection criteria
	Double_t fMaxDCAxy; // maximum dcaxy value
	Double_t fMaxDCAz; // maximum dcaz value
	Double_t fMaxNSigma; // maximum number of sigmas to primary vertex
	Double_t fMaxEta; // maximum pseudorapidity value
	Bool_t fTPCXRows; // enable cut on the number of crossed rows
	Int_t fMinTPCnClsOrXRows; // minimum number of TPC clusters or crossed rows
	Bool_t fTOFmatch; // if a TOF match signal is required
	
	ClassDef(AliLnAODtrackCuts,1)
};

#endif // ALILNAODTRACKCUTS_H
