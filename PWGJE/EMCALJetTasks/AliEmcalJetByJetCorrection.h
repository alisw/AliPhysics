#ifndef ALIEMCALJETBYJETCORRECTION_H
#define ALIEMCALJETBYJETCORRECTION_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TClonesArray.h>

class TH3;
class AliEmcalJet;

class AliEmcalJetByJetCorrection : public TNamed
{
 public:

  AliEmcalJetByJetCorrection();
  AliEmcalJetByJetCorrection(const char* name);
  AliEmcalJetByJetCorrection(const AliEmcalJetByJetCorrection &jet);
  AliEmcalJetByJetCorrection& operator=(const AliEmcalJetByJetCorrection &jet);
  ~AliEmcalJetByJetCorrection() {;}

  void SetTemplate(TH3 *h)          { fh3JetPtDRTrackPt = h; }
  void SetJetPtBinWidth(Double_t w) { fBinWidthJetPt    = w; }
  void SetJetPtRange(Double_t min, Double_t max) {fJetPtMin = min; fJetPtMax = max;};
  void SetFixedTrackEfficiency(Double_t eff) { fEfficiencyFixed = eff; }
  void SetEfficiencyHist(TH1 *h)             { fhEfficiency     = h  ; }

  Int_t        GetJetPtBin(const Double_t jetpt) const;
  Double_t     GetEfficiency(const Double_t pt) const;
  Double_t     GetMeanPtConstituents(const AliEmcalJet *jet, TClonesArray *fTracks) const;

  void         Init();
  AliEmcalJet *Eval(const AliEmcalJet *jet, TClonesArray *fTracks);
  
 protected:
  TH3      *fh3JetPtDRTrackPt;                 // 3D template
  Double_t  fBinWidthJetPt;                    // jet pt bin width in which to do correction
  Double_t  fJetPtMin;                         // minimum jet pt
  Double_t  fJetPtMax;                         // maximum jet pt
  TObjArray fCollTemplates;                    // templates (2D histos with track pT vs r)
  Bool_t    fInitialized;                      // status of initialization
  Double_t  fEfficiencyFixed;                  // fixed efficiency for all pT and all types of tracks
  TH1      *fhEfficiency;                      // single particle efficiency
 private:

  ClassDef(AliEmcalJetByJetCorrection, 1) // jet-by-jet correction class
};
#endif
