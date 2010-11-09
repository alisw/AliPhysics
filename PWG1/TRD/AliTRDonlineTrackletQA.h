#ifndef ALITRDONLINETRACKLETQA_H
#define ALITRDONLINETRACKLETQA_H

#include "AliAnalysisTask.h"

class TList;

class AliInputEventHandler;
class AliVEvent;
class AliAODEvent;
class AliMCEvent;

class AliTRDgeometry;
class AliTRDpadPlane;
class AliTRDtrackletMCM;
class AliTRDtrackletWord;

class AliTRDonlineTrackletQA : public AliAnalysisTask
{
 public:
  AliTRDonlineTrackletQA(const char *name);
  ~AliTRDonlineTrackletQA();

  void ConnectInputData(const Option_t *option);
  void CreateOutputObjects();
  void Exec(const Option_t *option);
  void LocalInit();
  void Terminate(const Option_t *option);

  void PlotMC(AliTRDtrackletMCM *trkl);
  void PlotESD(AliTRDtrackletMCM *trkl);

  void PlotESD(AliTRDtrackletWord *trkl);

  Int_t GetTrackletsForMC(Int_t label, Int_t idx[]);

 protected:
  AliESDEvent *fESD;                    //! current ESD event

  AliInputEventHandler *fInputHandler;  //! input handler
  AliVEvent            *fInputEvent;    //! input event
  AliAODEvent          *fOutputAOD;     //! output AOD
  AliMCEvent           *fMCEvent;       //! MC event

  TClonesArray         *fTrackletsRaw;  //! array of raw tracklets
  TClonesArray         *fTrackletsSim;  //! array of sim tracklets

  // ----- output objects -----
  TList                *fOutputList;	//! list of output objects
  TH1F                 *fHistYpos;      //! trkl y-position within chamber
  TH1F                 *fHistYres;      //! trkl y-residual to track reference
  TH2F                 *fHistYresDy;    //! trkl y-residual to track reference vs. deflection
  TH1F                 *fHistdY;	//! trkl dy 
  TH1F                 *fHistdYres;     //! trkl dy-residual to track reference
  TH2F                 *fHistYlocal[6]; //! trkl local y vs. MC (per layer)
  TH1F                 *fHistYresESD;   //! trkl y-residual to ESD track
  TH1F                 *fHistdYresESD;  //! trkl dy-residual to ESD track
  TH1F                 *fHistCanddY;    //! trkl cand at given dy
  TH1F                 *fHistFounddY;    //! trkl found at given dy
  TH1F                 *fHistTrklPerRef; //! no. of tracklets found per track reference
  TH2F                 *fHistdYdYref;    //! dy vs. dy from track reference
  TH1F                 *fHistYposRaw;    //! trkl y-position within chamber
  TH1F                 *fHistdYRaw;	//! trkl dy 
  TH2F                 *fHistYdYRaw;    //! y vs dy from raw
  TH1F                 *fHistZrow;      //! z-row
  TH1F                 *fHistZrowRaw;   //! z-row
  TH1F                 *fHistPid;       //! PID
  TH1F                 *fHistPidRaw;    //! PID
  TH1F                 *fHistYdiff;     //! difference in y between sim and raw
  TH1F                 *fHistdYdiff;    //! difference in dy between sim and raw
  TH2F                 *fHistdYdYraw;   //! dy sim vs. raw
  TTree                *fTreeTracklets; //! store tracklet information

  // ----- temporary storage -----
  Float_t fY;			// y-position
  Float_t fDY;			// deflection length
  Float_t fYdiff;		// difference in y-position
  Float_t fDYdiff;		// difference in deflection
  Int_t   fQ0;			// accumulated charge in first window
  Int_t   fQ1;			// accumulated charge in second window
  Int_t   fNHits;		// no. of hits

  // ----- configuration -----
  Float_t fMinPt;                       // minimal pt for tracks and track reference
					// which are used for comparison

  // ----- internal use -----
  AliTRDgeometry       *fGeo; //! TRD geometry

  Int_t fNevent;             // current event number

  TTree *fTrackletTree;  //! tracklets from simulation
  TTree *fTrackletTreeRaw; //! tracklets from raw data

 private:
  AliTRDonlineTrackletQA(const AliTRDonlineTrackletQA&); // not implemented
  AliTRDonlineTrackletQA& operator=(const AliTRDonlineTrackletQA&); // not implemented

  ClassDef(AliTRDonlineTrackletQA, 0);
};

#endif
