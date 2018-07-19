////////////////////////////////////////////////////////////////////////////////
//345678901234567890123456789012345678901234567890123456789012345678901234567890
//       1         2         3         4         5         6         7         8
//
// Class AliTwoTrackRes
// J. Mercado <mercado@physi.uni-heidelberg.de> Last modified: 20.01.2011
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALITWOTRACKRES_H
#define ALITWOTRACKRES_H

#include <TNtuple.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include "TMath.h"
#include "TBits.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "AliAnalysisTask.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"

//______________________________________________________________________________
class AliTwoTrackRes : public AliAnalysisTask {

public:

  AliTwoTrackRes() : AliAnalysisTask("",""), fChain(0), fESDEvent(0), fOutContainer(0),
    fTrackCuts(0), fNTuple1(0), fNTuple2(0), fP1(), fP2(), fPb1(), fPb2(), fP(),
    fQ(), fTpcEnt1(), fTpcEnt2(), fTpcDist(), fOutFilename() {}
  AliTwoTrackRes(const char *name);
  AliTwoTrackRes(const AliTwoTrackRes& aTwoTrackRes);
  virtual ~AliTwoTrackRes();

  AliTwoTrackRes& operator=(const AliTwoTrackRes& aTwoTrackRes);

  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *opt = "");
  virtual void Terminate(Option_t *opt = "");

  void   SetTr1(double pt1, double eta1, double phi1, double m);
  void   SetTr2(double pt2, double eta2, double phi2, double m);
  void   SetTpcEnt1(double x1, double y1, double z1);
  void   SetTpcEnt2(double x2, double y2, double z2);
  double Qinv2()         {fQ = fP2 - fP1; return -1.*fQ.M2();}
  double Qinv()          {return TMath::Sqrt(TMath::Abs(Qinv2()));}
  double MinDist(AliExternalTrackParam *trk1, AliExternalTrackParam *trk2);
  int    GetNSha(TBits cl, TBits sh);
  double Corr(TBits cl1, TBits cl2, TBits sh1, TBits sh2);
  double Qfac(TBits cl1, TBits cl2, TBits sh1, TBits sh2);
  double RotTr2Phi();
  double Dist()    {fTpcDist = fTpcEnt2 - fTpcEnt1; return fTpcDist.Mag();}
  double DEta()    const {return TMath::Abs(fP2.Eta()-fP1.Eta());}
  double DTheta()  const {return fP2.Theta()-fP1.Theta();}
  double DPhi()    const {return TVector2::Phi_mpi_pi(fP2.Phi()-fP1.Phi());}
  void   NoSwap()        {fPb1 = fP1; fPb2 = fP2;}
  void   Swap()          {fPb1 = fP2; fPb2 = fP1;}
  void   FillNTuple1(double minsep, double sep, double corr, double qf,
		     int ns1, int ns2);
  void   FillNTuple2(double minsep, double sep, double corr, double qf,
		     int ns1, int ns2);
  void   SetOutfile(const char *outfil) {fOutFilename = outfil;}

private:

  TTree           *fChain;        // Chain of ESD trees
  AliESDEvent     *fESDEvent;     // Leaves
  TObjArray       *fOutContainer; // Output data container
  AliESDtrackCuts *fTrackCuts;    // Track cuts
  TNtuple         *fNTuple1;      // True pairs
  TNtuple         *fNTuple2;      // Mixed pairs
  TLorentzVector  fP1;            // Four-momentum of track 1 (lab)
  TLorentzVector  fP2;            // Four-momentum of track 2 (lab)
  TLorentzVector  fPb1;           // Buffer single-track 1 for swapping
  TLorentzVector  fPb2;           // Buffer single-track 2 for swapping
  TLorentzVector  fP;             // Total four-momentum (lab)
  TLorentzVector  fQ;             // Four-momentum difference (lab)
  TVector3        fTpcEnt1;       // Nominal TPC entrance point track 1
  TVector3        fTpcEnt2;       // Nominal TPC entrance point track 2
  TVector3        fTpcDist;       // Nominal TPC entrance separation
  TString         fOutFilename;   // Output filename

  ClassDef(AliTwoTrackRes, 0);
};
#endif

//______________________________________________________________________________
// EOF


