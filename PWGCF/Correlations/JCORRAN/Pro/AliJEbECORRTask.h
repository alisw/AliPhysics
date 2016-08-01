#ifndef AliJEBECORR_cxx
#define AliJEBECORR_cxx

#include <TVector.h>
#include <TRandom.h>
#include <TString.h>
#include <TPRegexp.h>
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliJEfficiency.h"
#include "AliJFFlucAnalysis.h"
#include "AliMCParticle.h"

//#include <fstream>

class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliESDVertex;
class TClonesArray;
class AliJCard;
class AliJBaseTrack;
class AliJHistos;
class AliJEbeHistos;
class AliJCorrelations;
class AliJEbePercentile;
class AliJEfficiency;
class AliJEventPool;
class AliJRunTable;
class AliJFFlucAnalysis;

class AliJEbECORRTask : public AliAnalysisTaskSE {

	public:
		AliJEbECORRTask();
		AliJEbECORRTask(const char *name);
		AliJEbECORRTask(const AliJEbECORRTask& a); // not implemented
		AliJEbECORRTask& operator=(const AliJEbECORRTask& ap); // not implemented

		virtual ~AliJEbECORRTask(); 

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);
		virtual Bool_t UserNotify();

		bool IsGoodEvent(AliVEvent *event);
		double RunEbEFlowAnalysis(AliVEvent *event, TClonesArray* inputlist);
		void PlayCorrelation(TClonesArray *triggList, TClonesArray *assocList);
		void SetDebugMode( int debug) { fDebugMode = debug; };
		void ScaleNotEquidistantHisto(TH1D *hid, const double sc);
		AliJCard *GetCard() { return fCard; }
		void SetCard( AliJCard *c ) { fCard = c; }
		void SetEbePercentileInputFileName(TString name) { ebePercentileInputFileName = name; };
		void SetIsMC( Bool_t ismc){ IsMC = ismc; cout << "Settint IsMC = " << ismc << endl; };
		void SetEnableCORR( Bool_t runCORR ){ fenableCORR = runCORR; cout << "Settint enableCORR = " << fenableCORR << endl; };
		double GetCentralityFromImpactPar(double ip);
		double GetImpactParFromCentrality(double cent);
		void ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList);
		Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
		void SetKineOnly (Bool_t iskineonly) {IsKinematicOnly = iskineonly;};
		/*
		   TString MemoryStatus(){
		   int pid = getpid();
		   TString key;
		   TString v1;
		   ifstream ins(Form("/proc/%d/status",pid ));
		   double VmSwap= 0.;
		   double VmSize= 0.;
		   double VmPeak= 0.;
		   double nom = TMath::Power( 2, 20 );
		   TPMERegexp split("\\s+");
		   TString line;
		   while( ins.good() ){
		   line.ReadLine(ins);
		   split.Split(line);
		   if( split.NMatches() < 2 ) continue;
		   key = split[0];
		   v1 = split[1];
		   if( key == "VmSwap:" ) VmSwap = v1.Atof()/nom;
		   if( key == "VmSize:" ) VmSize = v1.Atof()/nom;
		   if( key.BeginsWith("VmPeak") ) VmPeak = v1.Atof()/nom;
		   }
		   TString res = Form("VmPeak:%10.2f VmSize:%10.2f VmSwap:%10.2f", VmPeak, VmSize, VmSwap);
		   ins.close();
		   return res;
		   };
		 */

	private:
		TDirectory           *fOutput;     // Output

		AliAnalysisUtils *fAnaUtils;
		AliJFFlucAnalysis *fFFlucAna;

		AliJCard * fCard;

		AliJHistos * fHistos;				//!
		AliJEbeHistos * fEbeHistos;			//!
		AliJEfficiency *fEfficiency;
		int fHadronSelectionCut;

		TClonesArray * fInputList;
		TClonesArray * fInputListSpectra;
		TClonesArray * ftriggList;
		TClonesArray * fassocList;

		AliJCorrelations *fcorrelations;	//!
		AliJEventPool *fassocPool;			//!
		AliJEbePercentile *fEbePercentile;	//!

		Bool_t fFirstEvent; //
		int cBin;
		int ebeBin;
		int zBin;
		double zVert;
		Int_t fevt; // event number
		int fDebugMode;
		Int_t trkfilterBit;
		TVector *fEbECentBinBorders;
		TString ebePercentileInputFileName;
		AliJRunTable *fRunTable; // 

		TRandom *fRandom;
		Bool_t IsMC;
		Bool_t IsKinematicOnly;
		Bool_t fenableCORR;

		ClassDef(AliJEbECORRTask, 1); // example of analysis
};

#endif
