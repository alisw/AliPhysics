#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <AliOADBContainer.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TString.h>
#include <TRandom3.h>
#include "createTree.C"
#endif

class TCalRun : public TObject {
public:
  TCalRun() : fLedM(-1), fLedR(-1), fMonM(-1), fMonR(-1), fSMT(-1), fLength(0), fRunNo(-1), fBad(0) {;}
  virtual ~TCalRun() {;}
  Double32_t fLedM;       //[0,0,16] led mean
  Double32_t fLedR;       //[0,0,16] led rms
  Double32_t fMonM;       //[0,0,16] mon mean
  Double32_t fMonR;       //[0,0,16] mon rms
  Double32_t fSMT;        //[0,0,16] sm T
  UInt_t     fLength;     // run length in seconds
  Int_t      fRunNo;      // run number
  Short_t    fBad;        // bad cell
  ClassDef(TCalRun, 1); // CalRun class
};

class TCalCellInfo : public TObject {
public:
  TCalCellInfo() : fCellId(-1),  fRuns("TCalRun") {}
  virtual ~TCalCellInfo() {;}
  Int_t         fCellId;        // cell IDs
  TClonesArray  fRuns;
  ClassDef(TCalCellInfo, 1); // CalCellInfo class
};

void rotateTree(  const char *ifile     ="treefile.root",
                  const char *ofile     ="treefile.root",
                  Int_t referenceRun    = -1
               ) {

  // Load EMCAL geometry for reference run
  AliEMCALGeometry *g = AliEMCALGeometry::GetInstanceFromRunNumber(referenceRun);
  const Int_t kSM     = g->GetNumberOfSuperModules();
  const Int_t kNcells = g->GetNCells();

  // initialize info from tree created by $ALICE_PHYSICS/PWGPP/EMCAL/TeCMacros/createTree.C
  TCalInfo *info = 0;
  TFile *in = TFile::Open(ifile,"read");
  TTree *treeRuns = (TTree*)in->Get("tcal");
  treeRuns->SetBranchAddress("event",&info);
  treeRuns->Branch("event", &info, 32000, 99);
  Int_t Nev=treeRuns->GetEntries();

  TCalCellInfo *cellinfo = new TCalCellInfo;
  TFile *out = TFile::Open(ofile,"recreate");
  out->SetCompressionLevel(9);
  TTree* treeCells = new TTree("tcalcell", "Temp calibration cell tree");
  treeCells->SetDirectory(out);
  treeCells->Branch("cells", &cellinfo, 32000, 99);
  TClonesArray &cRunArr    = cellinfo->fRuns;


  for (Int_t k = 0; k < kNcells; k++){
    if (k%50 == 0)
      cout << "starting with cell " << k << "/" << kNcells << endl;

    cellinfo->fCellId       = k;

    cRunArr.Clear();
    cRunArr.ExpandCreate(Nev);

    for (Int_t i=0;i<Nev;++i) {
      treeRuns->GetEvent(i);

      TCalRun *run = (TCalRun*)cRunArr.At(i);
      run->fLength            = ((ULong_t)info->fLastTime-(ULong_t)info->fFirstTime);     // run length in seconds
      run->fRunNo             = info->fRunNo;

      TClonesArray &cells     = info->fCells;
      TCalCell *cell          = static_cast<TCalCell*>(cells.At(k));
      if (cell->fId != k){
        cout << "missmatch for cell " << k << "\t" << cell->fId << endl;
        continue;
      }
      run->fSMT               = cell->fSMT;
      run->fBad               = cell->fBad;
      run->fLedM              = cell->fLedM;
      run->fLedR              = cell->fLedR;
      run->fMonM              = cell->fMonM;
      run->fMonR              = cell->fMonR;
    }
    treeCells->Fill();
  }
  out->cd();
  treeCells->Write();
  out->Close();
}
