#include "AliHFCutVarFDsubAnalysisManagerDplus.h"

#include <iostream>
#include <vector>

#include "TH1F.h"
#include "TList.h"
#include "TFile.h"
#include "TDatabasePDG.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubAnalysisManagerDplus);
/// \endcond

using std::cerr;
using std::endl;
using std::vector;

AliHFCutVarFDsubAnalysisManagerDplus::AliHFCutVarFDsubAnalysisManagerDplus()
  : AliHFCutVarFDsubAnalysisManager()
  , fNevents(0.)
{
  /// Default constructor
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerDplus::~AliHFCutVarFDsubAnalysisManagerDplus() {
  /// Destructor
}


//_________________________________________________________________________________________________
Int_t AliHFCutVarFDsubAnalysisManagerDplus::GetTHnSparses(const TString strFileData,
                                                          const TString strFileMC,
                                                          const TString strDirData,
                                                          const TString strDirMC,
                                                          const TString strListData,
                                                          const TString strListMC,
                                                          Bool_t MConly) {
  /// Retrieve the input data
  //------------------------------------------------------------------------------------------------
  //MC THnSparses
  TFile *infileMC = TFile::Open(strFileMC.Data(),"READ");//open MC file
  if (!infileMC) {cerr << "File " << strFileMC << " not found!" << endl; return 1;}
  TDirectoryFile *DplusDirMC = (TDirectoryFile*)infileMC->Get(strDirMC.Data());//open MC directory
  if (!DplusDirMC) {cerr << "Directory " << strDirMC << " not found!" << endl; return 2;}
  TList *ListMC = (TList*)DplusDirMC->Get(strListMC.Data());//get MC list
  if (!ListMC) {cerr << "List " << strListMC << " not found!" << endl; return 3;}

  //get THnSparses
  fMCafterCuts[kPrompt] = (THnSparseF*)ListMC->FindObject("hMassPtImpParPrompt");
  if (!fMCafterCuts[kPrompt]) { cerr << "THnSparseF hMassPtImpParPrompt not found!" << endl; return 4; }
  fMCafterCuts[kFD] = (THnSparseF*)ListMC->FindObject("hMassPtImpParBfeed");
  if (!fMCafterCuts[kFD]) { cerr << "THnSparseF hMassPtImpParBfeed not found!" << endl; return 5; }
  fMCgenLevel[kPrompt] = (THnSparseF*)ListMC->FindObject("hMCAccPrompt");
  if (!fMCgenLevel[kPrompt]) { cerr << "THnSparseF hMCAccPrompt not found!" << endl; return 6; }
  fMCgenLevel[kFD] = (THnSparseF*)ListMC->FindObject("hMCAccBFeed");
  if (!fMCgenLevel[kFD]) { cerr << "THnSparseF hMCAccBFeed not found!" << endl; return 7; }
  //-------------------------------------------------------------------------------------------------
  //Real Data
  TH1F *hNev = 0x0;
  if(MConly) {
    //get THnSparse from MC file
    fData = (THnSparseF*)ListMC->FindObject("hMassPtImpParAll");
    if (!fData) { cerr << "THnSparseF hMassPtImpParAll not found!" << endl; return 8; }
    hNev = (TH1F*)ListMC->FindObject("fHistNEvents");
    if(!hNev) { cerr << "THnSparseF fHistNEvents not found!" << endl; return 9; }
  }
  else {
    TFile *infileData = TFile::Open(strFileData.Data(),"READ");//open Data file
    if (!infileData) {cerr << "File " << strFileData << " not found!" << endl; return 10;}
    TDirectoryFile *DplusDirData = (TDirectoryFile*)infileData->Get(strDirData.Data());//open Data directory
    if (!DplusDirData) {cerr << "Directory " << strDirData << " not found!" << endl; return 11;}
    TList *ListData = (TList*)DplusDirData->Get(strListData.Data());//get Data list
    if (!ListData) {cerr << "List " << strListData << " not found!" << endl; return 12;}

    //get THnSparse from Data file
    fData = (THnSparseF*)ListData->FindObject("hMassPtImpParAll");
    if (!fData) { cerr << "THnSparseF hMassPtImpParAll not found!" << endl; return 8; }
    hNev = (TH1F*)ListData->FindObject("fHistNEvents");
    if(!hNev) { cerr << "THnSparseF fHistNEvents not found!" << endl; return 9; }
  }

  fNevents = hNev->GetBinContent(1);

  return 0;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerDplus::GetCuts() {
  /// Preparte the cuts

  Double_t massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();

  if (fCuts) {
    delete fCuts;
    fCuts = 0x0;
  }
  fCuts= new TList();
  fCuts->SetOwner(0);

  //pT bin [1,2] GeV/c
  TList *ptBin = new TList();

  //cut set 1 - Prompt enhanced
  AliHFCutVarFDsubCutSet *cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.007,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,1.,2.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.998,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.1));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,7,10));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  AliHFCutVarFDsubCutSet *cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.007,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,1.,2.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,7,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  AliHFCutVarFDsubCutSet *cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.007,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,1.,2.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.998));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.04,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,9,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //______________________________________________________________________________________________________________//

  //pT bin [2,3] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.008,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,2.,3.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.998,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.14));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,6,10));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.007,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,2.,3.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,7,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.009,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,2.,3.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.998));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.05,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [3,4] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.006,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,3.,4.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.14));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,10));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.008,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,3.,4.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.009,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,3.,4.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.998));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.05,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [4,5] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.008,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,4.,5.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.14));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,12));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.009,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,4.,5.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.0011,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,4.,5.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.999));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.06,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [5,6] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.009,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,5.,6.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.16));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,12));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.010,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,5.,6.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.014,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,5.,6.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.999));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.06,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [6,7] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.011,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,6.,7.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.16));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,12));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.010,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,6.,7.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.012,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,6.,7.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.999));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.06,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [7,8] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.013,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,7.,8.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.16));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,12));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.011,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,7.,8.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.011,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,7.,8.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,0.999));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.06,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,8,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [8,12] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.016,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,8.,12.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.04,0.16));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,12));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.013,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,8.,12.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.011,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,8.,12.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,1.));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.15,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,22,31));//normdlxy
  ptBin->Add((TObject*)cutset3);

  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [12,16] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.013,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,12.,16.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.25));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,14));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.014,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,12.,16.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.016,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,12.,16.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,1.));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.14,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,24,31));//normdlxy
  ptBin->Add((TObject*)cutset3);
  fCuts->Add((TObject*)ptBin);

  //_________________________________________________________________________________________________________//

  //pT bin [16,24] GeV/c
  ptBin = new TList();

  //cut set 1 - Prompt enhanced
  cutset1 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.011,4,0,0,0);
  cutset1->AddCut(new AliHFCutVarFDsubCut(0,16.,24.));
  cutset1->AddCut(new AliHFCutVarFDsubCut(1,0.999,1.));//cxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(2,0.1,0.3));//dlxy
  cutset1->AddCut(new AliHFCutVarFDsubCut(3,5,15));//normdlxy
  ptBin->Add((TObject*)cutset1);

  //cut set 2 - mixed
  cutset2 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.009,4,0,0,0);
  cutset2->AddCut(new AliHFCutVarFDsubCut(0,16.,24.));
  cutset2->AddCut(new AliHFCutVarFDsubCut(1,0.997,1.));//cxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(2,0.06,0.4));//dlxy
  cutset2->AddCut(new AliHFCutVarFDsubCut(3,5,28));//normdlxy
  ptBin->Add((TObject*)cutset2);

  //cut set 3 - FD enhanced
  cutset3 = new AliHFCutVarFDsubCutSet(1.72,2.05,massD,-0.014,4,0,0,0);
  cutset3->AddCut(new AliHFCutVarFDsubCut(0,16.,24.));
  cutset3->AddCut(new AliHFCutVarFDsubCut(1,0.95,1.));//cxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(2,0.12,1.1));//dlxy
  cutset3->AddCut(new AliHFCutVarFDsubCut(3,18,31));//normdlxy
  ptBin->Add((TObject*)cutset3);
  fCuts->Add((TObject*)ptBin);

}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerDplus::GetAxes(Int_t version) {
  if (fAxes) {
    delete fAxes;
    fAxes = 0x0;
  }
  fAxes = new TList();
  fAxes->SetOwner();
  switch (version) {
  case 0:
  default:
    // standard version
    // -1 means axis is not existent (or not used)
    //
    // On the generator level, one doesn't want to cut on the variation variables, but one would
    // like to select e.g. pt bins. In order to be to use the same cut for this, variables on
    // which one doesn't what to cut at generator level have the corresponding axes set to -1
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(1,0,1,"pt")); //pt axis
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(3,(UInt_t)-1,3,"cxy")); //cosP axis
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(4,(UInt_t)-1,4,"dlxy")); //decay length xy axis
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(5,(UInt_t)-1,5,"normdlxy")); //normalized decay length xy axis
    break;
  }
}


//_________________________________________________________________________________________________
TH1F* AliHFCutVarFDsubAnalysisManagerDplus::CalculateCrossSection(TString AccFilePath,
                                                                  TString GenLimAccHistoName,
                                                                  TString GenAccHistoName,
                                                                  TString system,
                                                                  TString PromptOrFD) {

  if(!fCorrYieldPrompt || !fCorrYieldFD) {
    cerr << "Minimisation needed before calculate cross section" << endl;
    return 0x0;
  }
  else if(AccFilePath=="") {
    cerr << "The Acceptance Factor needed to calculate the cross section" << endl;
    return 0x0;
  }
  else {
    TFile *AccFile = TFile::Open(AccFilePath.Data(),"READ");
    if (!AccFile) {
      cerr << "File " << AccFilePath << " not found!" << endl;
      return 0x0;
    }
    else {
      TH1F* hPtGenLimAcc = (TH1F*)AccFile->Get(GenLimAccHistoName.Data());
      TH1F* hPtGenAcc = (TH1F*)AccFile->Get(GenAccHistoName.Data());
      hPtGenAcc->SetDirectory(0);
      hPtGenLimAcc->SetDirectory(0);
      AccFile->Close();

      TAxis* PtAxis = (TAxis*)fCorrYieldPrompt->GetXaxis();
      TArrayD* PtBinsArray = (TArrayD*)PtAxis->GetXbins();
      Double_t* PtLims = (Double_t*)PtBinsArray->GetArray();
      Int_t nPtBins = PtBinsArray->GetSize()-1;

      TH1F* hPtGenLimAccReb = (TH1F*)hPtGenLimAcc->Rebin(nPtBins, "hPtGenLimAccReb", PtLims);
      TH1F* hPtGenAccReb = (TH1F*)hPtGenAcc->Rebin(nPtBins, "hPtGenAccReb", PtLims);
      hPtGenLimAccReb->Sumw2();
      hPtGenAccReb->Sumw2();

      TH1F* hPtAcc = new TH1F("hPtAcc","",nPtBins,PtLims);
      hPtAcc->Divide(hPtGenAccReb,hPtGenLimAccReb,1.,1.,"B");

      TH1F* hCrossSec = new TH1F(Form("hCrossSec%s",PromptOrFD.Data()),"",nPtBins,PtLims);

      Double_t BR=0.0913;
      Double_t ErrBR=0.0019;
      Double_t Sigma=0;
      Double_t ErrSigma=0;

      if(system=="pPb") {
        Sigma = 2.09;//barn
        ErrSigma = Sigma*3.5/100;
      }
      else if(system=="pp"){
        Sigma = 0.062;//barn
        ErrSigma = Sigma*3.5/100;
      }
      else {
        cerr << "only pPb and pp are implemented" << endl;
        return 0x0;
      }

      for(Int_t iBin=0; iBin<nPtBins; iBin++) {

        Double_t CrossSec;
        Double_t CrossSecError;

        if(PromptOrFD=="Prompt") {
          //cross section in mubarn
          CrossSec = fCorrYieldPrompt->GetBinContent(iBin+1)*Sigma*1000000./(2*(PtLims[iBin+1]-PtLims[iBin])*fNevents*BR);
          CrossSecError = TMath::Sqrt((fCorrYieldPrompt->GetBinError(iBin+1)/fCorrYieldPrompt->GetBinContent(iBin+1))*(fCorrYieldPrompt->GetBinError(iBin+1)/fCorrYieldPrompt->GetBinContent(iBin+1))+(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))*(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))+(ErrBR/BR)*(ErrBR/BR)+(ErrSigma/Sigma)*(ErrSigma/Sigma))*CrossSec;
        }
        else {
          CrossSec = fCorrYieldFD->GetBinContent(iBin+1)*Sigma*1000000./(2*(PtLims[iBin+1]-PtLims[iBin])*fNevents*BR);
          CrossSecError = TMath::Sqrt((fCorrYieldFD->GetBinError(iBin+1)/fCorrYieldFD->GetBinContent(iBin+1))*(fCorrYieldFD->GetBinError(iBin+1)/fCorrYieldFD->GetBinContent(iBin+1))+(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))*(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))+(ErrBR/BR)*(ErrBR/BR)+(ErrSigma/Sigma)*(ErrSigma/Sigma))*CrossSec;
        }

        hCrossSec->SetBinContent(iBin+1,CrossSec);
        hCrossSec->SetBinError(iBin+1,CrossSecError);
      }

      hCrossSec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      hCrossSec->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
      hCrossSec->GetYaxis()->SetTitleOffset(1.1);

      delete hPtAcc;
      hPtAcc = 0x0;

      return hCrossSec;
    }
  }
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerDplus::AliHFCutVarFDsubAnalysisManagerDplus(const AliHFCutVarFDsubAnalysisManagerDplus& am)
  : AliHFCutVarFDsubAnalysisManager(am)
  , fNevents(am.fNevents)
{
  /// Copy constructor

  // TODO: do proper deep copies
  Printf("Do not use this copy constructor!");
}
