#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TList.h"
#include "THashList.h"
#include "TMath.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliPHOSGeometry.h"
#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskPHOSTimeCalib.h"

// Author: Daiki Sekihata (Hiroshima University)
using namespace std;

ClassImp(AliAnalysisTaskPHOSTimeCalib)
//________________________________________________________________________
AliAnalysisTaskPHOSTimeCalib::AliAnalysisTaskPHOSTimeCalib(const char *name)
:	AliAnalysisTaskSE(name),
  fVEvent(0x0),
	fPHOSGeo(0),
	fOutputContainer(0),
	fRunNumber(0)
{
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());
}
//________________________________________________________________________
AliAnalysisTaskPHOSTimeCalib::~AliAnalysisTaskPHOSTimeCalib()
{

  //delete fPHOSGeo;


}
//________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

	if(fOutputContainer != NULL){
		delete fOutputContainer;
	}
 
	fOutputContainer = new THashList();
	fOutputContainer->SetOwner(kTRUE);


  TH1F *hEventSummary = new TH1F("hEventSummary","Event Summary",10,0.5,10.5);
  hEventSummary->GetXaxis()->SetBinLabel(1,"all");
  hEventSummary->GetXaxis()->SetBinLabel(2,"selected");
  hEventSummary->GetXaxis()->SetBinLabel(3,"kINT7");
  hEventSummary->GetXaxis()->SetBinLabel(4,"kPHI7 0PH0");
  hEventSummary->GetXaxis()->SetBinLabel(5,"kPHI7 1PHL");
  hEventSummary->GetXaxis()->SetBinLabel(6,"kPHI7 1PHM");
  hEventSummary->GetXaxis()->SetBinLabel(7,"kPHI7 1PHH");
  fOutputContainer->Add(hEventSummary);


	Char_t key[55];
	Char_t key1[55];
	Char_t key2[55];
  const Int_t Nmod=5;

  for(Int_t iddl=6;iddl<20;iddl++){
    for(Int_t igain=0;igain<2;igain++){
      fOutputContainer->Add(new TH2F(Form("hBC4vsRecTimeDDL%d_G%d",iddl,igain),Form("BC%%4 vs. Time DDL%d G%d",iddl,igain),400,-200,200,4,-0.5,3.5));
    }
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Post output data.
  fVEvent = InputEvent();
  if(!fVEvent) {
    AliError ("event not found. Nothing done!");
    return;
  }

  FillHistogram("hEventSummary",1);//all

  fRunNumber = fVEvent->GetRunNumber();

  if(!fPHOSGeo){

    if(fRunNumber < 209122)//Run1
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    else//Run2
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;

    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");

    for(Int_t mod=0; mod<6; mod++) {
      if(!matrixes->At(mod)) {
        if( fDebug )
          AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
        continue;
      }
      else {
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        if( fDebug >1 )
          AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
      }
    }//end of module loop
  }

  AliVCaloCells *cells = dynamic_cast<AliVCaloCells*>(fVEvent->GetPHOSCells());
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fVEvent);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fVEvent);

  UShort_t BC = fVEvent->GetBunchCrossNumber();
  //UShort_t BC = -999;
  //if(esd)      BC = esd->GetBunchCrossNumber();
  //else if(aod) BC = aod->GetBunchCrossNumber();


  CalibrateCellTime(cells,BC);

  PostData(1, fOutputContainer);
}      
//________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::Terminate(Option_t *) 
{
  // Called once at the end of the query

  fOutputContainer = dynamic_cast<THashList*>(GetOutputData(1));
  if (!fOutputContainer) {
    AliError ("Output list not available");
    return;
  }
 
}
//________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::CalibrateCellTime(AliVCaloCells *cells, UShort_t BC)
{
  Int_t relId[4]={};
  Double_t cellamp=0;
  Double_t celltime=0;
  Int_t module=0,cellx=0,cellz=0;
  Int_t cellAbsId=0;
  Short_t cellNumber=0;
  Int_t sru=-1, branch=-1, ddl=-1;
  Bool_t isHG=kFALSE;
  Double_t efrac=0;
  Int_t mclabel=-1;
  Double_t rawtime=0;
  Double_t corrtime=0;
  Int_t timeshift = -999;

  Int_t multCells = cells->GetNumberOfCells();
  //cout << "BC = " << BC << " , BC%4 = " << BC%4 << endl;

  Double_t value[4]={};

  for(Int_t iCell=0; iCell<multCells; iCell++){
    cells->GetCell(iCell, cellNumber, cellamp, celltime, mclabel, efrac);

    if(cellNumber<0) continue;//reject CPV clusters.

    isHG = cells->GetHighGain(iCell);
    cellAbsId = cells->GetCellNumber(iCell);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);

    module = relId[0];
    cellx  = relId[2];
    cellz  = relId[3];

    if(module < 1 || module > 4){
      AliError(Form("Wrong module number %d",module));
      return;
    }

    branch = WhichBranch(cellz);
    ddl = WhichDDL(module,cellx);

    if(cellamp>1.0){//above 1.0 GeV
      FillHistogram(Form("hBC4vsRecTimeDDL%d_G%d",ddl,isHG),celltime*1e+9,BC%4);
    }


  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSTimeCalib::WhichSRU(Int_t cellx)
{
  Int_t sru = -1;
  if(cellx<1 || 64<cellx){
    AliError ("cellx is wrong! sru=-1 will be returned");
    return -1;
  }
  else                    return (Int_t)((cellx-1)/16);
}
//______________________________________________________
Int_t AliAnalysisTaskPHOSTimeCalib::WhichBranch(Int_t cellz)
{
  Int_t branch = -1;
  if(cellz<1 || 56<cellz){
    AliError ("cellz is wrong! branch=-1 will be returned");
    return -1;
  }
  else                    return (Int_t)((cellz-1)/28);
}
//______________________________________________________
Int_t AliAnalysisTaskPHOSTimeCalib::WhichDDL(Int_t module, Int_t cellx)
{
  const Int_t Nmod=5;//totally, 5 PHOS modules are designed.
  Int_t ddl = -1;

  if(cellx<1 || 64<cellx) return -1;

  if(module<1 || 4<module){
    AliError ("module is wrong! ddl=-1 will be returned");
    return -1;
  }
  else{
    ddl = (Nmod-module) * 4 + (cellx-1)/16;//convert offline module numbering to online.
    return ddl;
  }
}
//______________________________________________________
void AliAnalysisTaskPHOSTimeCalib::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  TObject * obj = fOutputContainer->FindObject(key);
  TH1 * th1 = dynamic_cast<TH1*> (obj);
  TH2 * th2 = dynamic_cast<TH2*> (obj);
  if(th1){
    th1->Fill(x, y) ;
    return;
  }
  else if(th2){
    th2->Fill(x, y) ;
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskPHOSTimeCalib::FillHistogram(const char * key, Double_t x,Double_t y, Double_t z) const
{
  //FillHistogram
  TObject * obj = fOutputContainer->FindObject(key);

  TH2 * th2 = dynamic_cast<TH2*> (obj);
  if(th2) {
    th2->Fill(x, y, z) ;
    return;
  }
  TH3 * th3 = dynamic_cast<TH3*> (obj);
  if(th3) {
    th3->Fill(x, y, z) ;
    return;
  }

  AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}
//_______________________________________________________________________________
