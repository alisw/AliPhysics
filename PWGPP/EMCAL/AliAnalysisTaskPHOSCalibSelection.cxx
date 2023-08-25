#include "AliAnalysisTaskPHOSCalibSelection.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskPHOSCalibSelection)

AliAnalysisTaskPHOSCalibSelection::AliAnalysisTaskPHOSCalibSelection(): 
AliAnalysisTaskSE(),
fCellEmin(),
fInputEvent(nullptr),
fPHOSCells(0x0),
fSaveCells(true),
fSaveFullTree(true),
fOutputContainer(0x0),
fCellIDVector(),
fCellTimeVector(),
fCellEnergyVector(),
fCellTree(NULL)
{
    

};

AliAnalysisTaskPHOSCalibSelection::AliAnalysisTaskPHOSCalibSelection(const char * name):
AliAnalysisTaskSE(name),
fCellEmin(),
fInputEvent(nullptr),
fPHOSCells(0x0),
fSaveCells(true),
fSaveFullTree(true),
fOutputContainer(0x0),
fCellIDVector(),
fCellTimeVector(),
fCellEnergyVector(),
fCellTree(NULL)
{
   
   DefineOutput(1, TTree::Class());
};

AliAnalysisTaskPHOSCalibSelection::~AliAnalysisTaskPHOSCalibSelection()
{
  // default deconstructor
}

void AliAnalysisTaskPHOSCalibSelection::UserExec(Option_t* /* option */) {

  fInputEvent = InputEvent();

  if( fSaveCells ) ProcessCells();

  fCellTree->Fill();

  ResetBuffer();

  PostData(1, fCellTree);


}

void AliAnalysisTaskPHOSCalibSelection::UserCreateOutputObjects(){


    if(fSaveCells){
      
      fCellTree = new TTree(Form("PHOSCells"),"PHOSCells");

      fCellTree->Branch("Cell_ID", "std::vector<int>",      &fCellIDVector);
      fCellTree->Branch("Cell_E",  "std::vector<float>",      &fCellEnergyVector);
      fCellTree->Branch("Cell_time", "std::vector<float>",       &fCellTimeVector);

    }

  OpenFile(1);
  PostData(1, fCellTree);

}


void AliAnalysisTaskPHOSCalibSelection::ProcessCells(){

    fPHOSCells = fInputEvent->GetPHOSCells();
    int nCells = fPHOSCells->GetNumberOfCells();

    int nCellpass = 0;
    
    for(int i = 0; i < nCells; i++){
        if(fPHOSCells->GetCellAmplitude(i) > fCellEmin){
            fCellIDVector.push_back(fPHOSCells->GetCellNumber(i));
            fCellTimeVector.push_back(fPHOSCells->GetCellTime(i)*1e9);
            fCellEnergyVector.push_back(fPHOSCells->GetCellAmplitude(i));

        }
    }


    return;



};


void AliAnalysisTaskPHOSCalibSelection::ResetBuffer(){

  fCellIDVector.clear();
  fCellTimeVector.clear();
  fCellEnergyVector.clear();

}




