#include "TRDModule.h"
#include "TRDModuleImp.h"

#include "TGButton.h"
#include "TGColorSelect.h"
#include "TGNumberEntry.h"
#include "TG3DLine.h"
//#include "TColor.h"
//#include "TMath.h"

#include "Reve/RGTopFrame.h"
#include "Reve/Track.h"


using namespace Reve;
using namespace Alieve;
//using namespace std;

ClassImp(TRDModule)
ClassImp(TRDModuleEditor)

//________________________________________________________
TRDModule::TRDModule(const char *typ, const Int_t det) :
  TNamed(Form("%s%03d", typ, det), ""), fDet(det)
{
	fLoadHits = kFALSE;
	fRnrHits = kTRUE;
	
	fLoadDigits = kFALSE;
	fRnrDigits = kTRUE;
	fDigitsLog = kFALSE;
	fDigitsBox = kFALSE;

	fLoadRecPoints = kFALSE;
	fRnrRecPoints = kTRUE;

	fLoadTracklets = kFALSE;
	fRnrTracklets = kTRUE;

	kDigitsNeedRecompute = kTRUE;

	fDigitsThreshold = 15;
}



///////////////////////////////////////////////////////////
/////////////   TRDModuleEditor       /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDModuleEditor::TRDModuleEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options | kVerticalFrame, back)
{
	MakeTitle("TRD Detector");
	
	// "Hits" group frame
	TGGroupFrame *fGroupFrame615 = new TGGroupFrame(this,"Hits");
	fDisplayHits = new TGCheckButton(fGroupFrame615,"Display");
	fDisplayHits->SetToolTipText("Display hits attached to chamber");
	fDisplayHits->SetState(kButtonDown);
	fGroupFrame615->AddFrame(fDisplayHits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fDisplayHits->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "UpdateChamber()");

	// color select widget
	ULong_t ColPar1;
	gClient->GetColorByName("#ff00ff", ColPar1);
	fHitsColor = new TGColorSelect(fGroupFrame615, ColPar1, -1);

	fGroupFrame615->AddFrame(fHitsColor, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));
	fHitsColor->Connect("ColorSelected(Pixel_t)",
					"Alieve::TRDModuleEditor", this, "UpdateHits(Pixel_t)");
	
	fGroupFrame615->SetLayoutManager(new TGVerticalLayout(fGroupFrame615));
	AddFrame(fGroupFrame615, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	// "Digits" group frame
	TGGroupFrame *fGroupFrame636 = new TGGroupFrame(this,"Digits");
	TGVerticalFrame *fVerticalFrame538 = new TGVerticalFrame(fGroupFrame636,68,63,kVerticalFrame);
	fDisplayDigits = new TGCheckButton(fVerticalFrame538,"Display");
	fDisplayDigits->SetState(kButtonDown);
	fDisplayDigits->SetToolTipText("Toggle display digits attached to chamber");
	fVerticalFrame538->AddFrame(fDisplayDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fDisplayDigits->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "UpdateChamber()");
	fToggleLog = new TGCheckButton(fVerticalFrame538,"Log");
	fToggleLog->SetState(kButtonUp);
	fToggleLog->SetToolTipText("Toggle logarithmic representation of digits");
	fVerticalFrame538->AddFrame(fToggleLog, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fToggleLog->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "ModifyDigitsView()");

	fToggleBox = new TGCheckButton(fVerticalFrame538,"Box");
	fToggleBox->SetState(kButtonUp);
	fToggleBox->SetToolTipText("Toggle display digits as boxes");
	fVerticalFrame538->AddFrame(fToggleBox, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fToggleBox->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "ModifyDigitsView()");
	fGroupFrame636->AddFrame(fVerticalFrame538, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   
   TGVertical3DLine *fVertical3DLine773 = new TGVertical3DLine(fGroupFrame636,8,63);
   fGroupFrame636->AddFrame(fVertical3DLine773, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsExpandY,2,2,2,2));
	 
	 // vertical frame
	TGVerticalFrame *fVerticalFrame543 = new TGVerticalFrame(fGroupFrame636,81,47,kVerticalFrame);
	fThreshold = new TGCheckButton(fVerticalFrame543, "Threshold");
	fThreshold->SetToolTipText("Toggle threshold updating display");
	fVerticalFrame543->AddFrame(fThreshold, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fThresValue = new TGNumberEntry(fVerticalFrame543, (Int_t)10, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 0, 1024);
	fThresValue->Connect("ValueSet(Long_t)",
		  "Alieve::TRDModuleEditor", this, "SetThreshold(Long_t)");
	fVerticalFrame543->AddFrame(fThresValue, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	fGroupFrame636->AddFrame(fVerticalFrame543, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	fGroupFrame636->SetLayoutManager(new TGHorizontalLayout(fGroupFrame636));
	AddFrame(fGroupFrame636, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	// "Clusters" group frame
	TGGroupFrame *fGroupFrame649 = new TGGroupFrame(this,"Clusters");
	fDisplayClusters = new TGCheckButton(fGroupFrame649,"Display");
	fDisplayClusters->SetToolTipText("Display clusters attached to chamber");
	fGroupFrame649->AddFrame(fDisplayClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fDisplayClusters->SetState(kButtonDown);
	fDisplayClusters->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "UpdateChamber()");

	// color select widget
	ULong_t ColPar2;
	gClient->GetColorByName("#ff00ff", ColPar2);
	fClustersColor = new TGColorSelect(fGroupFrame649, ColPar2, -1);

	fGroupFrame649->AddFrame(fClustersColor, new TGLayoutHints(kLHintsLeft | kLHintsCenterX | kLHintsTop | kLHintsCenterY,2,2,2,2));
	fClustersColor->Connect("ColorSelected(Pixel_t)",
					"Alieve::TRDModuleEditor", this, "UpdateClusters(Pixel_t)");

	fGroupFrame649->SetLayoutManager(new TGVerticalLayout(fGroupFrame649));
	AddFrame(fGroupFrame649, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));

	// "Tracks" group frame
	TGGroupFrame *fGroupFrame676 = new TGGroupFrame(this,"Tracklets");
	fDisplayTracks = new TGCheckButton(fGroupFrame676,"Display");
	fDisplayTracks->SetToolTipText("Display MCM tracklets attached to chamber");
	fGroupFrame676->AddFrame(fDisplayTracks, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
	fDisplayTracks->SetState(kButtonDown);
	fDisplayTracks->Connect("Clicked()",
					"Alieve::TRDModuleEditor", this, "UpdateChamber()");

	fGroupFrame676->SetLayoutManager(new TGVerticalLayout(fGroupFrame676));
	AddFrame(fGroupFrame676, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY | kLHintsExpandX,2,2,2,2));
}

//________________________________________________________
TRDModuleEditor::~TRDModuleEditor()
{}

//_________________________________________________________
void TRDModuleEditor::SetModel(TObject* obj)
{
// Takes care to update UI state according to model state

	fM = dynamic_cast<TRDModule*>(obj);
	if(TRDNode *node = dynamic_cast<TRDNode*>(fM)) node->UpdateNode();
	
	fDisplayHits->SetEnabled(fM->fLoadHits);
	if(fM->fLoadHits) fDisplayHits->SetState( fM->fRnrHits ? kButtonDown : kButtonUp);
		
	fDisplayDigits->SetEnabled(fM->fLoadDigits);
	fToggleLog->SetEnabled(fM->fLoadDigits);
	fToggleBox->SetEnabled(fM->fLoadDigits);
	fThreshold->SetEnabled(fM->fLoadDigits);
	fThresValue->GetNumberEntry()->SetEnabled(fM->fLoadDigits);
	fThresValue->GetButtonDown()->SetEnabled(fM->fLoadDigits);
	fThresValue->GetButtonUp()->SetEnabled(fM->fLoadDigits);
	fThresValue->SetIntNumber(fM->fDigitsThreshold);
	if(fM->fLoadDigits){
		fDisplayDigits->SetState(fM->fRnrDigits ? kButtonDown : kButtonUp);
		fToggleLog->SetState(fM->fDigitsLog ? kButtonDown : kButtonUp);
		fToggleBox->SetState(fM->fDigitsBox ? kButtonDown : kButtonUp);
		fThreshold->SetState(fM->kDigitsNeedRecompute ? kButtonDown : kButtonUp);
	}
	
	if(fM->fLoadRecPoints) fDisplayClusters->SetState(fM->fRnrRecPoints ? kButtonDown : kButtonUp);
	else fDisplayClusters->SetEnabled(kFALSE);
	
	if(fM->fLoadTracklets) fDisplayTracks->SetState(fM->fRnrTracklets ? kButtonDown : kButtonUp);
	else fDisplayTracks->SetEnabled(kFALSE);
}

//_________________________________________________________
void TRDModuleEditor::ModifyDigitsView()
{
	fM->kDigitsNeedRecompute = kTRUE;
	UpdateChamber();
}

//________________________________________________________
void TRDModuleEditor::SetThreshold(Long_t tres)
{
	if(!fM->fLoadDigits) return;
	
	fM->fDigitsThreshold = (tres == 0) ? (int)fThresValue->GetNumber():tres;
	fM->kDigitsNeedRecompute = fThreshold->IsDown();
	UpdateChamber();
}


//_________________________________________________________
void TRDModuleEditor::UpdateChamber()
{
//	Info("UpdateChamber()", "");

	fM->fRnrHits = kFALSE;
	if(fDisplayHits->IsDown()) fM->fRnrHits = kTRUE;
	
	fM->fRnrDigits = kFALSE;
	if(fDisplayDigits->IsDown()){
		fM->fRnrDigits = kTRUE;
		if(fToggleLog->IsDown()) fM->fDigitsLog = kTRUE;
		else  fM->fDigitsLog = kFALSE;
		if(fToggleBox->IsDown()) fM->fDigitsBox = kTRUE;
		else  fM->fDigitsBox = kFALSE;
	}
	
	fM->fRnrRecPoints = kFALSE;
	if(fDisplayClusters->IsDown()) fM->fRnrRecPoints = kTRUE;
	
	fM->fRnrTracklets = kFALSE;
	if(fDisplayTracks->IsDown()) fM->fRnrTracklets = kTRUE;

	if(dynamic_cast<TRDNode*>(fM)) (dynamic_cast<TRDNode*>(fM))->UpdateLeaves();

	gReve->Redraw3D();
}

//_________________________________________________________
void TRDModuleEditor::UpdateClusters(Pixel_t /*col*/)
{
//	Info("UpdateClusters()","");
/*	
	if(fM->fRecPoints){
		fM->fRecPoints->SetMarkerColor(col ? fHitsColor->GetColor() : col);
		gReve->Redraw3D();
	}
*/
}

//_________________________________________________________
void TRDModuleEditor::UpdateHits(Pixel_t /*col*/)
{
//	Info("UpdateHits()","");
/*	
	if(fM->fHits){
		fM->fHits->SetMarkerColor(col ? fClustersColor->GetColor() : col);
		gReve->Redraw3D();
	}
*/
}
