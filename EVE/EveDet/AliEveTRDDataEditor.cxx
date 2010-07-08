#include "AliEveTRDDataEditor.h"
#include "AliEveTRDData.h"

ClassImp(AliEveTRDHitsEditor)
//ClassImp(AliEveTRDDigitsEditor)

///////////////////////////////////////////////////////////
////////////   AliEveTRDHitsEditor      ///////////////////
///////////////////////////////////////////////////////////
AliEveTRDHitsEditor::AliEveTRDHitsEditor(const TGWindow* p, Int_t width, Int_t height,
					 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options, back)
  ,fM(0x0)
{
  // Constructor.

  MakeTitle("TRD Points");
}

void AliEveTRDHitsEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveTRDHits*>(obj);

  // 	Float_t x, y, z;
  // 	for(int ihit=0; ihit<fM->GetN(); ihit++){
  // 		fM->GetPoint(ihit, x, y, z);
  // 		printf("%3d : x=%6.3f y=%6.3f z=%6.3f\n", ihit, x, y, z);
  // 	}
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDDigitsEditor /////////////////////
///////////////////////////////////////////////////////////
// AliEveTRDDigitsEditor::AliEveTRDDigitsEditor(const TGWindow* p, Int_t width, Int_t height,
// 					     UInt_t options, Pixel_t back) :
//   TGedFrame(p, width, height, options, back)
//   ,fM(0x0)
// {
//   // Constructor.
// 
//   MakeTitle("TRD Pixels");
// }
// 
// void AliEveTRDDigitsEditor::SetModel(TObject* obj)
// {
//   // Set model object.
// 
//   fM = dynamic_cast<AliEveTRDDigits*>(obj);
//   //fM->fParent->SpawnEditor();
// 
//   // 	printf("Chamber %d", fM->fParent->GetID());
//   // 	for (Int_t  row = 0;  row <  fM->fParent->GetRowMax();  row++)
//   // 		for (Int_t  col = 0;  col <  fM->fParent->GetColMax();  col++)
//   // 			for (Int_t time = 0; time < fM->fParent->GetTimeMax(); time++) {
//   // 				printf("\tA(%d %d %d) = %d\n", row, col, time, fM->fData.GetDataUnchecked(row, col, time));
//   // 			}
// }
// 
