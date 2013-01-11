#include "AliPSQAVisualization.h"

void PSQAV(){
  //
  // make nice plots for presentations of the Phys.Sel. QA
  //
  AliPSQAVisualization * QAVobject1 = new AliPSQAVisualization;
  //
  Bool_t useAutoScale    = kTRUE;
  Bool_t savePDFs        = kTRUE;
  Bool_t drawSelected    = kTRUE;  // draw all plots on the screen
  Bool_t drawOverPlot    = kTRUE;  // if true, draw to canvas to see
  Bool_t saveOverPlotPdf = kTRUE;  // if true, write to pdf file
  Bool_t saveOverPlotEps = kTRUE;  // if true, write to eps file
  //
  // configuration -- choose the input directory and file
  //
  QAVobject1->SetROOTInput("../QAoutputPerPeriod/12d_Pass1/LHC12d_PSQA.root");
  QAVobject1->InitializeSelectedPlots("../QAoutputPerPeriod/12d_Pass1/selectedPlots.list");
	//QAVobject1->InitializeSelectedPlots("plots.list");
  QAVobject1->SetOutDirectory("../QAoutputPerPeriod/12d_Pass1/plots");
  QAVobject1->SetOutPDFName("PlotsLHC12d_pass.pdf");
  QAVobject1->SetOutEPSName("PlotsLHC12d_pass.eps");
  QAVobject1->SetOverPlotTitle("");
  //
  // settings
  //
  QAVobject1->SetSavePDFs(savePDFs);
  QAVobject1->SetDrawSelected(drawSelected);
  QAVobject1->SetDrawOverPlot(drawOverPlot);
  QAVobject1->SetSaveOverPlotPDF(saveOverPlotPdf);
  QAVobject1->SetSaveOverPlotEPS(saveOverPlotEps); // if true, write to eps file
  QAVobject1->SetUseColorArray(kTRUE);
  //
  QAVobject1->InitializeColorArray("colorArray.list"); // Set the color enums for the overplot
    /// Use Auto Scale? See README.
  QAVobject1->SetScaleAuto(useAutoScale);
  if (QAVobject1->GetScaleAuto()) { 
    QAVobject1->SetScaleAutoDivMin(1.1);
    QAVobject1->SetScaleAutoMultMax(1.1);
  } else {
    QAVobject1->SetScaleManMin(0.0);
    QAVobject1->SetScaleManMax(1.1);
  }
  //
	cout << "Here we are?" << endl;
	
  QAVobject1->PostProcessQA();

 cout << "Done?" << endl;
};
