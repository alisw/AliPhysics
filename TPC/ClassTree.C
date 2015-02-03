/// \file ClassTree.C
/// \brief Macro generated from canvas: ClassTree
///
/// \date Tue Jun 1 17:01:38 1999
///
/// ROOT version 2.21/07

void ClassTree()
{

   TCanvas *ClassTree = new TCanvas("ClassTree", "",186,135,594,449);
   ClassTree->SetHighLightColor(2);
   ClassTree->Range(0,5,20,20);
   ClassTree->SetFillColor(33);
   ClassTree->SetBorderSize(2);
   TLine *line = new TLine(0.5,18.15,4.4,18.15);
   line->Draw();
   line = new TLine(4.4,17.725,4.4,18.575);
   line->Draw();
   
   TPaveLabel *pl = new TPaveLabel(1,17.895,4.205,18.405,"TArray","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(0.5,16.875,1,16.875);
   line->Draw();
   
   pl = new TPaveLabel(1,16.62,4.205,17.13,"TAttFill","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(0.5,16.025,1,16.025);
   line->Draw();
   
   pl = new TPaveLabel(1,15.77,4.205,16.28,"TAttLine","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(0.5,15.175,1,15.175);
   line->Draw();
   
   pl = new TPaveLabel(1,14.92,4.205,15.43,"TAttMarker","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(0.5,9.775,4.4,9.775);
   line->Draw();
   line = new TLine(4.4,7.65,4.4,12.325);
   line->Draw();
   
   pl = new TPaveLabel(1,9.52,4.205,10.03,"TObject","br");
   pl->SetFillColor(5);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,12.325,4.9,12.325);
   line->Draw();
   
   pl = new TPaveLabel(4.9,12.07,8.105,12.58,"AliArrayI","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,11.475,4.9,11.475);
   line->Draw();
   
   pl = new TPaveLabel(4.9,11.22,8.105,11.73,"AliArrayS","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,18.575,4.9,18.575);
   line->Draw();
   
   pl = new TPaveLabel(4.9,18.32,8.105,18.83,"TArrayI","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,17.725,4.9,17.725);
   line->Draw();
   
   pl = new TPaveLabel(4.9,17.47,8.105,17.98,"TArrayS","br");
   pl->SetFillColor(30);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,10.2,8.3,10.2);
   line->Draw();
   
   pl = new TPaveLabel(4.9,9.945,8.105,10.455,"TCollection","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(4.4,7.65,8.3,7.65);
   line->Draw();
   line = new TLine(8.3,6.8,8.3,8.075);
   line->Draw();
   
   pl = new TPaveLabel(4.9,7.395,8.105,7.905,"TNamed","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(8.3,8.075,12.2,8.075);
   line->Draw();
   
   pl = new TPaveLabel(8.8,7.82,12.005,8.33,"AliSegmentArray","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(8.3,10.2,12.2,10.2);
   line->Draw();
   
   pl = new TPaveLabel(8.8,9.945,12.005,10.455,"TSeqCollection","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(8.3,6.8,8.8,6.8);
   line->Draw();
   
   pl = new TPaveLabel(8.8,6.545,12.005,7.055,"TTree","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(12.2,8.075,12.7,8.075);
   line->Draw();
   
   pl = new TPaveLabel(12.7,7.82,15.905,8.33,"AliTPCClustersArray","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(12.2,10.2,16.1,10.2);
   line->Draw();
   
   pl = new TPaveLabel(12.7,9.945,15.905,10.455,"TObjArray","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   line = new TLine(16.1,10.2,16.6,10.2);
   line->Draw();
   
   pl = new TPaveLabel(16.6,9.945,19.805,10.455,"TClonesArray","br");
   pl->SetFillColor(18);
   pl->SetTextSize(0.9);
   pl->Draw();
   
   pl = new TPaveLabel(0.1,19.1,18.2,19.9,"*AliSegmet:*AliSegmentArray:*AliArrayI:*AliArrayS:TTree:*TObjArray","br");
   pl->SetFillColor(42);
   pl->SetTextSize(0.7);
   pl->Draw();
   line = new TLine(11.4041,6.8,14.3025,10.2);
   line->SetLineColor(6);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(11.4842,6.8,14.3025,10.2);
   line->SetLineColor(6);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(6.5025,12.325,6.5025,18.575);
   line->SetLineColor(4);
   line->SetLineStyle(2);
   line->Draw();
   line = new TLine(6.5025,11.475,6.5025,17.725);
   line->SetLineColor(4);
   line->SetLineStyle(2);
   line->Draw();
   line = new TLine(10.4025,6.8,2.6025,16.025);
   line->SetLineColor(4);
   line->SetLineStyle(2);
   line->Draw();
   line = new TLine(10.4025,6.8,2.6025,16.875);
   line->SetLineColor(4);
   line->SetLineStyle(2);
   line->Draw();
   line = new TLine(10.4025,6.8,2.6025,15.175);
   line->SetLineColor(4);
   line->SetLineStyle(2);
   line->Draw();
   TArrow *arrow = new TArrow(5.43417,10.2,6.5025,10.2,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(6.85861,10.2,2.6025,9.775,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(9.60125,8.075,14.3025,10.2,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(10.1354,8.075,6.5025,12.325,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(11.2037,8.075,10.4025,6.8,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(13.9019,10.2,2.6025,9.775,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   arrow = new TArrow(19.2708,10.2,14.3025,10.2,0.008,"|>");
   arrow->SetFillColor(2);
   arrow->SetFillStyle(1001);
   arrow->SetLineColor(2);
   arrow->Draw();
   ClassTree->Modified();
   ClassTree->cd();
}
