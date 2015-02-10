/// \file AnalyzeLaser.C


gSystem->Load("libSTAT");
TStatToolkit stat;
Int_t npoints;
Double_t chi2;
TVectorD vec;
TMatrixD mat;


TObjArray * array = AliTPCCalibViewerGUI::ShowGUI("laserTree.root");
AliTPCCalibViewerGUI * viewer = (AliTPCCalibViewerGUI*)array->At(0);
TTree * tree = viewer->GetViewer()->GetTree();
TFile fp("/data/calib/CalibTreePulser_run33834_Cside.root");
tree->AddFriend(treePulser,"P.");


tree->SetAlias("dt","(sector%36>30)*2+(sector<36)*0.3");
tree->SetAlias("T","T0_100_220.fElements-P..StandardTime0.fElements");
tree->SetAlias("Tm","T0_100_220_Median.fElements");
tree->SetAlias("Q","Q_100_220.fElements");
tree->SetAlias("Qm","Q_100_220_Median.fElements");

tree->SetAlias("Qcut","abs(Q/Qm-2)<1.4&&Q>6&&Q<200");
tree->SetAlias("Tcut","abs(T-Tm)<2");



TString strSector="";
{
  for (Int_t isec=54;isec<71;isec+=1){
    if (isec!=64) {
      strSector+="(sector==";
      strSector+=isec;
      strSector+=")++";
      strSector+="(lx.fElements-195.)*(sector==";
      strSector+=isec;
      strSector+=")++";
      strSector+="((lx.fElements-195)^2)*(sector==";
      strSector+=isec;
      strSector+=")++";
    }
  }
}


TCut cutA="Tcut&&Qcut&&sector%36>17";

TString *strFitG =stat.FitPlane(tree,"T+dt","gx.fElements++gy.fElements",cutA,chi2,npoints,vec,mat);

TString *strFitGL =stat.FitPlane(tree,"T+dt","lx.fElements++ly.fElements++gx.fElements++gy.fElements",cutA,chi2,npoints,vec,mat);

TString *strFitGL2 =stat.FitPlane(tree,"T+dt","lx.fElements++ly.fElements++gx.fElements++gy.fElements++lx.fElements^2++ly.fElements^2",cutA,chi2,npoints,vec,mat);

TString *strFitGLA =stat.FitPlane(tree,"T+dt",strSector+"lx.fElements++ly.fElements++gx.fElements++gy.fElements++lx.fElements^2++ly.fElements^2",cutA,chi2,npoints,vec,mat);


tree->SetAlias("tfitG",strFitG->Data())
tree->SetAlias("tfitGL",strFitGL->Data())
tree->SetAlias("tfitGL2",strFitGL2->Data())
tree->SetAlias("tfitGLA",strFitGLA->Data())


