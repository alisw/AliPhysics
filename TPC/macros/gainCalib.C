/*
 Small macro to show perforamance of the gain calibration
 tHe list of the  files with tpc calib tracks is supposed to be in cosmic.txt file 
 Supposing the 


 gSystem->Load("libSTAT.so")
 .x ~/rootlogon.C
 .L $ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+
 AliXRDPROOFtoolkit tool;
 //gSystem->Load("/usr/local/grid/XRootd/GSI/lib/libXrdClient.so");
 //TProof * proof = TProof::Open("miranov@lxgrid2.gsi.de");

*/

TChain * chain = 0;
TChain * chaing = 0;

TStatToolkit toolkit;

void Init(){

  AliCDBManager::Instance()->SetRun(0);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTPCClusterParam *clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  AliTPCClusterParam::SetInstance(clparam);  
  //
  AliXRDPROOFtoolkit tool;
  chain = tool.MakeChain("cosmic.txt","dEdxT",0,1000)
  chain->Lookup();
  chaing = tool.MakeChain("cosmic.txt","TrackG",0,1000)
  chaing->Lookup();
  //
  chain->SetAlias("dr","(250-abs(meanPos.fElements[4]))/250.");
  chain->SetAlias("tz","(0+abs(parZ.fElements[1]))");
  chain->SetAlias("ty","(0+abs(parY.fElements[1]))");
  chain->SetAlias("corrg","sqrt((1+ty^2)*(1+tz^2))");
}





void MakeFits(){
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD param;
  TMatrixD covar;
  Int_t npoints;
  TString *strq0 = toolkit.FitPlane(chain,"dedxQ.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==0",chi2,npoints,param,covar,0,20000);
  TString *strq1 = toolkit.FitPlane(chain,"dedxQ.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==1",chi2,npoints,param,covar,0,20000);
  TString *strq2 = toolkit.FitPlane(chain,"dedxQ.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==2",chi2,npoints,param,covar,0,20000);
  TString *strm0 = toolkit.FitPlane(chain,"dedxM.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==0",chi2,npoints,param,covar,0,20000);
  TString *strm1 = toolkit.FitPlane(chain,"dedxM.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==1",chi2,npoints,param,covar,0,20000);
  TString *strm2 = toolkit.FitPlane(chain,"dedxM.fElements[2]","dr++tz++ty++dr*tz++dr*ty++ty*tz","IPad==2",chi2,npoints,param,covar,0,20000);  
  chain->SetAlias("normqt0",strq0.Data());
  chain->SetAlias("normqt1",strq1.Data());
  chain->SetAlias("normqt2",strq2.Data());
  chain->SetAlias("normqm0",strm0.Data());
  chain->SetAlias("normqm1",strm1.Data());
  chain->SetAlias("normqm2",strm2.Data());
}

TFile fqplot("qplot.root","update");


void MakePlotsQS(){
  
  chain->Draw("dedxQ.fElements[2]:sector>>hisQ_sector(36,0,36)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/corrg:sector>>hisQ_sector_corrg(36,0,36)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/AliTPCClusterParam::SQnorm(0,0,dr,ty,tz):sector>>hisQ_sector_corrcal(36,0,36)","IPad==0&&P>1","prof*",1000000);

  TProfile * profq_sector = (TProfile*) gROOT->FindObject("hisQ_sector");
  TProfile * profq_sector_corrg = (TProfile*) gROOT->FindObject("hisQ_sector_corrg");
  TProfile * profq_sector_corrcal = (TProfile*) gROOT->FindObject("hisQ_sector_corrcal");
  profq_sector->SetMarkerStyle(22);
  profq_sector_corrg->SetMarkerStyle(24);
  profq_sector_corrcal->SetMarkerStyle(26);
  //
  profq_sector->SetXTitle("Sector number");
  profq_sector->SetYTitle("Mean amplitude");  
  profq_sector->SetMinimum(0);
  profq_sector->Draw();
  profq_sector_corrg->Draw("same");
  profq_sector_corrcal->Draw("same");
  fqplot.cd();
  profq_sector->Write("qt_sector_0");
  profq_sector_corrg->Write("qt_sector_1");
  profq_sector_corrcal->Write("qt_sector_2");
  gPad->Write("qt_sector");
}
  

void MakePlotsTY(){
  //
  //
  //
  chain->Draw("dedxQ.fElements[2]:ty>>hisQ_ty(20,0,1.5)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/corrg:ty>>hisQ_ty_corrg(20,0,1.5)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/AliTPCClusterParam::SQnorm(0,0,dr,ty,tz):ty>>hisQ_ty_corrcal(20,0,1.5)","IPad==0&&P>1","prof*",1000000);

  TProfile * profq_ty = (TProfile*) gROOT->FindObject("hisQ_ty");
  TProfile * profq_ty_corrg = (TProfile*) gROOT->FindObject("hisQ_ty_corrg");
  TProfile * profq_ty_corrcal = (TProfile*) gROOT->FindObject("hisQ_ty_corrcal");
  profq_ty->SetMarkerStyle(22);
  profq_ty_corrg->SetMarkerStyle(24);
  profq_ty_corrcal->SetMarkerStyle(26);
  //
  profq_ty->SetXTitle("tan(#phi)");
  profq_ty->SetYTitle("Mean amplitude");  
  profq_ty->SetMinimum(0);
  profq_ty->Draw();
  profq_ty_corrg->Draw("same");
  profq_ty_corrcal->Draw("same");
  fqplot.cd();
  profq_ty->Write("qt_ty_0");
  profq_ty_corrg->Write("qt_ty_1");
  profq_ty_corrcal->Write("qt_ty_2");
  gPad->Write("qt_ty");
}


void MakePlotsTZ(){
  //
  //
  //
  chain->Draw("dedxQ.fElements[2]:tz>>hisQ_tz(20,0,1.5)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/corrg:tz>>hisQ_tz_corrg(20,0,1.5)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/AliTPCClusterParam::SQnorm(0,0,dr,ty,tz):tz>>hisQ_tz_corrcal(20,0,1.5)","IPad==0&&P>1","prof*",1000000);

  TProfile * profq_tz = (TProfile*) gROOT->FindObject("hisQ_tz");
  TProfile * profq_tz_corrg = (TProfile*) gROOT->FindObject("hisQ_tz_corrg");
  TProfile * profq_tz_corrcal = (TProfile*) gROOT->FindObject("hisQ_tz_corrcal");
  profq_tz->SetMarkerStyle(22);
  profq_tz_corrg->SetMarkerStyle(24);
  profq_tz_corrcal->SetMarkerStyle(26);
  //
  profq_tz->SetXTitle("tan(#theta)");
  profq_tz->SetYTitle("Mean amplitude");  
  profq_tz->SetMinimum(0);
  profq_tz->Draw();
  profq_tz_corrg->Draw("same");
  profq_tz_corrcal->Draw("same");
  fqplot.cd();
  profq_tz->Write("qt_tz_0");
  profq_tz_corrg->Write("qt_tz_1");
  profq_tz_corrcal->Write("qt_tz_2");
  gPad->Write("qt_tz");
}

void MakePlotsDR(){
  //
  //
  //
  chain->Draw("dedxQ.fElements[2]:dr>>hisQ_dr(20,0,1.)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/corrg:dr>>hisQ_dr_corrg(20,0,1.)","IPad==0&&P>1","prof*",1000000);
  chain->Draw("dedxQ.fElements[2]/AliTPCClusterParam::SQnorm(0,0,dr,ty,dr):dr>>hisQ_dr_corrcal(20,0,1.)","IPad==0&&P>1","prof*",1000000);

  TProfile * profq_dr = (TProfile*) gROOT->FindObject("hisQ_dr");
  TProfile * profq_dr_corrg = (TProfile*) gROOT->FindObject("hisQ_dr_corrg");
  TProfile * profq_dr_corrcal = (TProfile*) gROOT->FindObject("hisQ_dr_corrcal");
  profq_dr->SetMarkerStyle(22);
  profq_dr_corrg->SetMarkerStyle(24);
  profq_dr_corrcal->SetMarkerStyle(26);
  //
  profq_dr->SetXTitle("drift length(unit)");
  profq_dr->SetYTitle("Mean amplitude");  
  profq_dr->SetMinimum(0);
  profq_dr->Draw();
  profq_dr_corrg->Draw("same");
  profq_dr_corrcal->Draw("same");
  fqplot.cd();
  profq_dr->Write("qt_dr_0");
  profq_dr_corrg->Write("qt_dr_1");
  profq_dr_corrcal->Write("qt_dr_2");
  gPad->Write("qt_dr");
}



void MakePlotsQ(){
  //
  //
  //
  //
  chaing->Draw("Track.fdEdx>>his0dedx(100,0,200)","Track.fN>100&&abs(Track.P()-15)<3","",20000);
  chaing->Draw("Track.CookdEdxNorm(0.02,0.6,0,0,160)/3.71>>his0dedxnormQ(100,0,200)","Track.fN>100&&abs(Track.P()-15)<3","",20000);
  TH1F * his0dedx = gROOT->FindObject("his0dedx");
  TH1F * his0dedxnormq = gROOT->FindObject("his0dedxnormQ");
  his0dedxnormq->SetXTitle("dEdx (rel. unit)")
  his0dedxnormq->Draw();
  his0dedx->Draw("same");
 

  chaing->Draw("Track.CookdEdxNorm(0.02,0.6,0,0,160)/3.71:Track.GetP()>>hispdedxnorm(10,5,100)","Track.fN>100","prof",20000);


  chaing->Draw("Track.fdEdx:sector0>>hisdedx(36,0,36)","Track.fN>60","prof",10000);
  chaing->Draw("Track.CookdEdxNorm(0.02,0.6,0,0,160)/3.71:sector0>>hisdedxnormQ(36,0,36)","Track.fN>60","prof",10000);

}












