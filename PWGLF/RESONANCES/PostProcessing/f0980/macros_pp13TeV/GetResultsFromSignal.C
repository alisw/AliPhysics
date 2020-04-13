double LevyTsallisF0(double *x, double *par){
// double mass = 0.970;
 double mass = par[3];
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / (par[1]*par[2]) ), -par[1] );
}
double LevyTsallisF01stMom(double *x, double *par){
// double mass = 0.970;
 double mass = par[3];
 return x[0] *
         x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / (par[1]*par[2]) ), -par[1] );
}


void GetResultsFromSignal(){
 TFile* fin = new TFile("./ResultFigures/SigOut_13TeV.root","read");

 const int nbins_mult = 6;
 double multmin[nbins_mult] = {
        0, 5, 10, 20, 50, 0 };
 double multmax[nbins_mult] = {
        5, 10, 20, 50, 100, 100 };

 int RainbowColor[nbins_mult]=
	{632, 801,402,412,434,602};
//	{632, 402,412,434,602 };

 TH1D* hYieldSysVtx[2][nbins_mult];
 TH1D* hYieldSysTrk[2][nbins_mult];
 TH1D* hYieldSysPID[2][nbins_mult];

 for(int i=0;i<2;i++){
	for(int j=0;j<nbins_mult;j++){
	}
 }

 TH1D* hY[nbins_mult];
 TH1D* hYSys[nbins_mult];

 for(int i=0;i<nbins_mult;i++){
	hY[i] = (TH1D*)fin->Get(Form("hy_%d",i));
	hY[i]->SetLineColor( RainbowColor[i] );
	hY[i]->SetMarkerColor( RainbowColor[i] );
	hY[i]->SetMarkerStyle( 20+i );
	hY[i]->Scale( pow( 2,nbins_mult-i-1 ) );

 	hYSys[i] = (TH1D*)hY[i]->Clone(0);
	hYSys[i]->SetFillStyle(0);
	hYSys[i]->SetFillColor( RainbowColor[i] );
 }


 TF1* fLevyTsallisF0[nbins_mult];
 TF1* fLevyTsallisF01stMom[nbins_mult];
 TFitResultPtr fLevyTsallisF0FitResults[nbins_mult];

 double yield[nbins_mult];
 double temp[nbins_mult];
 double expo[nbins_mult];
 double meanpt[nbins_mult];

 double yieldStat[nbins_mult];
 double tempStat[nbins_mult];
 double expoStat[nbins_mult];
 double meanptStat[nbins_mult];

 for(int i=0;i<nbins_mult;i++){
	fLevyTsallisF0[i] = new TF1("f1",LevyTsallisF0,0,13,4);
	fLevyTsallisF01stMom[i] = new TF1("f1",LevyTsallisF01stMom,0,13,4);

	fLevyTsallisF0[i]->SetLineColor( RainbowColor[i] );

	fLevyTsallisF0[i]->SetParLimits(0,0,100.0);
	fLevyTsallisF0[i]->SetParLimits(1,3,10);
	fLevyTsallisF0[i]->SetParLimits(2,0.1,0.5);

	fLevyTsallisF0[i]->FixParameter(3,0.98);

	fLevyTsallisF0FitResults[i] = hY[i]->Fit( fLevyTsallisF0[i], "s,i" ); 
 
	fLevyTsallisF01stMom[i]->SetParameters( fLevyTsallisF0[i]->GetParameters() );
	
//	yield[i] = fLevyTsallisF0[i]->GetParameter(0) / pow( 2,nbins_mult-i-1 );
	yield[i] = 0.0;
	yieldStat[i] = 0.0;

	meanpt[i] = 0.0;
	meanptStat[i] = 0.0;
	for(int j=0;j<hY[i]->GetNbinsX();j++){
		yield[i] += hY[i]->GetBinContent(j+1)*hY[i]->GetBinWidth(j+1);
		yieldStat[i] += pow( hY[i]->GetBinError(j+1)*hY[i]->GetBinWidth(j+1),2 );

		meanpt[i] += hY[i]->GetBinContent(j+1)*hY[i]->GetBinWidth(j+1)*hY[i]->GetBinCenter(j+1);
		meanptStat[i] += pow( hY[i]->GetBinError(j+1)*hY[i]->GetBinWidth(j+1)*hY[i]->GetBinCenter(j+1),2 );
	}
	yield[i] += fLevyTsallisF0[i]->Integral(0,0.3);
	yieldStat[i] += pow( fLevyTsallisF0[i]->IntegralError(0,0.3,
		fLevyTsallisF0[i]->GetParameters(),
		fLevyTsallisF0FitResults[i]->GetCovarianceMatrix().GetMatrixArray() ),2 );

	yieldStat[i] = sqrt( yieldStat[i] );

	yield[i] /= pow( 2,nbins_mult-i-1 );
	yieldStat[i] /= pow( 2,nbins_mult-i-1 );

	meanpt[i] += fLevyTsallisF01stMom[i]->Integral(0,0.3);
	meanptStat[i] += pow( fLevyTsallisF01stMom[i]->IntegralError(0,0.3,
		fLevyTsallisF01stMom[i]->GetParameters(),
		fLevyTsallisF0FitResults[i]->GetCovarianceMatrix().GetMatrixArray() ),2 );

	meanptStat[i] = sqrt( meanptStat[i] );

	meanpt[i] /= pow( 2,nbins_mult-i-1 );
	meanptStat[i] /= pow( 2,nbins_mult-i-1 );

	meanptStat[i] = meanpt[i]/yield[i] * sqrt( pow(meanptStat[i]/meanpt[i],2) + pow( yieldStat[i]/yield[i],2 ) );
	meanpt[i] /= yield[i];

	temp[i] = fLevyTsallisF0[i]->GetParameter(2);
	expo[i] = fLevyTsallisF0[i]->GetParameter(1);
//	meanpt[i] = fLevyTsallisF01stMom[i]->Integral(0,100) / yield[i] / pow( 2,nbins_mult-i-1 );
	
 
//        yieldStat[i] = fLevyTsallisF0[i]->GetParError(0) / pow( 2,nbins_mult-i-1 );
        tempStat[i] = fLevyTsallisF0[i]->GetParError(2);
        expoStat[i] = fLevyTsallisF0[i]->GetParError(1);
//	meanptStat[i] = fLevyTsallisF01stMom[i]->IntegralError(0,100,
//		fLevyTsallisF01stMom[i]->GetParameters(),
//		fLevyTsallisF0FitResults[i]->GetCovarianceMatrix().GetMatrixArray() ) / yield[i] / pow( 2,nbins_mult-i-1 );
 
 }

 TLegend* legAlldist = new TLegend(0.5,0.7,0.9,0.9);
 legAlldist->SetNColumns(2);
 TCanvas* c = new TCanvas("c","c",1200,800);
 gPad->SetLogy();
 for(int i=0;i<nbins_mult;i++){
	hY[i]->SetTitle("");
	hY[i]->SetMaximum(2e1);
	hY[i]->SetMinimum(1e-6);
	if(i==0) hY[i]->Draw("");
	else{ hY[i]->Draw("same"); }
	fLevyTsallisF0[i]->Draw("same");
//	hYSys[i]->Draw("same,e2");
	legAlldist->AddEntry( hY[i],Form("( %.0lf-%.0lf %% ) #times 2^{%d}",multmin[i],multmax[i], nbins_mult-i-1  ),"lp");
 }
 legAlldist->Draw();
 c->SaveAs("./GetResultsFromSignal/YieldAll.pdf");


 TLatex* latex = new TLatex();
 for(int i=0;i<nbins_mult;i++){
	c->cd();
	hY[i]->SetTitle(Form("(%.0lf-%.0lf)%% #times 2^{%d}",multmin[i],multmax[i],nbins_mult-i-1));
	hY[i]->Draw();
	fLevyTsallisF0[i]->Draw("same");

	latex->DrawLatexNDC(0.6,0.9,
		Form("#chi^{2} / NDF = %.3lf",fLevyTsallisF0[i]->GetChisquare()/fLevyTsallisF0[i]->GetNDF()));
	latex->DrawLatexNDC(0.6,0.8,
		Form("dN/dy = %.3lf #pm %.3lf",fLevyTsallisF0[i]->GetParameter(0)/pow( 2,nbins_mult-i-1 ),fLevyTsallisF0[i]->GetParError(0)/pow( 2,nbins_mult-i-1 )));
	latex->DrawLatexNDC(0.6,0.7,
		Form("T = %.3lf #pm %.3lf",fLevyTsallisF0[i]->GetParameter(2),fLevyTsallisF0[i]->GetParError(2)));
	latex->DrawLatexNDC(0.6,0.6,
		Form("n = %.3lf #pm %.3lf",fLevyTsallisF0[i]->GetParameter(1),fLevyTsallisF0[i]->GetParError(1)));

 	c->SaveAs(Form("./GetResultsFromSignal/Yield%d.pdf",i));
 }

 for(int i=0;i<nbins_mult;i++){
	cout << fLevyTsallisF0[i]->GetChisquare()/fLevyTsallisF0[i]->GetNDF() << ", " << yield[i] << ", " << yieldStat[i] << ", " << temp[i] << ", " << tempStat[i] << ", " << expo[i] << ", " << expoStat[i] << ", " << meanpt[i] << ", " << meanptStat[i] << endl;
 }

 const int RefMultBin = 10;
 const double RefMultPercMinBin[RefMultBin] = {
        1,      4,
        5,      5,
        5,      10,
        10,     10,
        20,     30 };


 const int index_ref2this[nbins_mult-1]={
        2, 3, 5, 8, 10};
//	2, 5, 8, 10};
//	5, 8, 10};

 double RealIncMult[RefMultBin]={
        26.02, 20.02, 16.17, 13.77, 12.04,
        10.02, 7.95, 6.32, 4.50, 2.55 };

 double RealIncMultStatErr[RefMultBin]={0};

 double RealIncMultSysUpErr[RefMultBin]={
        0.35, 0.27, 0.22, 0.19, 0.17,
        0.14, 0.11, 0.09, 0.07, 0.04 };

 double RealIncMultSysDownErr[RefMultBin]={
        0.29, 0.22, 0.18, 0.16, 0.14,
        0.11, 0.09, 0.07, 0.05, 0.03 };

 double RealIncMultSysErr[RefMultBin];
 for(int i=0;i<RefMultBin;i++){ RealIncMultSysErr[i] = ( RealIncMultSysUpErr[i]+RealIncMultSysDownErr[i] )/2.0; }

 double RealPi7TeV[RefMultBin]={
        10.035, 7.878, 6.459, 5.554, 4.892,
        4.138, 3.326, 2.699,1.989, 1.210 };

 double RealPi7TeVStatErr[RefMultBin]={
        0.007, 0.002, 0.001, 0.001, 0.001,
        0.001, 0.001, 0.001, 0.000, 0.000};

 double RealPi7TeVSysErr[RefMultBin]={
        0.519, 0.377, 0.306, 0.261, 0.228,
        0.192, 0.153, 0.125, 0.103, 0.086};

 double RealKaon7TeV[RefMultBin]={
        1.374, 1.062, 0.859, 0.731, 0.639,
        0.534, 0.423, 0.338, 0.243, 0.140 };

 double RealKaon7TeVStatErr[RefMultBin]={
        0.003, 0.001, 0.008, 0.006, 0.006,
        0.004, 0.003, 0.003, 0.002, 0.001};

 double RealKaon7TeVSysErr[RefMultBin]={
        0.093, 0.067, 0.054, 0.045, 0.039,
        0.033, 0.026, 0.021, 0.016, 0.012};

 double RealProton7TeV[RefMultBin]={
        5.488, 4.369, 3.599, 3.106, 2.741,
        2.316, 1.860, 1.491, 1.070, 0.586 }; //10^{-1}



 double RealPi13TeV[RefMultBin]={
	24.605, 19.016, 15.478, 13.242, 11.613,
	9.743, 7.780, 6.242, 4.531, 2.714};

 double RealPi13TeVStatErr[RefMultBin]={
	0.009, 0.004, 0.004, 0.003, 0.003,
	0.002, 0.002, 0.002, 0.001, 0.001};

 double RealPi13TeVSysErr[RefMultBin]={
	1.122, 0.856, 0.686, 0.583, 0.499,
	0.416, 0.327, 0.261, 0.186, 0.110};

for(int i=0;i<RefMultBin;i++){
	RealPi13TeV[i] /= 2.0;
	RealPi13TeVStatErr[i] /= 2.0;
	RealPi13TeVSysErr[i] /= 2.0;
 }

 double IncHad_ref2this[nbins_mult];
 double IncHad_ref2thisErr[nbins_mult];
 for(int i=0;i<nbins_mult-1;i++){
        IncHad_ref2this[i] = 0;
        IncHad_ref2thisErr[i] = 0;
        if(i==0){
                for(int j=0;j<index_ref2this[i];j++){
                        IncHad_ref2this[i] += RealIncMult[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        IncHad_ref2thisErr[i] += pow( RealIncMultSysErr[j] * RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        else{
                for(int j=index_ref2this[i-1];j<index_ref2this[i];j++){
                        IncHad_ref2this[i] += RealIncMult[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        IncHad_ref2thisErr[i] += pow( RealIncMultSysErr[j] * RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        IncHad_ref2thisErr[i] = sqrt( IncHad_ref2thisErr[i] );
 }
 IncHad_ref2this[nbins_mult-1] = 0;
 IncHad_ref2thisErr[nbins_mult-1] = 0;
 for(int j=0;j<RefMultBin;j++){
        IncHad_ref2this[nbins_mult-1] += RealIncMult[j]*RefMultPercMinBin[j]/100.0;
        IncHad_ref2thisErr[nbins_mult-1] += pow( RealIncMultSysErr[j]*RefMultPercMinBin[j]/100.0,2 );
 }
 IncHad_ref2thisErr[nbins_mult-1] = sqrt( IncHad_ref2thisErr[nbins_mult-1] );


/*
 double pi_ref2this[nbins_mult];
 double pi_ref2thisSys[nbins_mult];
 double pi_ref2thisStat[nbins_mult];
 for(int i=0;i<nbins_mult-1;i++){
        pi_ref2this[i] = 0;
        pi_ref2thisSys[i] = 0;
        pi_ref2thisStat[i] = 0;
        if(i==0){
                for(int j=0;j<index_ref2this[i];j++){
                        pi_ref2this[i] += RealPi7TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        pi_ref2thisSys[i] += pow( RealPi7TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        pi_ref2thisStat[i] += pow( RealPi7TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        else{
                for(int j=index_ref2this[i-1];j<index_ref2this[i];j++){
                        pi_ref2this[i] += RealPi7TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        pi_ref2thisSys[i] += pow( RealPi7TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        pi_ref2thisStat[i] += pow( RealPi7TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        pi_ref2thisSys[i] = sqrt( pi_ref2thisSys[i] );
        pi_ref2thisStat[i] = sqrt( pi_ref2thisStat[i] );
 }
 pi_ref2this[nbins_mult-1] = 0;
 pi_ref2thisSys[nbins_mult-1] = 0;
 pi_ref2thisStat[nbins_mult-1] = 0;
 for(int j=0;j<RefMultBin;j++){
        pi_ref2this[nbins_mult-1] += RealPi7TeV[j]*RefMultPercMinBin[j]/100.0;
        pi_ref2thisSys[nbins_mult-1] += pow( RealPi7TeVSysErr[j]*RefMultPercMinBin[j]/100.0,2 );
        pi_ref2thisStat[nbins_mult-1] += pow( RealPi7TeVStatErr[j]*RefMultPercMinBin[j]/100.0,2 );
 }
 pi_ref2thisSys[nbins_mult-1] = sqrt( pi_ref2thisSys[nbins_mult-1] );
 pi_ref2thisStat[nbins_mult-1] = sqrt( pi_ref2thisStat[nbins_mult-1] );
*/

 double pi_ref2this[nbins_mult];
 double pi_ref2thisSys[nbins_mult];
 double pi_ref2thisStat[nbins_mult];
 for(int i=0;i<nbins_mult-1;i++){
        pi_ref2this[i] = 0;
        pi_ref2thisSys[i] = 0;
        pi_ref2thisStat[i] = 0;
        if(i==0){
                for(int j=0;j<index_ref2this[i];j++){
                        pi_ref2this[i] += RealPi13TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        pi_ref2thisSys[i] += pow( RealPi13TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        pi_ref2thisStat[i] += pow( RealPi13TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        else{
                for(int j=index_ref2this[i-1];j<index_ref2this[i];j++){
                        pi_ref2this[i] += RealPi13TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        pi_ref2thisSys[i] += pow( RealPi13TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        pi_ref2thisStat[i] += pow( RealPi13TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        pi_ref2thisSys[i] = sqrt( pi_ref2thisSys[i] );
        pi_ref2thisStat[i] = sqrt( pi_ref2thisStat[i] );
 }
 pi_ref2this[nbins_mult-1] = 0;
 pi_ref2thisSys[nbins_mult-1] = 0;
 pi_ref2thisStat[nbins_mult-1] = 0;
 for(int j=0;j<RefMultBin;j++){
        pi_ref2this[nbins_mult-1] += RealPi13TeV[j]*RefMultPercMinBin[j]/100.0;
        pi_ref2thisSys[nbins_mult-1] += pow( RealPi13TeVSysErr[j]*RefMultPercMinBin[j]/100.0,2 );
        pi_ref2thisStat[nbins_mult-1] += pow( RealPi13TeVStatErr[j]*RefMultPercMinBin[j]/100.0,2 );
 }
 pi_ref2thisSys[nbins_mult-1] = sqrt( pi_ref2thisSys[nbins_mult-1] );
 pi_ref2thisStat[nbins_mult-1] = sqrt( pi_ref2thisStat[nbins_mult-1] );



 double kaon_ref2this[nbins_mult];
 double kaon_ref2thisSys[nbins_mult];
 double kaon_ref2thisStat[nbins_mult];
 for(int i=0;i<nbins_mult-1;i++){
        kaon_ref2this[i] = 0;
        if(i==0){
                for(int j=0;j<index_ref2this[i];j++){
                        kaon_ref2this[i] += RealKaon7TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        kaon_ref2thisSys[i] += pow( RealKaon7TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        kaon_ref2thisStat[i] += pow( RealKaon7TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        else{
                for(int j=index_ref2this[i-1];j<index_ref2this[i];j++){
                        kaon_ref2this[i] += RealKaon7TeV[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]);
                        kaon_ref2thisSys[i] += pow( RealKaon7TeVSysErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                        kaon_ref2thisStat[i] += pow( RealKaon7TeVStatErr[j]*RefMultPercMinBin[j] / (multmax[i]-multmin[i]),2 );
                }
        }
        kaon_ref2thisSys[i] = sqrt( kaon_ref2thisSys[i] );
        kaon_ref2thisStat[i] = sqrt( kaon_ref2thisStat[i] );
 }
 kaon_ref2this[nbins_mult-1] = 0;
 kaon_ref2thisSys[nbins_mult-1] = 0;
 kaon_ref2thisStat[nbins_mult-1] = 0;
 for(int j=0;j<RefMultBin;j++){
        kaon_ref2this[nbins_mult-1] += RealKaon7TeV[j]*RefMultPercMinBin[j]/100.0;
        kaon_ref2thisSys[nbins_mult-1] += pow( RealKaon7TeVSysErr[j]*RefMultPercMinBin[j]/100.0,2 );
        kaon_ref2thisStat[nbins_mult-1] += pow( RealKaon7TeVStatErr[j]*RefMultPercMinBin[j]/100.0,2 );
 }
 kaon_ref2thisSys[nbins_mult-1] = sqrt( kaon_ref2thisSys[nbins_mult-1] );
 kaon_ref2thisStat[nbins_mult-1] = sqrt( kaon_ref2thisStat[nbins_mult-1] );





 double GraphX[nbins_mult];
 double GraphXSys[nbins_mult];
 double GraphXStat[nbins_mult];
 double GraphY[nbins_mult];
 double GraphYSys[nbins_mult];
 double GraphYStat[nbins_mult];


//****************f0/f0_inc vs h/h_inc
 for(int i=0;i<nbins_mult;i++){
 	GraphX[i] = IncHad_ref2this[i] / IncHad_ref2this[nbins_mult-1];
	GraphXSys[i] = IncHad_ref2this[i] / IncHad_ref2this[nbins_mult-1] * sqrt(
		pow( IncHad_ref2thisErr[i]/IncHad_ref2this[i],2 )+
		pow( IncHad_ref2thisErr[nbins_mult-1]/IncHad_ref2this[nbins_mult-1],2 ) );
	GraphXStat[i] = 0.0;


	GraphY[i] = yield[i] / yield[nbins_mult-1];
//	GraphYSys[i] = Yield[1][i] / Yield[1][nbins_mult] * sqrt(
//		pow( YieldSys[1][i]/Yield[1][i],2 ) +
//		pow( YieldSys[1][nbins_mult]/Yield[1][nbins_mult],2 ) );
	GraphYSys[i] = 0.0;
	GraphYStat[i] = yield[i] / yield[nbins_mult-1] * sqrt(
		pow( yieldStat[i]/yield[i],2 ) +
		pow( yieldStat[nbins_mult-1]/yield[nbins_mult-1],2 ) );
 }


 TGraphErrors* gF0OverHadronSys = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXSys, GraphYSys );
 TGraphErrors* gF0OverHadronStat = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXStat, GraphYStat );



 gF0OverHadronStat->SetMarkerStyle(20);
 gF0OverHadronStat->SetMarkerColor(kBlack);
 gF0OverHadronStat->SetMarkerSize(1.5);

 gF0OverHadronStat->GetXaxis()->SetTitle("h / h^{inclusive}");
 gF0OverHadronStat->GetYaxis()->SetTitle("f_{0}(980) / f_{0}(980)^{inclusive}");

 gF0OverHadronStat->GetXaxis()->SetTitleSize(0.07);
 gF0OverHadronStat->GetYaxis()->SetTitleSize(0.07);
 gF0OverHadronStat->GetXaxis()->SetLabelSize(0.06);
 gF0OverHadronStat->GetYaxis()->SetLabelSize(0.06);

 gF0OverHadronStat->GetYaxis()->SetTitleOffset(0.85);
 gF0OverHadronStat->SetTitle("Yield Comparision in the pp at 13 TeV");

 gF0OverHadronSys->SetFillStyle(0);
 gF0OverHadronSys->SetFillColor(1);
 gF0OverHadronSys->SetMarkerColor(kBlack);

 TF1* fGuideLine = new TF1("f1","x",0,3.5);
 fGuideLine->SetLineColor(kBlue);
 fGuideLine->SetLineStyle(4);


 TCanvas* cf0_h_comp = new TCanvas("cf0_h_comp","cf0_h_comp",1200,800);
 cf0_h_comp->SetLeftMargin(0.15);
 cf0_h_comp->SetBottomMargin(0.15);
 gF0OverHadronStat->Draw("AP");
 gF0OverHadronSys->Draw("e2");
 fGuideLine->Draw("same");

 cf0_h_comp->SaveAs("./GetResultsFromSignal/f0_h_comp.pdf");
//***********************************






//************** Mean pT vs Mult
 for(int i=0;i<nbins_mult;i++){
        GraphX[i] = IncHad_ref2this[i];
        GraphXSys[i] = IncHad_ref2thisErr[i];
        GraphXStat[i] = 0.0;

        GraphY[i] = meanpt[i];
//        GraphYSys[i] = meanpTSys[i];
	GraphYSys[i] = 0.0;
	GraphYStat[i] = meanptStat[i];
 }

 TGraphErrors* gMeanpTOverHadronSys = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXSys, GraphYSys );
 TGraphErrors* gMeanpTOverHadronStat = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXStat, GraphYStat );

 gMeanpTOverHadronStat->SetMarkerStyle(20);
 gMeanpTOverHadronStat->SetMarkerColor(kBlack);
 gMeanpTOverHadronStat->SetMarkerSize(1.5);

 gMeanpTOverHadronStat->GetXaxis()->SetTitle("#frac{dN_{ch}}{dy}");
 gMeanpTOverHadronStat->GetYaxis()->SetTitle("<#font[12]{p}_{#font[22]{T}}> GeV/#font[12]{c}");

 gMeanpTOverHadronStat->GetXaxis()->SetTitleSize(0.07);
 gMeanpTOverHadronStat->GetYaxis()->SetTitleSize(0.07);
 gMeanpTOverHadronStat->GetXaxis()->SetLabelSize(0.06);
 gMeanpTOverHadronStat->GetYaxis()->SetLabelSize(0.06);

 gMeanpTOverHadronStat->GetYaxis()->SetTitleOffset(0.85);
 gMeanpTOverHadronStat->SetTitle("<#font[12]{p}_{#font[22]{T}}> of #font[22]{f}_{0}(980) in the pp at 13 TeV");

 gMeanpTOverHadronSys->SetFillStyle(0);
 gMeanpTOverHadronSys->SetFillColor(1);
 gMeanpTOverHadronSys->SetMarkerColor(kBlack);

 TCanvas* cpT_h_comp = new TCanvas("cpT_h_comp","cpT_h_comp",1200,800);
 cpT_h_comp->SetLeftMargin(0.15);
 cpT_h_comp->SetBottomMargin(0.15);
 gMeanpTOverHadronStat->Draw("AP");
 gMeanpTOverHadronSys->Draw("e2");

 cpT_h_comp->SaveAs("./GetResultsFromSignal/F0meanpT.pdf");
//*****************************

 
 
 
 
//**************** Ratio to pion
 for(int i=0;i<nbins_mult;i++){
        GraphX[i] = IncHad_ref2this[i];
        GraphXSys[i] = IncHad_ref2thisErr[i];
        GraphXStat[i] = 0.0;

        GraphY[i] = yield[i] / pi_ref2this[i];
//        GraphYSys[i] = Yield[1][i] / pi_ref2this[i] * sqrt(
//                pow( YieldSys[1][i] / Yield[1][i],2 )+
//                pow( pi_ref2thisSys[i] / pi_ref2this[i],2 ) );
	GraphYSys[i] = 0.0;
	GraphYStat[i] = yield[i] / pi_ref2this[i] * sqrt(
                pow( yieldStat[i] / yield[i],2 )+
                pow( pi_ref2thisStat[i] / pi_ref2this[i],2 ) );

 }



 TGraphErrors* gF0OverPion7TeVSys = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXSys, GraphYSys );
 TGraphErrors* gF0OverPion7TeVStat = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXStat, GraphYStat );

 gF0OverPion7TeVStat->SetMarkerStyle(20);
 gF0OverPion7TeVStat->SetMarkerColor(kBlack);
 gF0OverPion7TeVStat->SetMarkerSize(1.5);

 gF0OverPion7TeVStat->GetXaxis()->SetTitle("#frac{dN_{ch}}{dy}");
 gF0OverPion7TeVStat->GetYaxis()->SetTitle("#font[22]{f}_{0}(980) / (#pi^{+} + #pi^{-})");

 gF0OverPion7TeVStat->GetXaxis()->SetTitleSize(0.07);
 gF0OverPion7TeVStat->GetYaxis()->SetTitleSize(0.07);
 gF0OverPion7TeVStat->GetXaxis()->SetLabelSize(0.06);
 gF0OverPion7TeVStat->GetYaxis()->SetLabelSize(0.06);

 gF0OverPion7TeVStat->GetYaxis()->SetTitleOffset(0.85);
 gF0OverPion7TeVStat->SetTitle("Particle Ratio (#font[22]{f}_{0}(980) for 13 TeV and #pi for 7 TeV)");

 gF0OverPion7TeVSys->SetFillStyle(0);
 gF0OverPion7TeVSys->SetFillColor(1);
 gF0OverPion7TeVSys->SetMarkerColor(kBlack);

 TCanvas* cFPratio = new TCanvas("cFPratio","cFPratio",1200,800);
 cFPratio->SetLeftMargin(0.15);
 cFPratio->SetBottomMargin(0.15);
 gF0OverPion7TeVStat->Draw("AP");
 gF0OverPion7TeVSys->Draw("e2");
//*****************************




//**************** Ratio to kaon
 for(int i=0;i<nbins_mult;i++){
        GraphX[i] = IncHad_ref2this[i];
        GraphXSys[i] = IncHad_ref2thisErr[i];
        GraphXStat[i] = 0.0;

        GraphY[i] = yield[i] / kaon_ref2this[i];
/*
        GraphYSys[i] = Yield[1][i] / kaon_ref2this[i] * sqrt(
                pow( YieldSys[1][i] / Yield[1][i],2 )+
                pow( kaon_ref2thisSys[i] / kaon_ref2this[i],2 ) );
*/
	GraphYSys[i] = 0.0;
	GraphYStat[i] = yield[i] / kaon_ref2this[i] * sqrt(
                pow( yieldStat[i] / yield[i],2 )+
                pow( kaon_ref2thisStat[i] / kaon_ref2this[i],2 ) );
 }

 TGraphErrors* gF0OverKaon7TeVSys = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXSys, GraphYSys );
 TGraphErrors* gF0OverKaon7TeVStat = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXStat, GraphYStat );

 gF0OverKaon7TeVStat->SetMarkerStyle(20);
 gF0OverKaon7TeVStat->SetMarkerColor(kBlack);
 gF0OverKaon7TeVStat->SetMarkerSize(1.5);

 gF0OverKaon7TeVStat->GetXaxis()->SetTitle("#frac{dN_{ch}}{dy}");
 gF0OverKaon7TeVStat->GetYaxis()->SetTitle("#font[22]{f}_{0}(980) / (K^{+} + K^{-})");

 gF0OverKaon7TeVStat->GetXaxis()->SetTitleSize(0.07);
 gF0OverKaon7TeVStat->GetYaxis()->SetTitleSize(0.07);
 gF0OverKaon7TeVStat->GetXaxis()->SetLabelSize(0.06);
 gF0OverKaon7TeVStat->GetYaxis()->SetLabelSize(0.06);

 gF0OverKaon7TeVStat->GetYaxis()->SetTitleOffset(0.85);
 gF0OverKaon7TeVStat->SetTitle("Particle Ratio (#font[22]{f}_{0}(980) for 13 TeV and K for 7 TeV)");

 gF0OverKaon7TeVSys->SetFillStyle(0);
 gF0OverKaon7TeVSys->SetFillColor(1);
 gF0OverKaon7TeVSys->SetMarkerColor(kBlack);

 TCanvas* cFKratio = new TCanvas("cFKratio","cFKratio",1200,800);
 cFKratio->SetLeftMargin(0.15);
 cFKratio->SetBottomMargin(0.15);
 gF0OverKaon7TeVStat->Draw("AP");
 gF0OverKaon7TeVSys->Draw("e2");
//*****************************



//Double Ratio with pion***************
 for(int i=0;i<nbins_mult;i++){
        GraphX[i] = IncHad_ref2this[i];
        GraphXSys[i] = IncHad_ref2thisErr[i];
        GraphXStat[i] = 0.0;

	GraphY[i] = ( yield[i] / yield[nbins_mult-1] ) /
		( pi_ref2this[i] / pi_ref2this[nbins_mult-1] );
	GraphYSys[i] = 0.0;
	GraphYStat[i] = 0.1*GraphY[i];
 }

 TGraphErrors* gF0DoubleRatioWithPionSys = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXSys, GraphYSys );
 TGraphErrors* gF0DoubleRatioWithPionStat = new TGraphErrors( nbins_mult, GraphX, GraphY, GraphXStat, GraphYStat );

 gF0DoubleRatioWithPionStat->SetMarkerStyle(20);
 gF0DoubleRatioWithPionStat->SetMarkerColor(kBlack);
 gF0DoubleRatioWithPionStat->SetMarkerSize(1.5);

 gF0DoubleRatioWithPionStat->GetXaxis()->SetTitle("#frac{dN_{ch}}{dy}");
 gF0DoubleRatioWithPionStat->GetYaxis()->SetTitle("2#font[22]{f}_{0}(980) / (#pi^{+} + #pi^{-})");

 gF0DoubleRatioWithPionStat->GetXaxis()->SetTitleSize(0.07);
 gF0DoubleRatioWithPionStat->GetYaxis()->SetTitleSize(0.07);
 gF0DoubleRatioWithPionStat->GetXaxis()->SetLabelSize(0.06);
 gF0DoubleRatioWithPionStat->GetYaxis()->SetLabelSize(0.06);

 gF0DoubleRatioWithPionStat->GetYaxis()->SetTitleOffset(0.85);
 gF0DoubleRatioWithPionStat->SetTitle("Particle Ratio (#font[22]{f}_{0}(980) for 13 TeV and #pi for 7 TeV)");

 gF0DoubleRatioWithPionSys->SetFillStyle(0);
 gF0DoubleRatioWithPionSys->SetFillColor(1);
 gF0DoubleRatioWithPionSys->SetMarkerColor(kBlack);

 TCanvas* cDPion = new TCanvas("cDPion","cDPion",1200,800);
 cDPion->SetLeftMargin(0.15);
 cDPion->SetBottomMargin(0.15);
 gF0DoubleRatioWithPionStat->Draw("AP");
 gF0DoubleRatioWithPionSys->Draw("e2");
//****************************************************



}
