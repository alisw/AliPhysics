Int_t readPhojet() 
{
// read PHOJET INEL data from file

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read PHOJET INEL data from file" <<endl;
cout << "Number of bins: " << binsPhojet <<endl;
cout << "Filename:       " << filenamePhojet <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

TFile* filePhojet = new TFile(filenamePhojet,"OPEN");
Double_t nEventsPhojet = ((TH1F*)filePhojet->Get("eventsINEL"))->GetBinContent(1);

TH1F* histPhojet = (TH1F*)filePhojet->Get("ptINEL");

int bins = 0;
for (int i=0; i < histPhojet->GetNbinsX(); i++) {
    if (histPhojet->GetBinContent(i) > 0 ) {
        centerPtPhojet[bins]    = histPhojet->GetBinCenter(i);
        ptPhojet[bins]          = centerPtPhojet[bins];
        widthPtPhojet[bins]     = histPhojet->GetBinWidth(i);        
        errPtPhojet[bins]       = widthPtPhojet[bins] / 2.0;
        lowPtPhojet[bins]       = centerPtPhojet[bins] - errPtPhojet[bins];
        highPtPhojet[bins]      = centerPtPhojet[bins] + errPtPhojet[bins];
        
        inelPhojet[bins]        = histPhojet->GetBinContent(i) / (nEventsPhojet * etaRange * 2 * M_PI * ptPhojet[bins]);
        errInelPhojet[bins]     = histPhojet->GetBinError(i) / (nEventsPhojet * etaRange * 2 * M_PI * ptPhojet[bins]);
        lowErrInelPhojet[bins]  = inelPhojet[bins] - errInelPhojet[bins];
        highErrInelPhojet[bins] = inelPhojet[bins] + errInelPhojet[bins];
        relErrInelPhojet[bins]  = errInelPhojet[bins] / inelPhojet[bins];
     
        centerPt2PiPtPhojet[bins]    = histPhojet->GetBinCenter(i);
        pt2PiPtPhojet[bins]          = centerPt2PiPtPhojet[bins];
        widthPt2PiPtPhojet[bins]     = histPhojet->GetBinWidth(i);        
        errPt2PiPtPhojet[bins]       = widthPt2PiPtPhojet[bins] / 2.0;
        lowPt2PiPtPhojet[bins]       = centerPt2PiPtPhojet[bins] - errPt2PiPtPhojet[bins];
        highPt2PiPtPhojet[bins]      = centerPt2PiPtPhojet[bins] + errPt2PiPtPhojet[bins];        
        inel2PiPtPhojet[bins]        = histPhojet->GetBinContent(i) / (nEventsPhojet * etaRange);
        errInel2PiPtPhojet[bins]     = histPhojet->GetBinError(i) / (nEventsPhojet * etaRange);
        lowErrInel2PiPtPhojet[bins]  = inel2PiPtPhojet[bins] - errInel2PiPtPhojet[bins];
        highErrInel2PiPtPhojet[bins] = inel2PiPtPhojet[bins] + errInel2PiPtPhojet[bins];
        relErrInel2PiPtPhojet[bins]  = errInel2PiPtPhojet[bins] / inel2PiPtPhojet[bins];
        
        cout << "ptPhojet[" << bins << "]              = " << ptPhojet[bins] <<endl;
        cout << "   centerPtPhojet[" << bins << "]     = " << centerPtPhojet[bins] <<endl;
        cout << "   widthPtPhojet[" << bins << "]      = " << widthPtPhojet[bins] <<endl;
        cout << "   errPtPhojet[" << bins << "]        = " << errPtPhojet[bins] <<endl;
        cout << "   lowPtPhojet[" << bins << "]        = " << lowPtPhojet[bins] <<endl;
        cout << "   highPtPhojet[" << bins << "]       = " << highPtPhojet[bins] <<endl;
        cout << "inelPhojet[" << bins << "]            = " << inelPhojet[bins] <<endl;        
        cout << "errInelPhojet[" << bins << "]         = " << errInelPhojet[bins] <<endl;
        cout << "   lowErrInelPhojet[" << bins << "]   = " << lowErrInelPhojet[bins] <<endl;
        cout << "   highErrInelPhojet[" << bins << "]  = " << highErrInelPhojet[bins] <<endl;
        cout << "   relErrInelPhojet[" << bins << "]   = " << relErrInelPhojet[bins] <<endl;
        cout << endl;
        bins++;
    }
}

filePhojet->Close();

cout << "================================" << endl;
cout << endl;
cout << "Finished reading PHOJET INEL data" <<endl;
cout << "Number of bins read: " << bins <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return bins;
}