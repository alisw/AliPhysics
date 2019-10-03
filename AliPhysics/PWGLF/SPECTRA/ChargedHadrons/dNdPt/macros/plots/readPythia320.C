Int_t readPythia320() 
{
// read PYTHIA 320 INEL data from file

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read PHYTIA 320 INEL data from file" <<endl;
cout << "Number of bins: " << binsPythia320 <<endl;
cout << "Filename:       " << filenamePythia320 <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

TFile* filePythia320 = new TFile(filenamePythia320,"OPEN");
nEventsPythia320 = ((TH1F*)filePythia320->Get("eventsINEL"))->GetBinContent(1);

TH1F* histPythia320 = (TH1F*)filePythia320->Get("ptINEL");

int bins = 0;
for (int i=0; i < histPythia320->GetNbinsX(); i++) {
    if (histPythia320->GetBinContent(i) > 0 ) {
        centerPtPythia320[bins]    = histPythia320->GetBinCenter(i);
        ptPythia320[bins]          = centerPtPythia320[bins];
        widthPtPythia320[bins]     = histPythia320->GetBinWidth(i);        
        errPtPythia320[bins]       = widthPtPythia320[bins] / 2.0;
        lowPtPythia320[bins]       = centerPtPythia320[bins] - errPtPythia320[bins];
        highPtPythia320[bins]      = centerPtPythia320[bins] + errPtPythia320[bins];
        
        inelPythia320[bins]        = histPythia320->GetBinContent(i) / (nEventsPythia320 * etaRange * 2 * M_PI * ptPythia320[bins]);
        errInelPythia320[bins]     = histPythia320->GetBinError(i) / (nEventsPythia320 * etaRange * 2 * M_PI * ptPythia320[bins]);
        lowErrInelPythia320[bins]  = inelPythia320[bins] - errInelPythia320[bins];
        highErrInelPythia320[bins] = inelPythia320[bins] + errInelPythia320[bins];
        relErrInelPythia320[bins]  = errInelPythia320[bins] / inelPythia320[bins];
        
        
        
        
        
        
        centerPt2PiPtPythia320[bins]    = histPythia320->GetBinCenter(i);
        pt2PiPtPythia320[bins]          = centerPt2PiPtPythia320[bins];
        widthPt2PiPtPythia320[bins]     = histPythia320->GetBinWidth(i);        
        errPt2PiPtPythia320[bins]       = widthPt2PiPtPythia320[bins] / 2.0;
        lowPt2PiPtPythia320[bins]       = centerPt2PiPtPythia320[bins] - errPt2PiPtPythia320[bins];
        highPt2PiPtPythia320[bins]      = centerPt2PiPtPythia320[bins] + errPt2PiPtPythia320[bins];        
        inel2PiPtPythia320[bins]        = histPythia320->GetBinContent(i) / (nEventsPythia320 * etaRange);
        errInel2PiPtPythia320[bins]     = histPythia320->GetBinError(i) / (nEventsPythia320 * etaRange);
        lowErrInel2PiPtPythia320[bins]  = inel2PiPtPythia320[bins] - errInel2PiPtPythia320[bins];
        highErrInel2PiPtPythia320[bins] = inel2PiPtPythia320[bins] + errInel2PiPtPythia320[bins];
        relErrInel2PiPtPythia320[bins]  = errInel2PiPtPythia320[bins] / inel2PiPtPythia320[bins];         
        
        
        
        
        
        
        cout << "ptPythia320[" << bins << "]              = " << ptPythia320[bins] <<endl;
        cout << "   centerPtPythia320[" << bins << "]     = " << centerPtPythia320[bins] <<endl;
        cout << "   widthPtPythia320[" << bins << "]      = " << widthPtPythia320[bins] <<endl;
        cout << "   errPtPythia320[" << bins << "]        = " << errPtPythia320[bins] <<endl;
        cout << "   lowPtPythia320[" << bins << "]        = " << lowPtPythia320[bins] <<endl;
        cout << "   highPtPythia320[" << bins << "]       = " << highPtPythia320[bins] <<endl;
        cout << "inelPythia320[" << bins << "]            = " << inelPythia320[bins] <<endl;        
        cout << "errInelPythia320[" << bins << "]         = " << errInelPythia320[bins] <<endl;
        cout << "   lowErrInelPythia320[" << bins << "]   = " << lowErrInelPythia320[bins] <<endl;
        cout << "   highErrInelPythia320[" << bins << "]  = " << highErrInelPythia320[bins] <<endl;
        cout << "   relErrInelPythia320[" << bins << "]   = " << relErrInelPythia320[bins] <<endl;
        cout << endl;
        bins++;
    }
}

filePythia320->Close();

cout << "================================" << endl;
cout << endl;
cout << "Finished reading PYTHIA 320 INEL data" <<endl;
cout << "Number of bins read: " << bins <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return bins;
}