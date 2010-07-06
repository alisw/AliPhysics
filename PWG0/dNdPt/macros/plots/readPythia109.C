Int_t readPythia109() 
{
// read PYTHIA 109 INEL data from file

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read PHYTIA 109 INEL data from file" <<endl;
cout << "Number of bins: " << binsPythia109 <<endl;
cout << "Filename:       " << filenamePythia109 <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

TFile* filePythia109 = new TFile(filenamePythia109,"OPEN");
nEventsPythia109 = ((TH1F*)filePythia109->Get("eventsINEL"))->GetBinContent(1);

TH1F* histPythia109 = (TH1F*)filePythia109->Get("ptINEL");

int bins = 0;
for (int i=0; i < histPythia109->GetNbinsX(); i++) {
    if (histPythia109->GetBinContent(i) > 0 ) {
        centerPtPythia109[bins]    = histPythia109->GetBinCenter(i);
        ptPythia109[bins]          = centerPtPythia109[bins];
        widthPtPythia109[bins]     = histPythia109->GetBinWidth(i);        
        errPtPythia109[bins]       = widthPtPythia109[bins] / 2.0;
        lowPtPythia109[bins]       = centerPtPythia109[bins] - errPtPythia109[bins];
        highPtPythia109[bins]      = centerPtPythia109[bins] + errPtPythia109[bins];
        
        inelPythia109[bins]        = histPythia109->GetBinContent(i) / (nEventsPythia109 * etaRange * 2 * M_PI * ptPythia109[bins]);
        errInelPythia109[bins]     = histPythia109->GetBinError(i) / (nEventsPythia109 * etaRange * 2 * M_PI * ptPythia109[bins]);
        lowErrInelPythia109[bins]  = inelPythia109[bins] - errInelPythia109[bins];
        highErrInelPythia109[bins] = inelPythia109[bins] + errInelPythia109[bins];
        relErrInelPythia109[bins]  = errInelPythia109[bins] / inelPythia109[bins];
        
        
        
        centerPt2PiPtPythia109[bins]    = histPythia109->GetBinCenter(i);
        pt2PiPtPythia109[bins]          = centerPt2PiPtPythia109[bins];
        widthPt2PiPtPythia109[bins]     = histPythia109->GetBinWidth(i);        
        errPt2PiPtPythia109[bins]       = widthPt2PiPtPythia109[bins] / 2.0;
        lowPt2PiPtPythia109[bins]       = centerPt2PiPtPythia109[bins] - errPt2PiPtPythia109[bins];
        highPt2PiPtPythia109[bins]      = centerPt2PiPtPythia109[bins] + errPt2PiPtPythia109[bins];        
        inel2PiPtPythia109[bins]        = histPythia109->GetBinContent(i) / (nEventsPythia109 * etaRange);
        errInel2PiPtPythia109[bins]     = histPythia109->GetBinError(i) / (nEventsPythia109 * etaRange);
        lowErrInel2PiPtPythia109[bins]  = inel2PiPtPythia109[bins] - errInel2PiPtPythia109[bins];
        highErrInel2PiPtPythia109[bins] = inel2PiPtPythia109[bins] + errInel2PiPtPythia109[bins];
        relErrInel2PiPtPythia109[bins]  = errInel2PiPtPythia109[bins] / inel2PiPtPythia109[bins];         
        
        
        
        cout << "ptPythia109[" << bins << "]              = " << ptPythia109[bins] <<endl;
        cout << "   centerPtPythia109[" << bins << "]     = " << centerPtPythia109[bins] <<endl;
        cout << "   widthPtPythia109[" << bins << "]      = " << widthPtPythia109[bins] <<endl;
        cout << "   errPtPythia109[" << bins << "]        = " << errPtPythia109[bins] <<endl;
        cout << "   lowPtPythia109[" << bins << "]        = " << lowPtPythia109[bins] <<endl;
        cout << "   highPtPythia109[" << bins << "]       = " << highPtPythia109[bins] <<endl;
        cout << "inelPythia109[" << bins << "]            = " << inelPythia109[bins] <<endl;        
        cout << "errInelPythia109[" << bins << "]         = " << errInelPythia109[bins] <<endl;
        cout << "   lowErrInelPythia109[" << bins << "]   = " << lowErrInelPythia109[bins] <<endl;
        cout << "   highErrInelPythia109[" << bins << "]  = " << highErrInelPythia109[bins] <<endl;
        cout << "   relErrInelPythia109[" << bins << "]   = " << relErrInelPythia109[bins] <<endl;
        cout << endl;
        bins++;
    }
}

filePythia109->Close();

cout << "================================" << endl;
cout << endl;
cout << "Finished reading PYTHIA 109 INEL data" <<endl;
cout << "Number of bins read: " << bins <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return bins;
}