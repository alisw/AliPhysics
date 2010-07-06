Int_t readPythia306() 
{
// read PYTHIA 306 INEL data from file

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read PHYTIA 306 INEL data from file" <<endl;
cout << "Number of bins: " << binsPythia306 <<endl;
cout << "Filename:       " << filenamePythia306 <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

TFile* filePythia306 = new TFile(filenamePythia306,"OPEN");
nEventsPythia306 = ((TH1F*)filePythia306->Get("eventsINEL"))->GetBinContent(1);

TH1F* histPythia306 = (TH1F*)filePythia306->Get("ptINEL");

int bins = 0;
for (int i=0; i < histPythia306->GetNbinsX(); i++) {
    if (histPythia306->GetBinContent(i) > 0 ) {
        centerPtPythia306[bins]    = histPythia306->GetBinCenter(i);
        ptPythia306[bins]          = centerPtPythia306[bins];
        widthPtPythia306[bins]     = histPythia306->GetBinWidth(i);        
        errPtPythia306[bins]       = widthPtPythia306[bins] / 2.0;
        lowPtPythia306[bins]       = centerPtPythia306[bins] - errPtPythia306[bins];
        highPtPythia306[bins]      = centerPtPythia306[bins] + errPtPythia306[bins];
        
        inelPythia306[bins]        = histPythia306->GetBinContent(i) / (nEventsPythia306 * etaRange * 2 * M_PI * ptPythia306[bins]);
        errInelPythia306[bins]     = histPythia306->GetBinError(i) / (nEventsPythia306 * etaRange * 2 * M_PI * ptPythia306[bins]);
        lowErrInelPythia306[bins]  = inelPythia306[bins] - errInelPythia306[bins];
        highErrInelPythia306[bins] = inelPythia306[bins] + errInelPythia306[bins];
        relErrInelPythia306[bins]  = errInelPythia306[bins] / inelPythia306[bins];
        
        
        
        
        
        centerPt2PiPtPythia306[bins]    = histPythia306->GetBinCenter(i);
        pt2PiPtPythia306[bins]          = centerPt2PiPtPythia306[bins];
        widthPt2PiPtPythia306[bins]     = histPythia306->GetBinWidth(i);        
        errPt2PiPtPythia306[bins]       = widthPt2PiPtPythia306[bins] / 2.0;
        lowPt2PiPtPythia306[bins]       = centerPt2PiPtPythia306[bins] - errPt2PiPtPythia306[bins];
        highPt2PiPtPythia306[bins]      = centerPt2PiPtPythia306[bins] + errPt2PiPtPythia306[bins];        
        inel2PiPtPythia306[bins]        = histPythia306->GetBinContent(i) / (nEventsPythia306 * etaRange);
        errInel2PiPtPythia306[bins]     = histPythia306->GetBinError(i) / (nEventsPythia306 * etaRange);
        lowErrInel2PiPtPythia306[bins]  = inel2PiPtPythia306[bins] - errInel2PiPtPythia306[bins];
        highErrInel2PiPtPythia306[bins] = inel2PiPtPythia306[bins] + errInel2PiPtPythia306[bins];
        relErrInel2PiPtPythia306[bins]  = errInel2PiPtPythia306[bins] / inel2PiPtPythia306[bins];         
        
        
        
        
        
        
        
        
        cout << "ptPythia306[" << bins << "]              = " << ptPythia306[bins] <<endl;
        cout << "   centerPtPythia306[" << bins << "]     = " << centerPtPythia306[bins] <<endl;
        cout << "   widthPtPythia306[" << bins << "]      = " << widthPtPythia306[bins] <<endl;
        cout << "   errPtPythia306[" << bins << "]        = " << errPtPythia306[bins] <<endl;
        cout << "   lowPtPythia306[" << bins << "]        = " << lowPtPythia306[bins] <<endl;
        cout << "   highPtPythia306[" << bins << "]       = " << highPtPythia306[bins] <<endl;
        cout << "inelPythia306[" << bins << "]            = " << inelPythia306[bins] <<endl;        
        cout << "errInelPythia306[" << bins << "]         = " << errInelPythia306[bins] <<endl;
        cout << "   lowErrInelPythia306[" << bins << "]   = " << lowErrInelPythia306[bins] <<endl;
        cout << "   highErrInelPythia306[" << bins << "]  = " << highErrInelPythia306[bins] <<endl;
        cout << "   relErrInelPythia306[" << bins << "]   = " << relErrInelPythia306[bins] <<endl;
        cout << endl;
        bins++;
    }
}

filePythia306->Close();

cout << "================================" << endl;
cout << endl;
cout << "Finished reading PYTHIA 306 INEL data" <<endl;
cout << "Number of bins read: " << bins <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return bins;
}