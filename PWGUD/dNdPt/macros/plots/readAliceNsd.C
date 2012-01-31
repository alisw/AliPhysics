Int_t readAliceNsd() 
{
// read ALICE NSD data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read ALICE NSD data from file" <<endl;
cout << "Number of bins: " << binsNsdAlice <<endl;
cout << "Filename:       " << filenameNsdAlice <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileNsdAlice;
fileNsdAlice.open(filenameNsdAlice);

Int_t i = 0;
while(!fileNsdAlice.eof()) {
    if(i == binsNsdAlice) break;
    // textfile content: pt_bin_center:NSD_yield:errNSD_yield
    fileNsdAlice >> centerPtNsdAlice[i] >>  nsdAlice[i] >> statNsdAlice[i] >> systNsdAlice[i];
    //systNsdAlice[i] = 0.15*nsdAlice[i]+1e-7; // sys error has to be provided, currently 10%
    
    // the width of the pt bins, currently hardwired....
    widthPtNsdAlice[i] = 0.05;
    if (centerPtNsdAlice[i] > 1) widthPtNsdAlice[i] = 0.1;
    if (centerPtNsdAlice[i] > 2) widthPtNsdAlice[i] = 0.2;
    if (centerPtNsdAlice[i] > 4) widthPtNsdAlice[i] = 0.5;
    if (centerPtNsdAlice[i] > 7) widthPtNsdAlice[i] = 1.0;
    
    ptNsdAlice[i]       = centerPtNsdAlice[i];
    errPtNsdAlice[i]    = widthPtNsdAlice[i] / 2.0;
    lowPtNsdAlice[i]    = centerPtNsdAlice[i] - errPtNsdAlice[i];
    highPtNsdAlice[i]   = centerPtNsdAlice[i] + errPtNsdAlice[i];      
            
    lowStatNsdAlice[i]  = nsdAlice[i] - statNsdAlice[i];
    highStatNsdAlice[i] = nsdAlice[i] + statNsdAlice[i];
    relStatNsdAlice[i]  = statNsdAlice[i] / nsdAlice[i];    
    lowSystNsdAlice[i]  = nsdAlice[i] - systNsdAlice[i];
    highSystNsdAlice[i] = nsdAlice[i] + systNsdAlice[i];
    relSystNsdAlice[i]  = systNsdAlice[i] / nsdAlice[i];    
    errNsdAlice[i]      = systNsdAlice[i] + statNsdAlice[i];
    lowErrNsdAlice[i]   = nsdAlice[i] - errNsdAlice[i];
    highErrNsdAlice[i]  = nsdAlice[i] + errNsdAlice[i];
    relErrNsdAlice[i]   = errNsdAlice[i] / nsdAlice[i];
    err2NsdAlice[i]     = sqrt(systNsdAlice[i]*systNsdAlice[i] + statNsdAlice[i]*statNsdAlice[i]);
    lowErr2NsdAlice[i]  = nsdAlice[i] - err2NsdAlice[i];
    highErr2NsdAlice[i] = nsdAlice[i] + errNsdAlice[i];
    relErr2NsdAlice[i]  = err2NsdAlice[i] / nsdAlice[i];
        
        
        
        
    ptNsd2PiPtAlice[i]       = ptNsdAlice[i];
    centerPtNsd2PiPtAlice[i] = centerPtNsdAlice[i];
    widthPtNsd2PiPtAlice[i]  = widthPtNsdAlice[i];
    errPtNsd2PiPtAlice[i]    = errPtNsdAlice[i];
    lowPtNsd2PiPtAlice[i]    = lowPtNsdAlice[i];
    highPtNsd2PiPtAlice[i]   = highPtNsdAlice[i];

    nsd2PiPtAlice[i]         = nsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    statNsd2PiPtAlice[i]     = statNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    lowStatNsd2PiPtAlice[i]  = lowStatNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    highStatNsd2PiPtAlice[i] = highStatNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    relStatNsd2PiPtAlice[i]  = relStatNsdAlice[i];
    systNsd2PiPtAlice[i]     = systNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    lowSystNsd2PiPtAlice[i]  = lowSystNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    highSystNsd2PiPtAlice[i] = highSystNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    relSystNsd2PiPtAlice[i]  = relSystNsdAlice[i];
    errNsd2PiPtAlice[i]      = errNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    lowErrNsd2PiPtAlice[i]   = lowErrNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    highErrNsd2PiPtAlice[i]  = highErrNsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    relErrNsd2PiPtAlice[i]   = relErrNsdAlice[i];
    err2Nsd2PiPtAlice[i]     = err2NsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    lowErr2Nsd2PiPtAlice[i]  = lowErr2NsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    highErr2Nsd2PiPtAlice[i] = highErr2NsdAlice[i]*centerPtNsdAlice[i]*2*M_PI;
    relErr2Nsd2PiPtAlice[i]  = relErr2NsdAlice[i];        
        
        
        
        
        
        
        
    
    cout << "ptNsdAlice[" << i << "]          = " << ptNsdAlice[i] <<endl;
    cout << "   centerPtNsdAlice[" << i << "] = " << centerPtNsdAlice[i] <<endl;
    cout << "   widthPtNsdAlice[" << i << "]  = " << widthPtNsdAlice[i] <<endl;
    cout << "   errPtNsdAlice[" << i << "]    = " << errPtNsdAlice[i] <<endl;
    cout << "   lowPtNsdAlice[" << i << "]    = " << lowPtNsdAlice[i] <<endl;
    cout << "   highPtNsdAlice[" << i << "]   = " << highPtNsdAlice[i] <<endl;
    cout << "nsdAlice[" << i << "]            = " << nsdAlice[i] <<endl;
    cout << "   statNsdAlice[" << i << "]     = " << statNsdAlice[i] <<endl;
    cout << "   lowStatNsdAlice[" << i << "]  = " << lowStatNsdAlice[i] <<endl;
    cout << "   highStatNsdAlice[" << i << "] = " << highStatNsdAlice[i] <<endl;
    cout << "   relStatNsdAlice[" << i << "]  = " << relStatNsdAlice[i] <<endl;
    cout << "   systNsdAlice[" << i << "]     = " << systNsdAlice[i] <<endl;
    cout << "   lowSystNsdAlice[" << i << "]  = " << lowSystNsdAlice[i] <<endl;
    cout << "   highSystNsdAlice[" << i << "] = " << highSystNsdAlice[i] <<endl;
    cout << "   relSystNsdAlice[" << i << "]  = " << relSystNsdAlice[i] <<endl;
    cout << "errNsdAlice[" << i << "]         = " << errNsdAlice[i] <<endl;
    cout << "   lowErrNsdAlice[" << i << "]   = " << lowErrNsdAlice[i] <<endl;
    cout << "   highErrNsdAlice[" << i << "]  = " << highErrNsdAlice[i] <<endl;
    cout << "   relErrNsdAlice[" << i << "]   = " << relErrNsdAlice[i] <<endl;
    cout << "err2NsdAlice[" << i << "]        = " << err2NsdAlice[i] <<endl;
    cout << "   lowErr2NsdAlice[" << i << "]  = " << lowErr2NsdAlice[i] <<endl;
    cout << "   highErr2NsdAlice[" << i << "] = " << highErr2NsdAlice[i] <<endl;
    cout << "   relErr2NsdAlice[" << i << "]  = " << relErr2NsdAlice[i] <<endl;
    cout << endl;
    /*
    ptNsdAliceFit[i]      = ptNsdAlice[i];
    NsdNsdAliceFit[i]     = NsdNsdAlice[i]*ptNsdAlice[i];
    errNsdNsdAliceFit[i]  = errNsdNsdAlice[i];
    widthPtNsdAliceFit[i] = widthPtNsdAlice[i];
    errPtNsdAliceFit[i]  = 0; //errPtNsdAlice[i];
    lowPtNsdAliceFit[i]  = lowPtNsdAlice[i];
    highPtNsdAliceFit[i] = highPtNsdAlice[i];
    */
        
   i++;
} // while(!fileNsdAlice.eof())
fileNsdAlice.close();
//if (fileNsdAlice) { delete fileNsdAlice; }
//fileNsdAlice=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading ALICE NSD data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}