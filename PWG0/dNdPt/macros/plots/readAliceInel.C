Int_t readAliceInel() 
{
// read ALICE INEL data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read ALICE INEL data from file" <<endl;
cout << "Number of bins: " << binsInelAlice <<endl;
cout << "Filename:       " << filenameInelAlice <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileInelAlice;
fileInelAlice.open(filenameInelAlice);

Int_t i = 0;
while(!fileInelAlice.eof()) {
    if(i == binsInelAlice) break;    
    fileInelAlice >> centerPtInelAlice[i] >>  inelAlice[i] >> statInelAlice[i] >> systInelAlice[i];
    //systInelAlice[i] = 0.15*inelAlice[i]+1e-7; // sys error has to be provided, currently 10%
    
    // the width of the pt bins, currently hardwired....
    widthPtInelAlice[i] = 0.05;
    if (centerPtInelAlice[i] > 1) widthPtInelAlice[i] = 0.1;
    if (centerPtInelAlice[i] > 2) widthPtInelAlice[i] = 0.2;
    if (centerPtInelAlice[i] > 4) widthPtInelAlice[i] = 0.5;
    if (centerPtInelAlice[i] > 7) widthPtInelAlice[i] = 1.0;
    
    ptInelAlice[i]       = centerPtInelAlice[i];
    errPtInelAlice[i]    = widthPtInelAlice[i] / 2.0;
    lowPtInelAlice[i]    = centerPtInelAlice[i] - errPtInelAlice[i];
    highPtInelAlice[i]   = centerPtInelAlice[i] + errPtInelAlice[i];      
            
    lowStatInelAlice[i]  = inelAlice[i] - statInelAlice[i];
    highStatInelAlice[i] = inelAlice[i] + statInelAlice[i];
    relStatInelAlice[i]  = statInelAlice[i] / inelAlice[i];    
    lowSystInelAlice[i]  = inelAlice[i] - systInelAlice[i];
    highSystInelAlice[i] = inelAlice[i] + systInelAlice[i];
    relSystInelAlice[i]  = systInelAlice[i] / inelAlice[i];    
    errInelAlice[i]      = systInelAlice[i] + statInelAlice[i];
    lowErrInelAlice[i]   = inelAlice[i] - errInelAlice[i];
    highErrInelAlice[i]  = inelAlice[i] + errInelAlice[i];
    relErrInelAlice[i]   = errInelAlice[i] / inelAlice[i];
    err2InelAlice[i]     = sqrt(systInelAlice[i]*systInelAlice[i] + statInelAlice[i]*statInelAlice[i]);
    lowErr2InelAlice[i]  = inelAlice[i] - err2InelAlice[i];
    highErr2InelAlice[i] = inelAlice[i] + errInelAlice[i];
    relErr2InelAlice[i]  = err2InelAlice[i] / inelAlice[i];
    
    
    
    ptInel2PiPtAlice[i]       = ptInelAlice[i];
    centerPtInel2PiPtAlice[i] = centerPtInelAlice[i];
    widthPtInel2PiPtAlice[i]  = widthPtInelAlice[i];
    errPtInel2PiPtAlice[i]    = errPtInelAlice[i];
    lowPtInel2PiPtAlice[i]    = lowPtInelAlice[i];
    highPtInel2PiPtAlice[i]   = highPtInelAlice[i];

    inel2PiPtAlice[i]         = inelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    statInel2PiPtAlice[i]     = statInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    lowStatInel2PiPtAlice[i]  = lowStatInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    highStatInel2PiPtAlice[i] = highStatInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    relStatInel2PiPtAlice[i]  = relStatInelAlice[i];
    systInel2PiPtAlice[i]     = systInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    lowSystInel2PiPtAlice[i]  = lowSystInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    highSystInel2PiPtAlice[i] = highSystInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    relSystInel2PiPtAlice[i]  = relSystInelAlice[i];
    errInel2PiPtAlice[i]      = errInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    lowErrInel2PiPtAlice[i]   = lowErrInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    highErrInel2PiPtAlice[i]  = highErrInelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    relErrInel2PiPtAlice[i]   = relErrInelAlice[i];
    err2Inel2PiPtAlice[i]     = err2InelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    lowErr2Inel2PiPtAlice[i]  = lowErr2InelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    highErr2Inel2PiPtAlice[i] = highErr2InelAlice[i]*centerPtInelAlice[i]*2*M_PI;
    relErr2Inel2PiPtAlice[i]  = relErr2InelAlice[i];
    
        
    
    cout << "ptInelAlice[" << i << "]          = " << ptInelAlice[i] <<endl;
    cout << "   centerPtInelAlice[" << i << "] = " << centerPtInelAlice[i] <<endl;
    cout << "   widthPtInelAlice[" << i << "]  = " << widthPtInelAlice[i] <<endl;
    cout << "   errPtInelAlice[" << i << "]    = " << errPtInelAlice[i] <<endl;
    cout << "   lowPtInelAlice[" << i << "]    = " << lowPtInelAlice[i] <<endl;
    cout << "   highPtInelAlice[" << i << "]   = " << highPtInelAlice[i] <<endl;
    cout << "inelAlice[" << i << "]            = " << inelAlice[i] <<endl;
    cout << "   statInelAlice[" << i << "]     = " << statInelAlice[i] <<endl;
    cout << "   lowStatInelAlice[" << i << "]  = " << lowStatInelAlice[i] <<endl;
    cout << "   highStatInelAlice[" << i << "] = " << highStatInelAlice[i] <<endl;
    cout << "   relStatInelAlice[" << i << "]  = " << relStatInelAlice[i] <<endl;
    cout << "   systInelAlice[" << i << "]     = " << systInelAlice[i] <<endl;
    cout << "   lowSystInelAlice[" << i << "]  = " << lowSystInelAlice[i] <<endl;
    cout << "   highSystInelAlice[" << i << "] = " << highSystInelAlice[i] <<endl;
    cout << "   relSystInelAlice[" << i << "]  = " << relSystInelAlice[i] <<endl;
    cout << "errInelAlice[" << i << "]         = " << errInelAlice[i] <<endl;
    cout << "   lowErrInelAlice[" << i << "]   = " << lowErrInelAlice[i] <<endl;
    cout << "   highErrInelAlice[" << i << "]  = " << highErrInelAlice[i] <<endl;
    cout << "   relErrInelAlice[" << i << "]   = " << relErrInelAlice[i] <<endl;
    cout << "err2InelAlice[" << i << "]        = " << err2InelAlice[i] <<endl;
    cout << "   lowErr2InelAlice[" << i << "]  = " << lowErr2InelAlice[i] <<endl;
    cout << "   highErr2InelAlice[" << i << "] = " << highErr2InelAlice[i] <<endl;
    cout << "   relErr2InelAlice[" << i << "]  = " << relErr2InelAlice[i] <<endl;
    cout << endl;
    /*
    ptInelAliceFit[i]      = ptInelAlice[i];
    InelInelAliceFit[i]     = InelInelAlice[i]*ptInelAlice[i];
    errInelInelAliceFit[i]  = errInelInelAlice[i];
    widthPtInelAliceFit[i] = widthPtInelAlice[i];
    errPtInelAliceFit[i]  = 0; //errPtInelAlice[i];
    lowPtInelAliceFit[i]  = lowPtInelAlice[i];
    highPtInelAliceFit[i] = highPtInelAlice[i];
    */
        
   i++;
} // while(!fileInelAlice.eof())
fileInelAlice.close();
//if (fileInelAlice) { delete fileInelAlice; }
//fileInelAlice=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading ALICE INEL data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}