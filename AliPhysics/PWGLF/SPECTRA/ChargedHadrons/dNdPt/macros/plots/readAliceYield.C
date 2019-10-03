Int_t readAliceYield() 
{
// read ALICE INVARIANT YIELD data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read ALICE INVARIANT YIELD data from file" <<endl;
cout << "Number of bins: " << binsYieldAlice <<endl;
cout << "Filename:       " << filenameYieldAlice <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileYieldAlice;
fileYieldAlice.open(filenameYieldAlice);

Int_t i = 0;
while(!fileYieldAlice.eof()) {
    if(i == binsYieldAlice) break;    
    fileYieldAlice >> centerPtYieldAlice[i] >>  yieldAlice[i] >> statYieldAlice[i] >> systYieldAlice[i];
    //systYieldAlice[i] = 0.15*yieldAlice[i]+1e-7; // sys error has to be provided, currently 10%
    
    // the width of the pt bins, currently hardwired....
    widthPtYieldAlice[i] = 0.05;
    if (centerPtYieldAlice[i] > 1) widthPtYieldAlice[i] = 0.1;
    if (centerPtYieldAlice[i] > 2) widthPtYieldAlice[i] = 0.2;
    if (centerPtYieldAlice[i] > 4) widthPtYieldAlice[i] = 0.5;
    if (centerPtYieldAlice[i] > 7) widthPtYieldAlice[i] = 1.0;
    
    ptYieldAlice[i]       = centerPtYieldAlice[i];
    errPtYieldAlice[i]    = widthPtYieldAlice[i] / 2.0;
    lowPtYieldAlice[i]    = centerPtYieldAlice[i] - errPtYieldAlice[i];
    highPtYieldAlice[i]   = centerPtYieldAlice[i] + errPtYieldAlice[i];      
            
    lowStatYieldAlice[i]  = yieldAlice[i] - statYieldAlice[i];
    highStatYieldAlice[i] = yieldAlice[i] + statYieldAlice[i];
    relStatYieldAlice[i]  = statYieldAlice[i] / yieldAlice[i];    
    lowSystYieldAlice[i]  = yieldAlice[i] - systYieldAlice[i];
    highSystYieldAlice[i] = yieldAlice[i] + systYieldAlice[i];
    relSystYieldAlice[i]  = systYieldAlice[i] / yieldAlice[i];    
    errYieldAlice[i]      = systYieldAlice[i] + statYieldAlice[i];
    lowErrYieldAlice[i]   = yieldAlice[i] - errYieldAlice[i];
    highErrYieldAlice[i]  = yieldAlice[i] + errYieldAlice[i];
    relErrYieldAlice[i]   = errYieldAlice[i] / yieldAlice[i];
    err2YieldAlice[i]     = sqrt(systYieldAlice[i]*systYieldAlice[i] + statYieldAlice[i]*statYieldAlice[i]);
    lowErr2YieldAlice[i]  = yieldAlice[i] - err2YieldAlice[i];
    highErr2YieldAlice[i] = yieldAlice[i] + errYieldAlice[i];
    relErr2YieldAlice[i]  = err2YieldAlice[i] / yieldAlice[i];
        
        
        
        
    ptYield2PiPtAlice[i]       = ptYieldAlice[i];
    centerPtYield2PiPtAlice[i] = centerPtYieldAlice[i];
    widthPtYield2PiPtAlice[i]  = widthPtYieldAlice[i];
    errPtYield2PiPtAlice[i]    = errPtYieldAlice[i];
    lowPtYield2PiPtAlice[i]    = lowPtYieldAlice[i];
    highPtYield2PiPtAlice[i]   = highPtYieldAlice[i];

    yield2PiPtAlice[i]         = yieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    statYield2PiPtAlice[i]     = statYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    lowStatYield2PiPtAlice[i]  = lowStatYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    highStatYield2PiPtAlice[i] = highStatYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    relStatYield2PiPtAlice[i]  = relStatYieldAlice[i];
    systYield2PiPtAlice[i]     = systYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    lowSystYield2PiPtAlice[i]  = lowSystYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    highSystYield2PiPtAlice[i] = highSystYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    relSystYield2PiPtAlice[i]  = relSystYieldAlice[i];
    errYield2PiPtAlice[i]      = errYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    lowErrYield2PiPtAlice[i]   = lowErrYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    highErrYield2PiPtAlice[i]  = highErrYieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    relErrYield2PiPtAlice[i]   = relErrYieldAlice[i];
    err2Yield2PiPtAlice[i]     = err2YieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    lowErr2Yield2PiPtAlice[i]  = lowErr2YieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    highErr2Yield2PiPtAlice[i] = highErr2YieldAlice[i]*centerPtYieldAlice[i]*2*M_PI;
    relErr2Yield2PiPtAlice[i]  = relErr2YieldAlice[i];                
        
        
        
        
        
    
    cout << "ptYieldAlice[" << i << "]          = " << ptYieldAlice[i] <<endl;
    cout << "   centerPtYieldAlice[" << i << "] = " << centerPtYieldAlice[i] <<endl;
    cout << "   widthPtYieldAlice[" << i << "]  = " << widthPtYieldAlice[i] <<endl;
    cout << "   errPtYieldAlice[" << i << "]    = " << errPtYieldAlice[i] <<endl;
    cout << "   lowPtYieldAlice[" << i << "]    = " << lowPtYieldAlice[i] <<endl;
    cout << "   highPtYieldAlice[" << i << "]   = " << highPtYieldAlice[i] <<endl;
    cout << "yieldAlice[" << i << "]            = " << yieldAlice[i] <<endl;
    cout << "   statYieldAlice[" << i << "]     = " << statYieldAlice[i] <<endl;
    cout << "   lowStatYieldAlice[" << i << "]  = " << lowStatYieldAlice[i] <<endl;
    cout << "   highStatYieldAlice[" << i << "] = " << highStatYieldAlice[i] <<endl;
    cout << "   relStatYieldAlice[" << i << "]  = " << relStatYieldAlice[i] <<endl;
    cout << "   systYieldAlice[" << i << "]     = " << systYieldAlice[i] <<endl;
    cout << "   lowSystYieldAlice[" << i << "]  = " << lowSystYieldAlice[i] <<endl;
    cout << "   highSystYieldAlice[" << i << "] = " << highSystYieldAlice[i] <<endl;
    cout << "   relSystYieldAlice[" << i << "]  = " << relSystYieldAlice[i] <<endl;
    cout << "errYieldAlice[" << i << "]         = " << errYieldAlice[i] <<endl;
    cout << "   lowErrYieldAlice[" << i << "]   = " << lowErrYieldAlice[i] <<endl;
    cout << "   highErrYieldAlice[" << i << "]  = " << highErrYieldAlice[i] <<endl;
    cout << "   relErrYieldAlice[" << i << "]   = " << relErrYieldAlice[i] <<endl;
    cout << "err2YieldAlice[" << i << "]        = " << err2YieldAlice[i] <<endl;
    cout << "   lowErr2YieldAlice[" << i << "]  = " << lowErr2YieldAlice[i] <<endl;
    cout << "   highErr2YieldAlice[" << i << "] = " << highErr2YieldAlice[i] <<endl;
    cout << "   relErr2YieldAlice[" << i << "]  = " << relErr2YieldAlice[i] <<endl;
    cout << endl;
    /*
    ptYieldAliceFit[i]      = ptYieldAlice[i];
    YieldYieldAliceFit[i]     = YieldYieldAlice[i]*ptYieldAlice[i];
    errYieldYieldAliceFit[i]  = errYieldYieldAlice[i];
    widthPtYieldAliceFit[i] = widthPtYieldAlice[i];
    errPtYieldAliceFit[i]  = 0; //errPtYieldAlice[i];
    lowPtYieldAliceFit[i]  = lowPtYieldAlice[i];
    highPtYieldAliceFit[i] = highPtYieldAlice[i];
    */
        
   i++;
} // while(!fileYieldAlice.eof())
fileYieldAlice.close();
//if (fileYieldAlice) { delete fileYieldAlice; }
//fileYieldAlice=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading ALICE INVARIANT YIELD data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}