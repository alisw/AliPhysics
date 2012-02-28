Int_t readUa1() 
{
// read UA1 data from file
//

cout << endl;
cout << "================================" << endl;
cout << endl;
cout << "read UA1 YIELD data from file" <<endl;
cout << "Number of bins: " << binsUa1 <<endl;
cout << "Filename:       " << filenameUa1 <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

ifstream fileCrossUa1;
fileCrossUa1.open(filenameUa1);

Int_t i = 0;
while(!fileCrossUa1.eof()) {
    if(i == binsUa1) break;
    // textfile content: pt:cross_yield:(stat + syst error added linerarly)
    fileCrossUa1 >> centerPtUa1[i] >> crossUa1[i] >> errCrossUa1[i];
    
    ptUa1[i]          = centerPtUa1[i];
    
    // the width of the pt bins, currently hardwired....
    widthPtUa1[i] = 0.1;
    if (centerPtUa1[i] > 4) widthPtUa1[i] = 0.2;
    if (centerPtUa1[i] > 6) widthPtUa1[i] = 1.0;
    
    errPtUa1[i]       = widthPtUa1[i] / 2.0;              
    lowPtUa1[i]       = centerPtUa1[i] - errPtUa1[i];
    highPtUa1[i]      = centerPtUa1[i] + errPtUa1[i];
    
    /*
    lowStatCrossUa1[i]  = crossUa1[i] - statCrossUa1[i];
    highStatCrossUa1[i] = crossUa1[i] + statCrossUa1[i];
    relStatCrossUa1[i]  = statCrossUa1[i] / crossUa1[i];    
    lowSystCrossUa1[i]  = crossUa1[i] - systCrossUa1[i];
    highSystCrossUa1[i] = crossUa1[i] + systCrossUa1[i];
    relSystCrossUa1[i]  = systCrossUa1[i] / crossUa1[i];    
    */    
    lowErrCrossUa1[i]   = crossUa1[i] - errCrossUa1[i];
    highErrCrossUa1[i]  = crossUa1[i] + errCrossUa1[i];
    relErrCrossUa1[i]   = errCrossUa1[i] / crossUa1[i];
    err2CrossUa1[i]     = errCrossUa1[i];
    lowErr2CrossUa1[i]  = lowErrCrossUa1[i];
    highErr2CrossUa1[i] = highErrCrossUa1[i];
    relErr2CrossUa1[i]  = relErrCrossUa1[i];
    
    yieldUa1[i]         = crossUa1[i] * (avgToHadr / sigmaInelUa1);
    /*
    statYieldUa1[i]     = 
    lowStatYieldUa1[i]  = 
    highStatYieldUa1[i] = 
    relStatYieldUa1[i]  = 
    systYieldUa1[i]     = 
    lowSystYieldUa1[i]  = 
    highSystYieldUa1[i] = 
    relSystYieldUa1[i]  = 
    */
    errYieldUa1[i]      = errCrossUa1[i] * (avgToHadr / sigmaInelUa1);
    lowErrYieldUa1[i]   = yieldUa1[i] - errYieldUa1[i];
    highErrYieldUa1[i]  = yieldUa1[i] + errYieldUa1[i];
    relErrYieldUa1[i]   = errYieldUa1[i] / yieldUa1[i];
    err2YieldUa1[i]     = errYieldUa1[i];
    lowErr2YieldUa1[i]  = lowErrYieldUa1[i];
    highErr2YieldUa1[i] = highErrYieldUa1[i];
    relErr2YieldUa1[i]  = relErrYieldUa1[i];
        
    cout << "ptUa1[" << i << "]               = " << ptUa1[i] <<endl;
    cout << "   centerPtUa1[" << i << "]      = " << centerPtUa1[i] <<endl;
    cout << "   widthPtUa1[" << i << "]       = " << widthPtUa1[i] <<endl;
    cout << "   errPtUa1[" << i << "]         = " << errPtUa1[i] <<endl;
    cout << "   lowPtUa1[" << i << "]         = " << lowPtUa1[i] <<endl;
    cout << "   highPtUa1[" << i << "]        = " << highPtUa1[i] <<endl;
    cout << "crossUa1[" << i << "]            = " << crossUa1[i] <<endl;
    /*
    cout << "   statCrossUa1[" << i << "]     = " << statCrossUa1[i] <<endl;
    cout << "   lowStatCrossUa1[" << i << "]  = " << lowStatCrossUa1[i] <<endl;
    cout << "   highStatCrossUa1[" << i << "] = " << highStatCrossUa1[i] <<endl;
    cout << "   relStatCrossUa1[" << i << "]  = " << relStatCrossUa1[i] <<endl;
    cout << "   systCrossUa1[" << i << "]     = " << systCrossUa1[i] <<endl;
    cout << "   lowSystCrossUa1[" << i << "]  = " << lowSystCrossUa1[i] <<endl;
    cout << "   highSystCrossUa1[" << i << "] = " << highSystCrossUa1[i] <<endl;
    cout << "   relSystCrossUa1[" << i << "]  = " << relSystCrossUa1[i] <<endl;
    */
    cout << "errCrossUa1[" << i << "]         = " << errCrossUa1[i] <<endl;
    cout << "   lowErrCrossUa1[" << i << "]   = " << lowErrCrossUa1[i] <<endl;
    cout << "   highErrCrossUa1[" << i << "]  = " << highErrCrossUa1[i] <<endl;
    cout << "   relErrCrossUa1[" << i << "]   = " << relErrCrossUa1[i] <<endl;
    cout << "err2CrossUa1[" << i << "]        = " << err2CrossUa1[i] <<endl;
    cout << "   lowErr2CrossUa1[" << i << "]  = " << lowErr2CrossUa1[i] <<endl;
    cout << "   highErr2CrossUa1[" << i << "] = " << highErr2CrossUa1[i] <<endl;
    cout << "   relErr2CrossUa1[" << i << "]  = " << relErr2CrossUa1[i] <<endl;

    cout << "yieldUa1[" << i << "]            = " << yieldUa1[i] <<endl;
    /*
    cout << "   statYieldUa1[" << i << "]     = " << statYieldUa1[i] <<endl;
    cout << "   lowStatYieldUa1[" << i << "]  = " << lowStatYieldUa1[i] <<endl;
    cout << "   highStatYieldUa1[" << i << "] = " << highStatYieldUa1[i] <<endl;
    cout << "   relStatYieldUa1[" << i << "]  = " << relStatYieldUa1[i] <<endl;
    cout << "   systYieldUa1[" << i << "]     = " << systYieldUa1[i] <<endl;
    cout << "   lowSystYieldUa1[" << i << "]  = " << lowSystYieldUa1[i] <<endl;
    cout << "   highSystYieldUa1[" << i << "] = " << highSystYieldUa1[i] <<endl;
    cout << "   relSystYieldUa1[" << i << "]  = " << relSystYieldUa1[i] <<endl;
    */
    cout << "errYieldUa1[" << i << "]         = " << errYieldUa1[i] <<endl;
    cout << "   lowErrYieldUa1[" << i << "]   = " << lowErrYieldUa1[i] <<endl;
    cout << "   highErrYieldUa1[" << i << "]  = " << highErrYieldUa1[i] <<endl;
    cout << "   relErrYieldUa1[" << i << "]   = " << relErrYieldUa1[i] <<endl;
    cout << "err2YieldUa1[" << i << "]        = " << err2YieldUa1[i] <<endl;
    cout << "   lowErr2YieldUa1[" << i << "]  = " << lowErr2YieldUa1[i] <<endl;
    cout << "   highErr2YieldUa1[" << i << "] = " << highErr2YieldUa1[i] <<endl;
    cout << "   relErr2YieldUa1[" << i << "]  = " << relErr2YieldUa1[i] <<endl;
    cout << endl;
     
   i++;
} // while(!fileCrossUa1.eof())
fileCrossUa1.close();
//if (fileCrossUa1) { delete fileCrossUa1; }
//fileCrossUa1=0;

cout << "================================" << endl;
cout << endl;
cout << "Finished reading UA1 YIELD data" <<endl;
cout << "Number of bins read: " << i <<endl;
cout << endl;
cout << "================================" << endl;
cout << endl;

return i;
}