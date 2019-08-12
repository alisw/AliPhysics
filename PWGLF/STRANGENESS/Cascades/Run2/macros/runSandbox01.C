////////////////////////////////////////////////////////////////////////
//
// Code to replay V0 finding sequence for hypertritons (and others)
//
////////////////////////////////////////////////////////////////////////

int runSandbox01(){
    // Load common libraries
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libOADB");
    gSystem->Load("libPWGLFSTRANGENESS");
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    cout<<"Compiling sandbox..."<<endl;
    Int_t workedornot = gSystem->CompileMacro("sandbox-findable.C","-kfo");
    cout<<"----------------------------------------------------"<<endl;
    cout<<endl;
    if( workedornot == 0 ){
        cout<<"*************************************"<<endl;
        cout<<" sandbox class compilation failed! "<<endl;
        cout<<"*************************************"<<endl;
        return 0;
    }
    gROOT->ProcessLine("sandbox01();");
    return 0;
} 
