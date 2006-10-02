#include "AliITSdigit.h"
#include "AliITSLoader.h"
#include "AliRunLoader.h"
#include "AliLog.h"

///////////////////////////////////////////////////////////////////////////
// Loader for ITS
// it manages the I/O for:
// raw clusters, primary vertices
// V0 and cascade
// and tracks propagated to the origin
//////////////////////////////////////////////////////////////////////////
const TString AliITSLoader::fgkDefaultRawClustersContainerName = "TreeC";
const TString AliITSLoader::fgkDefaultBackTracksContainerName = "TreeB";
const TString AliITSLoader::fgkDefaultVerticesContainerName = "Vertex";
const TString AliITSLoader::fgkDefaultV0ContainerName = "V0";
const TString AliITSLoader::fgkDefaultCascadeContainerName = "Cascade";
ClassImp(AliITSLoader)

/**********************************************************************/
  AliITSLoader::AliITSLoader():AliLoader(),
fITSpid(0),
fGeom(0){
  // Default constructor
}
/*********************************************************************/
AliITSLoader::AliITSLoader(const Char_t *name,const Char_t *topfoldername):
AliLoader(name,topfoldername),
fITSpid(0),
fGeom(0){
  //Constructor   
    AliDataLoader* rawClustersDataLoader = new AliDataLoader(
        fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,
        "Raw Clusters");
    fDataLoaders->Add(rawClustersDataLoader);
    rawClustersDataLoader->SetEventFolder(fEventFolder);
    rawClustersDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* backTracksDataLoader =  new AliDataLoader(
        fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,
        "Back Propagated Tracks");
    fDataLoaders->Add(backTracksDataLoader);
    backTracksDataLoader->SetEventFolder(fEventFolder);
    backTracksDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* vertexDataLoader = new AliDataLoader(
        fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,
        "Primary Vertices","O");
    fDataLoaders->Add(vertexDataLoader);
    vertexDataLoader->SetEventFolder(fEventFolder);
    vertexDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* v0DataLoader = new AliDataLoader(
        fDetectorName + ".V0s.root",fgkDefaultV0ContainerName,"V0 Vertices");
    fDataLoaders->Add(v0DataLoader);
    v0DataLoader->SetEventFolder(fEventFolder);
    v0DataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* cascadeDataLoader = new AliDataLoader(
        fDetectorName + ".Cascades.root",fgkDefaultCascadeContainerName,
        "Cascades");
    fDataLoaders->Add(cascadeDataLoader);
    cascadeDataLoader->SetEventFolder(fEventFolder);
    cascadeDataLoader->SetFolder(GetDetectorDataFolder());
}
/**********************************************************************/
AliITSLoader::AliITSLoader(const Char_t *name,TFolder *topfolder): 
  AliLoader(name,topfolder),
fITSpid(0),
fGeom(0){
  //ctor  
    AliDataLoader*  rawClustersDataLoader = new AliDataLoader(
        fDetectorName + ".RawCl.root",fgkDefaultRawClustersContainerName,
        "Raw Clusters"); 
    fDataLoaders->Add(rawClustersDataLoader);
    rawClustersDataLoader->SetEventFolder(fEventFolder);
    rawClustersDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader*  backTracksDataLoader =  new AliDataLoader(
        fDetectorName + ".BackTracks.root",fgkDefaultBackTracksContainerName,
        "Back Propagated Tracks");
    fDataLoaders->Add(backTracksDataLoader);
    backTracksDataLoader->SetEventFolder(fEventFolder);
    backTracksDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* vertexDataLoader = new AliDataLoader(
        fDetectorName + ".Vertex.root",fgkDefaultVerticesContainerName,
        "Primary Vertices","O");
    fDataLoaders->Add(vertexDataLoader);
    vertexDataLoader->SetEventFolder(fEventFolder);
    vertexDataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* v0DataLoader = new AliDataLoader(
        fDetectorName + ".V0.root",fgkDefaultV0ContainerName,"V0 Vertices");
    fDataLoaders->Add(v0DataLoader);
    v0DataLoader->SetEventFolder(fEventFolder);
    v0DataLoader->SetFolder(GetDetectorDataFolder());

    AliDataLoader* cascadeDataLoader = new AliDataLoader(
        fDetectorName + ".Cascade.root",fgkDefaultCascadeContainerName,
        "Cascades");
    fDataLoaders->Add(cascadeDataLoader);
    cascadeDataLoader->SetEventFolder(fEventFolder);
    cascadeDataLoader->SetFolder(GetDetectorDataFolder());
}


/**********************************************************************/
AliITSLoader::~AliITSLoader(){
    //destructor
    AliDataLoader* dl = 0;
    UnloadRawClusters();
    dl = GetRawClLoader();
    fDataLoaders->Remove(dl);

    UnloadBackTracks();
    dl = GetBackTracksDataLoader();
    fDataLoaders->Remove(dl);

    UnloadVertices();
    dl = GetVertexDataLoader();
    fDataLoaders->Remove(dl);

    UnloadV0s();
    dl = GetV0DataLoader();
    fDataLoaders->Remove(dl);

    UnloadCascades();
    dl = GetCascadeDataLoader();
    fDataLoaders->Remove(dl);
  
    if(fGeom)delete fGeom;
    fGeom = 0;
}
/*
//----------------------------------------------------------------------
AliITS* AliITSLoader::GetITS(){
    // Returns the pointer to the ITS, kept on the file. A short cut metthod
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Returns:
    //    Returns a pointer to the ITS, if not found returns 0.
    AliITS *its;

    if(gAlice){
        its = dynamic_cast<AliITS*> (gAlice->GetDetector(
            GetDetectorName().Data()));
        if(its) return its;
    } // end if gAlice
    AliRunLoader *rl=0;
    rl = GetRunLoader();
    if(!rl) return 0;
    AliRun *ar=0;
    ar = rl->GetAliRun();
    if(!ar) return 0;
    its = dynamic_cast<AliITS*> (ar->GetDetector(GetDetectorName().Data()));
    return its;
}
//----------------------------------------------------------------------
void AliITSLoader::SetupDigits(AliITS *its){
    // Sets up to store ITS Digits in side AliITS::fDtype TObjArray
    // Inputs:
    //    AliITS *its  Pointer to the ITS
    // Outputs:
    //    none.
    // Return:
    //    none.

    its->SetTreeAddressD(TreeD());
}
*/
//----------------------------------------------------------------------
void AliITSLoader::SetupDigits(TObjArray *digPerDet,Int_t n,
                               const Char_t **digclass){
    // Sets up digPerDet to store ITS Digits.
    // Inputs:
    //    TObjArray *digPerDet   A pointer to a TObject Array size>=3.
    //    Int_t      n           The size of the TObjArray and digclass array
    //    Char_t     **digclass  Array of digit class names
    // Outputs:
    //    TObjArray *digPerDet   Setup and linked to the tree of digits
    // Return:
    //    none.
    Int_t i,m;
    TClonesArray *cl = 0;
    TTree *td = 0;
    TBranch *br = 0;
    Char_t branch[13];
    const Char_t *det[3] = {"SPD","SDD","SSD"};

    if(!digPerDet){
        Error("SetUpDigits","TObject Array digPerDet does not exist");
        return;
    } // end if
    m = digPerDet->GetSize();
    if(m<n){
        Error("SetUpDigits","TObject Array digPerDet=%p must have a size"
              " at least that of n=%d",digPerDet,n);
    } // end if
    if(m<3){
        Error("SetUpDigits","TObject Array digPerDet=%p must have a size >2",
              digPerDet);
        return;
    } // end if
    td = TreeD();
    for(i=0;i<n;i++){
        if(digPerDet->At(i)==0){ // set up TClones Array
            digPerDet->AddAt(new TClonesArray(digclass[i],1000),i);
            if(n==3) sprintf(branch,"ITSDigits%s",det[i]);
            else     sprintf(branch,"ITSDigits%d",i+1);
            br = td->GetBranch(branch);
            br->SetAddress(&((*digPerDet)[i]));
            continue; // do next one.
        } // end if
        cl =  dynamic_cast<TClonesArray*> (digPerDet->At(i));
        if(!cl && digPerDet->At(i)!=0){  // not a TClonesArray
            Error("SetUpDigits","TObject Array digPerDet-At(%d)=%p must be "
                  "zeroed or filled with TClonesArrays",i,digPerDet);
            return;
        } // end if
        if(!(cl->GetClass()->GetBaseClass(AliITSdigit::Class()))){
            Error("SetUPDigits","TClones array at digPerDet[%d}=%p must be"
                  "derived from AliITSdigit",i,digPerDet->At(i));
        } // end if
        cl->Clear();
        if(n==3) sprintf(branch,"ITSDigits%s",det[i]);
        else     sprintf(branch,"ITSDigits%d",i+1);
        br = td->GetBranch(branch);
        br->SetAddress(&((*digPerDet)[i]));
        continue;
    } // end for i
}
//---------------------------------------------------------------------
AliITSdigit * AliITSLoader::GetDigit(TObjArray *digPerDet,Int_t module,
                                     Int_t digit){
    // Gets the digit for for a specific detector type and module.
    // To be used in conjustion with Setupdigits(AliITS *its).
    // Inputs:
    //   TObjArray *digPereDet    Pointer to the Array of digits
    //   Int_t      module        Module number
    //   Int_t      digit         Digit number
    // Outputs:
    //   none.
    // Return:
    //   returns the pointer to the digit. if zero then last digit for that
    //   module.

    if(digPerDet==0){
        Error("GetDigit","digPerDet=%p, module=%d, digit=%d",
              digPerDet,module,digit);
        return 0;
    } // end if
    return 0;
}
/*
//---------------------------------------------------------------------
AliITSdigit * AliITSLoader::GetDigit(AliITS *its,Int_t module,Int_t digit){
    // Gets the digit for for a specific detector type and module.
    // To be used in conjustion with Setupdigits(AliITS *its).
    // Inputs:
    //   AliITS *its    Pointer to the ITS
    //   Int_t  module  Module number
    //   Int_t digit    Digit number
    // Outputs:
    //   none.
    // Return:
    //   returns the pointer to the digit. if zero then last digit for that
    //   module.
    //AliITSDetType *idtype;
    AliITSgeom *geom = its->GetITSgeom();
    Int_t idet = geom->GetModuleType(module);
    TClonesArray *digits;
 
    its->ResetDigits();
    TreeD()->GetEvent(module);
    digits = its->DigitsAddress(idet);
    if(digit>-1 && digit<digits->GetEntriesFast()){ // if in range.
        return (AliITSdigit*) digits->At(digit);
    } // end if
    return 0;
}
*/
//----------------------------------------------------------------------
void AliITSLoader::MakeTree(Option_t *opt){
    // invokes AliLoader::MakeTree + specific ITS tree(s)
    // Valid options: H,S,D,R,T and C (C=raw clusters)
    AliLoader::MakeTree(opt);
    const char *oC = strstr(opt,"C");
    if (oC) MakeRawClustersContainer();

    const char *oB = strstr(opt,"B");
    if (oB) MakeBackTracksContainer();

    const char *oV0 = strstr(opt,"V0");
    if (oV0) MakeV0Container();

    const char *oX = strstr(opt,"X");
    if (oX) MakeCascadeContainer();
}

//----------------------------------------------------------------------
AliITSgeom* AliITSLoader::GetITSgeom(Bool_t force) {
  // retrieves the ITS geometry from file
  if(fGeom && !force)return fGeom;
  if(fGeom && force){
    delete fGeom;
    fGeom = 0;
  }
  AliRunLoader *runLoader = GetRunLoader();
  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetITSgeom", "couldn't get AliRun object");
    return NULL;
  }
  
  TDirectory *curdir = gDirectory;
  runLoader->CdGAFile();
  fGeom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
  curdir->cd();
  if(!fGeom){
    Error("GetITSgeom","no ITS geometry available");
    return NULL;
  }
  AliWarning("AliITSgeom object has been fetched from galice.root file");
  return fGeom;
}
//______________________________________________________________________
void AliITSLoader::SetITSgeom(AliITSgeom *geom){
    // Replaces the AliITSgeom object read from file with the one
    // given.
    // Inputs:
    //   AliITSgeom *geom   The AliITSgeom object to replace the one
    //                      read from the file
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(fGeom==geom) return; // Same do nothing
    if(fGeom) {
	delete fGeom;
	fGeom=0;
    }// end if
    fGeom=geom;
}

