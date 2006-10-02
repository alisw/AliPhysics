#ifndef ALIITSLOADER_H
#define ALIITSLOADER_H
//////////////////////////////////////////////////////////
// Loader class for ITS                                 //
//////////////////////////////////////////////////////////
#include <AliLoader.h>
#include <AliESDVertex.h>
#include <AliITSgeom.h>
class AliITSpidESD;
class AliITSdigit;
class TObjArray;

class AliITSLoader: public AliLoader{
  public:
    AliITSLoader();
    AliITSLoader(const Char_t *name,const Char_t *topfoldername);
    AliITSLoader(const Char_t *name,TFolder *topfolder);

    virtual ~AliITSLoader();

    void           MakeTree(Option_t* opt);
    virtual void   SetupDigits(TObjArray *digPerDet,Int_t n,
                                const Char_t **digclass); // Sets up digits
    // Gets the AliITSdigit for a given module and a specific digit in that
    // module. Array of digits stored in AliITS (must use 
    // SetupDigits(AliITS *its)).
    // virtual AliITSdigit* GetDigit(AliITS *its,Int_t module,Int_t digit);
    // Gets the AliITSdigit for a given module and a specific digit in that
    // module. Array of digits stored in a user defined TObjArray digPerDet
    virtual AliITSdigit* GetDigit(TObjArray *digPerDet,Int_t module,Int_t digit);

    //Raw Clusters
    AliDataLoader* GetRawClLoader() {return GetDataLoader("Raw Clusters");}
    virtual void   CleanRawClusters() {
        GetRawClLoader()->GetBaseLoader(0)->Clean();}
    Int_t          LoadRawClusters(Option_t* opt=""){
        return GetRawClLoader()->GetBaseLoader(0)->Load(opt);}
    void           SetRawClustersFileName(const TString& fname){
        GetRawClLoader()->SetFileName(fname);}
    // returns a pointer to the tree of  RawClusters
    TTree*         TreeC(){ return GetRawClLoader()->Tree();} 
    void           UnloadRawClusters(){
        GetRawClLoader()->GetBaseLoader(0)->Unload();}
    virtual Int_t  WriteRawClusters(Option_t* opt=""){
        return GetRawClLoader()->GetBaseLoader(0)->WriteData(opt);}

    //Vertices
    AliDataLoader* GetVertexDataLoader() {
        return GetDataLoader("Primary Vertices");}
    virtual void   CleanVertices() {
        GetVertexDataLoader()->GetBaseLoader(0)->Clean();}
    Int_t          LoadVertices(Option_t* opt=""){
        return GetVertexDataLoader()->GetBaseLoader(0)->Load(opt);}
    void           SetVerticesFileName(const TString& fname){
        GetVertexDataLoader()->SetFileName(fname);}
    void           UnloadVertices(){
        GetVertexDataLoader()->GetBaseLoader(0)->Unload();}
    virtual Int_t  WriteVertices(Option_t* opt=""){
        return GetVertexDataLoader()->GetBaseLoader(0)->WriteData(opt);}
    virtual Int_t PostVertex(AliESDVertex *ptr){
        return GetVertexDataLoader()->GetBaseLoader(0)->Post(ptr);}
    //    virtual void SetVerticesContName(const char *name){
    //       GetVertexDataLoader()->GetBaseLoader(0)->SetName(name);}
    AliESDVertex *GetVertex(){
        return static_cast <AliESDVertex*>(GetVertexDataLoader()->
                                           GetBaseLoader(0)->Get());}
    //V0s
    AliDataLoader* GetV0DataLoader() {return GetDataLoader("V0 Vertices");}
    virtual void   CleanV0s() {GetV0DataLoader()->GetBaseLoader(0)->Clean();}
    Int_t          LoadV0s(Option_t* opt=""){
        return GetV0DataLoader()->GetBaseLoader(0)->Load(opt);}
    void           SetV0FileName(const TString& fname){
        GetV0DataLoader()->SetFileName(fname);}
    void           UnloadV0s(){GetV0DataLoader()->GetBaseLoader(0)->Unload();}
    virtual Int_t  WriteV0s(Option_t* opt=""){
        return GetV0DataLoader()->GetBaseLoader(0)->WriteData(opt);}
    TTree*         TreeV0(){ return GetV0DataLoader()->Tree();}

    //Cascades
    AliDataLoader* GetCascadeDataLoader() {return GetDataLoader("Cascades");}
    virtual void   CleanCascades() {
        GetCascadeDataLoader()->GetBaseLoader(0)->Clean();}
    Int_t          LoadCascades(Option_t* opt=""){
        return GetCascadeDataLoader()->GetBaseLoader(0)->Load(opt);}
    void           SetCascadeFileName(const TString& fname){
        GetCascadeDataLoader()->SetFileName(fname);}
    void           UnloadCascades(){
        GetCascadeDataLoader()->GetBaseLoader(0)->Unload();}
    virtual Int_t  WriteCascades(Option_t* opt=""){
        return GetCascadeDataLoader()->GetBaseLoader(0)->WriteData(opt);}
    TTree*         TreeX(){ return GetCascadeDataLoader()->Tree();}

    //Back Propagated Tracks
    AliDataLoader* GetBackTracksDataLoader() {
        return GetDataLoader("Back Propagated Tracks");}
    virtual void   CleanBackTracks() {
        GetBackTracksDataLoader()->GetBaseLoader(0)->Clean();}
    Int_t          LoadBackTracks(Option_t* opt=""){
        return GetBackTracksDataLoader()->GetBaseLoader(0)->Load(opt);}
    void           SetBackTracksFileName(const TString& fname){
        GetBackTracksDataLoader()->SetFileName(fname);}
     // returns a pointer to the tree of  BackTracks
    TTree*         TreeB(){ return GetBackTracksDataLoader()->Tree();}
    void           UnloadBackTracks(){
        GetBackTracksDataLoader()->GetBaseLoader(0)->Unload();}
    virtual Int_t  WriteBackTracks(Option_t* opt=""){
        return GetBackTracksDataLoader()->GetBaseLoader(0)->WriteData(opt);}

    // Geometry. Geom is read from file, unless already loaded
    // readout from file can be forced if force=kTRUE
    AliITSgeom* GetITSgeom(Bool_t force=kFALSE); 
    void SetITSgeom(AliITSgeom* g);
    // PID
    AliITSpidESD* GetITSpid() const {return fITSpid;}
    void  AdoptITSpid(AliITSpidESD* pid) {fITSpid=pid;}
  protected:

    AliITSLoader(const AliITSLoader &ob); // copy constructor
    AliITSLoader& operator=(const AliITSLoader & /* source */); // ass.

    // METHODS
    virtual void   MakeRawClustersContainer() {GetRawClLoader()->MakeTree();}
    Int_t          PostRawClusters(){
        return GetRawClLoader()->GetBaseLoader(0)->Post();}

    virtual void   MakeBackTracksContainer() {
        GetBackTracksDataLoader()->MakeTree();}
    Int_t          PostBackTracks(){
        return GetBackTracksDataLoader()->GetBaseLoader(0)->Post();}
    virtual void   MakeV0Container() {GetV0DataLoader()->MakeTree();}
    Int_t          PostV0s(){
        return GetV0DataLoader()->GetBaseLoader(0)->Post();}

    virtual void   MakeCascadeContainer() {GetCascadeDataLoader()->MakeTree();}
    Int_t          PostCascades(){
        return GetCascadeDataLoader()->GetBaseLoader(0)->Post();}

    // DATA
    static const TString fgkDefaultRawClustersContainerName;  //default for Raw Clusters container name
    static const TString fgkDefaultBackTracksContainerName;   //default for Back propag. tracks container name
    static const TString fgkDefaultVerticesContainerName;     //default for primary vertices container name
    static const TString fgkDefaultV0ContainerName;           //default for V0 container name
    static const TString fgkDefaultCascadeContainerName;      //default fo cascade container name
    AliITSpidESD* fITSpid; //! pointer for ITS pid
    AliITSgeom *fGeom;     //! pointer to the ITS geometry class


    ClassDef(AliITSLoader,5) // Loader for additional ITS specific trees.
};
 
#endif
