#include "TEveGeoShapeExtract.h"
#include "TBuffer3D.h"
#include "TColor.h"
#include <TXMLEngine.h>

class AliColladaBuffer
{
public:
    AliColladaBuffer();
    ~AliColladaBuffer();
    
    // method to add non-standard lambert material
    void AddNewLambertMaterial(const char* name, TColor *emission, TColor *ambient, TColor *diffuse, double opacity = 1.0);
    
    // add TEveGeoShapeExtract to Collada buffer
    void AddShape(const char *name,TEveGeoShapeExtract *shape,const char* materialName="default_lambert");
    void AddShape(const char* name,TGeoShape *shape, double translation[3], double rotation[3], double scale[3],const char* materialName="default_lambert");
    
    // add TBuffer3D to Collada buffer
    void AddShape(const char* name,TBuffer3D *buffer, double trans[16],const char* materialName="default_lambert");
    void AddShape(const char* name,TBuffer3D *buffer, double translation[3],double rotation[3], double scale[3],const char* materialName="default_lambert");
    
    // wrappers to add primitives to Collada buffer
    void AddBox(const char* name,double width, double height, double depth,double translation[3], double rotation[3], double scale[3],const char* materialName="default_lambert");
    
    // finally, save buffer to file
    void SaveAs(const char* filename);
    
private:
    TXMLEngine* fXml;             // XML engine
    XMLNodePointer_t fRootNode;   // root XML node
    
    // main COLLADA nodes
    XMLNodePointer_t fAsset;
    XMLNodePointer_t fLibraryMaterials;
    XMLNodePointer_t fLibraryEffects;
    XMLNodePointer_t fLibraryGeometries;
    XMLNodePointer_t fLibraryVisualScenes;
    XMLNodePointer_t fScene;
    
    // scenes
    XMLNodePointer_t fVisualScene;
    
    void SetupRootNode();
    void SetupMainNodes();
    void SetupAssetNode();
    void SetupScene();
    
    void GetTrianglesNormal(double p0[3],double p1[3], double p2[3], double *result);
    void GetTransMatrixFromParams(double *trans, double *translation, double *rotation, double *scale);
    
    void AddShapeToCollada(const char* name,
                           const char* verticesString,int verticesSize,
                           const char* normalString,
                           const char* trianglesString,int trianglesSize,
                           const char* transformationString,
                           const char* materialName);
};
