#include "AliColladaBuffer.h"

#include "TMath.h"
#include "TEveGeoShape.h"
#include "TEveGeoPolyShape.h"
#include "TGeoBBox.h"
#include "TGeoSphere.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

AliColladaBuffer::AliColladaBuffer()
{
    fXml = new TXMLEngine();
    SetupRootNode();
    SetupMainNodes();
    SetupAssetNode();
    SetupScene();
    
    TColor *emission = new TColor();
    TColor *ambient = new TColor();
    TColor *diffuse = new TColor();
    
    emission->SetRGB(0.0,0.0,0.0);
    emission->SetAlpha(1.0);
    ambient->SetRGB(0.0,0.0,0.0);
    ambient->SetAlpha(1.0);
    diffuse->SetRGB(1.0,0.4,0.4);
    diffuse->SetAlpha(1.0);
    
    // add one default material
    AddNewLambertMaterial("default_lambert",emission,ambient,diffuse,1.0);
    
    delete emission;
    delete ambient;
    delete diffuse;
}

AliColladaBuffer::~AliColladaBuffer()
{
    delete fXml;
}

void AliColladaBuffer::SetupRootNode()
{
    fRootNode = fXml->NewChild(0, 0, "COLLADA");
    fXml->NewAttr(fRootNode,0,"xmlns","http://www.collada.org/2005/11/COLLADASchema");
    fXml->NewAttr(fRootNode,0,"version","1.4.1");
}

void AliColladaBuffer::SetupMainNodes()
{
    fAsset = fXml->NewChild(fRootNode, 0, "asset");
    fLibraryMaterials = fXml->NewChild(fRootNode, 0, "library_materials");
    fLibraryEffects = fXml->NewChild(fRootNode, 0, "library_effects");
    fLibraryGeometries = fXml->NewChild(fRootNode, 0, "library_geometries");
    fLibraryVisualScenes = fXml->NewChild(fRootNode,0,"library_visual_scenes");
    fScene = fXml->NewChild(fRootNode,0,"scene");
}

void AliColladaBuffer::SetupAssetNode()
{
    fXml->NewChild(fAsset, 0, "contributor");
    fXml->NewChild(fAsset, 0, "created","2015-09-30T06:37:47Z");
    fXml->NewChild(fAsset, 0, "modified","2015-09-30T06:37:47Z");
    XMLNodePointer_t unit = fXml->NewChild(fAsset, 0, "unit","");
    fXml->NewAttr(unit,0,"meter","0.010000");
    fXml->NewAttr(unit,0,"name","centimeter");
}

void AliColladaBuffer::AddNewLambertMaterial(const char* name, TColor *emission, TColor *ambient, TColor *diffuse, double opacity)
{
    XMLNodePointer_t material = fXml->NewChild(fLibraryMaterials, 0, "material");
    fXml->NewAttr(material,0,"id",name);
    XMLNodePointer_t instanceEffect = fXml->NewChild(material,0,"instance_effect");
    fXml->NewAttr(instanceEffect,0,"url",Form("#%s-fx",name));
    
    XMLNodePointer_t effect = fXml->NewChild(fLibraryEffects, 0, "effect");
    fXml->NewAttr(effect,0,"id",Form("%s-fx",name));
    
    XMLNodePointer_t technique = fXml->NewChild(fXml->NewChild(effect,0,"profile_COMMON"),0,"technique");
    fXml->NewAttr(technique,0,"sid","standard");
    XMLNodePointer_t lambert = fXml->NewChild(technique,0,"lambert");
    
    XMLNodePointer_t emissionColor = fXml->NewChild(fXml->NewChild(lambert,0,"emission"),0,"color",Form("%f %f %f %f",
                                                                                                        emission->GetRed(),
                                                                                                        emission->GetGreen(),
                                                                                                        emission->GetBlue(),
                                                                                                        emission->GetAlpha()));
    XMLNodePointer_t ambientColor = fXml->NewChild(fXml->NewChild(lambert,0,"ambient"),0,"color",Form("%f %f %f %f",
                                                                                                      ambient->GetRed(),
                                                                                                      ambient->GetGreen(),
                                                                                                      ambient->GetBlue(),
                                                                                                      ambient->GetAlpha()));
    XMLNodePointer_t diffuseColor = fXml->NewChild(fXml->NewChild(lambert,0,"diffuse"),0,"color",Form("%f %f %f %f",
                                                                                                      diffuse->GetRed(),
                                                                                                      diffuse->GetGreen(),
                                                                                                      diffuse->GetBlue(),
                                                                                                      diffuse->GetAlpha()));
    
    
    fXml->NewAttr(emissionColor,0,"sid","emission");
    fXml->NewAttr(ambientColor,0,"sid","ambient");
    fXml->NewAttr(diffuseColor,0,"sid","diffuse");
    
    XMLNodePointer_t transparent = fXml->NewChild(lambert,0,"transparent");
    fXml->NewAttr(transparent,0,"opaque","RGB_ZERO");
    XMLNodePointer_t transparentColor = fXml->NewChild(transparent,0,"color",Form("%f %f %f 1.000000",1-opacity,1-opacity,1-opacity));
    fXml->NewAttr(transparentColor,0,"sid","transparent");
    
    XMLNodePointer_t transparency = fXml->NewChild(lambert,0,"transparency");
    XMLNodePointer_t transparencyFloatArray = fXml->NewChild(transparency,0,"float","1.0");
    fXml->NewAttr(transparencyFloatArray,0,"sid","transparency");
    
}

void AliColladaBuffer::SetupScene()
{
    fVisualScene = fXml->NewChild(fLibraryVisualScenes,0,"visual_scene");
    fXml->NewAttr(fVisualScene,0,"id","test");
    
    XMLNodePointer_t visualSceneInstance = fXml->NewChild(fScene,0,"instance_visual_scene","");
    fXml->NewAttr(visualSceneInstance,0,"url","#test");
}

XMLNodePointer_t AliColladaBuffer::AddVisualSceneGeometry(TString name, TString parent)
{
    XMLNodePointer_t parentNode = FindNodeById(parent);
    
    if(parentNode==0)
    {
        if(parent.EqualTo("")) // if parent was not specifed, add to main scene
        {
            cout<<"Parent was not specified. Adding shape to main scene."<<endl;
            parentNode = fVisualScene;
        }
        else // if parent was specified, but doesn't exist yet, create it
        {
            cout<<"As parent:"<<parent<<" was not found, I will create it for you."<<endl;
            parentNode = CreateNode(parent);
        }
    }
    XMLNodePointer_t geom = fXml->NewChild(parentNode,0,"node");
    fXml->NewAttr(geom,0,"id",name.Data());
    return geom;
}

XMLNodePointer_t AliColladaBuffer::FindNodeById(TString id)
{
  for(int i = 0;i < fVisualSceneNodes.size();i++){
        if(id.EqualTo(fXml->GetAttr(fVisualSceneNodes[i],"id"))){
            return fVisualSceneNodes[i];
        }
    }
    return 0;
}

XMLNodePointer_t AliColladaBuffer::CreateNode(TString name,TString parent)
{
    XMLNodePointer_t currentNode = FindNodeById(name);
    
    if(currentNode==0)  // if node with this id doesn't exists yet, create it
    {
        XMLNodePointer_t parentNode = FindNodeById(parent);
        
        if(parentNode==0 || parent.EqualTo(""))
        {
            cout<<"Parent of node to be created doesn't exist:"<<parent<<". Adding to main scene."<<endl;
            parentNode = fVisualScene;
        }
        
        XMLNodePointer_t node = fXml->NewChild(parentNode,0,"node");
        fXml->NewAttr(node,0,"id",name.Data());
        fVisualSceneNodes.push_back(node);
        
        currentNode = node;
    }
    else
    {
        cout<<"Node name is not unique, you should fix that. Skipping this shape."<<endl;
    }
    return currentNode;
}

void AliColladaBuffer::AddShape(const char *name,TEveGeoShapeExtract *shapeExtract,TString parent,const char* materialName)
{
    // get TEveGeoShape object
    TEveGeoShape *shape = TEveGeoShape::ImportShapeExtract(shapeExtract);
    // make buffer3D from TGeoShape and add it to Collada file
    AddShape(name,shape->GetShape()->MakeBuffer3D(),shapeExtract->GetTrans(),parent,materialName);
}

void AliColladaBuffer::AddBox(const char* name,double width, double height, double depth, double translation[3], double rotation[3], double scale[3],TString parent,const char* materialName)
{
    TGeoBBox *box = new TGeoBBox(width/2,height/2,depth/2);
    double trans[16];
    GetTransMatrixFromParams(trans,translation,rotation,scale);
    AddShape(name,box->MakeBuffer3D(),trans,parent,materialName);
    delete box;
}

void AliColladaBuffer::AddShape(const char* name,TGeoShape *shape, double translation[3], double rotation[3], double scale[3],TString parent,const char* materialName)
{
    double trans[16];
    GetTransMatrixFromParams(trans,translation,rotation,scale);
    AddShape(name,shape->MakeBuffer3D(),trans,parent,materialName);
}

void AliColladaBuffer::AddShape(const char* name,TBuffer3D *buffer, double translation[3],double rotation[3], double scale[3],TString parent, const char* materialName)
{
    double trans[16];
    GetTransMatrixFromParams(trans,translation,rotation,scale);
    AddShape(name,buffer,trans,parent,materialName);
}

void AliColladaBuffer::AddShape(const char* name,TBuffer3D *buffer, double trans[16],TString parent,const char* materialName)
{
    TString verticesString = "\n";
    TString normalString = "\n";
    TString trianglesString = "";
    TString transformationString = "";
    vector<double> verticesArray;
    vector<int> trianglesArray;
    set<int> triangle;
    int current=1;
    
    for(unsigned int i=0;i<buffer->NbPnts();i++)
    {
        verticesString += Form("%f %f %f\n",buffer->fPnts[3*i+0],buffer->fPnts[3*i+1],buffer->fPnts[3*i+2]);
        verticesArray.push_back(buffer->fPnts[3*i+0]);
        verticesArray.push_back(buffer->fPnts[3*i+1]);
        verticesArray.push_back(buffer->fPnts[3*i+2]);
    }
    
    for(unsigned int j=0;j<buffer->NbPols();j++) // loop over polygons
    {
        if(buffer->fPols[current] == 3) // if we have triangle
        {
            int a1,a2,b1,b2,c1,c2;
            
            a1 = buffer->fSegs[3*buffer->fPols[current+1]+1];
            a2 = buffer->fSegs[3*buffer->fPols[current+1]+2];
            
            b1 = buffer->fSegs[3*buffer->fPols[current+2]+1];
            b2 = buffer->fSegs[3*buffer->fPols[current+2]+2];
            
            c1 = buffer->fSegs[3*buffer->fPols[current+3]+1];
            c2 = buffer->fSegs[3*buffer->fPols[current+3]+2];
            
            int x,y,z;
            
            // disentangle segments
            if      (a2==b1 && b2 == c1 && c2 == a1){x = a1;y=b1;z=c1;}
            else if (a2==b2 && b1 == c1 && c2 == a1){x = a1;y=b2;z=c1;}
            else if (a2==b1 && b2 == c2 && c1 == a1){x = a1;y=b1;z=c2;}
            else if (a2==b2 && b1 == c2 && c1 == a1){x = a1;y=b2;z=c2;}
            else if (a1==b1 && b2 == c1 && c2 == a2){x = a2;y=b1;z=c1;}
            else if (a1==b2 && b1 == c1 && c2 == a2){x = a2;y=b2;z=c1;}
            else if (a1==b1 && b2 == c2 && c1 == a2){x = a2;y=b1;z=c2;}
            else if (a1==b2 && b1 == c2 && c1 == a2){x = a2;y=b2;z=c2;}
            
            
            trianglesArray.push_back(z);
            trianglesArray.push_back(y);
            trianglesArray.push_back(x);
        }
        else if(buffer->fPols[current] == 4) // if we have quad, make two triangles from it
        {
            int a1,a2,b1,b2,c1,c2,d1,d2;
            
            a1 = buffer->fSegs[3*buffer->fPols[current+1]+1];
            a2 = buffer->fSegs[3*buffer->fPols[current+1]+2];
            
            b1 = buffer->fSegs[3*buffer->fPols[current+2]+1];
            b2 = buffer->fSegs[3*buffer->fPols[current+2]+2];
            
            c1 = buffer->fSegs[3*buffer->fPols[current+3]+1];
            c2 = buffer->fSegs[3*buffer->fPols[current+3]+2];
            
            d1 = buffer->fSegs[3*buffer->fPols[current+4]+1];
            d2 = buffer->fSegs[3*buffer->fPols[current+4]+2];
            
            int x,y,z,u;
            
            // disentangle segments
            if     (a2==b1 && b2 == c1 && c2 == d1 && d2==a1){x=a1;y=b1;z=c1;u=d1;}
            else if(a2==b2 && b1 == c1 && c2 == d1 && d2==a1){x=a1;y=b2;z=c1;u=d1;}
            else if(a2==b1 && b2 == c2 && c1 == d1 && d2==a1){x=a1;y=b1;z=c2;u=d1;}
            else if(a2==b1 && b2 == c1 && c2 == d2 && d1==a1){x=a1;y=b1;z=c1;u=d2;}
            else if(a2==b2 && b1 == c2 && c1 == d1 && d2==a1){x=a1;y=b2;z=c2;u=d1;}
            else if(a2==b1 && b2 == c2 && c1 == d2 && d1==a1){x=a1;y=b1;z=c2;u=d2;}
            else if(a2==b2 && b1 == c2 && c1 == d2 && d1==a1){x=a1;y=b2;z=c2;u=d2;}
            else if(a1==b1 && b2 == c1 && c2 == d1 && d2==a2){x=a2;y=b1;z=c1;u=d1;}
            else if(a1==b2 && b1 == c1 && c2 == d1 && d2==a2){x=a2;y=b2;z=c1;u=d1;}
            else if(a1==b1 && b2 == c2 && c1 == d1 && d2==a2){x=a2;y=b1;z=c2;u=d1;}
            else if(a1==b1 && b2 == c1 && c2 == d2 && d1==a2){x=a2;y=b1;z=c1;u=d2;}
            else if(a1==b2 && b1 == c2 && c1 == d2 && d1==a2){x=a2;y=b2;z=c2;u=d2;}
            else if(a2==b2 && b1 == c1 && c2 == d2 && d1==a1){x=a1;y=b2;z=c1;u=d2;}
            else if(a1==b2 && b1 == c1 && c2 == d2 && d1==a2){x=a2;y=b2;z=c1;u=d2;}
            else if(a1==b2 && b1 == c2 && c1 == d1 && d2==a2){x=a2;y=b2;z=c2;u=d1;}
            else if(a1==b1 && b2 == c2 && c1 == d2 && d1==a2){x=a2;y=b1;z=c2;u=d2;}
            
            // 1st triangle
            trianglesArray.push_back(z);
            trianglesArray.push_back(y);
            trianglesArray.push_back(x);
            
            // 2nd triangle
            trianglesArray.push_back(x);
            trianglesArray.push_back(u);
            trianglesArray.push_back(z);
        }
        
        current = current + buffer->fPols[current] + 2;
        triangle.clear();
    }
    
    for(unsigned long i=0;i<trianglesArray.size();i+=3)
    {
        double normal[3];
        
        double p1[3] = {verticesArray[3*trianglesArray[i]+0],
            verticesArray[3*trianglesArray[i]+1],
            verticesArray[3*trianglesArray[i]+2]};
        
        double p2[3] = {verticesArray[3*trianglesArray[i+1]+0],
            verticesArray[3*trianglesArray[i+1]+1],
            verticesArray[3*trianglesArray[i+1]+2]};
        
        double p3[3] = {verticesArray[3*trianglesArray[i+2]+0],
            verticesArray[3*trianglesArray[i+2]+1],
            verticesArray[3*trianglesArray[i+2]+2]};
        
        GetTrianglesNormal(p1,p2,p3,normal);
        
        if(normal[0] == 0 && normal[1] == 0 && normal[2] == 0)
        {
//            cout<<name<<endl;
            //            return;
        }
        
        for(int j=0;j<3;j++){
            normalString += (float)normal[0];
            normalString += "\t";
            normalString += (float)normal[1];
            normalString += "\t";
            normalString += (float)normal[2];
            normalString += "\n";
        }
    }
    
    for (unsigned long i=0;i<trianglesArray.size(); i++)
    {
        trianglesString += trianglesArray[i];
        trianglesString += " ";
        trianglesString += i;
        trianglesString += " ";
    }
    trianglesString += "\n";
    
    for (int i=0; i<4; i++) {
        transformationString += trans[i];
        transformationString += " ";
        transformationString += trans[i+4];
        transformationString += " ";
        transformationString += trans[i+8];
        transformationString += " ";
        transformationString += trans[i+12];
        transformationString += " ";

    }
    transformationString+= "\n";
    
    AddShapeToCollada(name,verticesString,verticesArray.size(),normalString,trianglesString,trianglesArray.size(),transformationString,materialName,parent);
}


void AliColladaBuffer::AddShapeToCollada(const char* name,
                                       const char* verticesString,int verticesSize,
                                       const char* normalString,
                                       const char* trianglesString,int trianglesSize,
                                       const char* transformationString,
                                       const char* materialName,
                                       TString parent)
{
    XMLNodePointer_t geometry = fXml->NewChild(fLibraryGeometries, 0, "geometry");
    fXml->NewAttr(geometry,0,"id",Form("%s-lib",name));
    XMLNodePointer_t mesh_node = fXml->NewChild(geometry, 0, "mesh");
    
    // position
    XMLNodePointer_t sourcePosition = fXml->NewChild(mesh_node, 0, "source");
    fXml->NewAttr(sourcePosition,0,"id",Form("%s-POSITION",name));
    
    XMLNodePointer_t floatArrayPosition = fXml->NewChild(sourcePosition, 0, "float_array",verticesString);
    fXml->NewAttr(floatArrayPosition,0,"id",Form("%s-POSITION-array",name));
    fXml->NewAttr(floatArrayPosition,0,"count",Form("%d",verticesSize));
    XMLNodePointer_t techniquePosition = fXml->NewChild(sourcePosition,0,"technique_common");
    XMLNodePointer_t accessorPosition = fXml->NewChild(techniquePosition,0,"accessor");
    fXml->NewAttr(accessorPosition,0,"source",Form("#%s-POSITION-array",name));
    fXml->NewAttr(accessorPosition,0,"count",Form("%d",verticesSize/3));
    fXml->NewAttr(accessorPosition,0,"stride","3");
    fXml->NewAttr(fXml->NewChild(accessorPosition,0,"param"),0,"type","float");
    fXml->NewAttr(fXml->NewChild(accessorPosition,0,"param"),0,"type","float");
    fXml->NewAttr(fXml->NewChild(accessorPosition,0,"param"),0,"type","float");
    
    // normals
    XMLNodePointer_t sourceNormal = fXml->NewChild(mesh_node,0,"source");
    fXml->NewAttr(sourceNormal,0,"id",Form("%s-Normal0",name));
    
    XMLNodePointer_t floatArrayNormal = fXml->NewChild(sourceNormal,0,"float_array",normalString);
    fXml->NewAttr(floatArrayNormal,0,"id",Form("%s-Normal0-array",name));
    fXml->NewAttr(floatArrayNormal,0,"count",Form("%d",trianglesSize*3));
    
    XMLNodePointer_t techniqueNormal = fXml->NewChild(sourceNormal,0,"technique_common");
    XMLNodePointer_t accessorNormal = fXml->NewChild(techniqueNormal,0,"accessor");
    fXml->NewAttr(accessorNormal,0,"source",Form("#%s-Normal0-array",name));
    fXml->NewAttr(accessorNormal,0,"count",Form("%d",trianglesSize));
    fXml->NewAttr(accessorNormal,0,"stride","3");
    
    // vertices
    XMLNodePointer_t vertices = fXml->NewChild(mesh_node,0,"vertices");
    fXml->NewAttr(vertices,0,"id",Form("%s-VERTEX",name));
    XMLNodePointer_t inputVertices = fXml->NewChild(vertices,0,"input");
    fXml->NewAttr(inputVertices,0,"semantic","POSITION");
    fXml->NewAttr(inputVertices,0,"source",Form("#%s-POSITION",name));
    
    // triangles
    XMLNodePointer_t triangles = fXml->NewChild(mesh_node,0,"triangles");
    fXml->NewAttr(triangles,0,"count",Form("%d",trianglesSize/3));
    fXml->NewAttr(triangles,0,"material",materialName);
    
    XMLNodePointer_t inputTrianglesVertex = fXml->NewChild(triangles,0,"input");
    fXml->NewAttr(inputTrianglesVertex,0,"semantic","VERTEX");
    fXml->NewAttr(inputTrianglesVertex,0,"offset","0");
    fXml->NewAttr(inputTrianglesVertex,0,"source",Form("#%s-VERTEX",name));
    
    XMLNodePointer_t inputTrianglesNormal = fXml->NewChild(triangles,0,"input");
    fXml->NewAttr(inputTrianglesNormal,0,"semantic","NORMAL");
    fXml->NewAttr(inputTrianglesNormal,0,"offset","1");
    fXml->NewAttr(inputTrianglesNormal,0,"source",Form("#%s-Normal0",name));
    
    fXml->NewChild(triangles,0,"p",trianglesString);
    
    // add shape to scene
    XMLNodePointer_t geomNode = AddVisualSceneGeometry(name,parent);
    
    if(geomNode==0){
        cout<<"Couldn't add geometry to scene. Ignoring..."<<endl;
        return;
    }
    fXml->NewChild(geomNode,0,"matrix",transformationString);// transformation matrix
    
    XMLNodePointer_t instanceGeometry = fXml->NewChild(geomNode,0,"instance_geometry");
    fXml->NewAttr(instanceGeometry,0,"url",Form("#%s-lib",name));
    
    XMLNodePointer_t technique = fXml->NewChild(fXml->NewChild(instanceGeometry,0,"bind_material"),0,"technique_common");
    XMLNodePointer_t material = fXml->NewChild(technique,0,"instance_material");
    fXml->NewAttr(material,0,"symbol",materialName);
    fXml->NewAttr(material,0,"target",Form("#%s",materialName));
}

void AliColladaBuffer::SaveAs(const char* filename)
{
    XMLDocPointer_t xmldoc = fXml->NewDoc();
    fXml->DocSetRootElement(xmldoc, fRootNode);
    fXml->SaveDoc(xmldoc, filename);
    fXml->FreeDoc(xmldoc);
}

void AliColladaBuffer::GetTrianglesNormal(double p0[3],double p1[3], double p2[3], double *result)
{
    double u[3],v[3];
    
    for (int i=0; i<3; i++) {
        u[i] = p1[i]-p0[i];
        v[i] = p2[i]-p0[i];
    }
    double n[3];
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
    
    double den = TMath::Abs(n[0])+TMath::Abs(n[1])+TMath::Abs(n[2]);
    
    if(TMath::Abs(den)<0.00001){
//        cout<<"\nCannot calculate normal, it's not a triangle\n"<<endl;
//        cout<<p0[0]<<"\t"<<p0[1]<<"\t"<<p0[2]<<endl;
//        cout<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<endl;
//        cout<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<endl;
        
        result[0] = 0.;
        result[1] = 0.;
        result[2] = 0.;
        return;
    }
    
    result[0] = n[0]/den;
    result[1] = n[1]/den;
    result[2] = n[2]/den;
}

void AliColladaBuffer::GetTransMatrixFromParams(double *trans, double *translation, double *rotation, double *scale)
{
    trans[0] =  cos(rotation[1])*cos(rotation[2])*scale[0];
    trans[1] =  cos(rotation[1])*sin(rotation[2])*scale[0];
    trans[2] = -sin(rotation[1])*scale[0];
    trans[3]=  0;
    trans[4] = (cos(rotation[2])*sin(rotation[0])*sin(rotation[1])-cos(rotation[0])*sin(rotation[2]))*scale[1];
    trans[5] = (cos(rotation[0])*cos(rotation[2])+sin(rotation[0])*sin(rotation[1])*sin(rotation[2]))*scale[1];
    trans[6] =  cos(rotation[1])*sin(rotation[0])*scale[1];
    trans[7]=  0;
    trans[8] = (cos(rotation[0])*cos(rotation[2])*sin(rotation[1])+sin(rotation[0])*sin(rotation[2]))*scale[2];
    trans[9] = (cos(rotation[0])*sin(rotation[1])*sin(rotation[2])-cos(rotation[2])*sin(rotation[0]))*scale[2];
    trans[10]=  cos(rotation[0])*cos(rotation[1])*scale[2];
    trans[11]=  0;
    trans[12] =  translation[0];
    trans[13] =  translation[1];
    trans[14]=  translation[2];
    trans[15]=  1;
}

