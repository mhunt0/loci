#include <iostream>
#include <cstdio>
#include <stdlib.h>	
#include <vector>
#include <map>
#include <string>
#include "ADF.h"
#include "vogtools.h"
#include <Loci.h>
#include <sstream>
#include <vector>

using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;

typedef int int32;
typedef long int int64;

int readNodeis( double parentID,
                const char *name, int32** value );

void checkError( int err, double nodeID, string str )
{
  if(err!= -1){
    cerr << str << " (error " << err << ")";
    char nodeName[ADF_NAME_LENGTH]={'\0'};
    ADF_Get_Name(nodeID, nodeName, &err);
    cerr << " node(" << nodeName <<")"<<endl;
    exit(1);
  }
}


inline int type2size(string type){
  if(type=="MT") return 0;
  if(type=="R4" ||type=="I4") return 4;
  if(type=="C1")return 1;
  if(type=="R8" ||type=="I8")return 8;
  cerr<<"Bad type: " << type << endl;
  return 0;
}
void getDimensions( double node,
                    int64 *nDims, int64 *dims, bool ignoreExtended=false)
{
  int i;
  int32 nDims32, dims32[ADF_MAX_DIMENSIONS];
  int err;
  double nodeEx;

  {
     
    if (!nDims) return;
    *nDims = 0;
    ADF_Get_Node_ID(node, "ExtendedSize", &nodeEx, &err);
    if(err ==-1 && !ignoreExtended)
      { /* This node is an extendend node */
        /* Check that this node is really a 1D array like it should be */
        i = 0;
        ADF_Get_Number_of_Dimensions(nodeEx, &i, &err);
        if (err!=-1 || i != 1)
          return;

        /* Now find out how big it is */
        ADF_Get_Dimension_Values(nodeEx, dims32, &err);
        if (err!=-1)
          return;
        *nDims = int64(dims32[0]);
        
        /* Now read in the number of dimensions */
        if (dims)
          {
            ADF_Read_All_Data(nodeEx, (char *)dims, &err);
            if (err!=-1)
              return;
          }
      }
    else
      {   /* Just a regular node;  call the regular ADF functions */
        ADF_Get_Number_of_Dimensions(node, &nDims32, &err);
        if (err!=-1)
          return;
        *nDims = int64(nDims32);
        
        if (dims)
          {
            ADF_Get_Dimension_Values(node, dims32, &err);
            if (err == 32 /* No data */ ||
                err == 27 /* Dimension is zero */)
              {
                *nDims = int64(0);
                dims[0] = int64(0);
              }
            else if (err!=-1)
              return;

            for (i = 0;  i < *nDims;  ++i)
              dims[i] = int64(dims32[i]);
          }
      }

    return;
  } /* end CHECK_ERROR scope */
}

  

string getType(double nodeID){
  int err;
  char data_type[ADF_DATA_TYPE_LENGTH]={'\0'};
  ADF_Get_Data_Type(nodeID, data_type, &err);
  checkError(err, nodeID, "Error getting type");
  data_type[2]='\0';
  return string(data_type);
}

/*if it is a regular
  node this will just be the dimensions of the node;  if it is
  an extended node, these will be the dimensions of the extended
  data.) */
int64 getSize(double nodeID, bool ignoreExtended=false){
  int64 num_dims, dim_values[ADF_MAX_DIMENSIONS];
  string type = getType(nodeID);
  if(type2size(type)==0)return 0;
  getDimensions(nodeID, &num_dims, dim_values, ignoreExtended);
  if(num_dims==0)return 0;
  int64 size = 1; 
  for(int i = 0; i < num_dims; i++)size *=dim_values[i];
  return size;
}


void readAllData(double const nodeID, char* buffer)
{
  int err;
  ADF_Read_All_Data(nodeID, buffer, &err);
  checkError(err, nodeID, "Error reading data");


  double tmpID;
  ADF_Get_Node_ID(nodeID, "ExtendedSize", &tmpID, &err);
  // 29 - Specified child is not a child of the specified parent
  if(err == 29){
    return;
  }
 
  //is extended data
  
  // Read the rest of the extended data.
  long long size = getSize(nodeID, true);
  string type = getType(nodeID);
  char * local_buffer = buffer + type2size(type) * size;
  
  int num_children;
  ADF_Number_of_Children(nodeID, &num_children, &err);
  checkError(err, nodeID, "Error getting number of children");
  
 
  int inum_ret;
  char * names = new char[34*num_children];
  ADF_Children_Names (nodeID, 1, num_children, 33, &inum_ret, names,
                      &err);
  checkError(err, nodeID, "Error getting children names");

  char * cp;
  char * end = names + 34*inum_ret;
  for (cp=names; cp!=end; cp+=34)
    {
      static char const * extendedData = "ExtendedData-";
      if (strncmp(cp, extendedData, strlen(extendedData)))
        continue;

      double ID;
      ADF_Get_Node_ID(nodeID, cp, &ID, &err);
      checkError(err, nodeID, string("Error getting node ID ")+string(extendedData));

      ADF_Read_All_Data(ID, local_buffer, &err);
      checkError(err, nodeID, "Error reading data");

      size = getSize(ID, true);
      local_buffer += type2size(type) * size;
    }

  delete [] names;
}





string readNodestr( double parentID,
                    const char *name)
{
  double childID;
  int err;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  if(err!=-1)return string();//node doesn't exist
  
  int size = getSize(childID);
  if(size==0) return string();//node exist and size is 0
  
  char value[ADF_FILENAME_LENGTH];
  readAllData(childID, &value[0]);
  value[size]='\0';
 
  return string(value);
}

int readNodefs( double parentID,
                const char *name, float **value )
{
  double nodeID;
  int err;
   
  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID, string("Error in getting node ID of ")+string(name));

  int size = getSize(nodeID);
  if(size==0) return 0;
  *value = new float[size];
  readAllData(nodeID, (char*)(*value));
  return size;
}


int readNodeis( double parentID,
                const char *name, int32** value )
{
  double nodeID;
  int err;
 
  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID, string("Error in getting node ID of ")+string(name));
 
  int size = getSize(nodeID);
  if(size==0) return 0;
  *value = new int32[size];
  readAllData(nodeID, (char*)(*value));
  return size; 
  
}

int readNodeds( double parentID,
                const char *name, double **value )
{
  double nodeID;
  int err;
  
  ADF_Get_Node_ID(parentID, name, &nodeID, &err);
  checkError(err, parentID, string("Error in getting node ID of ")+string(name));

  int size = getSize(nodeID);
  if(size==0) return 0;
  *value = new double[size];
  readAllData(nodeID, (char*)(*value));
  
  return size;   
}





int32 readNodei32( double parentID,
                   const char *name)
{
  double childID;
  int err;
  int32 value;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID, string("Error in  getting node ID of ")+string(name));
  ADF_Read_All_Data(childID, (char*)&value, &err);
  checkError(err, childID, "Error reading data");
  return value;
}

float readNodef( double parentID,
                 const char *name)
{
  double childID;
  int err;
  float value;
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID, string("Error in getting  node ID of ")+string(name));
  ADF_Read_All_Data(childID, (char*)&value, &err);
  checkError(err, childID, "Error reading data");
  return value;
}

string getDataType( double parentID,
                    const char *name)
{
  double childID;
  int err;
  
  ADF_Get_Node_ID(parentID, name, &childID, &err);
  checkError(err, parentID, string("Error in  getting node ID of ")+string(name));
  return getType(childID);
}

int getNumNodes(double vertices){
  double nodeID;
  int err;
  int64 num_dims, dim_values[ADF_MAX_DIMENSIONS];
  
  ADF_Get_Node_ID(vertices, "Coordinates", &nodeID, &err);
  checkError(err, vertices, "Error in getting node ID of Coordinates");
  getDimensions(nodeID,&num_dims, dim_values);
  
  if(num_dims!=2) return 0;
  return dim_values[1];
}



int getNumFaces(double topoID){

  int numFaces = readNodei32(topoID, "InternalFaces/NumFaces");
  int err;
  int numChildren;
  ADF_Number_of_Children(topoID, &numChildren, &err);
  checkError(err, topoID, "Can not get number of children of the first child of States");
 
  char name[ADF_NAME_LENGTH]={'\0'};
  int num_ret = 0;
  
  for(int start = 1; start<=numChildren; start++){
    ADF_Children_Names (topoID,start, 1, ADF_NAME_LENGTH,&num_ret,name,&err);
    checkError(err, topoID, string("Can not get node ID of ")+string(name));
    double bfaceID;
    if(string(name).substr(0,14)=="BoundaryFaces-"){
      ADF_Get_Node_ID(topoID, name, &bfaceID, &err);
      checkError(err,topoID, string("Can not get node ID of ")+ string(name));
      numFaces += readNodei32(bfaceID,"NumFaces");
    }
  }
  
  return numFaces;
}


//read in pos  
void readVertices( double vertices,  store<vector3d<double> >& pos)
{
  
  string type = getDataType(vertices,"Coordinates");
  
  int numNodes=getNumNodes(vertices);
 
  float scale = readNodef(vertices, "ScaleFactor");
   
  int err;
  int32 mapIndex = readNodei32(vertices,"MapId"); 
  char nodeName[ADF_NAME_LENGTH]={'\0'};
  sprintf(nodeName, "/Maps/Map-%d/IdMap", mapIndex);
  double root;
  ADF_Get_Root_ID(vertices, &root, &err);
  checkError(err, vertices, "Error getting root ID");
  
  int32* mapData =0;
  readNodeis(root, nodeName, &mapData);
  if(mapData==0){
    cerr<< " Error reading vertex map " << endl;
    exit(1);
  }
  
  
  if(type=="R4"){
    float* verts =0;
    readNodefs(vertices, "Coordinates", &verts);
    //scale the coordiates and put into vdata
    if(verts==0){
      cerr<< " Error reading Coordinates " << endl;
      exit(1);
    }
    for ( int i = 0;  i < numNodes;  ++i)
      {
        vector3d<double> p;
        p.x = verts[3 * i    ] * scale;
        p.y = verts[3 * i + 1  ] * scale;
        p.z = verts[3 * i  + 2  ] * scale;
        pos[mapData[i]-1]  = p;
        
      }
    if(verts)delete [] verts;
  }else{
    double* verts =0;
    readNodeds(vertices, "Coordinates", &verts);
    for ( int i = 0;  i < numNodes;  ++i)
      {
        vector3d<double> p;
        p.x = verts[3 * i ] * (double)scale;
        p.y = verts[3 * i + 1 ] *(double) scale;
        p.z = verts[3 * i  + 2 ] * (double)scale;
        pos[mapData[i]-1]  = p;

      }
    if(verts)  delete [] verts;
  }

  if(mapData) delete [] mapData;
     
}
  

//read in facebased  topology 
void readMesh(  double topoID,
                vector<int32* > &allFaces,
                vector<int> &faceSizes,
                vector<int32* > &allFaceCells,
                vector<int> &faceCellSizes, 
                vector<pair<int,string> > &surf_ids)
{
  int err; 
  double root;
  bool isLocal = false;
  ADF_Get_Root_ID(topoID, &root, &err);
  checkError(err, topoID, "Error getting root ID");
  if(readNodestr(root, "/Meshes/Space")=="Local") isLocal = true;
  int32* vertexMapData = 0;
  int32* mapData = 0;

  int32* faceData = 0;
  int localFaceSize = readNodeis(topoID, "InternalFaces/Vertices", &faceData);
  faceSizes.push_back(localFaceSize);
 
  if(isLocal){
    char nodeName[ADF_NAME_LENGTH]={'\0'};
    ADF_Get_Name(topoID, nodeName, &err);
    checkError(err, topoID, "Error getting node name");
    int index = atoi(string(nodeName).substr(18).c_str());
    sprintf(nodeName, "/Meshes/Vertices-%d/MapId", index);
    int mapIndex = readNodei32(root, nodeName);
    
    //read vertices map data, assume FaceBasedTopology-i and  Vertices-i share the same vertex map
    sprintf(nodeName, "/Maps/Map-%d/IdMap", mapIndex);
   
    readNodeis(root, nodeName, &vertexMapData);
    if(vertexMapData==0){
      cerr<< " Error reading vertex map " << endl;
      exit(1);
    }
   
    //map face data
    int pointer = 0;
    while(pointer < localFaceSize){
      for(int  i = 1; i <= faceData[pointer]; i++){
        faceData[pointer+i]=vertexMapData[faceData[pointer+i]-1];
      
      } 
      
      pointer += faceData[pointer]+1;
    }
    
  }
  allFaces.push_back(faceData);

  
  int32* faceCellData = 0;
  faceCellSizes.push_back(readNodeis(topoID, "InternalFaces/Cells", &faceCellData));
  if(isLocal){
    
    char nodeName[ADF_NAME_LENGTH]={'\0'};  
    sprintf(nodeName, "Cells/MapId");
    int mapIndex = readNodei32(topoID, nodeName);
    
    sprintf(nodeName, "/Maps/Map-%d/IdMap", mapIndex);
    readNodeis(root, nodeName, &mapData);
    if(mapData==0){
      cerr<< " Error reading vertex map " << endl;
      exit(1);
    }
   
    //map the data
    for(int  i = 0; i < faceCellSizes.back(); i++){
      faceCellData[i]=mapData[faceCellData[i]-1];
    } 
    
  }

  allFaceCells.push_back(faceCellData);

  

  int numChildren;
  ADF_Number_of_Children(topoID, &numChildren, &err);
  checkError(err, topoID, "Can not get number of children of the first child of States");
 
  char name[ADF_NAME_LENGTH]={'\0'};
  int num_ret = 0;

  int numBoundaries = readNodei32(topoID, "NumBoundaryTypes");
  
  for(int start = 1; start<=numChildren; start++){
    ADF_Children_Names (topoID,start, 1, ADF_NAME_LENGTH,&num_ret,name,&err);
    checkError(err, topoID, string("Can not get node ID of ")+string(name));
    
    double bfaceID;
    if(string(name).substr(0,14)=="BoundaryFaces-"){
      ADF_Get_Node_ID(topoID, name, &bfaceID, &err);
      checkError(err,topoID, string("Can not get node ID of ")+ string(name));
      
      string index =string(name).substr(14);
      int btype= atoi(index.c_str());
      if(btype==0)btype = numBoundaries;//some index start with 0;
      int32* bfaceData = 0;
      faceSizes.push_back(readNodeis(bfaceID, "Vertices", &bfaceData));
      if(isLocal && vertexMapData !=0){
        int pointer = 0;
        int localFaceSize = faceSizes.back();
        while(pointer < localFaceSize){
          for(int  i = 1; i <= bfaceData[pointer]; i++){
            bfaceData[pointer+i]=vertexMapData[bfaceData[pointer+i]-1];
          } 
          
          pointer += bfaceData[pointer]+1;
        }
      }
      allFaces.push_back(bfaceData);
    
      int32* bfaceCellData = 0;
      int numFaces = readNodeis(bfaceID, "Cells", &bfaceCellData);
      faceCellSizes.push_back(2*numFaces);
      if(isLocal && mapData != 0){
        //map the data
        for(int  i = 0; i < numFaces; i++){
          bfaceCellData[i]=mapData[bfaceCellData[i]-1];
        } 
      }



      
      int32* tmpfacecells = new int32[2*numFaces];
      for( int i =0; i <numFaces; i++){
        tmpfacecells[2*i] = bfaceCellData[i];
        tmpfacecells[2*i+1] = -btype;
      }

      allFaceCells.push_back(tmpfacecells);
      delete [] bfaceCellData;
      
      //      string bcstr = readNodestr(bfaceID, "Label");
      //if(bcstr.size()==0){
      char bcname[512];
      if(btype >= 0)
        sprintf(bcname,"BC_%d",btype) ;
      else
        sprintf(bcname,"BC_m%d",-btype) ;
      
      string bcstr=string(bcname);
      // }
      
      surf_ids.push_back(pair<int, string>(btype, bcstr));
       
    }
  
  }
  if(vertexMapData) delete [] vertexMapData;
  if(mapData) delete [] mapData;
}
  
 

//for each processor, get its vertices and topology ID, might open another file
void getMeshID(double processorID, double* verticesID, double* topoID){
  int numChildren;
  int err;


  double root;
  ADF_Number_of_Children(processorID, &numChildren, &err);
  checkError(err, processorID,"Can not get number of children of the first child of States");


  //get vertices
  if(verticesID) {
   
    string verticesFileName = readNodestr(processorID, "VerticesFile"); 
    if(verticesFileName.size()>0){
      ADF_Database_Open(verticesFileName.c_str(), "READ_ONLY",
                        "NATIVE", &root, &err);
      checkError(err, root, "Error opening file");
    }else{
      ADF_Get_Root_ID(processorID, &root, &err);
      checkError(err, processorID, "Error get root ID");
    }
    
    int verticesIndex = readNodei32(processorID, "VerticesId"); 

    char nodeName[ADF_NAME_LENGTH]={'\0'};
    sprintf(nodeName, "/Meshes/Vertices-%d", verticesIndex);

    ADF_Get_Node_ID(root, nodeName, verticesID, &err);
    checkError(err, root, string("Error getting node ID ")+string(nodeName));
  }
  //get topoID
  if(topoID){
    
    string topoFileName = readNodestr(processorID, "TopologyFile"); 
    if(topoFileName.size()>0){
      ADF_Database_Open(topoFileName.c_str(), "READ_ONLY",
                        "NATIVE", &root, &err);
      checkError(err, root,"Error opening file");
    }else{
      ADF_Get_Root_ID(processorID, &root, &err);
      checkError(err, processorID, "Error get root ID");
    }
    
    
    int topoIndex = readNodei32(processorID, "TopologyId"); 
  
    char nodeName[ADF_NAME_LENGTH]={'\0'};
    sprintf(nodeName, "/Meshes/FaceBasedTopology-%d", topoIndex);
    
    ADF_Get_Node_ID(root, nodeName, topoID, &err);
    checkError(err, root, string("Error getting node ID ")+ string(nodeName));
  }

  
}


  
//facemap is not read, assume no duplicate faces between processors
int main( int argc, char *argv[])
{
  
  
  if (argc < 2 || argc > 3){
    cerr << "Usage: ccm2vog input.ccm[g] [output.vog]"<<endl;
    exit(1);
  }
  
  
  string infile=string(argv[1]);
  string outfile;
  string case_name;
  
  //only accept "*.ccm" or "*.ccmg" file as input
  size_t p = infile.rfind('.');
  if(p != string::npos){
    if(infile.substr(p,infile.size())!=".ccm" &&infile.substr(p,infile.size())!=".ccmg"){
      cerr << "Usage: ccm2vog input.ccm[g] [output.vog]"<<endl;
      exit(1);
    }
    case_name = infile.substr(0, p);
  }else{
    cerr << "Usage: ccm2vog input.ccm[g] [output.vog]"<<endl;
    
    exit(1); 
  }
  
    
  if(argc==2)outfile = case_name +".vog"; 
  else outfile=argv[2];
  

  using namespace Loci ;
  using namespace VOG ;
  bool optimize = true ;
  Loci::Init(&argc,&argv) ;//Loci initialize
  
  //Loci data structures
  store<vector3d<double> > pos;
  Map cl, cr;
  store<int> facecount;
  multiMap face2node;

  //open file
  double root;
  int err;
  ADF_Database_Open(argv[1], "READ_ONLY",
                    "NATIVE", &root, &err);
  checkError(err,  root,string("Can not open ")+string(argv[1])+ string(" for reading"));



  //get node States
  double statesID;
  ADF_Get_Node_ID(root,"States", &statesID, &err);
  checkError(err, root, "Can not get node ID of States");

  //get first State
  int num_ret;
  char name[ADF_NAME_LENGTH]={'\0'};
  double stateID;
  ADF_Children_Names (statesID,1,1,ADF_NAME_LENGTH,&num_ret, name, &err);
  checkError(err, statesID, "Can not children names of States");
  
  ADF_Get_Node_ID(statesID, name, &stateID, &err);
  checkError(err, stateID, "Can not get the node ID of the first child of States");
  
  //get all Processors
  int numChildren;
  ADF_Number_of_Children(stateID, &numChildren, &err);
  checkError(err, stateID, "Can not get number of children of the first child of States");
 
 
  vector<double> processorIDs;
  for(int start = 1; start<=numChildren; start++){
    ADF_Children_Names (stateID,start, 1,ADF_NAME_LENGTH,&num_ret,name,&err);
    checkError(err, stateID, "Can not get node ID of ");
  
    double processorID;
    string nodeName(name);
    if(string(name).substr(0,10)=="Processor-"){
      ADF_Get_Node_ID(stateID, name, &processorID, &err);
      checkError(err,stateID, string("Can not get node ID of ")+ nodeName);
      processorIDs.push_back(processorID);
    }
  }
  
  //go through all processors, compute total number of nodes
  int totalNumNode = 0;
  for(unsigned int i = 0; i < processorIDs.size(); i++){
    double verticesID;
    getMeshID(processorIDs[i], &verticesID, NULL);  
    totalNumNode += getNumNodes(verticesID);
  }
  
 
  //allocate pos read pos, pos always use global number
  entitySet nodes = interval(0, totalNumNode-1);
  pos.allocate(nodes);
  for(unsigned int i = 0; i < processorIDs.size(); i++){
    double verticesID;
    getMeshID(processorIDs[i], &verticesID, NULL);  
    readVertices(verticesID, pos);
  }
  
  
  //go through all processors, compute total number of faces
  int totalNumFace = 0;
  for(unsigned int i = 0; i < processorIDs.size(); i++){
    double topoID;
    getMeshID(processorIDs[i], NULL, &topoID);  
    totalNumFace += getNumFaces(topoID);
  }
  

  //allocate maps
  entitySet facesSet = interval(0, totalNumFace-1);
  cl.allocate(facesSet);
  cr.allocate(facesSet);
  facecount.allocate(facesSet);
  
  //read in all face info
  vector<int32*> allFaces;
  vector<int> faceSizes;
  vector<int32*> allFaceCells;
  vector<int32> faceCellSizes;
  vector<pair<int,string> > surf_ids;
  for(unsigned int i = 0; i < processorIDs.size(); i++){
    double topoID;
    getMeshID(processorIDs[i], NULL, &topoID);  
    readMesh( topoID,
              allFaces,
              faceSizes,
              allFaceCells,
              faceCellSizes, 
              surf_ids);
  }
  
  //close data base
  ADF_Database_Close(root, &err); 

 
  
 
  //put face info into maps
  entitySet::const_iterator ei = facesSet.begin();
  for(unsigned int fid = 0; fid < allFaces.size(); fid++){
    int pointer = 0;
   
    while(pointer < faceSizes[fid])
      { 
        
        facecount[*ei] = allFaces[fid][pointer];
        pointer += allFaces[fid][pointer]+1;
        ei++;
       
      }
  }
  face2node.allocate(facecount);
 
  ei = facesSet.begin();
  for(unsigned int fid = 0; fid < allFaces.size(); fid++){
    int pointer = 0;
    while(pointer < faceSizes[fid]){
      for(int  i = 0; i < facecount[*ei]; i++){
        face2node[*ei][i]=allFaces[fid][pointer+i+1]-1;
      } 
      
      pointer += allFaces[fid][pointer]+1;
      ei++;
    }
  }
 
  ei = facesSet.begin();
  for(unsigned int fid = 0; fid < allFaces.size(); fid++){
    int pointer = 0;
    while(pointer < faceCellSizes[fid]){
      cl[*ei] = allFaceCells[fid][pointer];
      cr[*ei] = allFaceCells[fid][pointer+1];

      pointer += 2;
      ei++;
    }
  }

  //release memory
  for(unsigned int fid = 0; fid < allFaces.size(); fid++)delete [] allFaces[fid];
  for(unsigned int fid = 0; fid < allFaceCells.size(); fid++)delete [] allFaceCells[fid];
   
 
  //Loci takes over
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;
  
  if(optimize) {
    if(MPI_rank == 0) 
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }
 
  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;
  Loci::writeVOG(outfile, pos, cl, cr, face2node ,surf_ids) ;
 
  Loci::Finalize() ;
  return 0 ;
    
}



