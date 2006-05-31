#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <algorithm> 

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;


//---------------------Array----------------------//
template <class T,size_t n> class Array {
    T x[n] ;
  public:
    typedef T * iterator ;
    typedef const T * const_iterator ;
    
    Array() {} ;
    Array(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; } 
    Array<T,n> &operator=(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 
    
    Array<T,n> &operator +=(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
    Array<T,n> &operator -=(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
    Array<T,n> &operator *=(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
    Array<T,n> &operator /=(const Array<T,n> &v)
        { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }
    
    T &operator[](size_t indx) { return x[indx]; }
    const T &operator[](size_t indx) const { return x[indx] ; }
    
    iterator begin() { return &x[0] ; }
    iterator end() { return begin()+n ; }
    const_iterator begin() const { return &x[0] ; }
    const_iterator end() const { return begin()+n ; }
    
    size_t size() const  { return n ; }
} ;

template <class T,size_t n> inline std::ostream &
operator<<(std::ostream &s, const Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
}

template <class T,size_t n> inline std::istream &
operator>>(std::istream &s, Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s >> v[i] ;
    return s ;
}


bool reverse_byteorder = false ;

void check_order() {
    static int test = 15 ;
    char *p = (char *)&test ;
    if(int(*p) == test) {
        reverse_byteorder = true ;
    }
}


void ug_io_reverse_byte_order
 (void * Data,
  size_t Size,
  int Number)

{
 
/*
 * Set file format and host to big or little endian byte ordering.
 * 
 */

  char *Data_Byte_Ptr;
  char Temp_Data_Byte;

  int Byte_Index, Index, Number_of_Bytes, Reverse_Byte_Index;

  Number_of_Bytes = int(Size);

  Data_Byte_Ptr = (char *) Data;

  for (Index = 0; Index < Number; ++Index)
  {
    Reverse_Byte_Index = Number_of_Bytes;

    for (Byte_Index = 0; Byte_Index < Number_of_Bytes/2; ++Byte_Index)
    {
      --Reverse_Byte_Index;

      Temp_Data_Byte = Data_Byte_Ptr[Byte_Index];

      Data_Byte_Ptr[Byte_Index] = Data_Byte_Ptr[Reverse_Byte_Index];

      Data_Byte_Ptr[Reverse_Byte_Index] = Temp_Data_Byte;
    }

    Data_Byte_Ptr += Number_of_Bytes;
  }

  return;
}

typedef Array<int,4> tri_info ;
typedef Array<int,5> quad_info ;

inline bool tri_sort(const tri_info &a1,const tri_info &a2) {
  return a1[0] < a2[0] || (a1[0]==a2[0] && a1[1]<a2[1]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2]<a2[2]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] < a2[3]) ;
}

inline bool quad_sort(const quad_info &a1,const quad_info &a2) {
  return a1[0] < a2[0] || (a1[0]==a2[0] && a1[1]<a2[1]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2]<a2[2]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] < a2[3]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] == a2[3]
     && a1[4] < a2[4]) ;
}

inline bool bnd_tri_sort(const Array<int,5> &a1,
                         const Array<int,5> &a2) {
  return a1[4] < a2[4] ;
}

inline bool bnd_quad_sort(const Array<int,6> &a1,
                          const Array<int,6> &a2) {
  return a1[5] < a2[5] ;
}

int main(int ac, char* av[]) {

  const char *filename ;
  std::string tmp_str ;
  bool binary = 0;
  if(ac == 3) {
    tmp_str.append(av[2]) ;
    if(!strcmp(av[1],"-b")) 
      binary = 1 ;
    else {
      cerr << "Right now the only option supported is '-b' which is the binary ugrid format" << endl ;
      exit(-1) ;
    }
  } else if(ac == 2) {
      tmp_str.append(av[1]) ;
  } else {
    cerr << "solidMesh2xdr requires one argument" << endl 
         << " (the -b flag may be specified for binary files)" << endl;
      exit(-1) ;
  }

  check_order() ;

  int loc = 0;
  loc = tmp_str.find('.') ;
  std::string new_str = tmp_str.substr(0, loc) ;
  filename = new_str.c_str() ;

  FILE* IFP ;
  char buf[512] ;
  if(!binary) {
    struct stat fstat ;
    sprintf(buf,"%s.ugrid",filename) ;
    if(stat(buf,&fstat)<0) {
      binary = true ;
    }
  }
    
  if(!binary)
    sprintf(buf,"%s.ugrid",filename) ;
  else
    sprintf(buf,"%s.b8.ugrid",filename) ;
  if(!binary)
    IFP = fopen(buf, "r") ;
  else
    IFP = fopen(buf, "rb") ;
  if(IFP == NULL) {
    cerr << "can't open '" << buf << "'" << endl ;
    exit(-1) ;
  }

  int num_nodes, num_sf_trias, num_sf_quads ;
  int num_vol_tets, num_vol_pents5, num_vol_pents6, num_vol_hexs ;
  
  if(!binary) {
    fscanf(IFP, "%d%d%d", &num_nodes, & num_sf_trias, & num_sf_quads) ;
    fscanf(IFP, "%d%d%d%d", &num_vol_tets, &num_vol_pents5, &num_vol_pents6, &num_vol_hexs) ;
  }
  else {
    fread(&num_nodes, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_nodes,sizeof(int),1) ;
    fread(&num_sf_trias, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_sf_trias,sizeof(int),1) ;
    fread(&num_sf_quads, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_sf_quads,sizeof(int),1) ;
    fread(&num_vol_tets, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_tets,sizeof(int),1) ;
    fread(&num_vol_pents5, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_pents5,sizeof(int),1) ;
    fread(&num_vol_pents6, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_pents6,sizeof(int),1) ;
    fread(&num_vol_hexs, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_hexs,sizeof(int),1) ;
  }

  cout << " Number of  nodes = " << num_nodes << endl ;
  cout << " Number of surface triangles = " << num_sf_trias << endl ;
  cout << " Number of surface quads = " << num_sf_quads << endl ;
  cout << " Number of volume tetrahedras = " << num_vol_tets << endl ;
  cout << " Number of volume pents_5 = " << num_vol_pents5 << endl ;
  cout << " Number of volume pents6 = " << num_vol_pents6 << endl ;
  cout << " Number of volume hexahedra = " << num_vol_hexs << endl ;

  int num_quad_faces = num_sf_quads + num_vol_pents5 + num_vol_pents6*3
    + num_vol_hexs*6 ;
  if((num_quad_faces & 1) == 1) {
    cerr << "not all quad faces can pair!" << endl ;
    exit(-1) ;
  }
  // We counted face pairs, this is twice the final number
  num_quad_faces = num_quad_faces >> 1 ;
  int num_tri_faces = num_sf_trias + num_vol_tets*4+num_vol_pents5*4+num_vol_pents6*2 ;
  if((num_tri_faces & 1) == 1) {
    cerr << "not all trianglular faces can pair!" << endl ;
    exit(-1) ;
  }
  num_tri_faces = num_tri_faces >> 1 ;
  
  int ndim = 3 ;
  int nzones = 1 ;
  int npatch = 0 ;
  int ncells = num_vol_tets + num_vol_pents5 + num_vol_pents6 + num_vol_hexs ;
  int nfaces = num_quad_faces+num_tri_faces ;
  int max_ppf = 3 ;
  int max_fpc = 4 ;
  if(num_vol_pents5>0 || num_vol_pents6>0 || num_vol_hexs>0)
    max_ppf = 4 ;
  if(num_vol_pents5>0 || num_vol_pents6>0)
    max_fpc = 5 ;
  if(num_vol_hexs>0)
    max_fpc = 6 ;
  
  char out_buf[512] ;
  sprintf(out_buf,"%s.xdr",filename) ;
  FILE *FP = fopen(out_buf, "w") ;
  if(FP == NULL) {
    cerr << "can't open " << out_buf <<  endl ;
    return(-1);
  }

  XDR xdr_handle ;
  xdrstdio_create(&xdr_handle, FP, XDR_ENCODE) ;

  xdr_int(&xdr_handle, &ndim) ;
  xdr_int(&xdr_handle, &nzones) ;
  xdr_int(&xdr_handle, &npatch) ;
  xdr_int(&xdr_handle, &num_nodes) ;
  xdr_int(&xdr_handle, &nfaces) ;
  xdr_int(&xdr_handle, &ncells) ;
  xdr_int(&xdr_handle, &max_ppf) ;
  xdr_int(&xdr_handle, &max_fpc) ;

  cout << "XDR grid contains " << num_nodes << " nodes, "
       << nfaces << " faces, and " << ncells << " cells" << endl ;

  cout << "copying node information..." << endl ;
  if(!binary) {
    double ptmp ;
    for(int i = 0; i < 3 * num_nodes; ++i) {
      fscanf(IFP, "%lf", &ptmp) ;
      xdr_double(&xdr_handle, &ptmp) ;
    }
  } else {
    double ptmp[3] ;
    if(reverse_byteorder)
      for(int i = 0; i < num_nodes; ++i) {
        fread(ptmp,sizeof(double),3,IFP) ;
        ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
        xdr_double(&xdr_handle,&ptmp[0]) ;
        xdr_double(&xdr_handle,&ptmp[1]) ;
        xdr_double(&xdr_handle,&ptmp[2]) ;
      }
    else
      for(int i = 0; i < num_nodes; ++i) {
        fread(ptmp,sizeof(double),3,IFP) ;
        xdr_double(&xdr_handle,&ptmp[0]) ;
        xdr_double(&xdr_handle,&ptmp[1]) ;
        xdr_double(&xdr_handle,&ptmp[2]) ;
      }
  }

  

  size_t tri_total = num_tri_faces*2 ;
  size_t quad_total = num_quad_faces*2 ;
  tri_info *tri_faces = new tri_info[tri_total] ;
  quad_info *quad_faces = new quad_info[quad_total] ;

  Array<int,5> *bnd_tri = new Array<int,5>[num_sf_trias] ;
  Array<int,6> *bnd_quad = new Array<int,6>[num_sf_quads] ;

  if(quad_faces == 0 || tri_faces == 0 || bnd_tri == 0 || bnd_quad == 0) {
      cerr << "unable to allocate space for faces" << endl ;
      exit(-1) ;
  }

  cout << "reading in boundary information..." << endl ;
  
  size_t tf = 0 ;// triangle and quad face pointers
  size_t qf = 0 ;

  // Read in boundary triangles
  for(int i = 0; i < num_sf_trias; ++i) {
    tri_info tmp_tria ;
    tmp_tria[3] = -100000 ;
    if(!binary)
      fscanf(IFP, "%d%d%d", &tmp_tria[0], &tmp_tria[1], &tmp_tria[2]) ;  
    else { 
        fread(&tmp_tria[0], sizeof(int), 3, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_tria[0],sizeof(int),3) ;

    }
    tri_faces[tf++] = tmp_tria  ;
  }
  // Read in boundary quads
  for(int i = 0; i < num_sf_quads; ++i) {
    quad_info tmp_quad ;
    tmp_quad[4] = -100000 ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]) ;  
    else {
        fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
    }
    quad_faces[qf++] = tmp_quad ;
  }

  // Read in boundary flags for surfaces triangles
  for(int i = 0; i < num_sf_trias; ++i) {	
    if(!binary)
      fscanf(IFP, "%d", &tri_faces[i][3]) ;
    else {
        fread(&tri_faces[i][3], sizeof(int), 1, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tri_faces[i][3],sizeof(int),1) ;
    }
    tri_faces[i][3] = -tri_faces[i][3] ;
  }

  // Read in boundary flags for surface quads
  for(int i = 0; i < num_sf_quads; ++i) {
    if(!binary)
      fscanf(IFP, "%d", &quad_faces[i][4]) ;
    else {
        fread(&quad_faces[i][4], sizeof(int), 1, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&quad_faces[i][4],sizeof(int),1) ;

    }
    quad_faces[i][4] = -quad_faces[i][4] ;
  }

  int cellnum = 1 ;

  cout << "reading volume elements..." << endl ;
  
  // Read in volume tets and convert to faces
  for(int i = 0; i < num_vol_tets; ++i) {
    int tmp_quad[4] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]) ;
    else {
        fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
    }

    tri_faces[tf][0] = tmp_quad[0] ;
    tri_faces[tf][1] = tmp_quad[1] ;
    tri_faces[tf][2] = tmp_quad[3] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_quad[1] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[3] ; 
    tri_faces[tf++][3] = cellnum ;
    
   
    tri_faces[tf][0] = tmp_quad[3] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[0] ;
    tri_faces[tf++][3] = cellnum ;
  
    tri_faces[tf][0] = tmp_quad[0] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[1] ;
    tri_faces[tf++][3] = cellnum ;

    cellnum++ ;
  }

  // read in volume pent 5's and convert to faces
  for(int i = 0; i < num_vol_pents5; ++i) {
    int tmp_pents5[5] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d", &tmp_pents5[0], &tmp_pents5[1], &tmp_pents5[2], &tmp_pents5[3], &tmp_pents5[4]) ;
    else {
        fread(&tmp_pents5[0], sizeof(int), 5, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_pents5[0],sizeof(int),5) ;
    }
    tri_faces[tf][0] = tmp_pents5[4] ;
    tri_faces[tf][1] = tmp_pents5[1] ;
    tri_faces[tf][2] = tmp_pents5[2] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[4] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[3] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[3] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[0] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[0] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[1] ;
    tri_faces[tf++][3] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents5[0] ;
    quad_faces[qf][1] = tmp_pents5[1] ;
    quad_faces[qf][2] = tmp_pents5[4] ;
    quad_faces[qf][3] = tmp_pents5[3] ;
    quad_faces[qf++][4] = cellnum ;
    cellnum++ ;
  }
  for(int i = 0; i < num_vol_pents6; ++i) {
    int tmp_pents6[6] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d%d", &tmp_pents6[0], &tmp_pents6[1], &tmp_pents6[2], &tmp_pents6[3], &tmp_pents6[4], &tmp_pents6[5]) ;
    else {
        fread(&tmp_pents6[0], sizeof(int), 6, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_pents6[0],sizeof(int),6) ;
    }
    
    tri_faces[tf][0] = tmp_pents6[3] ;
    tri_faces[tf][1] = tmp_pents6[4] ;
    tri_faces[tf][2] = tmp_pents6[5] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents6[0] ;
    tri_faces[tf][1] = tmp_pents6[2] ;
    tri_faces[tf][2] = tmp_pents6[1] ;
    tri_faces[tf++][3] = cellnum ;
   
    quad_faces[qf][0] = tmp_pents6[0] ;
    quad_faces[qf][1] = tmp_pents6[1] ;
    quad_faces[qf][2] = tmp_pents6[4] ;
    quad_faces[qf][3] = tmp_pents6[3] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[1] ;
    quad_faces[qf][1] = tmp_pents6[2] ;
    quad_faces[qf][2] = tmp_pents6[5] ;
    quad_faces[qf][3] = tmp_pents6[4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[3] ;
    quad_faces[qf][1] = tmp_pents6[5] ;
    quad_faces[qf][2] = tmp_pents6[2] ;
    quad_faces[qf][3] = tmp_pents6[0] ;
    quad_faces[qf++][4] = cellnum ;

    cellnum++ ;
  }

  for(int i = 0; i < num_vol_hexs; ++i) {
    int tmp_hexs[8] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d%d%d%d", &tmp_hexs[0], &tmp_hexs[1], &tmp_hexs[2], &tmp_hexs[3], &tmp_hexs[4], &tmp_hexs[5], &tmp_hexs[6], &tmp_hexs[7] ) ;
    else {
        fread(&tmp_hexs[0], sizeof(int), 8, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_hexs[0],sizeof(int),8) ;
    }
    
    quad_faces[qf][0] = tmp_hexs[0] ;
    quad_faces[qf][1] = tmp_hexs[1] ;
    quad_faces[qf][2] = tmp_hexs[5] ;
    quad_faces[qf][3] = tmp_hexs[4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[1] ;
    quad_faces[qf][1] = tmp_hexs[2] ;
    quad_faces[qf][2] = tmp_hexs[6] ;
    quad_faces[qf][3] = tmp_hexs[5] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[2] ;
    quad_faces[qf][1] = tmp_hexs[3] ;
    quad_faces[qf][2] = tmp_hexs[7] ;
    quad_faces[qf][3] = tmp_hexs[6] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[4] ;
    quad_faces[qf][1] = tmp_hexs[7] ;
    quad_faces[qf][2] = tmp_hexs[3] ;
    quad_faces[qf][3] = tmp_hexs[0] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[4] ;
    quad_faces[qf][1] = tmp_hexs[5] ;
    quad_faces[qf][2] = tmp_hexs[6] ;
    quad_faces[qf][3] = tmp_hexs[7] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[3] ;
    quad_faces[qf][1] = tmp_hexs[2] ;
    quad_faces[qf][2] = tmp_hexs[1] ;
    quad_faces[qf][3] = tmp_hexs[0] ;
    quad_faces[qf++][4] = cellnum ;
    
    cellnum++ ;
  }

  cout << "finished reading ugrid file, matching faces..." << endl ;
  fclose(IFP) ;

  if((qf != size_t(num_quad_faces)*2) || tf != size_t(num_tri_faces)*2) {
    cerr << "face numbers not consistent!" ;
    exit(-1) ;
  }

  // prepare triangle faces (sort them)
  for(int i=0;i<num_tri_faces*2;++i) {
    // xdr numbers nodes from zero
    tri_faces[i][0] -= 1 ;
    tri_faces[i][1] -= 1 ;
    tri_faces[i][2] -= 1 ;
    
    if(tri_faces[i][0] > tri_faces[i][1])
      std::swap(tri_faces[i][0],tri_faces[i][1]) ;
    if(tri_faces[i][0] > tri_faces[i][2])
      std::swap(tri_faces[i][0],tri_faces[i][2]) ;
    if(tri_faces[i][1] > tri_faces[i][2])
      std::swap(tri_faces[i][1],tri_faces[i][2]
) ;
  }

  // prepare quad faces (sort them, but be careful)
  for(int i=0;i<num_quad_faces*2;++i) {
    // xdr numbers nodes from zero
    quad_faces[i][0] -=1 ;
    quad_faces[i][1] -=1 ;
    quad_faces[i][2] -=1 ;
    quad_faces[i][3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad_faces[i][0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad_faces[i][j]) {
        vs = quad_faces[i][j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad_faces[i][(j+nv)&0x3] ;
    // next make orientation so that it will match other face 
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[j] ;
    else
      for(size_t j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[(4 - j) &0x3 ] ;
  }


  cout << "sorting faces..." << endl ;

  std::sort(tri_faces,tri_faces+num_tri_faces*2,tri_sort) ;
  std::sort(quad_faces,quad_faces+num_quad_faces*2,quad_sort) ;


  cout << "writing face information..." << endl ;
  
  size_t btf = 0 ;
  size_t bqf = 0 ;
  
  int off = 0 ;

  for(int i=0;i<num_tri_faces;++i) {
    if(tri_faces[i*2][0] != tri_faces[i*2+1][0] ||
       tri_faces[i*2][1] != tri_faces[i*2+1][1] ||
       tri_faces[i*2][2] != tri_faces[i*2+1][2] ||
       (tri_faces[i*2][3] < 0 && tri_faces[i*2+1][3] < 0)) {
      cerr << "trouble matching triangle faces! " << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
    
    if(tri_faces[i*2][3] < 0) {
      bnd_tri[btf][0] = tri_faces[i*2][0] ;
      bnd_tri[btf][1] = tri_faces[i*2][1] ;
      bnd_tri[btf][2] = tri_faces[i*2][2] ;
      bnd_tri[btf][3] = tri_faces[i*2+1][3] ;
      bnd_tri[btf++][4] = tri_faces[i*2][3] ;
    } else {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &tri_faces[i*2+1][3]) ;
      xdr_int(&xdr_handle, &tri_faces[i*2][3]) ;
      off += 3 ;
    }
  }

  for(int i=0;i<num_quad_faces;++i) {
    if(quad_faces[i*2][0] != quad_faces[i*2+1][0] ||
       quad_faces[i*2][1] != quad_faces[i*2+1][1] ||
       quad_faces[i*2][2] != quad_faces[i*2+1][2] ||
       quad_faces[i*2][3] != quad_faces[i*2+1][3] ||
       (quad_faces[i*2][4] < 0 && quad_faces[i*2+1][4] < 0)) {
      cerr << "trouble matching quad faces!" << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
    
    if(quad_faces[i*2][4] < 0) {
      bnd_quad[bqf][0] = quad_faces[i*2][0] ;
      bnd_quad[bqf][1] = quad_faces[i*2][1] ;
      bnd_quad[bqf][2] = quad_faces[i*2][2] ;
      bnd_quad[bqf][3] = quad_faces[i*2][3] ;
      bnd_quad[bqf][4] = quad_faces[i*2+1][4] ;
      bnd_quad[bqf++][5] = quad_faces[i*2][4] ;
    } else {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &quad_faces[i*2+1][4]) ;
      xdr_int(&xdr_handle, &quad_faces[i*2][4]) ;
      off += 4 ;
    }
  }

  std::sort(bnd_tri,bnd_tri+btf,bnd_tri_sort) ;
  std::sort(bnd_quad,bnd_quad+bqf,bnd_quad_sort) ;

  for(size_t i=0;i<btf;++i) {
    xdr_int(&xdr_handle, &off) ;
    xdr_int(&xdr_handle, &bnd_tri[i][3]) ;
    xdr_int(&xdr_handle, &bnd_tri[i][4]) ;
    off += 3 ;
  }

  for(size_t i=0;i<bqf;++i) {
    xdr_int(&xdr_handle, &off) ;
    xdr_int(&xdr_handle, &bnd_quad[i][4]) ;
    xdr_int(&xdr_handle, &bnd_quad[i][5]) ;
    off += 4 ;
  }

  
  xdr_int(&xdr_handle, &off) ;

  for(int i=0;i<num_tri_faces;++i) {
    if(tri_faces[i*2][3] >=0) {
      xdr_int(&xdr_handle, &tri_faces[i*2][0]) ;
      xdr_int(&xdr_handle, &tri_faces[i*2][1]) ;
      xdr_int(&xdr_handle, &tri_faces[i*2][2]) ;
    }
  }

  for(int i=0;i<num_quad_faces;++i) {
    if(quad_faces[i*2][4] >=0) {
      xdr_int(&xdr_handle, &quad_faces[i*2][0]) ;
      xdr_int(&xdr_handle, &quad_faces[i*2][1]) ;
      xdr_int(&xdr_handle, &quad_faces[i*2][2]) ;
      xdr_int(&xdr_handle, &quad_faces[i*2][3]) ;
    }
  }

  for(size_t i=0;i<btf;++i) {
    xdr_int(&xdr_handle, &bnd_tri[i][0]) ;
    xdr_int(&xdr_handle, &bnd_tri[i][1]) ;
    xdr_int(&xdr_handle, &bnd_tri[i][2]) ;
  }

  for(size_t i=0;i<bqf;++i) {
    xdr_int(&xdr_handle, &bnd_quad[i][0]) ;
    xdr_int(&xdr_handle, &bnd_quad[i][1]) ;
    xdr_int(&xdr_handle, &bnd_quad[i][2]) ;
    xdr_int(&xdr_handle, &bnd_quad[i][3]) ;
  }

  xdr_destroy(&xdr_handle) ;
  fclose(FP) ;
  return 0 ; 
}
