#include <istream>
#include <ostream>
#include <iostream>

#include <Map.h>
#include <DMultiMap.h>
#include <multiMap.h>
#include <Tools/hash_map.h>
using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;

namespace Loci 
{

  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;

  storeRepP dmultiMapRepI::thaw() {
    return getRep() ;
  }
  storeRepP dmultiMapRepI::freeze() {
    multiMap static_map ;
    store<int> count ;
    entitySet dom = domain() ;
    count.allocate(dom) ;
    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei)
      count[*ei] = attrib_data[*ei].size() ;
    static_map.allocate(count) ;
    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
      int i=0 ;
      for(std::vector<int,malloc_alloc<int> >::const_iterator vi=attrib_data[*ei].begin();
          vi!=attrib_data[*ei].end();++vi,++i)
        static_map[*ei][i] = *vi ;
    }
    return static_map.Rep() ;
  }
  
  storeRepP dmultiMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int,malloc_alloc<int> >::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int,malloc_alloc<int> > > copy(MPI_processes), send_clone(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = out_of_dom & ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    dmultiMap *map_entities = new dmultiMap[MPI_processes] ;

    for(int i = 0; i < MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi)
        (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet dom =map_entities[i].domain() ;
      send_count[i] = 2 * dom.size() ;
      for(entitySet::const_iterator hi = dom.begin(); hi != dom.end(); ++hi)
	send_count[i] += map_entities[i][*hi].size() ; 
    }
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += send_count[i] ;
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet dom = map_entities[i].domain() ;
      for(entitySet::const_iterator miv = dom.begin(); miv != dom.end(); ++miv) {
	send_map[size_send] = *miv ;
	++size_send ;
        int sz = map_entities[i][*miv].size() ;
	send_map[size_send] = sz ;
	++size_send ;
        for(int k=0;k<sz;++k) {
	  send_map[size_send] = map_entities[i][*miv][k] ;
	  ++size_send ;
	}
      }
    }
    delete[] map_entities ;
    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    dmultiMap hm ;
    std::vector<int,malloc_alloc<int> > ss ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	int count = recv_map[j+1] ;
	if(count)
	  for(int k = 0; k < count; ++k)
	    hm[recv_map[j]].push_back(recv_map[j+k+2]);
	else
	  hm[recv_map[j]] = ss ;
	j += count + 1 ;
      }
    }
    std::vector<int,malloc_alloc<int> > tmp_vec ;
    entitySet dom = hm.domain() ;
    for(entitySet::const_iterator hmi = dom.begin(); hmi != dom.end(); ++hmi)
      attrib_data[*hmi].swap(hm[*hmi]) ;
    dmultiMap dmul ;
    entitySet hdom = attrib_data.domain() ;
    for(entitySet::const_iterator hi = hdom.begin(); hi != hdom.end(); ++hi) {
      size_t sz = attrib_data.elem(*hi).size() ;
      vector<int,malloc_alloc<int> >(sz).swap(dmul[*hi]) ;
      dmul[*hi] = attrib_data.elem(*hi) ;
    }
    storeRepP sp = dmul.Rep() ;
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ; 
    return sp ;
  }
  //------------------------------------------------------------------
  
  void dmultiMapRepI::allocate(const entitySet &eset) 
  {

    entitySet redundant, newSet;
    entitySet :: const_iterator  ci;

    redundant = domain() -  eset;
    newSet    = eset - domain();

    attrib_data.erase_set(redundant) ;

    for( ci = newSet.begin(); ci != newSet.end(); ++ci)
      std::vector<int,malloc_alloc<int> >(0).swap(attrib_data[*ci]) ;

    dispatch_notify() ;
  }
  
  //------------------------------------------------------------------
  
  void dmultiMapRepI::allocate(const store<int> &sizes) 
  {
    entitySet ptn = sizes.domain() ;
    entitySet :: const_iterator  ci;
    for( ci = ptn.begin(); ci != ptn.end(); ++ci) {
      std::vector<int,malloc_alloc<int> >   newVec(sizes[*ci]) ;
      attrib_data[*ci] = newVec;
    }
    dispatch_notify() ;
  }
  
  //**************************************************************************/

  dmultiMapRepI::~dmultiMapRepI() 
  {
  }
  
  //**************************************************************************/
  
  storeRep *dmultiMapRepI::new_store(const entitySet &p) const 
  {
    return new dmultiMapRepI()  ;
  }
  storeRep *dmultiMapRepI::new_store(const entitySet &p, const int* cnt) const 
  {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new dmultiMapRepI(count)  ;
  }
  //**************************************************************************/
  
  storeRepP dmultiMapRepI::remap(const dMap &m) const {
    dmultiMap s ;
    
    //-------------------------------------------------------------------------
    // Select only those entities from domain of "m" (which is input) which 
    // are part of multimap. Since, we have mapping defined only on those
    // entities.
    //-------------------------------------------------------------------------
    
    entitySet newdomain = m.domain() & domain() ;
    
    //-------------------------------------------------------------------------
    // Get the preimage of entire map. We will get two entity set. The
    // first is the intersection, and the second  is union of entities.
    //-------------------------------------------------------------------------
    pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    
    //
    // Get all entities which are both in preimage and newdomain
    //
    newdomain &= mappimage.first ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(m,mapimage) ;
    
    multiMap   newmap; 
    newmap = MapRepP(s.Rep())->get_map() ;
    return newmap.Rep() ;  
   
    // return s.Rep() ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::compose(const dMap &m, const entitySet &context) 
  {
    vector<int,malloc_alloc<int> >    vec;
    
    //-------------------------------------------------------------------------
    // All the entities in the context should be present in the domain. ie. A->B
    // should be valid.
    //-------------------------------------------------------------------------
    fatal((context-domain()) != EMPTY) ;
    
    //-------------------------------------------------------------------------
    // All in the entities in B should be part of Map. Also. i.e. B->C should
    // be valid.
    //-------------------------------------------------------------------------
    fatal((image(context)-m.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      
      for(size_t j = 0; j < attrib_data.elem(i).size(); j++) {
        attrib_data.elem(i)[j] =   m[attrib_data.elem(i)[j]];
      }
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiMap s(st) ;
    vector<int,malloc_alloc<int> >    newVec;
    
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      size_t sz = s[i].size() ;
      std::vector<int,malloc_alloc<int> >(sz).swap(attrib_data[i]) ;
      for(size_t j = 0; j < sz; j++)
        attrib_data[i][j] = s[i][j] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::gather(const dMap &m, storeRepP &st, const entitySet  &context) 
  {
    const_dmultiMap s(st) ;
    vector<int,malloc_alloc<int> >    newVec;
    
    FORALL(context,i) {
      size_t sz = s[m[i]].size() ;
      std::vector<int,malloc_alloc<int> >(sz).swap(attrib_data[i]) ;
      for(size_t j = 0; j < sz; j++) 
        attrib_data[i][j] = s[m[i]][j] ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::scatter(const dMap &m, storeRepP &st, const entitySet  &context) 
  {
    const_dmultiMap s(st) ;
    vector<int,malloc_alloc<int> >    newVec;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    FORALL(context,i) {
      size_t sz = s[i].size() ;
      std::vector<int,malloc_alloc<int> >(sz).swap(attrib_data[m[i]]) ;
      for(size_t j = 0; j < sz; j++) 
        attrib_data[m[i]][j] = s[i][j] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  int dmultiMapRepI::pack_size(const  entitySet &e ) 
  {
    int size = 0 ;
    FORALL(e,i) {
      size  +=  attrib_data[i].size();
    } ENDFORALL ;
    
    return( size*sizeof(int) + e.size()*sizeof(int) ) ;
  }
  
  //**************************************************************************/
  
  void dmultiMapRepI::pack( void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
     int vsize;
    entitySet :: const_iterator ci;
    std::vector<int,malloc_alloc<int> >   newVec;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize  = attrib_data[*ci].size();
      MPI_Pack( &vsize, 1, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
      MPI_Pack( &attrib_data[*ci][0], vsize, MPI_INT, outbuf,outcount,
                &position, MPI_COMM_WORLD) ;
    }

  }
  
  
  //**************************************************************************/

  void dmultiMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) 
  {
    sequence:: const_iterator ci;
    std::vector<int,malloc_alloc<int> >   newVec;

    int vsize;
    for( ci = seq.begin(); ci != seq.end(); ++ci){
         MPI_Unpack( inbuf, insize, &position, &vsize,
                     1, MPI_INT, MPI_COMM_WORLD) ;
         std::vector<int,malloc_alloc<int> >(vsize).swap(attrib_data[*ci]) ;

         if(vsize != 0)
           MPI_Unpack( inbuf, insize, &position, &attrib_data[*ci][0],
                       vsize, MPI_INT, MPI_COMM_WORLD) ;
    }

  }   
      
  //**************************************************************************/
    
  entitySet dmultiMapRepI::domain() const 
  {
    return attrib_data.domain() ;
  }

  //**************************************************************************/
  
  entitySet dmultiMapRepI::image(const entitySet &iset) const 
  {
    entitySet :: const_iterator  ei;

    entitySet dom = attrib_data.domain() & iset ;
    
    vector<Entity> codomlist ;
    for( ei = dom.begin(); ei != dom.end(); ++ei){
      const vector<int,malloc_alloc<int> > &ref = attrib_data.elem(*ei) ;
      size_t sz = ref.size() ;
      for(size_t i = 0; i < sz; i++)
        codomlist.push_back(ref[i]) ;
    }
    return create_intervalSet(codomlist.begin(),codomlist.end()) ;
  }

  //**************************************************************************/
 
  pair<entitySet,entitySet>
  dmultiMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domaini,domainu ;
    FORALL(domain(),i) {
      bool vali = true ;
      const vector<int,malloc_alloc<int> > &ref = attrib_data.elem(i) ;

      bool valu = false;
      for(std::vector<int,malloc_alloc<int> >::const_iterator vi = ref.begin();
          vi != ref.end(); ++vi) {
        bool in_set = codomain.inSet(*vi) ;
	vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return make_pair(domaini,domainu) ;
  }

  //**************************************************************************/

  storeRepP dmultiMapRepI::get_map() 
  {
    multiMap   newmap;
    entitySet::const_iterator  ei;
    store<int> count ;
    entitySet dom = domain() ;
    count.allocate(dom) ;
    for(ei = dom.begin(); ei != dom.end(); ++ei)
      count[*ei] = attrib_data[*ei].size() ;
    newmap.allocate(count) ;
    for(ei = dom.begin(); ei != dom.end(); ++ei) 
      for(int i = 0; i < count[*ei]; i++)
	newmap[*ei][i] = attrib_data[*ei][i];
    
    return newmap.Rep();
  }
  
  //**************************************************************************/
    
  ostream &dmultiMapRepI::Print(ostream &s) const 
  {

    s << '{' << domain() << endl ;

    FORALL(domain(),ii) {
      const vector<int,malloc_alloc<int> > &ref = attrib_data.elem(ii) ;

      s << ref.size() << endl;
    } ENDFORALL ;
    
    FORALL(domain(),ii) {
      const vector<int,malloc_alloc<int> > &ref = attrib_data.elem(ii) ;
      for(size_t i = 0; i < ref.size(); i++)
        s << ref[i] << "    ";
      s << endl;
    } ENDFORALL ;
    
    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/

  istream &dmultiMapRepI::Input(istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;
    store<int> sizes ;

    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    allocate(sizes) ;
        
    FORALL(e,ii) {
      std::vector<int,malloc_alloc<int> >(sizes[ii]).swap(attrib_data[ii]) ;
      for( int i = 0; i < sizes[ii]; i++)
        s >> attrib_data[ii][i] ;
    } ENDFORALL ;
            
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }
    return s ;
  }
  DatatypeP dmultiMapRepI::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  
  frame_info dmultiMapRepI::read_frame_info(hid_t group_id) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  frame_info dmultiMapRepI::write_frame_info(hid_t group_id) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  //**************************************************************************/

  void dmultiMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
  }

  //**************************************************************************/

  void dmultiMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const
  {
    warn(true) ;
  } 

  //**************************************************************************/
  
  dmultiMap::~dmultiMap() {}

  //**************************************************************************/

  void dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  const_dmultiMap::~const_dmultiMap() { }

  //**************************************************************************/

  void const_dmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  store_instance::instance_type const_dmultiMap::access() const
  { return READ_ONLY ; }

  //**************************************************************************/


  void inverseMap(dmultiMap &result, const dMap &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) {
    entitySet preloop = input_preimage & input_map.domain() ;

    std::vector<int,malloc_alloc<int> > tmp_vec ;
    FORALL(input_image,i) {
      result[i] = tmp_vec ;
    } ENDFORALL ;

    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
	result[elem].push_back(i) ;
    } ENDFORALL ;
  }
  void inverseMap(dmultiMap &result, const Map &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) {
    entitySet preloop = input_preimage & input_map.domain() ;
    std::vector<int,malloc_alloc<int> > tmp_vec ;
    FORALL(input_image,i) {
      result[i] = tmp_vec ;
    } ENDFORALL ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
	result[elem].push_back(i) ;
    } ENDFORALL ;
  }

  void inverseMap(dmultiMap &result, const const_dMap &input_map,
                   const entitySet &input_image,
 		  const entitySet &input_preimage) {
     entitySet preloop = input_preimage & input_map.domain() ;
     std::vector<int,malloc_alloc<int> > tmp_vec ;
     FORALL(input_image,i) {
       result[i] = tmp_vec ;
     } ENDFORALL ;
     FORALL(preloop,i) {
       int elem = input_map[i] ;
       if(input_image.inSet(elem)) 
 	result[elem].push_back(i) ;
     } ENDFORALL ;
   }
   void inverseMap(dmultiMap &result, const const_Map &input_map,
                   const entitySet &input_image,
 		  const entitySet &input_preimage) {
     entitySet preloop = input_preimage & input_map.domain() ;
     std::vector<int,malloc_alloc<int> > tmp_vec ;
     FORALL(input_image,i) {
       result[i] = tmp_vec ;
     } ENDFORALL ;
     FORALL(preloop,i) {
       int elem = input_map[i] ;
       if(input_image.inSet(elem)) 
 	result[elem].push_back(i) ;
     } ENDFORALL ;
   }
 
   void inverseMap(dmultiMap &result, const dmultiMap &input_map,
                   const entitySet &input_image,
                   const entitySet &input_preimage) {
     entitySet preloop = input_preimage & input_map.domain() ;
     std::vector<int,malloc_alloc<int> > tmp_vec ;
     FORALL(input_image,i) {
       result[i] = tmp_vec ;
     } ENDFORALL ;
     FORALL(preloop,i) {
       std::vector<int,malloc_alloc<int> >::const_iterator vi ;
       for(vi=input_map[i].begin();vi!=input_map[i].end();++vi) {
         int elem = *vi ;
         if(input_image.inSet(elem)) 
           result[elem].push_back(i) ;
       }
     } ENDFORALL ;
     
   }
   
   void inverseMap(dmultiMap &result, const multiMap &input_map,
                   const entitySet &input_image,
                   const entitySet &input_preimage) {
     entitySet preloop = input_preimage & input_map.domain() ;
     std::vector<int,malloc_alloc<int> > tmp_vec ;
     FORALL(input_image,i) {
       result[i] = tmp_vec ;
     } ENDFORALL ;
     FORALL(preloop,i) {
       for(const Entity *ep = input_map.begin(i);ep!=input_map.end(i);++ep) {
         int elem = *ep ;
         if(input_image.inSet(elem)) 
           result[elem].push_back(i) ;
       }
     } ENDFORALL ;
     
   }
}

