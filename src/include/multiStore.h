#ifndef MULTISTORE_H
#define MULTISTORE_H

#include <istream>
#include <ostream>

extern "C" {
#include <hdf5.h>
}

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>

#include <Map.h>
#include <multiMap.h>
#include <DMultiStore.h>

namespace Loci {

  template<class T> class multiStoreRepI : public storeRep {
    entitySet store_domain ;
    T **index ;
    T *alloc_pointer ;
    T **base_ptr ;
    int size ;
    lmutex mutex ;
    bool istat ;

#ifdef ALLOW_DEFAULT_CONVERTER
    void  hdf5read( hid_t group_id, DEFAULT_CONVERTER c,      entitySet &en,
                    entitySet &usr);
    void  hdf5write( hid_t group_id, DEFAULT_CONVERTER c,      const entitySet &en) const;
    int   get_mpi_size( DEFAULT_CONVERTER c, const entitySet &eset);
    void  packdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(DEFAULT_CONVERTER c,      void *ptr, int &loc, int size,
                     const sequence &seq );
    void  StringVal( const int &entity, const int &ivec, std::string &memento);

#endif
    void  hdf5read( hid_t group_id, IDENTITY_CONVERTER c,     entitySet &en,
                    entitySet &usr);
    void  hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, entitySet &en,
                    entitySet &usr);
    void  hdf5write( hid_t group_id, IDENTITY_CONVERTER c,     const entitySet &en) const;
    void  hdf5write( hid_t group_id, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);

    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                    const sequence &seq );
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq);

  public:
    multiStoreRepI()
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; istat = 1; }

    multiStoreRepI(const entitySet &p)
    { index = 0; alloc_pointer = 0 ; base_ptr = 0 ; size=0; store_domain=p; istat = 1;}

    multiStoreRepI(const store<int> &sizes) {
      index = 0 ; alloc_pointer=0 ; base_ptr = 0; allocate(sizes) ;istat = 1; }

    void allocate(const store<int> &sizes) ;
    void multialloc(const store<int> &count, T ***index, T **alloc_pointer, T ***base_ptr) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~multiStoreRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( hid_t group_id, entitySet &en) ;
    virtual void writehdf5( hid_t group_id, entitySet& en) const ;

    bool is_static() { return istat ; } 
    T ** get_base_ptr() const { return base_ptr ; }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
  } ;
  
  //*************************************************************************/
  
  template<class T> class multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef Vect<T> containerType ;
    multiStore() {setRep(new storeType) ;}
    multiStore(multiStore<T> &var) {setRep(var.Rep()) ;}
    multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~multiStore() ;
    virtual void notification() ;
    
    multiStore<T> & operator=(multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    
    multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void setSizes(const const_multiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    const entitySet domain() const { return Rep()->domain() ; }

    Vect<T> elem(int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }

    Vect<T> operator[](int indx) 
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; 
    }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  //*************************************************************************/

  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const multiStore<T> &m)
  { return m.Print(s) ; }

  //************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, multiStore<T> &m)
  { return m.Input(s) ; }
 
  //************************************************************************/

  template<class T> 
  multiStore<T>::~multiStore() {}

  //************************************************************************/
  
  template<class T> 
  void multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //************************************************************************/
  
  template<class T> class const_multiStore : public store_instance {
    typedef multiStoreRepI<T> storeType ;
    T ** base_ptr ;
    int size ;
  public:
    typedef const_Vect<T> containerType ;
    const_multiStore() {setRep(new storeType) ;}
    const_multiStore(const_multiStore<T> &var) {setRep(var.Rep()) ;}
    const_multiStore(storeRepP rp) { setRep(rp) ;}
    
    virtual ~const_multiStore() ;
    virtual void notification() ;

    virtual instance_type access() const ;
    
    const_multiStore<T> & operator=(const multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(const const_multiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }

    const_multiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    const entitySet domain() const { return Rep()->domain() ; }

    containerType elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    containerType operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }

    const T *begin(int indx) const  { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  //*************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_multiStore<T>::access() const
  { return READ_ONLY ; }

  //*************************************************************************/
  
  template<class T> 
  const_multiStore<T>::~const_multiStore() {}

  //*************************************************************************/
  
  template<class T> 
  void const_multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::allocate(const store<int> &sizes) 
  {
    //-------------------------------------------------------------------------
    // Objective: Allocate memeory for multiStore data. This call reclaims 
    // all previously held memory
    //-------------------------------------------------------------------------
    // Assign new entitySet ...
    entitySet ptn = sizes.domain() ;
    store_domain  = ptn ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = 0 ;
    if(index) delete[] index ;

    index = 0 ;
    int sz = 0 ;
    if(ptn != EMPTY) {
      int top  = ptn.Min() ;
      int len  = ptn.Max() - top + 2 ;
      index    = new T *[len] ;
      base_ptr = index - top ;

      FORALL(ptn,i) {
        sz += sizes[i] ;
      } ENDFORALL ;

      alloc_pointer = new T[sz+1] ;
      sz = 0 ;
      for(int ivl=0;ivl< ptn.num_intervals(); ++ivl) {
        int i       = ptn[ivl].first ;
        base_ptr[i] = alloc_pointer + sz ;
        while(i<=ptn[ivl].second) {
          sz += sizes[i] ;
          ++i ;
          base_ptr[i] = alloc_pointer + sz ;
        }
      }

    }
    dispatch_notify();
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::multialloc(const store<int> &count, T ***index, 
                                     T **alloc_pointer, T ***base_ptr ) {
    entitySet ptn = count.domain() ;
    int top = ptn.Min() ;
    int len = ptn.Max() - top + 2 ;
    T **new_index = new T *[len] ;
    T **new_base_ptr = new_index - top ;
    int sz = 0 ;
    
    FORALL(ptn, i) {
      sz += count[i] ;
    } ENDFORALL ;
    
    T *new_alloc_pointer = new T[sz + 1] ;
    sz = 0 ;
    
    for(int ivl = 0; ivl < ptn.num_intervals(); ++ivl) {
      int i = ptn[ivl].first ;
      new_base_ptr[i] = new_alloc_pointer + sz ;
      while(i <= ptn[ivl].second) {
	sz += count[i] ;
	++i ;
	new_base_ptr[i] = new_alloc_pointer + sz ;
      }
    }
    
    *index = new_index ;
    *alloc_pointer = new_alloc_pointer ;
    *base_ptr = new_base_ptr ;
    
  }

  //*************************************************************************/
   
  template<class T> 
  void multiStoreRepI<T>::setSizes(const const_multiMap &mm)
  {
    //------------------------------------------------------------------------
    // Objective : Set the degree of each entity specified by the map..
    //------------------------------------------------------------------------

    mutex.lock() ;

    if(alloc_pointer != 0 && base_ptr[store_domain.Min()] == base_ptr[store_domain.Max()]) {
      delete[] index ;
      delete[] alloc_pointer ;
      index = 0 ;
      alloc_pointer = 0 ;
    }

    if(alloc_pointer != 0) {
      entitySet map_set = mm.domain() & store_domain ;
      entitySet problem ;
      FORALL(map_set,i) {
        if((end(i)-begin(i))<(mm.end(i)-mm.begin(i)))
          problem += i ;
      } ENDFORALL ;

      if(problem != EMPTY) {
        std::cerr << "reallocation of multiStore required for entities"
                  << problem << endl
                  << "Currently this reallocation isn't implemented."
                  << endl ;
      }
    } else {
      store<int> sizes ;
      sizes.allocate(store_domain) ;
      FORALL(store_domain,i) {
        sizes[i] = 0 ;
      } ENDFORALL ;
      entitySet map_set = mm.domain() & store_domain ;
      FORALL(map_set,i) {
        sizes[i] = (mm.end(i) - mm.begin(i)) ;
      } ENDFORALL ;
      allocate(sizes) ;
    }
    mutex.unlock() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::allocate(const entitySet &ptn) 
  {
    //------------------------------------------------------------------------
    // Objective : allocate memory specified by the entitySet. Allocation
    // doesn't resize the memory, therefore reclaims previously held memory.
    //------------------------------------------------------------------------
    
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
    
    alloc_pointer = 0 ;
    index         = 0 ;
    base_ptr      = 0 ;

    store_domain  = ptn ;

    //-------------------------------------------------------------------------
    // Initialize degree of each entity to zero 
    //-------------------------------------------------------------------------

    store<int> count ;
    count.allocate(ptn) ;
    
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    
    allocate(count) ;

    //-------------------------------------------------------------------------
    // Notify all observers ...
    //-------------------------------------------------------------------------

    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  multiStoreRepI<T>::~multiStoreRepI() 
  {
    if(alloc_pointer) delete[] alloc_pointer ;
    if(index) delete[] index ;
  }

  //***************************************************************************/

  template<class T> 
  storeRep *multiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    store<int> count ;
    count.allocate(p) ;
    entitySet ent = p - domain() ;
    
    for(entitySet::const_iterator ei = p.begin(); ei != p.end(); ++ei)
      count[*ei] = base_ptr[*ei+1] - base_ptr[*ei] ;
    
    for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei)
      count[*ei] = 0 ;
    
    return new multiStoreRepI<T>(count) ;
  }

  //***************************************************************************/

  template<class T> 
  storeRepP multiStoreRepI<T>::remap(const Map &m) const 
  {

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    multiStore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  //***************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_multiStore<T> s(st) ;
    fatal(alloc_pointer == 0) ;
    fatal((context - domain()) != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    store<int> count ;
    count.allocate(domain()) ;
    
    FORALL(domain() - context, i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context, i) {
      count[i] = s.end(i) - s.begin(i) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //***************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::gather(const Map &m, storeRepP &st,
                                 const entitySet &context) 
  {
    store<int> count ;
    const_multiStore<T> s(st) ;
    count.allocate(domain()) ;
    
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[i] = s.end(m[i])-s.begin(m[i]) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;

    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    
    FORALL(domain()-context,i) {
      for(int j = 0; j < count[i]; ++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j = 0; j < count[i]; ++j)
        new_base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;

    if(alloc_pointer) delete[] alloc_pointer ;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //***************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::scatter(const Map &m, storeRepP &st,
                                  const entitySet &context) 
  {
    
    store<int> count;
    
    const_multiStore<T> s(st) ;
    count.allocate(domain());
    
    FORALL(domain()-m.image(context),i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[m[i]] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ;
    
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
    
    FORALL(domain() - m.image(context),i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[m[i]];++j) {
        new_base_ptr[m[i]][j] = s[i][j] ;
      }
    } ENDFORALL ;
    
    if(alloc_pointer) delete[] alloc_pointer;
    alloc_pointer = new_alloc_pointer;
    if(index) delete[] index ;
    index = new_index ;
    base_ptr = new_base_ptr ;
    
    dispatch_notify() ;
  }

  //***************************************************************************/

  template<class T> 
  store_type multiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //***************************************************************************/
  
  template<class T> 
  entitySet multiStoreRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //***************************************************************************/
  
  template<class T> 
  std::ostream &multiStoreRepI<T>::Print(std::ostream &s) const 
  {
    //-------------------------------------------------------------------------
    // Objective : Print the multiStore data in the output stream.
    //-------------------------------------------------------------------------

    s << '{' << domain() << endl ;

    //-------------------------------------------------------------------------
    // Write the size of each entity
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Write the data of each entity in the domain 
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      for(const T *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << ' ' ;
      s << std::endl ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // close the bracket ...
    //-------------------------------------------------------------------------

    s << '}' << std::endl ;

    return s ;
  }

  //***************************************************************************/

  template<class T> 
  std::istream &multiStoreRepI<T>::Input(std::istream &s) 
  {
    //-------------------------------------------------------------------------
    // Objective : Read the multiStore data from the input stream
    //-------------------------------------------------------------------------

    entitySet e ;
    char ch ;

    // Look for opening bracket ....
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    //-------------------------------------------------------------------------
    // Read the entitySet intervals ....
    //-------------------------------------------------------------------------
    
    s >> e ;

    //-------------------------------------------------------------------------
    // Read the size of each entity in the set ...
    //-------------------------------------------------------------------------
    
    store<int> sizes ;
    sizes.allocate(e) ;

    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // read the data
    //-------------------------------------------------------------------------
    
    allocate(sizes) ;
    FORALL(e,ii) {
      for(T *ip = begin(ii);ip!=end(ii);++ip)
        s >> *ip  ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Look for closing brackets
    //-------------------------------------------------------------------------
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //**************************************************************************/

  template <class T> 
  int multiStoreRepI<T>::pack_size(const entitySet &eset ) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet ecommon;
    ecommon = eset & domain();

    return get_mpi_size( traits_type, ecommon );
  }
  //**************************************************************************/

#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T>
  void multiStoreRepI<T>::StringVal( const int &entity, const int &ivec, std::string &memento)
  {
    std::ostringstream oss;

    oss << base_ptr[entity][ivec] << endl;

    memento = oss.str();
  }
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(DEFAULT_CONVERTER c, const entitySet &eset ) 
  {
    int size;
    std::ostringstream oss;

    FORALL(eset,i) {
      size = base_ptr[i+1] - base_ptr[i] ;
      for(int j=0; j<size; j++)
        oss << base_ptr[i][j] << " ";
      oss << endl;
    } ENDFORALL ;

    std::string memento = oss.str();
    return(memento.length() );
  }
#endif
  //**************************************************************************/

  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset ) 
  {

    int size = 0 ;
    FORALL(eset,i) {
      size += base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;

    size *= sizeof(T) ;
    size += eset.size() * sizeof(int) ;
    return(size) ;
  }

  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset ) 
  {

    int        arraySize =0;
    entitySet  :: const_iterator ci;
    std::vector<T> avec;
    int       numContainers = 0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      size = base_ptr[(*ci)+1] - base_ptr[*ci] ;
      numContainers++;
      for( int ivec = 0; ivec < avec.size(); ivec++){
        Memento<T> memento( base_ptr[*ci][ivec] );
        arraySize += memento.getSize();
        numContainers++;
      }
    }

    typedef data_schema_converter_traits<T> converter_traits; 
    return( arraySize*sizeof(typename converter_traits::memento_type) +
            numContainers*sizeof(int) );
  }
  //**************************************************************************/

  template <class T> 
  void multiStoreRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &eset ) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    packdata( traits_type, ptr, loc, size, eset);

  }
  //**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
  template <class T> 
  void multiStoreRepI<T>::packdata( DEFAULT_CONVERTER c, void *outbuf,
                                    int &position,  int outcount,
                                    const entitySet &eset ) 
  {
    int  size;
    int  bufSize;
    std::ostringstream oss;
    std::string memento;

    FORALL(eset,ii){
      size = base_ptr[ii+1] - base_ptr[ii] ;
      MPI_Pack( &size, 1, MPI_INT, outbuf, outcount, &position, MPI_COMM_WORLD) ;
      for( int ivec = 0; ivec < size; ivec++){
        StringVal( ii, ivec, memento );
        bufSize = memento.length();
        MPI_Pack( &bufSize, 1, MPI_INT, outbuf, outcount, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack( &memento[0], bufSize, MPI_BYTE, outbuf, outcount, 
                  &position, MPI_COMM_WORLD) ;
      }
    }ENDFORALL ;
  }
#endif
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf,
                                    int &position,  int outcount,
                                    const entitySet &eset ) 
  {

    int  incount;
    store<int> count ;

    count.allocate(eset) ;

    FORALL(eset,i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;

    FORALL(eset,i) {
      incount = sizeof(int);
      MPI_Pack( &count[i], incount, MPI_BYTE, outbuf, outcount, &position, 
                MPI_COMM_WORLD) ; 
      incount = count[i]*sizeof(T);
      MPI_Pack( &base_ptr[i][0], incount, MPI_BYTE, outbuf, outcount, &position, 
                MPI_COMM_WORLD) ;
    } ENDFORALL ;

  }
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                    int &position, int outcount,
                                    const entitySet &eset ) 
  {
    int vecsize, stateSize, maxStateSize;

    //------------------------------------------------------------------------
    // Get the maximum size of container 
    //------------------------------------------------------------------------
    FORALL(eset,i) {
      vecsize = base_ptr[i+1] - base_ptr[i] ;
      for( int j = 0; j < vecsize; j++){
        Memento<T> memento( base_ptr[i][j]);
        stateSize = memento.getSize();
        maxStateSize = max( maxStateSize, stateSize);
      }
    } ENDFORALL ;

    typedef  data_schema_converter_traits<T> converter_traits; 
    typename converter_traits::memento_type *inbuf;

    int typesize = sizeof(typename converter_traits::memento_type);
    inbuf = new typename converter_traits::memento_type[maxStateSize];

    int incount;
    FORALL(eset,i) {
      vecsize = base_ptr[i+1] - base_ptr[i] ;

      incount =  sizeof(int);
      MPI_Pack(&vecsize, incount, MPI_BYTE, outbuf, outcount, &position,
               MPI_COMM_WORLD);

      for( int j = 0; j < vecsize; j++){
        Memento<T> memento( base_ptr[i][j] );
        memento.getState(inbuf, &stateSize);

        incount =  sizeof(int);
        MPI_Pack(&stateSize, incount, MPI_BYTE, outbuf, outcount, &position,
                 MPI_COMM_WORLD);

        incount =  stateSize*typesize;
        MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position, 
                 MPI_COMM_WORLD) ;
      }
    }
  }

  delete [] inbuf;

}

//**************************************************************************/

template <class T> 
void multiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, 
                               const sequence &seq) 
{

  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_type;

  unpackdata( traits_type, ptr, loc, size, seq); 

}

//**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T>::unpackdata( DEFAULT_CONVERTER c, void *inbuf, int &position, 
                                    int insize, const sequence &seq) 
{
  sequence:: const_iterator ci;
  char *outbuf;
  int   outcount, offset, size;
  entitySet eset(seq);

  for( ci = seq.begin(); ci != seq.end(); ++ci) {
    MPI_Unpack( inbuf, insize, &position, &size, 1, MPI_INT, 
                MPI_COMM_WORLD) ;
    for( int ivec = 0; ivec < size; ivec++){
      MPI_Unpack( inbuf, insize, &position, &outcount, 1, MPI_INT, 
                  MPI_COMM_WORLD) ;
      outbuf   = new char[outcount];

      MPI_Unpack( inbuf, insize, &position, outbuf, outcount, MPI_BYTE, 
                  MPI_COMM_WORLD) ;

      std::istringstream iss(outbuf);
      iss >> base_ptr[*ci][ivec];
      delete [] outbuf;
    }
  }
}
#endif

//**************************************************************************/
  
template <class T> 
void multiStoreRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *ptr, int &loc, int size, 
                                    const sequence &seq) {
  if(base_ptr == 0) return ;

  store<int> count ;
  bool conflict = 0 ;
  int temploc = loc ;
  entitySet new_dom = domain() | entitySet(seq) ;
  entitySet ent = domain() - entitySet(seq);
  count.allocate(new_dom) ;
    
  for(entitySet::const_iterator ei = domain().begin(); ei != domain().end(); ++ei)
    count[*ei] = base_ptr[*ei+1] - base_ptr[*ei] ;
    
  for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
    MPI_Unpack(ptr, size, &loc, &count[*si],sizeof(int), MPI_BYTE,
               MPI_COMM_WORLD) ;
    if(count[*si] != (base_ptr[*si+1] - base_ptr[*si])) conflict = 1 ;
    loc += count[*si] * sizeof(T) ;
  }
    
  if(conflict) {
    T **new_index ;
    T *new_alloc_pointer ;
    T **new_base_ptr ; 
    multialloc(count, &new_index, &new_alloc_pointer, &new_base_ptr) ;
      
    for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei) {
      for(int j = 0 ; j < count[*ei]; ++j) 
        new_base_ptr[*ei][j] = base_ptr[*ei][j] ;
    }
      
    if(alloc_pointer) delete [] alloc_pointer ;
      
    alloc_pointer = new_alloc_pointer;
      
    if(index) delete[] index ;
      
    index = new_index ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  loc = temploc ;
    
  for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
    loc += sizeof(int) ;
    MPI_Unpack(ptr, size, &loc, &base_ptr[*si][0],count[*si]*sizeof(T), MPI_BYTE,
               MPI_COMM_WORLD) ;
  }
}
  
//**************************************************************************
template <class T> 
void multiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                    int &position, int insize,
                                    const sequence &seq) 
{

  sequence :: const_iterator ci;
/*

  //-------------------------------------------------------------------------
  // Get the sum of each object size and maximum size of object in the
  // container for allocation purpose
  //-------------------------------------------------------------------------

  int  stateSize, outcount, vecsize, count_unpack = 0;

  typedef data_schema_converter_traits<T> converter_traits;
  typename converter_traits::memento_type *outbuf;
  int typesize = sizeof(typename converter_traits::memento_type);

  T  newObj;
  for( ci = seq.begin(); ci != seq.end(); ++ci) {
    if( count_unpack < insize ) {
      outcount      = sizeof(int);
      count_unpack += outcount;
      MPI_Unpack( inbuf, insize, &position, &vecsize, outcount, 
                  MPI_BYTE, MPI_COMM_WORLD) ;

      attrib_data[*ci].resize(vecsize);
      for( int i = 0; i < vecsize; i++) {
        outcount      = sizeof(int);
        count_unpack += outcount;
        MPI_Unpack(inbuf, insize, &position, &stateSize, outcount, 
                   MPI_BYTE, MPI_COMM_WORLD) ;

        outbuf = new typename converter_traits::memento_type[stateSize];

        outcount      = stateSize*typesize;
        count_unpack += outcount;
        MPI_Unpack(inbuf, insize, &position, outbuf, outcount, 
                   MPI_BYTE, MPI_COMM_WORLD) ;

        Memento<T> memento( attrib_data[*ci][i] );
        attrib_data[*ci][i] = memento.setState( outbuf, stateSize);
        delete [] outbuf;
      } 
    } else
      cout << "Warning: Insufficient buffer for multStore unpack " << *ci << endl;
  }
*/

}
//**************************************************************************/
  
template<class T> 
void multiStoreRepI<T>::readhdf5( hid_t group_id, entitySet &user_eset) 
{
  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_type;

  entitySet eset;
  HDF5_ReadDomain(group_id, eset);

  entitySet dom = eset & user_eset ;

  //-------------------------------------------------------------------------
  // Size of each main container....
  //--------------------------------------------------------------------------
  hsize_t  dimension;

  hid_t vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  hid_t vDataset   = H5Dopen(group_id,"ContainerSize");
  hid_t vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims (vDataspace, &dimension, NULL);

  int *ibuf = new int[dimension];
  H5Dread(vDataset,vDatatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

  store<int> container;
  container.allocate( eset );

  size_t indx      = 0;
  entitySet::const_iterator ci;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[*ci] = ibuf[indx++];

  allocate( eset );
  allocate( container);
  hdf5read( group_id, traits_type, eset, dom );

  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);
  delete ibuf;
}

//**************************************************************************/

template<class T> 
void multiStoreRepI<T>::writehdf5(hid_t group_id, entitySet &usr_eset) const 
{

  typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
  schema_converter traits_output_type;

  entitySet eset = domain()&usr_eset ;

  if( eset.size() < 1) return;

  //write out the domain
  HDF5_WriteDomain(group_id, eset);

  hdf5write(group_id, traits_output_type, eset);

}

//**************************************************************************/
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, DEFAULT_CONVERTER c, 
                                    entitySet &eset, entitySet &user_eset )
{
  hsize_t  dimension;
  hid_t vDatatype  = H5Tcopy(H5T_NATIVE_CHAR);
  hid_t vDataset   = H5Dopen(group_id,"VariableData");
  hid_t vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  char *ibuf = new char[dimension];

  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

  entitySet :: const_iterator ci;

  std::istringstream iss( ibuf );

  int count;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      iss >> base_ptr[*ci][j];
  }
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  delete ibuf;
}
#endif
//**************************************************************************/

template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, IDENTITY_CONVERTER c, 
                                    entitySet &eset, entitySet &usr_eset )
{

  hsize_t dimension;
  hid_t   vDatatype, vDataset, vDataspace, mDataspace;
  size_t indx = 0, arraySize;
  int    rank = 1;

  entitySet::const_iterator ci;

  //-------------------------------------------------------------------------
  // Size of each main container....
  //--------------------------------------------------------------------------
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dopen(group_id,"ContainerSize");
  vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  int *ibuf = new int[dimension];
  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

  store<int> container;
  container.allocate( eset );

  indx  = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[*ci] = ibuf[indx++];

  delete [] ibuf;
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  //---------------------------------------------------------------------------
  // Calculate the offset of each entity in file ....
  //---------------------------------------------------------------------------
  store<int>        bucket;
  bucket.allocate( usr_eset );

  for( ci = usr_eset.begin(); ci != usr_eset.end(); ++ci) 
    bucket[*ci] = container[*ci];
  allocate( bucket );

  store<unsigned>   offset;
  offset.allocate( eset );

  arraySize = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    offset[*ci] = arraySize;
    arraySize  += container[*ci];
  }

  //---------------------------------------------------------------------------
  // Read the data now ....
  //---------------------------------------------------------------------------
  int num_intervals = usr_eset.num_intervals();
  interval *it = new interval[num_intervals];

  for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

  dimension  = arraySize;
  mDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataset   = H5Dopen(group_id,"VariableData");

  AbstractDatatype  *dtype;
  typedef data_schema_traits<T> traits_type;
  dtype = traits_type::instance();
  vDatatype = dtype->get_hdf5_type();

  T *data;

  hssize_t  start_mem[] = {0};  // determines the starting coordinates.
  hsize_t   stride[]    = {1};  // which elements are to be selected.
  hsize_t   block[]     = {1};  // size of element block;
  hssize_t  foffset[]   = {0};  // location (in file) where data is read.
  hsize_t   count[]     = {0};  // how many positions to select from the dataspace

  for( int k = 0; k < num_intervals; k++) {
    count[0] = 0;
    for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  container[i];

    data = new T[count[0]];

    foffset[0] = offset[it[k].first];

    H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
    H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
    H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, data);

    indx = 0;
    for( int i = it[k].first; i <= it[k].second; i++) {
      for( int j = 0; j < container[i]; j++) 
        base_ptr[i][j] = data[indx++];
    }

    delete[] data;
  }
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);
  H5Sclose(mDataspace);
  delete dtype;
}

//**************************************************************************

template <class T> 
void multiStoreRepI<T> :: hdf5read( hid_t group_id, USER_DEFINED_CONVERTER c, 
                                    entitySet &eset, entitySet &usr_eset)
{

  hsize_t dimension;
  size_t indx = 0, arraySize;
  int    rank = 1;
  hid_t   vDatatype, vDataset, vDataspace, mDataspace;


  entitySet::const_iterator ci;

  typedef data_schema_converter_traits<T> converter_traits; 

  //-------------------------------------------------------------------------
  // Size of each main container....
  //--------------------------------------------------------------------------
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dopen(group_id,"ContainerSize");
  vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  int *ibuf = new int[dimension];
  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);

  store<int> container;
  container.allocate( eset );

  indx      = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[*ci] = ibuf[indx++];

  delete [] ibuf;
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  //---------------------------------------------------------------------------
  // Size of each sub-container ....
  //---------------------------------------------------------------------------
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dopen(group_id,"SubContainerSize");
  vDataspace = H5Dget_space(vDataset);
  H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

  ibuf = new int[dimension];
  H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, ibuf);
  H5Tclose(vDatatype);
  H5Dclose(vDataset);
  H5Sclose(vDataspace);

  int maxBucketSize = *max_element( ibuf, ibuf + (int)dimension );

  //---------------------------------------------------------------------------
  // Calculate the offset of each entity in file ....
  //---------------------------------------------------------------------------
  store<unsigned>   offset;
  store<int>        bucket;
  dmultiStore<int>  subcontainer;

  offset.allocate( eset );

  bucket.allocate( usr_eset );
  for( ci = usr_eset.begin(); ci != usr_eset.end(); ++ci) 
    bucket[*ci] = container[*ci];
  allocate( bucket );

  arraySize = 0;
  int indx2 = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    offset[*ci] = arraySize;
    for( int i = 0; i < container[*ci]; i++)  {
      size = ibuf[indx2];
      arraySize  += size;
      subcontainer[*ci].push_back( size );
      indx2++;
    }
  }

  delete [] ibuf;
  //---------------------------------------------------------------------------
  // Read the data now ....
  //---------------------------------------------------------------------------
  /*
    int num_intervals = usr_eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

    typedef data_schema_converter_traits<T> converter_traits; 
    typedef converter_traits::memento_type dtype;
    dtype *data, *buf;

    typedef data_schema_traits<dtype> traits_type;
    AtomicType *atom = traits_type::instance();
    vDatatype = atom->get_hdf5_type();


    dimension  = arraySize;
    vDataset   = H5Dopen(group_id,"VariableData");
    vDataspace = H5Dget_space(vDataset);
    mDataspace = H5Dget_space(vDataset);
    H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

    hssize_t  start_mem[] = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    buf  = new typename converter_traits::memento_type[maxBucketSize];

    for( int k = 0; k < num_intervals; k++) {
    count[0] = 0;
    for( int i = it[k].first; i <= it[k].second; i++){
    for( int j = 0; j < subcontainer[i].size(); j++)
    count[0] +=  subcontainer[i][j];
    }

    data = new typename converter_traits::memento_type[count[0]];
    foffset[0] = offset[it[k].first];

    H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem,  stride, 
    count, block);
    H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, 
    count, block);
    H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, 
    data);
 
    indx = 0;
    int size;
    for( int i = it[k].first; i <= it[k].second; i++) {
    for( int j = 0; j < subcontainer[i].size(); j++) {
    Memento<T> memento( base_ptr[i][j] );
    size = subcontainer[i][j];
    for( int m = 0; m < size; m++)
    buf[m] = data[indx++];
    base_ptr[i][j] = memento.setState( buf, size );
    }
    }

    delete[] data;
    }

    delete[] buf;
    H5Tclose(vDatatype);
    H5Dclose(vDataset);
    H5Sclose(vDataspace);
    H5Sclose(mDataspace);

    delete atom;
    delete buf;
  */

}; 

//**************************************************************************

template <class T> 
void multiStoreRepI<T>:: hdf5write( hid_t group_id, USER_DEFINED_CONVERTER g, 
                                    const entitySet &eset)  const
{   
  cout << " COMING TO MULTISTORE USER DEFINED " << endl;
  
  int rank = 1;
  hsize_t  dimension;
  hid_t    vDataspace, vDataset, vDatatype;
  hid_t cparms     = H5Pcreate (H5P_DATASET_CREATE);
    
  //-----------------------------------------------------------------------------
  // Get the sum of each object size and maximum size of object in the 
  // container for allocation purpose
  //-----------------------------------------------------------------------------

  entitySet :: const_iterator ci;
  size_t       arraySize= 0;
  int          count, stateSize, maxBucketSize;
  std::vector<int>  bucketSize;    //  Because we don't know in advance the size

  int *storeSize = new int[eset.size()];

  size_t indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    storeSize[indx++] = count;

    for( int j = 0; j < count; j++) {
      Memento<T> memento( base_ptr[*ci][j] );
      stateSize  = memento.getSize();
      arraySize += stateSize;
      bucketSize.push_back( stateSize );
    }
  }
  /*
    maxBucketSize = *max_element( bucketSize.begin(), bucketSize.end() );
    
    dimension=  eset.size();
    
    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
    cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, storeSize);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Tclose( vDatatype );
    delete [] storeSize;

    //-----------------------------------------------------------------------------
    // Write the size of each bucket...
    //-----------------------------------------------------------------------------
    dimension    =  bucketSize.size();
    int  *bucket = new int[bucketSize.size()];

    for( int i=0; i< bucketSize.size();i++)
    bucket[i] = bucketSize[i];
    bucketSize.clear();

    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDatatype  = H5Tcopy(H5T_NATIVE_INT);
    vDataset   = H5Dcreate( group_id, "SubContainerSize", vDatatype, vDataspace,
    cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, bucket);

    H5Tclose(vDatatype);
    H5Dclose(vDataset);
    H5Sclose(vDataspace);
    
    typedef data_schema_converter_traits<T> converter_traits; 
    typedef converter_traits::memento_type dtype;

    dtype *data, *buf;

    data =  new typename converter_traits::memento_type[arraySize];
    buf  =  new typename converter_traits::memento_type[maxBucketSize];
    //-----------------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-----------------------------------------------------------------------------

    indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) {
    Memento<T> memento( base_ptr[*ci][j] );
    memento.getState(buf, &stateSize);
    for( int i = 0; i < stateSize; i++)
    data[indx++] =  buf[i];
    }
    }

    //-----------------------------------------------------------------------------
    // Write (variable) Data into HDF5 format
    //-----------------------------------------------------------------------------
    dimension =  arraySize;

    typedef data_schema_traits<dtype> traits_type;
    AtomicType *atom = traits_type::instance();
    vDatatype = atom->get_hdf5_type();

    vDataspace = H5Screate_simple(rank, &dimension, NULL);
    vDataset   = H5Dcreate( group_id, "VariableData", vDatatype, vDataspace,
    cparms);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    //-----------------------------------------------------------------------------
    // Clean up
    //-----------------------------------------------------------------------------
    H5Tclose(vDatatype);
    H5Dclose(vDataset);
    H5Sclose(vDataspace);

    delete atom;
    delete [] data;
    delete [] buf;
  */
  cout << " COMING TO HDF5 USER DEFINED " << endl;

};

//*************************************************************************
#ifdef ALLOW_DEFAULT_CONVERTER
template <class T> 
void multiStoreRepI<T>::hdf5write( hid_t group_id, DEFAULT_CONVERTER g, 
                                   const entitySet &eset) const
{

  entitySet :: const_iterator ci;
  hid_t     vDataset, vDataspace, vDatatype;
    
  //------------------------------------------------------------------------
  // Write the container size:
  //------------------------------------------------------------------------
  int rank = 1;
  hsize_t  dimension;
  hid_t cparms     = H5Pcreate (H5P_DATASET_CREATE);

  dimension=  eset.size();
    
  int  *container = new int[eset.size()];

  size_t indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci)
    container[indx++] = end(*ci) - begin(*ci);

  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, container);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Tclose( vDatatype );
  delete [] container;
    

  std::ostringstream oss;
  int count;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      oss << base_ptr[*ci][j];
  }
  
  std::string memento = oss.str();

  hsize_t arraySize  =  memento.length();

  dimension  =  arraySize+1;

  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5Tcopy(H5T_NATIVE_CHAR);
  vDataset   = H5Dcreate( group_id, "VariableData", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, memento.c_str());

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Tclose( vDatatype );

}
#endif
//*************************************************************************

template <class T> 
void multiStoreRepI<T>::hdf5write( hid_t group_id, IDENTITY_CONVERTER g, 
                                   const entitySet &eset) const
{

  hid_t  vDatatype, vDataset, vDataspace;

  entitySet :: const_iterator ci;

  //-----------------------------------------------------------------------------
  // Get the sum of each object size and maximum size of object in the 
  // container for allocation purpose
  //-----------------------------------------------------------------------------
  int     count, newsize;

  int  *container = new int[eset.size()];

  size_t indx = 0, arraySize = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    newsize    = end(*ci) - begin(*ci);
    arraySize  += newsize;
    container[indx++] = newsize;
  }
         
  //------------------------------------------------------------------------
  // Write the Size of each multiStore ....
  //------------------------------------------------------------------------
  int rank = 1;
  hsize_t  dimension;

  dimension =  eset.size();
  hid_t cparms     = H5Pcreate (H5P_DATASET_CREATE);

  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDatatype  = H5Tcopy(H5T_NATIVE_INT);
  vDataset   = H5Dcreate( group_id, "ContainerSize", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, container);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Tclose( vDatatype );

  delete [] container;

  //-----------------------------------------------------------------------------
  // Collect state data from each object and put into 1D array
  //-----------------------------------------------------------------------------

  T  *data;
  data =  new T[arraySize];

  indx = 0;
  for( ci = eset.begin(); ci != eset.end(); ++ci) {
    count  = end(*ci) - begin(*ci);
    for( int j = 0; j < count; j++) 
      data[indx++] = base_ptr[*ci][j];
  }

  //-----------------------------------------------------------------------------
  // Write (variable) Data into HDF5 format
  //-----------------------------------------------------------------------------
  AbstractDatatype  *dtype;
  typedef data_schema_traits<T> traits_type;
  dtype = traits_type::instance();
  vDatatype = dtype->get_hdf5_type();

  dimension        = arraySize;
  vDataspace = H5Screate_simple(rank, &dimension, NULL);
  vDataset   = H5Dcreate( group_id, "VariableData", vDatatype, vDataspace,
                          cparms);
  H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose( vDataset  );
  H5Sclose( vDataspace);
  H5Tclose( vDatatype );
  
  delete [] data;

}; 

}

#endif

