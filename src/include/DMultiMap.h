#ifndef DMULTIMAP_H
#define DMULTIMAP_H

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <hdf5CC/H5cpp.h>

#include <vector>
#include <hash_map.h>

#include <Map.h>
#include <store.h>

namespace Loci {
class dmultiMapRepI : public MapRep {
    entitySet store_domain ;
    hash_map<int, vector<int> >  attrib_data;
  public:
    dmultiMapRepI() { }
    dmultiMapRepI(const store<int> &sizes) { allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dmultiMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRepP remap(const Map &m) const ;
    virtual void compose(const Map &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context)  ;
    virtual void gather(const Map &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const Map &m, storeRepP &st,
                         const entitySet &context) ;

    virtual int pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual multiMap get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group, entitySet &user_eset) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;

    hash_map<int,vector<int> > *get_attrib_data() {return &attrib_data;}

  } ;
      
  //***************************************************************************

  class dmultiMap : public store_instance {
    friend class const_dmultiMap ;
    typedef dmultiMapRepI MapType ;
    hash_map<int, vector<int> >   *attrib_data;
  public:
    dmultiMap() { setRep(new MapType) ; }
        
    dmultiMap(const store<int> &sizes) { setRep( new MapType(sizes) ); }
        
    dmultiMap(const dmultiMap &var) { setRep(var.Rep()) ; }

    dmultiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~dmultiMap() ;

    virtual void notification() ;

    dmultiMap & operator=(const dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

    dmultiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; 
    }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
         MapRepP p(Rep()) ;
         fatal(p==0) ;
         return p ; 
    }

    vector<int> &elem(int indx) {
      return (*attrib_data)[indx]; 
    }

    const vector<int> &const_elem(int indx)  const 
    {
      hash_map<int, vector<int> > :: const_iterator   ci;
     
      ci = attrib_data->find(indx);
      if( ci != attrib_data->end())
          return( ci->second );
      
      cerr << "Error: out of range entity " << endl;
      exit(0);
      return ci->second ;
    }

    vector<int> &operator[](int indx) { return elem(indx); }

    const vector<int> &operator[](int indx) const 
    { return const_elem(indx) ; }

    int num_elems(int indx) const 
    {
      hash_map<int, vector<int> > :: const_iterator   ci;
      vector<int>     newVec;
     
      ci = attrib_data->find(indx);
      if( ci != attrib_data->end())
          return( (ci->second).size() );

      return(0);

    }

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }

    std::istream &Input(std::istream &s) 
    { return Rep()->Input(s) ; }
  } ;
  
  //***************************************************************************

  inline std::ostream & operator<<(std::ostream &s, const dmultiMap &m)
    { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, dmultiMap &m)
    { return m.Input(s) ; }

  //***************************************************************************

  class const_dmultiMap : public store_instance {
    typedef dmultiMapRepI      MapType ;
    hash_map<int,vector<int> > *attrib_data;
  public:
    const_dmultiMap() { setRep(new MapType) ; }
    
    const_dmultiMap(const_dmultiMap &var) {  setRep(var.Rep()) ; }
    
    const_dmultiMap(dmultiMap &var) { setRep(var.Rep()) ; }
    
    const_dmultiMap(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dmultiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_dmultiMap & operator=(const_dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dmultiMap & operator=(dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dmultiMap & operator=(storeRepP p) 
    { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }

    operator MapRepP() 
    {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; 
    }

    const vector<int> const_elem(int indx)  const {
          hash_map<int,vector<int> > :: const_iterator   ci;
          ci = attrib_data->find(indx);
          if( ci != attrib_data->end() )
              return( ci->second );

          vector<int>  newVec;
          return( newVec );
    }

    const vector<int> operator[](int indx) const 
    { return const_elem(indx) ; }

    int num_elems(int indx) const 
    {
       hash_map<int,vector<int> > :: const_iterator   ci;
       ci = attrib_data->find(indx);
       if( ci != attrib_data->end() )
           return( (ci->second).size() );
       return(0);
    }

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }
  } ;

  //***************************************************************************

  inline std::ostream & operator<< (std::ostream &s, 
                                    const const_dmultiMap &m)
  { return m.Print(s) ; }

  //***************************************************************************

  void inverseMap(multiMap &result,
                  const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;

  //***************************************************************************

  void inverseMap(multiMap &result,
                  const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;

  //***************************************************************************

  void inverseMap( const dmultiMap &in, dmultiMap &out);
}

#endif
