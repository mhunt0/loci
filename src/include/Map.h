//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef LOCI_MAP_H_
#define LOCI_MAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <Tools/debug.h>
#include <Map_rep.h>

namespace Loci {

  // The integer version of image_section.
  entitySet image_section(const int *start, const int *end) ;

#ifdef ENTITY
  // The Entity version of image_section.
  entitySet image_section(const Entity *start, const Entity *end);
#endif
  
  class Map ;
  
  class MapRepI : public MapRep {
    entitySet store_domain ;
    Entity *base_ptr ;
  public:
    MapRepI() { base_ptr = 0 ; }
    MapRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~MapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    virtual void pack(void *ptr, int &loc,
                      int &size, const entitySet &e, const Map& remap) ;
    virtual void unpack(void *ptr, int &loc,
                        int &size, const sequence &seq, const dMap& remap) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name,  entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name,  entitySet& en, hid_t xfer_plist_id) const ;
#endif    

    Entity * get_base_ptr() const { Entity * p = 0 ; if(alloc_id>=0) p = ((Entity *)storeAllocateData[alloc_id].base_ptr) -storeAllocateData[alloc_id].base_offset ; return p ; }
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
#ifdef DYNAMICSCHEDULING
    virtual storeRepP freeze(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    virtual storeRepP thaw(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
#endif
  } ;
  
  class Map : public store_instance {
    friend class const_Map ;
    typedef MapRepI MapType ;
    Entity* base_ptr ;
    Map(const Map &var) { setRep(var.Rep()) ; }
    Map & operator=(const Map &str) { setRep(str.Rep()) ; return *this ;}
  public:
    Map() { setRep(new MapType) ;}
    Map(storeRepP rp) { setRep(rp) ; }

    virtual ~Map() ;

    virtual void notification() ;

    Map & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    Entity &elem(Entity indx) { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const Entity &const_elem(Entity indx)  const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    Entity &operator[](Entity indx) { return elem(indx); }
    const Entity &operator[](Entity indx) const { return const_elem(indx) ; }
    Entity &operator[](size_t indx) { return elem(indx); }
    const Entity &operator[](size_t indx) const { return const_elem(indx) ; }
    //        operator int*() { return base_ptr; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
    int getRangeKeySpace() const { return MapRepP(Rep())->getRangeKeySpace() ; }
    void setRangeKeySpace(int v) { MapRepP(Rep())->setRangeKeySpace(v) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const Map &m)
    { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, Map &m)
    { return m.Input(s) ; }

  class const_Map : public store_instance {
    typedef MapRepI MapType ;
    const Entity* base_ptr ;
    const_Map(const const_Map &var) {setRep(var.Rep()) ; }
    const_Map(const Map &var) {setRep(var.Rep()); }
    const_Map & operator=(const Map &str) { setRep(str.Rep()) ; return *this ;}
    const_Map & operator=(const const_Map &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    const_Map()
    { setRep(new MapType); }
    const_Map(storeRepP rp) { setRep(rp) ; }
    
    virtual ~const_Map() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_Map & operator=(storeRepP p) { setRep(p) ; return *this ;}

    entitySet domain() const { return Rep()->domain(); }
    operator MapRepP()
    { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    const Entity &const_elem(Entity indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return base_ptr[indx]; }
    const Entity &operator[](Entity indx) const { return const_elem(indx) ; }
    const Entity &operator[](size_t indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    int getRangeKeySpace() const { return MapRepP(Rep())->getRangeKeySpace() ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const const_Map &m)
  { return m.Print(s) ; }

  const size_t IMAGE_THRESHOLD = 4 ; 

  inline entitySet remapSet(entitySet s, const Map &m) {
    entitySet t ;
    std::vector<int> vec ;
    std::vector<int>::const_iterator vi ;
    FORALL(s,fc) {
      vec.push_back(m[fc]) ;
      //t += m[fc] ;
    } ENDFORALL ;
    std::sort(vec.begin(), vec.end()) ;
    for(vi = vec.begin(); vi != vec.end(); ++vi)
      t += *vi ;
    return t ;
  }

  inline void remapMap(Map &m, const Map &remap, Map &tmp) {
    FORALL(m.domain(),fc) {
      tmp[remap[fc]] = m[fc] ;
    } ENDFORALL ;
    FORALL(m.domain(),fc) {
      m[fc] = tmp[fc] ;
    } ENDFORALL ;
  }
}


#endif
