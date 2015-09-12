//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#ifndef GMULTIMAP_H_
#define GMULTIMAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <Tools/debug.h>
#include <gmap.h>
#include <gkeyspace.h>
namespace Loci {

  class gMap ;
  class gMultiMap ;
  
  class gMultiMapRepI : public gMapRep {
    bool sorted;
    std::vector<std::pair<gEntity,gEntity> > attrib_data;
    gEntitySet dom;
    gKeySpace* domain_space;
    gKeySpace* image_space;
  public:
    //convenient functions for vertor operations
    typedef  std::vector<std::pair<gEntity,gEntity> >::iterator iterator ;
    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    iterator begin() { return attrib_data.begin(); }
    iterator end() { return attrib_data.end(); }
    const_iterator begin()const { return attrib_data.begin(); }
    const_iterator end() const { return attrib_data.end(); }
    size_t size() const{return attrib_data.size();}
    void reserve (size_t n){attrib_data.reserve(n);}
    void clear(){attrib_data.clear();}
    void local_sort();
          
    gMultiMapRepI():sorted(true),dom(GEMPTY),domain_space(0),image_space(0) {}
    void set_domain_space(gKeySpace* space){domain_space = space;}
    gKeySpace* get_domain_space()const{return domain_space;}
    void set_image_space(gKeySpace* space){image_space = space;}
    gKeySpace* get_image_space()const{return image_space;}
   
    virtual ~gMultiMapRepI(){} 
    virtual gStoreRepP clone() const{return new gMultiMapRepI(*this); }  
    virtual gStoreRepP remap(const gMap &m) const ;
    virtual void inplace_remap(const gMap &m); 
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const ;

    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const ;
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_ptn,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const ;
   
    virtual void inplace_compose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD ) ;
    virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD ) ;

    virtual gStoreRepP local_inverse() const;
    virtual gStoreRepP distributed_inverse(const std::vector<gEntitySet> &init_ptn) const;
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;
    virtual void unpack(const void *ptr, int &loc, int size) ;
  
    virtual void insert(gEntity e, gEntity val){
      sorted = false;
      attrib_data.push_back(std::pair<gEntity, gEntity>(e, val));
      dom += e;
    }

    virtual void insert(const gEntitySet& seq,  const gEntity* vals){
      sorted = false;
      size_t idx= 0;
      for(gEntitySet::const_iterator itr = seq.begin();
          itr!= seq.end(); itr++){
        attrib_data.push_back(std::pair<gEntity, gEntity>(*itr, vals[idx++]));
      }
      dom += seq;
    }
  
    virtual gstore_type RepType() const {return GMULTIMAP;}
    virtual gEntitySet domain() const {return dom ;} 
    virtual bool isSorted() const {return sorted;}
    virtual gEntitySet image(const gEntitySet &domain) const ;
    virtual gEntitySet image( gEntity domain) const ;
    virtual gEntitySet image() const ;
    
    virtual std::pair<gEntitySet,gEntitySet>
    preimage(const gEntitySet &codomain) const ;
    virtual gStoreRepP get_map()const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    //  virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, gEntitySet &user_eset) ;
    //     virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, gEntitySet& en) const ;
    virtual gStoreRepP expand(gEntitySet &out_of_dom, std::vector<gEntitySet> &init_ptn,MPI_Comm comm=MPI_COMM_WORLD)const ;
  
    virtual void* get_attrib_data() {return &attrib_data; }
    virtual const void* get_attrib_data() const{return &attrib_data; }
  
    virtual DatatypeP getType()const ;
  
  } ;


  
  class gMultiMap : public gstore_instance {
    friend class const_gMultiMap ;
    typedef gMultiMapRepI MapType ;
    typedef std::vector<std::pair<gEntity,gEntity> > gRep ; 
  public:

    typedef  std::vector<std::pair<gEntity,gEntity> >::iterator iterator ;
    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    iterator begin()  { return static_cast<gRep*>(Rep()->get_attrib_data())->begin(); }
    iterator end(){ return static_cast<gRep*>(Rep()->get_attrib_data())->end(); }
    const_iterator begin()const  { return static_cast<const gRep*>(Rep()->get_attrib_data())->begin(); }
    const_iterator end() const  { return static_cast<const gRep*>(Rep()->get_attrib_data())->end(); }
    size_t size() const{return Rep()->size();}
    void reserve (size_t n){Rep()->reserve(n);}
    void clear(){Rep()->clear();}
    void local_sort(){Rep()->local_sort();}
    
    
    
    // These should be private, but too much code currently depends on it
    // being public.  This code is dangerous because it does a shallow
    // copy which means that results sometimes may be unpredictable when
    // these operations are used.
    gMultiMap(const gMultiMap &var) { setRep(var.Rep()) ; }
    gMultiMap & operator=(const gMultiMap &str) { setRep(str.Rep()) ;  return *this ;}
    gMultiMap & operator=(gStoreRepP p) { setRep(p) ;  return *this ;}
    
    gMultiMap() {setRep(new MapType) ; }
    gMultiMap(gStoreRepP &rp) { setRep(rp) ; }
    virtual ~gMultiMap(){}

    void set_domain_space(gKeySpace* space){static_cast<gMapRepP>(Rep())->set_domain_space(space);}
    gKeySpace* get_domain_space()const{return static_cast<gMapRepP>(Rep())->get_domain_space() ;}
    void set_image_space(gKeySpace* space){static_cast<gMapRepP>(Rep())->set_image_space(space);}
    gKeySpace* get_image_space()const{return static_cast<gMapRepP>(Rep())->get_image_space();}
    
    virtual gStoreRepP clone() const{return Rep()->clone() ;}
        
    virtual  void insert(gEntity e, gEntity val){
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        p->insert(e, val);
      warn(p==0) ;}
    virtual void insert(const gEntitySet& seq,  const gEntity* vals){
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        p->insert(seq, vals);
      warn(p==0);
    }
    virtual instance_type access() const{return READ_WRITE ; }
    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);}
    virtual gStoreRepP local_inverse() const{return static_cast<gMapRepP>(Rep())->local_inverse();}
    virtual gStoreRepP distributed_inverse(const std::vector<gEntitySet> &init_ptn) const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->distributed_inverse(init_ptn);
      else return gStoreRepP(0);
    }
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm) ;}
    
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_ptn, comm) ;}

    virtual gStoreRepP expand(gEntitySet &out_of_dom, std::vector<gEntitySet> &init_ptn,MPI_Comm comm=MPI_COMM_WORLD)const{return  gMapRepP(Rep())->expand(out_of_dom, init_ptn, comm);}

    storeRepP copy2store() const;
    gEntitySet domain() const { return Rep()->domain() ; }
    operator gMapRepP() {
      gMapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    gEntitySet image(const gEntitySet &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    gEntitySet image(const gEntity &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    
    std::pair<gEntitySet,gEntitySet> preimage(const gEntitySet &codomain) const {
      return gMapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const gMultiMap &m)
  { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, gMultiMap &m)
  { return m.Input(s) ; }

  class const_gMultiMap : public gstore_instance {
    typedef gMultiMapRepI MapType ;
    typedef std::vector<std::pair<gEntity,gEntity> > gRep ; 
    const_gMultiMap(const const_gMultiMap &var) {setRep(var.Rep()) ;}
    const_gMultiMap(const gMultiMap &var) {setRep(var.Rep());}
    const_gMultiMap & operator=(const gMultiMap &str) { setRep(str.Rep()) ; return *this ;}
    const_gMultiMap & operator=(const const_gMultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:

    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    const_iterator begin()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->begin(); }
    const_iterator end()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->end(); }
    size_t size() const{return Rep()->size();}

    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm) ;}
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split,remap, comm) ;}
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_ptn, comm) ;}

    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);}
    virtual gstore_type RepType() const {return Rep()->RepType();}
    virtual ~const_gMultiMap(){}
    
    const_gMultiMap(){ setRep(new MapType);  }
    const_gMultiMap(gStoreRepP &rp) { setRep(rp) ; }
    
    virtual instance_type access() const{return READ_ONLY ; }
        
    const_gMultiMap & operator=(gStoreRepP p) { setRep(p) ; return *this ;}

    gEntitySet domain() const { return Rep()->domain(); }
    
    gEntitySet image(const gEntitySet &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    gEntitySet image(const gEntity &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    std::pair<gEntitySet,gEntitySet> preimage(const gEntitySet &codomain) const {
      return gMapRepP(Rep())->preimage(codomain) ;
    }
    operator gMapRepP()
    { gMapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    virtual gStoreRepP distributed_inverse(const std::vector<gEntitySet> &init_ptn) const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->distributed_inverse(init_ptn);
      else return gStoreRepP(0);
    }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const_gMultiMap &m)
  { return m.Print(s) ; }
    
  // void inverseMap(gMultiMap &result, const gMap &input_map,
  //                 const gEntitySet &input_image,
  //                 const gEntitySet &input_preimage);
  // void inverseMap(gMultiMap &result, const gMultiMap &input_map,
  //                 const gEntitySet &input_image,
  //                 const gEntitySet &input_preimage);
}

#endif
