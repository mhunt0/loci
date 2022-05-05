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
#ifndef STORE_DEF_H
#define STORE_DEF_H

// This header file contains the class definition of
// storeRepI, store, and const_store.
// Their corresponding template implementation is in
// store_impl.h
// This separation is necessary to resolve some dependency
// problems in the class hierarchy.
// The same design applies to all other container classes.
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <DMap.h>
#include <store_rep.h>
#include <data_traits.h>
#include <sstream>
#include <hdf5_readwrite.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>

namespace Loci {
  extern int MPI_processes;
  extern int MPI_rank ;

  template<class T> class storeRepI : public storeRep {
    T *alloc_pointer ;
    T *base_ptr ;
    entitySet store_domain ;
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &usr);
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
#ifdef H5_HAVE_PARALLEL
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &usr, hid_t xfer_plist_id);
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en, hid_t xfer_plist_id) const;
#endif
    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int  get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void packdata(IDENTITY_CONVERTER, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;

    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en) const;
#ifdef H5_HAVE_PARALLEL
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en, hid_t xfer_plist_id) const;
#endif    
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    int  get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    storeRepI() { alloc_pointer = 0 ; base_ptr = 0; }
    storeRepI(const entitySet &p) { alloc_pointer=0 ; allocate(p) ;}
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~storeRepI()  ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual storeRepP thaw(const entitySet& es) const ;
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

    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    virtual entitySet domain() const ;
    T * get_base_ptr() const { return base_ptr ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;

  template<class T> class store : public store_instance {
    typedef storeRepI<T> storeType ;
    T* base_ptr ;
    store(const store &var) { setRep(var.Rep()) ; }
    store<T> & operator=(const store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    store() { setRep(new storeType); }
    store(storeRepP rp) { setRep(rp) ; }
    virtual ~store() {}
    virtual void notification() ;
    store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    T &elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return base_ptr[indx]; }
    const T &elem(Entity indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return base_ptr[indx]; }
    T &operator[](Entity indx) { return elem(indx);}
    const T&operator[](Entity indx) const { return elem(indx);}
    T &operator[](size_t indx) { return elem(indx);}
    const T&operator[](size_t indx) const { return elem(indx);}
  } ;

  template<class T> class const_store : public store_instance {
    typedef storeRepI<T> storeType ;
    const T * restrict base_ptr ;
    const_store(const store<T> &var) { setRep(var.Rep()) ; }
    const_store(const const_store &var) { setRep(var.Rep()) ; }
    const_store<T> & operator=(const const_store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_store<T> & operator=(const store<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    const_store() { setRep(new storeType) ; }
    const_store(storeRepP rp) { setRep(rp) ; }
    virtual ~const_store() {}
    virtual void notification() ;
    virtual instance_type access() const  ;
    const_store<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
#ifdef BOUNDS_CHECK
    const T &elem(Entity indx) const {
      fatal(base_ptr==NULL);
      fatal(!Rep()->domain().inSet(indx)) ;
      return base_ptr[indx]; }
    const T& operator[](Entity indx) const { return elem(indx); }
    const T& operator[](size_t indx) const { return elem(indx); }
#else
    const T &restrict elem(Entity indx) const restrict {
      return base_ptr[indx]; }
    const T& restrict operator[](Entity indx) const restrict { return elem(indx); }
    const T& restrict operator[](size_t indx) const restrict { return elem(indx); }
#endif

  } ;

  
} // end of namespace Loci

#endif
