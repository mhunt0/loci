#ifndef PARAMETER_H
#define PARAMETER_H

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <Tools/stream.h>
#include <hdf5CC/H5cpp.h>
#include <hdf5_traits.h>
#include <hdf5_write_template.h>

namespace Loci {
    
  template<class T> class paramRepI : public storeRep {
    entitySet store_domain ;
    T param ;
  public:
    paramRepI() { store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX) ; }
    paramRepI(const entitySet &p) { store_domain = p ;}
    virtual void allocate(const entitySet &p)  ;
    virtual ~paramRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual store_type RepType() const ;
    virtual const entitySet &domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5( H5::Group group) ;
    virtual void writehdf5( H5::Group group,entitySet& en) const ;
    T * get_param() { return &param ; }
  } ;

  template<class T> void paramRepI<T>::allocate(const entitySet &p) {
    store_domain = p ;
    dispatch_notify();
  }
  template<class T> paramRepI<T>::~paramRepI<T>() {}
    
  template<class T>
    storeRep *paramRepI<T>::new_store(const entitySet &p) const {
    return new paramRepI<T>(p) ;
  }

  template<class T> store_type paramRepI<T>::RepType() const {
    return PARAMETER ;
  }

  template<class T> const entitySet &paramRepI<T>::domain() const {
    return store_domain ;
  }
        
  template<class T> std::ostream &paramRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    s << param << std::endl ;
    s << '}' << std::endl ;
    return s ;
  }


  template<class T> std::istream &paramRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      s.putback(ch) ;
      e = ~EMPTY ;
      allocate(e) ;
      s>>param ;
      return s ;
    }
        
    s >> e ;
    allocate(e) ;
        
    s >> param ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading parameter" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  template<class T> void paramRepI<T>::readhdf5( H5::Group group){
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    entitySet en=param_hdf5read(group,traits_output_type,param);
    allocate(en);
  }

  template<class T> void paramRepI<T>::writehdf5( H5::Group group,entitySet& en) const{
    //entitySet en=domain();
    typedef typename hdf5_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    //schema_converter::hdf5write(group,param,en);
    param_hdf5write(group,traits_output_type,param,en);
  }

  template<class T> class param : public store_instance {
    typedef paramRepI<T> paramType ;
    T * data ;
  public:
    typedef T containerType ;
    param() { setRep(new paramType) ; }
    param(param &var) { setRep(var.Rep()) ; }
    param(storeRepP &rp) { setRep(rp); }

    virtual ~param() ;

    param & operator=(param &p) {setRep(p.Rep()) ; return *this ; }

    param & operator=(storeRepP p) {setRep(p) ; return *this ; }
    param & operator=(const T &v) { *data = v ; return *this ; }

    virtual void notification() ;

    
    T * operator->() { return data ; }
    const T * operator->() const { return data ; }
    
    T * operator&() { return data ; }
    const T * operator &() const { return data ; }

    T &operator*() { return *data ; }
    const T &operator*() const { return *data ; }

    T &operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }

    
    void set_entitySet(const entitySet &ptn) {Rep()->allocate(ptn); }

    operator storeRepP() { return Rep() ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  template<class T> param<T>::~param() {}
    
  template<class T> void param<T>::notification()
    {  NPTR<paramType> p(Rep());
       if(p!=0)
         data = p->get_param() ;
       warn(p==0);
    }

  template<class T> inline std::ostream &
    operator<<(std::ostream &s, const param<T> &t)
    { return t.Print(s) ; }

  template<class T> inline std::istream &
    operator>>(std::istream &s, param<T> &t)
    { return t.Input(s) ; }

  template<class T> class const_param : public store_instance {
    typedef T containerType ;
    typedef paramRepI<T> paramType ;
    const T * data ;
  public:
    const_param() { setRep(new paramType) ; }
    const_param(const_param<T> &var) { setRep(var.Rep()) ; }
    const_param(param<T> &var) { setRep(var.Rep()) ; }
    const_param(storeRepP &rp) { setRep(rp); }
    
    virtual ~const_param() ;

    const_param & operator=(const_param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(param<T> &p)
    { setRep(p.Rep) ; return *this ;}
    const_param & operator=(storeRepP p)
    { setRep(p) ; return *this ;}

    virtual void notification() ;
    virtual instance_type access() const ;
        
    const T * operator->() const { return data ; }
    
    const T * operator &() const { return data ; }

    const T &operator*() const { return *data ; }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL) ;
      fatal(!Rep()->domain().inSet(indx)) ;
#endif
      return *data ;
    }
    
    operator storeRepP() { return Rep() ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  template<class T> const_param<T>::~const_param() {}

  template<class T> void const_param<T>::notification() 
    {  NPTR<paramType> p(Rep());
       if(p!=0)
         data = p->get_param() ;
       warn(p==0);
    }
    
    
  template<class T> store_instance::instance_type
    const_param<T>::access() const
    { return READ_ONLY; }
    
  template<class T> inline std::ostream &
    operator<<(std::ostream &s, const const_param<T> &t)
    { return t.Print(s) ; }


}

#endif
    
