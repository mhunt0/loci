#ifndef LOCI_TYPES_H
#define LOCI_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <data_traits.h>

#include <Tools/options_list.h>

namespace Loci {

  //-----------STD pair-------------------------------//
  template<class T1,class T2> std::ostream &
    operator<<(std::ostream &s, const std::pair<T1,T2> &v) {
    s<<"["<<v.first<<","<<v.second<<"]";
    return s;
  }

  template<class T1,class T2> std::istream &
    operator>>(std::istream &s, std::pair<T1,T2> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a '[' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ',' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ']' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    return s;
  }

  //---------------------vector3d------------------//
  template <class T> 
    struct vector3d {
      T x,y,z ;
      vector3d() {} 
      vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
      vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
    {
      s >> v.x >> v.y >> v.z ;
      return s ;
    }

  template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
  }

  //  template <class T> inline T dot(const vector3d<T> &v1, const T ra2[]) {
  //    return v1.x*ra2[0] + v1.y*ra2[1] + v1.z*ra2[2] ;
  //  }

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
  }

  //  template <class T> inline T dot(const T ra1[], const vector3d<T> &v2) {
  //    return ra1[0]*v2.x + ra1[1]*v2.y + ra1[2]*v2.z ;
  //  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
    return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                       v1.z*ra2[0]-v1.x*ra2[2],
                       v1.x*ra2[1]-v1.y*ra2[0]) ;
  }
  template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
    return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                       ra1[2]*v2.x-ra1[0]*v2.z,
                       ra1[0]*v2.y-ra1[1]*v2.x) ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    target.z += val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    target.z -= val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(float r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(long double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template <class T>
    struct data_schema_traits< vector3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector3d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;

        LOCI_INSERT_TYPE(ct,vector3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };

  template<class T>  struct tensor3d : public vector3d<vector3d< T > > {
    tensor3d() {}
    tensor3d(vector3d<T> xx,vector3d<T> yy, vector3d<T> zz)
      : vector3d<vector3d< T> > (xx,yy,zz) {}
    tensor3d(const tensor3d &v) : vector3d<vector3d< T> >(v) {}
  } ;

  template<class T> inline vector3d<T> dot(const tensor3d<T> &t,
                                           const vector3d<T> &v) {
    return vector3d<T>(dot(t.x,v),dot(t.y,v),dot(t.z,v)) ;
  }

  template<class T> inline tensor3d<T> product(const tensor3d<T> &t1,
                                               const tensor3d<T> &t2) {
    tensor3d<T> temp ;
    temp.x.x = t1.x.x*t2.x.x+t1.x.y*t2.y.x+t1.x.z*t2.z.x ;
    temp.y.x = t1.y.x*t2.x.x+t1.y.y*t2.y.x+t1.y.z*t2.z.x ;
    temp.z.x = t1.z.x*t2.x.x+t1.z.y*t2.y.x+t1.z.z*t2.z.x ;

    temp.x.y = t1.x.x*t2.x.y+t1.x.y*t2.y.y+t1.x.z*t2.z.y ;
    temp.y.y = t1.y.x*t2.x.y+t1.y.y*t2.y.y+t1.y.z*t2.z.y ;
    temp.z.y = t1.z.x*t2.x.y+t1.z.y*t2.y.y+t1.z.z*t2.z.y ;

    temp.x.z = t1.x.x*t2.x.z+t1.x.y*t2.y.z+t1.x.z*t2.z.z ;
    temp.y.z = t1.y.x*t2.x.z+t1.y.y*t2.y.z+t1.y.z*t2.z.z ;
    temp.z.z = t1.z.x*t2.x.z+t1.z.y*t2.y.z+t1.z.z*t2.z.z ;

    return temp ;
  }

  template <class T>
    struct data_schema_traits< tensor3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector3d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;

        LOCI_INSERT_TYPE(ct,tensor3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };
  
    
  //---------------------vector2d------------------//
  template <class T> 
    struct vector2d {
      T x,y ;
      vector2d() {} 
      vector2d(T xx,T yy) : x(xx),y(yy) {}
      vector2d(const vector2d &v) {x=v.x;y=v.y;}
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector2d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector2d<T> &v)
    {
      s >> v.x >> v.y ;
      return s ;
    }

  template <class T> inline T dot(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y ;
  }

  template <class T> inline T dot(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[0] + v1.y*ra2[1] ;
  }

  template <class T> inline T norm(const vector2d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y) ;
  }

  //  template <class T> inline T dot(const T ra1[], const vector2d<T> &v2) {
  //    return ra1[0]*v2.x + ra1[1]*v2.y ;
  //  }

  template<class T> inline T cross(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.y-v1.y*v2.x ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[1]-v1.y*ra2[0] ;
  }

  //  template<class T> inline T cross(const T ra1[], const vector2d<T> &v2) {
  //    return ra1[0]*v2.y-ra1[1]*v2.x ;
  //  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator*=(vector2d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T> inline vector2d<T> &operator/=(vector2d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+=(vector2d<T> &target, const vector2d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator-=(vector2d<T> &target, const vector2d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x+v2.x,v1.y+v2.y) ;
  }

  template<class T> inline vector2d<T> operator-(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x-v2.x,v1.y-v2.y) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, float r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(float r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, float r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, double r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(double r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, double r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  template<class T> inline vector2d<T> operator*(const vector2d<T> &v1, long double r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T> inline vector2d<T> operator*(long double r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T> inline vector2d<T> operator/(const vector2d<T> &v1, long double r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  
  template <class T>
    struct  data_schema_traits< vector2d<T> > {
    public:
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector2d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,y) ;
        return DatatypeP(ct) ;
      }
    };
  
  
  
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
    for(size_t i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,size_t n> inline std::istream &
    operator>>(std::istream &s, Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }

  
  template <class T,size_t n> 
  class data_schema_traits< Array<T,n> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = n ;
      return new ArrayType(getLociType(T()),sizeof(Array<T,n>),1,&dim) ;
    }
  };

  template <> struct data_schema_traits<options_list> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<options_list> Converter_Type ;
  } ;

  // For allocating temporary arrays of small size
  const int tmp_array_internal_SIZE=25 ;
  template <class T> class tmp_array {
    int sz ;
    T data[tmp_array_internal_SIZE] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > tmp_array_internal_SIZE)
        p = new T[sz] ;
    }
    void free() {
      if(sz > tmp_array_internal_SIZE)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { return p[i] ; }
    T & operator[](int i) const { return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;
      
}

#endif
