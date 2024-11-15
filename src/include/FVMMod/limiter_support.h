/** ****************************************************************************
 * @file      limiter_support.h
 * @authors   Ed Luke
 *            Raymond Fontenot (CFDRC)
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Limiter support classes and limiter function defintions.
 * @attention This file is a part of the Loci Framework, a free software.
 * You can redistribute it and/or modify it under the terms of the Lesser
 * GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * The Loci Framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
 ******************************************************************************/

#include <Loci.h>

using std::cerr ;
using std::endl ;

namespace Loci {
  typedef real_t real ;
  typedef vector3d<real_t> vect3d ;

/**
 * @brief Class for storing max/min/norm for scalars
 * 
 */
  class stenMaxMinNorm {
    public:
      /**
       * @brief Public defintions for stenMaxMinNorm class
       * @param max       [-] max value
       * @param min       [-] min value
       * @param norm      [-] norm
       */
        real max;
        real min;
        real norm;
  };

/**
 * @brief Class for storing max/min/norm for vect3d
 * 
 */
  class stenMaxMinNormv3d {
    public:
      /**
       * @brief Public defintions for stenMaxMinNormv3d class
       * @param max       [-] max value
       * @param min       [-] min value
       * @param norm      [-] norm
       */
        vect3d max;
        vect3d min;
        real norm;
  };

/**
 * @brief Template operation for Max and Min of a scalar
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNorm
 */
    template <class T> struct MaxMin {
    void operator()(T &r, const T &s) {
      r.max = max(r.max,s.max) ;
      r.min = min(r.min,s.min) ;
    }
  } ;

/**
 * @brief Template operation for summing Max and Min of a scalar
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNorm
 */
    template <class T> struct SumMaxMin {
    void operator()(T &r, const T &s) {
      r.max += s.max;
      r.min += s.min;
    }
  } ;

/**
 * @brief Template operation for Max and Min of a vect3d
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNormv3d
 */
  template <class T> struct MaxMinV3D {
    void operator()(T &r, const T &s) {
      r.max.x = max(r.max.x,s.max.x);
      r.max.y = max(r.max.y,s.max.y);
      r.max.z = max(r.max.z,s.max.z);
      r.min.x = min(r.min.x,s.min.x);
      r.min.y = min(r.min.y,s.min.y);
      r.min.z = min(r.min.z,s.min.z);
    }
  } ;

/**
 * @brief Template operation for summing Max and Min of a vector
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNorm
 */
  template <class T> struct SumMaxMinV3D {
    void operator()(T &r, const T &s) {
      r.max.x += s.max.x;
      r.max.y += s.max.y; 
      r.max.z += s.max.z;
      r.min.x += s.min.x;
      r.min.y += s.min.y;
      r.min.z += s.min.z;
    }
  } ;

/**
 * @brief Template operation for Max and Min of a vector
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNorm
 */
  template <class T> struct MaxMinv {
    void operator()(T &r, const T &s) {
      // do nothing, as max/min for this is handled explicitly
    }
  } ;

/**
 * @brief Template operation for summing Max and Min of a vect3d
 * @param  r    in/out T
 * @param  s    in T for comparison
 * @tparam T    stenMaxMinNormv3d
 */
    template <class T> struct SumMaxMinv {
    void operator()(T &r, const T &s) {
      // do nothing, as max/min for this is handled explicitly
    }
  } ;

/**
 * @brief Venkatakrishnan Limiter Function
 * 
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param eps2     [-] cell volume function
 * @return real    [-] limiter
 */
real vlimit(real Xcc, real qmin, real qmax, real qdif, real eps2);


/**
 * @brief Barth Limiter function
 * 
 * @param Xcc       [-] Cell center value
 * @param qdif      [-] face - cell center value
 * @param qmax      [-] maximum value in local stencil
 * @param qmin      [-] minimum value in local stencil
 * @return real     [-] limiter value
 */
real barth_limit(real Xcc, real qdif, real qmax, real qmin);

/**
 * @brief Nishikawa limiter function. See () for reference
 * 
 * @param Xcc      [-] Cell center value
 * @param qmin     [-] min value in local stencil
 * @param qmax     [-] max value in local stencil
 * @param qdif     [-] face - cell center
 * @param epsp     [-] cell volume function
 * @param nisPow   [-] order of the limiter
 * @return real    [-] limiter
 */
real nis_limit(real Xcc, real qmin, real qmax, real qdif, real epsp, real nisPow);

  /**
   * @brief Register for stenMaxMinNorm class
   * 
   * @tparam  stenMaxMinNorm
   */
   template<> struct data_schema_traits<stenMaxMinNorm> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(stenMaxMinNorm()) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNorm,max) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNorm,min) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNorm,norm) ;
      return DatatypeP(ct) ;
    }
  } ;

  /**
   * @brief Register for stenMaxMinNormv3d class
   * 
   * @tparam  stenMaxMinNormv3d
   */
   template<> struct data_schema_traits<stenMaxMinNormv3d> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(stenMaxMinNormv3d()) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,max) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,min) ;
      LOCI_INSERT_TYPE(ct,stenMaxMinNormv3d,norm) ;
      return DatatypeP(ct) ;
    }
  } ;
}
