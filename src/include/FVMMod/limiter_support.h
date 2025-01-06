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
  * @brief Venkatakrishnan Limiter Function, 
  * written in standard form, see @cite nishikawa2022new for reference
  * 
  * @param Xcc      [-] Cell center value
  * @param qmin     [-] min value in local stencil
  * @param qmax     [-] max value in local stencil
  * @param qdif     [-] face - cell center
  * @param eps2     [-] cell volume function
  * @return real    [-] limiter
  */
inline real vlimit(real Xcc, real qmin, real qmax, real qdif, real eps2) 
  {
    /////////////////////////////////
    // delta +
    const real delp = (qdif>=0.0)?qmax-Xcc:qmin-Xcc;
    // delta -
    const real delm = (qdif > 0)?qdif+1e-30:qdif-1e-30;
    // numerator of limiter
    const real num = ((delp*delp+eps2)*delm+ 2.0*delm*delm*delp)  ;
    // denominator of limiter
    const real den = (delm*(delp*delp+2.0*delm*delm+delm*delp+eps2)) ;
    // make the limiting case of 0/0 work as expected
    //const real e = (den >= 0.0?1.0e-30:-1.0e-30) ;
    real lim = num/den; //(num+e)/(den+e);
    return max<real>(0.0,min<real>(1.0,lim)) ;
  }

  /**
  * @brief Barth Limiter function
  * 
  * @param Xcc       [-] Cell center value
  * @param qdif      [-] face - cell center value
  * @param qmax      [-] maximum value in local stencil
  * @param qmin      [-] minimum value in local stencil
  * @return real     [-] limiter value
  */
inline real barth_limit(real Xcc, real qdif, real qmax, real qmin)
  {
    real lim;
    if (qdif > 0){
      lim = (qmax-Xcc)/(qdif+1e-30);
    } else {
      lim = (qmin-Xcc)/(qdif-1e-30);
    }
    return max<real>(0.0,min<real>(1.0,lim));
  }

  /**
  * @brief Nishikawa limiter function. See @cite nishikawa2022new for reference
   * 
  * @param Xcc      [-] Cell center value
  * @param qmin     [-] min value in local stencil
  * @param qmax     [-] max value in local stencil
  * @param qdif     [-] face - cell center
  * @param epsp     [-] cell volume function
  * @param nisPow   [-] order for Nishikawa Limiter
  * @return real    [-] limiter
  */
inline real nis_limit(real Xcc, real qmin, real qmax, real qdif, real epsp, int nisPow) {
    // delta +
    const real delp = (qdif>=0.0)?qmax-Xcc:qmin-Xcc;
    // delta -
    const real delm = qdif ;
    const real a = abs(delp);
    const real b = abs(delm) + 1e-30;
    if (a > 2*b) return 1.0; // save on computation 
    // numerator of limiter
    const real fun1 = pow(a,nisPow) + epsp;
    real Sp = 2.0*b*b; // nisPow = 2
    if (nisPow == 3)
    {
      Sp = 4.0*b*b;
    } else if (nisPow == 4){
      Sp = 2.0*b*(a*a - 2.0*b*(a-2*b));
    } else if (nisPow == 5) {
      Sp = 8.0*b*b*(a*a - 2.0*b*(a-b));
    }
    const real num = fun1 + a*Sp ;
    // denominator of limiter
    const real den =  fun1 + b*(pow(delp,nisPow-1) + Sp);
    real lim = num/den;
    return max<real>(0.0,lim) ;
  }

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
