#ifndef OVERSET_H
#define OVERSET_H
#include "Loci.h"
#include <vector>
namespace Loci {
  typedef unsigned char byte_t ;
  struct Quaternion {
    real_t x,y,z,w ;
    Quaternion() {}
    Quaternion(real_t xi, real_t yi, real_t zi, real_t wi):x(xi),y(yi),z(zi),w(wi) {}
    Quaternion(vector3d<real_t> axis, real_t angle) {
      real_t sinAngle; 
      angle *= 0.5; 
      axis *= 1.0/(norm(axis)+1e-30) ; 
      sinAngle = sin(angle);
      x = (axis.x * sinAngle); 
      y = (axis.y * sinAngle); 
      z = (axis.z * sinAngle); 
      w = cos(angle);
    }
    Quaternion operator*(const Quaternion &q) const {
      vector3d<real_t> vector1(x,y,z), vector2(q.x,q.y,q.z); 

      const real_t angle = ((w * q.w) - (dot(vector1, vector2))); 
      const vector3d<real_t> across = cross(vector1, vector2);
      vector1 *= q.w ;
      vector2 *= w ;
      Quaternion result; 
      result.x = (vector1.x + vector2.x + across.x); 
      result.y = (vector1.y + vector2.y + across.y); 
      result.z = (vector1.z + vector2.z + across.z); 
      result.w = angle;
      return result ;
    }
    Quaternion &Normalize() {
      // reciprocal of the l2 norm 
      const real_t rl2 = 1.0 / sqrt((x*x) + (y*y) + (z*z) + (w*w)) ;
      x*=rl2 ;
      y*=rl2 ;
      z*=rl2 ;
      w*=rl2 ;
      return *this ;
    }
    
    Quaternion Inverse() const {
      Quaternion result = *this ;
      result.x *= -1 ; 
      result.y *= -1 ; 
      result.z *= -1 ; 
      return result ;
    }
    vector3d<real_t> operator*(const vector3d<real_t> &v) const {
      const Quaternion Q = *this ;
      const Quaternion Qinv = Inverse() ;
      const Quaternion vQ(v.x,v.y,v.z,0) ;
      Quaternion result = vQ*Qinv ;
      result = Q*result ;
      return vector3d<real_t>(result.x,result.y,result.z) ;
    }
  } ;

  struct componentXform {
    vector3d<real_t> cg, new_cg ;
    Quaternion q ;
    componentXform() {
      cg = vector3d<real_t>(0,0,0) ;
      new_cg = vector3d<real_t>(0,0,0) ;
      q.x = 0 ;
      q.y = 0 ;
      q.z = 0 ;
      q.w = 1 ;
    }
    vector3d<real_t> applyXform(vector3d<real_t> pos) const {
      return (q*(pos-cg)) + new_cg ;
    }
    vector3d<real_t> applyRotation(vector3d<real_t> vec) const {
      return q*vec ;
    }
  } ;

  inline std::ostream &operator<<(std::ostream &s, const componentXform x) {
    Loci::Abort() ;
    return s ;
  }
  inline std::istream &operator>>(std::istream &s, const componentXform x) {
    Loci::Abort() ;
    return s ;
  }

  class geometry_type : public Loci::CPTR_type {
  public:
    virtual geometry_type *applyXform(componentXform xform) const = 0 ;
    virtual bool inGeometry(vector3d<real_t> pt) const = 0 ;
    virtual real_t distToSurface(vector3d<real_t> pt) const = 0 ;
  } ;

  class motionType : public Loci::CPTR_type {
  public:
    // Compute initial state
    virtual void initialState(std::vector<real_t> &state) const = 0 ;
    // Advance state in time (usually just timestep), sometimes angle
    virtual std::vector<real_t>
      advanceStateInTime(const std::vector<real_t> &state,
			 const std::map<std::string,real_t> &vals) const = 0 ;
    // Get motion transform from current state
    virtual componentXform getXform(const std::vector<real_t> &state) const = 0 ;
  } ;
  
  // Find the time interval that we are splining
  int findt(const std::vector<real_t> &t, real_t tval) ;
  // Compute spline derivatives for hermite spline
  void splineD(std::vector<real_t> &xp, const std::vector<real_t> &x,
               const std::vector<real_t> &t) ;
  // cubic spline
  real_t spline(int ind,real_t tval,const std::vector<real_t> &t,
                const std::vector<real_t> &x,
                const std::vector<real_t> &xp) ;
  
  struct motionSplines {
    std::vector<real_t> t,x,y,z,q0,q1,q2,q3 ;
    std::vector<real_t> xp,yp,zp,q0p,q1p,q2p,q3p ;
    void initialize(std::vector<real_t> ti,std::vector<real_t> xi,
                    std::vector<real_t> yi,std::vector<real_t> zi,
                    std::vector<real_t> q0i,std::vector<real_t> q1i,
                    std::vector<real_t> q2i,std::vector<real_t> q3i) {
      t = ti ;
      x = xi ;
      splineD(xp,x,t) ;
      y = yi ;
      splineD(yp,y,t) ;
      z = zi ;
      splineD(zp,z,t) ;
      q0 = q0i ;
      splineD(q0p,q0,t) ;
      q1 = q1i ;
      splineD(q1p,q1,t) ;
      q2 = q2i ;
      splineD(q2p,q2,t) ;
      q3 = q3i ;
      splineD(q3p,q3,t) ;
    }
    void getMotion(vector3d<real_t> &cg,Quaternion &q,real_t tv) const {
      int ind = findt(t,tv) ;
      cg.x = spline(ind,tv,t,x,xp) ;
      cg.y = spline(ind,tv,t,y,yp) ;
      cg.z = spline(ind,tv,t,z,zp) ;
      q.x = spline(ind,tv,t,q0,q0p) ;
      q.y = spline(ind,tv,t,q1,q1p) ;
      q.z = spline(ind,tv,t,q2,q2p) ;
      q.w = spline(ind,tv,t,q3,q3p) ;
      q.Normalize() ;
    }
  } ;

  template<> struct data_schema_traits<Loci::Quaternion> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::Quaternion()) ;
      LOCI_INSERT_TYPE(ct,Loci::Quaternion,x) ;
      LOCI_INSERT_TYPE(ct,Loci::Quaternion,y) ;
      LOCI_INSERT_TYPE(ct,Loci::Quaternion,z) ;
      LOCI_INSERT_TYPE(ct,Loci::Quaternion,w) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<Loci::componentXform> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::componentXform()) ;
      LOCI_INSERT_TYPE(ct,Loci::componentXform,cg) ;
      LOCI_INSERT_TYPE(ct,Loci::componentXform,new_cg) ;
      LOCI_INSERT_TYPE(ct,Loci::componentXform,q) ;

      return DatatypeP(ct) ;
    }
  } ;
  class interpolate_points {
  public:
    Loci::kdTree::KDTree<float> *kd ;
    store<vector3d<real_t> > pos ;
    store<int> posid ;
    std::vector<int> distribution ;
    interpolate_points &operator=(const interpolate_points &in) {
      std::cerr << "interpolate_data shouldn't be copied!" << std::endl ;
      kd = 0 ;
      pos.setRep(in.pos.Rep()) ;
      return *this ;
    }
    ~interpolate_points() {
      if(kd !=0)
        delete kd ;
      kd = 0 ;
      pos.allocate(EMPTY) ;
      posid.allocate(EMPTY) ;
    }
    interpolate_points() {
      kd = 0 ;
    }
  } ;

  typedef vector3d<float> coord3df ;
  class StencilSizer {
  public:
    Loci::kdTree::KDTree<float> *kd ;
    float ref_size ;
    float getSpacing(coord3df pt) const {
      int id = kd->find_closest(pt) ;
      if(id <= 0)
        return ref_size ;
      return ref_size/float(id) ;
    }


    StencilSizer &operator=(const StencilSizer &in) {
      std::cerr << "interpolate_data shouldn't be copied!" << std::endl ;
      kd = 0 ;
      return *this ;
    }
    ~StencilSizer() {
      if(kd !=0)
        delete kd ;
      kd = 0 ;
    }
    StencilSizer() {
      kd = 0 ;
    }
  } ;

  struct stencil_info {
    std::vector<Loci::Array<real_t,4> > weights ;
    std::vector<Loci::Array<int ,4> > stencils ;
    std::vector<int> send_info, req_sizes, snd_sizes ;
    Loci::storeRepP slookup ;
  } ;
}
#endif
