#ifndef LOCI_H
#define LOCI_H

#include <Loci>

using Loci::entitySet ;
using Loci::create_entitySet ;
using Loci::EMPTY ;
using Loci::interval ;

// Here we use a #define to define Map because bastard MPI implementations fail
// to keep Map in the MPI namespace !@#$^!@^^$%
//using Loci::Map ;
#ifdef Map
#undef Map
#endif
#define Map Loci::Map

using Loci::const_Map ;
using Loci::MapVec ;
using Loci::const_MapVec ;
using Loci::multiMap ;
using Loci::const_multiMap ;
using Loci::inverseMap ;
using Loci::distributed_inverseMap ;

using Loci::store ;
using Loci::const_store ;
using Loci::multiStore ;
using Loci::const_multiStore ;
using Loci::storeVec ;
using Loci::const_storeVec ;
using Loci::Vect ;
using Loci::Mat ;
using Loci::const_Vect ;
using Loci::const_Mat ;
using Loci::Scalar ;
using Loci::mk_Scalar ;
using Loci::pivot_type ;

using Loci::param ;
using Loci::const_param ;

using Loci::constraint ;
using Loci::const_constraint ;
using Loci::Constraint ;

using Loci::unit_rule ;
using Loci::apply_rule ;
using Loci::pointwise_rule ;
using Loci::singleton_rule ;
using Loci::default_rule ;
using Loci::constraint_rule ;
using Loci::map_rule ;
using Loci::optional_rule ;
using Loci::register_rule ;


using Loci::storeMat ;
using Loci::const_storeMat ;

using Loci::dstore ;
using Loci::const_dstore ;

using Loci::dstoreVec;
using Loci::const_dstoreVec;

//using Loci::dstoreMat;
//using Loci::const_dstoreMat;

using Loci::dmultiStore;
using Loci::const_dmultiStore;

using Loci:: dmultiMap;
using Loci:: const_dmultiMap;

using Loci:: dMap;
using Loci:: const_dMap;

using Loci:: dMapVec;
using Loci:: const_dMapVec;

using Loci::sequence ;
using Loci::create_sequence ;
using Loci::do_loop ;
using Loci::Entity ;

using Loci::global_rule_list ;
using Loci::register_rule_list ;
using Loci::rule_db ;

using Loci::fact_db ;
using Loci::executeP ;
using Loci::create_execution_schedule ; 

using Loci::set_fpe_abort ;

using Loci::options_list ;

using Loci::data_schema_traits ;
using Loci::data_schema_converter_traits ;

using Loci::vector3d ;
using Loci::vector2d ;
using Loci::norm ;
using Loci::dot ;
using Loci::cross ;

using Loci::Array ;
using Loci::rule_impl_list ;

using Loci::blackbox;
using Loci::const_blackbox;

#endif





