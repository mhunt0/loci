#include <string>

#ifndef H5_NO_NAMESPACE
using namespace std;
#endif

#include "Exception.h"
#include "Atomtype.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

// Default constructor
AtomType::AtomType() : DataType() {}

// Constructor that takes an existing id 
AtomType::AtomType( const hid_t existing_id ) : DataType( existing_id ) {}

// Copy constructor: makes a copy of the original AtomType object.
AtomType::AtomType( const AtomType& original ) : DataType( original ) {}

// Sets the total size for an atomic datatype. 
void AtomType::setSize( size_t size ) const
{
   // Call C routine H5Tset_size to set the total size
   herr_t ret_value = H5Tset_size( id, size );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Returns the byte order of an atomic datatype.  Inheritance class???
H5T_order_t AtomType::getOrder( string& order_string ) const
{
   // Call C routine to get the byte ordering
   H5T_order_t type_order = H5Tget_order( id );

   // return a byte order constant if successful
   if( type_order == H5T_ORDER_ERROR )
   {
      throw HDF5DatatypeInterfaceException();
   }
   if( type_order == H5T_ORDER_LE )
      order_string = "Little endian byte ordering (0)";
   else if( type_order == H5T_ORDER_BE )
      order_string = "Big endian byte ordering (1)";
   else if( type_order == H5T_ORDER_VAX )
      order_string = "VAX mixed byte ordering (2)";
   return( type_order );
}

// Sets the byte ordering of an atomic datatype.  Inheritance class???
void AtomType::setOrder( H5T_order_t order ) const
{
   // Call C routine to set the byte ordering
   herr_t ret_value = H5Tset_order( id, order );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Returns the precision of an atomic datatype. 
size_t AtomType::getPrecision() const
{
   size_t num_signi_bits = H5Tget_precision( id );  // C routine

   // returns number of significant bits if successful
   if( num_signi_bits == 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( num_signi_bits );
}

// Sets the precision of an atomic datatype. 
void AtomType::setPrecision( size_t precision ) const
{
   // Call C routine to set the datatype precision
   herr_t ret_value = H5Tset_precision( id, precision );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Retrieves the bit offset of the first significant bit. 
size_t AtomType::getOffset() const
{
   size_t offset = H5Tget_offset( id );  // C routine

   // returns a positive offset value if successful
   if( offset == 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( offset );
}

// Sets the bit offset of the first significant bit. 
void AtomType::setOffset( size_t offset ) const
{
   // Call C routine to set the bit offset
   herr_t ret_value = H5Tset_offset( id, offset );
   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Retrieves the padding type of the least and most-significant bit padding. 
// these two are for Opaque type
//void AtomType::getPad( H5T_pad_t * lsb, H5T_pad_t * msb ) const
//{
   // Call C routine to get the padding type
   //herr_t ret_value = H5Tget_pad( id, lsb, msb );
   //if( ret_value < 0 )
   //{
      //throw HDF5DatatypeInterfaceException();
   //}
//}

// Sets the least and most-significant bits padding types
//void AtomType::setPad( H5T_pad_t lsb, H5T_pad_t msb ) const
//{
   // Call C routine to set the padding type
   //herr_t ret_value = H5Tset_pad( id, lsb, msb );
   //if( ret_value < 0 )
   //{
      //throw HDF5DatatypeInterfaceException();
   //}
//}

// This destructor terminates access to the datatype; it calls ~DataType
AtomType::~AtomType() {}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif
