#include <string>

#ifndef H5_NO_NAMESPACE
using namespace std;
#endif

#include "Exception.h"
#include "Strtype.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

// Default constructor
StrType::StrType() : AtomType() {}

// Creates a string type using a predefined type
StrType::StrType( const PredType& pred_type ) : AtomType()
{
   // use DataType::copy to make a copy of this predefined type
   copy( pred_type );
}

// Creates a string datatype using an existing id
StrType::StrType( const hid_t existing_id ) : AtomType( existing_id ) {}

// Copy constructor: makes copy of the original StrType object
StrType::StrType( const StrType& original ) : AtomType ( original ) {}

// Gets the string datatype of the specified dataset - will reimplement
StrType::StrType( const DataSet& dataset ) : AtomType ()
{
   // Calls C function H5Dget_type to get the id of the datatype
   id = H5Dget_type( dataset.getId() );

   if( id <= 0 )
   {
      throw HDF5DatasetInterfaceException();
   }
}

// Retrieves the character set type of a string datatype. 
H5T_cset_t StrType::getCset() const
{
   H5T_cset_t cset = H5Tget_cset( id );

   // Returns a valid character set type if successful
   if( cset == H5T_CSET_ERROR )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( cset );
}

// Sets character set to be used. 
void StrType::setCset( H5T_cset_t cset ) const
{
   herr_t ret_value = H5Tset_cset( id, cset );

   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// Retrieves the string padding method for a string datatype. 
H5T_str_t StrType::getStrpad() const
{
   H5T_str_t strpad = H5Tget_strpad( id );

   // Returns a valid string padding type if successful
   if( strpad == H5T_STR_ERROR )
   {
      throw HDF5DatatypeInterfaceException();
   }
   return( strpad );
}

// Defines the storage mechanism for character strings. 
void StrType::setStrpad( H5T_str_t strpad ) const
{
   herr_t ret_value = H5Tset_strpad( id, strpad );

   if( ret_value < 0 )
   {
      throw HDF5DatatypeInterfaceException();
   }
}

// This destructor terminates access to the datatype
StrType::~StrType() {}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif
