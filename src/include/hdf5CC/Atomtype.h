#ifndef _H5AtomType_H
#define _H5AtomType_H


#include "Tools/stream.h"
extern "C"{
#include <hdf5.h>
}
#include "Datatype.h"


#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif
class AtomType : public DataType {
   protected:
	// Default constructor
	AtomType();

	// Constructor that takes an existing id - for predefined type
	AtomType( const hid_t existing_id );

   public:
	// Copy constructor - makes copy of the original object
	AtomType( const AtomType& original );

	// Sets the total size for an atomic datatype. 
	void setSize( size_t size ) const;

	// Returns the byte order of an atomic datatype. 
	H5T_order_t getOrder( string& order_string ) const;

	// Sets the byte ordering of an atomic datatype. 
	void setOrder( H5T_order_t order ) const;

	// Returns the precision of an atomic datatype. 
	size_t getPrecision() const;

	// Sets the precision of an atomic datatype. 
	void setPrecision( size_t precision ) const;

	// Retrieves the bit offset of the first significant bit. 
	size_t getOffset() const;

	// Sets the bit offset of the first significant bit. 
	void setOffset( size_t offset ) const;

	// The followings will go into Opaque type when completed
	// Retrieves the padding type of the least and most-significant bit padding. 
	// void getPad( H5T_pad_t * lsb, H5T_pad_t * msb ) const;

	// Sets the least and most-significant bits padding types
	// void setPad( H5T_pad_t lsb, H5T_pad_t msb ) const;

	virtual ~AtomType();
};
#ifndef H5_NO_NAMESPACE
}
#endif
#endif
