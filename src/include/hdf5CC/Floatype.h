#ifndef _H5FloatType_H
#define _H5FloatType_H

#include "Tools/stream.h"
extern "C"{
#include <hdf5.h>
}
#include "Predtype.h"
#include "Dataset.h"
#include "Atomtype.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#endif

class FloatType : public AtomType {
   public:
	// default constructor
	FloatType();

	// Creates a floating-point datatype using an existing id
	FloatType( hid_t existing_id );

        // Creates a floating-point type using a predefined type
        FloatType( const PredType& pred_type );

	// Copy constructor: makes a copy of the original FloatType object.
	FloatType( const FloatType& original );

	// Gets the floating-point datatype of the specified dataset        
	FloatType( const DataSet& dataset );

	// Retrieves floating point datatype bit field information. 
	void getFields( size_t* spos, size_t * epos, size_t * esize, size_t * mpos, size_t * msize ) const;

	// Sets locations and sizes of floating point bit fields. 
	void setFields( size_t spos, size_t epos, size_t esize, size_t mpos, size_t msize ) const;

	// Retrieves the exponent bias of a floating-point type. 
	size_t getEbias() const;

	// Sets the exponent bias of a floating-point type. 
	void setEbias( size_t ebias ) const;

	// Retrieves mantissa normalization of a floating-point datatype. 
	H5T_norm_t getNorm( string& norm_string ) const;

	// Sets the mantissa normalization of a floating-point datatype. 
	void setNorm( H5T_norm_t norm ) const;

	// Retrieves the internal padding type for unused bits in floating-point datatypes. 
	H5T_pad_t getInpad( string& pad_string ) const;
	
	// Fills unused internal floating point bits. 
	void setInpad( H5T_pad_t inpad ) const;

	virtual ~FloatType();
};
#ifndef H5_NO_NAMESPACE
}
#endif
#endif
