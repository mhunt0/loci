#include "dist_internal.h"
#include <vector>
#include <mpi.h>

namespace Loci {

  void read_vector_int(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size) {
    std::vector<int> vec_size = all_collect_sizes(dom_size) ;
    hsize_t dimension = 0 ;
    hid_t dataset = 0;
    hid_t dataspace = 0;
    if(Loci::MPI_rank == 0) {
      dataset = H5Dopen(group_id, name) ;
      dataspace = H5Dget_space(dataset) ;
      H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    }
    int dim = dimension ;
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
    int rank = 1 ;
    hid_t datatype = H5T_NATIVE_INT ;
    int total_size = dim / MPI_processes ;
    std::vector<int> sizes = all_collect_sizes(dom_size) ;
    total_size = *std::max_element(sizes.begin(), sizes.end() );
    int *tmp_int = new int[total_size] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      MPI_Recv(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      for(int i = 0; i < sizes[MPI_rank]; ++i)
	vint.push_back(tmp_int[i]) ;
    } else { 
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	dimension = sizes[p] ;
	count = dimension ;
	if(dimension != 0) {
	  hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	  hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
			      H5P_DEFAULT, tmp_int) ;
	  if(err < 0) {
	    cerr << "H5Dread() failed" << endl ;
	  }
	  H5Sclose(memspace) ;
	}
	start += count ;
	if(p == 0) {
	  for(int i = 0; i < sizes[p]; ++i) 
	    vint.push_back(tmp_int[i]) ;
	} else 
	  MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, MPI_COMM_WORLD) ;
      }
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
    }
    delete [] tmp_int ; 
  }

  void read_multi_vector_int(hid_t group_id, const char* name, int dim,  std::vector<int>& vint) {
    hsize_t dimension = 0 ;
    hid_t dataset = 0;
    hid_t dataspace = 0;
    if(Loci::MPI_rank == 0) {
      dataset = H5Dopen(group_id, name) ;
      dataspace = H5Dget_space(dataset) ;
      H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
    }
    int rank = 1 ;
    hid_t datatype = H5T_NATIVE_INT ;
    int total_size = 0 ;
    std::vector<int> sizes = all_collect_sizes(dim) ;
    total_size = *std::max_element(sizes.begin(), sizes.end());
    int *tmp_int = new int[total_size] ;
    MPI_Status status ;
    if(Loci::MPI_rank != 0) {
      MPI_Recv(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12,
	       MPI_COMM_WORLD, &status) ;  
      for(int i = 0; i < sizes[MPI_rank]; ++i)
	vint.push_back(tmp_int[i]) ;
    } else { 
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      for(int p = 0; p < Loci::MPI_processes; ++p) {
	dimension = sizes[p] ;
	count = dimension ;
	if(dimension != 0) {
	  hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	  hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
			      H5P_DEFAULT, tmp_int) ;
	  if(err < 0) {
	    cerr << "H5Dread() failed" << endl ;
	  }
	  H5Sclose(memspace) ;
	}
	start += count ;
	if(p == 0) {
	  for(int i = 0; i < sizes[p]; ++i) 
	    vint.push_back(tmp_int[i]) ;
	} else 
	  MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, MPI_COMM_WORLD) ;
      }
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
    }
    delete [] tmp_int ; 
  }
  
  void write_vector_int(hid_t group_id, const char* name, std::vector<int>& vint) {
    std::vector<int> sort_max(MPI_processes) ;
    int tot_entities = vint.size() ;
    sort_max = all_collect_sizes(tot_entities) ;
    std::vector<int> sizes = sort_max ;
    std::sort(sort_max.begin(), sort_max.end()) ;
    tot_entities = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      tot_entities += sort_max[i] ;
    int *tmp_int = new int[sort_max[MPI_processes-1]] ;
    int tmp = 0 ;
    for(std::vector<int>::iterator vi = vint.begin(); vi != vint.end(); ++vi)
      tmp_int[tmp++] = *vi ;
    if(Loci::MPI_rank != 0) {
      MPI_Status status ;
      int flag = 0 ;
      MPI_Recv(&flag,1, MPI_INT, 0, 11, MPI_COMM_WORLD, &status) ;
      if(flag)
	MPI_Send(tmp_int, sizes[MPI_rank], MPI_INT, 0, 12, MPI_COMM_WORLD) ;
    } else {
      hid_t datatype = H5T_NATIVE_INT ;
      int rank = 1 ;
      hssize_t start = 0 ; 
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      hsize_t dimension = tot_entities ;
      if(dimension != 0) {
	hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
	count = sizes[0] ;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
	dimension = sizes[0] ;
	start += dimension ;
	if(dimension != 0) {
	  hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	  hid_t dataset = H5Dcreate(group_id, name , datatype, dataspace,H5P_DEFAULT) ;
	  H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
	  H5Dclose(dataset) ;
	  H5Sclose(memspace) ;
	}
	for(int i = 1; i < Loci::MPI_processes; ++i) {
	  MPI_Status status ;
	  int flag = 1 ;
	  MPI_Send(&flag, 1, MPI_INT, i, 11, MPI_COMM_WORLD) ;
	  MPI_Recv(tmp_int, sizes[i], MPI_INT, i, 12, MPI_COMM_WORLD, &status) ;
	  dimension = sizes[i] ;
	  count = dimension ;
	  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
	  start += count ;
	  if(dimension != 0) {
	    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
	    hid_t dataset = H5Dopen(group_id, name) ;
	    H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
	    H5Dclose(dataset) ;
	    H5Sclose(memspace) ;
	  }
	}
	H5Sclose(dataspace) ;
      }
    }
    delete [] tmp_int ;
  }

  std::vector<int> all_collect_sizes(int size) {
    std::vector<int> vset( MPI_processes) ;
    if(MPI_processes > 1) {
      int *recv_count = new int[ MPI_processes] ;
      int *send_count = new int[ MPI_processes] ;
      for(int i = 0; i <  MPI_processes; ++i) 
	send_count[i] = size ;
      
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ; 
      for(int i = 0; i <  MPI_processes; ++i)
	vset[i] = recv_count[i] ;
      
      delete [] send_count ;
      delete [] recv_count ;
    }
    else
      vset[0] = size ;
    return vset ;
  }

}
