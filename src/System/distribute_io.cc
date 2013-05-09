//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;

#include <algorithm>
using std::sort;
using std::pair ;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>

namespace Loci {

  entitySet BcastEntitySet(entitySet set, int root, MPI_Comm comm) {

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    // Now lets share the domain with all other processors ;
    int sz = set.num_intervals() ;
    MPI_Bcast(&sz,1,MPI_INT,root,comm) ;
    vector<interval> vlist(sz) ;
    if(rank == 0) {
      for(int i=0;i<sz;++i)
        vlist[i] = set[i] ;
    }
    MPI_Bcast(&vlist[0],sz*2,MPI_INT,root,comm) ;
    set = EMPTY ;
    for(int i = 0;i<sz;++i)
      set += vlist[i] ;
    return set ;
  }

  vector<int> simplePartitionVec(int mn, int mx, int p) {
    vector<int> nums(p+1) ;
    int n = mx-mn+1 ;
    int dn = n/p ; // divisor
    int rn = n%p ; // remainder
    int start = mn ;
    nums[0] = start ;
    for(int i=0;i<p;++i) {
      start += dn+((i<rn)?1:0) ;
      nums[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return nums ;
  }

  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    vector<int> pl = simplePartitionVec(mn,mx,p) ;
    vector<entitySet> ptn(p) ;
    for(int i=0;i<p;++i)
      ptn[i] = interval(pl[i],pl[i+1]-1) ;
    return ptn ;
  }

  void write_frame_info_param(hid_t group_id, frame_info &fi, MPI_Comm comm) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Write out is_stat and vector size
    if(prank==0) {
      hsize_t dimension = 1 ;
      int rank = 1 ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      if(fi.is_stat != 0) {
        rank = 1 ;
        dimension = 1 ;
        dataspace = H5Screate_simple(rank,&dimension,NULL) ;
        hid_t dataset = H5Dcreate(group_id,"second_level",H5T_NATIVE_INT,
                                  dataspace,H5P_DEFAULT) ;
        hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
        H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &fi.second_level[0]) ;
        H5Sclose(memspace) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
    }
  }

  void write_frame_info(hid_t group_id, frame_info &fi, MPI_Comm comm) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Write out is_stat and vector size
    if(prank==0) {
      hsize_t dimension = 1 ;
      int rank = 1 ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,H5P_DEFAULT) ;
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
    int is_stat = fi.is_stat ;
    int size = fi.size ;
    MPI_Bcast(&is_stat,1,MPI_INT,0,comm) ;
    MPI_Bcast(&size,1,MPI_INT,0,comm) ;
    if(size == 0) { // Two level framing
      write_vector_int(group_id, "first_level", fi.first_level,comm) ;
      if(is_stat != 0) {
        write_vector_int(group_id, "second_level", fi.second_level,comm) ;
      }
    } else {
      if(is_stat != 0) {
        write_vector_int(group_id,"second_level",fi.second_level,comm) ;
      }
    }
  }

  frame_info read_frame_info_param(hid_t group_id,  int dom_size, MPI_Comm comm) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    int is_stat = 0 ;
    int sz = 0 ;
    // Write out is_stat and vector size
    frame_info fi ;
    if(prank == 0) {
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dopen(group_id, "is_stat") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dopen(group_id, "vec_size") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
      H5Dclose(dataset) ;
      fi.is_stat = is_stat ;
      fi.size = sz ;
      if(is_stat != 0) {
        hid_t dataset = H5Dopen(group_id, "second_level") ;
        hid_t dataspace = H5Dget_space(dataset) ;
        hsize_t dimension = 1 ;
        H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
        std::vector<int> vint ;
        int tmp ;
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT, &tmp) ;
        vint.push_back(tmp) ;
        fi.second_level = vint ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
    }
    return fi ;
  }

  frame_info read_frame_info(hid_t group_id,  int dom_size, MPI_Comm comm) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    int is_stat = 0 ;
    int sz = 0 ;
    // Write out is_stat and vector size
    if(prank == 0) {
      hid_t datatype = H5T_NATIVE_INT ;
      hid_t dataset = H5Dopen(group_id, "is_stat") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
      H5Dclose(dataset) ;
      dataset = H5Dopen(group_id, "vec_size") ;
      H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
      H5Dclose(dataset) ;
    }
    int dim[2] ;
    dim[0] = is_stat ;
    dim[1] = sz ;
    MPI_Bcast(&dim, 2, MPI_INT, 0, comm) ;
    frame_info fi = frame_info(dim[0], dim[1]) ;

    if(dim[1] == 0) { // level 1 framing
      read_vector_int(group_id, "first_level", fi.first_level, dom_size,comm) ;
      int total_size = fi.first_level.size() ;
      if(dim[0] != 0) {
        int dims = 0 ;
        for(int i = 0; i < total_size; ++i)
          dims += (fi.first_level)[i] ;
        read_multi_vector_int(group_id, "second_level", dims, fi.second_level,comm) ;
      }
    } else { // level 2 framing only
      if(dim[0] !=0) {
        std::vector<int> vint ;
        read_vector_int(group_id, "second_level", vint, dom_size,comm) ;
        fi.second_level = vint ;
      }
    }
    return fi ;
  }

  void write_parameter(hid_t group_id, storeRepP qrep, MPI_Comm comm) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    frame_info fi = qrep->get_frame_info() ;
    write_frame_info_param(group_id,fi,comm) ;

    if(prank==0) {
      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          array_size = fi.size  ;
      else
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;

      if(array_size == 0)
        array_size = 1 ;
      hsize_t dimension = array_size ;
      int rank = 1 ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = array_size ;

      hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      dimension = count ;
      start += dimension ;
      hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace, H5P_DEFAULT) ;
      entitySet dom = ~EMPTY ;
      qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
  }

  void write_store(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm) {
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    // Shift domain by offset
    entitySet dom_file = dom >> offset ;

    // Compute overall domain across processors
    std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
    entitySet q_dom;
    for(int i = 0; i < np; i++)
      q_dom += dom_vector[i];

    if(prank == 0)
      HDF5_WriteDomain(group_id, q_dom);

    // If nothing to write, don't proceed
    if(q_dom == EMPTY)
      return ;

    // Allocate buffer for largest processor buffer size
    std::vector<int> sort_max ;
    int local_size = qrep->pack_size(dom) ;
    sort_max = all_collect_sizes(local_size,comm) ;
    int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
    vector<unsigned char> tmp_send_buf(total_size) ;


    frame_info fi = qrep->get_frame_info() ;
    write_frame_info(group_id,fi,comm) ;


    int array_size = 0 ;
    if(fi.size)
      if(fi.is_stat)
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	array_size = fi.size * dom.size() ;
    else
      if(fi.is_stat)
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
      else
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	  array_size += *fvi ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
    size_t tot_arr_size = 0 ;
    for(int i = 0; i < np; ++i)
      tot_arr_size += size_t(max(0,arr_sizes[i])) ;


    if(prank != 0) {
      // Client processor code, pack data into send buffer
      MPI_Status status ;
      int send_size_buf ;
      send_size_buf = qrep->pack_size(dom) ;
      int tot_size = send_size_buf ;
      int loc_pack = 0 ;
      qrep->pack(&tmp_send_buf[0], loc_pack, total_size, dom) ;
      // Wait for signal to send message to root processor
      int flag = 0 ;
      MPI_Recv(&flag,1, MPI_INT, 0, 10, comm, &status) ;
      // received token to send, so send message
      if(flag) {
	MPI_Send(&tot_size, 1, MPI_INT, 0, 11, comm) ;
	MPI_Send(&tmp_send_buf[0], tot_size, MPI_PACKED, 0, 12, comm) ;
      }
    } else {
      // Begin writing array
      int rank = 1 ;
      hsize_t dimension = 1 ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = arr_sizes[0] ;
      dimension =  tot_arr_size ;
      if(dimension != 0) {
        // First write local data
        hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
	DatatypeP dp = qrep->getType() ;
        hid_t datatype = dp->get_hdf5_type() ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        dimension = count ;
        start += dimension ;
        hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace, H5P_DEFAULT) ;
        qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
        H5Dclose(dataset) ;

        // Now write remaining vectors

	for(int i = 1; i < np; ++i) {
	  MPI_Status status ;
	  int recv_total_size ;
          // Allocate over 0-size-1, this allows for greater scalability when
          // sets data exceeds 2gig
	  entitySet tmpset = interval(0,dom_vector[i].size()-1);

          storeRepP t_qrep = qrep->new_store(tmpset) ;

	  int loc_unpack = 0 ;
	  int flag = 1 ;
	  MPI_Send(&flag, 1, MPI_INT, i, 10, comm) ;
	  MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, comm, &status) ;
	  MPI_Recv(&tmp_send_buf[0], recv_total_size, MPI_PACKED, i, 12, comm, &status) ;

	  sequence tmp_seq = sequence(tmpset) ;
          t_qrep->unpack(&tmp_send_buf[0], loc_unpack, total_size, tmp_seq) ;
	  dimension = arr_sizes[i] ;
	  count = dimension ;

          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
	  start += count ;

          dataset = H5Dopen(group_id, "data") ;
          t_qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", tmpset) ;
          t_qrep->allocate(EMPTY) ;

          H5Dclose(dataset) ;
	}
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
      }
      //add else part by Qiuhan to avoid MPI communication get stuck
      else{
	for(int i = 1; i < np; ++i) {
	  MPI_Status status ;
	  int recv_total_size ;
          int flag = 1 ;
          MPI_Send(&flag, 1, MPI_INT, i, 10, comm) ;
	  MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, comm, &status) ;
	  MPI_Recv(&tmp_send_buf[0], recv_total_size, MPI_PACKED, i, 12, comm, &status) ;
        }
      }
    }
  }

  void write_container(hid_t group_id, storeRepP qrep) {
    int offset = 0 ;
    if(qrep->RepType() == PARAMETER)
      write_parameter(group_id,qrep,MPI_COMM_WORLD) ;
    else
      write_store(group_id,qrep,qrep->domain(),offset,MPI_COMM_WORLD) ;
  }

  void read_parameter(hid_t group_id, storeRepP qrep, MPI_Comm comm) {
    int pack_size = 0 ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    frame_info fi = read_frame_info_param(group_id,1,comm) ;
    if(prank==0) {
      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
          array_size = fi.size ;
        } else {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
        }
      }

      hid_t dimension = array_size ;
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
      entitySet dom = ~EMPTY ;
      qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;

      pack_size = qrep->pack_size(dom) ;

    }

    // Now broadcast the result to other processors
    if(np > 1) {
      MPI_Bcast(&pack_size,1,MPI_INT,0,comm) ;
      unsigned char *pack_buf = new unsigned char[pack_size] ;
      if(prank == 0) {
        int loc_pack = 0 ;
        int sz = pack_size ;
        entitySet dom = ~EMPTY ;
        qrep->pack(pack_buf,loc_pack,sz,dom) ;
      }
      MPI_Bcast(pack_buf,pack_size,MPI_PACKED,0,comm) ;
      if(prank != 0) {
        int loc_pack = 0 ;
        int sz = pack_size ;
        entitySet dom = ~EMPTY ;
        qrep->unpack(pack_buf,loc_pack,sz,dom) ;
      }
      delete[] pack_buf ;
    }
  }

  void read_store(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm) {
    offset = 0 ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    // Here we read in a store container.  First lets read in the domain
    entitySet q_dom ;
    if(prank == 0)
      HDF5_ReadDomain(group_id, q_dom) ;

    // Now lets share the domain with all other processors ;
    q_dom = BcastEntitySet(q_dom,0,comm) ;

    unsigned char* tmp_buf = 0;
    std::vector<int> interval_sizes ;
    entitySet dom ;
    if(q_dom != EMPTY) {
      vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
      for(int i=0;i<np;++i) {
        entitySet qset = ptn[i] &q_dom ;
        interval_sizes.push_back(qset.size()) ;
      }
      dom = ptn[prank] &q_dom ;
    } else
      for(int i=0;i<np;++i)
        interval_sizes.push_back(0) ;

    offset = dom.Min() ;
    dom <<= offset ;
    qrep->allocate(dom) ;

    frame_info fi = read_frame_info(group_id,dom.size(),comm) ;
    int array_size = 0 ;
    int vec_size = 0 ;

    if(fi.size) {
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
	vec_size = fi.second_level.size() ;
      } else {
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	array_size = fi.size * dom.size() ;
      }
    } else {
      if(fi.is_stat) {
	for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	  array_size += *vi ;
	vec_size = fi.second_level.size() + dom.size() ;
      } else {
	for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	  array_size += *fvi ;
	vec_size = dom.size() ;
      }
    }

    std::vector<int> tmp_sizes = all_collect_sizes(vec_size,comm) ;
    int max_tmp_size = *std::max_element(tmp_sizes.begin(), tmp_sizes.end()) ;
    int max_eset_size = *std::max_element(interval_sizes.begin(), interval_sizes.end()) ;
    int* tmp_int  ;
    tmp_int = new int[max_tmp_size] ;
    std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
    size_t tot_arr_size = 0 ;
    for(int i = 0; i < np; ++i)
      tot_arr_size += size_t(max(0,arr_sizes[i])) ;
    MPI_Status status ;
    if(prank != 0) {
      int t = 0 ;
      if(fi.size) {
	if(fi.size > 1)
	  qrep->set_elem_size(fi.size) ;
	if(fi.is_stat)
	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	    tmp_int[t++] = *vi ;
      } else {
	if(fi.is_stat) {
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;

	  for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
	    tmp_int[t++] = *vi ;
	} else
	  for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
	    tmp_int[t++] = *fvi ;
      }
      if(tmp_sizes[prank])
	MPI_Send(tmp_int, tmp_sizes[prank], MPI_INT, 0, 10, comm) ;
      int total_size = 0 ;
      MPI_Recv(&total_size, 1, MPI_INT, 0, 11,comm, &status) ;
      tmp_buf = new unsigned char[total_size] ;
      MPI_Recv(tmp_buf, total_size, MPI_PACKED, 0, 12, comm, &status) ;
      sequence tmp_seq = sequence(dom) ;
      int loc_unpack = 0 ;
      qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
    } else {
      // processor zero
      hid_t dataset =  H5Dopen(group_id, "data") ;
      hid_t dataspace = H5Dget_space(dataset) ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
#else
      hssize_t start = 0 ;
#endif
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      int curr_indx = 0 ;
      int total_size = 0 ;
      int tmp_total_size = 0 ;
      entitySet max_set;
      if(max_eset_size > 0)
	max_set = interval(0, max_eset_size-1) ;
      if(np == 1) {
        qrep->allocate(dom) ;
        if(fi.size)
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
        hsize_t dimension = arr_sizes[0] ;
        count = dimension ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
      } else {
        storeRepP tmp_sp ;
        if(fi.size)
          tmp_sp = qrep->new_store(max_set) ;
        for(int p = 0; p < np; ++p) {
          entitySet local_set;
          if(interval_sizes[p] > 0)
            local_set = entitySet(interval(curr_indx, interval_sizes[p]+curr_indx-1)) ;
          curr_indx += interval_sizes[p] ;
          hsize_t dimension = arr_sizes[p] ;
          count = dimension ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          entitySet tmp_set;
          if(local_set.size())
            tmp_set = interval(0, local_set.size()-1) ;
          if(p && tmp_sizes[p]) {
            MPI_Recv(tmp_int, tmp_sizes[p], MPI_INT, p, 10, comm, &status) ;
            std::vector<int> vint, fvint ;
            int t = 0 ;
            if(fi.size) {
              if(fi.is_stat) {
                for(int i = 0; i < tmp_sizes[p]; ++i)
                  vint.push_back(tmp_int[t++]) ;
                fi.second_level = vint ;
              }
            } else {
              for(int i = 0; i < local_set.size(); ++i)
                fvint.push_back(tmp_int[t++]) ;
              for(int i = 0; i < tmp_sizes[p]-local_set.size(); ++i)
                vint.push_back(tmp_int[t++]) ;
              fi.first_level = fvint ;
              fi.second_level = vint ;
            }
          }
          storeRepP t_sp ;
          int t = 0 ;
          if(p == 0)
            if(!fi.size)
              for(std::vector<int>::const_iterator vi = fi.first_level.begin(); vi != fi.first_level.end(); ++vi)
                tmp_int[t++] = *vi ;
          if(fi.size) {
            tmp_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ;
            tmp_total_size = tmp_sp->pack_size(tmp_set) ;
          } else {
            t_sp = qrep->new_store(tmp_set, tmp_int) ;
            t_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ;
            tmp_total_size = t_sp->pack_size(tmp_set) ;
          }

          if(tmp_total_size > total_size) {
            total_size = tmp_total_size ;
            if(p)
              delete [] tmp_buf ;
            tmp_buf = new unsigned char[total_size] ;
          }
          start += count ;
          int loc = 0 ;
          if(fi.size)
            tmp_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
          else
            t_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
          if(p == 0) {
            int loc_unpack = 0 ;

            if(fi.size)
              if(fi.size > 1)
                qrep->set_elem_size(fi.size) ;
            
            sequence tmp_seq = sequence(dom) ;
            qrep->allocate(dom) ;
            qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
          } else {
            MPI_Send(&total_size, 1, MPI_INT, p, 11, comm) ;
            MPI_Send(tmp_buf, total_size, MPI_PACKED, p, 12, comm) ;
          }
        }
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    }
    delete [] tmp_buf ;
    delete [] tmp_int ;
  }

  entitySet findBoundingSet(entitySet in) {
    Entity max_val = in.Max() ;
    Entity min_val = in.Min() ;
    Entity gmin_val = min_val;
    Entity gmax_val = max_val ;
    MPI_Allreduce(&min_val,&gmin_val,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
    MPI_Allreduce(&max_val,&gmax_val,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

    return entitySet(interval(gmin_val,gmax_val)) ;
  }

  vector<sequence> transposeSeq(const vector<sequence> sv) {
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = sv[i].num_intervals()*2 ;
    vector<int> recv_sz(MPI_processes) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }
    //    outRep->allocate(new_alloc) ;
    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<sv[i].num_intervals();++j) {
        send_store[send_displacement[i]+j*2] = sv[i][j].first ;
        send_store[send_displacement[i]+j*2+1] = sv[i][j].second ;
      }


    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;

    vector<sequence> sv_t(MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        sv_t[i] += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return sv_t ;
  }


  // Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm) {


    // Get local numbering of entities owned by this processor, only write
    // out these entities.
    constraint my_entities ;
    my_entities = dist->my_entities ;
    dom = *my_entities & dom ;

    // Get mapping from local to global numbering
    Map l2g ;
    l2g = dist->l2g.Rep() ;
    // Compute domain in global numbering
    entitySet dom_global = l2g.image(dom) ;

    // This shouldn't happen
    FATAL(dom.size() != dom_global.size()) ;

    // Now get global to file numbering
    dMap g2f ;
    g2f = dist->g2f.Rep() ;

    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[l2g[i]] ;
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    // Find bounds in file numbering from this processor
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;

    // Find overall bounds
    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Get partitioning of file numbers across processors
    vector<entitySet> out_ptn = simplePartition(imn,imx,comm) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;

    // Loop over processors and compute sets of entities to send
    // To efficiently compute this mapping, first sort the transpose
    // of the newnum map to quickly find the set of entities to send
    // without searching entire newnum map for each processor
    vector<pair<int,int> > file2num(dom.size()) ;
    size_t cnt = 0 ;
    FORALL(dom,ii) {
      file2num[cnt].first = newnum[ii] ;
      file2num[cnt].second = ii ;
      cnt++ ;
    } ENDFORALL ;
    sort(file2num.begin(),file2num.end()) ;

    // Check each processor, find out which sets to send
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = out_ptn[i].Max() ;
      while(cnt < file2num.size() && file2num[cnt].first <= mxi) {
        send_sets[i] += file2num[cnt].second ;
        cnt++ ;
      }
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;


    // shift by the offset
    offset = out_ptn[prank].Min() ;
    for(int i=0;i<p;++i)
      recv_seqs[i] <<= offset ;

    // Compute allocation domain
    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    // allocate store over shifted domain
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;

    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
		  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    return qcol_rep ;
  }

  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm) {
    Map newnum ;
    newnum.allocate(resultSet) ;

    if(dist !=0 ) {
      dMap g2f ;
      g2f = dist->g2f.Rep() ;
      Map l2g ;
      l2g = dist->l2g.Rep() ;
      FORALL(resultSet,i) {
        newnum[i] = g2f[l2g[i]] ;
      } ENDFORALL ;
    } else {
      result->copy(input,resultSet) ;
      return ;
    }

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int mn = input->domain().Min() ;
    int mx = input->domain().Max() ;
    if(input->domain() != EMPTY) {
      mn += offset ;
      mx += offset ;
    }
    vector<int> allmx(p) ;
    vector<int> allmn(p) ;
    MPI_Allgather(&mx,1,MPI_INT,&allmx[0],1,MPI_INT,comm) ;
    MPI_Allgather(&mn,1,MPI_INT,&allmn[0],1,MPI_INT,comm) ;

    vector<pair<int,int> > file_requests ;
    FORALL(resultSet,i) {
      file_requests.push_back(pair<int,int>(newnum[i],i)) ;
    } ENDFORALL ;
    sort(file_requests.begin(),file_requests.end()) ;
    // Get distribution plan
    vector<vector<pair<int,int> > > dist_plan(p) ;

    int proc = 0 ;
    for(size_t i=0;i<file_requests.size();++i) {
      int fn = file_requests[i].first ;
      while(proc < p && (fn < allmn[proc] || fn > allmx[proc]))
        proc++ ;
      if(fn < allmn[proc] || fn > allmx[proc]) {
        cerr << "Unable to find processor that contains entity!" << endl ;
        Loci::Abort() ;
      }
      dist_plan[proc].push_back(pair<int,int>(fn,file_requests[i].second)) ;
    }

    // Compute recv requests from distribution plan
    vector<sequence> recv_seq(p),send_req(p) ;
    for(int i=0;i<p;++i) {
      sort(dist_plan[i].begin(),dist_plan[i].end()) ;
      sequence s1,s2 ;
      int psz = dist_plan[i].size() ;
      for(int j=0;j<psz;++j) {
        s1 +=dist_plan[i][j].first ;
        s2 +=dist_plan[i][j].second ;
      }
      send_req[i] = s1 ;
      recv_seq[i] = s2 ;
    }

    // Transpose the send requests to get the sending sequences
    // from this processor
    vector<sequence> send_seq = transposeSeq(send_req) ;
    vector<entitySet> send_sets(p) ;
    for(int i=0;i<p;++i) {
      send_seq[i] <<= offset ;
      send_sets[i] = entitySet(send_seq[i]) ;
    }

    vector<int> send_sizes(p), recv_sizes(p) ;


    for(int i=0;i<p;++i)
      send_sizes[i] = input->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,&recv_sizes[0],1,MPI_INT, comm) ;

    vector<int> send_dspl(p), recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz), recv_store(recv_sz) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      input->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
                  send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
		  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      result->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                     recv_seq[i]) ;
    }
  }


  void writeContainerRAW(hid_t file_id, std::string vname,
                         storeRepP var, MPI_Comm comm) {
    hid_t group_id = 0 ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    if(prank == 0)
      group_id = H5Gcreate(file_id, vname.c_str(), 0) ;

    if(var->RepType() != PARAMETER) {
      int offset = 0 ;
      write_store(group_id,var,var->domain(),offset,comm) ;
    } else {
      write_parameter(group_id, var,comm) ;
    }
    if(prank == 0)
      H5Gclose(group_id) ;
  }

  void readContainerRAW(hid_t file_id, std::string vname,
                        storeRepP var,
                        MPI_Comm comm ) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    hid_t group_id = 0;
    if(prank == 0)
      group_id = H5Gopen(file_id, vname.c_str()) ;
    if(var->RepType() == PARAMETER) {
      read_parameter(group_id, var, comm) ;
      if(prank == 0)
        H5Gclose(group_id) ;
      return ;
    }

    // Read in store in file numbering
    int offset = 0 ;
    read_store( group_id, var,offset,comm) ;
    var->shift(offset) ;

    if(prank == 0)
      H5Gclose(group_id) ;
  }

  void redistribute_write_container(hid_t file_id, std::string vname,
                                    storeRepP var, fact_db &facts) {
    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist == 0) {
      writeContainerRAW(file_id,vname,var,MPI_COMM_WORLD) ;
      return ;
    }

    hid_t group_id = 0 ;
    if(MPI_rank == 0)
      group_id = H5Gcreate(file_id, vname.c_str(), 0) ;

    // Redistribute container to map from local to global numbering
    if(var->RepType() != PARAMETER && MPI_processes != 1) {
      // parallel store write.. reorder to file numbering then write out
      // reorder from local to file numbering
      int offset = 0 ;
      entitySet dom = var->domain() ;
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      storeRepP vardist = Local2FileOrder(var,dom,offset,dist,MPI_COMM_WORLD) ;
      // Write out container that has been distributed in the file numbering
      write_store(group_id,vardist,vardist->domain(),offset,MPI_COMM_WORLD) ;
    } else {
      // No need to reorder container if parameter or only one
      // processor, so just write out the container.
      write_container(group_id, var) ;
    }

    if(MPI_rank == 0)
      H5Gclose(group_id) ;
  }

  void read_container_redistribute(hid_t file_id, std::string vname,
                                   storeRepP var, entitySet read_set,
                                   fact_db &facts) {
    hid_t group_id = 0;
    if(MPI_rank == 0)
      group_id = H5Gopen(file_id, vname.c_str()) ;
    if(var->RepType() == PARAMETER) {
      read_parameter(group_id, var, MPI_COMM_WORLD) ;
      if(MPI_rank == 0)
        H5Gclose(group_id) ;
      return ;
    }

    // Read in store in file numbering
    int offset = 0 ;
    storeRepP new_store = var->new_store(EMPTY) ;
    read_store( group_id, new_store,offset,MPI_COMM_WORLD) ;
      
    
    // map from file number to local numbering
    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist != 0) {
      // Allocate space for reordered container
      storeRepP result = var->new_store(read_set) ;
      File2LocalOrder(result,read_set,new_store,offset,dist,MPI_COMM_WORLD) ;
      // Copy results into container
      if(read_set == EMPTY) {
        read_set = result->domain() ;
        var->allocate(read_set) ;
      }
      var->copy(result,read_set) ;
    } else {
      if(offset != 0) {
        // shift new store by offset to correct alignment
        new_store->shift(offset) ;
      }
      if(read_set == EMPTY) {
        read_set = new_store->domain() ;
        var->allocate(read_set) ;
      }
      var->copy(new_store,read_set) ;
    }

    if(MPI_rank == 0)
      H5Gclose(group_id) ;
  }

  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) {
    vector<int> ids(local_set.size()) ;

    int c = 0 ;
    if(MPI_processes > 1) {
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      dMap g2f ;
      g2f = df->g2f.Rep() ;
      FORALL(local_set,ii) {
        ids[c++] = g2f[l2g[ii]] ;
      } ENDFORALL ;
    } else {
      // Note, this means that a single processor run will not be compatible
      // with a parallel processor run
      FORALL(local_set,ii) {
        ids[c++] = ii ;
      } ENDFORALL ;
    }
    writeUnorderedVector(file_id,"entityIds",ids) ;
  }

  hid_t createUnorderedFile(const char * filename, entitySet set, fact_db &facts) {
    hid_t file_id = 0;
    hid_t group_id = 0 ;
    if(MPI_rank == 0) {
      file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
      group_id = H5Gcreate(file_id,"dataInfo",0) ;
    }
    writeSetIds(group_id,set,facts) ;
    if(MPI_rank == 0)
      H5Gclose(group_id) ;
    return file_id ;
  }

  void closeUnorderedFile(hid_t file_id) {
    if(MPI_rank == 0)
      H5Fclose(file_id) ;
  }

}
