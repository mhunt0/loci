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
#include <constraint.h>
#include "dist_tools.h"
#include <Tools/debug.h>
#include <entitySet.h>

#include "Tools/debugger.h"
#include "Tools/hash_map.h"

#include <vector>
using std::vector ;

#include <set>
using std::set ;

#include <mpi.h>

#include <iostream>
using std::endl ;

#include <algorithm>
using std::sort ;

namespace Loci {
  
  entitySet collect_entitySet(entitySet e, fact_db &facts) {
    if(!facts.isDistributed())
      return e ;
    entitySet re ;
    if(facts.isDistributed()) {  
      Map l2g ;
      entitySet::const_iterator ti ;
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      d = facts.get_distribute_info() ;
      l2g = d->l2g.Rep() ;
      if(MPI_processes == 1) {
	entitySet temp = e ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  re += l2g[*ti] ;
        return re ;
      }
      if(d->myid == 0) {
	MPI_Status *status, *size_status ;
	MPI_Request *recv_request, *size_request ;
	int **recv_buffer ;
	int *recv_size ;
	int k = 0 ;
	entitySet temp = e & d->my_entities ;
	for(ti = temp.begin(); ti != temp.end(); ++ti)
	  re += l2g[*ti] ;
        
	recv_size = new int[MPI_processes-1] ;
	size_request = new MPI_Request[MPI_processes-1] ;
	size_status = new MPI_Status[MPI_processes-1] ;
	for(k = 0; k < MPI_processes-1; k++) {
	  MPI_Irecv(&recv_size[k],1,MPI_INT, k+1,1, MPI_COMM_WORLD, &size_request[k]);
        }
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(MPI_processes-1, size_request, size_status) ;
	FATAL(err != MPI_SUCCESS) ;
	recv_buffer = new int*[MPI_processes-1] ;
        int total_size = 0 ;
        for(int i=0;i<MPI_processes-1;++i)
          total_size += recv_size[i] ;
        recv_buffer[0] = new int[total_size] ;
        
	for(int i = 1; i < MPI_processes-1; ++i)
	  recv_buffer[i] = recv_buffer[i-1]+recv_size[i-1] ;
        
	recv_request = new MPI_Request[MPI_processes-1] ;
	status = new MPI_Status[MPI_processes-1] ;
	
	for(k = 0; k < MPI_processes-1; k++) 
	  MPI_Irecv(&recv_buffer[k][0], recv_size[k],MPI_INT, k+1,2, MPI_COMM_WORLD, &recv_request[k] );  

#ifdef DEBUG
	err =
#endif
          MPI_Waitall(MPI_processes-1, recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
	for(k = 0; k < MPI_processes-1; ++k)       
	  for(int i = 0 ; i < recv_size[k]; ++i) {
	    FATAL(re.inSet(recv_buffer[k][i])) ;
	    re += recv_buffer[k][i] ;
	  }
        delete [] recv_buffer[0] ;
	delete [] recv_buffer ;
	delete [] recv_size ;
        delete[] size_request ;
        delete[] size_status ;
        
        delete[] recv_request ;
        delete[] status ;
      } else {
	int *send_buffer;
	int send_size ;
	entitySet temp = e & d->my_entities ;
	send_size = temp.size() ;
	send_buffer = new int[send_size] ;
	int j = 0 ;
	for(ti = temp.begin(); ti != temp.end(); ++ti) {
	  send_buffer[j] = l2g[*ti] ; 
	  re += send_buffer[j] ;
	  ++j ;
	}
	MPI_Send(&send_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD) ;
	MPI_Send(&send_buffer[0], send_size, MPI_INT, 0, 2, MPI_COMM_WORLD) ;
	delete [] send_buffer ;
      }
    }
    return re ;
  }      

  using std::vector ;
  void fill_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {

    vector<entitySet> recv_req(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      if(i!=MPI_rank) 
        recv_req[i] = out_of_dom & init_ptn[i] ;

    // send the recieve requests
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    for(int i=0;i<MPI_processes;++i)
      send_count[i] = recv_req[i].num_intervals() * 2 ;
    MPI_Alltoall(send_count,1,MPI_INT, recv_count, 1, MPI_INT,MPI_COMM_WORLD) ;

    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }
    int mp = MPI_processes-1 ;
    int send_sizes = send_displacement[mp]+send_count[mp] ;
    int recv_sizes = recv_displacement[mp]+recv_count[mp] ;
    int * send_set_buf = new int[send_sizes] ;
    int * recv_set_buf = new int[recv_sizes] ;
    for(int i=0;i<MPI_processes;++i) {
      for(int j=0;j<recv_req[i].num_intervals();++j) {
        send_set_buf[send_displacement[i]+j*2  ] = recv_req[i][j].first ;
        send_set_buf[send_displacement[i]+j*2+1] = recv_req[i][j].second ;
      }
    }

    MPI_Alltoallv(send_set_buf, send_count, send_displacement , MPI_INT,
		  recv_set_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;

    vector<entitySet> send_set(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i) {
      for(int j=0;j<recv_count[i]/2;++j) {
        int i1 = recv_set_buf[recv_displacement[i]+j*2  ] ;
        int i2 = recv_set_buf[recv_displacement[i]+j*2+1] ;
        send_set[i] += interval(i1,i2) ;
      }
    }
    delete[] recv_set_buf ;
    delete[] send_set_buf ;
    
#ifdef DEBUG
    // Sanity check, no send set should be outside of entities we own
    entitySet mydom = init_ptn[MPI_rank] ;
    for(int i=0;i<MPI_processes;++i) {
      if((send_set[i] & mydom) != send_set[i]) {
        cerr << "problem with partitioning in fill_clone!" ;
        debugout << "send_set["<< i << "] = " << send_set[i]
                 << "not owned = " << (send_set[i]-mydom) << endl ;
      }
    }
#endif
    
    // Now that we know what we are sending and receiving (from recv_req and
    // send_set) we can communicate the actual information...

    // Compute sizes of sending buffers
    for(int i=0;i<MPI_processes;++i) 
      send_count[i] =  sp->pack_size(send_set[i]) ;

    // Get sizes needed for receiving buffers
    MPI_Alltoall(send_count,1,MPI_INT, recv_count, 1, MPI_INT,MPI_COMM_WORLD) ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }

    send_sizes = send_displacement[mp]+send_count[mp] ;
    recv_sizes = recv_displacement[mp]+recv_count[mp] ;

    unsigned char *send_store = new unsigned char[send_sizes] ;
    unsigned char *recv_store = new unsigned char[recv_sizes] ;

    for(int i=0;i<send_sizes;++i)
      send_store[i] = 0 ;
    

    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_displacement[i]], loc_pack, send_count[i],
               send_set[i]) ;
    }
    
    MPI_Alltoallv(send_store, send_count, send_displacement, MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;

    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;
      sp->unpack(&recv_store[recv_displacement[i]], loc_pack,recv_count[i],
                 sequence(recv_req[i])) ; 
    }
    delete[] recv_store ;
    delete[] send_store ;

    delete[] recv_count ;
    delete[] send_count ;
    delete[] send_displacement ;
    delete[] recv_displacement ;
  }
  
  storeRepP send_clone_non( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    for(int i = 0; i <  MPI_processes; ++i)
      tmp_dom += entitySet(recv_dom[i]) ;
    storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] =  tmp_sp->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }

  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec( MPI_processes) ;
    std::vector< storeRepP> tmp_sp( MPI_processes) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      recv_count[i] =  tmp_sp[i]->pack_size(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      tmp_sp[i]->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }

  // This only used by the code for the min2noslip computation!  Is it really
  // needed?
#define MINNOSLIP
#ifdef MINNOSLIP
  
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	hm[recv_map[j]] = recv_map[j+1];
	j++ ;
      }
    }
    dMap tmp_dm ;
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) {
      tmp_dm[hmi->first] = hmi->second ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_dm ;
  }
  
  std::vector<dMap> send_global_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    std::vector<HASH_MAP(int, int) > hm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm ;
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	tmp_hm[recv_map[j]] = recv_map[j+1];
	j++ ;
      }
      hm[i] = tmp_hm ;
    }
    std::vector<dMap> v_dm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm = hm[i] ; 
      dMap tmp_dm ;
      for(HASH_MAP(int, int)::const_iterator hmi = tmp_hm.begin(); hmi != tmp_hm.end(); ++hmi) 
	tmp_dm[hmi->first] = hmi->second ;
      v_dm[i] = tmp_dm.Rep() ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return v_dm ;
  }
#endif

  // Collect entitities to a unified entitySet that is distributed across
  // processors according to the partition ptn.
  entitySet dist_collect_entitySet(entitySet inSet, const vector<entitySet> &ptn) {
    const int p = MPI_processes ;
    vector<entitySet> distSet(p) ;
    for(int i=0;i<p;++i) 
      distSet[i] = inSet & ptn[i] ;

    vector<int> send_sz(p) ;
    for(int i=0;i<p;++i)
      if(i!=MPI_rank)
        send_sz[i] = distSet[i].num_intervals()*2 ;
      else
        send_sz[i] = 0 ;
    
    vector<int> recv_sz(p) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<p;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }

    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  p; ++i)
      if(i != MPI_rank)
        for(int j=0;j<distSet[i].num_intervals();++j) {
          send_store[send_displacement[i]+j*2] = distSet[i][j].first ;
          send_store[send_displacement[i]+j*2+1] = distSet[i][j].second ;
        }
    
    
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  

    entitySet retval = distSet[MPI_rank] ;
    for(int i = 0; i <  p; ++i) 
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        retval += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return retval ;
  }
  
  entitySet dist_expand_entitySet(entitySet inSet, entitySet copy,
                                  const vector<entitySet> &ptn) {
    vector<int> send_req(ptn.size()) ;
    for(size_t i=0;i < ptn.size();++i)
      if((copy&ptn[i]) != EMPTY)
        send_req[i] = 1 ;
      else
        send_req[i] = 0 ;
    vector<int> recv_req(ptn.size()) ;
    MPI_Alltoall(&send_req[0],1,MPI_INT, &recv_req[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    vector<int> send_sizes(ptn.size()) ;
    int send_intervals = inSet.num_intervals() ;
    
    for(size_t i=0;i < ptn.size();++i)
      if(recv_req[i]!=0)
        send_sizes[i] = send_intervals ;
      else
        send_sizes[i] = 0 ;
    vector<int> recv_sizes(ptn.size()) ;
    MPI_Alltoall(&send_sizes[0],1,MPI_INT, &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    vector<vector<int> > recv_buffers ;
    vector<int> buf_size ;
    vector<int> proc ;
    vector<MPI_Request> recv_Requests ;
    for(size_t i=0;i < ptn.size();++i) {
      if(recv_sizes[i] != 0) {
        int recv_size = recv_sizes[i]*2 ;
        recv_buffers.push_back(vector<int>(recv_size))  ;
        buf_size.push_back(recv_sizes[i]) ;
        recv_Requests.push_back(MPI_Request()) ;
        proc.push_back(i) ;
      }
    }
    for(size_t i=0;i < recv_buffers.size();++i) {
      int recv_size = buf_size[i]*2 ;
      
      MPI_Irecv(&recv_buffers[i][0],recv_size,MPI_INT, proc[i],3,
                MPI_COMM_WORLD,&recv_Requests[i]) ;
    }

    vector<int> send_buf(send_intervals*2) ;
    for(int j=0;j < send_intervals;++j) {
      send_buf[j*2] = inSet[j].first ;
      send_buf[j*2+1] = inSet[j].second ;
    }

    for(size_t i=0;i < ptn.size();++i) {
      if(send_sizes[i] != 0) {
        MPI_Send(&send_buf[0],send_intervals*2,MPI_INT,i,3,
                 MPI_COMM_WORLD) ;
      }
    }
    for(size_t i=0;i<recv_Requests.size();++i) {
      MPI_Status stat ;
      MPI_Wait(&recv_Requests[i], &stat ) ;
    }
    entitySet recvSet ;
    for(size_t i=0;i<recv_buffers.size();++i)
      for(int j=0;j<buf_size[i];++j) {
        recvSet += interval(recv_buffers[i][j*2],recv_buffers[i][j*2+1]) ;
      }

    inSet += recvSet ;
    return inSet ;
  }
    
    
  entitySet all_collect_entitySet(const entitySet &e) {
    entitySet collect ;
    if(MPI_processes > 1) {
      int *recv_count = new int[MPI_processes] ;
      int *send_count = new int[MPI_processes] ;
      int *send_displacement = new int[MPI_processes];
      int *recv_displacement = new int[MPI_processes];
      int size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i) {
	send_count[i] = 2 * e.num_intervals() ;
	size_send += send_count[i] ; 
      }
      int *send_buf = new int[size_send] ;
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ; 
      int size_recv = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	size_recv += recv_count[i] ;
      int *recv_buf = new int[size_recv] ;
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i)
	for(int j = 0; j < e.num_intervals(); ++j) {
	  send_buf[size_send] = e[j].first ;
	  ++size_send ;
	  send_buf[size_send] = e[j].second ;
	  ++size_send ;
	}
      send_displacement[0] = 0 ;
      recv_displacement[0] = 0 ;
      for(int i = 1; i < MPI_processes; ++i) {
	send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
	recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
      }
      MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		    recv_buf, recv_count, recv_displacement, MPI_INT,
		    MPI_COMM_WORLD) ;

      // Sort intervals to make union operation fast
      vector<interval> ivl_list ;
      for(int i=0;i<size_recv;i+=2) 
        ivl_list.push_back(interval(recv_buf[i],recv_buf[i+1])) ;
      std::sort(ivl_list.begin(),ivl_list.end()) ;
      // Union all of the intervals into the final entitySet
      for(size_t i=0;i<ivl_list.size();++i)
        collect += ivl_list[i] ;

      delete [] send_count ;
      delete [] recv_count ;
      delete [] recv_displacement ;
      delete [] send_displacement ;
      delete [] send_buf ;
      delete [] recv_buf ;
    }
    else
      return e ;
    return collect ;
  }


  std::vector<entitySet> all_collect_vectors(entitySet &e,MPI_Comm comm) {
    int np = 1 ;
    MPI_Comm_size(comm,&np) ;
    int rnk = 0 ;
    MPI_Comm_rank(comm,&rnk) ;
    int *recv_count = new int[ np ] ;
    int *send_count = new int[ np ] ;
    int *send_displacement = new int[ np ];
    int *recv_displacement = new int[ np ];
    int size_send = 0 ;

    for(int i = 0; i <  np; ++i) {
      send_count[i] = e.num_intervals()*2 ;
      size_send += send_count[i] ; 
    }  
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  np ; ++i)
      size_send += recv_count[i] ;
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  np ; ++i) {
      for(int j=0;j<e.num_intervals();++j) {
        send_buf[size_send++] = e[j].first ;
        send_buf[size_send++] = e[j].second ;
      }
    }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  np ; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    } 
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    std::vector<entitySet> vset( np ) ;
    
    for(int i = 0; i < np ; ++i) {
      int *buf = recv_buf+recv_displacement[i] ;
      for(int j=0;j<recv_count[i];j+=2) {
        vset[i] += interval(buf[j],buf[j+1]) ;
      }
    } 
    delete [] send_count ;
    delete [] recv_count ;
    delete [] recv_displacement ;
    delete [] send_displacement ;
    delete [] send_buf ;
    delete [] recv_buf ;
    return vset ;
  }
  std::vector<entitySet> all_collect_vectors(entitySet &e) {
    return all_collect_vectors(e,MPI_COMM_WORLD) ;
  }

  int GLOBAL_OR(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LOR,MPI_COMM_WORLD) ;
    return result ;
  }
  
  int GLOBAL_AND(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LAND,MPI_COMM_WORLD) ;
    return result ;
  }
  
  int GLOBAL_MAX(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD) ;
    return result ;
  }

  int GLOBAL_MIN(int b) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD) ;
    return result ;
  }

} 
