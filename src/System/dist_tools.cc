//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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

#include <set>
using std::set;

#include <string>
using std::string ;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <algorithm>
using std::sort;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>
#include "loci_globs.h"
#include <execute.h>

namespace Loci {
  namespace {
    /*This routine loops over all the rules in the database and extracts
      all the variables associated with the mappings in the head, body and
      the constraints of the rules. */
    //rule_part = 0: Output and Input of the rule mapping
    //          = 1: Only Outupt
    //          = 2: Only Input
    void get_mappings(const rule_db &rdb, fact_db &facts,
		      set<vector<variableSet> > &maps_ret, unsigned int rule_part) {
      ruleSet rules = rdb.all_rules("main") ;
      set<vector<variableSet> > maps ;

      for(ruleSet::const_iterator ri = rules.begin(); ri != rules.end(); ++ri) {
	set<vmap_info>::const_iterator vmsi ;
	if(rule_part != 2) {
	  for(vmsi = ri->get_info().desc.targets.begin();
	      vmsi != ri->get_info().desc.targets.end();
	      ++vmsi) {
	    if(vmsi->mapping.size() != 0) {
	      vector<variableSet> vvs ;
	      for(size_t i = 0; i < vmsi->mapping.size(); ++i) {
		variableSet v ;
		for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
		    vi != vmsi->mapping[i].end();
		    ++vi) {
		  v += variable(*vi,time_ident()) ;
		}
		vvs.push_back(v) ;
	      }
	      maps.insert(vvs) ;
	    }
	  }
	}
	if(rule_part != 1) {
	  for(vmsi = ri->get_info().desc.sources.begin();
	      vmsi != ri->get_info().desc.sources.end();
	      ++vmsi) {
	    //	  if(vmsi->mapping.size() != 0) {
	    vector<variableSet> vvs ;
	    for(size_t i = 0; i < vmsi->mapping.size(); ++i) {
	      variableSet v ;
	      for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
		  vi != vmsi->mapping[i].end();
		  ++vi) {
		v += variable(*vi,time_ident()) ;
	      }
	      vvs.push_back(v) ;
	    }
            // we need to also check to see if there's
            // any maps in "var" field. This is currently just
            // a hack to make the input of a chain of maps
            // work correctly. e.g., we have a rule like this:
            // A <- B->C->D, E. "B, C, and D" are all Maps,
            // currently, when building the rule, only
            // "B->C" is recorded in the vmap_info's mapping
            // field and the "D" is recorded in the "var"
            // field in the vmap_info record. In parallel run,
            // this will cause that the Map "D" not being
            // expanded and will lead to later problems.
            // Therefore, here we also need to check into
            // the "var" field to see if there's any Map
            // inside, in case of yes, we need to collect
            // them as Maps so that they are also expanded.
            //
            // NOTE: this is only a hack, to fix this
            // elegantly and correctly, we need to modify the
            // rule building procedures. Also, recording a
            // Map in the "var" field in the vmap_info structure
            // is probably not a very good idea and may not
            // be entirely safe to other analyses. So we
            // definitely need to fix this in the future
            // if we intend to support a chain of Maps
            // as rule input.
            variableSet vs ;
            for(variableSet::const_iterator vi=vmsi->var.begin();
                vi!=vmsi->var.end();++vi) {
              storeRepP srp = facts.get_variable(*vi) ;
              if(srp == 0)
                continue ;
              if(srp->RepType() == MAP)
                vs += variable(*vi,time_ident()) ;
            }
            if(vs != EMPTY)
              vvs.push_back(vs) ;

            if(!vvs.empty())
              maps.insert(vvs) ;
            //	  }
	  }

	  for(vmsi = ri->get_info().desc.constraints.begin();
	      vmsi != ri->get_info().desc.constraints.end();
	      ++vmsi) {
	    vector<variableSet> vvs ;
	    if(vmsi->mapping.size() != 0) {
	      for(size_t i = 0; i < vmsi->mapping.size(); i++) {
		variableSet v ;
		for(variableSet::const_iterator vi = vmsi->mapping[i].begin();
		    vi != vmsi->mapping[i].end();
		    ++vi) {
		  v += variable(*vi,time_ident()) ;
		}
		vvs.push_back(v) ;
	      }
	      maps.insert(vvs) ;
	    }
	  }
	}
      }

      set<vector<variableSet> >::const_iterator mi ;

      variableSet vars = facts.get_typed_variables() ;
      maps_ret.clear() ;
      for(mi=maps.begin();mi!=maps.end();++mi) {
	const vector<variableSet> &vss = *mi ;
	// If the maps aren't in the fact database then exclude it
	variableSet notvars ;
	for(size_t i=0;i<vss.size();++i)
	  notvars += vss[i]-vars ;
	if(notvars != EMPTY)
	  continue ;

	// Now expand mapping of length up two 2
	if(vss.size() < 2) {
	  for(variableSet::const_iterator vi=vss[0].begin();vi!=vss[0].end();++vi) {
	    vector<variableSet> vs(1) ;
	    vs[0] += *vi ;
	    maps_ret.insert(vs) ;
	  }
	} else {
	  for(variableSet::const_iterator vi=vss[0].begin();
	      vi!=vss[0].end();++vi)
	    for(variableSet::const_iterator vvi=vss[1].begin();
		vvi!=vss[1].end();++vvi) {
	      vector<variableSet> vs(vss.size()) ;
	      vs[0] += *vi ;
	      vs[1] += *vvi ;
	      for(size_t i=2;i<vss.size();++i)
		vs[i] = vss[i] ;
	      maps_ret.insert(vs) ;
	    }
	}

      }

    }

    //This routine  is similar to the expand map but it works for maps
    //which are distributed across processors.
    vector<entitySet> dist_expand_map(vector<entitySet> domain, fact_db &facts,
			      const std::set<std::vector<variableSet> > &maps) {

      vector<entitySet> dom = domain ;
      variableSet vars = facts.get_typed_variables() ;
      std::set<std::vector<variableSet> >::const_iterator smi ;
      for(smi = maps.begin(); smi != maps.end(); ++smi) {
	vector<entitySet> locdom = domain ;
	const vector<variableSet> &mv = *smi ;
	      
	for(size_t i = 0; i < mv.size(); ++i) {
	  variableSet v = mv[i] ;
	  v &= vars ;
	  vector<entitySet> image(domain.size()) ;

	  for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	    storeRepP p = facts.get_variable(*vi) ;
	    if(p->RepType() ==  MAP) {      
	      int kd = p->getDomainKeySpace() ;
	      std::vector<entitySet> ptn = facts.get_init_ptn(kd) ;

	      entitySet tmp_dom = p->domain() ;
	      MapRepP mp =  MapRepP(p->getRep()) ;
	      entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	      entitySet tmp_out = (glob_dom & locdom[kd]) - tmp_dom ;
	      storeRepP sp = mp->expand(tmp_out, ptn) ;
	      if(sp->domain() != tmp_dom) {
		//facts.update_fact(variable(*vi), sp) ;
		sp->setDomainKeySpace(kd) ;
		MapRepP(sp)->setRangeKeySpace(mp->getRangeKeySpace()) ;
		facts.replace_fact(*vi,sp) ;
	      }
	      int ki = mp->getRangeKeySpace() ;
	      entitySet tmp = MapRepP(sp)->image((sp->domain()) & locdom[kd]) ;
	      image[ki] +=  tmp ;

	    }
	  }
	  for(size_t i=0;i<domain.size();++i)
	    dom[i] += image[i] ;
	  locdom = image ;
	}
      }
      return dom ;
    }

    // This is only used by the experimental work replication option
    std::set<std::vector<variableSet> > classify_moderate_maps(fact_db &facts, const std::set<std::vector<variableSet> > &maps) {
      variableSet vars = facts.get_typed_variables() ;
      std::set<std::vector<variableSet> >::const_iterator smi ;
      std::set<std::vector<variableSet> > return_maps;
      for(smi = maps.begin(); smi != maps.end(); ++smi) {
	entitySet domain;
	entitySet image;
	entitySet tmp;
	const vector<variableSet> &mv = *smi ;
	for(unsigned int i = 0; i < mv.size(); i++) {
	  variableSet v = mv[i] ;
	  v &= vars ;
	  if(i == 0) {
	    for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	      storeRepP p = facts.get_variable(*vi) ;
	      if(p->RepType() ==  MAP) {
		MapRepP mp =  MapRepP(p->getRep()) ;
		domain += mp->domain();
	      }
	    }
	    domain = all_collect_entitySet(domain);
	    tmp = domain;
	  }

	  for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	    storeRepP p = facts.get_variable(*vi) ;
	    if(p->RepType() ==  MAP) {
	      MapRepP mp =  MapRepP(p->getRep()) ;
	      image += mp->image(tmp);
	    }
	  }

	  tmp = all_collect_entitySet(image);
	  image = EMPTY;
	}

	image = tmp;
	debugout << "map ratio: " << 1.0*domain.size()/image.size() << endl;
	for(unsigned int i = 0; i < mv.size(); i++) {
	  debugout << "(" << mv[i] << ")  ->  ";
	}
	debugout << endl;
	if(domain.size()*1.0 > image.size() *  10.0) {
	  facts.intensive_output_maps.insert(*smi);
	}
	else {
	  return_maps.insert(*smi);
	}
      }
      return return_maps;
    }


    // This only used for experimental work replication feature.
    // Find sthe context for maps whose final image will be in provided domain.
    entitySet context_for_map_output(entitySet domain, fact_db &facts,
				     const std::set<std::vector<variableSet> > &maps) {
      std::vector<entitySet> ptn = facts.get_init_ptn(0) ;// FIX THIS
      entitySet context;
      variableSet vars = facts.get_typed_variables() ;
      std::set<std::vector<variableSet> >::const_iterator smi ;
      for(smi = maps.begin(); smi != maps.end(); ++smi) {
	std::vector<entitySet>  preimage_vec = all_collect_vectors(domain);
	const vector<variableSet> &mv = *smi ;
	for(int i = mv.size() -1; i >= 0; --i) {
	  variableSet v = mv[i] ;
	  v &= vars ;
	  std::vector<entitySet>  tmp_preimage_vec(MPI_processes);
	  for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {
	    storeRepP p = facts.get_variable(*vi) ;
	    if(p->RepType() ==  MAP) {
	      MapRepP mp =  MapRepP(p->getRep()) ;
	      for(int j = 0; j < MPI_processes; j++) {
		tmp_preimage_vec[j] += mp->preimage(preimage_vec[j]).second;
	      }
	    }
	  }
	  for(int j = 0; j < MPI_processes; j++) {
	    preimage_vec[j] = all_collect_entitySet(tmp_preimage_vec[j]);
	  }

	  if(i == 0) {
	    context += preimage_vec[MPI_rank];
	  }
	}
      }
      return context ;
    }

    // This is only used by the experimental work replication
    // It works just like dist_expand_map except the return entities only contain
    // the image of the last map in the chain. the difference of the returned
    // value is only 1 or 2 entities
    entitySet dist_special_expand_map(entitySet domain, fact_db &facts,
				      const std::set<std::vector<variableSet> > &maps) {
      entitySet special_return;

      // FIX THIS
      std::vector<entitySet> ptn = facts.get_init_ptn(0) ;

      variableSet vars = facts.get_typed_variables() ;
      std::set<std::vector<variableSet> >::const_iterator smi ;
      for(smi = maps.begin(); smi != maps.end(); ++smi) {
	entitySet locdom = domain ;
	const vector<variableSet> &mv = *smi ;
	for(size_t i = 0; i < mv.size(); ++i) {
	  variableSet v = mv[i] ;
	  v &= vars ;
	  entitySet image ;
	  for(variableSet::const_iterator vi = v.begin(); vi != v.end(); ++vi) {

	    storeRepP p = facts.get_variable(*vi) ;
	    if(p->RepType() ==  MAP) {
	      entitySet tmp_dom = p->domain() ;
	      MapRepP mp =  MapRepP(p->getRep()) ;
	      entitySet glob_dom = all_collect_entitySet(tmp_dom) ;
	      entitySet tmp_out = (glob_dom & locdom) - tmp_dom ;
	      storeRepP sp = mp->expand(tmp_out, ptn) ;
	      if(sp->domain() != tmp_dom) {
		//facts.update_fact(variable(*vi), sp) ;
		facts.replace_fact(*vi,sp) ;
	      }
	      image +=  MapRepP(sp)->image((sp->domain()) & locdom) ;

	    }
	  }

	  locdom = image ;
	  if(i == mv.size() -1)
	    special_return += image;
	}
      }
      return special_return ;
    }

  }

  void get_clone(fact_db &facts, const rule_db &rdb) {
    fact_db::distribute_infoP df = facts.get_distribute_info()  ;
    // reallocate constraints on whole domain(instead of local domain) on each process
    variableSet tmp_vars = facts.get_typed_variables();
    for(variableSet::const_iterator vi = tmp_vars.begin(); vi != tmp_vars.end(); ++vi) {
      storeRepP tmp_sp = facts.get_variable(*vi) ;
      if(tmp_sp->RepType() == CONSTRAINT) {
        entitySet tmp_dom = tmp_sp->domain() ;
        if(GLOBAL_OR(tmp_dom != ~EMPTY,MPI_COMM_WORLD)) {
          entitySet global_tmp_dom = all_collect_entitySet(tmp_dom) ;
          constraint tmp ;
          *tmp = global_tmp_dom ;
	  tmp.Rep()->setDomainKeySpace(tmp_sp->getDomainKeySpace()) ;
          facts.replace_fact(*vi,tmp) ;
        }
      }
      if(tmp_sp->RepType() == MAP) {
	MapRepP tmp_mp = MapRepP(tmp_sp->getRep()) ;
	int kd_domain = tmp_sp->getDomainKeySpace() ;
	int kd_range = tmp_mp->getRangeKeySpace() ;
        storeRepP map_sp = MapRepP(tmp_sp->getRep())->thaw() ;
	map_sp->setDomainKeySpace(kd_domain) ;
	tmp_mp = MapRepP(map_sp->getRep()) ;
	tmp_mp->setRangeKeySpace(kd_range) ;
        facts.replace_fact(*vi, map_sp) ;
      }
    }

    //find context. with duplication, it's the preimage of my_entites of output mapping .
    //(if extended_duplication, the mapping includes both input and output)
    //(if multilevel_duplication, the context includes mapping back and forth n times)
    //without duplication, the context is my_entities
    int nkd = facts.getNumKeyDomains() ;
    vector<entitySet> tmp_set(nkd);
    // This part will be broken if the keyDomains feature is used

    if(duplicate_work) {
      std::vector<entitySet> &ptn = facts.get_init_ptn(0) ; // FIX THIS
      std::set<std::vector<variableSet> > context_maps ;

      if(extended_duplication)
	get_mappings(rdb, facts, context_maps, 0);
      else
	get_mappings(rdb, facts, context_maps, 1);

      context_maps = classify_moderate_maps(facts, context_maps);

      //Find out entities that can produce target variable entities owned by a processor
      facts.global_comp_entities += context_for_map_output(ptn[MPI_rank], facts, context_maps);

      //Add entities so that maximum depth of duplication can be achieved
      if(multilevel_duplication) {
	entitySet mySet = ptn[MPI_rank];
	bool continue_adding = true;
	int num_levels = 1;
	do{
	  mySet += dist_special_expand_map(facts.global_comp_entities,
					   facts, context_maps) ;
	  entitySet added_entities = context_for_map_output(mySet,  facts, context_maps);
	  added_entities -= facts.global_comp_entities;
	  if(all_collect_entitySet(added_entities) == EMPTY)
	    continue_adding = false;
	  else {
	    facts.global_comp_entities += added_entities;
	    num_levels++;
	  }
	}while(continue_adding);
	facts.global_comp_entities += mySet;
#ifdef VERBOSE
	debugout << "Number of Duplication Levels: " << num_levels << endl;
#endif
      }
      //the union of context and my entities is comp_entities
      tmp_set[0] = facts.global_comp_entities;
    }
    else {
      for(int i=0;i<nkd;++i) {
	tmp_set[i] = facts.get_init_ptn(i)[MPI_rank] ;
      }
    }
    //The clone region is the image and domain of each map in both input and output on comp_entities
    //minus my_entities
    std::set<std::vector<variableSet> > dist_maps ;
    get_mappings(rdb,facts,dist_maps, 0) ;


    double clone_time_start = MPI_Wtime() ;

    vector<entitySet> image =dist_expand_map(tmp_set, facts, dist_maps) ;

    vector<vector<entitySet> > copyv(nkd), send_clonev(nkd) ;
    for(int kd=0;kd<nkd;++kd) {
      vector<entitySet>  copy(MPI_processes), send_clone(MPI_processes) ;
      entitySet tmp_copy ;
      std::vector<entitySet> &ptn = facts.get_init_ptn(kd) ; 
      tmp_copy =  image[kd] - ptn[MPI_rank] ;
      
      int *recv_count = new int[MPI_processes] ;
      int *send_count = new int[MPI_processes] ;
      int *send_displacement = new int[MPI_processes];
      int *recv_displacement = new int[MPI_processes];
      int size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i) {
	copy[i] = tmp_copy & ptn[i] ;
	send_count[i] = copy[i].size() ;
	size_send += send_count[i] ;
      }
      
      int *send_buf = new int[size_send] ;
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		   MPI_COMM_WORLD) ;
      size_send = 0 ;
      for(int i = 0; i < MPI_processes; ++i) {
	size_send += recv_count[i] ;
      }
      int *recv_buf = new int[size_send] ;
      size_send = 0 ;
      entitySet::const_iterator ei ;
      for(int i = 0; i < MPI_processes; ++i)
	for(ei = copy[i].begin(); ei != copy[i].end(); ++ei) {
	  send_buf[size_send] = *ei ;
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
      std::vector<entitySet> add(MPI_processes) ;
      for(int i = 0; i < MPI_processes; ++i) {
	for(int j = recv_displacement[i]; j <
	      recv_displacement[i]+recv_count[i]; ++j)
	  send_clone[i] += recv_buf[j] ;
      }
      delete [] recv_count ;
      delete [] send_count ;
      delete [] send_displacement ;
      delete [] recv_displacement ;
      delete [] send_buf ;
      delete [] recv_buf ;
      copyv[kd] = copy ;
      send_clonev[kd] = send_clone ;
    }
    double clone_time_end = MPI_Wtime() ;
    debugout << "  Time taken for creating clone info =  "
             << clone_time_end - clone_time_start << endl ;

    variableSet vars = facts.get_typed_variables() ;
    double start = MPI_Wtime() ;
    //    timer_token creating_initial_info_timer = new timer_token;
    //	if(collect_perf_data)
    //		creating_initial_info_timer = perfAnalysis->start_timer("Create initial info");
    int myid = MPI_rank ;
    int size = 0 ;
    Map l2g ;
    int isDistributed ;

    vector<vector<entitySet> > proc_entities(nkd) ;
    
    for(int kd=0;kd<nkd;++kd) {
      stopWatch s ;
      s.start() ;
      std::vector<entitySet> iv ;
      categories(facts,iv,kd) ; // FIX THIS
      std::vector<entitySet> &ptn = facts.get_init_ptn(kd) ; 
      debugout << "finding categories time = " << s.stop() << endl ;
      s.start() ;
      for(size_t i=0;i < iv.size(); ++i) {
	//      iv[i] = all_collect_entitySet(iv[i]) ;
	debugout << "iv["<< i << "] size before expand = " << iv[i].num_intervals()<< endl ;
	entitySet tmp_copy =  image[kd] - ptn[MPI_rank] ;
	iv[i] = dist_expand_entitySet(iv[i],tmp_copy,ptn) ;
	debugout << "iv["<< i << "] size after expand = " << iv[i].num_intervals() <<endl;
      }
      debugout << "time expanding categories = " << s.stop() << endl ;
      s.start() ;
      entitySet::const_iterator ti ;
      entitySet e ;

      debugout << " initial_categories =  " << iv.size() << endl ;

      for(size_t i = 0; i < iv.size(); ++i) {
	//    debugout << " iv[ " << i << " ] = " << iv[i] << endl ;
	// Within each category:
	// 1) Number local processor entities first
	e = ptn[myid] & iv[i] ;
	if(e != EMPTY) {
	  proc_entities[kd].push_back(e) ;
	  size += e.size() ;
	}
	// 2) Number clone region entities next
	for(int j = 0; j < MPI_processes; ++j)
	  if(myid != j) {
	    e = copyv[kd][j] & iv[i];
	    if(e != EMPTY) {
	      proc_entities[kd].push_back(e) ;
	      size += e.size() ;
	    }
	  }
      }
      iv.clear() ;
    }

    int j = 0 ;
    entitySet e = interval(0, size-1) ;
    l2g.allocate(e) ;
    store<unsigned char> key_domain ;
    key_domain.allocate(e) ;

    stopWatch s ;
    s.start() ;
    vector<entitySet> glist ;
    for(int kd=0;kd<nkd;++kd) {
      entitySet g ;
      for(size_t i = 0; i < proc_entities[kd].size(); ++i) {
	g += proc_entities[kd][i] ;
	entitySet::const_iterator ei ;
	for(ei = proc_entities[kd][i].begin(); ei != proc_entities[kd][i].end(); ++ei ) {
	  l2g[j] = *ei ;
	  key_domain[j] = kd ;
	  ++j ;
	}
      }
      proc_entities[kd].clear() ;
      glist.push_back(g) ;
    }
    debugout << "time setting up l2g = " << s.stop() ;
    s.start() ;
    df->l2g = l2g.Rep() ;
    df->key_domain = key_domain.Rep() ; // FIX THIS

    vector<dMap> tmp2(nkd) ;
    df->g2lv=tmp2 ;
    for(int kd=0;kd<nkd;++kd) {
      dMap tmp ;
      df->g2lv[kd].setRep(tmp.Rep()) ;
      df->g2lv[kd].allocate(glist[kd]) ;
    }
    entitySet ldom = l2g.domain() ;
    entitySet::const_iterator ei ;
    for(entitySet::const_iterator ei=ldom.begin();ei!=ldom.end();++ei) {
      int kd = key_domain[*ei] ;
      df->g2lv[kd][l2g[*ei]] = *ei ;
    }


    entitySet send_neighbour ;
    entitySet recv_neighbour ;
    store<entitySet> send_entities ;
    store<entitySet> recv_entities ;
    for(int i = 0 ; i < MPI_processes; ++i)
      if(myid != i )
	for(int kd=0;kd<nkd;++kd) 
	  if(copyv[kd][i] != EMPTY)
	    recv_neighbour += i ;

    for(int i = 0; i < MPI_processes; ++i)
      if(myid != i)
	for(int kd=0;kd<nkd;++kd) 
	  if(send_clonev[kd][i] != EMPTY)
	    send_neighbour += i ;

    send_entities.allocate(send_neighbour) ;
    recv_entities.allocate(recv_neighbour) ;
    entitySet::const_iterator ti ;
    for(ei = recv_neighbour.begin(); ei != recv_neighbour.end(); ++ei) {
      for(int kd=0;kd<nkd;++kd) 
	for(ti =  copyv[kd][*ei].begin(); ti != copyv[kd][*ei].end(); ++ti) {
	  recv_entities[*ei] += df->g2lv[kd][*ti] ;
	}
    }

    for(ei = send_neighbour.begin(); ei!= send_neighbour.end(); ++ei) {
      for(int kd=0;kd<nkd;++kd) 
	for(ti =  send_clonev[kd][*ei].begin(); ti != send_clonev[kd][*ei].end(); ++ti)
	  send_entities[*ei] +=  df->g2lv[kd][*ti] ;
    }


    debugout << "time setting up send and recieve info = " << s.stop() << endl ;
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(creating_initial_info_timer);
    double end_time =  MPI_Wtime() ;
    debugout << "  Time taken for creating initial info =  " << end_time - start << endl ;

    start = MPI_Wtime() ;
    //    timer_token reordering_timer = new timer_token;
    //	if(collect_perf_data)
    //		reordering_timer = perfAnalysis->start_timer("Reordering");
    //----------------------------------------------------------------------
    // reorder facts to local numbering
    reorder_facts(facts, df->g2lv[0]) ; 

    //----------------------------------------------------------------------
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(reordering_timer);
    end_time =  MPI_Wtime() ;
    debugout << "  Time taken for reordering =  " << end_time - start << endl ;

    isDistributed = 1 ;
    df->isDistributed = isDistributed ;
    entitySet mySetLocal ;
    for(int kd=0;kd<nkd;++kd)  {
      entitySet g ;
      std::vector<entitySet> &ptn = facts.get_init_ptn(kd) ; 

#ifdef DEBUG
      entitySet g2ldom = df->g2lv[kd].domain() ;
      if(ptn[myid]-g2ldom != EMPTY) {
	cerr << "problem with g2lv " << ptn[myid]-g2ldom << endl ;
      }
#endif
      for(ei = ptn[myid].begin(); ei != ptn[myid].end(); ++ei) {
	int gv = df->g2lv[kd][*ei] ;
	if(g.inSet(gv)) 
	  cerr << "repeated values in gv!" << endl ;
	g += gv ;
      }
      if((g&mySetLocal) != EMPTY)
	cerr << "intersection between mySetLocals from keyspaces!" << endl ;
      mySetLocal += g ;
    }
    //Add comp_entities
    if(duplicate_work) {
      entitySet g = EMPTY;
      for(ei = facts.global_comp_entities.begin();
	  ei != facts.global_comp_entities.end(); ++ei)
	g += df->g2lv[0][*ei] ;
      df->comp_entities = g;
    }
    /*xmit data structure contains the information as to what
      entities are to be send to what processor . The copy data
      structure contains the entities that are to be received from a
      particular processor(information regarding the clone region
      entities). All the entities are stored in their local
      numbering. A local to global numbering l2g  is provided to send the
      entities in their original global numbering.*/
    for(ei=send_neighbour.begin(); ei != send_neighbour.end();++ei)
      df->xmit.push_back
        (fact_db::distribute_info::dist_data(*ei,send_entities[*ei])) ;
    for(ei=recv_neighbour.begin(); ei != recv_neighbour.end();++ei)
      df->copy.push_back
        (fact_db::distribute_info::dist_data(*ei,recv_entities[*ei])) ;

    df->myid = myid ;
    constraint my_entities ;
    my_entities = mySetLocal ;
    df->my_entities = mySetLocal ;

    int total = 0 ;

    for(size_t i=0;i<df->xmit.size();++i)
      total += df->xmit[i].size ;
    df->xmit_total_size = total ;

    total = 0 ;
    for(size_t i=0;i<df->copy.size();++i)
      total += df->copy[i].size ;
    df->copy_total_size = total ;
    facts.put_distribute_info(df) ;
    // this needs to be an intensional fact
    facts.create_intensional_fact("my_entities",my_entities);

  }

  /*! The fill_entitySet routine fills in the clone region entities
    . The send_buffer and the recv_buffer are allocated only once to
    contain the maximum clone region

    arguments: e: entitySet that if it's in my xmit region, I need send it to others
    returned: entitySet I received
  */

  entitySet fill_entitySet(const entitySet& e, fact_db &facts) {

    entitySet re ;

    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;
        recv_buffer[0] = new int[2*d->copy_total_size] ;
        recv_size[0] = 2*d->copy[0].size ;
        for(size_t i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+2*d->copy[i-1].size ;
	  recv_size[i] = 2*d->copy[i].size ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[2*d->xmit_total_size] ;
        for(size_t i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+2*d->xmit[i-1].size ;
      }


      Map l2g ;
      l2g = d->l2g.Rep() ;
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;
      /*The recv_size is the maximum possible, so that even if we
	receive a shorter message there won't be any problem */
      for(size_t i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(size_t i=0;i<d->xmit.size();++i) {
        entitySet temp = e & d->xmit[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	  send_buffer[i][j++] = key_domain[*ei] ;
          send_buffer[i][j++] = l2g[*ei] ;
	}

	int send_size = 2*temp.size() ;
        MPI_Send(send_buffer[i], send_size, MPI_INT, d->xmit[i].proc,
                 1, MPI_COMM_WORLD) ;
      }


      if(d->copy.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }
      entitySet tst ;
      for(size_t i = 0; i < d->copy.size(); ++i) {
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        for(int j = 0 ; j < recieved; ++j) {
	  int kd = recv_buffer[i][j++] ;
          re += d->g2lv[kd][recv_buffer[i][j]] ;
	}
      }


      if(d->copy.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->xmit.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] status ;
      delete [] recv_request ;

    }
    return re ;
  }

  /*This is an optimization to the fill_entitySet routine to which we
    passed only an entitySet. In this case we pass in a vector of
    entitySet so that we can group the communication of entities. This
    avoids the additional start up cost incurred when we send the
    entities corresponding to an entitySet
    each time . ie with one startup cost ts we can send all the
    entities required to a particular processor. */
  /*! ev: the entitySets I send if they are in my xmit region
    return: the entitySets I recieve*/
  vector<entitySet> fill_entitySet(const vector<entitySet>& ev,
                                   fact_db &facts) {

    vector<entitySet> re(ev.size()) ;

    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      const int evsz = ev.size() ;
      int *evtmp = new int[evsz] ;
      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->copy.size() > 0) {
        recv_buffer = new int*[d->copy.size()] ;
        recv_size = new int[d->copy.size()] ;

        recv_buffer[0] = new int[2*d->copy_total_size*evsz+evsz*d->copy.size()] ;
        recv_size[0] = 2*d->copy[0].size*evsz+evsz ;
        for(size_t i=1;i<d->copy.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(2*d->copy[i-1].size*evsz+evsz) ;
          recv_size[i] = 2*d->copy[i].size*evsz+evsz ;
        }
      }

      if(d->xmit.size() > 0) {
        send_buffer = new int*[d->xmit.size()] ;

        send_buffer[0] = new int[2*d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        for(size_t i=1;i<d->xmit.size();++i)
          send_buffer[i] = send_buffer[i-1]+(2*d->xmit[i-1].size*evsz+evsz) ;
      }


      Map l2g ;
      l2g = d->l2g.Rep() ;
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->copy.size()] ;
      MPI_Status *status = new MPI_Status[d->copy.size()] ;

      for(size_t i=0;i<d->copy.size();++i)
        MPI_Irecv(recv_buffer[i],recv_size[i],MPI_INT,d->copy[i].proc,1,
                  MPI_COMM_WORLD, &recv_request[i]) ;

      for(size_t i=0;i<d->xmit.size();++i) {
        int j=evsz ;
        for(int k=0;k<evsz;++k) {
          entitySet temp = ev[k] & d->xmit[i].entities ;
          send_buffer[i][k] = temp.size() ;

          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	    send_buffer[i][j++] = key_domain[*ei] ;
            send_buffer[i][j++] = l2g[*ei] ;
	  }
        }
        int send_size = j ;
        MPI_Send(send_buffer[i], send_size, MPI_INT, d->xmit[i].proc,
                 1, MPI_COMM_WORLD) ;
      }


      if(d->copy.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->copy.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }

      for(size_t i = 0; i < d->copy.size(); ++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=evsz ;
        WARN(recieved < evsz) ;
        for(int k=0;k<evsz;++k) {
          for(int l=0;l<recv_buffer[i][k];++l) {
	    int kd = recv_buffer[i][j++] ;
       	    re[k] += d->g2lv[kd][recv_buffer[i][j++]] ;
	  }
        }
        WARN(j!=recieved) ;
      }


      if(d->copy.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->xmit.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] evtmp ;
      delete [] status ;
      delete [] recv_request ;

    }
    return re ;
  }

  //===================================================================
  //
  // The send_entitySet routine is used to handle cases when there are
  // mapping in the output. Sometimes we might compute entities in the
  // clone region. Since these entities are not owned by the processor
  // it needs to be send to the processor that actually owns them.
  //! e: the entitySet I send if it's in my clone region
  //    return: the entitySet I  receive
  // Note, this assumes that the entities are in local number

  entitySet send_entitySet(const entitySet &e, fact_db &facts) {
    vector<entitySet> in ;
    in.push_back(e) ;
    vector<entitySet> r = fill_entitySet(in,facts) ;
    return r[0] ;
  }


  //Similar to send_entitySet defined above.
  //Instead of returning set of entities, it returns
  //vector of entitySet describing exactly which entities
  //are received from which processor.
  vector<entitySet> send_entitySetv(const entitySet& e, fact_db &facts) {
    vector<entitySet> re(MPI_processes);
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[2*d->xmit_total_size] ;
        recv_size[0] = 2*d->xmit[0].size ;

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+2*d->xmit[i-1].size ;
          recv_size[i] = 2*d->xmit[i].size ;
        }
      }

      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[2*d->copy_total_size] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+2*d->copy[i-1].size ;
      }
      Map l2g ;
      l2g = d->l2g.Rep() ;
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;

      /*By intersecting the given entitySet with the clone region
	entities we can find out which entities are to be sent */
      for(size_t i=0;i<d->copy.size();++i) {
        entitySet temp = e & d->copy[i].entities ;

        int j=0 ;
        for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	  send_buffer[i][j++] = key_domain[*ei] ;
          send_buffer[i][j++] = l2g[*ei] ;
	}

        int send_size = temp.size() ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }

      if(d->xmit.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
	FATAL(err != MPI_SUCCESS) ;
      }

      for(size_t i=0;i<d->xmit.size();++i) {
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
        for(int j=0;j<recieved;++j) {
	  int kd = recv_buffer[i][j++] ;
          re[d->xmit[i].proc] += d->g2lv[kd][recv_buffer[i][j]] ;
	}
      }

      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;

    }
    return re ;
  }
  /*! ev: the entittySets that I send if they are in my clone region
    return: the entitySets I receive; its index is the same as that of ev*/
  vector<entitySet> send_entitySet(const vector<entitySet>& ev,
                                   fact_db &facts) {
    vector<entitySet> re(ev.size()) ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;

      const int evsz = ev.size() ;
      int **send_buffer = 0 ;
      int **recv_buffer = 0 ;
      int *recv_size = 0 ;

      if(d->xmit.size() > 0) {
        recv_buffer = new int*[d->xmit.size()] ;
        recv_size = new int[d->xmit.size()] ;

        recv_buffer[0] = new int[2*d->xmit_total_size*evsz+evsz*d->xmit.size()] ;
        recv_size[0] = 2*d->xmit[0].size*evsz + evsz ;

        for(size_t i=1;i<d->xmit.size();++i) {
          recv_buffer[i] = recv_buffer[i-1]+(2*d->xmit[i-1].size*evsz+evsz) ;
          recv_size[i] = 2*d->xmit[i].size*evsz+evsz ;
        }
      }

      if(d->copy.size() > 0 ) {
        send_buffer = new int*[d->copy.size()] ;
        send_buffer[0] = new int[2*d->copy_total_size*evsz+evsz*d->copy.size()] ;
        for(size_t i=1;i<d->copy.size();++i)
          send_buffer[i] = send_buffer[i-1]+2*d->copy[i-1].size*evsz+evsz ;
      }
      Map l2g ;
      l2g = d->l2g.Rep();
      store<unsigned char> key_domain ;
      key_domain = d->key_domain.Rep() ;

      MPI_Request *recv_request = new MPI_Request[d->xmit.size()] ;
      MPI_Status *status = new MPI_Status[d->xmit.size()] ;

      for(size_t i=0;i<d->xmit.size();++i)
	MPI_Irecv(recv_buffer[i], recv_size[i], MPI_INT, d->xmit[i].proc, 1,
                  MPI_COMM_WORLD, &recv_request[i] ) ;

      for(size_t i=0;i<d->copy.size();++i) {
        int j=evsz ;
        for(int k=0;k<evsz;++k) {
          entitySet temp = ev[k] & d->copy[i].entities ;
          send_buffer[i][k] = temp.size() ;

          for(entitySet::const_iterator ei=temp.begin();ei!=temp.end();++ei) {
	    send_buffer[i][j++] = key_domain[*ei] ;
            send_buffer[i][j++] = l2g[*ei] ;
	  }
        }
        int send_size = j ;
        MPI_Send(send_buffer[i],send_size, MPI_INT, d->copy[i].proc,
                 1,MPI_COMM_WORLD) ;
      }

      if(d->xmit.size() > 0) {
#ifdef DEBUG
	int err =
#endif
          MPI_Waitall(d->xmit.size(), recv_request, status) ;
      	FATAL(err != MPI_SUCCESS) ;
      }
      for(size_t i=0;i<d->xmit.size();++i) {
#ifdef DEBUG
        int recieved ;
	MPI_Get_count(&status[i], MPI_INT, &recieved) ;
#endif
        int j=evsz ;
        for(int k=0;k<evsz;++k)
          for(int l=0;l<recv_buffer[i][k];++l) {
	    int kd = recv_buffer[i][j++] ;
            re[k] += d->g2lv[kd][recv_buffer[i][j++]] ;
	  }
        WARN(j!=recieved) ;
      }

      if(d->xmit.size() > 0) {
        delete [] recv_size ;
        delete [] recv_buffer[0] ;
        delete [] recv_buffer ;
      }
      if(d->copy.size() > 0) {
        delete [] send_buffer[0] ;
        delete [] send_buffer ;
      }
      delete [] recv_request ;
      delete [] status ;

    }
    return re ;
  }


  // The motivation of this routine is that in the parallel code,
  // a Map used as constraint will often not be expanded enough to
  // include the clone region, thus causing problems in the rule
  // execution schedule. Since the context of the rule including
  // the clone region
  // will often exceeds the domain of the constraint, this will either
  // cause the clone region not being computed properly, or in the case
  // of unit/apply rule, cause conflicts in the existential analysis.
  //
  // We could expand the Maps used in rule contraints to include the
  // clone region. However doing so would often require duplicating
  // the Map on all processes, thus incurring a memory cost, or else,
  // we would be allocating the Map on domains that do not have
  // meaningful values for the Map.
  rule_db
  replace_map_constraints(fact_db& facts, const rule_db& rdb) {

    rule_db new_rdb ;
    ruleSet all_rules = rdb.all_rules() ;

    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri) {

      rule_implP rp = ri->get_rule_implP() ;
      if(rp == 0) {
        new_rdb.add_rule(*ri) ;
        continue ;
      }

      // save the old rule name
      string old_name = ri->get_info().name() ;

      rp->replace_map_constraints(facts) ;

      rule nr(rp) ;
      // rename the rule to its old name gives the rule
      // signature to still have the original Map name
      // in the constraint field, which is easier for
      // the output routine and for users to read
      nr.rename(old_name) ;

      new_rdb.add_rule(nr) ;
    }
    return new_rdb ;
  }

}

