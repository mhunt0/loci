//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#include <mod_db.h>
#include "loci_globs.h"
#include <iostream>
#include <distribute.h>
#include <vector>
#include <string>
//#define VERBOSE
namespace Loci {
  using std::vector ;
  using std::string ;

  vector<string> ModuleDirectoryPath ;

  void AddModuleSearchDir(string dirname) {
    ModuleDirectoryPath.push_back(dirname) ;
  }

  using std::cout ;
  std::string remove_space(const std::string &str){
    std::string tmp_str ;
    for(std::string::const_iterator si = str.begin(); si != str.end(); ++si) {
      if(*si != ' ')
	tmp_str += *si ;
    }
    return tmp_str ;
  }
  //This is a hack for specifying hierarchy of namespaces. We specify
  // the name to be loaded from as say for eg: chem_heat_.... When we
  // call the load_module and pass in this string it loads all the
  // variables into the namespace chem@heat@ .... etc.   
  void parse_str(const std::string& str, std::vector<std::string> &str_vec) {
    if(!str.empty()) {
      size_t tmp = 0 ;
      std::string sub_str = remove_space(str) ;
      while(tmp != std::string::npos) {
	tmp = sub_str.find(",") ;
	str_vec.push_back(sub_str.substr(0,tmp)) ;
        sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
      }
    }
  }

  void *dlopen_helper(string module) {
    void *tmp ;
    for(size_t i=0;i<ModuleDirectoryPath.size();++i) {
      string test = ModuleDirectoryPath[i] + "/" + module ;
      tmp = dlopen(test.c_str(),RTLD_GLOBAL|RTLD_NOW) ;
      if(tmp != 0) {
        Loci::debugout << "reading modules in from file '" << test <<"'" << endl ;
        return tmp ;
      }
    }
    tmp = dlopen(module.c_str(),RTLD_GLOBAL|RTLD_NOW) ;
    if(tmp == 0) {
      if(Loci::MPI_rank == 0)
        cerr << "unable to open " << module << endl ;
      const char *error = dlerror() ;
      if(error)
        if(Loci::MPI_rank == 0)
          cerr << "reason for failure is " << error << endl ;
      Loci::Abort() ;
    }
    Loci::debugout << "reading modules in from file '" << module <<"'" << endl ;
    return tmp ;
  }
  
  //The rules loaded in from a module are stored in a data-structure
  //called mod_info. There is a static struct mod_db which stores the
  //mod_infos' so that we load in the rules from a module only
  //once. When we pass in a string to get_info it searches the mod_db
  //to see if the rules are already loaded or else it opens the shared 
  //module and stores the rule_list in the mod_info. 
  mod::mod_info& mod::mod_db::get_info(const std::string &str) { 
    MI msi ;
    mod_info md ;
    std::string tmp_str ;
    tmp_str.append(str) ;
    tmp_str.append("_m.so") ;
    std::vector<std::string> default_ns_vec ;
    
    if((msi = mod_map.find(tmp_str)) == mod_map.end()) {
      default_ns_vec.push_back(str) ;
      if(Loci::MPI_rank == 0)
	cout << "Loading  in rules from  " << tmp_str << endl ;
      md.m_library = dlopen_helper(tmp_str) ;
      if(md.m_library == 0) {
	if(Loci::MPI_rank == 0)
	  cerr << "unable to open " << tmp_str.c_str() << endl ;
	const char *error = dlerror() ;
	if(error)
	  if(Loci::MPI_rank == 0)
	    cerr << "reason for failure is " << error << endl ;
	Loci::Abort() ;
      }
      //Copy the loaded rule_list to the mod_info rule_list. 
      md.loaded_rule_list.copy_rule_list(register_rule_list) ;
      // copy the loaded keyspace list to the mod_info keyspace list
      md.loaded_keyspace_list.copy_space_list(register_key_space_list) ;
      // We  have to clear the static register_rule_list as it is
      // going to be used to push in the next loaded rule list. 
      register_rule_list.clear() ;
      register_key_space_list.clear() ;
      md.mod_name = tmp_str ;

      put_info(md) ;
      return mod_map[tmp_str] ;
    }
    return msi->second ;
  }
  //This is similar to the above get_info but for the init_model
  //part. Some modules need an init_model. In that case we need to
  //pass in the problem name(usually the grid/vars file to be read)
  //along with the fact_database.    
  
  mod::mod_info& mod::mod_db::get_info(const std::string &str, const std::string &to_str, const char* problem_name, gfact_db &facts) { 
    MI msi ;
    mod_info md ;
    std::string tmp_str ;
    tmp_str.append(str) ;
    tmp_str.append("_m.so") ;
    std::vector<std::string> default_ns_vec ;

#ifdef VERBOSE
    debugout << "get_info("<< str << ","<<to_str<<","<< problem_name<<")"
             << endl ;
#endif
    if((msi = mod_map.find(tmp_str)) == mod_map.end()) {
      default_ns_vec.push_back(str) ;
      if(Loci::MPI_rank == 0)
	cout << "Loading  in rules from  " << tmp_str << endl ;
      md.m_library = dlopen_helper(tmp_str) ;
      if(md.m_library == 0) {
	if(Loci::MPI_rank == 0)
	  cerr << "unable to open " << tmp_str.c_str() << endl ;
	const char *error = dlerror() ;
	if(error)
	  if(Loci::MPI_rank == 0)
	    cerr << "reason for failure is " << error << endl ;
	Loci::Abort() ;
      }
      md.m_init_model = (void (*)(gfact_db &, rule_db &, const char *))
	dlsym(md.m_library,"init_model") ;
      md.loaded_rule_list.copy_rule_list(register_rule_list) ;
      md.loaded_keyspace_list.copy_space_list(register_key_space_list) ;
      register_rule_list.clear() ;
      register_key_space_list.clear() ;
      md.mod_name = tmp_str ;

      put_info(md) ;
      return mod_map[tmp_str] ; 
    }
    return msi->second ;
  } 
  
  mod::mod_db *mod::mdb = 0 ;
  vector<variableSet> unnamedVarList ;

  void load_module(const std::string from_str, const std::string to_str, rule_db& rdb, std::set<std::string> &str_set) {
#ifdef VERBOSE
    debugout << "Calling load_module with " << from_str
             << "," << to_str << endl ;
#endif
    std::vector<std::string> using_ns_vec ;
    str_set.insert(from_str) ;
    mod md(from_str) ;
    mod::mod_info m = md.get_info(from_str) ;
    variableSet nonamespace_vars  ;
    if(!unnamedVarList.empty())
      nonamespace_vars = unnamedVarList.back() ;

    for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
#ifdef VERBOSE
      debugout << "iterating over *gi = " << *gi << endl ;
#endif
      if(!(gi.get_p())->rr->is_module_rule()) {
	if(!to_str.empty()) {
	  size_t tmp = 0 ;
	  std::string sub_str = to_str ;
	  while(tmp != std::string::npos) {
	    tmp = sub_str.find("_") ;
	    using_ns_vec.push_back(sub_str.substr(0,tmp)) ;
	    sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
	  }
	  rule_implP rp = *gi;
	  variableSet vars = rp->get_var_list() ;
	  std::map<variable,variable> new_vars;
	  for(variableSet::variableSetIterator i=vars.begin();i!=vars.end();++i) {
	    if(nonamespace_vars.inSet(*i) || 
	       i->is_time_variable() || 
	       i->get_info().name == "OUTPUT" ||
	       i->get_info().name == "UNIVERSE" ||
	       i->get_info().name == "EMPTY")
	      new_vars[*i] = *i ;
	    else {
	      variable tmp_var = *i ;
	      for(int j = using_ns_vec.size()-1; j >= 0 ; --j)
		tmp_var =  tmp_var.add_namespace(using_ns_vec[j]); 
	      new_vars[*i] = tmp_var ;
	    }
	  }
	  rp->rename_vars(new_vars) ;
	  rdb.add_rule(Loci::rule(rp)) ; 
	} else {
#ifdef VERBOSE
          debugout << "adding rule " << *gi << endl ;
#endif
	  rdb.add_rule(Loci::rule(*gi)) ;
        }
      } else {
	if(Loci::MPI_rank == 0)
	  cerr << "Module rule found in " << from_str << endl ;
	std::string load  =  ((Loci::register_module*)(gi.get_p()->rr))->using_nspace() ;
	std::vector<std::string> str_vec ;
	parse_str(load, str_vec) ;
	for(size_t i = 0; i < str_vec.size(); ++i) 
	  if(str_set.find(str_vec[i]) == str_set.end()) {
	    if(Loci::MPI_rank == 0)
	      cout << "loading in rules from " << str_vec[i] <<"  for module " << to_str << endl ; 
	    load_module(str_vec[i], to_str, rdb, str_set) ;
	  }
      }
    }
    // finally add the keyspace list to the global one
    global_key_space_list.copy_space_list(m.loaded_keyspace_list) ;
  }
  
  void load_module(const std::string from_str, const std::string to_str, const char* problem_name, gfact_db &facts, rule_db& rdb, std::set<std::string> &str_set) {
#ifdef VERBOSE
    debugout << "load_module using " << from_str << "," << to_str
             << "," << problem_name << endl ;
#endif
    str_set.insert(from_str) ;
    mod md(from_str, to_str, problem_name, facts) ;
    mod::mod_info m = md.get_info(from_str, to_str, problem_name, facts) ;
#ifdef VERBOSE
    debugout << "after get_info" << endl ;
#endif
    variableSet input_vars, output_vars ;
    std::vector<std::string> using_ns_vec ;
    
    variableSet nonamespace_vars  ;
    if(!unnamedVarList.empty())
      nonamespace_vars = unnamedVarList.back() ;

    if(!to_str.empty()) {
      size_t tmp = 0 ;
      std::string sub_str = to_str ;
      while(tmp != std::string::npos) {
	tmp = sub_str.find("_") ;
	using_ns_vec.push_back(sub_str.substr(0,tmp)) ;
	sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
      }
      for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
	if((gi.get_p())->rr->is_module_rule()) {
	  if(Loci::MPI_rank == 0)
	    cerr << "Module rule found in " << from_str << endl ;
	  std::string load  =  ((Loci::register_module*)(gi.get_p()->rr))->using_nspace() ;
	  std::vector<std::string> str_vec ;
	  parse_str(load, str_vec) ;
	  for(size_t i = 0; i < str_vec.size(); ++i) 
	    if(str_set.find(str_vec[i]) == str_set.end()) {
	      if(Loci::MPI_rank == 0)
		cout << "loading in rules from " << str_vec[i] <<"  for module " << to_str << endl ; 
	      load_module(str_vec[i], to_str, problem_name, facts, rdb, str_set) ;
	    }
	  load = ((Loci::register_module*)(gi.get_p()->rr))->input_vars() ;
	  input_vars = variableSet(expression::create(load)) ;
	  load = ((Loci::register_module*)(gi.get_p()->rr))->output_vars() ;
	  output_vars = variableSet(expression::create(load)) ;
	  break ;
	}
      }
      nonamespace_vars += input_vars + output_vars ;

      for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
	if(!(gi.get_p())->rr->is_module_rule()) {
	  rule_implP rp = *gi ;
	  variableSet vars = rp->get_var_list() ;
	  if(!rp->is_specialized()) {
	    if(rp->is_parametric_provided())
	      vars += rp->get_parametric_variable() ;
	  }
	  std::map<variable,variable> new_vars;
	  for(variableSet::variableSetIterator i=vars.begin();i!=vars.end();++i)
	    if(nonamespace_vars.inSet(*i) || 
	       i->is_time_variable() || 
	       i->get_info().name == "OUTPUT" ||
	       i->get_info().name == "UNIVERSE" ||
	       i->get_info().name == "EMPTY")
	      new_vars[*i] = *i ;
	    else {
	      variable tmp_var = *i ;
	      for(int j = using_ns_vec.size()-1; j >= 0; --j) 
		tmp_var = tmp_var.add_namespace(using_ns_vec[j]); 
	      new_vars[*i] = tmp_var ;
	    }
	  rp->rename_vars(new_vars) ;
	  rdb.add_rule(Loci::rule(rp)) ;
	}
      }
    } else {
      for(rule_impl_list::iterator gi = m.loaded_rule_list.begin(); gi !=m.loaded_rule_list.end(); ++gi) {
	if(!(gi.get_p())->rr->is_module_rule()) {
	  rule_implP rp = *gi ;
	  rdb.add_rule(Loci::rule(rp)) ;
	}
      }
    }
    if(m.m_init_model != 0) {
#ifdef VERBOSE
      debugout << "init model with to_str = " << to_str << endl ;
#endif
      unnamedVarList.push_back(nonamespace_vars) ;
      int cnt = 0 ;
      if(!to_str.empty()) {
        size_t tmp = 0 ;
        std::string sub_str = to_str ;
        while(tmp != std::string::npos) {
          tmp = sub_str.find("_") ;
#ifdef VERBOSE
	  debugout << "adding namespace " << sub_str.substr(0,tmp) << endl ;
#endif
          facts.set_namespace(sub_str.substr(0,tmp)) ;
          sub_str = sub_str.substr(tmp+1, sub_str.size()) ;
	  cnt++ ;
        }
      }
      m.m_init_model(facts,rdb,problem_name) ;
      
#ifdef VERBOSE
      debugout << "removing namespace cnt = " << cnt << endl ;
#endif
      for(int i=0;i<cnt;++i)
	facts.remove_namespace() ;
      unnamedVarList.pop_back() ;
    }
    // finally add the keyspace list to the global one
    global_key_space_list.copy_space_list(m.loaded_keyspace_list) ;
  }	 
}
