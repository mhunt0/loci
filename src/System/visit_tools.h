#ifndef VISIT_TOOLS_H
#define VISIT_TOOLS_H

#include "visitor.h"
#include "comp_tools.h"

#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <list>

namespace Loci {

  // some inlined small functions
  
  inline rule create_rule(variable sv, variable tv, std::string qualifier) {
    std::ostringstream oss ;
    oss << "source(" << sv << ')' ;
    oss << ",target(" << tv << ')' ;
    oss << ",qualifier(" << qualifier << ')' ;
    std::string sig = oss.str() ;
    rule r(sig) ;
    return r ;
  }
  
  inline rule create_rule(variableSet source, variableSet target,
                          std::string qualifier) {
    std::ostringstream oss ;
    oss << "source(" << source << ')' ;
    oss << ",target(" << target << ')' ;
    oss << ",qualifier(" << qualifier << ")" ;
    std::string sig = oss.str() ;
    rule r(sig) ;
    return r ;
  }

  inline bool is_super_node(int rid)
  {
    if(rid >= 0) // a variable
      return false ;
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;
    
    return (rqualifier.substr(0,2) == "SN") ;
  }
  
  inline bool is_super_node(const ruleSet::const_iterator& ruleIter) {
    std::string rqualifier = ruleIter->get_info().qualifier() ;

    return (rqualifier.substr(0,2) == "SN") ;
  }

  inline bool is_super_node(const rule& r) {
    std::string rqualifier = r.get_info().qualifier() ;

    return (rqualifier.substr(0,2) == "SN") ;
  }

  inline int get_supernode_num(const rule& r) {
    std::string rqualifier = r.get_info().qualifier() ;
    std::string head = rqualifier.substr(0,2) ;
    if(head != "SN") {
      std::cerr << "get_supernode_num error! pass in rule is not a super node!"
           << std::endl ;
      exit(-1) ;
    }
      
    std::string number = rqualifier.substr(2,rqualifier.size()-2) ;
    std::stringstream ss ;
    ss << number ;
    int ret ;
    ss >> ret ;
      
    return ret ;
  }

  template<typename T>
  inline bool inSet(const std::set<T>& s, const T& elem) {
    typename std::set<T>::const_iterator si ;
    si = s.find(elem) ;
    if(si == s.end())
      return false ;
    return true ;
  }

  // is there a path from source vertex to target vertex in gr?
  inline bool has_path(const digraph& gr, int source, int target) {
    bool path = false ;
    digraph::vertexSet working = gr[source] ;
    while(working != EMPTY) {
      if(working.inSet(target)) {
        path = true ;
        break ;
      }
      digraph::vertexSet children ;
      for(digraph::vertexSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        children += gr[*vi] ;
      }
      working = children ;
    }

    return path ;
  }

  // given a variableSet, convert them to vertexSet
  inline digraph::vertexSet get_vertexSet(const variableSet& vars) {
    digraph::vertexSet ret ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
      ret += vi->ident() ;

    return ret ;
  }

  // given a ruleSet, convert them to vertexSet
  inline digraph::vertexSet get_vertexSet(const ruleSet& rules) {
    digraph::vertexSet ret ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      ret += ri->ident() ;

    return ret ;
  }

  // get the init source variableSet from a target variable and
  // the recurrence target to source table, or
  // get the final target variableSet from a source variable and
  // the recurrence source to target table
  inline
  variableSet get_leaf_recur_vars(const std::map<variable,variableSet>& t,
                                  const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret ;
    std::list<variable> working ;

    // first check if v is in the table
    found = t.find(v) ;
    // if v is not recurrence variable, return empty set
    if(found == t.end())
      return variableSet(EMPTY) ;

    // do a depth first search
    working.push_front(v) ;
    while(!working.empty()) {
      variable cur = working.front() ;
      working.pop_front() ;
        
      variableSet tmp ;
      found = t.find(cur) ;
        
      if(found != t.end())
        tmp = found->second ;
      else
        ret += cur ;
        
      for(variableSet::const_iterator vi=tmp.begin();
          vi!=tmp.end();++vi)
        working.push_front(*vi) ;
    }

    return ret ;
  }

  // get all the source variables from a target variable and
  // the recurrence target to source table, or
  // get all the target variables from a source variable and
  // the recurrence source to target table
  inline
  variableSet get_all_recur_vars(const std::map<variable,variableSet>& t,
                                 const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret, working ;
      
    // first check if v is in the table
    found = t.find(v) ;
    // if v is not recurrence variable, return empty set
    if(found == t.end())
      return variableSet(EMPTY) ;

    // do a breadth first search, and add up all the results
    working += v ;
    while(working != EMPTY) {
      variableSet tmp ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        found = t.find(*vi) ;
        if(found != t.end())
          tmp += found->second ;
      }
      ret += tmp ;
      working = tmp ;
    }

    return ret ;
  }

  // given a variable, if it is a recurrence variable, then get
  // all other leaf target variables that all its source can reach
  inline
  variableSet get_all_leaf_target(const std::map<variable,variableSet>& t2s,
                                  const std::map<variable,variableSet>& s2t,
                                  const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret, working ;

    working += v ;
    while(working != EMPTY) {
      variableSet tmp ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        found = t2s.find(*vi) ;
        if(found != t2s.end()) {
          variableSet cur = found->second ;
          tmp += cur ;
          for(variableSet::const_iterator vi2=cur.begin();
              vi2!=cur.end();++vi2)
            ret += get_leaf_recur_vars(s2t,*vi2) ;
        }
      }
      working = tmp ;
    }

    ret -= v ;
    return ret ;
  }

  inline bool is_dg_empty(const digraph& gr) {
    digraph::vertexSet allv = gr.get_all_vertices() ;
    return (allv == EMPTY) ;
  }

  // return if a rule is a generalize, promote or priority rule
  inline bool is_recur_rule(const ruleSet::const_iterator& ruleIter) {
    if(ruleIter->type() == rule::INTERNAL) {
      if( (ruleIter->get_info().qualifier() == "promote") ||
          (ruleIter->get_info().qualifier() == "generalize") ||
          (ruleIter->get_info().qualifier() == "priority")
          )
        return true ;
    }
    return false ;
  }
  
  inline bool is_internal_rule(const ruleSet::const_iterator& ruleIter) {
    return (ruleIter->type() == rule::INTERNAL) ;
  }

  inline bool is_internal_rule(const rule& r) {
    return (r.type() == rule::INTERNAL) ;
  }
  
  inline bool is_internal_rule(int i) {
    if(i>=0) // is a variable
      return false ;
    rule r(i) ;
    return (r.type() == rule::INTERNAL) ;
  }

  inline bool is_virtual_rule(const rule& r) {
    if(r.type() == rule::INTERNAL) {
      return (
              is_super_node(r) ||
              (r.get_info().qualifier() == "generalize") ||
              (r.get_info().qualifier() == "promote") ||
              (r.get_info().qualifier() == "priority")
              ) ;
    }else
      return false ;
  }

  inline bool time_before(const variable& v1, const variable& v2) {
    return v1.time().before(v2.time()) ;
  }
  
  inline bool time_equal(const variable& v1, const variable& v2) {
    return (!time_before(v1,v2) && !time_before(v2,v1)) ;
  }
  
  inline bool time_after(variable v1, variable v2) {
    return (!time_before(v1,v2) && !time_equal(v1,v2)) ;
  }

} // end of namespace Loci

#endif
