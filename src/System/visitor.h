#ifndef VISITOR_H
#define VISITOR_H

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <Tools/intervalSet.h>
#include <Tools/digraph.h>
#include "sched_tools.h"
#include "visitorabs.h"
#include <algorithm>

namespace Loci {

  class orderVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    std::vector<digraph::vertexSet>
    order_dag(const digraph &g,
              digraph::vertexSet start_vertices = EMPTY,
              digraph::vertexSet only_vertices =
              interval(UNIVERSE_MIN,UNIVERSE_MAX)
              ) ;
  } ;

  class assembleVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    void compile_dag_sched(std::vector<rule_compilerP> &dag_comp,
                           const std::vector<digraph::vertexSet> &dag_sched,
                           const rulecomp_map& rcm,
                           const digraph& dag) ;
  } ;

  class graphVisualizeVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  } ;
  
  // generate allocation information table
  // used in top - down order
  class allocInfoVisitor: public visitor {
  public:
    // need graph_sn information
    allocInfoVisitor(const std::set<int>& gsn,
                     const variableSet& rtv,
                     const std::set<int>& lsn,
                     const std::map<int,variableSet>& rot_vt,
                     const std::map<int,variableSet>& lcommont,
                     const variableSet& untyped_vars)
      :graph_sn(gsn),recur_target_vars(rtv),
       loop_sn(lsn),rotate_vtable(rot_vt),
       loop_common_table(lcommont),
       allocated_vars(untyped_vars){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<int,variableSet> get_alloc_table() const
    {return alloc_table ;}
    std::map<int,variableSet> get_loop_alloc_table() const
    {return loop_alloc_table ;}
  protected:
    // return the all the variables that allocated
    // in the graph
    variableSet gather_info(const digraph& gr,int id) ;
    std::map<int,variableSet> alloc_table ;
    // holds info about the allocation of every loop super node
    std::map<int,variableSet> loop_alloc_table ;
    // variables that have been allocated up to now
    // or reserved not to be allocated
    variableSet allocated_vars ; 
    // super nodes that have a graph inside it
    std::set<int> graph_sn ;
    // the set of all the recurrence target variable
    variableSet recur_target_vars ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // used to decorate the graph to include allocation rules
  class allocGraphVisitor: public visitor {
  public:
    allocGraphVisitor(const std::map<int,variableSet>& t,
                      const std::set<int>& lsn,
                      const std::map<int,variableSet>& rot_vt,
                      const std::map<int,variableSet>& lcommont)
      :alloc_table(t),loop_sn(lsn),rotate_vtable(rot_vt),
       loop_common_table(lcommont){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  protected:
    void edit_gr(digraph& gr,rulecomp_map& rcm,int id) ;
    std::map<int,variableSet> alloc_table ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // generate delete information table
  // used in top - down order
  class deleteInfoVisitor: public visitor {
  public:
    // need loop allocation info and recurrence variable info
    deleteInfoVisitor(const std::map<int,variableSet>& lat,
                      const std::map<variable,variableSet>& rvt2s,
                      const std::map<variable,variableSet>& rvs2t,
                      const std::set<int>& gsn,
                      const std::map<int,int>& pnt,
                      const std::map<int,int>& lct,
                      const std::set<int>& lsn,
                      const std::map<int,variableSet>& rot_vt,
                      const std::map<int,variableSet>& lcommont,
                      const variableSet& reserved_vars) ;
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<int,variableSet> get_delete_table() const
    {return delete_table ;}
    std::map<variable,ruleSet> get_recur_source_other_rules() const
    {return recur_source_other_rules ;}
  protected:
    // determine whether variables in working_vars
    // should be deleted in this graph

    // return all the variables that deleted
    // in the graph
    variableSet gather_info(const digraph& gr,
                            const variableSet& working_vars) ;

    // gather information for gather_info function use
    // get working_vars in the graph
    variableSet get_start_info(const digraph& gr,int id) ;
    // constrain deletions of graphs that are inside n number of loops
    int only_loop_alloc(variableSet& working_vars,int id,int loop_num) ;
    // the looping algorithm for schedule deletion of variables
    // in the multilevel graph
    void looping_algr(const variableSet& working_vars,const digraph& dg,
                      int id,int start_loop_num) ;
    // given a rule set and a graph
    // return true if there is only one super node and
    // all other rules have path to the super node, or if
    // it is only one super node
    // return false otherwise
    bool let_it_go(const digraph& gr, ruleSet rules, const variable& v) ;
    // the delete information table
    std::map<int,variableSet> delete_table ;
    // variables that have been deleted up to now or reserved
    // not to be deleted
    variableSet deleted_vars ;
    // info about recurrence variables (target to source map)
    std::map<variable,variableSet> recur_vars_t2s ;
    // source to target map
    std::map<variable,variableSet> recur_vars_s2t ;
    // if a variable is a recurrence target variable,
    // and there are other rules (other than the recurrence rule)
    // that its recurrence source variables reach in the same graph,
    // we then record such rules in this table
    std::map<variable,ruleSet> recur_source_other_rules ;
    // info about every loop's allocation
    std::map<int,variableSet> loop_alloc_table ;
    // super nodes that have a graph inside it
    std::set<int> graph_sn ;
    // table that records the parent node 
    std::map<int,int> pnode_table ;
    // table that holds the collapse node of each loop
    std::map<int,int> loop_ctable ;
    // set that holds all the loop node id
    std::set<int> loop_sn ;
    // table that holds rotate list variables in each loop
    std::map<int,variableSet> rotate_vtable ;
    // table that holds the common variables
    // between adv & col part of each loop
    std::map<int,variableSet> loop_common_table ;
    // all the collapse node id
    std::set<int> col_sn ;
    // all the loop rotate variables
    variableSet all_rot_vars ;
  } ;

  // visitor that get all the recurrence variables in the
  // multilevel graph
  class recurInfoVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<variable,variableSet> get_recur_vars_t2s() const
    {return recur_vars_t2s ;}
    std::map<variable,variableSet> get_recur_vars_s2t() const
    {return recur_vars_s2t ;}
    variableSet get_recur_source_vars() const
    {return recur_source_vars ;}
    variableSet get_recur_target_vars() const
    {return recur_target_vars ;}
    
    std::map<variable,variableSet> get_generalize_t2s() const
    {return generalize_t2s ;}
    std::map<variable,variableSet> get_generalize_s2t() const
    {return generalize_s2t ;}
    variableSet get_generalize_source_vars() const
    {return generalize_source_vars ;}
    variableSet get_generalize_target_vars() const
    {return generalize_target_vars ;}
    
    std::map<variable,variableSet> get_promote_t2s() const
    {return promote_t2s ;}
    std::map<variable,variableSet> get_promote_s2t() const
    {return promote_s2t ;}
    variableSet get_promote_source_vars() const
    {return promote_source_vars ;}
    variableSet get_promote_target_vars() const
    {return promote_target_vars ;}
    
    std::map<variable,variableSet> get_priority_t2s() const
    {return priority_t2s ;}
    std::map<variable,variableSet> get_priority_s2t() const
    {return priority_s2t ;}
    variableSet get_priority_source_vars() const
    {return priority_source_vars ;}
    variableSet get_priority_target_vars() const
    {return priority_target_vars ;}
    
    std::map<variable,variableSet> get_rename_t2s() const
    {return rename_t2s ;}
    std::map<variable,variableSet> get_rename_s2t() const
    {return rename_s2t ;}
    variableSet get_rename_source_vars() const
    {return rename_source_vars ;}
    variableSet get_rename_target_vars() const
    {return rename_target_vars ;}
  protected:
    void gather_info(const digraph& gr) ;
    // from x{n} -> x, i.e. from target -> source
    std::map<variable,variableSet> recur_vars_t2s ; 
    // from source -> target, e.g. x -> x{n}
    std::map<variable,variableSet> recur_vars_s2t ;

    std::map<variable,variableSet> generalize_t2s ;
    std::map<variable,variableSet> generalize_s2t ;

    std::map<variable,variableSet> promote_t2s ;
    std::map<variable,variableSet> promote_s2t ;

    std::map<variable,variableSet> priority_t2s ;
    std::map<variable,variableSet> priority_s2t ;

    std::map<variable,variableSet> rename_t2s ;
    std::map<variable,variableSet> rename_s2t ;
    
    // set of recurrence source and target variables
    variableSet recur_source_vars ;
    variableSet recur_target_vars ;

    variableSet generalize_source_vars ;
    variableSet generalize_target_vars ;

    variableSet promote_source_vars ;
    variableSet promote_target_vars ;

    variableSet priority_source_vars ;
    variableSet priority_target_vars ;

    variableSet rename_source_vars ;
    variableSet rename_target_vars ;
  } ;

  // used to decorate the graph to include deletion rules
  class deleteGraphVisitor: public visitor {
  public:
    deleteGraphVisitor(const std::map<int,variableSet>& t,
                       const std::map<variable,ruleSet>& rsor)
      :delete_table(t),recur_source_other_rules(rsor){}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  protected:
    void edit_gr(digraph& gr,rulecomp_map& rcm,int id) ;
    std::map<int,variableSet> delete_table ;
    // if a variable is a recurrence target variable,
    // and there are other rules (other than the recurrence rule)
    // that its recurrence source variables reach in the same graph,
    // this table holds such rules
    std::map<variable,ruleSet> recur_source_other_rules ;
  } ;

  // visitor to get some inter- super node information
  class snInfoVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::set<int> get_graph_sn() const
    {return graph_sn ;}
    std::set<int> get_loop_sn() const
    {return loop_sn ;}
    std::map<int,std::set<int> > get_subnode_table() const
    {return subnode_table ;}
    std::map<int,int> get_loop_col_table() const
    {return loop_col_table ;}
  protected:
    void fill_subnode_table(const digraph& gr, int id) ;
    // super node that has a graph in it
    // i.e. the loop, dag, and conditional super node
    std::set<int> graph_sn ;
    // all the loop super node id
    std::set<int> loop_sn ;
    // table that holds all the super nodes' id inside a super node
    std::map<int,std::set<int> > subnode_table ;
    // table that holds the collapse node of each loop
    std::map<int,int> loop_col_table ;
  } ;

  // check the if the allocation and deletion table are
  // consistent
  class allocDeleteStat {
    variableSet allocated_vars ;
    variableSet deleted_vars ;
    std::map<variable,variableSet> recur_vars_t2s ;
    std::map<variable,variableSet> recur_vars_s2t ;
    variableSet recur_source_vars ;
    variableSet recur_target_vars ;
  public:
    allocDeleteStat(const std::map<int,variableSet>& alloc_table,
                    const std::map<int,variableSet>& delete_table,
                    const std::map<variable,variableSet>& t2s,
                    const std::map<variable,variableSet>& s2t,
                    const variableSet& rsv,const variableSet& rtv) ;
    std::ostream& report(std::ostream& s) const ;
  } ;

  // function to get the multilevel graph parent hierarchy table
  std::map<int,int>
  get_parentnode_table(const std::map<int,std::set<int> >& subnode_table) ;
  // function to get the allocation of each loop
  std::map<int,variableSet>
  get_loop_alloc_table(const std::map<int,variableSet>& alloc_table,
                       const std::map<int,std::set<int> >& subnode_table,
                       const std::set<int>& loop_sn,
                       const std::set<int>& graph_sn,
                       const std::map<variable,variableSet>& rvs2t) ;

  // visitor to compute the loop_rotate lists
  class rotateListVisitor: public visitor {
  public:
    rotateListVisitor(const sched_db& sd):scheds(sd) {}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) {}
    virtual void visit(conditional_compiler& cc) {}
    std::map<int,variableSet> get_rotate_vars_table() const
    {return rotate_vars_table ;}
    std::map<int,variableSet> get_loop_common_table() const
    {return loop_common_table ;}
  private:
    // reference to the schedule database
    const sched_db& scheds ;
    // table that holds variables in each loop's
    // rotate list
    std::map<int,variableSet> rotate_vars_table ;
    // table holds the common varibles between adv & col part of loop
    std::map<int,variableSet> loop_common_table ;
  } ;

  // visitor that checks if a graph has cycle
  class dagCheckVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
  private:
    bool check_dag(digraph gr) ;
  } ;

  // visitor that discover all the un-typed variables
  // in the multilevel graph
  class unTypedVarVisitor: public visitor {
  public:
    unTypedVarVisitor(const std::map<variable,variableSet>& s2t,
                      const std::map<variable,variableSet>& t2s,
                      const variableSet& input) ;

    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    variableSet get_untyped_vars() const
    {return untyped_vars ;}
  private:
    void discover(const digraph& gr) ;
    variableSet untyped_vars ;
    variableSet typed_vars ;
    // recurrence variable mapping table
    // source to target
    std::map<variable,variableSet> recur_vars_s2t ;
    // target to source
    std::map<variable,variableSet> recur_vars_t2s ;
  } ;

  // visitor that collects all the variables
  // in the multilevel graph
  class getAllVarVisitor: public visitor {
  public:
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    variableSet get_all_vars_ingraph() const
    {return all_vars ;}
  private:
    void collect_vars(const digraph& gr) ;
    variableSet all_vars ;
  } ;

  // visitor that does preprocess for promoted variables
  // determine which promoted variable should be deleted.
  // e.g. there are x -> x{n} -> x{n,it}, we can only delete
  // one of them. For generalize,priority and renamed variables,
  // we should always delete and process the final one.
  class promotePPVisitor: public visitor {
  public:
    promotePPVisitor(const std::map<variable,variableSet>& pt2s,
                     const std::map<variable,variableSet>& ps2t,
                     const variableSet& psource,
                     const variableSet& ptarget,
                     const std::set<int>& gsn,
                     variableSet& input) ;
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    variableSet get_rep() const
    {return rep ;}
    variableSet get_remaining() const
    {return remaining ;}
  private:
    // set that holds the representitive promote variable
    variableSet rep ;
    // set that holds the remaining promote variables
    // that we want to reserve
    variableSet remaining ;
    variableSet processed ;
    bool is_rep(const digraph& gr, ruleSet rules) ;
    void pick_rep(const digraph& gr) ;
    variableSet promote_source_vars, promote_target_vars ;
    std::map<variable,variableSet> promote_t2s, promote_s2t ;
    std::set<int> graph_sn ;
    variableSet reserved_vars ;
  } ;

  // function to analyze the renamed target variable cluster
  // and find out a representitive variable for each cluster,
  // only this variable is allowed to be deleted
  variableSet pick_rename_target(const std::map<variable,variableSet>& s2t,
                                 const std::map<variable,variableSet>& t2s,
                                 const variableSet& allvars) ;

  // function that checks some preconditions that all the recurrence
  // variables should meet
  void check_recur_precondition(const recurInfoVisitor& v,
                                const variableSet& input) ;

  // visitor that finds rule chains that are suitable for chopping
  typedef std::pair<digraph,variableSet> chop_chain ;
  class chopRuleVisitor: public visitor {
  public:
    chopRuleVisitor(fact_db& fd):facts(fd) {}
    virtual void visit(loop_compiler& lc) ;
    virtual void visit(dag_compiler& dc) ;
    virtual void visit(conditional_compiler& cc) ;
    std::map<int,std::list<chop_chain> > get_all_chains() const
    {return all_chains ;}
    std::ostream& visualize(std::ostream& s) const ;
  private:
    void find_chain(const digraph& gr, int id) ;
    std::map<int,std::list<chop_chain> > all_chains ;
    // reference to the fact database
    fact_db& facts ;
  } ;
  
  // overload "<<" to print out an std::set
  template<typename T>
  inline std::ostream& operator<<(std::ostream& s, const std::set<T>& ss) {
    typename std::set<T>::const_iterator si ;
    si=ss.begin() ;
    s << *si ;
    ++si ;
    for(;si!=ss.end();++si) {
      s << " " ;
      s << *si ;
    }

    return s ;
  }
  
  // overload "<<" to print out an std::map
  template<typename T1, typename T2>
  inline std::ostream& operator<<(std::ostream& s, const std::map<T1,T2>& m) {
    typename std::map<T1,T2>::const_iterator mi ;
    for(mi=m.begin();mi!=m.end();++mi) {
      s << mi->first
        << ": " << mi->second
        << std::endl ;
    }

    return s ;
  }

  // pretty printing of a rule's signature
  // i.e. remove namespace info, if any
  inline std::string pretty_sig(const rule& r) {
    //the following line is the simplest, but it does not
    //include the qualify, such like "SN1:" in a super node
    //return r.get_info().desc.rule_identifier() ;
    std::string name = r.get_info().name() ;
    if(r.type() == rule::INTERNAL)
      return name ;
    else {
      std::string::iterator pos ;
      pos = std::find(name.begin(),name.end(),'#') ;
      return std::string( (pos==name.end()?name.begin():pos+1),name.end()) ;
    }
  }

  // function to get all the recur targets of a variableSet
  // from a given recur mapping table
  variableSet
  get_recur_target_for_vars(const variableSet& vars,
                            const std::map<variable,variableSet>& t) ;

}

#endif
