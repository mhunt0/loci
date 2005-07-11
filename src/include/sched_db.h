#ifndef SCHED_DB_H
#define SCHED_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <fact_db.h>
#include <rule.h>

namespace Loci {
  class sched_db {
    
    struct sched_data {
      variableSet aliases,rotations,antialiases ;
      bool ismap ;
      MapRepP minfo ;
      std::map<entitySet,entitySet> imageMap ;
      std::map<entitySet,std::pair<entitySet,entitySet> > preimageMap
      ;
      sched_data() {} 
      sched_data(variable v, storeRepP &st)
	{ aliases += v ; 
	ismap = (st->RepType() == Loci::MAP);
	if(ismap) minfo = MapRepP(st->getRep()) ; }
    } ;
    
    struct existential_info {
      variable v ;
      entitySet exists ;
      existential_info() {} 
      existential_info(const variable &vin,const entitySet &e) :
        v(vin), exists(e) {}
    } ;
    struct sched_info ;
    friend struct sched_info ;
    struct sched_info {
      int sched_info_ref ;
      entitySet fact_installed ;
      variableSet synonyms ;
      std::map<rule,existential_info> exist_map ;
      entitySet existence ;
      entitySet requested ;
      entitySet shadow ; // Used by distributed memory apply rules
      
      //Apply rules may have target variables both on input and output side.
      //Sometimes apply rules may need to access entities of target variable
      //which are not requested.  In that case, UNIT rule also allocates these 
      //entities in addition to entities being requested.
      entitySet extra_unit_request;

      //////////////////////////Duplication Related:////////////// 
      //Defines maximum which target variable entities a rule can compute 
      std::map<rule, entitySet> proc_able_map;

      //Defines which target variable entities a rule can compute using context
      //that is subset of my_entities
      std::map<rule, entitySet> my_proc_able_map;

      unsigned int policy;  //Each bit defines a policy 
      bool duplicate_variable; //Defines if a variable is duplicated
      
      //Applies only to reduce variable: considering all the rules of a variable, 
      //it defines which entities can definitely computed on a processor
      entitySet reduce_proc_able_entities; 

      bool reduction_outputmap; //If any of the apply or unit rule has mapping in output

      double original_computation_time, duplication_computation_time;
      double original_communication_time, duplication_communication_time;
      //////////////////////////////////////////////////////////////

      sched_info(int ref = -1) {
	sched_info_ref = ref ;
	policy = 0;
	duplicate_variable = false;
	reduce_proc_able_entities = ~EMPTY;
	reduction_outputmap = false;

	original_computation_time = 0;
	duplication_computation_time = 0;
	original_communication_time = 0;
	duplication_communication_time = 0;
      }
    } ;

    struct model {
      double ts, tw;
      static const double INVALID_TS;
      model(double t0, double tc) {
	ts = t0;
	tw = tc;
      }
      model() { ts = INVALID_TS; tw = 0;}
      void get_parameters(double &t0, double &tc) { t0 = ts; tc = tw; }
      void set_parameters(double t0, double tc) { ts = t0; tw = tc; }
    };

    variableSet possible_duplicate_vars;
    std::map<rule, model> comp_model;
    model comm_model;
    void register_variable(variable v) ;
  
    variableSet all_vars ;
    std::map<variable,variable> synonyms ;
    typedef std::map<variable, sched_info> vmap_type ;
    vmap_type vmap ;
    std::vector<sched_data> sched_infov ;
    intervalSet free_set ;

    bool detected_errors ;
    
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      if((mi=synonyms.find(v)) != synonyms.end()) 
	return mi->second ;
      return v ;
    }
  public:
    enum  duplicate_policy{NEVER, ALWAYS, MODEL_BASED};
    sched_db() ;
    ~sched_db() ;
    sched_db(fact_db &facts) ;

    bool errors_found() {return detected_errors ;}
    void clear_errors() {detected_errors = false ;}
    void set_error() { detected_errors = true ; }
    
    void install_sched_data(variable v, sched_data data) ;
    void install_sched_info(variable v, sched_info info) ;
    
    variableSet get_synonyms(variable v) const 
      { return get_sched_info(v).synonyms ; }

    variableSet get_aliases(variable v) const
    { return get_sched_data(v).aliases ; }

      variableSet get_antialiases(variable v) const
    { return get_sched_data(v).antialiases ; }

    variableSet get_rotations(variable v) const
    { return get_sched_data(v).rotations ; }
    
    void set_variable_type(variable v, storeRepP st, fact_db &facts) ;
    void set_variable_type(std::string vname,const storeRepP st, fact_db &facts)
      { set_variable_type(variable(vname),st, facts) ; }
    
    void set_variable_type(variable v, store_instance &si, fact_db &facts)
      { set_variable_type(v,si.Rep(), facts) ; }
    void set_variable_type(std::string vname, store_instance &si, fact_db &facts)
      { set_variable_type(variable(vname),si, facts) ; }
    
    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname,const storeRepP st)
      { set_variable_type(variable(vname),st) ; }
    
    void set_variable_type(variable v, store_instance &si)
      { set_variable_type(v,si.Rep()) ; }
    void set_variable_type(std::string vname, store_instance &si)
      { set_variable_type(variable(vname),si) ; } 
    
    void variable_is_fact_at(variable v,entitySet s, fact_db &facts) ;
    void variable_is_fact_at(std::string vname, const entitySet s, fact_db &facts)
      { variable_is_fact_at(variable(vname),s, facts) ; }
    
    void variable_is_fact_at(variable v,entitySet s) ;
    void variable_is_fact_at(std::string vname, const entitySet s)
      { variable_is_fact_at(variable(vname),s) ; }

    void set_variable_rotations(variableSet rot) ;
     
    void alias_variable(variable v, variable alias, fact_db &facts) ;
    void alias_variable(std::string vname, std::string alias, fact_db &facts)
      { alias_variable(variable(vname),variable(alias), facts) ; }
   
    void alias_variable(variable v, variable alias) ;
    void alias_variable(std::string vname, std::string alias)
      { alias_variable(variable(vname),variable(alias)) ; }
    
    void synonym_variable(variable v, variable synonym, fact_db &facts) ;
    void synonym_variable(std::string vname, std::string synonym, fact_db &facts)
      { synonym_variable(variable(vname),variable(synonym), facts) ; }
    
    void synonym_variable(variable v, variable synonym) ;
    void synonym_variable(std::string vname, std::string synonym)
      { synonym_variable(variable(vname),variable(synonym)) ; }
    
    void set_existential_info(variable v,rule f,entitySet x) ;
    ruleSet get_existential_rules(variable v) ;
    const sched_info & get_sched_info(variable v) const {
      vmap_type::const_iterator mi = vmap.find(remove_synonym(v)) ;
      if(mi == vmap.end()) {
        cerr << "const get_sched_info: variable " << v << " does not exist" << endl ;
        abort() ;
      }
      return mi->second ;
    }
    const sched_data & get_sched_data(variable v) const 
  { return sched_infov[get_sched_info(v).sched_info_ref] ; }
    
    sched_info & get_sched_info(variable v) ;
    sched_data & get_sched_data(variable v) 
      { return sched_infov[get_sched_info(v).sched_info_ref] ; }

    bool is_a_Map(variable v) {
      return get_sched_data(v).ismap ;
    }
    entitySet get_existential_info(variable v, rule f) {
      sched_info &finfo = get_sched_info(v) ;
      std::map<rule,existential_info>::const_iterator mi ;
      mi = finfo.exist_map.find(f) ;
      if(mi!=finfo.exist_map.end()) {
        return mi->second.exists ;
      } else 
        return EMPTY ;
    }

    entitySet variable_existence(variable v) ;
    void variable_request(variable v, entitySet e) ;
    void set_variable_shadow(variable v, entitySet e)
    { get_sched_info(v).shadow = e ; }
    void variable_shadow(variable v, entitySet e)
    { get_sched_info(v).shadow += e ; }
    entitySet get_variable_shadow(variable v) const
    { return get_sched_info(v).shadow ; }
    entitySet get_variable_request(rule f, variable v) ;
    entitySet get_variable_requests(variable v) const 
    { return  get_sched_info(v).requested ; }
    variableSet get_typed_variables() const { return all_vars ; }
    entitySet image(variable v, entitySet e) ;
    std::pair<entitySet,entitySet> preimage(variable v, entitySet e) ;

    void set_policy(variable v, unsigned int p) { get_sched_info(v).policy = p;}
    unsigned int get_policy(variable v) { return get_sched_info(v).policy; }
    
    void add_policy(variable v, duplicate_policy p);
    
    bool is_policy(variable v, duplicate_policy p);
    
    bool is_duplicate_variable(variable v) { return get_sched_info(v).duplicate_variable; }
    
    void set_duplicate_variable(variable v, bool p) {
      variableSet synonyms = get_synonyms(v);
      for(variableSet::const_iterator vi = synonyms.begin(); vi != synonyms.end(); vi++)
	get_sched_info(*vi).duplicate_variable = p;
    } 

    entitySet get_proc_able_entities(variable v, rule f) {
      sched_info &finfo = get_sched_info(v);
      std::map<rule, entitySet>::const_iterator mi;
      mi = finfo.proc_able_map.find(f);
      if(mi != finfo.proc_able_map.end())
	return mi->second;
      else
	return EMPTY;
    }

    void set_proc_able_entities(variable v, rule f, entitySet x) {
      sched_info &finfo = get_sched_info(v);
      finfo.proc_able_map[f] += x;
    }

    entitySet get_my_proc_able_entities(variable v, rule f) {
      sched_info &finfo = get_sched_info(v);
      std::map<rule, entitySet>::const_iterator mi;
      mi = finfo.my_proc_able_map.find(f);
      if(mi != finfo.my_proc_able_map.end())
	return mi->second;
      else
	return EMPTY;
    }

    void set_my_proc_able_entities(variable v, rule f, entitySet x) {
      sched_info &finfo = get_sched_info(v);
      finfo.my_proc_able_map[f] += x;
    }

    entitySet get_reduce_proc_able_entities(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.reduce_proc_able_entities);
    }

    void set_reduce_proc_able_entities(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v);
      finfo.reduce_proc_able_entities = x;
    }

    bool is_reduction_outputmap(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.reduction_outputmap);
    }

    void set_reduction_outputmap(variable v, bool x) {
      sched_info &finfo = get_sched_info(v);
      finfo.reduction_outputmap = x;
    }

    void add_extra_unit_request(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v);
      finfo.extra_unit_request += x;
    }

    entitySet get_extra_unit_request(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.extra_unit_request);
    }

    //For map variables, it is necessary to change existence without
    //supllying any rules because Loci does not have rules
    //that create maps
    void set_variable_existence(variable v, entitySet x) {
      sched_info &finfo = get_sched_info(v);
      finfo.existence += x;
    }
    void add_model_info(double comm_ts, double comm_tw, const std::map<rule, std::pair<double, double> > &comp_info);

    void add_original_computation_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v);
      finfo.original_computation_time += add;
    }
    
    void add_duplication_computation_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v);
      finfo.duplication_computation_time += add;
    }
    
    void add_original_communication_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v);
      finfo.original_communication_time += add;
    }

    void add_duplication_communication_time(variable v, double add) {
      sched_info &finfo = get_sched_info(v);
      finfo.duplication_communication_time += add;
    }

    double get_precalculated_original_computation_time(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.original_computation_time);
    }

    double get_precalculated_duplication_computation_time(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.duplication_computation_time);
    }

    double get_precalculated_original_communication_time(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.original_communication_time);
    }

    double get_precalculated_duplication_communication_time(variable v) {
      sched_info &finfo = get_sched_info(v);
      return(finfo.duplication_communication_time);
    }

    model get_comp_model(rule r) const {
      std::map<rule, model>::const_iterator mi = comp_model.find(r);
      if(mi != comp_model.end())
	return mi->second;
      else {
	return model(model::INVALID_TS, 0);
      }
    }
    
    model get_comm_model() const { return comm_model; }

    variableSet get_possible_duplicate_vars() const { return possible_duplicate_vars; }
    void add_possible_duplicate_vars(variableSet vars) {
      for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); vi++) {
	variableSet synonyms = get_synonyms(*vi);
	possible_duplicate_vars += synonyms;
      }
    }
    
    std::ostream &print_summary(fact_db &facts, std::ostream &s) ;
  } ;
}

#endif
