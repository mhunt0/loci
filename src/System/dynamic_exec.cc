#include "sched_tools.h"
using std::cout ;
using std::cerr ;
using std::endl ;
#include "distribute.h"

namespace Loci {

void Allocate_func(fact_db &local_facts,entitySet &exec_set,variableSet inputs,variableSet outputs);

void Deallocate_func(fact_db &local_facts,variableSet inputs,variableSet outputs);

void SendInput (int tStart, int tSize, int dest, int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,fact_db &facts) ;

void ReceiveInput (int rcvStart,int &rcvSize,int src,int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,fact_db &local_facts);

void SendOutput (int tStart, int tSize, int dest,int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &outputs,fact_db &local_facts) ;

void ReceiveOutput (int rcvStart, int rcvSize, int src, int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,variableSet &outputs,fact_db &facts);

void ExecuteLoop (rule_implP rp1,entitySet &e,int method,int *yMap,fact_db &facts,fact_db &local_facts1,fact_db &local_facts2,rule_implP local_compute1,rule_implP local_compute2,variableSet &inputs,variableSet &outputs,double *stats) ;

}

namespace Loci {

dynamic_schedule_rule::dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds)  {
  
  rp = fi.get_rule_implP() ; //get rule from rule database 
  rule_tag = fi ;            //store rule tag in rule_tag
  local_compute1 = rp->new_rule_impl() ; //another instance of rule 
  local_compute2 = rp->new_rule_impl() ; //another instance of rule 
  entitySet in = rule_tag.sources() ; //inputs from rule rhs
  outputs = rule_tag.targets() ;      //outputs as in rhs  
  exec_set = eset ;    //for second half
  

  

  // Setup local facts input variables (types only no allocation)
  for(variableSet::const_iterator vi=in.begin();vi!=in.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    if((store_ptr != 0) && store_ptr->RepType() == Loci::STORE) {	  
      inputs += *vi ; 
      local_facts[0].create_fact(*vi,store_ptr->new_store(EMPTY)) ; 
      local_facts[1].create_fact(*vi,store_ptr->new_store(EMPTY)) ;             
    } else {
      local_facts[0].create_fact(*vi,facts.get_variable(*vi)) ;
      local_facts[1].create_fact(*vi,facts.get_variable(*vi)) ;
    }
  }
  // Setup local facts output variables 
  for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
    storeRepP store_ptr = rp->get_store(*vi) ;    
    local_facts[0].create_fact(*vi,store_ptr->new_store(EMPTY)) ;
    local_facts[1].create_fact(*vi,store_ptr->new_store(EMPTY)) ;
  }
  
  // Initialize both functions for remote and local execution.
  local_compute1->initialize(local_facts[0]) ;
  local_compute2->initialize(local_facts[1]) ;
  rp->initialize(facts) ;   

}


void dynamic_schedule_rule::execute(fact_db &facts) { 

  int *yMap=new int[2*Loci::MPI_processes]; 
  for(int i = 0; i < Loci::MPI_processes; i++) {
    yMap[2*i] = 0;
    yMap[2*i+1] = 0;
  }
   //Broadcast start,size to all procs
  int sendstart=0;
  sendstart = exec_set.Min();
  int *rstart = new int[Loci::MPI_processes];
  for(int q = 0; q < Loci::MPI_processes; q++) {
    rstart[q] = 0;
  }
  MPI_Allgather(&sendstart, 1, MPI_INT, rstart, 1, MPI_INT, MPI_COMM_WORLD);
  int sendsize=0;
  sendsize = exec_set.Max()-sendstart+1;
  int *rsize = new int[Loci::MPI_processes];
  for(int n = 0; n < Loci::MPI_processes; n++) {
    rsize[n] = 0;
  }
  MPI_Allgather(&sendsize, 1, MPI_INT, rsize, 1, MPI_INT, MPI_COMM_WORLD);
 
  for(int i = 0; i < Loci::MPI_processes; i++) {
    yMap[2*i] = rstart[i];
    yMap[2*i+1] = rsize[i];
  }
  delete [] rsize;
  delete [] rstart;
 
  double *stats=new double[4*MPI_processes];
  for(int i = 0; i < Loci::MPI_processes; i++) {
   stats[4*i+0] = 0;
   stats[4*i+1] =0;
   stats[4*i+2] = 0;
   stats[4*i+3] =0;
  }
  double t1=0.0;
  double t2=0.0;
  double t3=0.0;
  t1 = MPI_Wtime();
  extern int method; 
Loci::ExecuteLoop(rp,exec_set,method,yMap,facts,local_facts[0],local_facts[1],local_compute1,local_compute2,inputs,outputs,stats);
  t2 = MPI_Wtime();
  delete [] yMap;        
  t3=t2-t1;
  std::cout<<Loci::MPI_rank<<'\t'<<t3<<'\t'<<stats[0]<<std::endl;
}

void dynamic_schedule_rule::Print(std::ostream &s) const {

 s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl ; 
  
}


}
