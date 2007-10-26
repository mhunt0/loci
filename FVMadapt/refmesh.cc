//#include <fstream>
#include <iostream>
#include <string>
//#include <utility>
//#include <vector>
//#include <list>
//#include <algorithm>
#include <Loci.h>
//#include "hexcell.h"
#include "defines.h"
//#include <rpc/xdr.h>
//#include <rpc/rpc.h>
using std::string;
using std::cout;
using std::endl;
using std::cerr;
//using std::ofstream;
//using std::ifstream;
using std::ios;

namespace Loci{
  void parallelClassifyCell(fact_db &facts) ;
}


int main(int argc, char ** argv) {
  // Let Loci initialize itself.
  // This also gives Loci first dibs on the command line arguments.
  Loci::Init(&argc, &argv);
  

  // This is the name of the mesh file that you want to read in.
  // This may be overridden by the command line argument "-g file.xdr"
  //  string meshfile = "testGrid.xdr";
  string meshfile = "testGrid.vog";
  //This is the name of the refinement plan file
  string planFile = "out.plan";
  //This is the name of the output gridfile
  string outFile  = "out.xdr";
  
  // Here's where we parse out the command line arguments that are
  // relevant to this program.
  int j=1;
  //string pathname = "../gridgen/restore/mixedcell/";
  // string pathname = "/var/tmp/qxue/grid/";
   string pathname = "";
  bool restart = false;
  bool xdr = true;


  
    //print out help info
  if( (argc == 1)||(argc==2) ){
    if(Loci::MPI_rank == 0){
      cout<<"command line:" << endl;
      cout <<"refmesh <options> <filename> <options> <filename> ... "<< endl;
      cout << endl;
      cout << "options:" << endl;
      cout <<"-g <file> -- original grid file, refinement plans are based on this grid" << endl;
      cout <<"-r <file> -- input refinement plan file" <<endl;
      cout <<"-o <file> -- output grid file" << endl;
      cout << "-xdr -- output grid is in xdr format" << endl;
    }
    return 0;
  }
  
  for (int i=1; i<argc; i++) {
      // Let's look at the i'th argument
    string arg(argv[i]);

     if (arg == "-g" && (i+1) < argc) {
      // Replace the mesh filename with the next argument
      meshfile =  argv[++i];      
      }
    else if(arg == "-o" && (i+1) < argc){
      //replace the output filename with the next argument
      outFile =  argv[++i];
      }
    
    else if(arg == "-r" && (i+1) < argc){
      //replace the input refinement plan filename with the next argument
      planFile =  argv[++i];
      restart = true;
    }
    else if(arg == "-xdr" ){
      //output .xdr file
      
        xdr = true;
      }
      
      else{
        // Anything that we don't recognize gets recycled back into argv
        argv[j++] = argv[i];
      }
  }
  // Set the new size of the argument list to however many arguments
  // we put back in argv.
  argc = j;
  meshfile = pathname + meshfile;
  outFile = pathname + outFile;
  planFile = pathname + planFile;
  



 // Setup the fact database.
  fact_db facts;
  
  param<bool> restart_par;
  *restart_par = restart;
  facts.create_fact("restart_par", restart_par);
  
  
  param<std::string> planfile_par ;
  *planfile_par = planFile;
  facts.create_fact("planfile_par",planfile_par) ;

  param<std::string> outfile_par ;
  *outfile_par = outFile;
  facts.create_fact("outfile_par",outfile_par) ;

  
  // Setup the rule database.
  // Add all registered rules.
  rule_db rules;
  rules.add_rules(global_rule_list);
  
 
  std::cout <<"reading in meshfile" << std::endl;
  // Read in the mesh file.  Setup Loci datastructures
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    std::cerr << "unable to read grid file '" << meshfile << "'" << std::endl ;
    Loci::Abort() ;
  }
  
  createLowerUpper(facts) ;
  createEdgesPar(facts) ;
   Loci:: parallelClassifyCell(facts);
 
 cerr << currentMem() << " after read in grid " <<  Loci::MPI_rank << endl;

  
 //  std::cout << "Making query." << std::endl;
//   if(!Loci::makeQuery(rules, facts, "num_inner_nodes,num_fine_cells")) {
//     std::cerr << "query failed!" << std::endl;
//     Loci::Abort();
//   }

 

//   Loci::storeRepP num_inner_nodes;
//   num_inner_nodes = facts.get_variable("num_inner_nodes");
//   facts.create_fact("num_inner_nodes_copy", num_inner_nodes);
  
//   Loci::storeRepP num_fine_cells;
//   num_fine_cells = facts.get_variable("num_fine_cells");
//   facts.create_fact("num_fine_cells_copy", num_fine_cells);
  
//     cerr << currentMem() << " after num_inner_nodes " <<  Loci::MPI_rank << endl;

//     std::cout << "Making query." << std::endl;
//   if(!Loci::makeQuery(rules, facts, "inner_nodes,fine_faces")) {
//     std::cerr << "query failed!" << std::endl;
//     Loci::Abort();
//   }

//   Loci::storeRepP inner_nodes;
//   inner_nodes = facts.get_variable("inner_nodes");
//   facts.create_fact("inner_nodes_copy", inner_nodes);
  
//   Loci::storeRepP fine_faces;
//   fine_faces = facts.get_variable("fine_faces");
//   facts.create_fact("fine_faces_copy", fine_faces);
  
//   cerr << currentMem() << " after inner_nodes " <<  Loci::MPI_rank << endl;
  
//   std::cout << "Making query." << std::endl;
//   if(!Loci::makeQuery(rules, facts, "npnts")) {
//     std::cerr << "query failed!" << std::endl;
//     Loci::Abort();
//   }
 

//     Loci::storeRepP npnts;
  
//     npnts = facts.get_variable("npnts");
  
//    facts.create_fact("npnts_copy", npnts);

//    // cout << "npnts: " << num_points <<" myID : " << Loci::MPI_rank << std::endl;
  
//   std::cout << "Making query." << std::endl;
//   if(!Loci::makeQuery(rules, facts, "ncells")) {
//     std::cerr << "query failed!" << std::endl;
//     Loci::Abort();
//   }
 
//   Loci::storeRepP ncells;
//   ncells = facts.get_variable("ncells");
//   facts.create_fact("ncells_copy", ncells);
 
  
//  std::cout << "Making query." << std::endl;
//   if(!Loci::makeQuery(rules, facts, "nfaces")) {
//     std::cerr << "query failed!" << std::endl;
//     Loci::Abort();
//   }

//    Loci::storeRepP nfaces;
//    //  param<int> nfaces;
//   nfaces = facts.get_variable("nfaces");
//   facts.create_fact("nfaces_copy", nfaces);
 
  
//   cerr << currentMem() << " after nfaces " <<  Loci::MPI_rank << endl;
  
  if(!Loci::makeQuery(rules, facts, "pos_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }
 

   if(!Loci::makeQuery(rules, facts, "node_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
   }

   if(!Loci::makeQuery(rules, facts, "face_output")) {
     std::cerr << "query failed!" << std::endl;
     Loci::Abort();
    }
  // Tell Loci to cleanup after itself and exit gracefully.
  Loci::Finalize();
  
}
