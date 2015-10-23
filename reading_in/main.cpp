  #include <boost/graph/adjacency_list.hpp>
  #include <boost/graph/visitors.hpp>
  #include <boost/graph/breadth_first_search.hpp>
  #include <cstdlib>
  #include <boost/graph/bipartite.hpp>
  
  #include<string>
  #include<iostream>
  #include<fstream>
  #include<utility>
  #include<algorithm>

  #include"reaction.h"




int main (int argc, char* argv[]){
  using namespace boost;
  
  
  std::cout<<"Tests begin."<<std::endl;
  
  std::vector<reaction> reacVector;
  
  
  
           std::vector<Vertex> reacVList,compoundVList;
           ReactionNetwork lofasz;
  
  
           //reaction::readCompounds("compounds_list__4C_v3_2_2_ext_100.dat",lofasz,compoundVList);
           reaction::readCompounds("newsortedcompounds.txt",lofasz,compoundVList);
  
           std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
           //reaction::readReactions("reactions__4C_v3_2_2_ext_100.dat", reacVector,lofasz,reacVList,compoundVList);
           reaction::readReactions("newreactions.txt", reacVector,lofasz,reacVList,compoundVList);
           std::cout<<"length of the vector is: "<<reacVector.size()<<std::endl;
           reacVector[13].printReaction();
           reacVector[299].printReaction();
  
  
  
           graph_traits<ReactionNetwork>::vertex_iterator vi, vi_end;
           int count=0;
           for (boost::tie(vi, vi_end)=vertices(lofasz); vi!=vi_end; ++vi){
           count++;
           }
  
           std::cout<<"The graph has "<<count<<" vertices."<<std::endl;
  
  
           typedef graph_traits <ReactionNetwork> traits;
           typename traits::vertex_iterator vertex_iter, vertex_end;
  
             bool bipartiteee = boost::is_bipartite(lofasz);
  
           if (bipartiteee) {std::cout<<"The graph is bipartite"<<std::endl;}
           else {std::cout<<"The graph is not bipartite."<<std::endl;}
  
           reaction newreaction;
  
           reaction evennewer;
           evennewer.printReaction();
           evennewer=newreaction;
  
           newreaction.printReaction();
           evennewer.printReaction();
           std::cin.ignore();
  
  
  
  std::cout<<"Tests completed."<<std::endl;
  }

