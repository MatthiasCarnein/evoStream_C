#include <sstream>
#include <fstream>
#include <vector>

#include "evoStream.hpp"

int main(){


  // Main Interface:
  EvoStream evo = EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 250); // init
  std::vector<double> observation { 10.0, 20.0, 30.0 }; // read observation
  evo.cluster(observation); // cluster new observation
  evo.get_microclusters();
  evo.get_microweights();
  evo.get_macroclusters();
  evo.get_macroweights();
  evo.recluster(100); // evaluate 100 more macro solutions
  evo.microToMacro();


  // Full Example: Read CSV file (here: comma-separated, numeric values)
  evo = EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 250);

  std::ifstream in("data.csv");
  std::string line;
  while (std::getline(in, line)) {
      std::stringstream sep(line);
      std::string field;
      std::vector<double> fields;
      while (std::getline(sep, field, ',')) {
          fields.push_back(stod(field));
      }
      evo.cluster(fields); // insert observation
      evo.recluster(1); // evaluate 1 generation after every observation. This can be adapted to the available time
  }


  // Get micro-clusters
  std::vector< std::vector<double> > micro = evo.get_microclusters();
  std::vector<double> microWeights = evo.get_microweights();

  std::cout << "Micro Clusters" << std::endl;
  for(unsigned int row=0; row < micro.size(); row++){
    for(unsigned int col=0; col < micro[0].size(); col++){
      std::cout << micro[row][col] << " ";
    }
    std::cout << "(weight: " << microWeights[row] << ")" <<std::endl;
  }


  // Get macro-clusters (here: performs an additional 250 reclustering steps, see parameter)
  std::vector< std::vector<double> > macro = evo.get_macroclusters(); // reclustering 
  std::vector<double> macroWeights = evo.get_macroweights(); 

  std::cout << "\n\nMacro Clusters" << std::endl;
  for(unsigned int row=0; row < macro.size(); row++){
    for(unsigned int col=0; col < macro[0].size(); col++){
      std::cout << macro[row][col] << " ";
    }
    std::cout << "(weight: " << macroWeights[row] << ")" <<std::endl;
  }


  std::cout << "\n\nAssignment of Micro Clusters to Macro Clusters" << std::endl;
  std::vector<int> microToMacro = evo.microToMacro();
  for(unsigned int i=0; i < microToMacro.size(); i++){
    std::cout << "Micro " << i << " -> " << "Macro " << microToMacro[i] << std::endl;
  }


  return 0;

}