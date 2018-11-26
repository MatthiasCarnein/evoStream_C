#include <sstream>
#include <fstream>
#include <vector>

#include "evoStream.hpp"

int main(){

  // Init evoStream 
  EvoStream evo = EvoStream(0.05, 0.001, 100, 4, .8, .001, 100, 2*4, 1000);

  // Read CSV file (here comma-separated)
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
      evo.recluster(1); // evaluate 1 generations after every observation. This can be adapted to the available time
  }

  std::vector< std::vector<double> > res = evo.get_macroclusters(); // reclustering evaluates 1000 additional generations (parameter)

  // Print macro-clusters
  for(unsigned int row=0; row<res.size(); row++){
    for(unsigned int col=0; col<res[0].size(); col++){
      std::cout << res[row][col] << ", ";
    }
    std::cout << std::endl;
  }
  return 0;

}