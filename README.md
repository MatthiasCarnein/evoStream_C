# evoStream - Evolutionary Stream Clustering Utilizing Idle Times

This is the implementation of an evolutionary stream clustering algorithm as proposed in our article in the Journal of Big Data Research.
The online component uses a simplified version of \code{DBSTREAM} to generate micro-clusters.
The micro-clusters are then incrementally reclustered using an evloutionary algorithm.
Evolutionary algorithms create slight variations by combining and randomly modifying existing solutions.
By iteratively selecting better solutions, an evolutionary pressure is created which improves the clustering over time.
Since the evolutionary algorithm is incremental, it is possible to apply it between observations, e.g. in the idle time of the stream.
Whenever there is idle time, we can call the \code{recluster} function of the reference class to improve the macro-clusters (see example).
The evolutionary algorithm can also be applied as a traditional reclustering step, or a combination of both.
In addition, this implementation also allows to evaluate a fixed number of generations after each observation.

## Usage

This is a port of the original code from the R-Package, available here: [https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream](https://wiwi-gitlab.uni-muenster.de/m_carn01/evoStream)
This implementation uses plain C++ without the interfaces to R. This should allow easier interfacing for other languages. It also makes use of C++11 features.

The main interfaces are the same as in the R-Package. An example how to read a comma-separated file and cluster the data points is shown below:

```cpp
int main()
{

  // Init evoStream
  EvoStream evo = EvoStream(0.05, 0.001, 100, 4, .8, .5, 100, 2*4, 1000);

  // Read CSV file (here comma-separated)
  std::ifstream in("data.csv");
  string line;
  while (getline(in, line)) {
      stringstream sep(line);
      string field;
      vector<double> fields;
      while (getline(sep, field, ',')) {
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
```
