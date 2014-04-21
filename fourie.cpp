
#include <stdio.h>  /* for sprintf */
#include <stdlib.h> /* for strtod */
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

const float pi = acos(-1);

int numTerms;
int numVars;
int Order;
vector<vector<double> > multipliers;
vector<vector<double> > obsRanges;


/********** GDK's Fourier Code (Ported from Java) **************/
/**                                                                                                                                                                                                                                         * This method iterates through a coefficient vector                                                                                                                                                                                        * up to a given degree. (like counting in a base of                                                                                                                                                                                        * that degree).                                                                                                                                                                                                                            *                                                                                                                                                                                                                                          * @param cthe coefficient vector.                                                                                                                                                                                                          * @param NVariablesthe number of variables in c.                                                                                                                                                                                           * @param Degreethe degree up to which to increment.                                                                                                                                                                                        */
void Iterate(int *c, int NVariables, int Degree)
{
  c[NVariables - 1] = c[NVariables-1] + 1;

  if(c[NVariables - 1] > Degree)
    {
      if(NVariables > 1)
        {
          c[NVariables - 1]  = 0;
          Iterate(c, NVariables - 1, Degree);
        }
    }
}

/*
 * Compute the full Fourier Basis coefficient matrix for 
 * a given number of variables up to a given order.
 *
 * @param nvarsthe number of variables.
 * @param orderthe highest coefficient value for any individual variable.  
 * @returna two dimensional array of doubles. The first dimension length is
 * the number of basis functions, and the second is the number of state variables. 
 */
void computeFourierCoefficients(int nvars, int order) {
  int nterms = (int)pow(order + 1.0, nvars);
  numTerms = nterms;
  numVars = nvars;
  Order = order;

  int pos = 0;
  multipliers.resize(nterms);
  for (int i=0; i<nterms; i++)
    multipliers[i].resize(nvars);
  int *c = new int[nvars];
  for(int j = 0; j < nvars; j++)
    c[j] = 0;

  do
    {
      for(int k = 0; k < nvars; k++)
        {
          multipliers[pos][k] = c[k];
        }

      pos++;
      // Iterate c                                                                                                                                                                                                                         
      Iterate(c, nvars, order);
    }
  while(c[0] <= order);
}

/**
 * Scale a state variable to between 0 and 1. 
 * (this is required for the Fourier Basis). 
 *
 * @param valthe state variable. 
 * @param posthe state variable number.
 * @returnthe normalized state variable. 
 */
double scale(double val, int pos)
{
  //  cout << pos << "," << val << " " << obsRanges[pos][0] << "," << obsRanges[pos][1] << " : " << (val - obsRanges[pos][0])/(obsRanges[pos][1] - obsRanges[pos][0]) << endl;
  return (val - obsRanges[pos][0]) / (obsRanges[pos][1] - obsRanges[pos][0]);
}


/**
 * Compute the feature vector for a given state. 
 * This is achieved by evaluating each Fourier Basis function
 * at that state.
 *
 * @param sthe state in question.
 * @return a vector of doubles representing each basis function evaluated at s.
 */
vector<double> computeFeatures(vector<double> features)
{
  vector<double> phi(numTerms);
  for(int pos = 0; pos < numTerms; pos++)
    {
      double dsum = 0;
      for(int j = 0; j < numVars; j++)
        {
          double sval = scale(features[j], j);
          dsum += sval*multipliers[pos][j];
        }

      phi[pos] = cos((pi) * dsum);
    }

  return phi;
}

  typedef boost::unordered_map<std::string, int> map;  

int main(int argc, char **argv) {
  char input[4096];
  char_separator<char> sep(", "); // lets start with csv and move to VW input later
  char_separator<char> vwsep(":");
  double order = lexical_cast<int>(*++argv);
  //  vector<string> features;
  int vw_offset = 0;
  bool skip_normal = false;
  boost::unordered_map<string, int> features;
  // Default to CSV conversion
  // Take parameter flag: --vw 
  // to switch to vw input based
  while(*++argv) {
      if(strcmp(*argv, "--vw") == 0)
          vw_offset = 1;
      else if(strcmp(*argv, "--nonorm") == 0)
  	  skip_normal = true;
      else {
	cout << "Unknown argument: " << *argv << endl;
	exit(1);
      }
  }
  
  //cout << "Fourier Order: " << order << endl;
  int line_count = 0;
  bool translating = false;

  while(!cin.eof()) {
    bool line_has_data = false;
    cin.getline(input,4096);

    string line(input);
    tokenizer< char_separator<char> > tokens(line, sep);
    vector<double> vanilla_features;
    int feature_counter = 0;

    if(vw_offset > 0 && line_count > 0)  {
      vanilla_features.resize(features.size());
    }
    BOOST_FOREACH (const string& t, tokens) {
      try {
	if(line_count == 0) {
	  obsRanges.push_back(vector<double>(2));
	  obsRanges[feature_counter][0] = 0.0;
	  obsRanges[feature_counter][1] = 1.0;
	} 
	// min then max then features
	if(line_count < vw_offset) { // first line, and we ARE doing VW inputs
	  features[t] = feature_counter;
	} else if(line_count == vw_offset && !skip_normal) { // ready for the min-max values
	  obsRanges[feature_counter][0] = lexical_cast<double>(t);
	} else if(line_count == vw_offset+1 && !skip_normal) {
	  obsRanges[feature_counter][1] = lexical_cast<double>(t);
	} else {
	  translating = true;
	  if(vw_offset > 0) {
	    boost::tokenizer< char_separator<char> > vwtoken(t, vwsep);
	    boost::tokenizer< char_separator<char> >::iterator beg=vwtoken.begin();
	    string fname(*beg);
	    if(features.find(*beg) != features.end()) {
	      int findex = features.at(*beg);
	      vanilla_features[findex] = lexical_cast<double>(*++beg);
	    } else {
	      cout << t << " ";
	      feature_counter--;
	      line_has_data = true;
	    }
	  } else {
	    vanilla_features.push_back(lexical_cast<double>(t));      
	  }
	  if(feature_counter >= 0 && (vanilla_features[feature_counter] < obsRanges[feature_counter][0] || vanilla_features[feature_counter] > obsRanges[feature_counter][1])) {
	    cout << "ERROR: feature " << feature_counter << " is out of supplied range"<<endl;
	    exit(1);
	  } //else cout << vanilla_features[feature_counter] << endl;
	}
      }
      catch(bad_lexical_cast&)
	{}

      feature_counter++;
    }
    //    cout << "Counts: " << feature_counter << " " << line_count << " " << vanilla_features.size() << endl;
    //head ../uk3day/train.vw | sed '1d' | ./fourie 3 --vw
    if(line_count == 0) {
      if(vw_offset > 0 && features.size() == 0) {
	cout << "ERROR: feature labels not provided" << endl;
	exit(1);
      }
      //      cout << feature_counter << endl;
      //cout << "Num Features: " << feature_counter << endl;
      computeFourierCoefficients(feature_counter, order);
      //      cout << "Num Fourier Terms: " << numTerms << endl;
    }
      
    if(translating && feature_counter > 0) {
      //cout << "Computing features... " <<endl;
      vector<double> fourie_features = computeFeatures(vanilla_features);
      if(vw_offset > 0) {
	cout <<"FOURIER0:" << fourie_features[0] << " ";
	for(int i=1; i<(int)fourie_features.size(); i++) {
	  cout << "FOURIER" << i << ":" << fourie_features[i] << " ";
	}
	cout << endl;
      } else {
	cout << fourie_features[0];
	for(int i=1; i<(int)fourie_features.size(); i++) {
	  cout << "," << fourie_features[i];
	}
	cout << endl;
      }
    } else if(line_has_data) {
      cout << endl;
    }
    line_count++;
  }
  return 0;
}
