#ifndef MLREGRESSOR_H__
#define MLREGRESSOR_H__

#include <stdio.h>
#include <cmath>
#include <stdexcept>


// Add dlib library
#include <dlib/svm.h>

#define MaxPartition 100
#define MaxTau 65
#define MaxQueries 20000

using namespace dlib;

typedef matrix<double> feature_type;

//typedef radial_basis_kernel<feature_type> kernel_type;
typedef linear_kernel<feature_type> kernel_type;

// Combinations of bits
template<class iterator_type>
class combination_generator {
    iterator_type first, last;
    std::vector<bool> use;
    unsigned r;
    typedef typename std::iterator_traits<iterator_type>::value_type element_type;
public:
    combination_generator(iterator_type first_, iterator_type last_, unsigned r_) : first(first_), last(last_) , r(r_)
    {
        use.resize(std::distance(first, last), false);
        if (r > use.size())
            throw std::domain_error("can't select more elements than exist for combination");
        std::fill(use.end()-r, use.end(), true);
    }
    template<class output_iterator>
    bool operator()(output_iterator result) 
    {
        iterator_type c=first;
        for (unsigned i = 0; i<use.size(); ++i,++c) {
            if (use[i]) 
                *result++ = *c;
        }
        return std::next_permutation(use.begin(), use.end());
    }
};

template<class iterator_type>
combination_generator<iterator_type> make_combination_generator(iterator_type first, iterator_type last, unsigned r)
{return combination_generator<iterator_type>(first, last, r);}



/* ml class for predicting candidate size
 * by giving a chunk and a threshold
 */
class mlregressor {
private:
	// Calcualte the value of selecting K from N
	int cal_combination(int n, int k);


public:
	// number of regressor, which is also the 
	// maximum threshold currently.
	int numreg;
	// combination of bits
	int bit_num;
	// dimension of feature vector
	int D;
	// number of bis in chunk
	int d;
	// decision function of dlib
	decision_function<kernel_type> df[MaxPartition][MaxTau];
	// feature vectors
	std::vector<feature_type> ft[MaxQueries];
	int M;

	// temp records to store
	double temp[MaxQueries][MaxPartition][MaxTau];

	/* maxtau:	the maximum threshold
	 * prefixfile:	the prefix file of regressors
	 * d:		the number of bits in single chunk
	 * M:		the number of chunks
	 */
	mlregressor(int maxtau, char* prefixfile, int d, int M, int bit_num) {
		this->numreg = std::min(maxtau + 1, d);
		/*
		// Malloc the decision_function
		df = new decision_function<kernel_type>*[M];
		for(int i = 0; i < M ; i++)
		{
			df[i] = new decision_function<kernel_type>[this->numreg];
		}
		*/
		// Load the decision_functions from file
		char regfile[256];
		for(int j = 0; j < M; j++) {
			for(int i = 0; i < this->numreg; i++) {
				sprintf(regfile, "%s_P=%d_t=%d.reg", prefixfile, j, i);
				deserialize(regfile) >> this->df[j][i];
			}
		}
		this->bit_num = bit_num;
		this->d = d;
		this->M = M;

		// dimension of the feature vector
		this->D = d + this->cal_combination(d, bit_num) * (uint32_t)pow(2.0, bit_num);
		// Malloc the feature vectors
		//this->ft = new feature_type[M];
		// Assign the dimension of each feature vector
		
		/* useless for now
		for(int i = 0; i < M; i++)
		{
			feature_type t(this->D, 1);
			this->ft[i] = t;

			//this->ft[i](this->D, 1);
		}
		*/
		// Initialize temp
		for(int i = 0; i < MaxQueries; i++)
		 for(int j = 0 ; j < MaxPartition; j++)
		  for(int z = 0; z < MaxTau; z++)
			temp[i][j][z] = -1.0;
	}

	~mlregressor() {
		if (this->df != NULL) {
			//delete[] df;
		}
	}

	void transfer_to_high_dimen(std::vector<uint8_t>& feature_vec, int pid, int queryid);
	void gen_feature_vec(uint64_t chunk, int pid, int queryid);
	void gen_feature_vecs(uint64_t* chunks, int queryid);
	double search(int script, int error, int queryid);
	void fillallpredicts(int NQ, int M, int maxTau);

};

#endif
