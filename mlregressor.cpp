#include "mlregressor.h"

#include <bitset>


uint32_t transfer_bits_to_string(std::vector<uint8_t>& vec, uint32_t start, uint32_t end)
{
	uint32_t result = 0;
	int p = end - start - 1;
	for(uint32_t i = start; i < end; i++)
	{
		result += (uint32_t)(pow(2.0, p--)) * vec[i];				
		//result += (uint32_t)(pow(2.0, i - start)) * vec[i];
	}
	return result;
}


// Transfer one feature vector to a high-dimensions vector
void mlregressor::transfer_to_high_dimen(std::vector<uint8_t>& feature_vec, int pid, int queryid)
{
	uint8_t* p = feature_vec.data();
	uint64_t size = feature_vec.size();
	auto generator = make_combination_generator(p, p + size, bit_num);

	std::vector<uint8_t> t;
	// Allocate space 
	t.resize((int)(this->D - this->d) / (uint32_t)pow(2.0, bit_num) * this->bit_num);

	std::vector<uint8_t>::iterator it = t.begin();
	while(generator(it))
	{
		it += bit_num;
	}
	
	feature_type ft_single(this->D, 1);
	//printf("Dimension of a single feature vector : %d\n", this->D);
	// the dimension of new feature vector
	uint32_t ft_count = 0;
	// Insert the old features
	for(int i = 0; i < this->d; i++)
	{
		ft_single(ft_count++) = feature_vec[i];
	}

	// Add the new features
	uint32_t base = (uint32_t)pow(2.0, bit_num);
	for(uint32_t i = 0; i < t.size(); i += bit_num)
	{
		uint32_t r = transfer_bits_to_string(t, i, i + bit_num);	
		// Mark the corresponding feature
		for(uint32_t j = 0; j < base; j++)
		{
			ft_single(ft_count++) = (r == j ? 1 : 0);			

		}
	}
	//printf("count : %d\n", ft_count);
	//this->ft[pid] = ft_single;
	ft[queryid].push_back(ft_single);
}


void mlregressor::gen_feature_vec(uint64_t chunk, int pid, int queryid) {
	std::bitset<64> bit_partition(chunk);
	std::vector<uint8_t> ft_old;
	for(int i = this->d - 1; i >= 0; i--)
	{
		ft_old.push_back(bit_partition[i]);
	}
	this->transfer_to_high_dimen(ft_old, pid, queryid);

}


void mlregressor::gen_feature_vecs(uint64_t* chunks, int queryid)
{
	for(int i = 0; i < this->M; i++)
	{
		gen_feature_vec(chunks[i], i, queryid);
	}	
}




// Search
double mlregressor::search(int script, int error, int queryid) {
	if(error < 0)
		return 0;
        if(error > this->d) error = this->d;
	if(temp[queryid][script][error] < 0) {
		temp[queryid][script][error] =  this->df[script][error](this->ft[queryid][script]);
	}
	return exp(temp[queryid][script][error]); 
	//return this->df[script][error](this->ft[queryid][script]);
}



// Calcualte the value of selecting K from N
int mlregressor::cal_combination(int n, int k)
{
	if(k == 0)
	       	return 1;
	return (n * cal_combination(n - 1, k - 1)) / k;
}


// Just for test
void mlregressor::fillallpredicts(int NQ, int M, int maxTau)
{
	for(int i = 0; i < NQ; i++)
	  for(int j = 0; j < M; j++)
	   for(int z = 0; z < maxTau; z++)
		this->temp[i][j][z] = this->df[j][z](this->ft[i][j]);
}
