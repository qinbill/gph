void linscan_query(UINT32 *counter, UINT32 *res, UINT8 *codes, UINT8 *queries, int N, UINT32 NQ, int B, int R, \
		   int dim1codes, int dim1queries);


void linscan_rangequery(UINT32 *results, UINT32 *numres, UINT8 *codes, UINT8 *queries, int N, UINT32 NQ, int B,
                        int dim1codes, int dim1queries, int tau);
