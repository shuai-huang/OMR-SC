#ifndef PGD_H
#define PGD_H

#include "UVT.h"
#include <thread>
#include <cmath>

// objective function should contain the term corresponding to zero values
// when computing gradient, even if p(r) = 0, we still need to compute it, right?
class PGD:public UVT
{
	private:
	
	public:
		PGD();

        virtual void ComputeCorObjFunMuti( vector<DtPair> lm_partition_block, VectorXd* est_cor_feat_seq_pt );
        virtual double ComputeCorObjFun();
        
        virtual void ComputeRadialObjFunMuti( vector<int> radial_partition_block );
        virtual double ComputeRadialObjFun();

        virtual void ComputeProjObjFunMuti( vector<DtPair> proj_partition_block );
        virtual double ComputeProjObjFun();
        
        virtual void ComputeCorGradientMuti( vector<DtPair> lm_partition_block, VectorXd* gradient_val_seq_pt );
        virtual void ComputeCorGradient();
        
        virtual void ComputeRadialGradientMuti( vector<int> radial_partition_block, VectorXd* gradient_val_seq_pt );
        virtual void ComputeRadialGradient();

        virtual void ComputeProjGradientMuti( vector<DtPair> proj_partition_block, VectorXd* gradient_val_seq_pt );
        virtual void ComputeProjGradient();
		
		virtual ~PGD();
};

#endif // PGD_H
