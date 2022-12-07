#include "PGD.h"

// Could remove the columns that are all zeros and speed up the computation

PGD::PGD(){}



void PGD::ComputeCorObjFunMuti( vector<DtPair> lm_partition_block, VectorXd* est_cor_feat_seq_pt ) {

    VectorXd est_cor_feat_seq_tmp = VectorXd::Zero(num_cor_feat);
    for (int i=0; i<lm_partition_block.size(); i++) {
        DtPair lm_pair = lm_partition_block[i];
        int l=lm_pair[0];
        int m=lm_pair[1];

        //cout<<"m1"<<" "<<l<<" "<<m<<endl;

        //if (l==0) {
        //    continue;
        //}
        
        map<int, double> density_alpha_sh_val; // the radial index, the sum of multiplication of density, alpha, sh values
        // note that only r>=1 has alpha values
        //////////////////
        // I should save this so that i don't have to compute it again while computing the gradients !!!
        //////////////////
        for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {
            if (r_candidate.find(r_idx)!=r_candidate.end()) {
                double density_alpha_sh_val_tmp = (radial_density_val[r_idx].array() * radial_alpha_sh_seq[{r_idx, l, m}].array() ).sum();
                density_alpha_sh_val[r_idx]=density_alpha_sh_val_tmp;
            } else {
                density_alpha_sh_val[r_idx]=0;
            }
            //if ((l==1)&&(m==-1)) {
            //    cout<<r<<" "<<density_alpha_sh_val[r]<<endl;
            //}
        }

        //cout<<"m2"<<" "<<radial_map.size()<<endl;
        //for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {
        //    cout<<r_idx<<" ";
        //}
        //cout<<endl;
        
        for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {
            //cout<<r_idx<<" ";
            if (r_candidate.find(r_idx)!=r_candidate.end()) {
            if (cor_feat_all_loc.find({l,m,r_idx})!=cor_feat_all_loc.end()) {
                int cor_feat_all_loc_val = cor_feat_all_loc[{l,m,r_idx}];
                est_cor_feat_seq_tmp(cor_feat_all_loc_val) = est_cor_feat_seq_tmp(cor_feat_all_loc_val) + density_alpha_sh_val[r_idx];  // does not multiply by (M_PI/2) here
            }
            }
        }
        //cout<<endl;
        //cout<<"m3"<<endl;

    }
    
    (*est_cor_feat_seq_pt) = est_cor_feat_seq_tmp;
}



double PGD::ComputeCorObjFun() {
    
    vector<VectorXd> est_cor_feat_seq;
    for (int i=0; i<num_lm_thread_assign; i++) {
        VectorXd est_cor_feat_seq_tmp = VectorXd::Zero(num_cor_feat);
        est_cor_feat_seq.push_back(est_cor_feat_seq_tmp);
    }
    
    if (num_lm_thread_assign>1) {
        thread *multi_thread = new thread[num_lm_thread_assign-1];
        for (int i=0; i<num_lm_thread_assign-1; i++) {
            vector<DtPair> lm_partition_idx_tmp = lm_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeCorObjFunMuti, this,  lm_partition_idx_tmp, &est_cor_feat_seq[i]);
        }
        
        ComputeCorObjFunMuti( lm_partition_idx[num_lm_thread_assign-1], &est_cor_feat_seq[num_lm_thread_assign-1]);
        
        for (int i=0; i<num_lm_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeCorObjFunMuti( lm_partition_idx[num_lm_thread_assign-1], &est_cor_feat_seq[num_lm_thread_assign-1]);
    }
    
    est_cor_feat_all = VectorXd::Zero(num_cor_feat);
    for (int i=0; i<num_lm_thread_assign; i++) {
        est_cor_feat_all += est_cor_feat_seq[i];
    }
   
    // rearrange est_cor_feat_all to est_cor_feat_mat_all
    for (int l=0; l<=L_max; l++) {
        MatrixXd cor_feat_mat_all_tmp = MatrixXd::Zero(num_radial_cor_feat,2*l+1);
        for (int m=-l; m<=l; m++) {
            for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {
                DtTriplet cor_feat_all_loc_key = {l,m,r_idx};
                int cor_feat_all_loc_val = cor_feat_all_loc[cor_feat_all_loc_key];
                cor_feat_mat_all_tmp(r_idx,m+l) = est_cor_feat_all(cor_feat_all_loc_val);
            }
        }
        est_cor_feat_mat_all[l] = cor_feat_mat_all_tmp;
    }
    
    // compute est_fourier_cor_feat_mat_all
    for (int l=0; l<=L_max; l++) {
        MatrixXd est_cor_feat_mat_all_tmp = est_cor_feat_mat_all[l];
        MatrixXd radial_lgwt_mat_tmp = radial_lgwt_mat[l];

        MatrixXd est_fourier_cor_feat_mat_all_tmp = radial_lgwt_mat_tmp.adjoint() * est_cor_feat_mat_all_tmp;
        est_fourier_cor_feat_mat_all[l] = est_fourier_cor_feat_mat_all_tmp;
    }
    
    // rearrange est_fourier_cor_feat_mat_all to est_fourier_cor_feat_all
    est_fourier_cor_feat_all = VectorXd::Zero(num_fourier_cor_feat);
    for (int l=0; l<=L_max; l++) {
        MatrixXd est_fourier_cor_feat_mat_all_tmp = est_fourier_cor_feat_mat_all[l];
        for (int m=-l; m<=l; m++) {
            for (int r_idx=0; r_idx<fourier_radial_map.size(); r_idx++) {
                DtTriplet fourier_cor_feat_all_loc_key = {l,m,r_idx};
                int fourier_cor_feat_all_loc_val = fourier_cor_feat_all_loc[fourier_cor_feat_all_loc_key];
                
                est_fourier_cor_feat_all(fourier_cor_feat_all_loc_val) = est_fourier_cor_feat_mat_all_tmp(r_idx,m+l);
            }
        }
    }
    
    
    VectorXd fourier_cor_feat_diff = fourier_cor_feat_all_transform-est_fourier_cor_feat_all;
    
    double fourier_cor_obj = (fourier_cor_feat_diff.array() * fourier_cor_feat_diff.array()).sum();

    return fourier_cor_obj;
}



void PGD::ComputeRadialObjFunMuti( vector<int> radial_partition_block ) {

    for (int i=0; i<radial_partition_block.size(); i++) {
        int r_idx=radial_partition_block[i];
        
        if (r_candidate.find(r_idx)!=r_candidate.end()) {
            est_radial_feat_all(radial_feat_all_loc[r_idx]) = (radial_density_val[r_idx].array() * radial_beta_seq[r_idx].array() ).sum();
        }
    }
    
}



double PGD::ComputeRadialObjFun() {
    
    est_radial_feat_all = VectorXd::Zero(num_radial_feat);
    
    if (num_radial_thread_assign>1) {
        thread *multi_thread = new thread[num_radial_thread_assign-1];
        for (int i=0; i<num_radial_thread_assign-1; i++) {
            vector<int> radial_partition_idx_tmp = radial_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeRadialObjFunMuti, this, radial_partition_idx_tmp);
        }
        
        ComputeRadialObjFunMuti( radial_partition_idx[num_radial_thread_assign-1] );
        
        for (int i=0; i<num_radial_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeRadialObjFunMuti( radial_partition_idx[num_radial_thread_assign-1] );
    }
    
    //VectorXd radial_feat_diff = radial_feat_all-est_radial_feat_all;
    
    //double radial_obj = radial_feat_diff.array().square().sum();
    double radial_obj = 0;
    for (int r_idx=0; r_idx<radial_feat_idx.size(); r_idx++) {
        int radial_feat_loc_tmp = radial_feat_all_loc[r_idx];
        radial_obj = radial_obj + pow(radial_lgwt_weights(r_idx)*(radial_feat_all(radial_feat_loc_tmp)-est_radial_feat_all(radial_feat_loc_tmp)), 2);
    }



    return radial_obj;
}



void PGD::ComputeProjObjFunMuti( vector<DtPair> proj_partition_block ) { 

    for (int i=0; i<proj_partition_block.size(); i++) {
        DtPair proj_idx=proj_partition_block[i];
    
        if (voxel_proj_idx.find(proj_idx)!=voxel_proj_idx.end()) {
            vector<int> voxel_group = voxel_proj_idx[proj_idx];
            int num_voxel_group = voxel_group.size();
            double est_proj_feat_all_tmp = 0;
            for (int i=0; i<num_voxel_group; i++) {
                est_proj_feat_all_tmp = est_proj_feat_all_tmp + density_val(voxel_group[i]);
            }
            est_proj_feat_all(proj_feat_all_loc[proj_idx]) = est_proj_feat_all_tmp;
        }
    }   
    
}



double PGD::ComputeProjObjFun() {
    
    est_proj_feat_all = VectorXd::Zero(num_proj_feat);
    
    if (num_proj_thread_assign>1) {
        thread *multi_thread = new thread[num_proj_thread_assign-1];
        for (int i=0; i<num_proj_thread_assign-1; i++) {
            vector<DtPair> proj_partition_idx_tmp = proj_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeProjObjFunMuti, this, proj_partition_idx_tmp);
        }
    
        ComputeProjObjFunMuti( proj_partition_idx[num_proj_thread_assign-1] );
    
        for (int i=0; i<num_proj_thread_assign-1; i++) {
            multi_thread[i].join();
        }
    
        delete [] multi_thread;
    
    } else {
        ComputeProjObjFunMuti( proj_partition_idx[num_proj_thread_assign-1] );
    }   
    
    VectorXd proj_feat_diff = proj_feat_all-est_proj_feat_all;
    
    double proj_obj = proj_feat_diff.array().square().sum();

    return proj_obj;
}




void PGD::ComputeCorGradientMuti( vector<DtPair> lm_partition_block, VectorXd* gradient_val_seq_pt ) {

    VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
    for (int i=0; i<lm_partition_block.size(); i++) {
        DtPair lm_pair = lm_partition_block[i];
        int l=lm_pair[0];
        int m=lm_pair[1];

        //if (l==0) {
        //    continue;
        //}
        
        MatrixXd fourier_radial_gradient_mut = MatrixXd::Zero(num_fourier_radial_cor_feat,1);
        for (int r_idx=0; r_idx<fourier_radial_map.size(); r_idx++) {
            int fourier_cor_feat_all_loc_val = fourier_cor_feat_all_loc[{l,m,r_idx}];
            fourier_radial_gradient_mut(r_idx,0) = (est_fourier_cor_feat_all(fourier_cor_feat_all_loc_val) - fourier_cor_feat_all_transform(fourier_cor_feat_all_loc_val));
        }
        fourier_radial_gradient_mut = radial_lgwt_mat[l] * fourier_radial_gradient_mut;
       
        // change to radial_map,  because fourier_radial_gradient_mut is num_radial_cor_feat by 1 now
        for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {
            if (r_candidate.find(r_idx)!=r_candidate.end()) {
        
                double fourier_radial_gradient_mut_tmp = fourier_radial_gradient_mut(r_idx,0);
                
                VectorXd radial_gradient_val_seq_tmp_r = radial_alpha_sh_seq[{r_idx, l, m}] * fourier_radial_gradient_mut_tmp;
                
                // add the gradients corresponding to r
                vector<int> radial_idx_r = voxel_radial_idx[r_idx];
                for (int r_idx2=0; r_idx2<radial_idx_r.size(); r_idx2++) {
                    int r_loc = radial_idx_r[r_idx2];
                    gradient_val_seq_tmp(r_loc) = gradient_val_seq_tmp(r_loc) + radial_gradient_val_seq_tmp_r[r_idx2];
                }
            }

        }

    }
    
    (*gradient_val_seq_pt) = gradient_val_seq_tmp;
    
}

void PGD::ComputeCorGradient() {

    vector<VectorXd> gradient_val_seq;
    for (int i=0; i<num_lm_thread_assign; i++) {
        VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
        gradient_val_seq.push_back(gradient_val_seq_tmp);
    }
    
    if (num_lm_thread_assign>1) {
        thread *multi_thread = new thread[num_lm_thread_assign-1];
        for (int i=0; i<num_lm_thread_assign-1; i++) {
            vector<DtPair> lm_partition_idx_tmp = lm_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeCorGradientMuti, this, lm_partition_idx_tmp, &gradient_val_seq[i]);
        }
        
        ComputeCorGradientMuti( lm_partition_idx[num_lm_thread_assign-1], &gradient_val_seq[num_lm_thread_assign-1]);
        
        for (int i=0; i<num_lm_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeCorGradientMuti( lm_partition_idx[num_lm_thread_assign-1], &gradient_val_seq[num_lm_thread_assign-1]);
    }
    
    cor_gradient_val = VectorXd::Zero(num_loc);
    for (int i=0; i<num_lm_thread_assign; i++) {
        cor_gradient_val += gradient_val_seq[i];
    }

}



void PGD::ComputeRadialGradientMuti( vector<int> radial_partition_block, VectorXd* gradient_val_seq_pt ) {

    VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
    for (int i=0; i<radial_partition_block.size(); i++) {
        int r_idx=radial_partition_block[i];
        
        if (r_candidate.find(r_idx)!=r_candidate.end()) {
            double radial_gradient_mut_tmp = radial_lgwt_weights(r_idx) * ( est_radial_feat_all(radial_feat_all_loc[r_idx]) - radial_feat_all(radial_feat_all_loc[r_idx]) );
            VectorXd radial_gradient_val_seq_tmp = radial_beta_seq[r_idx]*radial_gradient_mut_tmp;
            
            vector<int> radial_idx_r = voxel_radial_idx[r_idx];
            for (int r_idx2=0; r_idx2<radial_idx_r.size(); r_idx2++) {
                int r_loc = radial_idx_r[r_idx2];
                gradient_val_seq_tmp(r_loc) = gradient_val_seq_tmp(r_loc) + radial_gradient_val_seq_tmp[r_idx2];
            }
        }
    }
    
    (*gradient_val_seq_pt) = gradient_val_seq_tmp;

}



void PGD::ComputeRadialGradient() {

    vector<VectorXd> gradient_val_seq;
    for (int i=0; i<num_radial_thread_assign; i++) {
        VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
        gradient_val_seq.push_back(gradient_val_seq_tmp);
    }
    
    if (num_radial_thread_assign>1) {
        thread *multi_thread = new thread[num_radial_thread_assign-1];
        for (int i=0; i<num_radial_thread_assign-1; i++) {
            vector<int> radial_partition_idx_tmp = radial_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeRadialGradientMuti, this, radial_partition_idx_tmp, &gradient_val_seq[i]);
        }
        
        ComputeRadialGradientMuti( radial_partition_idx[num_radial_thread_assign-1], &gradient_val_seq[num_radial_thread_assign-1]);
        
        for (int i=0; i<num_radial_thread_assign-1; i++) {
            multi_thread[i].join();
        }
        
        delete [] multi_thread;
        
    } else {
        ComputeRadialGradientMuti( radial_partition_idx[num_radial_thread_assign-1], &gradient_val_seq[num_radial_thread_assign-1]);
    }
    
    radial_gradient_val = VectorXd::Zero(num_loc);
    for (int i=0; i<num_radial_thread_assign; i++) {
        radial_gradient_val += gradient_val_seq[i];
    }

}


void PGD::ComputeProjGradientMuti( vector<DtPair> proj_partition_block, VectorXd* gradient_val_seq_pt ) {

    VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
    for (int i=0; i<proj_partition_block.size(); i++) {
        DtPair proj_idx=proj_partition_block[i];

        if (voxel_proj_idx.find(proj_idx)!=voxel_proj_idx.end()) {
            
        
            double proj_gradient_mut_tmp = ( est_proj_feat_all(proj_feat_all_loc[proj_idx]) -  proj_feat_all(proj_feat_all_loc[proj_idx]) );
            
            vector<int> voxel_group = voxel_proj_idx[proj_idx];
            int num_voxel_group = voxel_group.size();
            for (int i=0; i<num_voxel_group; i++) {
                gradient_val_seq_tmp(voxel_group[i]) = gradient_val_seq_tmp(voxel_group[i]) + proj_gradient_mut_tmp;
            }

        }
    }

    (*gradient_val_seq_pt) = gradient_val_seq_tmp;

}



void PGD::ComputeProjGradient() {

    vector<VectorXd> gradient_val_seq;
    for (int i=0; i<num_proj_thread_assign; i++) {
        VectorXd gradient_val_seq_tmp = VectorXd::Zero(num_loc);
        gradient_val_seq.push_back(gradient_val_seq_tmp);
    }

    if (num_proj_thread_assign>1) {
        thread *multi_thread = new thread[num_proj_thread_assign-1];
        for (int i=0; i<num_proj_thread_assign-1; i++) {
            vector<DtPair> proj_partition_idx_tmp = proj_partition_idx[i];
            multi_thread[i] = thread(&PGD::ComputeProjGradientMuti, this, proj_partition_idx_tmp, &gradient_val_seq[i]);
        }

        ComputeProjGradientMuti( proj_partition_idx[num_proj_thread_assign-1], &gradient_val_seq[num_proj_thread_assign-1]);

        for (int i=0; i<num_proj_thread_assign-1; i++) {
            multi_thread[i].join();
        }

        delete [] multi_thread;

    } else {
        ComputeProjGradientMuti( proj_partition_idx[num_proj_thread_assign-1], &gradient_val_seq[num_proj_thread_assign-1]);
    }

    proj_gradient_val = VectorXd::Zero(num_loc);
    for (int i=0; i<num_proj_thread_assign; i++) {
        proj_gradient_val += gradient_val_seq[i];
    }

}


PGD::~PGD() {}
