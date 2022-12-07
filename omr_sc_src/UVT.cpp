#include<chrono>
#include "UVT.h"
#include<zbessel.hh>

UVT::UVT() {

}

void UVT::SetInitFile(char* init_file) {
    density_init_file = init_file;
}

void UVT::SetData(DataReader* data_reader) {

    // add radial feat
    vector<vector<double>> raw_radial_feat = data_reader->getRadialData();
    for (int i=0; i<num_radial; i++) {
        int radial_idx = (int)raw_radial_feat[i][0];
        double raw_radial_feat_val = raw_radial_feat[i][1];
        
        radial_feat_idx.push_back(radial_idx);
        radial_feat[radial_idx] = raw_radial_feat_val;
    }

    // add radial map
    vector<vector<double>> raw_radial_map = data_reader->getRadialMapData();
    // the radial_idx must start from 0 consequtively
    for (int i=0; i<num_radial; i++) {
        int radial_idx = (int)raw_radial_map[i][0];
        double raw_radial_map_val = raw_radial_map[i][1];

        radial_map[radial_idx] = raw_radial_map_val;
    }
    num_radial_cor_feat = radial_map.size();    // for computing the 

    // add projection data
    vector<vector<double>> raw_proj_feat = data_reader->getProjData();
    for (int i=0; i<num_proj; i++) {

        int x_loc = (int)raw_proj_feat[i][0];
        int y_loc = (int)raw_proj_feat[i][1];
    
        DtPair xy_loc = {x_loc, y_loc};
        proj_feat_idx.push_back(xy_loc);
        proj_feat[xy_loc] = raw_proj_feat[i][2];
    }


    // add fourier radial map
    vector<vector<double>> raw_fourier_radial_map = data_reader->getFourierRadialMapData();
    for (int i=0; i<num_fourier_radial; i++) {
        int fourier_radial_idx = (int)raw_fourier_radial_map[i][0];
        double raw_fourier_radial_map_val = raw_fourier_radial_map[i][1];

        fourier_radial_map[fourier_radial_idx] = raw_fourier_radial_map_val;
    }

    // add fourier correlation feat
    vector<vector<double>> raw_fourier_cor_feat = data_reader->getFourierCorData();

    for (int i=0; i<num_fourier_cor_seq; i++) {
        vector<double> raw_fourier_cor_feat_val = raw_fourier_cor_feat[i];

        int l_idx = (int) raw_fourier_cor_feat_val[0];
        int m_idx = (int) raw_fourier_cor_feat_val[1];
        // double check if the fourier correlation features are symmetric
        raw_fourier_cor_feat_val.erase(raw_fourier_cor_feat_val.begin(), raw_fourier_cor_feat_val.begin()+2);
        DtPair fourier_cor_feat_idx = {l_idx, m_idx};
        
        num_fourier_radial_cor_feat = raw_fourier_cor_feat_val.size();
        
        VectorXd raw_fourier_cor_feat_seq(raw_fourier_cor_feat_val.size());
        for (int j=0; j<raw_fourier_cor_feat_val.size(); j++) {
            raw_fourier_cor_feat_seq(j) = raw_fourier_cor_feat_val[j];
        }
        fourier_cor_feat[fourier_cor_feat_idx] = raw_fourier_cor_feat_seq;
        
        // fourier correlation feature pair sequence
        // fourier_cor_feat_lm_unique should be the same for both fourier correlation and spatial correlation
        // but we will do it anyway, error 
        if (fourier_cor_feat_lm_unique.find(fourier_cor_feat_idx)==fourier_cor_feat_lm_unique.end()) {
            fourier_cor_feat_lm_unique[fourier_cor_feat_idx]=1;
        }

    }
    
    if (num_fourier_radial_cor_feat!=fourier_radial_map.size()) {
        cout<<"num_fourier_radial_cor_feat is not consistent with fourier_radial_map size! Error!"<<endl;
        cout<<num_fourier_radial_cor_feat<<" "<<fourier_radial_map.size()<<endl;
        abort();
    }

    // add fourier l_idx data
    // the order of the l_idx data should be consistent with cor_feat
    // the radius DOES NOT start with radius = 0
    vector<vector<double>> raw_fourier_radial_l_idx = data_reader->getFourierLIdxData();
    for (int i=0; i<num_fourier_l_idx_seq; i++) {
        vector<double> raw_fourier_radial_l_idx_val = raw_fourier_radial_l_idx[i];
        int l_idx = (int) raw_fourier_radial_l_idx_val[0];
        int m_idx = (int) raw_fourier_radial_l_idx_val[1];

        raw_fourier_radial_l_idx_val.erase(raw_fourier_radial_l_idx_val.begin(), raw_fourier_radial_l_idx_val.begin()+2);
        DtPair fourier_cor_feat_idx = {l_idx, m_idx};
        int raw_fourier_radial_l_idx_num = raw_fourier_radial_l_idx_val.size();
        
        vector<int> raw_fourier_radial_l_idx_seq;
        for (int j=0; j<raw_fourier_radial_l_idx_num; j++) {
            int raw_fourier_radial_l_idx_int_val = (int) raw_fourier_radial_l_idx_val[j];
            raw_fourier_radial_l_idx_seq.push_back(raw_fourier_radial_l_idx_int_val);
        }
        fourier_radial_l_idx_seq[fourier_cor_feat_idx] = raw_fourier_radial_l_idx_seq;

    }
    
    
    // add radial lgwt data
    radial_lgwt_range = VectorXd::Zero(num_radial);
    radial_lgwt_weights = VectorXd::Zero(num_radial);
    vector<vector<double>> raw_radial_lgwt_data = data_reader->getRadialLGWTWTSData();
    for (int i=0; i<num_radial; i++) {
        radial_lgwt_range(i) = raw_radial_lgwt_data[i][0];
        radial_lgwt_weights(i) = raw_radial_lgwt_data[i][1];
    }
    
    //// set cor_feat_mat_all
    //for (int l=0; l<=L_max; l++) {
    //    MatrixXd cor_feat_mat_all_tmp = MatrixXd::Zero(num_radial_cor_feat,2*l+1);
    //    cor_feat_mat_all[l] = cor_feat_mat_all_tmp;
    //}    
    
    // set est_cor_feat_mat_all and also calculate num_cor_feat
    num_cor_feat = 0;
    for (int l=0; l<=L_max; l++) {
        MatrixXd cor_feat_mat_all_tmp = MatrixXd::Zero(num_radial_cor_feat,2*l+1);
        est_cor_feat_mat_all[l] = cor_feat_mat_all_tmp;
        num_cor_feat += num_radial_cor_feat*(2*l+1);
    }    


    // set fourier_cor_feat_mat_all, note this is real matrix from decomposition
    for (int l=0; l<=L_max; l++) {
        MatrixXd fourier_cor_feat_mat_all_tmp = MatrixXd::Zero(num_fourier_radial_cor_feat,2*l+1);
        for (int m=-l; m<=l; m++) {
            DtPair lm_pair = {l,m};
            fourier_cor_feat_mat_all_tmp.col(m+l) = fourier_cor_feat[lm_pair];
        }
        fourier_cor_feat_mat_all[l] = fourier_cor_feat_mat_all_tmp;
    }
    
    // set fourier_est_cor_feat_mat_all, note this is complex
    for (int l=0; l<=L_max; l++) {
        MatrixXd fourier_cor_feat_mat_all_tmp = MatrixXd::Zero(num_fourier_radial_cor_feat,2*l+1);
        est_fourier_cor_feat_mat_all[l] = fourier_cor_feat_mat_all_tmp;
    }
    
    // set fourier orthogonal matrices
    for (int l=0; l<=L_max; l++) {
        MatrixXd fourier_orth_mat_all_tmp = MatrixXd::Zero(2*l+1,2*l+1);
        for (int i=0; i<2*l+1; i++) {
            fourier_orth_mat_all_tmp(i,i)=1;
        }
        fourier_orth_mat_all[l] = fourier_orth_mat_all_tmp;
    }

    int cor_feat_all_loc_val = 0;   // the location of the correlation feature in cor_feat_all
    // The l,m pair is the same as the fourier_cor_feat
    for (map<DtPair, int>::iterator ite=fourier_cor_feat_lm_unique.begin(); ite!=fourier_cor_feat_lm_unique.end(); ite++) {
    
        DtPair cor_feat_idx = ite->first;
        for (int i=0; i<radial_map.size(); i++) {
            DtTriplet cor_feat_all_loc_key = {cor_feat_idx[0], cor_feat_idx[1], i};
            cor_feat_all_loc[cor_feat_all_loc_key] = cor_feat_all_loc_val;
            cor_feat_all_loc_val += 1;
        }
    }    


    
    // set fourier_cor_feat_all 
    // note that this is different from before, we are using l m here
    num_fourier_cor_feat = 0;
    for (map<DtPair, int>::iterator ite=fourier_cor_feat_lm_unique.begin(); ite!=fourier_cor_feat_lm_unique.end(); ite++) {
        num_fourier_cor_feat = num_fourier_cor_feat + fourier_cor_feat[ite->first].size();
    }
    
    fourier_cor_feat_all = VectorXd::Zero(num_fourier_cor_feat);
    int fourier_cor_feat_all_idx = 0;
    int fourier_cor_feat_all_loc_val = 0;   // the location of the fourier correlation feature in cor_feat_all
    for (map<DtPair, int>::iterator ite=fourier_cor_feat_lm_unique.begin(); ite!=fourier_cor_feat_lm_unique.end(); ite++) {
        VectorXd fourier_cor_feat_tmp = fourier_cor_feat[ite->first];
        vector<int> fourier_radial_l_idx_seq_tmp = fourier_radial_l_idx_seq[ite->first];
        fourier_cor_feat_all.segment(fourier_cor_feat_all_idx, fourier_cor_feat_tmp.size()) = fourier_cor_feat_tmp;
        fourier_cor_feat_all_idx += fourier_cor_feat_tmp.size();
    
        // check to make sure the lengths of cor_feat_tmp and radial_l_idx_seq_tmp are the same
        if (fourier_cor_feat_tmp.size()!=fourier_radial_l_idx_seq_tmp.size()) {
            cout<<"The lengths of fourier_cor_feat and fourier_l_idx are not the same!"<<endl;
            abort();
        }
    
        DtPair fourier_cor_feat_idx = ite->first;
        for (int i=0; i<fourier_radial_l_idx_seq_tmp.size(); i++) {
            DtTriplet fourier_cor_feat_all_loc_key = {fourier_cor_feat_idx[0], fourier_cor_feat_idx[1], fourier_radial_l_idx_seq_tmp[i]};
            fourier_cor_feat_all_loc[fourier_cor_feat_all_loc_key] = fourier_cor_feat_all_loc_val;
            fourier_cor_feat_all_loc_val += 1;
        }
    }

    
    // set est_fourier_cor_feat_all
    est_fourier_cor_feat_all = VectorXd::Zero(num_fourier_cor_feat);
    
    
    // set transformed fourier_cor_feat_all using fourier_orth_mat_all
    // compute fourier_cor_feat_all_transform
    map<int, MatrixXd> fourier_cor_feat_mat_all_transform;
    for (int l=0; l<=L_max; l++) {
        fourier_cor_feat_mat_all_transform[l] = fourier_cor_feat_mat_all[l] * fourier_orth_mat_all[l];
    }
    
    fourier_cor_feat_all_transform = VectorXd::Zero(num_fourier_cor_feat);
    for (int l=0; l<=L_max; l++) {
        MatrixXd fourier_cor_feat_mat_all_transform_tmp = fourier_cor_feat_mat_all_transform[l];
        for (int m=-l; m<=l; m++) {
            for (int r=0; r<fourier_radial_map.size(); r++) {
                DtTriplet fourier_cor_feat_all_loc_key = {l,m,r};
                int fourier_cor_feat_all_loc_val = fourier_cor_feat_all_loc[fourier_cor_feat_all_loc_key];
                fourier_cor_feat_all_transform(fourier_cor_feat_all_loc_val) = fourier_cor_feat_mat_all_transform_tmp(r,m+l);
            }
        }
    }

    // set radial_feat
    num_radial_feat = radial_feat.size();
    radial_feat_all = VectorXd::Zero(num_radial_feat);
    for (int i=0; i<radial_feat_idx.size(); i++) {
        radial_feat_all(i) = radial_feat[radial_feat_idx[i]];
        radial_feat_all_loc[radial_feat_idx[i]] = i;
    }

    // set proj_feat
    num_proj_feat = proj_feat.size();
    proj_feat_all = VectorXd::Zero(num_proj_feat);
    max_proj_feat = 0;
    for (int i=0; i<proj_feat_idx.size(); i++) {
        proj_feat_all(i) = proj_feat[proj_feat_idx[i]];
        proj_feat_all_loc[proj_feat_idx[i]] = i;
        if (proj_feat_all(i) > max_proj_feat) {
            max_proj_feat = proj_feat_all(i);
        }
    }

    cout<<"Max proj feat: "<<max_proj_feat<<endl;
    
    // set est_radial_feat_all
    est_radial_feat_all = VectorXd::Zero(num_radial_feat);

    // set est_proj_feat_all
    est_proj_feat_all = VectorXd::Zero(num_proj_feat);
    
    cout<<"Setting data finished!"<<endl;
    
}



void UVT::Census3DSpace() {

    // the location of the centroid is (0,0,0)
    // starting from (-N_half,-N_half,-N_half), we check every location

    int idx_pos = 0;    //index of the voxel locations, starts from 0
    
    N_half = (N-1)/2;
    for (int i=-N_half; i<=N_half; i++) {
        for (int j=-N_half; j<=N_half; j++) {
            for (int k=-N_half; k<=N_half; k++) {
            
                double radial_db = sqrt(i*i*1.0+j*j*1.0+k*k*1.0);
                // find the radial_idx

                // just incase we want to shrink the space even further
                if (radial_db>r_max*1.0) {
                    continue;
                }

                // prune the locations where proj_feat are not recorded (i.e. zero)
                if (proj_feat.find({i,j})==proj_feat.end()) {
                    continue;
                }

                
                DtTriplet loc_tmp = {0,0,0};
                loc_tmp[0]=i; loc_tmp[1]=j; loc_tmp[2]=k;

                if (voxel_proj_idx.find({i,j})==voxel_proj_idx.end()) {
                    vector<int> voxel_proj_idx_tmp;
                    voxel_proj_idx_tmp.push_back(idx_pos);
                    voxel_proj_idx[{i,j}] = voxel_proj_idx_tmp;
                } else {
                    voxel_proj_idx[{i,j}].push_back(idx_pos);
                }

                
                voxel voxel_map_tmp;
                voxel_map_tmp.coordinate = loc_tmp;
                voxel_map_tmp.radius = radial_db;                    

                voxel_map[idx_pos] = voxel_map_tmp;
                voxel_map_reverse[loc_tmp] = idx_pos;

                // find the comass_idx
                if (radial_db==0) {
                    comass_idx = idx_pos;
                }
                
                idx_pos++;
            }
        }
    }
    
    // the total number of voxels
    num_loc = idx_pos;
    cout<<"Total number candidate locations: "<<num_loc<<endl;
    
    // save the location values
    loc_val = MatrixXd::Zero(num_loc, 3);
    loc_sph_val = MatrixXd::Zero(num_loc, 3);
    radius_val = VectorXd::Zero(num_loc);
    for (int i=0; i<num_loc; i++) {
        DtTriplet coordinate_tmp = voxel_map[i].coordinate;

        loc_val(i,0) = (double)coordinate_tmp[0];
        loc_val(i,1) = (double)coordinate_tmp[1];
        loc_val(i,2) = (double)coordinate_tmp[2];

        loc_val_x.push_back(coordinate_tmp[0]);
        loc_val_y.push_back(coordinate_tmp[1]);
        loc_val_z.push_back(coordinate_tmp[2]);

        radius_val(i) = voxel_map[i].radius;
    
        if (i==comass_idx) {
            loc_sph_val(i,0) = 0;
            loc_sph_val(i,1) = 0;
            loc_sph_val(i,2) = radius_val(i);
        } else {
            loc_sph_val(i,0) = atan2(loc_val(i,1),loc_val(i,0));
            loc_sph_val(i,1) = acos(loc_val(i,2)/radius_val(i));    // between 0 and pi
            loc_sph_val(i,2) = radius_val(i);
        }

        // this is the coordinates of the first location in the test_script
        //if (coordinate_tmp[0]==-27) {
        //    if (coordinate_tmp[1]==4) {
        //        if (coordinate_tmp[2]==-9) {
        //            cout<<i<<" "<<loc_sph_val(i,0)<<" "<<loc_sph_val(i,1)<<" "<<loc_sph_val(i,2)<<endl;
        //        }
        //    }
        //}
    }


    map<DtPair, VectorXd> sh_val;               // the (l,m) pair, and the spherical harmonic values according to voxel indices
    for (int l=0; l<=L_max; l++) {
        for (int m=-l; m<=l; m++) {
            DtPair lm_pair = {l, m};
            VectorXd sh_val_tmp = VectorXd::Zero(num_loc);
            for (int i=0; i<num_loc; i++) {
                if (i==comass_idx) {
                    if (l==0) {
                        sh_val_tmp(i) = sqrt(1/(4*M_PI));
                    } else {
                        sh_val_tmp(i) = 0;
                    }
                } else {
                    if (m==0) {
                        sh_val_tmp(i) = sph_legendre(l,m,loc_sph_val(i,1));
                    } else if (m>0) {
                        sh_val_tmp(i) = sqrt(2)*sph_legendre(l,m,loc_sph_val(i,1))*cos(m*loc_sph_val(i,0));
                    } else {
                        // there are different definitions when m<0 between wikipedia and nist handbook, which one is correct??? the following one is based on wikipedia
                        sh_val_tmp(i) = sqrt(2)*sph_legendre(l,abs(m),loc_sph_val(i,1))*sin(abs(m)*loc_sph_val(i,0));
                    }
                }
            }
            
            //if (abs(m) % 2 == 1) {
            //    sh_val_tmp = (-1)*sh_val_tmp;
            //}
            //cout<<l<<" "<<m<<" "<<sh_val_tmp(32)<<endl;
            sh_val[lm_pair] = sh_val_tmp;
        }
    }

    
    for (int r_idx=0; r_idx<radial_map.size(); r_idx++) {

        double rval = radial_map[r_idx];
   
        //cout<<r<<endl;
        VectorXd radius_diff_val = abs(radius_val.array() - 1.0*rval);
        vector<int> voxel_radial_idx_tmp;
        double voxel_std_3times = 3.0*voxel_std;
        vector<double> radius_val_vect_tmp; // the ||mu||_2 that falls within the range of r
        // check every candidate voxel 
        //int spot_32 = 0;
        //int spot_32_idx = -1;
        for (int i=0; i<num_loc; i++) {
            if (radius_diff_val(i)<=voxel_std_3times) {
                voxel_radial_idx_tmp.push_back(i);
                radius_val_vect_tmp.push_back(radius_val(i));
                //if (i==32) {
                //    spot_32 = 1;
                //    spot_32_idx = radius_val_vect_tmp.size()-1;
                //}
            }
        }

        if (radius_val_vect_tmp.size()==0) {
            continue;
        }
        
        int radius_val_len_tmp = radius_val_vect_tmp.size();
        
        r_candidate[r_idx] = 1;
        voxel_radial_idx[r_idx] = voxel_radial_idx_tmp;
        
        VectorXd radius_val_tmp = VectorXd::Zero(radius_val_len_tmp);   // a VectorXd version of radius_val_vect_tmp
        for (int i=0; i<radius_val_len_tmp; i++) {
            radius_val_tmp(i) = radius_val_vect_tmp[i];
        }

        //// Note there is no 4.0*M_PI here!!!! just r*r
        
        if (rval==0) {
            // compute radial_beta_seq
            radial_beta_seq[r_idx] = 1/pow(sqrt(2*M_PI)*voxel_std,3) * (-0.5/voxel_var * (radius_val_tmp.array()).pow(2)).exp();
        } else {
            // compute radial_alpha_seq
            // does not multiply by r_val here!!!
            VectorXd mut_part_1 = (1/voxel_var) * (-0.5/voxel_var * (rval-radius_val_tmp.array()).pow(2)).exp() / (rval*radius_val_tmp.array()).sqrt(); // what if radius_val_tmp==0???

            VectorXd bessel_input = (rval/voxel_var)*radius_val_tmp;
            vector<double> bessel_input_vect_real(radius_val_len_tmp);
            vector<double> bessel_input_vect_imag(radius_val_len_tmp);
            for (int i=0; i<radius_val_len_tmp; i++) {
                bessel_input_vect_real[i] = bessel_input(i);
                bessel_input_vect_imag[i] = 0;
            }
            
            MatrixXd radial_alpha_seq_tmp = MatrixXd::Zero(L_max+1,radius_val_len_tmp);
            
            VectorXd mut_part_2 = rval*rval *  1/(voxel_var) * (-0.5/voxel_var * (rval-radius_val_tmp.array()).pow(2)).exp() / (rval*radius_val_tmp.array()).sqrt();    // what if radius_val_tmp==0???
            VectorXd radial_beta_seq_tmp = VectorXd::Zero(radius_val_len_tmp);
            
            double fnu = 0.5;
            int kode = 2;
            int l_seq_num = L_max+1;
            double zbesi_val_real[l_seq_num] {}, zbesi_val_imag[l_seq_num] {};
            int nz {}, ierr {};
                
            for (int i=0; i<radius_val_len_tmp; i++) {
                ierr = zbessel::zbesi(bessel_input_vect_real[i], bessel_input_vect_imag[i], fnu, kode, l_seq_num, zbesi_val_real, zbesi_val_imag, &nz);
                if (ierr!=0) {
                    cout<<"Computation of scaled bessel function error: "<<ierr<<endl;
                    abort();
                }

                if (radius_val_tmp(i)==0) {
                    for (int l=0; l<=L_max; l++) {
                        if (l==0) {
                            radial_alpha_seq_tmp(l,i) = rval * 1/pow(sqrt(2*M_PI)*voxel_std,3) * exp(-0.5/voxel_var*rval*rval) * sqrt(4*M_PI);
                        } else {
                            radial_alpha_seq_tmp(l,i) = 0;
                        }
                    }
                    radial_beta_seq_tmp(i) = 1/pow(sqrt(2*M_PI)*voxel_std,3) * exp(-0.5/voxel_var*rval*rval) * (4*M_PI*rval*rval);
                } else {
                    for (int l=0; l<=L_max; l++) {
                        radial_alpha_seq_tmp(l,i) = mut_part_1(i) * zbesi_val_real[l];
                    }
                    radial_beta_seq_tmp(i) = mut_part_2(i) * zbesi_val_real[0];
                }
            }
            
            for (int l=0; l<=L_max; l++) {
                for (int m=-l; m<=l; m++) {
                    VectorXd sh_val_tmp = sh_val[{l,m}];
                    VectorXd radial_sh_seq_tmp = VectorXd::Zero(radius_val_len_tmp);
                    for (int i=0; i<radius_val_len_tmp; i++) {
                        radial_sh_seq_tmp(i) = sh_val_tmp(voxel_radial_idx_tmp[i]);
                    }
                    // it is a simple transpose here
                    VectorXd radial_alpha_sh_seq_tmp = radial_alpha_seq_tmp.row(l).transpose().array() * radial_sh_seq_tmp.array();
                    radial_alpha_sh_seq[{r_idx, l, m}] = radial_alpha_sh_seq_tmp;

                }
            }
            
            // compute radial_beta_seq
            radial_beta_seq[r_idx] = radial_beta_seq_tmp;

            //if (spot_32==1) {
            //    cout<<"r: "<<r<<" beta_seq: "<<radial_beta_seq_tmp(spot_32_idx)<<endl;
            //}
            //
            //if (spot_32==1) {
            //    for (int l=0; l<=L_max; l++) {
            //        cout<<"r: "<<r<<" "<<l<<" alpha_seq: "<<radial_alpha_seq_tmp(l,spot_32_idx)<<endl;
            //    }
            //}


        }
        
    }
   
    
    // coompute radial_lgwt_mat
    for (int l=0; l<=L_max; l++) {
        MatrixXd radial_lgwt_mat_tmp = MatrixXd::Zero(num_radial_cor_feat,num_fourier_radial_cor_feat);
        for (int i=0; i<num_radial_cor_feat; i++) {
            double radial_lgwt_pre_mut = radial_lgwt_weights(i)*pow(radial_lgwt_range(i),2);
            for (int j=0; j<num_fourier_radial_cor_feat; j++) {
                radial_lgwt_mat_tmp(i,j) = radial_lgwt_pre_mut*sph_bessel(l,radial_lgwt_range(i)*fourier_radial_map[j]);
            }
        }
        radial_lgwt_mat_tmp = (4.0*M_PI)*radial_lgwt_mat_tmp;
        //radial_lgwt_mat[l] = pow(1i,l)*radial_lgwt_mat_tmp;
        radial_lgwt_mat[l] = radial_lgwt_mat_tmp; // remove pow(1i,l) since the coefficients are purely real or imaginary
    }
    
    
    int num_lm_total = (1+L_max)*L_max + L_max + 1;
    num_lm_thread_assign = num_thread > num_lm_total ? num_lm_total : num_thread;
    
    vector<int> lm_partition_sz(num_lm_thread_assign,0);
    for (int l=0; l<=L_max; l++) {
        for (int m=-l; m<=l; m++) {
            int idx_min = distance( lm_partition_sz.begin(), min_element(lm_partition_sz.begin(), lm_partition_sz.end()) );
            DtPair lm_idx_tmp = {l,m};
            if (lm_partition_idx.find(idx_min)==lm_partition_idx.end()) {
                vector<DtPair> lm_partition_idx_tmp;
                lm_partition_idx_tmp.push_back(lm_idx_tmp);
                lm_partition_idx[idx_min] = lm_partition_idx_tmp;
            } else {
                lm_partition_idx[idx_min].push_back(lm_idx_tmp);
            }
            lm_partition_sz[idx_min] += 1;
        }
    }
    
    num_radial_thread_assign = num_thread > num_radial_feat ? num_radial_feat : num_thread;
    vector<int> radial_partition_sz(num_radial_thread_assign,0);
    for (int i=0; i<radial_feat_idx.size(); i++) {
        int idx_min = distance( radial_partition_sz.begin(), min_element(radial_partition_sz.begin(), radial_partition_sz.end()) );
        if (radial_partition_idx.find(idx_min)==radial_partition_idx.end()) {
            vector<int> radial_partition_idx_tmp;   // push radii indices
            radial_partition_idx_tmp.push_back(radial_feat_idx[i]);
            radial_partition_idx[idx_min] = radial_partition_idx_tmp;
        } else {
            radial_partition_idx[idx_min].push_back(radial_feat_idx[i]);
        }
        radial_partition_sz[idx_min] += 1;
    }

    num_proj_thread_assign = num_thread > num_proj_feat ? num_proj_feat : num_thread;
    vector<int> proj_partition_sz(num_proj_thread_assign,0);
    for (int i=0; i<proj_feat_idx.size(); i++) {
        int idx_min = distance( proj_partition_sz.begin(), min_element(proj_partition_sz.begin(), proj_partition_sz.end()) );
        if (proj_partition_idx.find(idx_min)==proj_partition_idx.end()) {
            vector<DtPair> proj_partition_idx_tmp;   // push radii indices
            proj_partition_idx_tmp.push_back(proj_feat_idx[i]);
            proj_partition_idx[idx_min] = proj_partition_idx_tmp;
        } else {
            proj_partition_idx[idx_min].push_back(proj_feat_idx[i]);
        }
        proj_partition_sz[idx_min] += 1;
    }    

    cout<<"Census3DSpace finished"<<endl;
    
}



void UVT::Initialization() {

    if (init_type==1) { // random initialization
    
        // density_val ordered according to the voxel indices
        density_val = VectorXd::Zero(num_loc);
        unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int i=0; i<num_loc; i++) {
            density_val(i) = abs(distribution(generator));
        }
        
    } else if (init_type==2) {  // read from initialization file
        
        density_val = VectorXd::Zero(num_loc);
        
        ifstream read_data;
        istringstream istr;
        string str; 
        
        read_data.open(density_init_file);

        if (!read_data) {
            cout<<"Initialization does not exist!!!"<<endl;
            abort();
        }

        vector<vector<double>> raw_density_init;
        int line_num=0;

        while(getline(read_data, str)) {
            istr.str(str);
            double tmp;
            vector<double> tmpvec;
            while(istr>>tmp) {
                tmpvec.push_back(tmp);
            }
            raw_density_init.push_back(tmpvec);

            tmpvec.clear();
            istr.clear();
            line_num++;
        }   

        read_data.close();
        
        for (int i=0; i<raw_density_init.size(); i++) {
            int loc_x = (int)raw_density_init[i][0];
            int loc_y = (int)raw_density_init[i][1];
            int loc_z = (int)raw_density_init[i][2];
            double raw_density_init_tmp = raw_density_init[i][3];
            DtTriplet loc_tmp = {loc_x, loc_y, loc_z};
            if (voxel_map_reverse.find(loc_tmp)!=voxel_map_reverse.end()) {
                density_val(voxel_map_reverse[loc_tmp]) = raw_density_init_tmp;
            }
        }
    
    } else{ // Error Message

        cout<<"Unknown initialization type."<<endl;
        abort();
    }
    //cout<<"density: "<<density_val(32)<<endl;
    
    // enforce projection constraint and make the three types of densities consistent with one another
    for (map<DtPair, vector<int>>::iterator ite=voxel_proj_idx.begin(); ite!=voxel_proj_idx.end(); ite++) {
        DtPair proj_tmp = ite->first;
        vector<int> voxel_group = ite->second;
        int num_voxel_group = voxel_group.size();

        VectorXd density_val_tmp = VectorXd::Zero(num_voxel_group);
        for (int i=0; i<num_voxel_group; i++) {
            density_val_tmp(i) = density_val(voxel_group[i]);
        }

        // enforce projection constraint
        proj_density_val[proj_tmp] = density_val_tmp;

    }


    // projection of density_val
    density_val = ProjectOntoCvxSet(density_val/max_proj_feat, density_all/max_proj_feat);
    density_val = density_val*max_proj_feat;
    
    // since we are relying on radial_density_val to compute the gradients, so we are going to update the radial_desntity_val. 
    // we can actually do it iteratively until the radial and projection constraints are both satisfied, 
    // here we are going to just do it once
    
    for (map<int, vector<int>>::iterator ite=voxel_radial_idx.begin(); ite!=voxel_radial_idx.end(); ite++) {
        int radial_tmp = ite->first;
        vector<int> voxel_group = ite->second;
        int num_voxel_group = voxel_group.size();
        
        VectorXd density_val_tmp = VectorXd::Zero(num_voxel_group);
        for (int i=0; i<num_voxel_group; i++) {
            density_val_tmp(i) = density_val(voxel_group[i]);
        }
        
        // enforce radial constraint
        radial_density_val[radial_tmp] = density_val_tmp;
    }
    //cout<<"density: "<<density_val(32)<<endl;

    // Compute objective function
    cor_obj_val = ComputeCorObjFun();
    radial_obj_val = ComputeRadialObjFun();
    proj_obj_val = ComputeProjObjFun();
    obj_val = cor_obj_val + radial_weight * radial_obj_val + proj_weight * proj_obj_val;
    cout<<"Objective function: "<<cor_obj_val<<"\t"<<radial_obj_val<<"\t"<<proj_obj_val<<"\t"<<obj_val<<endl;
    
    cout<<"Initialization finished."<<endl;
    

}


VectorXd UVT::ProjectOntoCvxSet(VectorXd smp_vec, double num_smp_proj) {
    
    // make sure all the entries of v are above the convex set
    // Note that the vector index starts from 0
    double offset = -smp_vec.array().minCoeff() + num_smp_proj;
    smp_vec = smp_vec.array() + offset;
    
    VectorXd smp_sort_vec = smp_vec;
    sort(smp_sort_vec.data(), smp_sort_vec.data()+smp_sort_vec.size(), greater<double>());
    int smp_vec_len=smp_sort_vec.size();
    
    double theta = 0;
    int check_status = 0;
    for (int r=1; r<=num_smp_proj; r++) {
        VectorXd smp_sort_vec_new = smp_sort_vec.segment(r-1, smp_vec_len-r+1);
        VectorXd smp_sort_vec_new_cumsum = smp_sort_vec_new;
        for (int j=1; j<smp_sort_vec_new_cumsum.size(); j++) {
            smp_sort_vec_new_cumsum(j) += smp_sort_vec_new_cumsum(j-1);
        }
        VectorXd smp_sort_vec_new_thr = smp_sort_vec_new_cumsum;
        for (int j=0; j<smp_sort_vec_new_thr.size(); j++) {
            smp_sort_vec_new_thr(j) = smp_sort_vec_new(j) - (smp_sort_vec_new_thr(j)-(num_smp_proj-r+1))/(j+1);
        }
        int rho_new = -1;
        for (int j=smp_sort_vec_new_thr.size()-1; j>=0; j--) {
            if (smp_sort_vec_new_thr(j)>0) {
                rho_new = j+1;  // Note that the vector index starts from 0, we need to add 1 here
                break;
            }
        }
        if (rho_new==-1) {break;}
        
        // check that rho = rho_new+r-1 is larger than num_smp_proj
        if (rho_new+r-1<=num_smp_proj) {continue;}
        // Compute the threshold
        double theta_new = (smp_sort_vec_new_cumsum(rho_new-1)-(num_smp_proj-r+1))/rho_new;
        
        int break_marker = 0;
        double w_r = smp_sort_vec(r-1)-theta_new;
        if ( (w_r>0) && (w_r<1) ) {
            if (r==1) {
                break_marker = 1;
                check_status = 1;
            } else {
                double w_rm1 = smp_sort_vec(r-2)-theta_new ;
                if (w_rm1>=1) {
                    break_marker=1;
                    check_status=1;
                }
            }
        } else {
            continue;
        }
        
        theta = theta_new;
        if (break_marker==1) {
            break;
        }
    }
    
    VectorXd smp_vec_proj = smp_vec;
    
    if (check_status==1) {
        smp_vec_proj = smp_vec.array() - theta;
        for (int i=0; i<smp_vec_len; i++) {
            if (smp_vec_proj(i)<0) {
                smp_vec_proj(i)=0;
            }
            if (smp_vec_proj(i)>1) {
                smp_vec_proj(i)=1;
            }
        }
    } else {
        // set the top N entries to 1 and the rest to 0
        double thd_tmp = smp_sort_vec(num_smp_proj-1);
        for (int i=0; i<smp_vec_len; i++) {
            if (smp_vec_proj(i)>=thd_tmp) {
                smp_vec_proj(i)=1;
            } else {
                smp_vec_proj(i)=0;
            }
        }
    }
    
    return smp_vec_proj;
    
}

void UVT::OrthoMatrixSync() {

    // when l=0; the orthognal matrix becomes 1
    // we only update orthogonal matrices for l>=1
    
    // we need a better way to coordinate the computation of objective function to avoid duplicate computations between orthomatrixsync and gradient descent
    // Compute objective function to get the estimated correlation features
    //cor_obj_val = ComputeCorObjFun();
    
    
    

    
    // multiply conjungate transpose(est_cor_feat_mat_all) with cor_feat_mat_all
    // is it transpose or conjungate transpose?    
    map<int, MatrixXd> mut_mat_all;
    for (int l=0; l<=L_max; l++) {
        mut_mat_all[l] = est_fourier_cor_feat_mat_all[l].adjoint() * fourier_cor_feat_mat_all[l];
    }
    
    // perform svd to compute the orthognal matrix
    for (int l=0; l<=L_max; l++) {
        MatrixXd mut_mat_all_tmp = mut_mat_all[l];
        JacobiSVD<MatrixXd> svd_tmp(mut_mat_all_tmp, ComputeFullU | ComputeFullV);
        fourier_orth_mat_all[l] = svd_tmp.matrixV() * svd_tmp.matrixU().adjoint();
    }

    // compute fourier_cor_feat_all_transform
    map<int, MatrixXd> fourier_cor_feat_mat_all_transform;
    for (int l=0; l<=L_max; l++) {
        fourier_cor_feat_mat_all_transform[l] = fourier_cor_feat_mat_all[l] * fourier_orth_mat_all[l];
    }
    
    fourier_cor_feat_all_transform = VectorXd::Zero(num_fourier_cor_feat);
    for (int l=0; l<=L_max; l++) {
        MatrixXd fourier_cor_feat_mat_all_transform_tmp = fourier_cor_feat_mat_all_transform[l];
        for (int m=-l; m<=l; m++) {
            for (int r_idx=0; r_idx<fourier_radial_map.size(); r_idx++) {
                DtTriplet fourier_cor_feat_all_loc_key = {l,m,r_idx};
                int fourier_cor_feat_all_loc_val = fourier_cor_feat_all_loc[fourier_cor_feat_all_loc_key];
                fourier_cor_feat_all_transform(fourier_cor_feat_all_loc_val) = fourier_cor_feat_mat_all_transform_tmp(r_idx,m+l);
            }
        }
    }
    
    // Compute objective function
    cor_obj_val = ComputeCorObjFun();
    radial_obj_val = ComputeRadialObjFun();
    proj_obj_val = ComputeProjObjFun();
    obj_val = cor_obj_val + radial_weight * radial_obj_val + proj_weight * proj_obj_val;
    cout<<"OrthoMatrixSync Objective function: "<<cor_obj_val<<"\t"<<radial_obj_val<<"\t"<<proj_obj_val<<"\t"<<obj_val<<endl;
    
    
}

void UVT::GradientDescent() {

    step = step_ori;

    double obj_val_pre;
    VectorXd density_val_pre;
    map<int, VectorXd> radial_density_val_pre;
    map<DtPair, VectorXd> proj_density_val_pre;

    int ite;    // count the total iteration number
    
    int num_cov = 0;    // the number of times the method reaches convergence
    for (ite =1; ite<max_ite; ite++) {

        obj_val_pre = obj_val;
        density_val_pre = density_val;
        radial_density_val_pre = radial_density_val;
        proj_density_val_pre = proj_density_val;
      
        ComputeRadialGradient();
        ComputeProjGradient();
        ComputeCorGradient();
        gradient_val = cor_gradient_val + radial_weight * radial_gradient_val + proj_weight * proj_gradient_val;
        
        VectorXd density_val_bkt = VectorXd::Zero(density_val.size());
        map<int, VectorXd> radial_density_val_bkt;
        map<DtPair, VectorXd> proj_density_val_bkt;
        while (step>step_thd) {
        
            density_val_bkt = density_val - gradient_val*step;
            

            // projection of density_val
            density_val_bkt = ProjectOntoCvxSet(density_val_bkt/max_proj_feat, density_all/max_proj_feat);
            density_val_bkt = density_val_bkt*max_proj_feat;

            // make radial_density_val_bkt consistent
            // we use radial_density_val_bkt and density_val_bkt to compute the objective function
            // so these two must be consistent all the time!!!
            for (map<int, vector<int>>::iterator ite=voxel_radial_idx.begin(); ite!=voxel_radial_idx.end(); ite++) {
                int radial_tmp = ite->first;
                vector<int> voxel_group = ite->second;
                int num_voxel_group = voxel_group.size();
                
                VectorXd density_val_bkt_tmp = VectorXd::Zero(num_voxel_group);
                for (int i=0; i<num_voxel_group; i++) {
                    density_val_bkt_tmp(i) = density_val_bkt(voxel_group[i]);
                }

                radial_density_val_bkt[radial_tmp] = density_val_bkt_tmp;
            }

            // make proj_density_val_bkt consistent
            for (map<DtPair, vector<int>>::iterator ite=voxel_proj_idx.begin(); ite!=voxel_proj_idx.end(); ite++) {
                DtPair proj_tmp = ite->first;
                vector<int> voxel_group = ite->second;
                int num_voxel_group = voxel_group.size();

                VectorXd density_val_bkt_tmp = VectorXd::Zero(num_voxel_group);
                for (int i=0; i<num_voxel_group; i++) {
                    density_val_bkt_tmp(i) = density_val_bkt(voxel_group[i]);
                }

                // enforce projection constraint
                proj_density_val_bkt[proj_tmp] = density_val_bkt_tmp;

            }

            density_val = density_val_bkt;
            radial_density_val = radial_density_val_bkt;
            proj_density_val = proj_density_val_bkt;

            cor_obj_val = ComputeCorObjFun();
            radial_obj_val = ComputeRadialObjFun();
            proj_obj_val = ComputeProjObjFun();
            obj_val = cor_obj_val + radial_weight * radial_obj_val + proj_weight * proj_obj_val;


            if ((obj_val<obj_val_pre)||(obj_val==0)) {
                if (step<step_max) {
                    step = step/bkt_rate;
                }
                break;
            } else {
                step = step*bkt_rate;
                density_val = density_val_pre;
                radial_density_val = radial_density_val_pre;
            }
        }
        
        if (step<=step_thd) {
            density_val = density_val_pre;
            radial_density_val = radial_density_val_pre;

            cor_obj_val = ComputeCorObjFun();
            radial_obj_val = ComputeRadialObjFun();
            obj_val = cor_obj_val + radial_weight * radial_obj_val + proj_weight * proj_obj_val;;

            cout<<"Step size too small!"<<"\t"<<step<<endl;
            break;
        }

        cout<<"Ite: "<<ite<<"\t"<<"Step size: "<<step<<" "<<density_val.array().sum()<<endl;
        cout<<"Objective function: "<<cor_obj_val<<"\t"<<radial_obj_val<<"\t"<<proj_obj_val<<"\t"<<obj_val<<endl;
        //for (int i=0; i<est_radial_feat_all.size(); i++) {
        //    cout<<est_radial_feat_all(i)<<" ";
        //}
        //cout<<endl;
        
        VectorXd density_val_diff = density_val - density_val_pre;
        double cvg_val = density_val_diff.norm() / density_val.norm();
        
        if (cvg_val<cvg_thd) {
            num_cov +=1;
        }
        
        if (num_cov>=10) {
            cout<<"Convergence reached!"<<endl;
            break;
        }
        
    }
    
    if (ite>=max_ite) {
        cout<<"Max iteration reached!"<<endl;
    }
    
}


void UVT::WriteFile(char* density_output_file, char* obj_output_file) {

    ofstream write_result(density_output_file, ios_base::trunc);

    for (int i=0; i<num_loc; i++) {
        write_result<<loc_val_x[i]<<" "<<loc_val_y[i]<<" "<<loc_val_z[i]<<" "<<density_val(i)<<"\n";
    }
    
    write_result.close();

    ofstream write_result2(obj_output_file, ios_base::trunc);
    write_result2<<cor_obj_val<<" "<<radial_obj_val<<" "<<proj_obj_val<<" "<<obj_val;
    write_result2.close();

}


void UVT::ComputeCorObjFunMuti( vector<DtPair> lm_partition_block, VectorXd* est_cor_feat_seq_pt ) {}
double UVT::ComputeCorObjFun() { return 0.0; }

void UVT::ComputeRadialObjFunMuti( vector<int> radial_partition_block ) {}
double UVT::ComputeRadialObjFun() { return 0.0; }

void UVT::ComputeProjObjFunMuti( vector<DtPair> proj_partition_block ) {}
double UVT::ComputeProjObjFun() { return 0.0; }

void UVT::ComputeCorGradientMuti( vector<DtPair> lm_partition_block, VectorXd* gradient_val_seq_pt ) {}
void UVT::ComputeCorGradient() {}

void UVT::ComputeRadialGradientMuti( vector<int> radial_partition_block, VectorXd* gradient_val_seq_pt ) {}
void UVT::ComputeRadialGradient() {}

void UVT::ComputeProjGradientMuti( vector<DtPair> proj_partition_block, VectorXd* gradient_val_seq_pt ) {}
void UVT::ComputeProjGradient() {}

UVT::~UVT() {}
