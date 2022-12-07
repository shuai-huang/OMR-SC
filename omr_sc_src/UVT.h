#ifndef UVT_H
#define UVT_H

#include <iostream>
#include <map>
#include <vector>
#include <random>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "DataReader.h"
#include "global.h"

#include "math.h"

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef std::array<int,3> DtTriplet;
typedef std::array<int,2> DtPair;

struct voxel {
    DtTriplet coordinate;
    double radius;
};

using namespace std;

class UVT
{    
    protected:
    
        map<int, vector<int>> voxel_radial_idx;     // the radial index, the group of voxel indices
        
        map<int, double> radial_map;                // the radial index (start from 0 consequitively), the actual radius value
        map<int, double> fourier_radial_map;        // the fourier radial index (start from 0), the actual fourier radius value
        
        map<int, vector<DtPair>> lm_partition_idx;  // the partition index, the partition for multi-threading of (l,m) pair
        int num_lm_thread_assign;                   // the number of assigned multi_threads for (l,m) pair
       
        // should i use integer (radial index) here or double (radius value) here?
        // use radial_index here
        map<int, vector<int>> radial_partition_idx; // the partition idx, the partition for multi-threading of the actual radii values, not indices
        int num_radial_thread_assign;               // the number of assigned multi threads for radii
        map<int, vector<DtPair>> proj_partition_idx;   // the partition idx, the partition for multi-threading of the proj locations
        int num_proj_thread_assign;

        map<int, voxel> voxel_map;                  // the voxel indices, the voxel info, the orders are important
        map<DtTriplet, int> voxel_map_reverse;      // the coordinates, the voxel indices
        map<int, int> r_candidate;                  // the radii that holds candidate voxels, "1" indicating it is true
        
        map<DtPair, vector<int>> fourier_radial_l_idx_seq;  // the (l,m) pair, the fourier radial indices associated with each (l,m) pair
        
        map<DtTriplet, VectorXd> radial_alpha_sh_seq;     // the (radial index, l, m), the multiplication of alpha and sh values corresponding to group of voxel indices
        map<int, VectorXd> radial_beta_seq;         // the radial index, the beta values corresponding to group of voxel indices
       
        map<DtTriplet, int> cor_feat_all_loc;       // the (r1, r2, l) triplet, the location/index at the vectorized correlation feature cor_feat_all
        map<DtTriplet, int> fourier_cor_feat_all_loc;   // the (k1, k2, l) triplet, the location/index at the vectorized fourier correlation feature fourier_cor_feat_all
        
        
        VectorXd density_val;                       // the density values ordered according to voxel indices
        MatrixXd loc_val;                           // num_loc by 3 matrix, corresponding to ordered cartesian locations
        vector<int> loc_val_x;
        vector<int> loc_val_y;
        vector<int> loc_val_z;
        MatrixXd loc_sph_val;                       // num_loc by 3 matrix, corresponding to ordered spherical locations
        VectorXd radius_val;                        // the radius values ordered according to voxel indices

        
        map<int, VectorXd> radial_density_val;      // the radial index, the density values of the corresponding group of voxels

        map<DtPair, VectorXd> proj_density_val;     // the projection index, the desntiy values of the voxels
        map<DtPair, vector<int>> voxel_proj_idx;   // the projection index, the group of voxel indices

        map<DtPair, double> proj_feat;      // the (x,y) pair, the projection density
        vector<DtPair> proj_feat_idx;          // a sequence of all proj indices corresponding to proj_feat
        VectorXd proj_feat_all;             // to make it easier for block computation
        int num_proj_feat;
        map<DtPair, int> proj_feat_all_loc;    // the proj index, the location/index at the vectorized proj_feat_all
        double max_proj_feat;

        
        VectorXd cor_gradient_val;            // the radial index, the gradients corresponding to each voxel
        VectorXd radial_gradient_val;
        VectorXd proj_gradient_val;
        VectorXd gradient_val;
        
        
        VectorXd radial_lgwt_range;
        VectorXd radial_lgwt_weights;
        map<int, MatrixXd> radial_lgwt_mat;     // the l index, the integration matrix via lgwt

        int num_cor_feat;                 // the number of single corelation features

        map<DtPair, VectorXd> fourier_cor_feat; // the (l,m) pair, the fourier autocorrelation features ordered according to fourier radius values                                  
        map<DtPair,int> fourier_cor_feat_lm_unique; // this should be the same for both fourier correlation and spatial correlation
        

        map<int, MatrixXd> fourier_cor_feat_mat_all;              // the fourier correlation feature ordered according to l. "int" corresponds to l value; "MatrixXd" corresponds to the fourier_cor_features, the rows of "MatrixXd" correspond to fourier radial values, the columns of "MatrixXd" corresponds to m values
        // note this is real from decomposition, but we set is to complex format anyway
        map<int, MatrixXd> fourier_orth_mat_all;        // the fourier orthogonal matrix ordered according to l, "int" corresponds to l value; "MatrixXd" corresponds to the (2l+1) by (2l+1) fourier orthogonal matrix, so for l=0, "MatrixXd" is always 1
        // note this is complex from estimation
        VectorXd fourier_cor_feat_all;              // the fourier correlation feature ordered according to cor_feat_lm_unique
                                                    // to make it easier for block computation
        VectorXd fourier_cor_feat_all_transform;    // transformed fourier_cor_feat_all with fourier_orth_mat_all
        int num_fourier_cor_feat;                 // the number of single fourier corelation features
        
        map<int, double> radial_feat;       // the radial index, the radial feature
        vector<int> radial_feat_idx;    // a sequence of all radii indices according to radial_feat, not indices, some radial indices in the 3D census might be missing
        VectorXd radial_feat_all;   // to make it easier for block computation
        int num_radial_feat;
        map<int, int> radial_feat_all_loc;  // the radial index, the location/index at the vectorized radial_feat_all
        
        map<int, double> fourier_radial_feat;       // the fourier radial index, the fourier radial feature
        vector<int> fourier_radial_feat_idx;    // a sequence of all fourier radii indices according to radial_feat, not indices, some radial indices in the 3D census might be missing
        VectorXd fourier_radial_feat_all;   // to make it easier for block computation
        int num_fourier_radial_feat;
        map<int, int> fourier_radial_feat_all_loc;  // the radial index, the location/index at the vectorized radial_feat_all
        
       
        map<int, MatrixXd> est_cor_feat_mat_all;
        VectorXd est_cor_feat_all;
        VectorXd est_radial_feat_all;

        VectorXd est_proj_feat_all;
        
        map<int, MatrixXd> est_fourier_cor_feat_mat_all;
        VectorXd est_fourier_cor_feat_all;
        
        int N_half;                         // N is an odd number, this is (N-1)/2
        int num_loc;                        // the total number of candidate voxels
        int num_radial_cor_feat;            // the number of radii in cor_feat;
        int num_fourier_radial_cor_feat;            // the number of fourier radii in cor_feat;

        int comass_idx;                      // the index of the location corresponding to the center of mass
        
        double step;
        double cor_obj_val;
        double radial_obj_val;
        double proj_obj_val;
        double obj_val;
        
        char* density_init_file;
        
    public:
        UVT();
        
        void SetInitFile(char* init_file);
        void SetData(DataReader* data_reader);
        void Census3DSpace();

        void Initialization();
        void OrthoMatrixSync();
        void GradientDescent();
        void WriteFile(char* density_output_file, char* obj_output_file);
        
        VectorXd ProjectOntoCvxSet(VectorXd smp_vec, double num_smp_proj);
        
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
        
        virtual ~UVT();
        
};

#endif // UVT_H
