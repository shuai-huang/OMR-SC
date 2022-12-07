#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

#include "global.h"
#include "UVT.h"
#include "PGD.h"


using namespace std;
using namespace Eigen;

using std::ifstream;
using std::string;

void process_cmdline(char *argv[], char **radial_file, char **radial_map_file, char **radial_lgwt_weight_file, char **proj_file, char **fourier_cor_file, char **fourier_radial_map_file, char **fourier_l_idx_file, char **option_file, char **init_file, char **output_file, char **obj_file) {
	
    *radial_file                = argv[1];
    *radial_map_file            = argv[2];
    *radial_lgwt_weight_file    = argv[3];
    *proj_file                  = argv[4];
	*fourier_cor_file	        = argv[5];
    *fourier_radial_map_file    = argv[6];
    *fourier_l_idx_file         = argv[7];
	*option_file	            = argv[8];
	*init_file                  = argv[9];
	*output_file	            = argv[10];
    *obj_file                   = argv[11];
}

int main(int argc, char *argv[]){

	char *radial_file, *radial_map_file, *radial_lgwt_weight_file, *proj_file, *fourier_cor_file, *fourier_radial_map_file, *fourier_l_idx_file, *option_file, *init_file, *output_file, *obj_file;
	if (argc==12) {
		process_cmdline(argv, &radial_file, &radial_map_file, &radial_lgwt_weight_file, &proj_file, &fourier_cor_file, &fourier_radial_map_file, &fourier_l_idx_file, &option_file, &init_file, &output_file, &obj_file);
	} else {
		cout<<"Command line error!"<<endl;
		abort();
	}
	
	DataReader data_input(option_file);
	data_input.SetParameters();
    data_input.ReadRadial(radial_file);
    data_input.ReadRadialMap(radial_map_file);
    data_input.ReadProj(proj_file);
    data_input.ReadFourierCor(fourier_cor_file);
    data_input.ReadFourierRadialMap(fourier_radial_map_file);
    data_input.ReadFourierLIdx(fourier_l_idx_file);
    data_input.ReadRadialLGWTWTS(radial_lgwt_weight_file);

	UVT *uvt_pgd;

    uvt_pgd = new PGD();

    uvt_pgd->SetInitFile(init_file);
    uvt_pgd->SetData(&data_input);
    uvt_pgd->Census3DSpace();
    uvt_pgd->Initialization();
    
    for (int alt_ite=1; alt_ite<=max_alt_ite; alt_ite++) {
        uvt_pgd->OrthoMatrixSync();
        uvt_pgd->GradientDescent();
    }
    uvt_pgd->WriteFile(output_file, obj_file);

    delete uvt_pgd;

	exit(0);
}
