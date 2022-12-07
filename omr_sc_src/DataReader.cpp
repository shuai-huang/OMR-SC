#include "DataReader.h"
#include "global.h"
#include "omp.h"

DataReader::DataReader(char* option_file) {
	// read option files and reset the parameter values
	ifstream read_option;
	read_option.open(option_file);
	string opt_tmp, opt_name;
	double opt_val;
	int opt_idx=0;
	while(read_option>>opt_tmp){
		if (opt_idx==0) {opt_name=opt_tmp; opt_idx++;}
		else {opt_val=atof(opt_tmp.c_str()); opt_idx=0;}
		
		if(opt_idx==0) {options[opt_name]=opt_val;}
	}
	read_option.close();
}

void DataReader::SetParameters() {

	// reset the parameters
	

	if (options.find("r_max")!=options.end())   {r_max=(int)options["r_max"];}
	if (options.find("N")!=options.end())   {N=(int)options["N"];}  
	if (options.find("L_max")!=options.end())   {L_max=(int)options["L_max"];}
    if (options.find("num_radial")!=options.end())   {num_radial=(int)options["num_radial"];}
		
    if (options.find("num_proj")!=options.end())   {num_proj=(int)options["num_proj"];}

	if (options.find("num_fourier_l_idx_seq")!=options.end())   {num_fourier_l_idx_seq=(int)options["num_fourier_l_idx_seq"];}
	if (options.find("num_fourier_cor_seq")!=options.end())   {num_fourier_cor_seq=(int)options["num_fourier_cor_seq"];}
	if (options.find("num_fourier_radial")!=options.end())   {num_fourier_radial=(int)options["num_fourier_radial"];}
	
	if (options.find("num_thread")!=options.end())   {num_thread=(int)options["num_thread"];}
	if (options.find("init_type")!=options.end()) {init_type=(int)options["init_type"];}
    if (options.find("max_alt_ite")!=options.end())   {max_alt_ite=(int)options["max_alt_ite"];}
	if (options.find("max_ite")!=options.end())   {max_ite=(int)options["max_ite"];}
    if (options.find("density_all")!=options.end()) {density_all=options["density_all"];}
    if (options.find("radial_weight")!=options.end()) {radial_weight=options["radial_weight"];}
    if (options.find("proj_weight")!=options.end()) {proj_weight=options["proj_weight"];}
    if (options.find("voxel_std")!=options.end())  {voxel_std=options["voxel_std"];}
	if (options.find("bkt_rate")!=options.end())  {bkt_rate=options["bkt_rate"];}
	if (options.find("step_ori")!=options.end())  {step_ori=options["step_ori"];}
	if (options.find("step_thd")!=options.end())  {step_thd=options["step_thd"];}
    if (options.find("step_max")!=options.end())  {step_max=options["step_max"];}
	if (options.find("cvg_thd")!=options.end())   {cvg_thd=options["cvg_thd"];}

    // compute voxel_var
    voxel_var = voxel_std*voxel_std;
    // just reset the num_thread value if the pre-specified value is larger than the actual number of threads in the machine
    int num_thread_machine = omp_get_max_threads();
    if (num_thread>num_thread_machine) {
        num_thread = num_thread_machine;
    }
    
    cout<<"No. of threads:          "<<num_thread<<endl;
    cout<<"No. of Raidal feature:   "<<num_radial<<endl;
    cout<<"No. of Proj feature:     "<<num_proj<<endl;
    cout<<"No. of Fourier Cor feat Seq:     "<<num_fourier_cor_seq<<endl;
    cout<<"No. of Fourier Raidal feature:   "<<num_fourier_radial<<endl;
    cout<<"No. of Fourier l_idx sequence:   "<<num_fourier_l_idx_seq<<endl;
    cout<<"L max:                   "<<L_max<<endl;
    cout<<"Domain size:             "<<N<<endl;
    cout<<"Initialization type:     "<<init_type<<endl;
    cout<<"Max alt iteration:       "<<max_alt_ite<<endl;
    cout<<"Max iteration:           "<<max_ite<<endl;
    cout<<"Density all:             "<<density_all<<endl;
    cout<<"Radial weight:           "<<radial_weight<<endl;
    cout<<"Voxel std:               "<<voxel_std<<endl;
    cout<<"Backtracking rate:       "<<bkt_rate<<endl;
    cout<<"Stepsize init, min, max: "<<step_ori<<"\t"<<step_thd<<"\t"<<step_max<<endl;
    cout<<"Convergence threshold:   "<<cvg_thd<<endl;
    cout<<endl;
}


void DataReader::ReadRadial(char* radial_file) {

    // read radial data files and save it in a 2-dimensional vector mat
    ifstream read_radial;

    read_radial.open(radial_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_radial, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        radial_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }   

    read_radial.close();
    
    if (line_num != num_radial) {
        cout<<"Radial data dimension mismatch!!!"<<endl;
        cout<<line_num<<" "<<num_radial<<endl;
        abort();
    }

    cout<<"Reading radial data finished."<<endl;
}

void DataReader::ReadRadialMap(char* radial_map_file) {

    // read radial map data files and save it in a 2-dimensional vector mat
    ifstream read_radial_map;

    read_radial_map.open(radial_map_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_radial_map, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        radial_map_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }   

    read_radial_map.close();
    
    if (line_num != num_radial) {
        cout<<"radial map data dimension mismatch!!!"<<endl;
        abort();
    }   

    cout<<"Reading radial map data finished."<<endl;

}


void DataReader::ReadProj(char* proj_file) {

    // read proj data files and save it in a 2-dimensional vector mat
    // the points locations that have a nonzero value
    ifstream read_proj;

    read_proj.open(proj_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_proj, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        proj_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }

    read_proj.close();

    if (line_num != num_proj) {
        cout<<"Projection data dimension mismatch!!!"<<endl;
        cout<<line_num<<"\t"<<num_proj<<endl;
        abort();
    }

    cout<<"Reading proj data finished."<<endl;
}


void DataReader::ReadFourierCor(char* fourier_cor_file) {
	
	// read data files and save it in a 2-dimensional vector mat
	ifstream read_fourier_cor;
	
    read_fourier_cor.open(fourier_cor_file);
	istringstream istr;
	string str;
	int line_num=0;
    
    while(getline(read_fourier_cor, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        fourier_cor_data.push_back(tmpvec);

	    tmpvec.clear();
	    istr.clear();
	    line_num++;
    }
	
	read_fourier_cor.close();
	
	if (line_num != num_fourier_cor_seq) {
	    cout<<"Fourier Correlation data dimension mismatch!!!"<<endl;
        cout<<line_num<<" "<<num_fourier_cor_seq<<endl;
	    abort();
    }

	cout<<"Reading fourier correlation data finished."<<endl;
}


void DataReader::ReadFourierRadialMap(char* fourier_radial_map_file) {

    // read radial map data files and save it in a 2-dimensional vector mat
    ifstream read_fourier_radial_map;

    read_fourier_radial_map.open(fourier_radial_map_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_fourier_radial_map, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        fourier_radial_map_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }   

    read_fourier_radial_map.close();
    
    if (line_num != num_fourier_radial) {
        cout<<"Fourier radial map data dimension mismatch!!!"<<endl;
        abort();
    }   

    cout<<"Reading fourier radial map data finished."<<endl;

}


void DataReader::ReadFourierLIdx(char* fourier_l_idx_file) {

    // read l_idxs and save it in a vector
    ifstream read_fourier_l_idx;

    read_fourier_l_idx.open(fourier_l_idx_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_fourier_l_idx, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        fourier_l_idx_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }   

    read_fourier_l_idx.close();
    
    if (line_num != num_fourier_l_idx_seq) {
        cout<<"Fourier L idx data dimension mismatch!!!"<<endl;
        cout<<line_num<<" "<<num_fourier_l_idx_seq<<endl;
        abort();
    }

    cout<<"Reading fourier l_idx data finished."<<endl;
}


void DataReader::ReadRadialLGWTWTS(char* radial_lgwt_weight_file) {

    // read lgwt and save it in a vector
    ifstream read_radial_lgwt_weight;

    read_radial_lgwt_weight.open(radial_lgwt_weight_file);
    istringstream istr;
    string str;
    int line_num=0;

    while(getline(read_radial_lgwt_weight, str)) {
        istr.str(str);
        double tmp;
        vector<double> tmpvec;
        while(istr>>tmp) {
            tmpvec.push_back(tmp);
        }
        radial_lgwt_weight_data.push_back(tmpvec);

        tmpvec.clear();
        istr.clear();
        line_num++;
    }   

    read_radial_lgwt_weight.close();
    
    if (line_num != num_radial) {
        cout<<"Radial Lgwt data dimension mismatch!!!"<<endl;
        abort();
    }   

    cout<<"Reading radial lgwt weights data finished."<<endl;
}



vector< vector<double> > DataReader::getRadialData() {return radial_data;}
vector< vector<double> > DataReader::getRadialMapData() {return radial_map_data;}
vector< vector<double> > DataReader::getProjData() {return proj_data;}

vector< vector<double> > DataReader::getFourierCorData() {return fourier_cor_data;}
vector< vector<double> > DataReader::getFourierRadialMapData() {return fourier_radial_map_data;}
vector< vector<double> > DataReader::getFourierLIdxData() {return fourier_l_idx_data;}

vector< vector<double> > DataReader::getRadialLGWTWTSData() {return radial_lgwt_weight_data;}

void DataReader::RemoveRadialData() {radial_data.clear();}
void DataReader::RemoveRadialMapData() {radial_map_data.clear();}
void DataReader::RemoveProjData() {proj_data.clear();}

void DataReader::RemoveFourierCorData() {fourier_cor_data.clear();}
void DataReader::RemoveFourierRadialMapData() {fourier_radial_map_data.clear();}
void DataReader::RemoveFourierLIdxdata() {fourier_l_idx_data.clear();}

void DataReader::RemoveRadialLGWTWTSData() {radial_lgwt_weight_data.clear();}
