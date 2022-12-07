#include <map>
#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>

using namespace std;
using std::string;
using std::cerr;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


class DataReader {
	
	private:
		map<string, double> options;
        vector< vector<double> > radial_data;
        vector< vector<double> > radial_map_data;
        vector< vector<double> > proj_data;
        
		vector< vector<double> > fourier_cor_data;
        vector< vector<double> > fourier_radial_map_data;
        vector< vector<double> > fourier_l_idx_data;

        vector< vector<double> > radial_lgwt_weight_data;
        
	public:
		DataReader(char* option_file);
        vector< vector<double> > getRadialData();
        vector< vector<double> > getRadialMapData();
        vector< vector<double> > getProjData();
        
		vector< vector<double> > getFourierCorData();  
        vector< vector<double> > getFourierRadialMapData();
        vector< vector<double> > getFourierLIdxData();
        
        vector< vector<double> > getRadialLGWTWTSData();
		
		void SetParameters();
                                            // the order of the correlation feature should be consistent with the l_idx data
        void ReadRadial(char* radial_file); // in each line, the first is radial index, the second is the radial feature
        void ReadRadialMap(char* radial_map_file); // in each line the first is radial index, the second is the radius value
                                            // make sure radial feature reaches radial_max, we need to this to determine candidate locations
        void ReadProj(char* proj_file);     // in each line, the first is the x coordinate, the second is the y coordinate, the third is the projection value
                                            
        void ReadFourierCor(char* fourier_cor_file);       // in each line, the first is ||r_1||_2, the second is ||r_2||_2, the rest are ordered correlation features corresponding to l_indices ranging from 0 to L_max, but of course, the l_indices can be randomly sampled as well
                                            // the order of the correlation feature should be consistent with the l_idx data
        void ReadFourierRadialMap(char* fourier_radial_map_file); // in each line the first is radial index, the second is the radius value
                                            // make sure radial feature reaches radial_max, we need to this to determine candidate locations
        void ReadFourierLIdx(char* fourier_l_idx_file);   // in each line, the first is ||r_1||_2, the second is ||r_2||_2, the rest are the l_indices ranging from 0 to L_max, but of course, the l_indices can be randomly sampled as well
                                            // the order of the l_idx data should be consistent with cor_feat
                                            
        void ReadRadialLGWTWTS(char* radial_lgwt_weight_file);
        
		void RemoveRadialData();
        void RemoveRadialMapData();
		void RemoveProjData();
		
		void RemoveFourierCorData();
        void RemoveFourierRadialMapData();
		void RemoveFourierLIdxdata();
		
		void RemoveRadialLGWTWTSData();
		
};
