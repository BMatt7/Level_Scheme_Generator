#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <stdlib.h>
#include "Riostream.h"
#include <math.h>       /* ceil */


/* How To Run //

Compile with: g++ -std=c++11 -o <executable_name> Level_Scheme_Latex_Generator.C


./<executable_name> -l <Levels_List_Data_Filename.txt> -t <Transitions_List_Data_Filename.txt> -max <Maximum Energy of Level Scheme> -of <Output_Filename.txt>

Example Run Line:

./generate -l 208Po_Levels.txt -t 208Po_Transitions.txt -max 2000.0 -of 208Po_Latex.txt


Additional Notes:
* Maximum Energy is needed to scale the diagram according to the maximum energy
* shift changes spacing between transitions
* flat_part changes length of kink for levels that need it (may affect placement of labels if changed)
* line_length changes width of the level scheme

* Energy states should be in increasing energy order
* Transitions can be in any order but will be placed in the order given, so may affect placing

*/

using namespace std;

int main(int argc, char* argv[])
{

	string max_energy_read = "empty";
	string levels_filename = "empty";
	string transitions_filename = "default.txt";
	string out_filename = "default.txt";
	bool thick_arrows = false;
	string thicc;
	
	for(int i = 0; i < argc; i++)
	{

		if((string)argv[i] == "-max" || (string)argv[i] == "--maxenergy") max_energy_read = (string)argv[i+1];
		if((string)argv[i] == "-l" || (string)argv[i] == "--levels") levels_filename = (string)argv[i+1];
		if((string)argv[i] == "-t" || (string)argv[i] == "--transitions") transitions_filename = (string)argv[i+1];
		if((string)argv[i] == "-of" || (string)argv[i] == "--outputfile") out_filename = (string)argv[i+1];
		if((string)argv[i] == "-thicc" || (string)argv[i] == "--thiccarrows"){
		
			thick_arrows = true;
			thicc = (string)argv[i+1];
			
		}
		
	}

	double max_energy = stof(max_energy_read);
	
	max_energy = max_energy/25.0;

    ifstream lfs (levels_filename.c_str());
    ifstream tfs (transitions_filename.c_str());
    ofstream ofs (out_filename.c_str());

    
    if(lfs.fail()){
        cout << "Couldn't open " << levels_filename << " for reading"<<endl;
        exit(0);
    }
    if(tfs.fail()){
        cout << "Couldn't open " << transitions_filename << " for reading"<<endl;
        exit(0);
    }
    
	ofs<<"\\documentclass[11pt]{article}"<<endl;
	ofs<<"\\usepackage[T1]{fontenc}"<<endl;
	ofs<<"\\usepackage[version=3]{mhchem}"<<endl;
	ofs<<"\\usepackage{siunitx}"<<endl;
	ofs<<"\\usepackage{tikz}"<<endl;
	ofs<<"\\usepackage[margin=0.5cm]{geometry}"<<endl;
	ofs<<"\\usetikzlibrary{arrows}"<<endl;
	ofs<<"\\usetikzlibrary{shapes.geometric}"<<endl;
	ofs<<"\\usetikzlibrary{arrows.meta,arrows}"<<endl;

	ofs<<"\\begin{document}"<<endl;

	ofs<<"\\begin{tikzpicture}"<<endl;
	
	ofs<<"%%Energy States"<<endl;
	
	
	const char* format = "%lf %lf %lf %lf";
	const char* format_2 = "%lf %lf %lf";
	
	int num_states = 0;
	int num_transitions = 0;
	bool new_scheme = false;
	
    string line;
    
    double flat_part = 1.5;
    double line_length = 16; // max 18
   	double shift = 0.5; // in cm
   	double max_thickness;
   	if(thick_arrows) max_thickness = stof(thicc); // in mm 

    vector<vector<double>> a;
    vector<double> energy_state_vec;
    vector<double> spin_vec;
    vector<double> parity_vec;
    vector<double> line_type_vec;
    
    vector<vector<double>> b;
    vector<double> initial_state_vec;
    vector<double> final_state_vec;
    vector<double> intensity_vec;
    
    
	int n = 0;
	double prev_energy = -10000.0;
	double height = 0.0;
	
	double correction = 0.55;
	
	vector<double> levels;
    vector<double> heights_vec;
    	
    double prev_level_height = -10000.0;
    double prev_kink_height = -10000.0;
    
    
    
    double energy_of_state;
    double spin;
    double parity;
    double line_type;
    
    double initial_state;
    double final_state;
    double intensity;

    

    while(lfs.good()){
    
        getline(lfs,line,'\n');
        if (lfs.eof()) break;
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format, &energy_of_state, &spin, &parity, &line_type);
		
		energy_state_vec.push_back(energy_of_state);
		spin_vec.push_back(spin);
		parity_vec.push_back(parity);
		line_type_vec.push_back(line_type);
		
		num_states++;
		
	}
	
	a.push_back(energy_state_vec);
	a.push_back(spin_vec);
	a.push_back(parity_vec);
	a.push_back(line_type_vec);
	
	
	while(tfs.good()){
    
        getline(tfs,line,'\n');
        if (tfs.eof()) break;
        if(line[0] == '#') continue;
        sscanf(line.c_str(),format_2, &initial_state, &final_state, &intensity);
		
		initial_state_vec.push_back(initial_state);
		final_state_vec.push_back(final_state);
		intensity_vec.push_back(intensity);
		
		num_transitions++;
	
	}
	
	b.push_back(initial_state_vec);
	b.push_back(final_state_vec);
	b.push_back(intensity_vec);
	
	
	for(int i = 0; i < num_states; ++i){

		cout<<a[0][i]<<" "<<a[1][i]<<" "<<a[2][i]<<" "<<a[3][i]<<" "<<endl;

	}
		
	for(int i = 0; i < num_transitions; ++i){
	
		cout<<endl;

		cout<<b[0][i]<<" "<<b[1][i]<<" "<<b[2][i]<<" "<<endl;

	}
	
	int shift_checker[num_states];
	
	for(int i = 0; i < num_states; ++i) shift_checker[i] = 0;
    
    bool energy_found;
    int value = 0;
    int num_shifts = 0;
    
    int start_energy = -1;
    int end_energy = -1;
    
    int shift_amount;
	
	bool transitions_unplaced = true;
	
	int last_placed = 0;
	
	
	while(transitions_unplaced){
	
	prev_level_height = -10000.0;
    prev_kink_height = -10000.0;
    n = 0;
	
	for(int i = 0; i < num_states; ++i){
	
	    if(a[3][i] == 0) ofs<<"\\draw ";
        else if(a[3][i] == 1) ofs<<"\\draw[very thick] ";
        else if(a[3][i] == 2) ofs<<"\\draw[densely dashed] ";
        
        if(a[0][i]/max_energy <= prev_level_height + 0.55){
        
        	height = prev_level_height + 0.55;
            
            if (a[0][i]/max_energy <= (prev_kink_height + 0.1)) correction = 0.4;
       		else correction = height - a[0][i]/max_energy;
       		
        }
        else{
        
        	height = a[0][i]/max_energy;
        	correction = 0.0;
        }
        

        
        ofs<<"(1,"<<height<<") -- ("<<1+flat_part<<","<<height<<") node[pos=0.2,above]{"<<a[1][i]<<"\\textsuperscript{";
        
        if(a[2][i] == 1) ofs<<"+}} -- (";
        else if(a[2][i] == 0) ofs<<"-}} -- (";
        
        
        ofs<<flat_part+1.5<<","<<height-correction<<")  -- ("<<line_length-flat_part-0.5<<","<<height-correction<<") -- ("<<line_length-flat_part<<","<<height<<") -- ("<<line_length<<","<<height<<") node[pos=0.5,above] {\\ce{"<<a[0][i]<<"}};";
        
        ofs<<endl;
        ofs<<endl;
        
		prev_energy = a[0][i];
		prev_level_height = height; 
		prev_kink_height = height - correction;
		
		levels.push_back(prev_energy);
		heights_vec.push_back(prev_kink_height);
		
        n++;
        
	}
	
		
	for(int j = last_placed; j < num_transitions; ++j){
	
				cout<<b[0][j]<<" "<<b[1][j]<<" "<<b[2][j]<<" "<<endl;
			
				energy_found = false;
				value = 0;
				
				num_shifts = 0;
				
				start_energy = -1;
				
				if(thick_arrows && b[2][j] >= 15) shift_amount = 1 + ceil((max_thickness*(b[2][j]/100.0))/(shift*10));
				else if(thick_arrows && b[2][j-1] >= 15) shift_amount = 1 + ceil((max_thickness*(b[2][j-1]/200.0))/(shift*10));
				else shift_amount = 1;
				
				if(thick_arrows && b[2][j] >= 15 && j == 0) shift_amount = 1 + ceil((max_thickness*(b[2][j-1]/50.0))/(shift*10));
				
				while(!energy_found){

					cout<<"Up En: "<<b[0][j]<<"   Low En: "<<b[1][j]<<" Level: "<<levels[value]<<"  "<<value<<endl;

					if(levels[value] < b[1][j]){

					 value++;
					 
					 continue;
					 
					}
					else if (levels[value] >= b[1][j] && levels[value] <= b[0][j]){
			
						if (start_energy == -1) start_energy = value;
						end_energy = value;
				
						shift_checker[value] += shift_amount;
				
						cout<<"Shift Checker Increased: "<<value<<" "<<shift_checker[value]<<endl;
					 
						if (num_shifts < shift_checker[value]){
						 
						 
							num_shifts = shift_checker[value];
			
							cout<<"num_shifts = "<<num_shifts<<endl;
			
						}
						else if (num_shifts > shift_checker[value]) shift_checker[value] = num_shifts;
				
					}
					 
					if(b[0][j] == levels[value]) energy_found = true; 
			
					value++;

				}
				
				if((line_length-flat_part-(num_shifts*shift)-1.0) <= flat_part) new_scheme = true;
				//else if(b[2][j] >= 15 && line_length-flat_part-(num_shifts*shift)-(1.4*(max_thickness*(b[2][j]/100.0))) <= flat_part) new_scheme = true;
				if(new_scheme){
				
					cout<<"NEW LEVEL SCHEME WOO!!!!!"<<endl;
					last_placed = j;
					break;
					
				}
		
				if(!thick_arrows){

					if(b[2][j] >= 15) ofs<<"\\draw[red,->] ";
					else if(b[2][j] >= 5) ofs<<"\\draw[blue,->] ";
					else ofs<<"\\draw[->] ";
			
					ofs<<"("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[end_energy]<<") -- ("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[start_energy]<<") node[pos=0.0,right,rotate=60,black,fill=white]{"<<b[0][j]-b[1][j]<<" "<<b[2][j]<<"};";
	
		
				}

				else if(thick_arrows){
				
					if(b[2][j] >= 15){
				
						ofs<<"\\draw[>={Triangle[width="<<1.4*(max_thickness*(b[2][j]/100.0))<<"mm,length="<<2.0*(heights_vec[end_energy] - heights_vec[start_energy])<<"mm]}, line width="<<max_thickness*(b[2][j]/100.0)<<"mm,->] ";
				
						ofs<<"("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[end_energy]<<") -- ("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[start_energy]<<") node[pos=0.5,black,fill=white]{"<<b[0][j]-b[1][j]<<"};";
	
					}	
					else{
				
						ofs<<"\\draw[->] ";
				
						ofs<<"("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[end_energy]<<") -- ("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[start_energy]<<") node[pos=0.0,right,rotate=60,black,fill=white]{"<<b[0][j]-b[1][j]<<" "<<b[2][j]<<"};";
	
				
					}
		
				}
		
		
				//ofs<<"("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[end_energy]<<") -- ("<<line_length-flat_part-(num_shifts*shift)<<","<<heights_vec[start_energy]<<") node[pos=0.0,right,rotate=60,black,fill=white]{"<<b[0]-b[1]<<" "<<b[2]<<"};";
	
				ofs<<endl;
				ofs<<endl;
			
			
			}
			
			if(!new_scheme) transitions_unplaced = false;
			
			if(new_scheme){
			
				ofs<<"\\end{tikzpicture}"<<endl;
				ofs<<endl;
				ofs<<"\\pagebreak"<<endl;
				ofs<<endl;
				ofs<<"\\begin{tikzpicture}"<<endl;
				
				ofs<<endl;
				ofs<<endl;
				ofs<<endl;
				ofs<<endl;
				ofs<<endl;
				ofs<<endl;
				
				for(int i = 0; i < n; ++i) shift_checker[i] = 0;
    
    			value = 0;
    			num_shifts = 0;
    
   				start_energy = -1;
    			end_energy = -1;
    			
    			new_scheme = false;
							
			}
	
	}
    
        

	
	ofs<<"\\end{tikzpicture}"<<endl;
	ofs<<"\\end{document}"<<endl;
			
			
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"File complete!"<<endl;
	cout<<""<<endl;
    cout<<"Data written to file: "<<out_filename<<endl;
   	cout<<"----------------------------------------------------------"<<endl;
   	
   	

}
