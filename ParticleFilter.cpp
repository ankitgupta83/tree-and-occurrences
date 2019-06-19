//Particle filter approximation of the likelihood: Implements the method from: Vaughan, T.G., Leventhal, G.E., Rasmussen, D.A., Drummond, A.J., Welch, D., Stadler, T., 2018. Estimating epidemic incidence and prevalence from genomic data. Molecular Biology and Evolution, msz106, //https://doi.org/10.1093/molbev/msz106

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <random>
#include "BirthDeath.h"
using namespace std;


std::default_random_engine generator;
std::uniform_real_distribution<double> Random(0.0,1.0);

/* Global variables */
int NextreactionIndex;
double TimeIncrement;
double PropensitySum;
double DecayFactor;
double beta;
double mu;
double psi;
double omega;
double rho;
double PF_Log_LikelihoodSquared;




double nchoosek(int a , int b){

    if(a < b){
        return 0;
    }
    else if(a==b || b == 0 ) {
        return 1;
    }
    else{
        return nchoosek(a-1,b-1) + nchoosek(a-1,b);
    }
}


double Myfactorial(int a){

    if(a < 2){
        return 1;
    }
    else{
        return a*Myfactorial(a-1);
    }
}

void FindNextReaction(int * state){
    
    PropensitySum = CalculatePropensities(state);
    
    if(PropensitySum == 0){
        TimeIncrement = MAX_TIME_VALUE;
        NextreactionIndex = -1;
        return;
    }
    
    TimeIncrement = -log(Random(generator))/PropensitySum;
    double sum = 0;
    int i;
    double unif = Random(generator);
    for(i = 0; i < NUMBER_OF_REACTIONS ; i++ ){                
        sum = sum + propensity[i]/PropensitySum;        
        if(unif < sum){
            break;
        }                
    }
    NextreactionIndex = i;
}



int GenerateTrajectory(double Finaltime, int initialstate, int NumObservedLineages){

    int state[NUMBER_OF_SPECIES];
    double t = 0;

    double Num;
    double temp_state;
    bool identifier = true;

    for(int i = 0; i < NUMBER_OF_SPECIES ; i++){
        state[i] = initialstate;
    }
    while(1){
        FindNextReaction(state);
        t = t + TimeIncrement;
        if(t > Finaltime){
            break;
        }
        else if(NextreactionIndex >=0){
            for( int j = 0 ; j < NUMBER_OF_SPECIES ; j++){
                state[j] = state[j] + StoichiometryMatrix[j][NextreactionIndex];
            }

            if((NextreactionIndex == 2) || (NextreactionIndex == 3) ||(state[0] < NumObservedLineages)){ // no sampling events or omega events or insufficient number of lineages
                identifier = false;
                break;
            }
            else if((NextreactionIndex == 0) && (NumObservedLineages>=2)){
                Num = (double) NumObservedLineages;
                temp_state = (double) state[0];
                DecayFactor = DecayFactor * (1 - (Num * (Num - 1)) / (temp_state * (temp_state - 1)));

            }
        }
    }

    if(!identifier){
        return -1;
    }
    else {
        return state[0];
    }
}



double EstimateLikelihood(int NumEvents, double * ObservedEventTimes, int * ObservedEventCodes){

    double PF_Log_Likelihood = 0;
    double Weights[MAX_NUM_TOTAL_LINEAGES];
    double weight;
    double unif;
    double temp_state, NumLineages;
    int InitialStates[NUMBER_OF_PARTICLES];
    int final_state;
    double weightsum;
    
    int NumObservedLineages = 1;
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
        InitialStates[i] = NumObservedLineages;
    }

    double time =  ObservedEventTimes[(NumEvents-1)]; // This is the time of origin.

    for(int J = (NumEvents - 2) ; J >=0; J--) {

        //initialize the weight vector
        for (int i = 0; i < MAX_NUM_TOTAL_LINEAGES; i++) {
            Weights[i] = 0;
        }
        weightsum = 0;
        

        for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
            DecayFactor = 1;
    
            final_state = GenerateTrajectory(time - ObservedEventTimes[J] , InitialStates[i], NumObservedLineages);
            temp_state = (double) final_state;
            NumLineages = (double) NumObservedLineages;
            
            switch(ObservedEventCodes[J]) {
                case 0  :
                    //"Branching (beta event)";
                    if (final_state >= NumObservedLineages && final_state< (MAX_NUM_TOTAL_LINEAGES-1)) {
                        weight = DecayFactor * beta/(temp_state+1);
                        Weights[final_state+1] = Weights[final_state+1] + weight;
                        weightsum = weightsum + weight;
                    }
                    break;
                case 1  :
                    //"Sampling with Removal Event (psi event)";
                    if (final_state >= NumObservedLineages && final_state < (MAX_NUM_TOTAL_LINEAGES)) {
                        weight = DecayFactor*psi*temp_state;
                        Weights[final_state-1] = Weights[final_state-1] + weight;
                        weightsum = weightsum + weight;
                    }
                    break;
                case 2  :
                    //"Occurrence Event (omega event)";
                    if (final_state >= NumObservedLineages && final_state < (MAX_NUM_TOTAL_LINEAGES+1)) {
                        weight = DecayFactor*omega*temp_state;
                        Weights[final_state - 1] = Weights[final_state - 1] + weight;
                        weightsum = weightsum + weight;
                    }
                    break;
            }
        }
        PF_Log_Likelihood = PF_Log_Likelihood + log(weightsum) - log(NUMBER_OF_PARTICLES);
        //cout<<"log-likelihood after event "<<(J+1)<<" "<< PF_Log_Likelihood<<endl;
        //random sampling  of initial states...
        for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
            unif = Random(generator) * weightsum;
            double subsum = 0;
            for (int j = 0; j < MAX_NUM_TOTAL_LINEAGES; j++) {

                if (Weights[j] == 0) {
                    continue;
                } else {
                    subsum = subsum + Weights[j];
                    if (unif < subsum) {
                        InitialStates[i] = j;
                        break;
                    }
                }

            }
        }
        time = ObservedEventTimes[J];

        switch(ObservedEventCodes[J]) {
            case 0  :
                NumObservedLineages++;
                break;
            case 1:
                NumObservedLineages--;
                break;
        }

    }
    
    // last time-period
    for (int i = 0; i < MAX_NUM_TOTAL_LINEAGES; i++) {
        Weights[i] = 0;
    }
    weightsum = 0;
    
    for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
        DecayFactor = 1;
        final_state = GenerateTrajectory(time, InitialStates[i], NumObservedLineages);
        temp_state = (double) final_state;
        NumLineages = (double) NumObservedLineages;
        if (final_state >= NumObservedLineages && final_state<MAX_NUM_TOTAL_LINEAGES) {
            weight = DecayFactor*pow(rho,NumLineages)*pow(1 - rho, temp_state - NumLineages)*nchoosek(final_state,NumObservedLineages)*Myfactorial(NumObservedLineages);
            Weights[final_state] = Weights[final_state] + weight;
            weightsum = weightsum + weight;
            
        }
    }
    PF_Log_Likelihood = PF_Log_Likelihood + log(weightsum) - log(NUMBER_OF_PARTICLES);
    
    
    

    return PF_Log_Likelihood;
}



int main(int argc, char* argv[]){
    
	if(argc < 5){
	cout<<" Insufficient Number of Arguments "<<endl;
		return -1;
	}
    	beta = atof(argv[2]);
    	psi = atof(argv[3]);
    	mu = atof(argv[4]);
    	omega = atof(argv[5]);
    	rho = atof(argv[6]);

    string TreeInput = argv[7];
    string OutputFile= argv[8];

    SetStoichiometryMatrix();
    double ObservedEventTimes[MAX_NUM_OBSERVED_EVENTS] ={0};
    int  ObservedEventCodes[MAX_NUM_OBSERVED_EVENTS] = {0};
    int NumObservedEvents = 0;


    string line;
    ifstream ObservedTree(TreeInput);

    if (ObservedTree.is_open()) {

       // cout<<endl<<"Reading Observed Tree from file..."<<TreeInput<<endl<<endl;

        while ( getline (ObservedTree,line)) {
            if(line.length()==0){
                continue;
            }
            NumObservedEvents++;
            if(NumObservedEvents > 1){
                //cout << line << '\n';
                sscanf(line.c_str(),"%lf %d",&ObservedEventTimes[NumObservedEvents-2],&ObservedEventCodes[NumObservedEvents-2]);
            }
        }
        NumObservedEvents--;

        ObservedTree.close();

    }
    else{
        cout<<"Unable to open the tree file"<<endl;
        return -1;
    }
    /*
    cout<<endl<<"Number of Observed Events "<<NumObservedEvents<<endl;
    cout<<" The list of events is below "<<endl;
 
    
    for(int i = 0; i < NumObservedEvents;i++){
        cout<<"Event No: "<<(i+1)<<" Time: "<<ObservedEventTimes[i]<<" "<<" Event Type: ";

        switch(ObservedEventCodes[i]) {
            case 0  :
                cout<<"Branching (beta event)";
                break;
            case 1  :
                cout<<"Sampling with Removal Event (psi event)";
                break;
            case 2  :
                cout<<"Occurrence Event (omega event)";
                break;
            case 3  :
                cout<<"Time of origin";
                break;
            default:
                break;
        }
        cout<<endl<<endl;

    }
    */
 
    double EstimateLogLikelihood = EstimateLikelihood(NumObservedEvents,ObservedEventTimes,ObservedEventCodes);

  //cout<<" The estimated likelihood is "<<EstimateLogLikelihood<<endl;

   

    ofstream Output(OutputFile, std::ofstream::out | std::ofstream::app);
    cout.precision(10);
    Output<<omega<<" "<<EstimateLogLikelihood<<endl;
   	Output.close();
    
    return 0;
}
