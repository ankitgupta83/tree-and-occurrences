
#ifndef _BirthDeathExample_h
#define _BirthDeathExample_h

#define MAX_TIME_VALUE 10000  // proxy for infinity
#define NUMBER_OF_PARTICLES 100000
#define MAX_NUM_TOTAL_LINEAGES 100
#define MAX_NUM_OBSERVED_EVENTS 100
#define NUMBER_OF_SPECIES 1 
#define NUMBER_OF_REACTIONS 4
int StoichiometryMatrix[NUMBER_OF_SPECIES][NUMBER_OF_REACTIONS];
double propensity[NUMBER_OF_REACTIONS];




extern double beta, mu, psi, rho, omega;

// propensity functions

double CalculatePropensities(int * state){
    int n = state[0]; // total number of lineages
    
 
    propensity[0] = beta*n; // infection event..
    propensity[1] = mu*n;  // death without sampling
    propensity[2] = psi*n; // death with sampling
    propensity[3] = omega*n; // occurrence with removal

    double sum = 0;
    
    for(int i = 0; i < NUMBER_OF_REACTIONS ; i++){
        sum = sum + propensity[i];
    }
    return sum;
}


void SetStoichiometryMatrix(){
    for(int i = 0 ; i < NUMBER_OF_SPECIES ; i++){        
        for(int j = 0 ; j < NUMBER_OF_REACTIONS ; j++ ){            
            StoichiometryMatrix[i][j] = 0;
        }
    }
    // set the stoichiometry matrix. Only specify the non-zero entries   
    StoichiometryMatrix[0][0] = 1;
    StoichiometryMatrix[0][1] = -1;
    StoichiometryMatrix[0][2] = -1;
    StoichiometryMatrix[0][3] = -1;
}




#endif
