% clear screen and variables
clc
clearvars

% enter cpp source file for Particle Filter based likelihood estimation
Cpp_Source = './ParticleFilter.cpp';
Cpp_Binary = './ParticleFilter.o';


% enter input and output files
Tree_input = './ExampleTrees/Tree4.txt'; % this is the file from which the tree will be read
Output_File = './OutputSimulationBasedLikelihoods.txt'; % this is a temporary file where estimated likelihoods will be written




% set the parameter values 
beta =1.25;
mu = 0;
psi = 0.1;
rho = 0.7;
rhoc = 1 - rho;
omega = 0.1;


% read the input file
delimiterIn = ' ';
headerlinesIn = 2;
ObservedTreeData = importdata(Tree_input,delimiterIn,headerlinesIn);
ObservedEventTimes = ObservedTreeData.data(:,1);
ObservedEventCodes = ObservedTreeData.data(:,2);
clear ObservedTreeData;




% compute total number of events
NumEvents = length(ObservedEventTimes);
    

% set the number of observed lineages at t = t_or
NumObservedLineages = 1;
   
% compute number of observed lineages at the present t = 0
FinalNumObservedLineages = NumObservedLineages;
for J=(NumEvents-1):-1:1
    switch ObservedEventCodes(J)
        case 0
        %Branching Event
            FinalNumObservedLineages = FinalNumObservedLineages + 1;
        case 1
            %Sampling with removal;
            FinalNumObservedLineages = FinalNumObservedLineages - 1;
        otherwise
            % observed lineages do not change
    end
end
    
    




LogLikelihoodValues = zeros(20,1);

 % we now estimate with simulations   
 str1 = sprintf('g++ %s -o %s',Cpp_Source,Cpp_Binary);
 system(str1);

for inum = 1:1:20
    beta = inum*0.1;
    %omega = inum*0.05;
    % first compute using simulations


    fprintf('\n Estimating likelihood with parameters beta = %d, psi = %d, omega = %d, mu =%d and rho = %d\n',beta,psi,omega,mu,rho);
   
    str2 = sprintf('%s %i %d %d %d %d %d %s %s',Cpp_Binary,5,beta,psi,mu,omega,rho,Tree_input,Output_File);
    system(str2);
    disp('Particle Filter Estimation completed!');

    
    % Numerical estimation starts
    
    NumObservedLineages = FinalNumObservedLineages;
    p0 = rhoc;
    v = 1;
    logfactor = 0;
    for J=1:1:NumEvents
        
        oldp0 = p0;
        p0 = GetP0(ObservedEventTimes(J),beta,rhoc,mu,psi+omega,0);
        if J > 1
            deltat = ObservedEventTimes(J)-ObservedEventTimes(J-1);
        else
            deltat = ObservedEventTimes(J);
        end
        
        v = TransformDerivativeContrVec(deltat,beta,oldp0,mu,psi+omega,0,NumObservedLineages,v);
        switch ObservedEventCodes(J)
            case 0
                % Branching Event; 
               logfactor = logfactor +log(beta);
               NumObservedLineages = NumObservedLineages - 1;
               
            case 1
                % Sampling with removal
                logfactor = logfactor +log(psi);
                NumObservedLineages = NumObservedLineages + 1;
                
            case 2
                % Occurrence Event   
                logfactor = logfactor +log(omega);
                v= [0;v];
            otherwise
        end
    end
    
    LogLikelihoodValues(inum) = FinalNumObservedLineages*log(rho)+logfactor+log(v(1));
    disp('Analytical Estimation completed!');
    
end



%read estimated likelihoods from the file
SimLik = importdata(Output_File);

%clear unnecessary files
system(sprintf('rm %s',Output_File));
system(sprintf('rm %s',Cpp_Binary));

figure('Name','log-probability');
plot(0.05:0.05:1,LogLikelihoodValues,'-','LineWidth',1.5,'Color','red'); hold on
plot(0.05:0.05:1,SimLik(:,2),'*','LineWidth',1.5,'Color','blue'); hold on
ylabel('log-probability','FontSize',15); hold on
xlabel('\beta','FontSize',15); hold off
legend('Analytical Likelihood', 'Simulated Likelihood');


% Below are the helper functions we need for analytical computation

function q = GetQ(t,beta,rhoc,mu,psi,omega)

    c1 = sqrt( (beta - mu - omega - psi)*(beta - mu - omega - psi)+4*beta*psi);
    c2 = (beta + mu + omega + psi - 2*beta*rhoc)/c1;
    q = 2*(1-c2*c2) + exp(-c1*t)*(1-c2)*(1-c2) + exp(c1*t)*(1+c2)*(1+c2);
end


function q = GetDerivativeQ(t,beta,rhoc,mu,psi,omega,n)
    c1 = sqrt( (beta - mu - omega - psi)*(beta - mu - omega - psi)+4*beta*psi);
    c2 = (beta + mu + omega + psi - 2*beta*rhoc)/c1;
    if(n == 1)    
        q = -(4*beta/c1)*( exp(c1*t) -exp(-c1*t) + c2*( exp(c1*t) + exp(-c1*t) - 2 )    );
    elseif ( n==2 )
        q = 8*( (beta/c1)^2)*(exp(c1*t) + exp(-c1*t) - 2);
    else
        q = 0;
    end
end



function q = GetMultiDerivativeRecQ(t,beta,rhoc,mu,psi,omega,n,NumObservedLineages)
    q0 = 1/GetQ(t,beta,rhoc,mu,psi,omega);
    if (n == 0)
        q = (q0)^NumObservedLineages;
    elseif (NumObservedLineages == 0)
        q = 0;
    else
        q = 0;
        Qder1 = GetDerivativeQ(t,beta,rhoc,mu,psi,omega,1);
        Qder2 = GetDerivativeQ(t,beta,rhoc,mu,psi,omega,2);  
        for m1=0:1:n
            if(mod(n-m1,2) == 0)
                m2 = (n-m1)/2;
                q = q + (-1)^(m1+m2)*(2^(-m2))*(q0^(m1 +m2+NumObservedLineages))*nchoosek(m1+m2+NumObservedLineages-1,m1+m2)*nchoosek(m1+m2,m1)*( Qder1^m1 )*(Qder2^m2);
            end
        end
    end
    q = q*(GetQ(0,beta,rhoc,mu,psi,omega)^NumObservedLineages)*factorial(n);
end

function p0 = GetP0(t,beta,rhoc,mu,psi,omega)
    c1 = sqrt( (beta - mu - omega - psi)*(beta - mu - omega - psi)+4*beta*psi);
    c2 = (beta + mu + omega + psi - 2*beta*rhoc)/c1;
    factor = ( exp(c1*t)*(1+c2)^2 - exp(-c1*t)*(1-c2)^2)/GetQ(t,beta,rhoc,mu,psi,omega);
    p0 = (1/(2*beta))*( beta + mu +omega + psi)  - (c1/(2*beta))*factor;

end


function derp0 = GetDerP0(t,beta,rhoc,mu,psi,omega,n)
    c1 = sqrt( (beta - mu - omega - psi)*(beta - mu - omega - psi)+4*beta*psi);
    c2 = (beta + mu + omega + psi - 2*beta*rhoc)/c1;

    derfac = -2*beta/c1;
    factor = zeros(n+1,1);
    factor(1) = ( exp(c1*t)*(1+c2)^2 - exp(-c1*t)*(1-c2)^2)/4; 
    factor(2) =  derfac*( exp(c1*t)*(1+c2) + exp(-c1*t)*(1-c2))/2;
    factor(3) = (derfac^2)*( exp(c1*t) - exp(-c1*t))/2;
    derp0 = 0;
    for i = 0:1:n
        derp0 = derp0 + nchoosek(n,i)*factor(i+1)*GetMultiDerivativeRecQ(t,beta,rhoc,mu,psi,omega,n-i,1);
    end
    derp0 = -derp0*(c1/(2*beta));
end

function v1 = TransformDerivativeContrVec(t,beta,rhoc,mu,psi,omega,NumObservedLineages,v)
    n = length(v);
    if(n > 1)
        %requires computation of (n-1) derivatives%
        derivativeVector = zeros(n-1,1);
        for i=1:1:(n-1)
            derivativeVector(i) = GetDerP0(t,beta,rhoc,mu,psi,omega,i);
        end
        B = IncompleteBellPolynomial(n-1,n-1,derivativeVector);
       A = zeros(n,n);

        for i=1:1:n
            for j=1:1:i
                A(i,j) = nchoosek(i-1,j-1)*GetMultiDerivativeRecQ(t,beta,rhoc,mu,psi,omega,i-j,NumObservedLineages);
            end
        end
        C = A*B;
        v1 = (C')*v;
    else
        v1 = GetMultiDerivativeRecQ(t,beta,rhoc,mu,psi,omega,0,NumObservedLineages)*v;
    end
end


function OutputMatrix = IncompleteBellPolynomial(N,K,Vector)
A = zeros(N,K);
A(1:N,1) = Vector(1:N);
 for n=2:N 
  for k=2:K 
   for m=1:n-k+1 %recursive relation
    A(n,k) =A(n,k) + nchoosek(n-1,m-1)*Vector(m)*A(n-m,k-1); 
   end
  end
 end
OutputMatrix = eye(N+1,K+1);
OutputMatrix(2:N+1,2:K+1) =A;
end 


