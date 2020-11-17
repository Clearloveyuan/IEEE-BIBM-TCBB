%Code by Dongyuan Li
%BIBM-19
clear;
clc;

options = [];
options.maxIter = 200;
options.error = 1e-4;
options.nRepeat = 30;
options.minIter = 50;
options.meanFitRatio = 0.1;
options.rounds = 30;
options.kmeans = 1;


%import my dataset��maybe it called AdjMtx; I have four layers and i will
%import four layers dataset including AdjMtx1 AdjMtx2 AdjMtx3 and AdjMtx4 .
AdjMtx1=importdata('H:\net1');
AdjMtx2=importdata('H:\net2');
AdjMtx3=importdata('H:\net3');
AdjMtx4=importdata('H:\net4');
K=271; 
S=importdata('H:\attri');

NumOfNodes = size(AdjMtx1,1);%AdjMtx1 has how many samples.
ClusteringList1=zeros(NumOfNodes,1);%clustering for node
U1=zeros(NumOfNodes,K);%U1 to store the result
%To get U1 and the first 10 times clustering result 
k=8;    
  A1 = createSPPMIMtx(AdjMtx1,k); 
        U_ = [];
        V_ = [];
        % to guarantee the property that U is nonnegative
        [U, V] = NMF(A1, K, options, U_, V_);        
        %for our method
        C = max(pinv(U)*S,0); 
        m = size(S, 2);

        %for tesing the numbers of negative samples, we set alpha = belta = 1
        %referencing 'Experiment-Parameter Sensitivity Analysis' section 
        belta = 1;
        alpha = 1;

        disp('Computing......');
        [Ucontent,C] = iterationNew(A1,S,U,C,belta,alpha);
        [Umax, Uindex] = max(Ucontent,[],2);
        U1 = Ucontent;
        ClusteringList1 = Uindex;

U2=zeros(NumOfNodes,K);
ClusteringList2=zeros(NumOfNodes,1);
    A2 = createSPPMIMtx(AdjMtx2,k); 
        U_ = [];
        V_ = [];
        % to guarantee the property that U is nonnegative
        [U, V] = NMF(A2, K, options, U_, V_);
        

        %for our method
        C = max(pinv(U)*S,0); 
        m = size(S, 2);

        %for tesing the numbers of negative samples, we set alpha = belta = 1
        %referencing 'Experiment-Parameter Sensitivity Analysis' section 
        belta = 1;
        alpha = 1;
        gama = 1;

        disp('Computing......');
        [Ucontent,C] = iteration(A2,S,U,C,U1,belta,alpha,gama);
        [Umax, Uindex] = max(Ucontent,[],2);
        ClusteringList2 = Uindex;
        U2 = Ucontent;



U3=zeros(NumOfNodes,K);
ClusteringList3 = zeros(NumOfNodes,1); 
    A3 = createSPPMIMtx(AdjMtx3,k); 
        U_ = [];
        V_ = [];
        % to guarantee the property that U is nonnegative
        [U, V] = NMF(A3, K, options, U_, V_);
        

        %for our method
        C = max(pinv(U)*S,0); 
        m = size(S, 2);

        %for tesing the numbers of negative samples, we set alpha = belta = 1
        %referencing 'Experiment-Parameter Sensitivity Analysis' section 
        belta = 1;
        alpha = 1;
        gama = 1;

        disp('Computing......');
        [Ucontent,C] = iteration(A3,S,U,C,U2,belta,alpha,gama);
        [Umax, Uindex] = max(Ucontent,[],2);
        ClusteringList3 = Uindex;
        U3 = Ucontent;



U4 = zeros(NumOfNodes,K);
ClusteringList4 = zeros(NumOfNodes,1);
% parfor k = 1:10    
    A4 = createSPPMIMtx(AdjMtx4,k); 
        U_ = [];
        V_ = [];
        % to guarantee the property that U is nonnegative
        [U, V] = NMF(A4, K, options, U_, V_);
        

        %for our method
        C = max(pinv(U)*S,0); 
        m = size(S, 2);

        %for tesing the numbers of negative samples, we set alpha = belta = 1
        %referencing 'Experiment-Parameter Sensitivity Analysis' section 
        belta = 1;
        alpha = 1;
        gama = 1;

        disp('Computing......');
        [Ucontent,C] = iteration(A4,S,U,C,U3,belta,alpha,gama);
        [Umax, Uindex] = max(Ucontent,[],2);
        ClusteringList4 = Uindex;
        U4 = Ucontent;


