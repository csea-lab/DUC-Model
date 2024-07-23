%% fits for the grand mean
% compmodelfit
clear 
clc
close all
load('workspace.mat')
taxis = 1:6000;

RDKtemp = rdkmat(3, 201:12200).*10^7;
RSVPtemp = RSVPmat(3,201:12200).*10^7;

RDKtemp = bslcorr(RDKtemp, 1:6000);
RSVPtemp = bslcorr(RSVPtemp, 1:6000);

Input_empirical = [RSVPtemp RDKtemp];

% first, fit the DUC model
 [BETA_DUC,R,J,COVB,MSE_DUC] = nlinfit(taxis,Input_empirical,@DUCmodel,[1 1 2 -.5 2]);
 BETA_DUC
subplot(2,1,1)
 [ssveptotal_DUC] = DUCmodel(BETA_DUC,taxis); 
 plot(ssveptotal_DUC(2001:12000), 'LineWidth',2)
 hold on
 plot(ssveptotal_DUC(14001:end), 'LineWidth',2)
 plot(Input_empirical(2001:12000),  'LineWidth',2)
  plot(Input_empirical(14001:end),  'LineWidth',2)

  % Now the equivalent biased competition model
 [BETA_BIAS,R,J,COVB,MSE_BIAS] = nlinfit(taxis,Input_empirical,@BiasedCompmodel,[1 1 1 1 .2]);
 BETA_BIAS
subplot(2,1,2)
 [ssveptotal_BIAS] = BiasedCompmodel(BETA_BIAS,taxis); 
 plot(ssveptotal_BIAS(2001:12000),  'LineWidth',2)
 hold on
 plot(ssveptotal_BIAS(14001:end), 'LineWidth',2)
 plot(Input_empirical(2001:12000),  'LineWidth',2)
  plot(Input_empirical(14001:end),  'LineWidth',2)


  %% this cell computes fits based on bootstrapped means and resultin stats fro model fit comparison
  cd '/Users/andreaskeil/Dropbox (UFL)/ak_own/Ducmodel'
  load('RSVP_RDK_bysubjects.mat')
  taxis = 1:6000;

warning('off')
for bootstrapInd = 1:5000
    % make a new grand mean for each stimulus
    subjectvec = randi(36, 36, 1); 
    RDKtempU = mean(RDK_unpleasants(subjectvec, 201:12200).*10^7); 
    RSVPtempU = mean(RSVP_unpleasants(subjectvec, 201:12200).*10^7); 

    RDKtempP = mean(RDK_pleasants(subjectvec, 201:12200).*10^7); 
    RSVPtempP = mean(RSVP_pleasants(subjectvec, 201:12200).*10^7); 

    RDKtempN = mean(RDK_neutrals(subjectvec, 201:12200).*10^7); 
    RSVPtempN = mean(RSVP_neutrals(subjectvec, 201:12200).*10^7); 
    
    % baseline correction of bootstrap means
    RDKtempU = bslcorr(RDKtempU, 1:6000);
    RSVPtempU = bslcorr(RSVPtempU, 1:6000);

    RDKtempP = bslcorr(RDKtempP, 1:6000);
    RSVPtempP = bslcorr(RSVPtempP, 1:6000);

     RDKtempN = bslcorr(RDKtempN, 1:6000);
    RSVPtempN = bslcorr(RSVPtempN, 1:6000);
    
    % concatenate so that they serve as one input to nlinfit
    Input_empiricalU = [RSVPtempU RDKtempU];
    Input_empiricalP = [RSVPtempP RDKtempP];
    Input_empiricalN= [RSVPtempN RDKtempN];

     [BETA_DUCU(bootstrapInd, :),~,~,~,MSE_DUCU(bootstrapInd)] = nlinfit(taxis,Input_empiricalU,@DUCmodel,[1 1 3 -.5 2]);
     [BETA_BIASU(bootstrapInd, :),~,~,~,MSE_BIASU(bootstrapInd)] = nlinfit(taxis,Input_empiricalU,@BiasedCompmodel,[1 1 1 1 1]);

      [BETA_DUCP(bootstrapInd, :),~,~,~,MSE_DUCP(bootstrapInd)] = nlinfit(taxis,Input_empiricalP,@DUCmodel,[1 1 3 -.5 2]);
     [BETA_BIASP(bootstrapInd, :),~,~,~,MSE_BIASP(bootstrapInd)] = nlinfit(taxis,Input_empiricalP,@BiasedCompmodel,[1 1 1 1 1]);
      
     [BETA_DUCN(bootstrapInd, :),~,~,~,MSE_DUCN(bootstrapInd)] = nlinfit(taxis,Input_empiricalN,@DUCmodel,[1 1 3 -.5 2]);
     [BETA_BIASN(bootstrapInd, :),~,~,~,MSE_BIASN(bootstrapInd)] = nlinfit(taxis,Input_empiricalN,@BiasedCompmodel,[1 1 1 1 1]);

     if bootstrapInd/100 == round(bootstrapInd./100), disp(bootstrapInd), end

end

[BFU] = bootstrap2BF_z(MSE_BIASU,MSE_DUCU, 1);
log10BFU = log10(BFU);
[BFP] = bootstrap2BF_z(MSE_BIASP,MSE_DUCP, 1);
log10BFP = log10(BFP);
[BFN] = bootstrap2BF_z(MSE_BIASN,MSE_DUCN, 1);
log10BFN = log10(BFN);

[BF_UN] = bootstrap2BF_z(MSE_DUCN,MSE_DUCU, 1);
log10BFUN = log10(BF_UN);
[BF_UP] = bootstrap2BF_z(MSE_DUCP,MSE_DUCU, 1);
log10BFUP = log10(BF_UP);
[BF_NP] = bootstrap2BF_z(MSE_DUCP,MSE_DUCN, 1);
log10BFNP = log10(BF_NP);

[BF_para3DUCUN] = bootstrap2BF_z(BETA_DUCU(:, 3) ,BETA_DUCN(:,3), 1);
[BF_para3DUCUP] = bootstrap2BF_z(BETA_DUCU(:, 3) ,BETA_DUCP(:,3), 1);
[BF_para3DUCNP] = bootstrap2BF_z(BETA_DUCN(:, 3) ,BETA_DUCP(:,3), 1);

[BF_para4DUCUN] = bootstrap2BF_z(BETA_DUCU(:, 4) ,BETA_DUCN(:,4), 1);
[BF_para4DUCUP] = bootstrap2BF_z(BETA_DUCU(:, 4) ,BETA_DUCP(:,4), 1);
[BF_para4DUCNP] = bootstrap2BF_z(BETA_DUCN(:, 4) ,BETA_DUCP(:,4), 1);

[BF_para5DUCUN] = bootstrap2BF_z(BETA_DUCU(:, 5) ,BETA_DUCN(:,5), 1);
[BF_para5DUCUP] = bootstrap2BF_z(BETA_DUCU(:, 5) ,BETA_DUCP(:,5), 1);
[BF_para5DUCNP] = bootstrap2BF_z(BETA_DUCN(:, 5) ,BETA_DUCP(:,5), 1);

