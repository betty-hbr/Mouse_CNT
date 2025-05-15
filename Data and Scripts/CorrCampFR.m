%%%
%%% This script is to evaluate the correlation between C amplitude and
%%% neural average firing rate.
%%%
%%% Input: "mouse_XX_anald.mat"
%%%        "mouse_XX_syn.mat"
%%%
%%% Output: "FiringRate_mean" the trial-averaged mean firing rate of each
%%%         locus of each subject.
%%%         "Coef_mean" the trial-averaged mean coefficient amplitude of
%%%         each locus of each subject.
%%%
%%% 
%%% Create by Borong 2022/Oct @CUHK

clc
clear
close all

subjects = [{'17_10'},{'19_12'},{'20_13'},...
                                {'21_14'},{'22_15'},{'23_16'},{'24_17'}];
stim_con = 'SS';
inten_con = 'T';

for s = 1:length(subjects)
    ind = str2double(subjects{s}(4:5));
    filename_anald = ['mouse_',subjects{s},'_anald_',stim_con,'_',inten_con,'.mat'];
    filename_syn = ['mouse_',subjects{s},'_syn_',stim_con,'_',inten_con,'.mat'];
    if isfile(filename_anald) && isfile(filename_syn)
        mat_anald = matfile(filename_anald);
        anald = mat_anald.anald;
        mat_syn = matfile(filename_syn);
        syn = mat_syn.syn;
        clear mat_anald mat_syn
        
        % -----------------------------------------------------------------
        %  EXTRACT SPIKES FROM RECORDINGS OF EACH INDIVIDUAL LOCUS (BLOCKS)
        %  ----------------------------------------------------------------

        % Set up variables
        u = 1; % specify which unit to be considered (single-unit analysis), 
        % "u = 1" means considering only the first unit.
        Location = anald.blocks;
        fs = anald.fs;
        if isfield(anald,'LClr')
            TimePoints = 1:length(anald.LClr);
        elseif isfield(anald,'LC1r')
            TimePoints = 1:length(anald.LC1r);
        end
        TrialOnset = syn.C_temp.sp;
        StimOnset = syn.C_temp.sp+100;
        StimOffset = syn.C_temp.ep-100;
        TrialOffset = syn.C_temp.ep;
        Unit = anald.units{u};
        
        % Calling a constructor of `SpikeTrainClass`
        % Inside the 'SpikeTrainClass' constructor:
            % 'SpikeTrainClass.Spikes' dataset with dimension of (Num of blocks, 
            % Spike train for each trial 400ms, 5 trials that were run for each 
            % block).
            % 'SpikeTrainClass.FiringRate' trial-averaged firing rate of all the 
            % blocks (PSTH).
        MySpikeData = SpikeTrainClass(filename_anald,Location,TimePoints,fs,TrialOnset,...
                        StimOnset,StimOffset,TrialOffset,Unit);                                      
        clear u Unit filename_anald
        % Get trial-averaged firing rate for each block.
        NumTrials = 5; % number of trials performed at each locus
        FiringRate_mean{ind,1} = zeros(length(MySpikeData.Location),1);
        for b = 1:length(MySpikeData.Location)
            fr_temp = sum(MySpikeData.Spikes(b,:,:),3)/NumTrials;
            FiringRate_mean{ind,1}(b) = mean(smooth(1000*fr_temp,5));
            clear fr_temp
        end
        clear b
        
        %  ----------------------------------------------------------------
        %  GET THE TRIAL-AVERAGED C AMPLITUDE FOR EACH BLOCK.
        %  ----------------------------------------------------------------
        
        % Calling a constructor of "CoefficientClass"
        % Inside the "CoefficientClass" constructor:
            % "CoefficientTrain"
            % "Coefficients" is a cell structure with the length of number
            % of coefficients. "Coefficients{n}" dataset with dimension of 
            % (Num of blocks, coefficient train for each trial 400ms, 5 
            % trials that were run for each block).
            % 'CoefficientsAve' is the trial-averaged coefficients of all 
            % the blocks.
        MyCoefficientData = CoefficientClass(filename_syn,Location,TimePoints,...
                            fs,TrialOnset,StimOnset,StimOffset,TrialOffset,syn); 
        clear filename_syn Location TimePoints fs TrialOnset StimOnset StimOffset TrialOffset syn
        
        for c = 1:length(MyCoefficientData.Coefficients)
            for b = 1:length(MyCoefficientData.Location)
                Coef_temp = sum(MyCoefficientData.Coefficients{c}(b,:,:),3)/NumTrials;
                Coef_mean{ind,1}{c,1}(b,1) = mean(smooth(Coef_temp,5));
                Coef_max{ind,1}{c,1}(b,1) = max(smooth(Coef_temp,5));
            end
            clear b Coef_temp
        end
        clear c
    end
    clear ind anald
end
clear s

%% ------------------------------------------------------------------------
%  STACK THE MEAN FIRING RATE AND SELECTED MEAN COEFFICIENT AMPLITUDE
%  TOGETHER.
%   If the "flag_all = 0" and "flag_preference = 0":
%   The C would be selected by the highest correlation coefficient with
%   averaged neural firing rate, denoted as "C_preferred". 
%   If the "flag_all = 0" and "flag_preference = 1":
%   The corresponding C would be selected according to the W preference 
%   pair then correlate with the averaged neural firing rate. 
%  ------------------------------------------------------------------------
flag_all = 0;
flag_preference = 1;
TimeSeries_mean_total = [];
C_preference_id = [0;0;0;0;0;0;0;0;0;1;0;4;3;4;2;3;4]; % corresponding to the preference index 
% of synergies.
for s = 1:length(FiringRate_mean)
    if ~isempty(FiringRate_mean{s}) && flag_all
        for c = 1:length(Coef_mean{s}) % interate all the Cs
            % Stack the firing rate and all the C mean amplitudes
            TimeSeries_mean_total = [TimeSeries_mean_total ; [FiringRate_mean{s},Coef_max{s}{c}]];
        end
        clear c
    elseif ~isempty(FiringRate_mean{s}) && ~flag_all && ~flag_preference
        for c = 1:length(Coef_mean{s})
            Corr_temp = corrcoef(FiringRate_mean{s},Coef_mean{s}{c});
            Corr(c) = Corr_temp(1,2);
        end
        clear c Corr_temp
        % Stack the firing rate and the C mean amplitudes with the highset
        % correlation coefficient value.
        [~,C_id] = max(Corr);
        TimeSeries_mean_total = [TimeSeries_mean_total ; [FiringRate_mean{s},Coef_max{s}{C_id}]];
        clear Corr
    elseif ~isempty(FiringRate_mean{s}) && ~flag_all && flag_preference
        TimeSeries_mean_total = [TimeSeries_mean_total ; [FiringRate_mean{s},Coef_max{s}{C_preference_id(s)}]];
    end
end
clear s

x = TimeSeries_mean_total(:,1); % mean neural firing rate
y = TimeSeries_mean_total(:,2); % mean coefficient amplitude
[fitresult, gof] = createFit(x, y);


%% ------------------------------------------------------------------------
%  An example of mean firing rate and activation coefficient amplitude
%  ------------------------------------------------------------------------

% From "mouse_19_12"

figure,
subplot(5,1,1), 
plot(smooth(FiringRate_mean{12}),'r-o','LineWidth',1.2); 
xticklabels([]); 
yticklabels([]);
ylabel('F.R.'); 
grid off
set(gca,'FontSize',20);

subplot(5,1,2), 
plot(smooth(Coef_mean{12}{1}),'b-o','LineWidth',1.2); 
xticklabels([]); 
yticklabels([]);
ylabel('C1'); 
yyaxis right 
yticks([])
ylabel(num2str(round(corr(FiringRate_mean{12},Coef_mean{12}{1}),2)));
grid off
set(gca,'FontSize',20);

subplot(5,1,3), 
plot(smooth(Coef_mean{12}{2}),'b-o','LineWidth',1.2); 
xticklabels([]); 
yticklabels([]);
ylabel('C2'); 
yyaxis right 
yticks([])
ylabel(num2str(round(corr(FiringRate_mean{12},Coef_mean{12}{2}),2)));
grid off
set(gca,'FontSize',20);

subplot(5,1,4), 
plot(smooth(Coef_mean{12}{3}),'b-o','LineWidth',1.2); 
xticklabels([]); 
yticklabels([]);
ylabel('C3'); 
yyaxis right 
yticks([])
ylabel(num2str(round(corr(FiringRate_mean{12},Coef_mean{12}{3}),2)));
grid off
set(gca,'FontSize',20);

subplot(5,1,5), 
plot(smooth(Coef_mean{12}{4}),'b-o','LineWidth',1.2); 
xticklabels([]); 
yticklabels([]);
ylabel('C4'); 
yyaxis right 
yticks([])
ylabel(num2str(round(corr(FiringRate_mean{12},Coef_mean{12}{4}),2)));
grid off
set(gca,'FontSize',20);

xlabel('Rostral <--> Caudal');
set(gca,'FontSize',20);
set(gcf,'Position',[400 400 500 700]);



function [fitresult, gof] = createFit(x, y)
    %CREATEFIT(X,Y)
    %  Create a fit.
    %
    %  Data for 'untitled fit 1' fit:
    %      X Input : x
    %      Y Output: y
    %  Output:
    %      fitresult : a fit object representing the fit.
    %      gof : structure with goodness-of fit info.
    %
   
    [xData, yData] = prepareCurveData( x, y );
    
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    
    % Plot fit with data.
    figure( 'Name', 'Linear fit' );
    h = plot( fitresult, xData, yData);
    h(2).LineWidth = 1.2;
    h(1).MarkerSize = 10;
    legend( h, 'Coefficient amplitude vs. Neural firing rate', 'Linear fit',...
                        'Location', 'NorthEast', 'Interpreter', 'none' );
    text(25, 0.006, ['Rsquare is ', num2str(gof.rsquare)],'FontSize',15);
    % Label axes
    xlabel( 'Neural firing rate', 'Interpreter', 'none' );
    ylabel( 'Coefficient amplitude', 'Interpreter', 'none' );
    set(gca, 'FontSize', 20);
    set(gcf, 'Position', [400 400 800 400])
    grid on
end
















