
clear all,clc,close all;
% 
maindir = '/Users/ece.kaya/MATLAB/Githubbed/Misc/HWE/jatos_resultfiles_HWE';
maindir_contents = dir(maindir);

maindir_contents = maindir_contents(3:end); % for macos users. check unused rows & change if needed

df = table();
mITIs = [];

for subj = 1:length(maindir_contents)
    
    currdir = maindir_contents(subj).name;
    subfile = [maindir filesep currdir filesep 'data.csv'];

    data = importfile(subfile);

    dfs = table();

    trial = 1;
    for i = 1:4:height(data)-4

        stim1 = str2double(extractBefore(data{i+1,1},'_ms.wav'));
        stim2 = str2double(extractBefore(data{i+2,1},'_ms.wav'));

        tap_times = str2double(split(data{i+3,4}, ','));

        mITI = median(rmoutliers(diff(tap_times)));

        dfs.tNr(trial) = trial;
        dfs.stim1(trial) = stim1;
        dfs.stim2(trial) = stim2;
        dfs.miti(trial) = mITI;
        trial = trial +1;
    end
    dfs = sortrows(dfs,'stim1','ascend');
    mITIs = horzcat(mITIs, dfs{:,4});

end

df = dfs(:,2:3);
df.mitis = mITIs;

%%
figure;
k = 1;
checklist = zeros(height(df),1);

for i = 1:height(df)
    
    if checklist(i) == 0
        
    s_val = df.stim1(i,1);              %standard 
    pair_sc = double(df.mitis(i,:));

    idx_cs = df.stim2==s_val & df.stim1==df.stim2(i,1);
    c_val = df.stim1(idx_cs,:);

    pair_cs = double(df.mitis(idx_cs,:));

    subplot(2,5,k)
    scatter(pair_sc, pair_cs, 100, 'filled'); hold on;
    xlim([300 900]); ylim([300 900]);
    xline(s_val, '--', 'LineWidth', 3, 'Color', [0.4667    0.6745    0.1882]);
    xline(c_val, ':', 'LineWidth', 3, 'Color', [0.4667    0.6745    0.1882]);
    yline(c_val, '--', 'LineWidth', 3, 'Color', [0.4667    0.6745    0.1882]);
    yline(s_val, ':', 'LineWidth', 3, 'Color', [0.4667    0.6745    0.1882]);
    
    geo_mean = geomean([s_val, c_val]);
    arit_mean = mean([s_val, c_val]);
    
%     xline(arit_mean, '-', 'LineWidth', 2, 'Color', [0.9294    0.6941    0.1255]); %arit-mean
%     yline(arit_mean, '-', 'LineWidth', 2, 'Color', [0.9294    0.6941    0.1255]); %arit-mean
%     
%     xline(geo_mean, '-', 'LineWidth', 2, 'Color', [0.8510    0.3255    0.0980]); %-geo_mean
%     yline(geo_mean, '-', 'LineWidth', 2, 'Color', [0.8510    0.3255    0.0980]); % geo_mean  
    
    scatter(mean(pair_sc), mean(pair_cs), 200, 'filled', 'MarkerFaceAlpha',.5);
%     tittxt = ['mean result=' num2str(mean([mean(pair_sc),mean(pair_cs)])) 'ms'];
%     EK_plotlabels([num2str(s_val) 'ms & ' num2str(c_val) 'ms'], ...
%         [num2str(c_val) 'ms & ' num2str(s_val) 'ms'],tittxt,13);
 
    tittxt = [num2str(s_val) 'ms & ' num2str(c_val) 'ms --> ' num2str(abs(c_val-s_val)) 'ms diff'];
    EK_plotlabels('st-comp', 'comp-st', tittxt, 13);    

    k = k+1;

    checklist(idx_cs) = 1;
    end

    
end
legend('indiv. data','stim. heard first', 'stim. heard second');
% legend('','stim. heard first', 'stim. heard second','','','ar. mean', 'geo. mean','particip. mean');

%% 
function data = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DATA = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  DATA = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  data = importfile("/Users/ece.kaya/MATLAB/Githubbed/Misc/HWE/jatos_resultfiles_HWE/comp-result_6834/data.csv", [2, 82]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 03-May-2021 11:27:16

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "stimulus", "Var3", "Var4", "trial_index", "time_elapsed", "Var7", "Var8", "tap_times"];
opts.SelectedVariableNames = ["stimulus", "trial_index", "time_elapsed", "tap_times"];
opts.VariableTypes = ["string", "string", "string", "string", "double", "double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "stimulus", "Var3", "Var4", "Var7", "Var8", "tap_times"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "stimulus", "Var3", "Var4", "Var7", "Var8", "tap_times"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["trial_index", "time_elapsed"], "ThousandsSeparator", ",");

% Import the data
data = readtable(filename, opts);

end