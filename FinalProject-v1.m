%% Initialization
clearvars;
clc;

%% Loading the relevant files
% gfunda, gsecd, gexrtdly

% Read & Screen gfunda data
gfundaFile = "/Users/zacharyzhu/Desktop/Interns/Round 2/John's Data/comp/Global/gfunda/gfunda_20190611.csv";
opts = detectImportOptions(gfundaFile);
openvar opts.VariableNames
opts.SelectedVariableNames = {'gvkey','datadate','curcd','indfmt','seq','txditc','pstk'}; 
gfundaData = readtable(gfundaFile, opts);
tokeep = strcmp(gfundaData.indfmt,'INDL');
gfundaData = gfundaData(tokeep, :);

gfundaData.gvkey(isnan(gfundaData.gvkey)) = -9999;
gfundaData.datadate(isnan(gfundaData.datadate)) = -9999;
gfundaData.seq(isnan(gfundaData.seq)) = -9999;
gfundaData.txditc(isnan(gfundaData.txditc)) = -9999;
gfundaData.pstk(isnan(gfundaData.pstk)) = -9999;
gfundaData = unique(gfundaData, 'rows');
gfundaData.gvkey(gfundaData.gvkey == -9999) = nan;
gfundaData.datadate(gfundaData.datadate == -9999) = nan;
gfundaData.seq(gfundaData.seq == -9999) = nan;
gfundaData.txditc(gfundaData.txditc == -9999) = nan;
gfundaData.seq(gfundaData.seq == -9999) = nan;

gsecdFile = "/Users/zacharyzhu/Desktop/Interns/Round 2/John's Data/comp/Global/gsecd/gsecd_20190611.csv";
gsecdInfo = detectImportOptions(gsecdFile);
gsecdData = readtable(gsecdFile, gsecdInfo);
gexrtdlyFile = readtable("gexrtdly_20190611.txt");
