%% Initialization
clearvars;
clc;

% gfunda, gsecd, gexrtdly

%% Read and screen the gfunda data
gfundaFile = "/Users/zacharyzhu/Desktop/Interns/Round 2/John's Data/comp/Global/gfunda/gfunda_20190611.csv";
opts = detectImportOptions(gfundaFile);
opts.SelectedVariableNames = {'gvkey','datadate','curcd','indfmt','seq','txditc','pstk','fic','loc'}; 
gfundaData = readtable(gfundaFile,opts);
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

%% Divide gfunda data into regions
gfundaData.region = zeros(height(gfundaData),1);
tochange = strcmp(gfundaData.loc,'USA') | strcmp(gfundaData.loc,'CAN') | strcmp(gfundaData.fic,'USA') | strcmp(gfundaData.fic,'CAN');
gfundaData.region(tochange) = 1;
tochange = strcmp(gfundaData.loc,'JAP') | strcmp(gfundaData.fic,'JAP');
gfundaData.region(tochange) = 2;
tochange = strcmp(gfundaData.loc, 'AUS') | strcmp(gfundaData.loc, 'NZL') | strcmp(gfundaData.loc, 'HKG')| strcmp(gfundaData.loc,'SGP')...
    | strcmp(gfundaData.fic, 'AUS') | strcmp(gfundaData.fic, 'NZL') | strcmp(gfundaData.fic, 'HKG') | strcmp(gfundaData.fic,'SGP');
gfundaData.region(tochange) = 3;
tochange =  strcmp(gfundaData.loc, 'AUT') | strcmp(gfundaData.loc, 'BEL') | strcmp(gfundaData.loc, 'DNK') | strcmp(gfundaData.loc, 'FIN') | strcmp(gfundaData.loc, 'FRA') |...
    strcmp(gfundaData.loc, 'DEU') | strcmp(gfundaData.loc, 'GRC') | strcmp(gfundaData.loc, 'IRL') | strcmp(gfundaData.loc, 'ITA') | strcmp(gfundaData.loc, 'NLD') |...
    strcmp(gfundaData.loc, 'NOR') | strcmp(gfundaData.loc, 'PRT') | strcmp(gfundaData.loc, 'ESP') | strcmp(gfundaData.loc, 'SWE') | strcmp(gfundaData.loc, 'CHE') | strcmp(gfundaData.loc, 'GBR')|...
    strcmp(gfundaData.fic, 'AUT') | strcmp(gfundaData.fic, 'BEL') | strcmp(gfundaData.fic, 'DNK') | strcmp(gfundaData.fic, 'FIN') | strcmp(gfundaData.fic, 'FRA') |...
    strcmp(gfundaData.fic, 'DEU') | strcmp(gfundaData.fic, 'GRC') | strcmp(gfundaData.fic, 'IRL') | strcmp(gfundaData.fic, 'ITA') | strcmp(gfundaData.fic, 'NLD') |...
    strcmp(gfundaData.fic, 'NOR') | strcmp(gfundaData.fic, 'PRT') | strcmp(gfundaData.fic, 'ESP') | strcmp(gfundaData.fic, 'SWE') | strcmp(gfundaData.fic, 'CHE') | strcmp(gfundaData.fic, 'GBR');
gfundaData.region(tochange) = 4;
tokeep = gfundaData.region ~= 0;
gfundaData = gfundaData(tokeep, :);


%% Load & screen the gsecd data --figure out how to load greater quantities
gsecdFile = "/Users/zacharyzhu/Desktop/Interns/Round 2/John's Data/comp/Global/gsecd/gsecd_20190611.csv";
opts = detectImportOptions(gsecdFile);
opts.DataLines = [2,100001];
opts.SelectedVariableNames = {'gvkey','iid','datadate','loc','cshoc','curcdd','prccd'};
gsecdData = readtable(gsecdFile,opts);
tokeep = strcmp(gsecdData.iid,'01W');
gsecdData = gsecdData(tokeep, :);

gsecdData.gvkey(isnan(gsecdData.gvkey)) = -9999;
gsecdData.datadate(isnan(gsecdData.datadate)) = -9999;
gsecdData.seq(isnan(gsecdData.cshoc)) = -9999;
gsecdData.prccd(isnan(gsecdData.prccd)) = -9999;
gsecdData = unique(gsecdData, 'rows');
gsecdData.gvkey(gsecdData.gvkey == -9999) = nan;
gsecdData.datadate(gsecdData.datadate == -9999) = nan;
gsecdData.cshoc(gsecdData.cshoc == -9999) = nan;
gsecdData.prccd(gsecdData.prccd  == -9999) = nan;

%% Divide gsecd data into regions
gfundaData.region = zeros(height(gfundaData),1);
tochange = strcmp(gfundaData.loc,'USA') | strcmp(gfundaData.loc,'CAN');
gfundaData.region(tochange) = 1;
tochange = strcmp(gfundaData.loc,'JAP') | strcmp(gfundaData.fic,'JAP');
gfundaData.region(tochange) = 2;
tochange = strcmp(gfundaData.loc, 'AUS') | strcmp(gfundaData.loc, 'NZL') | strcmp(gfundaData.loc, 'HKG')| strcmp(gfundaData.loc,'SGP');
gfundaData.region(tochange) = 3;
tochange =  strcmp(gfundaData.loc, 'AUT') | strcmp(gfundaData.loc, 'BEL') | strcmp(gfundaData.loc, 'DNK') | strcmp(gfundaData.loc, 'FIN') | strcmp(gfundaData.loc, 'FRA') |...
    strcmp(gfundaData.loc, 'DEU') | strcmp(gfundaData.loc, 'GRC') | strcmp(gfundaData.loc, 'IRL') | strcmp(gfundaData.loc, 'ITA') | strcmp(gfundaData.loc, 'NLD') |...
    strcmp(gfundaData.loc, 'NOR') | strcmp(gfundaData.loc, 'PRT') | strcmp(gfundaData.loc, 'ESP') | strcmp(gfundaData.loc, 'SWE') | strcmp(gfundaData.loc, 'CHE') | strcmp(gfundaData.loc, 'GBR')
gfundaData.region(tochange) = 4;
tokeep = gfundaData.region ~= 0;
gfundaData = gfundaData(tokeep, :);

% Do I need to recheck for duplicates?

%% Minor adjustments to nan entries and negative prices

% txditc
tochange = isnan(gfundaData.txditc);
gfundaData.txditc(tochange) = 0;

% pstk
tochange = isnan(gfundaData.pstk);
gfundaData.pstk(tochange) = 0;

% prccd
tokeep = 0 < gsecdData.prccd;
gsecdData = gsecdData(tokeep, :);


%% Load in exchange rates & normalize things in terms of the dollar
data_curr1 = readtable("gexrtdly_20190611.txt");
data_curr1.(1)=[]; % if col 1 is python?s row number index
data_curr1.datadate=yyyymmdd(data_curr1.datadate);
tokeep=strcmp(data_curr1.tocurd,'USD');
data_curr2=data_curr1(tokeep,:);

% join curr and compute fxrates
data_curr1.Properties.VariableNames{'tocurd'}='curcdd';
gsecdData=outerjoin(gsecdData,data_curr1,'keys',{'datadate' 'curcdd'},'type','left','mergekeys',true);
gsecdData=outerjoin(gsecdData,data_curr2,'keys',{'datadate' 'fromcurd'},'type','left','mergekeys',true);
gsecdData.fxrates=gsecdData.(9)./gsecdData.(12);
gsecdData(:,9:12)=[];
gsecdData = sortrows(gsecdData,{'gvkey','datadate'});
gsecdData.usprice = gsecdData.(7)./gsecdData.(9);


%% Book value
gfundaData.BE = (gfundaData.seq + gfundaData.txditc - gfundaData.pstk).*1000000;
tokeep = gfundaData.BE > 0;
gfundaData = gfundaData(tokeep, :);

%% Market Value
gsecdData.ME = gsecdData.cshoc.*gsecdData.usprice.*1000;
tokeep = gsecdData.ME > 0;
gsecdData = gsecdData(tokeep, :);

%% Calculating the daily returns from price
returnsArray = zeros([height(gsecdData),1]);
tempMatrix = [gsecdData.gvkey, gsecdData.usprice];
for i = 2:length(tempMatrix)
    if tempMatrix(i,1) == tempMatrix(i-1,1)
        returnsArray(i) = tempMatrix(i,2)/tempMatrix(i-1,2)-1;
    end
end
gsecdData.returns = returnsArray;

%% Obtaining the monthly returns from the daily returns


% Get the # of monthly returns
numReturns = 0;
tempMatrix = [gsecdData.gvkey, gsecdData.datadate, gsecdData.usprice];
for i = 2:length(tempMatrix)
    if tempMatrix(i, 1) ~= tempMatrix(i-1,1) || floor(tempMatrix(i,2)/100) ~= floor(tempMatrix(i-1,2)/100)
        numReturns = numReturns + 1;
    end
end

% Generate a list of last entry for date/month
% [return index]
% Ask what to do if month skipped/years skipped
returnPoints = zeros([numReturns+1,1]);
returnPoints(length(returnPoints)) = length(returnsArray);
rpCounter = 1;
for i = 2:length(tempMatrix)
    % The starting indices should be different depending on if its new
    % gvkey or new month
    if tempMatrix(i,1) == tempMatrix(i-1,1) && floor(tempMatrix(i,2)/100) ~= floor(tempMatrix(i-1,2)/100)
        returnPoints(rpCounter) = i-1;
        rpCounter = rpCounter + 1;
    elseif tempMatrix(i,1) ~= tempMatrix(i-1,1)
        returnPoints(rpCounter,1) = i;
        rpCounter = rpCounter + 1;
    end
end

% Calculate the monthly returns
monthlyReturns = nan([numReturns,1]);
mrCounter = 1;
for i = 2:length(returnPoints)
    monthlyReturns(mrCounter) = tempMatrix(returnPoints(i),3)/tempMatrix(returnPoints(i-1),3)-1;
    mrCounter = mrCounter + 1;
end

returnPoints = returnPoints(2:length(returnPoints),:);

% Put monthly returns + date points into an array
monthlyReturnsFinal = nan([length(returnsArray),1]);
tempMatrix = [returnPoints, monthlyReturns];
for i = 1:length(returnPoints)
    monthlyReturnsFinal(tempMatrix(i,1)) = tempMatrix(i,2);
end

gsecdData.mreturns = monthlyReturnsFinal;

%% Annual Returns

% Calculate the number of annual returns
numReturns = 0;
tempMatrix = [gsecdData.gvkey, gsecdData.datadate, gsecdData.usprice];
for i = 2:length(tempMatrix)
    if tempMatrix(i, 1) ~= tempMatrix(i-1,1) || floor(tempMatrix(i,2)/10000) ~= floor(tempMatrix(i-1,2)/10000)
        numReturns = numReturns + 1;
    end
end

% Generate a list of last entry for date/year
% [return index]
% Ask what to do if month skipped/years skipped
returnPoints = zeros([numReturns+2,1]);
returnPoints(1) = 1;
returnPoints(length(returnPoints)) = length(returnsArray);
rpCounter = 2;
for i = 2:length(tempMatrix)
    % The starting indices should be different depending on if its new
    % gvkey or new year
    if tempMatrix(i,1) == tempMatrix(i-1,1) && floor(tempMatrix(i,2)/10000) ~= floor(tempMatrix(i-1,2)/10000)
        returnPoints(rpCounter) = i-1;
        rpCounter = rpCounter + 1;
    elseif tempMatrix(i,1) ~= tempMatrix(i-1,1)
        returnPoints(rpCounter,1) = i;
        rpCounter = rpCounter + 1;
    end
end

% Calculate the annual returns
annualReturns = nan([numReturns+1,1]);
arCounter = 1;
for i = 2:length(returnPoints)
    annualReturns(arCounter) = tempMatrix(returnPoints(i),3)/tempMatrix(returnPoints(i-1),3)-1;
    arCounter = arCounter + 1;
end

returnPoints = returnPoints(2:length(returnPoints),:);

% Put monthly returns + date points into an array
annualReturnsFinal = nan([length(returnsArray),1]);
tempMatrix = [returnPoints, annualReturns];
for i = 1:length(returnPoints)
    annualReturnsFinal(tempMatrix(i,1)) = tempMatrix(i,2);
end

gsecdData.areturns = annualReturnsFinal;




% This is more or less Robert's method

%{
numReturns = 0;
tempMatrix = [gsecdData.gvkey, gsecdData.datadate, gsecdData.usprice];
for i = 2:length(tempMatrix)
    if tempMatrix(i, 1) ~= tempMatrix(i-1,1) || floor(tempMatrix(i,2)/100) ~= floor(tempMatrix(i-1,2)/100)
        numReturns = numReturns + 1;
    end
end

% Generate a list of last entry for date/month
% [return index, isgvkeychange]
% Ask what to do if month skipped/years skipped
returnPoints = zeros([numReturns+1,1]);
returnPoints(length(returnPoints)) = length(returnsArray);
rpCounter = 1;
for i = 2:length(tempMatrix)
    % The starting indices should be different depending on if its new
    % gvkey or new month
    if tempMatrix(i,1) == tempMatrix(i-1,1) && floor(tempMatrix(i,2)/100) ~= floor(tempMatrix(i-1,2)/100)
        returnPoints(rpCounter) = i-1;
        rpCounter = rpCounter + 1;
    elseif tempMatrix(i,1) ~= tempMatrix(i-1,1)
        returnPoints(rpCounter,1) = i;
        rpCounter = rpCounter + 1;
    end
end

monthlyReturns = ones([numReturns+1,1]);
mrCounter = 1;
for i = 1:length(returnPoints)-1
    tempReturn = 1;
    for j = returnPoints(i,1):returnPoints(i+1,1)-1
        tempReturn = tempReturn * (1+returnsArray(i));
    end
    monthlyReturns(mrCounter) = tempReturn-1;
    mrCounter = mrCounter + 1;
end

monthlyReturnsFinal = nan([length(returnsArray),1]);
tempMatrix = [returnPoints, monthlyReturns];
for i = 1:length(returnPoints)
    monthlyReturnsFinal(tempMatrix(i,1)) = tempMatrix(i,2);
end

gsecdData.mreturns = monthlyReturnsFinal;

%}

%% Annual Returns
    
numReturns = 0;
tempMatrix = [gsecdData.gvkey, gsecdData.datadate, gsecdData.usprice];
for i = 2:length(tempMatrix)
    if tempMatrix(i, 1) ~= tempMatrix(i-1,1) || floor(tempMatrix(i,2)/10000) ~= floor(tempMatrix(i-1,2)/10000)
        numReturns = numReturns + 1;
    end
end

% Generate a list of last entry for date/year
% [return index, isgvkeychange]
% Ask what to do if month skipped/years skipped
returnPoints = zeros([numReturns+1,1]);
returnPoints(length(returnPoints)) = length(returnsArray);
rpCounter = 1;
for i = 2:length(tempMatrix)
    % The starting indices should be different depending on if its new
    % gvkey or new month
    if tempMatrix(i,1) == tempMatrix(i-1,1) && floor(tempMatrix(i,2)/10000) ~= floor(tempMatrix(i-1,2)/10000)
        returnPoints(rpCounter) = i-1;
        rpCounter = rpCounter + 1;
    elseif tempMatrix(i,1) ~= tempMatrix(i-1,1)
        returnPoints(rpCounter,1) = i;
        rpCounter = rpCounter + 1;
    end
end

annualReturns = ones([numReturns+1,1]);
arCounter = 1;
for i = 1:length(returnPoints)-1
    tempReturn = 1;
    for j = returnPoints(i,1):returnPoints(i+1,1)-1
        tempReturn = tempReturn * (1+returnsArray(i));
    end
    annualReturns(arCounter) = tempReturn-1;
    arCounter = arCounter + 1;
end

annualReturnsFinal = nan([length(returnsArray),1]);
tempMatrix = [returnPoints, annualReturns];
for i = 1:length(returnPoints)
    annualReturnsFinal(tempMatrix(i,1)) = tempMatrix(i,2);
end

gsecdData.areturns = annualReturnsFinal;
        
        


x = [0,0.0134890222888007,-0.00581950421316124,-0.142308732008472,-0.138504345874530,0.0904029974994609,0.0444914222493644,0.0160957530993862,-0.00212976824400979,-0.0257532210406607,0.0214611165061034,-0.0238732880430901,0.0382331165363263,0.0196063390686252,0.0238182115650614,0.0455823091954308,0.00459641069353389,0.0103618935740224,0.00421692030908361];
prod = 1;
for i = length(x)
    prod = prod * (1+x(i));
end
disp(prod);


