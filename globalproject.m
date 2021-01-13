%% Clear
clearvars;
clc;

%% Load tables
file_in = "/Users/zacharyzhu/Desktop/Interns2/Round 2/John's Data/comp/Global/gfunda/gfunda_20190628.csv";
opts = detectImportOptions(file_in);
opts.SelectedVariableNames = {'gvkey' 'datadate' 'loc' 'fic' 'curcd' 'indfmt' 'seq' 'txditc' 'pstk'}; %'ceq' 'at' 'lt' 'txdb'
gfunda_table = readtable(file_in, opts);
tokeep=strcmp(gfunda_table.indfmt,'INDL');
gfunda_table=gfunda_table(tokeep,:);
gfunda_table.indfmt=[];
gvkey_list = unique(gfunda_table.gvkey);

dates = str2num(datestr(datetime(1990, 06, 01):calmonths:datetime(2019, 07, 01), 'yyyymm'));
filein = "/Users/zacharyzhu/Desktop/Interns2/Round 2/John's Data/comp/Global/gsecd_ver2_pt2_bydate/gsecd_"+dates+".csv" ;
gsecm_cell = cell(350, 1);
for ii = 1:length(dates)
    gsecm_cell{ii, 1} = readtable(filein(ii));
end
gsecm_table = gsecm_cell{1, 1};
for ii = 2:350
     temp_table = gsecm_cell{ii, 1};
     gsecm_table = [gsecm_table; temp_table];
end
gsecm_table = sortrows(gsecm_table, {'gvkey', 'datadate'});
gsecm_198912 = readtable("/Users/zacharyzhu/Desktop/Interns2/Round 2/John's Data/comp/Global/gsecd_ver2_pt2_bydate/gsecd_198912.csv");

gexrtdly_table = readtable("/Users/zacharyzhu/Desktop/Interns2/Round 2/John's Data/comp/Global/gexrtdly/gexrtdly_20190611.txt");

%% Group gfunda by region
gfunda_table.region_id = nan(height(gfunda_table), 1);
tochange = strcmp(gfunda_table.loc, 'USA') | strcmp(gfunda_table.loc, 'CAN') | strcmp(gfunda_table.fic, 'USA') | strcmp(gfunda_table.fic, 'CAN');
gfunda_table.region_id(tochange) = 1;
tochange = strcmp(gfunda_table.loc, 'JPN') | strcmp(gfunda_table.fic, 'JPN');
gfunda_table.region_id(tochange) = 2;
tochange = strcmp(gfunda_table.loc, 'AUS') | strcmp(gfunda_table.loc, 'NZL') | strcmp(gfunda_table.loc, 'HKG') | strcmp(gfunda_table.loc, 'SGP') | strcmp(gfunda_table.fic, 'AUS') | strcmp(gfunda_table.fic, 'NZL') | strcmp(gfunda_table.fic, 'HKG') | strcmp(gfunda_table.fic, 'SGP');
gfunda_table.region_id(tochange) = 3;
tochange = strcmp(gfunda_table.loc, 'AUT') | strcmp(gfunda_table.loc, 'BEL') | strcmp(gfunda_table.loc, 'DNK') | strcmp(gfunda_table.loc, 'FIN') | strcmp(gfunda_table.loc, 'FRA') | strcmp(gfunda_table.loc, 'DEU') | strcmp(gfunda_table.loc, 'GRC') | strcmp(gfunda_table.loc, 'IRL') | strcmp(gfunda_table.loc, 'ITA') | strcmp(gfunda_table.loc, 'NLD') | strcmp(gfunda_table.loc, 'NOR') | strcmp(gfunda_table.loc, 'PRT') | strcmp(gfunda_table.loc, 'ESP') | strcmp(gfunda_table.loc, 'SWE') | strcmp(gfunda_table.loc, 'CHE') | strcmp(gfunda_table.loc, 'GBR') | strcmp(gfunda_table.fic, 'AUT') | strcmp(gfunda_table.fic, 'BEL') | strcmp(gfunda_table.fic, 'DNK') | strcmp(gfunda_table.fic, 'FIN') | strcmp(gfunda_table.fic, 'FRA') | strcmp(gfunda_table.fic, 'DEU') | strcmp(gfunda_table.fic, 'GRC') | strcmp(gfunda_table.fic, 'IRL') | strcmp(gfunda_table.fic, 'ITA') | strcmp(gfunda_table.fic, 'NLD') | strcmp(gfunda_table.fic, 'NOR') | strcmp(gfunda_table.fic, 'PRT') | strcmp(gfunda_table.fic, 'ESP') | strcmp(gfunda_table.fic, 'SWE') | strcmp(gfunda_table.fic, 'CHE') | strcmp(gfunda_table.fic, 'GBR');
gfunda_table.region_id(tochange) = 4;
tokeep = ~isnan(gfunda_table.region_id);
gfunda_table = gfunda_table(tokeep, :);

%% Screen data
% remove duplicates from gfunda
gfunda_table.gvkey(isnan(gfunda_table.gvkey)) = -9999;
gfunda_table.datadate(isnan(gfunda_table.datadate)) = -9999;
gfunda_table.seq(isnan(gfunda_table.seq)) = -9999;
gfunda_table.txditc(isnan(gfunda_table.txditc)) = -9999;
gfunda_table.pstk(isnan(gfunda_table.pstk)) = -9999;
gfunda_table = unique(gfunda_table, 'rows');
gfunda_table.gvkey(gfunda_table.gvkey == -9999) = nan;
gfunda_table.datadate(gfunda_table.datadate == -9999) = nan;
gfunda_table.seq(gfunda_table.seq == -9999) = nan;
gfunda_table.txditc(gfunda_table.txditc == -9999) = nan;
gfunda_table.pstk(gfunda_table.pstk == -9999) = nan;

% remove duplicates from gsecm
gsecm_table.gvkey(isnan(gsecm_table.gvkey)) = -9999;
gsecm_table.datadate(isnan(gsecm_table.datadate)) = -9999;
gsecm_table.date_m(isnan(gsecm_table.date_m)) = -9999;
gsecm_table.exchg(isnan(gsecm_table.exchg)) = -9999;
gsecm_table.prccd_j(isnan(gsecm_table.prccd_j)) = -9999;
gsecm_table.cshoc_j(isnan(gsecm_table.cshoc_j)) = -9999;
gsecm_table.vol_m(isnan(gsecm_table.vol_m)) = -9999;
gsecm_table.ret_m(isnan(gsecm_table.ret_m)) = -9999;
gsecm_table = unique(gsecm_table, 'rows');
gsecm_table.gvkey(gsecm_table.gvkey == -9999) = nan;
gsecm_table.datadate(gsecm_table.datadate == -9999) = nan;
gsecm_table.date_m(gsecm_table.date_m == -9999) = nan;
gsecm_table.exchg(gsecm_table.exchg == -9999) = nan;
gsecm_table.prccd_j(gsecm_table.prccd_j == -9999) = nan;
gsecm_table.cshoc_j(gsecm_table.cshoc_j == -9999) = nan;
gsecm_table.vol_m(gsecm_table.vol_m == -9999) = nan;
gsecm_table.ret_m(gsecm_table.ret_m == -9999) = nan;

% add adjdate
gfunda_table.adjdate = dateshift(datetime(gfunda_table.datadate, 'convertfrom', 'yyyyMMdd'), 'end', 'year');
gsecm_table.adjdate = dateshift(datetime(gsecm_table.datadate, 'convertfrom', 'yyyyMMdd'), 'end', 'month');
gsecm_198912.adjdate = dateshift(datetime(gsecm_198912.datadate, 'convertfrom', 'yyyyMMdd'), 'end', 'month');

% txditc
tochange = isnan(gfunda_table.txditc);
gfunda_table.txditc(tochange) = 0;

% pstk
tochange = isnan(gfunda_table.pstk);
gfunda_table.pstk(tochange) = 0;

% prccd_j
tokeep = 0 < gsecm_table.prccd_j;
gsecm_table = gsecm_table(tokeep, :);

%% Currency conversion
tokeep = strcmp(gexrtdly_table.tocurd, 'USD');
gexrtdly_table.(1) = [];
gexrtdly_usdexrt = gexrtdly_table(tokeep, [4 2]);
gexrtdly_table = outerjoin(gexrtdly_table, gexrtdly_usdexrt, 'keys', {'datadate'}, 'type', 'left', 'mergekeys', true);
gexrtdly_table.Properties.VariableNames = {'curcd', 'exratd', 'fromcurd', 'datadate', 'exratd_usd'};

gfunda_table.datadate = datetime(gfunda_table.datadate, 'convertfrom', 'yyyyMMdd');
gfunda_table = outerjoin(gfunda_table, gexrtdly_table, 'keys', {'datadate', 'curcd'}, 'type', 'left', 'mergekeys', true);
gfunda_table.datadate = yyyymmdd(gfunda_table.datadate);
gfunda_table.fxrates = gfunda_table.exratd./gfunda_table.exratd_usd;

%% Calculate book value
gfunda_table.BE = gfunda_table.seq + gfunda_table.txditc - gfunda_table.pstk;
gfunda_table.BE = gfunda_table.BE./gfunda_table.fxrates.*1000000; % make sure this is still supposed to be 1000000
tokeep = gfunda_table.BE > 0;
gfunda_table = gfunda_table(tokeep, :);

tokeep = year(gfunda_table.adjdate) == 1989;
gfunda_1989 = gfunda_table(tokeep, :);
gfunda_1989 = sortrows(gfunda_1989, {'gvkey', 'adjdate'});

%% Calculate market value
gsecm_table.ME = gsecm_table.cshoc_j.*gsecm_table.prccd_j;
tokeep = gsecm_table.ME > 0;
gsecm_table = gsecm_table(tokeep, :);

%% Merge the tables
BEME_table = gsecm_table(:, {'gvkey', 'adjdate', 'ret_m', 'ME'});
BEME_table.Properties.VariableNames = {'gvkey', 'adjdate', 'ret', 'ME'};

BE_table = gfunda_table(:, {'gvkey', 'adjdate', 'BE', 'region_id'});
BEME_table = outerjoin(BEME_table, BE_table, 'keys', {'gvkey', 'adjdate'}, 'type', 'left', 'mergekeys', true);
BEME_table.BEME = BEME_table.BE./BEME_table.ME;

gsecm_198912 = gsecm_198912(:, {'gvkey', 'adjdate', 'ME'});
gfunda_1989 = gfunda_1989(:, {'gvkey', 'adjdate', 'region_id', 'BE'});
BEME_table_1989 = outerjoin(gsecm_198912, gfunda_1989, 'keys', {'gvkey', 'adjdate'}, 'type', 'left', 'mergekeys', true);
tokeep = ones(size(BEME_table_1989.gvkey));
for ii = 1:height(BEME_table_1989) - 1
    if BEME_table_1989.gvkey(ii) == BEME_table_1989.gvkey(ii+1)
        tokeep(ii+1) = 0;
    end
end
tokeep = logical(tokeep);
BEME_table_1989 = BEME_table_1989(tokeep, :);
BEME_table_1989.BEME = BEME_table_1989.BE./BEME_table_1989.ME;

%% Quantiling
month_list = unique(dateshift(BEME_table.adjdate, 'end', 'month'));
month_len = length(month_list);
bm_vars = BEME_table_1989(:, [1 6 4]);

BEME_cell = cell(month_len,1);
for ii = 1:month_len
    tokeep = BEME_table.adjdate == month_list(ii);
    temp = BEME_table(tokeep, :);
    BEME_cell{ii, 2} = temp;
    BEME_cell{ii, 1} = yyyymmdd(month_list(ii));
end

SG_array_NA = nan(29, 1); SN_array_NA = nan(29, 1); SV_array_NA = nan(29, 1);
BG_array_NA = nan(29, 1); BN_array_NA = nan(29, 1); BV_array_NA = nan(29, 1);
SG_array_JP = nan(29, 1); SN_array_JP = nan(29, 1); SV_array_JP = nan(29, 1); 
BG_array_JP = nan(29, 1); BN_array_JP = nan(29, 1); BV_array_JP = nan(29, 1);
SG_array_AP = nan(29, 1); SN_array_AP = nan(29, 1); SV_array_AP = nan(29, 1);
BG_array_AP = nan(29, 1); BN_array_AP = nan(29, 1); BV_array_AP = nan(29, 1);
SG_array_EU = nan(29, 1); SN_array_EU = nan(29, 1); SV_array_EU = nan(29, 1);
BG_array_EU = nan(29, 1); BN_array_EU = nan(29, 1); BV_array_EU = nan(29, 1);
for ii = 1:25
    ret_array = BEME_cell{12*(ii-1)+1, 2}(:, [1 4 7 6]);
    
    if ii ~= 1
        bm_vars = BEME_cell{12*(ii-1)+1-6, 2}(:, [1 7 6]);
    end
    % get tp1:tp12 ret
    for jj = 1:12
        temp2 = BEME_cell{12*(ii-1)+1+jj, 2};
        temp2 = temp2(:, {'gvkey' 'ret'});
        temp2.Properties.VariableNames = ["gvkey" "ret"+"_tp"+jj];
        ret_array = outerjoin(ret_array, temp2, 'keys', 'gvkey', 'type', 'left', 'mergekeys', true);
    end
    
    % get ret for 12 months
    ret_array.ret = prod(ret_array{:, 5:16}+1, 2, 'omitnan')-1;
    
    % join ret_array with bm_stuff
    tokeep = ~isnan(bm_vars.BEME) & ~isnan(bm_vars.region_id);
    bm_vars = bm_vars(tokeep, :);
    ret_array = outerjoin(ret_array, bm_vars, 'keys', {'gvkey'}, 'type', 'left', 'mergekeys', true);
    ret_array(:, 3:4) = ret_array(:, 18:19);
    ret_array(:, 18:19) = [];
    ret_array.Properties.VariableNames{3} = 'BEME';
    ret_array.Properties.VariableNames{4} = 'region_id';
    
    % separate by region
    tokeep = ret_array.region_id == 1;
    ret_array_NA = ret_array(tokeep, :);
    tokeep = ret_array.region_id == 2;
    ret_array_JP = ret_array(tokeep, :);
    tokeep = ret_array.region_id == 3;
    ret_array_AP = ret_array(tokeep, :);
    tokeep = ret_array.region_id == 4;
    ret_array_EU = ret_array(tokeep, :);
    
    % breakpoints
    % size
    bp_sz_NA = quantile(ret_array_NA.ME, [0.1 0.9]);
    bp_sz_JP = quantile(ret_array_JP.ME, [0.1 0.9]);
    bp_sz_AP = quantile(ret_array_AP.ME, [0.1 0.9]);
    bp_sz_EU = quantile(ret_array_EU.ME, [0.1 0.9]);
    % b/m
    tokeep = ret_array_NA.ME > bp_sz_NA(2);
    NA_val = ret_array_NA(tokeep, :);
    bp_bm_NA = quantile(NA_val.BEME, [0.3 0.7]);
    tokeep = ret_array_JP.ME > bp_sz_JP(2);
    JP_val = ret_array_JP(tokeep, :);
    bp_bm_JP = quantile(JP_val.BEME, [0.3 0.7]);
    tokeep = ret_array_AP.ME > bp_sz_AP(2);
    AP_val = ret_array_AP(tokeep, :);
    bp_bm_AP = quantile(AP_val.BEME, [0.3 0.7]);
    tokeep = ret_array_EU.ME > bp_sz_EU(2);
    EU_val = ret_array_EU(tokeep, :);
    bp_bm_EU = quantile(EU_val.BEME, [0.3 0.7]);
    
    % group by bp
    % NA
    % size
    ret_array_NA{:, 'sz_id'} = nan(size(ret_array_NA, 1), 1);
    ret_array_NA{ret_array_NA{:, 'ME'} < bp_sz_NA(1), 'sz_id'} = 1;
    ret_array_NA{ret_array_NA{:, 'ME'} > bp_sz_NA(2), 'sz_id'} = 2;
    % b/m
    ret_array_NA{:, 'bm_id'} = nan(size(ret_array_NA, 1), 1);
    ret_array_NA{ret_array_NA{:, 'BEME'} < bp_bm_NA(1), 'bm_id'} = 1;
    ret_array_NA{bp_bm_NA(1) < ret_array_NA{:, 'BEME'} & ret_array_NA{:, 'BEME'} < bp_bm_NA(2), 'bm_id'} = 2;
    ret_array_NA{bp_bm_NA(2) < ret_array_NA{:, 'BEME'}, 'bm_id'} = 3;
    
    % JP
    % size
    ret_array_JP{:, 'sz_id'} = nan(size(ret_array_JP, 1), 1);
    ret_array_JP{ret_array_JP{:, 'ME'} < bp_sz_JP(1), 'sz_id'} = 1;
    ret_array_JP{ret_array_JP{:, 'ME'} > bp_sz_JP(2), 'sz_id'} = 2;
    % b/m
    ret_array_JP{:, 'bm_id'} = nan(size(ret_array_JP, 1), 1);
    ret_array_JP{ret_array_JP{:, 'BEME'} < bp_bm_JP(1), 'bm_id'} = 1;
    ret_array_JP{bp_bm_JP(1) < ret_array_JP{:, 'BEME'} & ret_array_JP{:, 'BEME'} < bp_bm_JP(2), 'bm_id'} = 2;
    ret_array_JP{bp_bm_JP(2) < ret_array_JP{:, 'BEME'}, 'bm_id'} = 3;
    
    % AP
    % size
    ret_array_AP{:, 'sz_id'} = nan(size(ret_array_AP, 1), 1);
    ret_array_AP{ret_array_AP{:, 'ME'} < bp_sz_AP(1), 'sz_id'} = 1;
    ret_array_AP{ret_array_AP{:, 'ME'} > bp_sz_AP(2), 'sz_id'} = 2;
    % b/m
    ret_array_AP{:, 'bm_id'} = nan(size(ret_array_AP, 1), 1);
    ret_array_AP{ret_array_AP{:, 'BEME'} < bp_bm_AP(1), 'bm_id'} = 1;
    ret_array_AP{bp_bm_AP(1) < ret_array_AP{:, 'BEME'} & ret_array_AP{:, 'BEME'} < bp_bm_AP(2), 'bm_id'} = 2;
    ret_array_AP{bp_bm_AP(2) < ret_array_AP{:, 'BEME'}, 'bm_id'} = 3;
    
    % EU
    % size
    ret_array_EU{:, 'sz_id'} = nan(size(ret_array_EU, 1), 1);
    ret_array_EU{ret_array_EU{:, 'ME'} < bp_sz_EU(1), 'sz_id'} = 1;
    ret_array_EU{ret_array_EU{:, 'ME'} > bp_sz_EU(2), 'sz_id'} = 2;
    % b/m
    ret_array_EU{:, 'bm_id'} = nan(size(ret_array_EU, 1), 1);
    ret_array_EU{ret_array_EU{:, 'BEME'} < bp_bm_EU(1), 'bm_id'} = 1;
    ret_array_EU{bp_bm_EU(1) < ret_array_EU{:, 'BEME'} & ret_array_EU{:, 'BEME'} < bp_bm_EU(2), 'bm_id'} = 2;
    ret_array_EU{bp_bm_EU(2) < ret_array_EU{:, 'BEME'}, 'bm_id'} = 3;
    
    % returns
        % NA
    for jj = 1:2
        tokeep = ret_array_NA{:, 'sz_id'} == jj;
        temp = ret_array_NA(tokeep, :);
        for kk = 1:3
            tokeep = temp{:, 'bm_id'} == kk;
            temp2 = temp(tokeep, {'ME', 'ret'});
            temp2.weight = (temp2.ME./sum(temp2.ME))+1;
            temp2.weighted_ret = temp2.ret.*temp2.weight;
            tempr = temp2.weighted_ret;
            if jj == 1 && kk == 1
                SG_array_NA(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 2
                SN_array_NA(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 3
                SV_array_NA(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 1
                BG_array_NA(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 2
                BN_array_NA(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 3
                BV_array_NA(ii, 1) = nanmean(tempr);
            end
        end
    end
    
    % JP
    for jj = 1:2
        tokeep = ret_array_JP{:, 'sz_id'} == jj;
        temp = ret_array_JP(tokeep, :);
        for kk = 1:3
            tokeep = temp{:, 'bm_id'} == kk;
            temp2 = temp(tokeep, {'ME', 'ret'});
            temp2.weight = (temp2.ME./sum(temp2.ME))+1;
            temp2.weighted_ret = temp2.ret.*temp2.weight;
            tempr = temp2.weighted_ret;
            if jj == 1 && kk == 1
                SG_array_JP(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 2
                SN_array_JP(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 3
                SV_array_JP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 1
                BG_array_JP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 2
                BN_array_JP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 3
                BV_array_JP(ii, 1) = nanmean(tempr);
            end
        end
    end
    
        % AP
    for jj = 1:2
        tokeep = ret_array_AP{:, 'sz_id'} == jj;
        temp = ret_array_AP(tokeep, :);
        for kk = 1:3
            tokeep = temp{:, 'bm_id'} == kk;
            temp2 = temp(tokeep, {'ME', 'ret'});
            temp2.weight = (temp2.ME./sum(temp2.ME))+1;
            temp2.weighted_ret = temp2.ret.*temp2.weight;
            tempr = temp2.weighted_ret;
            if jj == 1 && kk == 1
                SG_array_AP(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 2
                SN_array_AP(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 3
                SV_array_AP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 1
                BG_array_AP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 2
                BN_array_AP(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 3
                BV_array_AP(ii, 1) = nanmean(tempr);
            end
        end
    end
    
        % EU
    for jj = 1:2
        tokeep = ret_array_EU{:, 'sz_id'} == jj;
        temp = ret_array_EU(tokeep, :);
        for kk = 1:3
            tokeep = temp{:, 'bm_id'} == kk;
            temp2 = temp(tokeep, {'ME', 'ret'});
            temp2.weight = (temp2.ME./sum(temp2.ME))+1;
            temp2.weighted_ret = temp2.ret.*temp2.weight;
            tempr = temp2.weighted_ret;
            if jj == 1 && kk == 1
                SG_array_EU(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 2
                SN_array_EU(ii, 1) = nanmean(tempr);
            elseif jj == 1 && kk == 3
                SV_array_EU(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 1
                BG_array_EU(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 2
                BN_array_EU(ii, 1) = nanmean(tempr);
            elseif jj == 2 && kk == 3
                BV_array_EU(ii, 1) = nanmean(tempr);
            end
        end
    end
    
end

%% Return series
% HMLs
HMLs_NA = SV_array_NA - SG_array_NA;
HMLs_JP = SV_array_JP - SG_array_JP;
HMLs_AP = SV_array_AP - SG_array_AP;
HMLs_EU = SV_array_EU - SG_array_EU;

% HMLb
HMLb_NA = BV_array_NA - BG_array_NA;
HMLb_JP = BV_array_JP - BG_array_JP;
HMLb_AP = BV_array_AP - BG_array_AP;
HMLb_EU = BV_array_EU - BG_array_EU;

% HML
HML_array_NA = mean([HMLs_NA, HMLb_NA], 2);
HML_array_JP = mean([HMLs_JP, HMLb_JP], 2);
HML_array_AP = mean([HMLs_AP, HMLb_AP], 2);
HML_array_EU = mean([HMLs_EU, HMLb_EU], 2);

HML_NA = nanmean(HML_array_NA);
HML_JP = nanmean(HML_array_JP);
HML_AP = nanmean(HML_array_AP);
HML_EU = nanmean(HML_array_EU);

% plot(cumprod(HML_array_EU+1))

% SMB
`

%% Notes
% making 2x3(5x5?) sort based on size and book/market
% size breakpoints - in paper
% b/m breakpoints - in paper

% june months become part of portfolio formed the following june (why??)
% 2018630 or 20180601--> 20190630
% calculate monthly returns last day of previous month to last day of month
% (next/previous)^

% BEME needs to be valid in Dec, prc needs to be valid in june
% rebalance BEME monthly or find it once in Dec for whole year


% 19870630-20190531, raw
% 19891231-20171231, BEME
% 19900630-20180630, portfolio formation
% 19910630-20190630, returns



% 3 part presentation: 
% 1 - paper and data
% 2 - code and fixes
% 3 - results, experiences, possible improvements
