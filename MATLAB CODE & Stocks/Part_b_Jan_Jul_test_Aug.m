
% List of CSV files for different assets
csvFiles = {'F.csv', 'CAT.csv', 'DIS.csv', 'MCD.csv', 'KO.csv', 'PEP.csv', 'WMT.csv', 'C.csv', 'WFC.csv', 'JPM.csv', 'AAPL.csv', 'IBM.csv', 'PFE.csv', 'JNJ.csv', 'XOM.csv', 'MRO.csv', 'ED.csv', 'T.csv', 'VZ.csv', 'NEM.csv'}; % Add more filenames as needed

% Initialize cell arrays to store adjusted closing prices for July and August
adjCloseJuly2011 = zeros(1, length(csvFiles));
adjCloseAugust2011 = zeros(1, length(csvFiles));

% Iterate over each asset
for fileIndex = 1:length(csvFiles)
    
    % Load Data_Stock_Features from the current CSV file
    Stock_Features = detectImportOptions(csvFiles{fileIndex});
    Data_Stock_Features = readtable(csvFiles{fileIndex}, Stock_Features);
    Data_Stock_Features.Date = datetime(Data_Stock_Features.Date);

    % Filter out dates for July 2011
    Data_Stock_Features_July2011 = Data_Stock_Features(month(Data_Stock_Features.Date) == 7 & year(Data_Stock_Features.Date) == 2011, :);
    
    % Extract adjusted closing price for the last trading day of July 2011
    adjCloseJuly2011(fileIndex) = Data_Stock_Features_July2011.AdjClose(end);

    % Filter out dates for August 2011
    Data_Stock_Features_August2011 = Data_Stock_Features(month(Data_Stock_Features.Date) == 8 & year(Data_Stock_Features.Date) == 2011, :);
    
    % Extract adjusted closing price for the last trading day of August 2011
    adjCloseAugust2011(fileIndex) = Data_Stock_Features_August2011.AdjClose(end);
end

% Calculate monthly returns for August 2011
monthlyReturnsAugust2011 = (adjCloseAugust2011 ./ adjCloseJuly2011) - 1;

% Display monthly returns for each asset in August 2011
for assetIndex = 1:length(csvFiles)
    disp(['Asset: ', csvFiles{assetIndex}]);
    disp(['Monthly Return (August 2011): ', num2str(monthlyReturnsAugust2011(assetIndex))]);
    fprintf('\n');
end


%Portfolio return for just Aug
%MVO - Aug

Portfolio_return_Aug_MVO = monthlyReturnsAugust2011 * x; 

sharpeRatio = Portfolio_return_Aug_MVO / sqrt(std_devi);

disp('Portfolio return for Aug MVO');
disp(Portfolio_return_Aug_MVO);
disp('Variance of the portfolio');
disp(std_devi);
disp('Standard deviation of the portfolio');
disp(sqrt(std_devi));
disp('Sharpe Ratio');
disp(sharpeRatio);

%Robust MVO - Aug

Portfolio_return_Aug_RMVO = monthlyReturnsAugust2011 * x_RMVO; 

sharpeRatio_1 = Portfolio_return_Aug_RMVO / sqrt(std_devi_RMVO);

disp('Portfolio return for Aug Robust MVO (CL=95%)');
disp(Portfolio_return_Aug_RMVO);
disp('Variance of the portfolio');
disp(std_devi_RMVO);
disp('Standard deviation of the portfolio');
disp(sqrt(std_devi_RMVO));
disp('Sharpe Ratio');
disp(sharpeRatio_1);

Portfolio_return_Aug_RMVO_1 = monthlyReturnsAugust2011 * x_RMVO_1; 

sharpeRatio_2 = Portfolio_return_Aug_RMVO_1 / sqrt(std_devi_RMVO_1);

disp('Portfolio return for Aug Robust MVO (CL=90%)');
disp(Portfolio_return_Aug_RMVO_1);
disp('Variance of the portfolio');
disp(std_devi_RMVO_1);
disp('Standard deviation of the portfolio');
disp(sqrt(std_devi_RMVO_1));
disp('Sharpe Ratio');
disp(sharpeRatio_2);
