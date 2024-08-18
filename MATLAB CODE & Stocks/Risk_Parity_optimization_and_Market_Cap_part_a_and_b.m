%READ ME
%BEFORE YOU RUN THIS MATLAB CODE, MAKE SURE YOU SELECT THE SPAN OF DATES
%AND CORRESPONDING MARKET CAP WEIGHTS YOU ARE AFTER. SEE IN LINES 11-16 AND
%159-167

% Define the list of CSV file names
fileNames = {'F.csv', 'CAT.csv', 'DIS.csv','MCD.csv', 'KO.csv', 'PEP.csv','WMT.csv', 'C.csv', 'WFC.csv','JPM.csv', 'AAPL.csv', 'IBM.csv','PFE.csv', 'JNJ.csv', 'XOM.csv','MRO.csv', 'ED.csv', 'T.csv','VZ.csv', 'NEM.csv'}; % Add more filenames as needed

% Define the rolling periods
periods = {
    {'2007-12-31', '2011-06-30', '2011-07-29'};
    %{'2011-01-31', '2011-06-30', '2011-07-29'};
    %{'2011-01-31', '2011-07-29', '2011-08-31'};
    %{'2011-01-31', '2011-08-31', '2011-09-29'};
    %{'2011-01-31', '2011-09-29', '2011-10-31'};
    %{'2011-01-31', '2011-10-31', '2011-11-29'};
};

% Loop through each rolling period
for p = 1:length(periods)
    trainStart = periods{p}{1};
    trainEnd = periods{p}{2};
    testEnd = periods{p}{3};

    % Initialize arrays to store training and testing prices
    prices_train = [];
    prices_test = zeros(length(fileNames), 1);
    monthlyReturnsAllAssets = cell(1, length(fileNames));

    % Loop through each file and extract the adjusted closing prices
    for i = 1:length(fileNames)
        % Read the table from the CSV file
        T = readtable(fileNames{i});

        % Convert 'Date' to datetime
        T.Date = datetime(T.Date, 'InputFormat', 'yyyy-MM-dd');

        % Extract training period prices
        trainPrices = T.AdjClose(T.Date >= datetime(trainStart) & T.Date <= datetime(trainEnd));
        prices_train = [prices_train, trainPrices];

        % Extract testing period price
        price_test = T.AdjClose(T.Date == datetime(testEnd));
        
        % Store the testing price in the array
        prices_test(i) = price_test;
    end

    % Initialize a cell array to store monthly returns for each asset
    for fileIndex = 1:length(fileNames)
        % Load data from the current CSV file
        opts = detectImportOptions(fileNames{fileIndex});
        data = readtable(fileNames{fileIndex}, opts);

        % Convert date strings to datetime format
        data.Date = datetime(data.Date);

        % Filter data for the training period
        trainData = data(data.Date >= datetime(trainStart) & data.Date <= datetime(trainEnd), :);

        % Find the unique years and months in the training data
        [years, months, ~] = ymd(trainData.Date);
        uniqueYearMonths = unique(years*100 + months);

        % Preallocate the output variables
        endOfMonthDates = NaT(size(uniqueYearMonths));
        adjCloseLastDayOfMonth = zeros(size(uniqueYearMonths));

        % Loop through each unique year-month combination
        for i = 1:length(uniqueYearMonths)
            year = floor(uniqueYearMonths(i) / 100);
            month = uniqueYearMonths(i) - year * 100;

            % Find indices of the current month and year
            ind = find(years == year & months == month);

            % Find the last day of the current month in the data
            [~, lastDayIndex] = max(trainData.Date(ind));
            actualIndex = ind(lastDayIndex);

            % Store the results
            endOfMonthDates(i) = trainData.Date(actualIndex);
            adjCloseLastDayOfMonth(i) = trainData.AdjClose(actualIndex);
        end

        % Calculate monthly returns for the current asset
        numMonths = length(adjCloseLastDayOfMonth);
        monthlyReturns = zeros(numMonths-1, 1);

        for i = 2:numMonths
            monthlyReturns(i-1) = ((adjCloseLastDayOfMonth(i) / adjCloseLastDayOfMonth(i-1)) - 1);
        end

        % Handle NaN values in monthly returns
        monthlyReturns = monthlyReturns(~isnan(monthlyReturns));

        % Store monthly returns of the current asset for later analysis
        monthlyReturnsAllAssets{fileIndex} = monthlyReturns;
    end

    % Analysis of all assets
    num_assets = length(monthlyReturnsAllAssets);
    averages = zeros(1, num_assets);
    expected_returns = zeros(1, num_assets);
    var = zeros(1, num_assets);

    % Concatenate all monthly returns for covariance calculation
    allReturns = horzcat(monthlyReturnsAllAssets{:});

    % Calculate the covariance matrix for all assets
    cov_matrix = cov(allReturns);

    % Calculate realized returns for the training period
    returns_train = diff(prices_train) ./ prices_train(1:end-1, :);

    % Risk Parity Optimization
    % Number of assets (companies)
    n = length(fileNames);

    % Combine all monthly returns into a single matrix
    allMonthlyReturns = cell2mat(monthlyReturnsAllAssets');

    % Compute the covariance matrix of the monthly returns
    Q = cov(allMonthlyReturns);

    % Initial guess for portfolio weights and theta
    x0 = ones(n, 1) / n;
    theta0 = 0.01;
    x_theta0 = [x0; theta0];

    % Define the objective function for fmincon
    objective = @(x_theta) sum((x_theta(1:n) .* (Q * x_theta(1:n)) - x_theta(n+1)).^2);

    % Define the equality constraint
    Aeq = [ones(1, n), 0];
    beq = 1;

    % Define bounds for weights and theta
    lb = [zeros(n, 1); -Inf];  % Restrict weights to be non-negative
    ub = [ones(n, 1); Inf];

    % Set options for fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

    % Solve the optimization problem using fmincon
    [x_theta_opt, fval, exitflag, output] = fmincon(objective, x_theta0, [], [], Aeq, beq, lb, ub, [], options);

    % Extract optimal weights and theta
    x_opt = x_theta_opt(1:n);
    theta_opt = x_theta_opt(n+1);
    
    % Display risk parity results
    disp('Optimal weights for risk parity:');
    disp(x_opt);
    disp('Optimal theta:');
    disp(theta_opt);

    %%Market Cap Weights for June 2011
    weights_market_cap = [0.0187; 0.0245; 0.0264; 0.0312; 0.0550; 0.0398; 0.0659; 0.0434; 0.0530; 0.0581; 0.1109; 0.0742; 0.0581; 0.0651; 0.1432; 0.0134; 0.0056; 0.0664; 0.0376; 0.0094];
    %%Market Cap Weights for July 2011
    %weights_market_cap = [0.0167; 0.0230; 0.0263; 0.0324; 0.0561; 0.0366; 0.0660; 0.0404; 0.0533; 0.0580; 0.1306; 0.0783; 0.0548; 0.0640; 0.1417; 0.0080; 0.0055; 0.0625; 0.0360; 0.0098];
    %%Market Cap Weights for August 2011
    %weights_market_cap = [0.0158; 0.0220; 0.0237; 0.0350; 0.0606; 0.0382; 0.0692; 0.0340; 0.0516; 0.0549; 0.1337; 0.0769; 0.0555; 0.0676; 0.1349; 0.0072; 0.0062; 0.0632; 0.0384; 0.0114];
    %%Market Cap Weights for September 2011
    %weights_market_cap = [0.0144; 0.0187; 0.0220; 0.0356; 0.0609; 0.0385; 0.0702; 0.0294; 0.0500; 0.0461; 0.1388; 0.0820; 0.0542; 0.0685; 0.1387; 0.0060; 0.0066; 0.0664; 0.0409; 0.0121];
    %%Market Cap Weights for October 2011
    %weights_market_cap = [0.0163; 0.0224; 0.0237; 0.0351; 0.0569; 0.0361; 0.0717; 0.0338; 0.0502; 0.0497; 0.1380; 0.0798; 0.0551; 0.0647; 0.1393; 0.0068; 0.0062; 0.0637; 0.0384; 0.0120];


    % Calculate realized returns for the testing month
    realizedReturns = (prices_test - prices_train(end, :)') ./ prices_train(end, :)';

    cov_matrix = cov(returns_train);

    % Calculate portfolio return for the testing period
    R_p_risk_parity = sum(x_opt .* realizedReturns);
    R_p_market_cap = sum(weights_market_cap .* realizedReturns);

    % Calculate portfolio variance using the covariance matrix
    sigma_p2_risk_parity = x_opt' * cov_matrix * x_opt;
    sigma_p2_market_cap = weights_market_cap' * cov_matrix * weights_market_cap;

    % Calculate portfolio standard deviation
    sigma_p_risk_parity = sqrt(sigma_p2_risk_parity);
    sigma_p_market_cap = sqrt(sigma_p2_market_cap);

    monthlyRiskFreeRate = 0;

    % Calculate Sharpe Ratio
    sharpeRatio_risk_parity = (R_p_risk_parity - monthlyRiskFreeRate) / sigma_p_risk_parity;
    sharpeRatio_market_cap = (R_p_market_cap - monthlyRiskFreeRate) / sigma_p_market_cap;


    % Display the results for the current period
    fprintf('Results for testing month ending %s:\n', testEnd);
    fprintf('Risk Parity: Return = %.4f, Variance = %.4f, StdDev = %.4f, Sharpe Ratio = %.4f\n', R_p_risk_parity, sigma_p2_risk_parity, sigma_p_risk_parity, sharpeRatio_risk_parity);
    fprintf('Market Cap: Return = %.4f, Variance = %.4f, StdDev = %.4f, Sharpe Ratio = %.4f\n\n', R_p_market_cap, sigma_p2_market_cap, sigma_p_market_cap, sharpeRatio_market_cap);

end
