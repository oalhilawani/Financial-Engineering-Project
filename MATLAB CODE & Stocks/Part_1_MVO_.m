% List of CSV files for different assets
csvFiles = {'F.csv', 'CAT.csv', 'DIS.csv', 'MCD.csv', 'KO.csv', 'PEP.csv', 'WMT.csv', 'C.csv', 'WFC.csv', 'JPM.csv', 'AAPL.csv', 'IBM.csv', 'PFE.csv', 'JNJ.csv', 'XOM.csv', 'MRO.csv', 'ED.csv', 'T.csv', 'VZ.csv', 'NEM.csv'}; % Add more filenames as needed

% Initialize a cell array to store monthly returns for each asset
monthlyReturnsAllAssets = cell(1, length(csvFiles));
sumMonthlyReturns = zeros(1, length(csvFiles));
sampleMeans = zeros(1, length(csvFiles));
Variances = zeros(1, length(csvFiles));
STD = zeros(1, length(csvFiles));

for fileIndex = 1:length(csvFiles)
    
    % Load Data_Stock_Features from the current CSV file
    Stock_Features = detectImportOptions(csvFiles{fileIndex});
    Data_Stock_Features = readtable(csvFiles{fileIndex}, Stock_Features);
    Data_Stock_Features.Date = datetime(Data_Stock_Features.Date);

    % Filter out dates after June 2011
    Data_Stock_Features = Data_Stock_Features(Data_Stock_Features.Date <= datetime(2011, 6, 30), :);

    % Extract years and months
    [years, months] = ymd(Data_Stock_Features.Date);
    uniqueYearMonths = unique(years * 100 + months);
    numMonths = length(uniqueYearMonths);
    
    % Preallocate arrays
    adjCloseLastDayOfMonth = zeros(1, numMonths);
    
    % Find the last adjusted close price of each month
    for i = 1:numMonths
        ym = uniqueYearMonths(i);
        ind = find(years * 100 + months == ym);
        adjCloseLastDayOfMonth(i) = Data_Stock_Features.AdjClose(ind(end));
    end
    
    % Calculate monthly returns
    monthlyReturns = (adjCloseLastDayOfMonth(2:end) ./ adjCloseLastDayOfMonth(1:end-1)) - 1;
    
    % Store results
    monthlyReturnsAllAssets{fileIndex} = monthlyReturns;
    sumMonthlyReturns(fileIndex) = sum(monthlyReturns);
    sampleMeans(fileIndex) = mean(monthlyReturns);
    Variances(fileIndex) = var(monthlyReturns, 1); %n-1 since using "var" by defult divides by "n". This corrects the bias.
    STD(fileIndex) = std(monthlyReturns, 1);
end

% Compute the sample covariance matrix
allMonthlyReturns = cell2mat(monthlyReturnsAllAssets')'; % Convert cell array to matrix
Covariances = cov(allMonthlyReturns, 1); %n-1 since using "cov" by defult divides by "n". This corrects the bias.
STD_1 = std(allMonthlyReturns, 1);

% Display monthly returns for each asset
for assetIndex = 1:length(csvFiles)
    disp(['Asset: ', csvFiles{assetIndex}]);
    disp(['Sample Mean: ', num2str(sampleMeans(assetIndex))]);
    disp(['Variance: ', num2str(Variances(assetIndex))]);
    fprintf('\n');
end

assetNames = csvFiles;

% Create table for variance
Variance_Table = table(assetNames', Variances', ...
          'VariableNames', {'assetNamess', 'Variance'});

% Display table for variance
disp('Asset Names and Variance');
disp(Variance_Table);

% Create table for covariance
Table_names = ['Asset', csvFiles];
Covariance_data = [csvFiles', num2cell(Covariances)];
Covariances_Table = cell2table(Covariance_data, 'VariableNames', Table_names);

% Display the table for covariances
disp('Asset Names and Covariance');
disp(Covariances_Table);


% MVO 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In the provided MVO, need to calculate lambda. To do so we need the market
%risk free rate, and variance of the market. The approach used was the CAPM
%model to calculate the expected return of the market.

% %6 month US treasury bill rates which act as the risk free rate. 
%Spans 2008 to 2011. For more accuracy, the 6-month US treasury bill was
%used which was found online.
six_month_rates = [1.71, 0.29, 0.20, 0.10];

% Market capitalization for all 20 assets in billions from 2008 to 2011.
% For more accuracy the monthly market cap should have been used since the
% data only spans till June-2011 however that data could not be found
% online for free as it required payment.
market_cap = [35.71, 45.22, 60.13, 80.01, 136.41, 96.81, 205.07, 86.67, 136.55, 143.54, 235.39, 170.12, 143.21, 173.16, 374.68, 22.07, 14.01, 171.56, 103.41, 24.98];

% Load the monthly returns data from the variable
returns = allMonthlyReturns;

% Calculate the risk-free rate (rf) by taking the average of the six_month_rates
rf = mean(six_month_rates);

% Calculate total market capitalization
total_market_cap = sum(market_cap);

% Calculate market weights (w_mkt) using market capitalization
x_market = (market_cap / total_market_cap)'; %vector of 20 x 1

% Calculate monthly excess returns for each asset by subtracting the risk-free rate
excess_returns = returns - rf/(12*100);

% Calculate the expected monthly excess returns (mu_i) for each asset
%mu_i = mean(excess_returns, 1);  % Take mean across rows

% Calculate the covariance matrix (Q) of the excess returns for all assets
Q1 = cov(excess_returns, 1);  

% Calculate the parameter lambda
%lambda = ((w_mkt' * mu_i') - (rf/(12*100))) / (w_mkt' * Q1 * w_mkt);

% Calculate the expected return of the market portfolio
expected_market_return = sum(x_market' .* (mean(returns, 1)));

% Calculate the variance of market returns
variance_market_returns = x_market' * Q1 * x_market;

% Calculate lambda as (expected market return - risk-free rate) / variance of market returns
lambda = (expected_market_return - (rf/(12*100))) / variance_market_returns;

%disp(lambda)
disp(lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solving MVO
n = 20;
Q = Covariances;
mu = sampleMeans;
c = -mu'; 
Aeq = [ones(1,n)]; %e'x = 1
beq = [1];
lb = [-inf(1, n)]; % with short selling
ra = lambda;
H = 2 * ra * Q;

x = quadprog(H, c, Aeq, beq,[],[], lb, []);

std_devi = sqrt(x' * Q * x); % standard deviation = sqrt(x'*Q*x)

variance_values = std_devi;

% Optimal weights
MVO_Weight_Table = table(assetNames', x, ...
          'VariableNames', {'Assets', 'Weights'});

% Display table for variance
disp('Asset Names and Optimal weights for MVO');
disp(MVO_Weight_Table);


disp('Variance of the portfolio:');
disp(variance_values);


expected_return_MVO = (mu * x);

disp('Expected return of MVO =');
disp(expected_return_MVO);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solving robust MVO
Q = Covariances;
mu = sampleMeans;
T = size(allMonthlyReturns, 1); %number of data points for the 20 assets. 

%Error_sigma = sigma / sqrt(dp);
ep = 1.96; %95% CL

Theta = (1/T)*diag(diag(Q));

delta = ep *sqrt(diag(Theta));

c = -mu'; 
Aeq = [ones(1,n)]; %e'x = 1
beq = [1];
lb = []; % with short selling
ra = lambda;
H = 2 * ra * Q;

c_1 = [c; delta];
H_1 = blkdiag(H, zeros(n, n)); % Extend H to account for artificial variables

% Modify the equality constraints to include artificial variables
Aeq_new = [Aeq, zeros(1, n)];
beq_new = beq;

% Inequality constraints for box uncertainty
A_new = [eye(n), -eye(n); -eye(n), -eye(n)];
b_new = zeros(2 * n, 1);

% Solve the robust quadratic programming problem
x_new = quadprog(H_1, c_1, A_new, b_new, Aeq_new, beq_new, lb, []);

% Extract the optimal weights for the original variables
x_RMVO = x_new(1:n);

% Display the result
disp('Optimal weights with robust MVO (box uncertainty) with 95% CL:');
disp(x_RMVO);

std_devi_RMVO = sqrt(x_RMVO' * Q * x_RMVO); % standard deviation = sqrt(x'*Q*x)


expected_return_RMVO = (mu * x_RMVO);

disp('Expected return of RMVO with 95% CL =');
disp(expected_return_RMVO);


%For 90% CL level
ep_1 = 1.645; %90% CL
delta_1 = ep_1 * sqrt(diag(Theta));
c_new = -mu'; 
c_new_2 = [c_new; delta_1];
x_new_1 = quadprog(H_1, c_new_2, A_new, b_new, Aeq_new, beq_new, lb, []);
% Extract the optimal weights for the original variables
x_RMVO_1 = x_new_1(1:n);
% Display the result
disp('Optimal weights with robust MVO (box uncertainty) with 90% CL:');
disp(x_RMVO_1);


std_devi_RMVO_1 = sqrt(x_RMVO_1' * Q * x_RMVO_1); % standard deviation = sqrt(x'*Q*x)

expected_return_RMVO_1 = (mu * x_RMVO_1);

disp('Expected return of RMVO with 90% CL =');
disp(expected_return_RMVO_1);















