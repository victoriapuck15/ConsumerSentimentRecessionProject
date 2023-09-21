%-------------------------------------------------------------------------
% ESTIMATION CODE FOR HOUSEHOLDS' MACROECONOMIC CONDITIONS INDEX (H-MCI)
% 
% H-MCI summarizes the information contained in various economic 
% aggregates, households’ disposable income, labor market indicators, 
% asset prices, interest rates, and external environment indicators. 
% These variables form a comprehensive set of indicators reflecting 
% the overall macroeconomic conditions faced by households. 
% The index database provide four version of the index for 22 OECD 
% countries with quaterly frequency from 2002 to 2018. 
%
% We provide the code and the resulting time series in good faight. 
% The views expressed are those of the authors and do not necessarily 
% reflect the official views of the CNB.
%
% The use of the resulting time series is free. However, we kindly ask 
% the user to reference the following methodological paper:
%   Hodula, M., S. Malovaná, and J. Frait (2019). Introducing a New Index 
%   of Households’ Macroeconomic Conditions. Czech National Bank Working 
%   Paper No. 10/2019.
%
% .........................................................................
% The code and the dataset are maintained by the authors. It is an ongoing 
% research and we welcome suggestions to enhance coverage or 
% the methodlogical approach. For questions and feedback please 
% contact simona.malovana@cnb.cz or martin.hodula@cnb.cz. 
% .........................................................................
%
% File:        DFM.m
%
% Description: The main file for CNB WP 10/2019
%              (1) Prepares data to be used in estimation.
%              (2) Runs the dynamic factor model, estimates H-MCI and 
%                  saves output.
%              (3) Plots figures.
%
% Dynamic factor model; state space representation:
%   y = C*Z + e (measurement)
%   Z = A*Z(-1) + v (transition)
%-------------------------------------------------------------------------

% IRIS required for seasonal adjustment (if indicated, x13) and time series
addpath K:\Dynare-IrisToolBox\IRIS-Toolbox-master
iris.startup

clear
clc

warning('off','all')
RootDir = cd;
FunDir  = [RootDir,'\functions\'];
addpath(genpath(FunDir));

% Preliminarities
YearEnd = 2018;
YearStart = 2000;

% Load data
DataFile =  [RootDir,'\data\data.xlsx'];  
[a,b] = xlsread(DataFile);
data = a(:,1:11);
ListC = b(5:end,1); % list of countries
Countries = unique(ListC); % list of countries unique
T = size(ListC,1)/size(Countries,1); % no. of time periods 

Variables = b(1,3:13); % variable names
Transformation = b(2,3:13); % transformation indication
SA = b(3,3:13); % seasonal adjsutment indication
Reciprocal = b(4,3:13); % reciprocal value indication
Dates = datenum(b(5:end,2),'dd.m.yyyy');
 
opt_fact = 3; % no. of factors used to construct HMCI
p = 1; % no. of lags
r = size(data,2);
 
[~, ind] = ismember(ListC,Countries(1));
X = data(find(ind),:);
[tt,~] = size(X);    
 
% Create empty structure
fields = Countries;
c = cell(length(fields),1);
Res = cell2struct(c,fields);

for j = 1:length(Countries) % for each country j
    [~, ind] = ismember(ListC,Countries(j));
    Xnsa = data(find(ind),:);
    
    % Seasonal adjustment
    X = repmat(NaN, size(Xnsa));
    [~, ind] = ismember(SA,'y');  
    X(:,find(ind)) = x13(tseries(qq(YearStart,1):qq(YearEnd,4), ...
        Xnsa(:,find(ind))));
    X(:,find(abs(ind-repmat(1,size(ind))))) = ...
        Xnsa(:,find(abs(ind-repmat(1,size(ind)))));
    
    % Transformation to growth rates if indicated
    [~, ind] = ismember(Transformation,'gr');   
    X(:,find(ind)) = [repmat(NaN,4,length(find(ind)));...
        (X(5:end,find(ind))./X(1:end-4,find(ind))-1)*100];
     
    % Reciprocal values
    [~, ind] = ismember(Reciprocal,'y');   
    X(:,find(ind)) = X(:,find(ind))*(-1);
    
    % NaNs out
    X(any(isnan(X), 2), :) = [];
    YearStart0 = YearEnd-floor(size(X,1)/4);
    if size(X,1)/4 > floor(size(X,1)/4)
        QuarterStart0 = 4-(size(X,1)/4-floor(size(X,1)/4))/0.25+1;
    else
        QuarterStart0 = 1;
        YearStart0 = YearStart0+1;
    end
   
    % Initial Conditions
    X = tseries(qq(YearStart0,QuarterStart0):qq(YearEnd,4), X);
    [~, ind] = ismember(Transformation,'gap');
    if sum(ind)>0 % if transformation to gaps required
        [A, C, Q, R, Z_0, V_0] = InitCondA(X,r,p,Transformation);
    else % if transformation to growth rates required
        [A, C, Q, R, Z_0, V_0] = InitCondB(X,r,p);
    end    

    % Standardise
    X = X(:,1:end);
    [T,N] = size(X);
    Mx = nanmean(X);
    Wx = (nanstd(X));
    X_st = (X-repmat(Mx,T,1))./repmat(Wx,T,1);
    X_st = tseries(qq(YearStart0,QuarterStart0):qq(YearEnd,4), X_st);    
    y = X_st(:,1:end)';

    % Run Kalman filter
    [Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(y, A, C, Q, R, Z_0, V_0);
    x_sm = Zsmooth(:,2:end)'*C';
    
    % Standard errors
    for i = 2:size(Vsmooth,3)
        se(:,:,i-1) = C*Vsmooth(:,:,i)*C' + R'*R;
    end
    
    % Correlation matrix, eigenvalues and HMCI
    HMCI = repmat(NaN,size(Zsmooth,2)-1,1);
    [~, ind] = ismember(Transformation,'gap');
    if sum(ind)>0
        factors = Zsmooth(2*N+1:2*N+opt_fact,2:end)';
        gaps = Zsmooth(1:N,2:end)';
        tempCorr = corr(factors,gaps);        
        for ij = 1:opt_fact
            maxAbsCorr = sort(abs(tempCorr(ij,:)),'descend');
            if maxAbsCorr(1) > max(tempCorr(ij,:))
                factors(:,ij) = -1*factors(:,ij);
            end
        end
        for ii = 1:(size(Vsmooth,3)-1)
            eigen_real = eig(corrcoef(Q(2*N+1:end,2*N+1:end)));
            eigen = sort(eigen_real(end-opt_fact+1:end),'descend');            
            HMCI(ii) = sum(factors(ii,:).*(eigen./sum(eigen))');
        end
        HMCI = (HMCI-mean(HMCI))./std(HMCI);
    else
        factors = Zsmooth(1:opt_fact,2:end)';
        tempCorr = corr(factors,X_st(:,1:end));
        for ij = 1:opt_fact
            if max(abs(tempCorr(ij,:))) > max(tempCorr(ij,:))
                factors(:,ij) = -1*factors(:,ij);
            end
        end            
        for ii = 1:(size(Vsmooth,3)-1)
            eigen_real = eig(corrcoef(Q));
            eigen = sort(eigen_real(end-opt_fact+1:end),'descend');
            HMCI(ii) = sum(factors(ii,:).*(eigen./sum(eigen))');
        end
    end
    
    % Save results
    Res.(string(Countries(j))).eigen = eigen;
    Res.(string(Countries(j))).X_st = X_st;
    Res.(string(Countries(j))).xx_sm = x_sm;
    Res.(string(Countries(j))).X_sm = repmat(Wx,T,1).*x_sm+repmat(Mx,T,1);
    Res.(string(Countries(j))).se = se;    
    Res.(string(Countries(j))).hmci = HMCI;
    [~, ind] = ismember(Transformation,'gap');
    if sum(ind)>0  
        Res.(string(Countries(j))).gaps = Zsmooth(1:N,2:end)';
        Res.(string(Countries(j))).trends = Zsmooth(N+1:2*N,2:end)';
        Res.(string(Countries(j))).F = Zsmooth(2*N+1:end,2:end)';
    else
        Res.(string(Countries(j))).F = Zsmooth(:,2:end)';
    end
end

%% Save estimated HMCI
hmci = repmat(NaN,T,size(Countries,1));
for j = 1:length(Countries)
    hmci(:,j) = [repmat(NaN,76-size(Res.(string(Countries(j))).hmci,1),1); ...
         Res.(string(Countries(j))).hmci];
end
xlswrite('hmci.xlsx',hmci,'hmci','B2')
xlswrite('hmci.xlsx',Countries','hmci','B1')
xlswrite('hmci.xlsx',b(5:4+T,2),'hmci','A2')

%% Load CCI (OECD)
DataFile =  [RootDir,'\data\cci'];  
[cci,b] = xlsread(DataFile);
cciList = b(1,2:end);

% Figures
figure(1)
for j = 1:length(Countries)
    subplot(8,3,j)
    xx = Res.(string(Countries(j))).hmci;
    [~, ind] = ismember(cciList,Countries(j));
    yy = cci(:,find(ind));
    YearStart0 = YearEnd-floor(size(xx,1)/4);
    if size(xx,1)/4 > floor(size(xx,1)/4)
        QuarterStart0 = 4-(size(xx,1)/4-floor(size(xx,1)/4))/0.25+1;
    else
        QuarterStart0 = 1;
        YearStart0 = YearStart0+1;
    end
    plot(tseries(qq(YearStart0,QuarterStart0):qq(YearEnd,4), ...
        xx))
    hold on
    if size(yy,2)>0
        YearEnd0 = 2018;
        YearStart0 = YearEnd0-floor(size(yy,1)/4);
        if size(yy,1)/4 > floor(size(yy,1)/4)
            QuarterStart0 = 4-(size(yy,1)/4-floor(size(yy,1)/4))/0.25+1;
        else
            QuarterStart0 = 1;
            YearStart0 = YearStart0+1;
        end
        plot(tseries(qq(YearStart0,QuarterStart0):qq(YearEnd0,4), ...
            yy))
    end
    title(string(Countries(j)))
end
legend('HMCI (growth)','CCI (OECD)')
