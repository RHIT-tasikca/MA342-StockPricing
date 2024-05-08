close all hidden
clear all
tic
% tickers = ["AAPL",'BILI','CLSK','COIN','DIS','F','NKE','NVDA','PDD','RIOT','SE','TROX','TSLA','WMT','CROX',...
%     'TGT','NTDOY','MCD','VTI','HON','GOOG','MSFT','EBAY','MAR','SBUX','UBER','FUJHY','QCOM','DUOL','HLT',...
%     'AMZN','AVGO','META','COST','TXN','RRGB','GS','JPM','PINS','BIG','LOW','HD','WEN','APD','WBD'];
tickers = ["AAPL","F","NVDA","DIS","TSLA"];

for i = 1:size(tickers,2)
    P = hist_stock_data('30042023','30042024', tickers(i));
    Price(:,i) = P.AdjClose(1:end);
end

Return(:,:) = (Price(1:(size(Price,1)-1),:) - Price(2:end,:))./Price(2:end,:);

xbar(1,:) = mean(Return);

CR = cov(Return);
CP = cov(Price);

numRowPrice = size(Price,1);
numDayToPredict = 30*4;
delt = 1;
for i = numRowPrice:numRowPrice+numDayToPredict
    for j = 1:size(tickers,2)

        Price(i+1,j) = Price(i,j) + (xbar(1,j)*Price(i,j)*delt) + (sqrt(CR(j,j))*Price(i,j)*sqrt(delt)*randn(1));
        if Price(i+1,j) <= 0.001
            Price(i+1,j) = 10^-8;
        end
    end
end

hold on
figure(1)
for k = 1:size(tickers,2)
    plot([1:1:size(Price,1)],Price(:,k))
end
%legend(tickers)
xlabel('Time (days)')
ylabel('Price ($)')
grid on
hold off

clear Return
Return(:,:) = (Price(1:(size(Price,1)-1),:) - Price(2:end,:))./Price(2:end,:);

alpha = 0.95;
n = size(tickers,2);
dayInMonth = 30;
numMonthToPredict = numDayToPredict/dayInMonth;
for i = 1:numMonthToPredict
    MonthsBefore = size(Return,1)-((numMonthToPredict-i)*dayInMonth)-(6*dayInMonth);
    pastDates(i) = MonthsBefore;
    currentDate = size(Return,1)-((numMonthToPredict-i)*dayInMonth);
    curDates(i) = currentDate;
    xbar(i+1,:) = mean(Return(MonthsBefore:currentDate,:))

    w(:,i) = quadprog((1-alpha)*2*cov(Return(MonthsBefore:currentDate,:)), -alpha*mean(Return(MonthsBefore:currentDate,:)), [], [], ones(1,n), ...
    [1], zeros(n,1), 0.2*ones(n,1))
end

startMoney = 10000;
lastMoney = startMoney;


for i = 1:size(w,2)
    moneyPerShare = lastMoney*w(:,i)';
    Price(pastDates(i),:);
    numShares = moneyPerShare./Price(pastDates(i),:);
    Price(curDates(i),:)
    afterMoney = Price(curDates(i),:).*numShares;
    curMoney = sum(afterMoney);
    monthReturn(i) = curMoney - lastMoney;
    lastMoney = curMoney
end

lastMoney

toc
