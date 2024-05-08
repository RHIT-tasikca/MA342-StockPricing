close all hidden
clear all
tic
tickers = ["AAPL",'BILI','CLSK','COIN','DIS','F','NKE','NVDA','PDD','RIOT','SE','TROX','TSLA','WMT','CROX',...
    'TGT','NTDOY','MCD','VTI','HON','GOOG','MSFT','EBAY','MAR','SBUX','UBER','FUJHY','QCOM','DUOL','HLT',...
    'AMZN','AVGO','META','COST','TXN','RRGB','GS','JPM','PINS','BIG','LOW','HD','WEN','APD','WBD'];

for i = 1:size(tickers,2)
    P = hist_stock_data('09032023','09092023', tickers(i));
    Price(:,i) = P.AdjClose(1:end);
    Return(:,i) = (P.AdjClose(1:(size(P.AdjClose,1)-1),1) - P.AdjClose(2:end,1))./P.AdjClose(2:end,1);
end
xbar = mean(Return);

CR = cov(Return);
CP = cov(Price);

furPrice(1,:) = Price(end,:);
numDayToPredict = 250;
delt = 1;
for i = 1:numDayToPredict
    for j = 1:size(tickers,2)
        furPrice(i+1,j) = furPrice(i,j) + (xbar(j)*furPrice(i,j)*delt) + (sqrt(CR(j,j))*furPrice(i,j)*sqrt(delt)*randn(1));
        if furPrice(i+1,j) <= 0.001
            furPrice(i+1,j) = 10^-8;
        end
    end
end

hold on
figure(1)
for k = 1:size(tickers,2)
    plot([-(size(Price,1)-1):0],Price(:,k))
    plot([0:numDayToPredict],furPrice(:,k))
end
legend(tickers)
grid on
hold off

for i= 2:size(furPrice,1)  %skip first row b/c that is intial Price
    R(i,:) = (furPrice(i,:) - furPrice(i-1,:))./furPrice(i-1,:);
end

figure(2)
hold on
counter = 1;
for alpha =  [0:0.01:1].^2
n = size(tickers,2);
[ww(:,counter), optVal] = quadprog((1-alpha)*2*CR, -alpha*xbar, [], [], ones(1,n), ...
    [1], zeros(n,1), 0.2*ones(n,1));

plot(ww(:,counter)'*CR*ww(:,counter), xbar*ww(:,counter),".")
xlabel('Risk')
ylabel('Expected Return')
counter = counter + 1;
end

for i = 1:1:length(ww)
    x_vec(i) = ww(:,i)'*CR*ww(:,i);
    y_vec(i) = xbar*ww(:,i);
end

alpha_vec = [0:0.01:1].^2;

for i = 1:1:(length(x_vec)-1)
    slopes(i) = (y_vec(i+1) - y_vec(i))/(x_vec(i+1)-x_vec(i));
end

for i = length(slopes):-1:1
    if (slopes(i)-1) > 0
        optimum_alpha_index = i;
        break
    end
end

plot(ww(:,72)'*CR*ww(:,72), xbar*ww(:,72),'diamond','MarkerSize',10)
axis equal

hold off

toc
