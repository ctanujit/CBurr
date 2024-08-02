% matlab plot actual and estimated frequencies of the proposed CBurr and Other distributions  

clear all;


data=csvread('compared_output_bio-SC-HT_CBurr.csv',1);
xDeg=data(:,1);
xAct=data(:,2);
xburr=data(:,3);
xburr_exponiented=data(:,7);
xburr_mo=data(:,8);
xCburr=data(:,4);
xPow=data(:,11);
xPar=data(:,12);
xLog=data(:,10);
xPoC=data(:,6);
xLomax=data(:,9);
%xExp=data(:,10);
%xPoi=data(:,11);



figure

loglog(xDeg,xAct,'.','MarkerSize',7,'MarkerEdgeColor','b')
%semilogy(xDeg,xAct,'.','MarkerSize',7)
%plot(xDeg,xAct,'.','MarkerSize',7)

%txt1 = 'Actual';
%text(xDeg(10),xAct(10),txt1)

hold on
loglog(xDeg,xPar,'linewidth',1,'color',[0, 0.75, 0.75])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xPow,'linewidth',1,'color',[0.4940, 0.1840, 0.5560])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xLog,'linewidth',1,'color',[0.75,0,0.75])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)


hold on
loglog(xDeg,xPoC,'linewidth',1,'color',[0.75,0.75,0])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xLomax,'linewidth',1,'color',[0.25,0.25,0.25])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xburr,'linewidth',1,'color',[0, 0.4470, 0.7410])
%semilogy(xDeg,xPar,'linewidth',1.5)
%plot(xDeg,xPar,'linewidth',1.5)


hold on
loglog(xDeg,xburr_exponiented,'linewidth',1.2,'color',[0.8500,0.3250, 0.0980])
%semilogy(xDeg,xPow,'linewidth',1.5)
%plot(xDeg,xPow,'linewidth',1.5)


hold on
loglog(xDeg,xburr_mo,'linewidth',1.4,'color',[0.9290, 0.6940, 0.1250])
%semilogy(xDeg,xLog,'linewidth',1.5)
%plot(xDeg,xLog,'linewidth',1.5)

hold on
loglog(xDeg,xCburr,'linewidth',1.6,'color',[1, 0, 0])
%semilogy(xDeg,xPoC,'linewidth',1.5)
%plot(xDeg,xPoC,'linewidth',1.5)


% x1=xDeg(34);
% y1=xEst(34);
% %stem(x1,0.31,'*','MarkerSize',15,'color','red')
% %loglog(x1,y1,'square','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor',[0.5,0.5,0.5])
% loglog(x1,0.31,'*','MarkerSize',12,'MarkerEdgeColor','r','MarkerFaceColor',[0.5,0.5,0.5])
% %xticks(x1)
% %xticklabels('TD')
% %xlabel(x1,0.31,'TD');
% text(x1,0.31,'TD','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold');

ylim([0.3 1000]);
xlim([1 2000]);
set(gca,'fontweight','bold','fontsize',12);
xlabel('Node Degree');
ylabel('Frequency');

L=legend('bio-SC-HT Network','Pareto Type-I','Power-law','Log-normal','Power-law cutoff','Lomax','Burr','Exponentiated Burr','MO Burr','Proposed CBurr','Location',[0.5, 0.5, .25, .25]);

%L=legend('ego-Twitter Network','Exponential','Pareto','Power-law','Lognormal','Power-law Exp.Cutoff','Double Power-law','Location',[0.5, 0.5, .25, .25]);
%L=legend('Artist Facebook Network','Double Power-law','Power-law','Pareto','Lognormal','Power-law Exp. Cutoff','Exponential','Poisson','Location',[0.5, 0.5, .25, .25]);

%L.FontWeight='bold';
%print(gcf,'foo.png','-dpng','-r300');  

